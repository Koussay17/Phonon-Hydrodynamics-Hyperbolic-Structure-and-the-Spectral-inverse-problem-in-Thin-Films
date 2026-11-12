"""RS2b: FDTR-visible first-order time of the slabs computed by rs2_slab.py (same extraction as the audit's RTA slab).

For each tensor and operator (full, rta) and thickness (300, 500 nm):
  * two-port at the 80 FDTR frequencies: exact (freqset 'all') or AAA rational interpolation in p of each entry of Y
    and Th from the exact 'sub' frequencies (validated on REF against 'all', see rs2c_interp_check.py);
  * static: G0 = Re Y11(100 Hz), lam_eff = G0 d, R_slab = 1/G0; first-order two-port coefficients a1 = Im Y11/omega,
    b1 = Im Y12/omega (Richardson 1 and 2 kHz), delta_2port = (C d/6 - (a1 - b1)) / (2 G0) (homogeneous GK film with
    the slab's static conductance, no separate wall resistance);
  * FDTR data: 1D stack Au 80 nm / R1 / slab / R2 / sapphire (au6_slab.Z_stack_slab, port 'wall' = R1, R2 on the
    wall emission temperatures, 'energy' = on the energy temperatures at the walls), note-16 (1e-8, 2e-8) and
    favourable (2e-9, 2e-9) stacks; fitted by the twin's 1D GK film (au5_fisher1d.lnH, local convention) with
    tau_R fixed at the base value and tau_l free (delta_fit = tau_R,base - tau_l,fit), nuisances N2lam (lam_z free,
    C known) and N3 (lam_z, C free); also N2fix / N1fix (lam_z, and for N1 also R1, R2, fixed at the REF fit),
    the GK film with the regularised tau_l (no free target), Cattaneo (tau_l = 0, tau_R free) and Fourier films.
  Base point (as au8): the REF regularised film of the same operator family - full: C + QDQ kinetic times and K;
  rta: RTA V1 film times and lam_M - for every tensor.
usage: python -B rs2b_fits.py <mesh> <group> <freqset> <name> [name ...]     Output: out/rs2b_<mesh>_<group>_<freqset>.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import AAA

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
AUD = HERE.parent / "R_A_audit"
sys.path.insert(0, str(AUD))
import au5_fisher1d as a5  # noqa: E402
import au6_slab as a6  # noqa: E402

STACKS = {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}
SETS = {"N2lam": a5.SETS["N2"] + [2], "N3": a5.SETS["N3"], "N2fix": a5.SETS["N2"], "N1fix": a5.SETS["N1"]}


def wall_tag(w):
    return "black" if w[0] == "black" else f"{w[0]}{w[1]}"


def load(mesh, group, op, ftag, nm, wall=0):
    """Two-ports of one wall model (index into the stored wall list): f, d, Y[d, f, 2, 2], Th[...], diagnostics."""
    Z = np.load(OUT / f"rs2_{mesh}_{group}_{op}_{ftag}_{nm}.npz")
    Y, Th = Z["Y"], Z["Th"]
    if Y.ndim == 5:
        Y, Th = Y[wall], Th[wall]
    return Z["f"], Z["d"], Y, Th, {k: Z[k] for k in Z.files if k.startswith(("dg_", "ei_"))}


def walls_of(mesh, group, op, ftag, nm):
    Z = np.load(OUT / f"rs2_{mesh}_{group}_{op}_{ftag}_{nm}.npz")
    return [tuple(w) for w in json.loads(str(Z["walls"]))] if "walls" in Z.files else [("black",)]


def at_fdtr(f, Y, Th, ftag):
    """Y, Th at the 80 FDTR frequencies: exact for 'all', AAA in p otherwise; returns arrays [d, 80, 2, 2]."""
    F = a5.FREQ
    if ftag == "all":
        idx = [int(np.argmin(np.abs(f - x))) for x in F]
        assert np.allclose(f[idx], F, rtol=1e-12)
        return Y[:, idx], Th[:, idx], 0.0
    pz = 2j * np.pi * f
    pF = 2j * np.pi * F
    Yo, To = np.zeros((Y.shape[0], F.size, 2, 2), complex), np.zeros((Y.shape[0], F.size, 2, 2), complex)
    worst = 0.0
    for i in range(Y.shape[0]):
        for a in range(2):
            for b in range(2):
                for src, dst in ((Y, Yo), (Th, To)):
                    r = AAA(pz, src[i, :, a, b], rtol=1e-14)
                    dst[i, :, a, b] = r(pF)
                    worst = max(worst, float(np.max(np.abs(r(pz) - src[i, :, a, b])) / np.max(np.abs(src[i, :, a, b]))))
    return Yo, To, worst


def static_numbers(f, Y, d, C):
    i100 = int(np.argmin(np.abs(f - 1e2)))
    i1, i2 = int(np.argmin(np.abs(f - 1e3))), int(np.argmin(np.abs(f - 2e3)))
    out = {}
    for k, dd in enumerate(d):
        G0 = float(Y[k, i100, 0, 0].real)
        w1, w2 = 2 * np.pi * f[i1], 2 * np.pi * f[i2]
        a1 = (4 * Y[k, i1, 0, 0].imag / w1 - Y[k, i2, 0, 0].imag / w2) / 3
        b1 = (4 * Y[k, i1, 0, 1].imag / w1 - Y[k, i2, 0, 1].imag / w2) / 3
        out[f"{dd * 1e9:.0f}"] = {"G0": G0, "R_slab": 1 / G0, "lam_eff": G0 * dd, "a1": float(a1), "b1": float(b1),
                                  "a1_plus_b1_over_Cd2": float((a1 + b1) / (C * dd / 2)),
                                  "delta_2port_ps": float((C * dd / 6 - (a1 - b1)) / (2 * G0) * 1e12)}
    return out


def refit_free_taul(y, bb, nuis, sp):
    """GK-film fit with tau_l free from four starting points delta_start = f * delta_base (f = 1, 0.4, 0.2, 0);
    returns the lowest misfit, its parameters, the base it was started from, and the spread of the converged deltas
    among starts that reach the minimum to within 1e-6 (identifiability check)."""
    d0 = bb["tauR"] - bb["taul"]
    runs = []
    for frac in (1.0, 0.4, 0.2, 0.0):
        b2 = dict(bb, taul=bb["tauR"] - frac * d0)
        chi, th = a5.gn_refit(y, b2, nuis, sp, free=(1,))
        runs.append((chi, th, b2, (b2["tauR"] - b2["taul"] - th[1] * 1e-12) * 1e12))
    chi, th, b2, dl = min(runs, key=lambda r: r[0])
    same = [r[3] for r in runs if r[0] <= chi + 1e-6]
    return chi, th, b2, {"n_starts_at_min": len(same), "delta_spread_ps": float(max(same) - min(same)),
                         "worst_start_dchi2": float(max(r[0] for r in runs))}


def _absolute(b2, th):
    return {"taul": b2["taul"] + th[1] * 1e-12, "lamz": b2["lamz"] * float(np.exp(th[2])),
            "C": b2["C"] * float(np.exp(th[3])), "R1": b2["R1"] * float(np.exp(th[4])), "R2": b2["R2"] * float(np.exp(th[5]))}


def fits_one(Ys, Ths, d, base, port, stack, sets=("N2lam", "N3"), sps=(0.01, 0.1), extra=True, fix=None):
    """GK-film fits of the slab data (tau_l free) for each nuisance set; N2lam and N3 are also restarted from each
    other's converged point (warm start; the set's fixed parameters keep their base values)."""
    p = 2j * np.pi * a5.FREQ
    R1, R2 = STACKS[stack]
    y = np.log(a6.Z_stack_slab(list(Ys), list(Ths), p, R1, R2, port=port))
    b0 = dict(base, R1=R1, R2=R2, d=d)
    best = {}
    for N in sets:
        for sp in sps:
            bb = dict(b0)
            if N in ("N2fix", "N1fix"):
                if fix is None:
                    continue
                f_ = fix[f"{port}|{stack}|{sp}"]
                bb["lamz"] = f_["lamz"]
                if N == "N1fix":
                    bb["R1"], bb["R2"] = f_["R1"], f_["R2"]
            chi, th, b2, ms = refit_free_taul(y, bb, SETS[N], sp)
            best[(N, sp)] = [chi, th, b2, ms, bb, "base"]
    for N, other in (("N2lam", "N3"), ("N3", "N2lam")):
        for sp in sps:
            if (N, sp) not in best or (other, sp) not in best:
                continue
            sol = _absolute(best[(other, sp)][2], best[(other, sp)][1])
            bb = best[(N, sp)][4]
            bw = dict(bb, taul=sol["taul"], lamz=sol["lamz"], R1=sol["R1"], R2=sol["R2"])
            if N == "N3":
                bw["C"] = sol["C"]
            chi, th, b2, ms = refit_free_taul(y, bw, SETS[N], sp)
            if chi < best[(N, sp)][0] - 1e-9:
                best[(N, sp)] = [chi, th, b2, ms, bb, f"warm from {other}"]
    out = {}
    for (N, sp), (chi, th, b2, ms, bb, origin) in best.items():
        a = _absolute(b2, th)
        ent = {"delta_fit_ps": (b2["tauR"] - a["taul"]) * 1e12, "dchi2": chi, "multistart": dict(ms, origin=origin),
               "lamz_fit": a["lamz"], "C_fit": a["C"], "R1_fit": a["R1"], "R2_fit": a["R2"]}
        if extra and N in ("N2lam", "N3"):
            chiA, _ = a5.gn_refit(y, bb, SETS[N], sp)                       # GK with the regularised tau_l
            bC = dict(bb, taul=0.0)
            chiC, thC = a5.gn_refit(y, bC, SETS[N], sp, free=(0,))           # Cattaneo, tau_R free
            bF = dict(bb, tauR=0.0, taul=0.0)
            chiF, _ = a5.gn_refit(y, bF, SETS[N], sp)                       # Fourier
            ent.update({"dchi2_gk_taul_reg": chiA, "dchi2_cattaneo": chiC,
                        "tauR_fit_cattaneo_ps": (bb["tauR"] + thC[0] * 1e-12) * 1e12, "dchi2_fourier": chiF})
        out[f"{N}|{sp}"] = ent
    return out, y


def main():
    mesh, group, ftag = sys.argv[1], sys.argv[2], sys.argv[3]
    names = sys.argv[4:]
    t0 = time.perf_counter()
    resp = OUT / f"rs2b_{mesh}_{group}_{ftag}.json"
    res = json.loads(resp.read_text()) if resp.exists() else {}
    stat = {nm: json.loads((OUT / f"rs2_{mesh}_{group}_static_{nm}.json").read_text()) for nm in set(names) | {"REF"}}
    C = stat["REF"]["C_tot"]
    for op in ("full", "rta"):
        if not (OUT / f"rs2_{mesh}_{group}_{op}_{ftag}_REF.npz").exists():
            continue
        walls = walls_of(mesh, group, op, ftag, "REF")
        fixes = {}
        for nm in ["REF"] + [x for x in names if x != "REF"]:
            if not (OUT / f"rs2_{mesh}_{group}_{op}_{ftag}_{nm}.npz").exists():
                continue
            ent = res.get(op, {}).get(nm, {})
            for iw, w in enumerate(walls):
                wt = wall_tag(w)
                sfx = "" if wt == "black" else f"|{wt}"
                f, d, Y, Th, dg = load(mesh, group, op, ftag, nm, iw)
                Yf, Tf, aaa_err = at_fdtr(f, Y, Th, ftag)
                ent["static" + sfx] = static_numbers(f, Y, d, stat[nm]["C_tot"])
                ent["aaa_support_residual" + sfx] = aaa_err
                ent["max_energy_balance" + sfx] = float(dg["dg_energy_balance_rel"][iw].max())                     if dg["dg_energy_balance_rel"].ndim == 3 else float(dg["dg_energy_balance_rel"].max())
                ent["max_eig_backward"] = float(dg["ei_eig_backward_max"].max())
                for k, dd in enumerate(d):
                    L = f"{dd * 1e9:.0f}"
                    sr = stat["REF"]
                    if op == "full":
                        kr = sr[f"reg_full_QDQ_{L}"]
                        base = {"tauR": kr["tauR"], "taul": kr["taul"], "lamz": kr["K"], "C": C}
                        kw = stat[nm][f"reg_full_QDQ_{L}"]
                        ent[f"delta_reg_{L}_ps"] = kw["delta"] * 1e12
                        ent[f"delta_reg_phys_{L}_ps"] = stat[nm][f"reg_full_phys_{L}"]["delta"] * 1e12
                        ent[f"kappa_reg_{L}"] = kw["K"]
                    else:
                        kr = sr[f"reg_rta_V1_{L}"]
                        base = {"tauR": kr["tauR"], "taul": kr["taul"], "lamz": kr["lam"], "C": C}
                        kw = stat[nm][f"reg_rta_V1_{L}"]
                        ent[f"delta_reg_{L}_ps"] = kw["delta"] * 1e12
                        ent[f"kappa_reg_{L}"] = kw["lam"]
                    ent[f"base_{L}"] = {kk: (vv * 1e12 if kk in ("tauR", "taul") else vv) for kk, vv in base.items()}
                    fk = (L, wt)
                    if nm == "REF" and wt == "black":
                        for port in ("wall", "energy"):
                            for stack in STACKS:
                                r_, _ = fits_one(Yf[k], Tf[k], dd, base, port, stack)
                                ent[f"fit_{L}|{port}|{stack}"] = r_
                                for sp in (0.01, 0.1):
                                    fixes.setdefault(fk, {})[f"{port}|{stack}|{sp}"] = {
                                        "lamz": r_[f"N2lam|{sp}"]["lamz_fit"], "R1": r_[f"N2lam|{sp}"]["R1_fit"],
                                        "R2": r_[f"N2lam|{sp}"]["R2_fit"]}
                                r2, _ = fits_one(Yf[k], Tf[k], dd, base, port, stack, sets=("N2fix", "N1fix"),
                                                 extra=False, fix=fixes[fk])
                                ent[f"fit_{L}|{port}|{stack}"].update(r2)
                        sens = {}
                        for fac in (0.5, 2.0):
                            b2 = dict(base, tauR=base["tauR"] * fac)
                            r_, _ = fits_one(Yf[k], Tf[k], dd, b2, "wall", "note16", sps=(0.01,), extra=False)
                            sens[str(fac)] = {kk: vv["delta_fit_ps"] for kk, vv in r_.items()}
                        ent[f"tauR_sensitivity_{L}"] = sens
                    else:
                        if nm == "REF":
                            r_, _ = fits_one(Yf[k], Tf[k], dd, base, "wall", "note16", extra=True)
                            for sp in (0.01, 0.1):
                                fixes.setdefault(fk, {})[f"wall|note16|{sp}"] = {
                                    "lamz": r_[f"N2lam|{sp}"]["lamz_fit"], "R1": r_[f"N2lam|{sp}"]["R1_fit"],
                                    "R2": r_[f"N2lam|{sp}"]["R2_fit"]}
                        else:
                            r_, _ = fits_one(Yf[k], Tf[k], dd, base, "wall", "note16", extra=False)
                        r2, _ = fits_one(Yf[k], Tf[k], dd, base, "wall", "note16", sets=("N2fix", "N1fix"), extra=False,
                                         fix=fixes[fk])
                        r_.update(r2)
                        ent[f"fit_{L}|wall|note16" + sfx] = r_
                        if wt == "black":
                            r3, _ = fits_one(Yf[k], Tf[k], dd, base, "energy", "note16", sps=(0.01,), extra=False)
                            ent[f"fit_{L}|energy|note16"] = r3
            res.setdefault(op, {})[nm] = ent
            fr = ent["fit_500|wall|note16"]
            print(op, nm, "delta_reg 300/500", round(ent["delta_reg_300_ps"], 2), round(ent["delta_reg_500_ps"], 2),
                  "| fit 500 N2lam/N3", round(fr["N2lam|0.01"]["delta_fit_ps"], 2), round(fr["N3|0.01"]["delta_fit_ps"], 2),
                  "chi", round(fr["N2lam|0.01"]["dchi2"], 4), "| fit 300 N3",
                  round(ent["fit_300|wall|note16"]["N3|0.01"]["delta_fit_ps"], 2),
                  "| walls 500 N3:", {wall_tag(w): round(ent[f"fit_500|wall|note16" + ("" if wall_tag(w) == "black" else "|" + wall_tag(w))]["N3|0.01"]["delta_fit_ps"], 2) for w in walls},
                  "| aaa", f"{ent['aaa_support_residual']:.1e}", f"[{time.perf_counter() - t0:.0f}s]", flush=True)
            resp.write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
