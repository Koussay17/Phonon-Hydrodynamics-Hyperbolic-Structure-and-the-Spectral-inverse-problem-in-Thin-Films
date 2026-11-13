"""SA2: audit of the extraction of the FDTR-visible delta and of the resolvability margins (7x7x5, full operator,
black walls, wall port, note-16 stack), with the independent FDTR side sa_fdtr.py.

Stages:
  ref     REF data from the stored exact 80-frequency two-ports (rs2_m775_film2_full_all_REF.npz): multistart GK fits
          (N2lam, N3; 0.01, 0.1 deg), misfit profile over delta (global-minimum check), frequency-subset and band
          sensitivity; cross-check of the stack against au6_slab.Z_stack_slab (check only).
  interp  interpolation check on REF: AAA (scipy) from the sub4 frequencies vs exact; ratio-spline interpolation
          (Y_W / Y_REF,exact splined in ln f) validated on the exact REF wall variants.
  wit     REF + 12 witnesses at 300 and 500 nm (sub4 two-ports): delta_fit with AAA and with the ratio-spline
          interpolation; half-ranges (campaign rule); Fisher sigma(delta) (Cauchy-integral derivatives) at the REF
          slab-fit base point; margins; direct test m = 1/sqrt(max Delta chi^2).
Output: out/sa2_fits.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
from scipy.interpolate import AAA, CubicSpline  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sa_fdtr as fd  # noqa: E402

RS = HERE.parent / "R_slab" / "out"
RES = HERE / "out" / "sa2_fits.json"
T0 = time.perf_counter()
F = fd.FREQ
P = 2j * np.pi * F
R1N, R2N = 1e-8, 2e-8
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
STAT = json.loads((RS / "rs2_m775_film2_static_REF.json").read_text())


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


def save(res):
    RES.write_text(json.dumps(res, indent=1, default=float))


def base_par(L):
    kr = STAT[f"reg_full_QDQ_{L}"]
    return np.array([kr["tauR"], kr["taul"], np.log(kr["K"]), np.log(STAT["C_tot"]), np.log(R1N), np.log(R2N)])


def load_exact(op="full", wall=0):
    Z = np.load(RS / f"rs2_m775_film2_{op}_all_REF.npz")
    f, d, Y, Th = Z["f"], Z["d"], Z["Y"], Z["Th"]
    if Y.ndim == 5:
        Y, Th = Y[wall], Th[wall]
    idx = [int(np.argmin(np.abs(f - x))) for x in F]
    assert np.allclose(f[idx], F, rtol=1e-12)
    return d, Y[:, idx], Th[:, idx]


def load_sub(nm, op="full", wall=0, ftag="sub4"):
    Z = np.load(RS / f"rs2_m775_film2_{op}_{ftag}_{nm}.npz")
    f, d, Y, Th = Z["f"], Z["d"], Z["Y"], Z["Th"]
    if Y.ndim == 5:
        Y, Th = Y[wall], Th[wall]
    return f, d, Y, Th


def aaa_interp(fsub, Ysub):
    """AAA in p of each entry (support = all stored frequencies) evaluated at the 80 FDTR frequencies."""
    out = np.zeros((F.size, 2, 2), complex)
    for a in range(2):
        for b in range(2):
            out[:, a, b] = AAA(2j * np.pi * fsub, Ysub[:, a, b], rtol=1e-14)(P)
    return out


def ratio_spline(fsub, Ysub, Yref_sub, Yref_exact):
    """Y_W(f) = Y_REF,exact(f) * spline_{ln f}(Y_W/Y_REF) built on the in-band stored frequencies."""
    m = fsub >= 9e3
    x = np.log(fsub[m])
    out = np.zeros((F.size, 2, 2), complex)
    for a in range(2):
        for b in range(2):
            r = Ysub[m, a, b] / Yref_sub[m, a, b]
            sr = CubicSpline(x, r.real)(np.log(F))
            si = CubicSpline(x, r.imag)(np.log(F))
            out[:, a, b] = Yref_exact[:, a, b] * (sr + 1j * si)
    return out


def data(Y):
    return np.log(fd.Z_slab(Y, R1N, R2N, P, port="wall"))


def stage_ref(res):
    d, Ye, The = load_exact()
    out = {}
    # cross-check of the stack against the audited code (check only)
    sys.path.insert(0, str(HERE.parent / "R_A_audit"))
    import au6_slab as a6
    ychk = np.log(a6.Z_stack_slab(list(Ye[1]), list(The[1]), P, R1N, R2N, port="wall"))
    out["stack_crosscheck_max_abs_dlnH"] = float(np.max(np.abs(ychk - data(Ye[1]))))
    log("stack cross-check", out["stack_crosscheck_max_abs_dlnH"])
    for i, L in enumerate(("300", "500")):
        dd = float(L) * 1e-9
        y = data(Ye[i])
        b = base_par(L)
        ent = {}
        for N in ("N2lam", "N3"):
            for sp in (0.01, 0.1):
                r = fd.fit_delta(y, dd, b, sp, N)
                ent[f"{N}|{sp}"] = r
                log(L, N, sp, "delta", round(r["delta_ps"], 4), "chi2", f"{r['chi2']:.3e}", "starts at min",
                    r["n_at_min"], "/", r["n_starts"], "spread", f"{r['spread_ps']:.1e}", "worst", f"{r['worst_chi2']:.2e}")
        # profile over delta (0.01 deg)
        for N in ("N2lam", "N3"):
            grid = np.concatenate([np.arange(-400, -150, 10.0), np.arange(-150, 60, 2.5), np.arange(60, 160, 10.0)])
            prof = fd.profile_delta(y, dd, b, 0.01, N, grid, warm=ent[f"{N}|0.01"]["par"])
            k = int(np.argmin(prof))
            # local minima of the profile
            locmin = [float(grid[j]) for j in range(1, grid.size - 1) if prof[j] < prof[j - 1] and prof[j] < prof[j + 1]]
            # curvature near the minimum -> sigma_profile (delta direction = -tau_l at fixed tau_R)
            fine = ent[f"{N}|0.01"]["delta_ps"] + np.linspace(-6, 6, 13)
            pf = fd.profile_delta(y, dd, b, 0.01, N, fine, warm=ent[f"{N}|0.01"]["par"])
            cq = np.polyfit(fine - ent[f"{N}|0.01"]["delta_ps"], pf, 2)
            ent[f"profile|{N}|0.01"] = {"grid_ps": grid.tolist(), "chi2": prof.tolist(), "argmin_grid_ps": float(grid[k]),
                                        "min_grid_chi2": float(prof[k]), "local_minima_ps": locmin,
                                        "chi2_at_ends": [float(prof[0]), float(prof[-1])],
                                        "sigma_profile_ps": float(1 / np.sqrt(cq[0])), "fine_chi2": pf.tolist()}
            log(L, "profile", N, "argmin", grid[k], "local minima", locmin, "sigma_prof", round(1 / np.sqrt(cq[0]), 3),
                "chi2 at -400/+150", round(prof[0], 1), round(prof[-1], 1))
        # frequency subsets and bands (0.01 deg)
        sub = {}
        for tag, sel in (("even40", np.arange(0, 80, 2)), ("odd40", np.arange(1, 80, 2)), ("every4_20", np.arange(0, 80, 4))):
            for N in ("N2lam", "N3"):
                r = fd.fit_delta(y[sel], dd, b, 0.01, N, freq=F[sel])
                sub[f"{tag}|{N}"] = {"delta_ps": r["delta_ps"], "chi2": r["chi2"]}
        for tag, (lo, hi) in (("10k_50M", (1e4, 5e7)), ("1M_200M", (1e6, 2e8)), ("10M_200M", (1e7, 2e8)), ("10k_10M", (1e4, 1e7))):
            sel = np.where((F >= lo * 0.999) & (F <= hi * 1.001))[0]
            for N in ("N2lam", "N3"):
                r = fd.fit_delta(y[sel], dd, b, 0.01, N, freq=F[sel])
                sub[f"band_{tag}|{N}"] = {"delta_ps": r["delta_ps"], "chi2": r["chi2"], "n": int(sel.size)}
        ent["subsets"] = sub
        log(L, "subsets", {k: round(v["delta_ps"], 2) for k, v in sub.items()})
        # sensitivity of delta_fit to a smooth relative perturbation of the two-port (size eps)
        pert = {}
        lf = np.log(F / F[0]) / np.log(F[-1] / F[0])
        shapes = {"cos": np.cos(np.pi * lf), "ramp": lf, "hf": lf ** 4}
        for sh_name, g in shapes.items():
            for eps in (1e-6, 1e-5):
                Yp = Ye[i] * (1 + eps * (1 + 0.5j) * g[:, None, None])
                rp = fd.fit_delta(data(Yp), dd, b, 0.01, "N2lam")
                pert[f"{sh_name}|{eps:.0e}"] = rp["delta_ps"] - ent["N2lam|0.01"]["delta_ps"]
        ent["perturbation_sensitivity_N2lam"] = pert
        log(L, "perturbation sensitivity (ps)", pert)
        # tau_R sensitivity (base tau_R x 0.5, x 2)
        sens = {}
        for fac in (0.5, 2.0):
            b2 = b.copy()
            b2[0] *= fac
            for N in ("N2lam", "N3"):
                sens[f"{fac}|{N}"] = fd.fit_delta(y, dd, b2, 0.01, N)["delta_ps"]
        ent["tauR_sensitivity"] = sens
        log(L, "tau_R sensitivity", sens)
        out[L] = ent
        res["ref"] = out
        save(res)


def stage_interp(res):
    d, Ye, The = load_exact()
    fsub, dsub, Ysub, Thsub = load_sub("REF")
    out = {}
    for i, L in enumerate(("300", "500")):
        dd = float(L) * 1e-9
        b = base_par(L)
        Ya = aaa_interp(fsub, Ysub[i])
        ent = {"aaa_max_rel_err_Y": float(np.max(np.abs(Ya - Ye[i])) / np.max(np.abs(Ye[i]))),
               "aaa_max_abs_dlnH": float(np.max(np.abs(data(Ya) - data(Ye[i]))))}
        # stored sub4 REF equals the exact all-run at the shared frequencies?
        idx = [int(np.argmin(np.abs(F - x))) for x in fsub if x >= 9e3]
        m = fsub >= 9e3
        ent["sub4_vs_all_same_freqs"] = float(np.max(np.abs(Ysub[i][m] - Ye[i][idx])) / np.max(np.abs(Ye[i])))
        for N in ("N2lam", "N3"):
            r1 = fd.fit_delta(data(Ye[i]), dd, b, 0.01, N)["delta_ps"]
            r2 = fd.fit_delta(data(Ya), dd, b, 0.01, N)["delta_ps"]
            ent[f"delta_exact_vs_aaa|{N}"] = [r1, r2, r2 - r1]
        # ratio-spline validated on the wall variants (exact on 80 points): variant / black
        for w in (1, 3, 4):
            _, Yw, _ = load_exact(wall=w)
            Ysw = Yw[i][idx]
            Ybs = Ye[i][idx]
            fin = fsub[m]
            # build ratio spline on the in-band sub4 points only
            x = np.log(fin)
            Yr = np.zeros_like(Ye[i])
            for a in range(2):
                for bb in range(2):
                    rr = Ysw[:, a, bb] / Ybs[:, a, bb]
                    Yr[:, a, bb] = Ye[i][:, a, bb] * (CubicSpline(x, rr.real)(np.log(F)) + 1j * CubicSpline(x, rr.imag)(np.log(F)))
            ent[f"ratio_spline_wall{w}_max_rel_err_Y"] = float(np.max(np.abs(Yr - Yw[i])) / np.max(np.abs(Yw[i])))
            ent[f"ratio_spline_wall{w}_max_abs_dlnH"] = float(np.max(np.abs(data(Yr) - data(Yw[i]))))
            for N in ("N2lam", "N3"):
                r1 = fd.fit_delta(data(Yw[i]), dd, b, 0.01, N)["delta_ps"]
                r2 = fd.fit_delta(data(Yr), dd, b, 0.01, N)["delta_ps"]
                ent[f"ratio_spline_wall{w}_delta|{N}"] = [r1, r2, r2 - r1]
        out[L] = ent
        log(L, "interp", ent)
    res["interp"] = out
    save(res)


def stage_wit(res):
    d, Ye, The = load_exact()
    fsub, dsub, Yref_sub, _ = load_sub("REF")
    out = res.get("wit", {})
    for i, L in enumerate(("300", "500")):
        dd = float(L) * 1e-9
        b = base_par(L)
        ent = out.get(L, {"delta": {}})
        for nm in ["REF"] + WIT:
            if nm in ent["delta"]:
                continue
            fs_, ds_, Ys_, _ = load_sub(nm)
            Ya = aaa_interp(fs_, Ys_[i])
            Yr = Ye[i] if nm == "REF" else ratio_spline(fs_, Ys_[i], Yref_sub[i], Ye[i])
            row = {}
            for N in ("N2lam", "N3"):
                for sp in (0.01, 0.1):
                    ra = fd.fit_delta(data(Ya), dd, b, sp, N)
                    row[f"{N}|{sp}|aaa"] = {"delta_ps": ra["delta_ps"], "chi2": ra["chi2"], "par": ra["par"],
                                            "spread_ps": ra["spread_ps"], "n_at_min": ra["n_at_min"]}
                    if sp == 0.01:
                        rr = fd.fit_delta(data(Yr), dd, b, sp, N)
                        row[f"{N}|{sp}|ratio"] = {"delta_ps": rr["delta_ps"], "chi2": rr["chi2"]}
            ent["delta"][nm] = row
            log(L, nm, {k: round(v["delta_ps"], 3) for k, v in row.items()})
            out[L] = ent
            res["wit"] = out
            save(res)
        # ranges, Fisher, margins
        summ = {}
        for N in ("N2lam", "N3"):
            for sp in (0.01, 0.1):
                for meth in ("aaa", "ratio") if sp == 0.01 else ("aaa",):
                    key = f"{N}|{sp}|{meth}"
                    ref = ent["delta"]["REF"][key]["delta_ps"]
                    ex = {w: ent["delta"][w][key]["delta_ps"] - ref for w in WIT}
                    up, low = max(list(ex.values()) + [0.0]), min(list(ex.values()) + [0.0])
                    summ[key] = {"ref_ps": ref, "up": up, "low": low, "half_range_ps": 0.5 * (up - low),
                                 "argmax": max(ex, key=ex.get), "argmin": min(ex, key=ex.get), "excursions": ex}
        # Fisher sigma at the REF slab-fit base point of each set (0.01-deg fit), all four sets read at that point
        sig = {}
        for Nb in ("N2lam", "N3"):
            par = np.array(ent["delta"]["REF"][f"{Nb}|0.01|aaa"]["par"])
            J16 = fd.cauchy_jac(par, dd, M=16, rel=0.05)
            J24 = fd.cauchy_jac(par, dd, M=24, rel=0.025)
            sig[f"jac_check_base{Nb}"] = float(np.max(np.abs(J16 - J24)) / np.max(np.abs(J24)))
            for N in ("N1", "N2", "N2lam", "N3"):
                for sp in (0.01, 0.03, 0.1):
                    s_tauR = fd.sigma_direction(J16, 0, fd.SETCOLS[N], sp) * 1e12
                    s_taul = fd.sigma_direction(J16, 1, fd.SETCOLS[N], sp) * 1e12
                    sig[f"base{Nb}|{N}|{sp}"] = {"sigma_via_tauR_ps": s_tauR, "sigma_via_taul_ps": s_taul}
        marg = {}
        for N in ("N1", "N2", "N2lam", "N3"):
            Nb = N if N in ("N2lam", "N3") else "N2lam"
            hr = summ[f"{Nb}|0.01|aaa"]["half_range_ps"]
            for sp in (0.01, 0.03, 0.1):
                marg[f"{N}|{sp}"] = sig[f"base{Nb}|{N}|{sp}"]["sigma_via_tauR_ps"] / hr
        # direct test: witness data against the GK film at the REF fit (tau_l fixed), the set's nuisances refitted
        direct = {}
        for N in ("N2lam", "N3", "N2", "N1"):
            Nb = N if N in ("N2lam", "N3") else "N2lam"
            par = np.array(ent["delta"]["REF"][f"{Nb}|0.01|aaa"]["par"])
            for sp in (0.01, 0.1):
                chis = {}
                for nm in WIT:
                    fs_, ds_, Ys_, _ = load_sub(nm)
                    y = data(aaa_interp(fs_, Ys_[i]))
                    Fi = fd.Fitter(y, dd, par, sp)
                    nl = fd.NL_SETS[N]
                    chi, _, _ = Fi.fit(nl, par, fd.nuisance_starts(par, nl))
                    chis[nm] = chi
                wm = max(chis, key=chis.get)
                direct[f"{N}|{sp}"] = {"max_dchi2": chis[wm], "at": wm, "m_direct": float(1 / np.sqrt(chis[wm])),
                                       "per_witness": chis}
        ent.update({"ranges": summ, "sigma": sig, "margins": marg, "direct": direct})
        out[L] = ent
        res["wit"] = out
        save(res)
        log(L, "half-ranges", {k: round(v["half_range_ps"], 3) for k, v in summ.items()})
        log(L, "sigma", {k: (round(v["sigma_via_tauR_ps"], 2) if isinstance(v, dict) else v) for k, v in sig.items()})
        log(L, "margins", {k: round(v, 3) for k, v in marg.items()})
        log(L, "direct", {k: (round(v["m_direct"], 3), v["at"]) for k, v in direct.items()})


WITW = ["cal1_fL100_max", "cal1_fL100_min"] + [f"orb0.1_fL{L}_{s}" for L in (30, 100, 300, 1000) for s in ("max", "min")]


def stage_mesh(res):
    """Same extraction and margins on 5x5x3 (REF + ten film witnesses, admissible on 5x5x3) and on 7x7x5 for the two
    film witnesses run there (same-tensor probe of the mesh trend)."""
    out = {}
    for mesh, names in (("m553", WITW), ("m775", ["cal1_fL100_min", "orb0.1_fL1000_max"])):
        st = json.loads((RS / f"rs2_{mesh}_filmw_static_REF.json").read_text())
        ZR = np.load(RS / f"rs2_{mesh}_filmw_full_sub4_REF.npz")
        fR, dR, YR = ZR["f"], ZR["d"], ZR["Y"]
        YR = YR[0] if YR.ndim == 5 else YR
        ent = {}
        for i, L in enumerate(("300", "500")):
            dd = float(L) * 1e-9
            kr = st[f"reg_full_QDQ_{L}"]
            b = np.array([kr["tauR"], kr["taul"], np.log(kr["K"]), np.log(st["C_tot"]), np.log(R1N), np.log(R2N)])
            dl = {}
            for nm in ["REF"] + names:
                Z = np.load(RS / f"rs2_{mesh}_filmw_full_sub4_{nm}.npz")
                Y = Z["Y"][0] if Z["Y"].ndim == 5 else Z["Y"]
                y = data(aaa_interp(Z["f"], Y[i]))
                dl[nm] = {N: fd.fit_delta(y, dd, b, 0.01, N) for N in ("N2lam", "N3")}
            row = {"delta": {nm: {N: v[N]["delta_ps"] for N in v} for nm, v in dl.items()},
                   "chi2_max": max(v[N]["chi2"] for v in dl.values() for N in v)}
            for N in ("N2lam", "N3"):
                ex = [dl[w][N]["delta_ps"] - dl["REF"][N]["delta_ps"] for w in names]
                up, low = max(ex + [0.0]), min(ex + [0.0])
                hr = 0.5 * (up - low)
                J = fd.cauchy_jac(np.array(dl["REF"][N]["par"]), dd)
                sig = fd.sigma_direction(J, 0, fd.SETCOLS[N], 0.01) * 1e12
                row[N] = {"half_range_ps": hr, "up": up, "low": low, "sigma_ps": sig, "margin": sig / hr,
                          "ref_delta_ps": dl["REF"][N]["delta_ps"], "lamz_fit": float(np.exp(dl["REF"][N]["par"][2]))}
            ent[L] = row
            log(mesh, L, {N: {k: round(v, 3) for k, v in row[N].items()} for N in ("N2lam", "N3")})
        out[mesh] = ent
        res["mesh_trend"] = out
        save(res)


def main():
    stages = sys.argv[1:] or ["ref", "interp", "wit"]
    res = json.loads(RES.read_text()) if RES.exists() else {}
    for st in stages:
        {"ref": stage_ref, "interp": stage_interp, "wit": stage_wit, "mesh": stage_mesh}[st](res)
    res["seconds"] = time.perf_counter() - T0
    save(res)


if __name__ == "__main__":
    main()
