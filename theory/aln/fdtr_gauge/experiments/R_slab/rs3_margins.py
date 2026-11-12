"""RS3: static-data range of the FDTR-visible delta and its resolvability (7x7x5, REF + 12 witnesses).

Three readings on the same mesh and witnesses:
  REG  : first-order kinetic closure of the regularised full operator (C + QDQ, |v_z|/L): delta_reg = tau_R - tau_l;
  RTA  : kinetic slab, diffuse black walls, RTA diagonal of each tensor (unbinned): delta_fit;
  FULL : kinetic slab, diffuse black walls, full conserving operator of each tensor: delta_fit.
delta_fit = tau_R,base - tau_l,fit from the GK-film fit of the slab's FDTR data (rs2b; wall port, note-16 stack, 0.01 deg),
witnesses fitted with the REF base (as au8). Half-range by the campaign rule: ends = largest excursion of each sign
(0 if none), half-range = (up - low)/2.
sigma(delta) = marginal standard deviation along tau_R at fixed tau_l (= d/d delta) of the 1D first-order GK film
(au5_fisher1d.jac, sigma_dir) at the reading's REF base point: REG - (tau_R, tau_l, K)_reg with nominal R1, R2;
slab readings - (tau_R,base, tau_R,base - delta_fit, lam_z,fit, C, R1,fit, R2,fit) of the REF slab fit (N2lam), and as a
check the audit's convention (lam_z = K_reg, nominal R1, R2). Sets N1, N2 (lam_z and C known; in the slab readings at
the REF-fitted lam_z), N2lam, N3; sigma_phi 0.01, 0.03, 0.1 deg; note-16 and favourable stacks.
Margin m = sigma / half-range (m <= 1: resolvable at 1 sigma). Direct test: each witness slab's data against the GK film
at the REF-fitted point with tau_l fixed and the set's nuisances refitted; m_direct = 1/sqrt(max Delta chi^2).
Required prior on lam_z (au9 method) for the FULL reading. 9x9x6 projection: campaign half-range (au4) x the 7x7x5
range transfer FULL/REG, with and without the audit's geometric mesh continuation (labelled projections).
usage: python -B rs3_margins.py [freqset]      Output: out/rs3_margins.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
AUD = HERE.parent / "R_A_audit"
sys.path.insert(0, str(AUD))
sys.path.insert(0, str(HERE))
import au5_fisher1d as a5  # noqa: E402
import au6_slab as a6  # noqa: E402
import rs2b_fits as fb  # noqa: E402
import rs_core as rc  # noqa: E402

FT = sys.argv[1] if len(sys.argv) > 1 else "sub"
R2B = json.loads((OUT / f"rs2b_m775_film2_{FT}.json").read_text())
AU4 = json.loads((AUD / "out" / "au4_delta_range.json").read_text())
ALLSETS = {**a5.SETS, **a5.SETS_EXTRA}
SIGS = (0.01, 0.03, 0.1)
WIT = rc.WIT12


def halfrange(dv):
    up, low = max(list(dv) + [0.0]), min(list(dv) + [0.0])
    return {"up": up, "low": low, "half_range_ps": 0.5 * (up - low)}


def sigma_prior(J, nuis, sp, s_rel):
    Jw = a5.whiten(J, sp)
    v = Jw[:, 0]
    N = Jw[:, nuis]
    col = nuis.index(2)
    prow = np.zeros(len(nuis))
    prow[col] = 1.0 / s_rel
    A = np.vstack([np.column_stack([v, N]), np.concatenate([[0.0], prow])[None, :]])
    return float(np.sqrt(np.linalg.inv(A.T @ A)[0, 0]))


def main():
    res = {"freqset": FT}
    for L in ("300", "500"):
        d = float(L) * 1e-9
        rl = {}
        # ---- ranges
        reg = R2B["full"]
        dREG = {nm: reg[nm][f"delta_reg_{L}_ps"] for nm in ["REF"] + WIT}
        dREGp = {nm: reg[nm][f"delta_reg_phys_{L}_ps"] for nm in ["REF"] + WIT}
        rng_ = {"REG": halfrange([dREG[w] - dREG["REF"] for w in WIT]),
                "REG_phys": halfrange([dREGp[w] - dREGp["REF"] for w in WIT])}
        rng_["REG"]["ref_ps"] = dREG["REF"]
        rng_["REG_phys"]["ref_ps"] = dREGp["REF"]
        dRTAreg = {nm: R2B["rta"][nm][f"delta_reg_{L}_ps"] for nm in ["REF"] + WIT}
        rng_["RTA_reg"] = dict(halfrange([dRTAreg[w] - dRTAreg["REF"] for w in WIT]), ref_ps=dRTAreg["REF"])
        for op, tag in (("rta", "RTA"), ("full", "FULL")):
            for N in ("N2lam", "N3", "N2fix", "N1fix"):
                for sp in ("0.01", "0.1"):
                    key = f"{N}|{sp}"
                    vals = {nm: R2B[op][nm][f"fit_{L}|wall|note16"][key]["delta_fit_ps"] for nm in ["REF"] + WIT}
                    chis = {nm: R2B[op][nm][f"fit_{L}|wall|note16"][key]["dchi2"] for nm in ["REF"] + WIT}
                    rng_[f"{tag}|{key}"] = dict(halfrange([vals[w] - vals["REF"] for w in WIT]), ref_ps=vals["REF"],
                                                per_witness={w: vals[w] - vals["REF"] for w in WIT},
                                                max_dchi2=max(chis.values()))
            # energy-port convention check (0.01 deg)
            for N in ("N2lam", "N3"):
                vals = {nm: R2B[op][nm][f"fit_{L}|energy|note16"][f"{N}|0.01"]["delta_fit_ps"] for nm in ["REF"] + WIT}
                rng_[f"{tag}|energyport|{N}|0.01"] = dict(halfrange([vals[w] - vals["REF"] for w in WIT]), ref_ps=vals["REF"])
        # transfers (7x7x5): slab-visible half-range / regularised half-range of the same operator family
        rng_["transfer"] = {
            "FULL_over_REG|N2lam": rng_["FULL|N2lam|0.01"]["half_range_ps"] / rng_["REG"]["half_range_ps"],
            "FULL_over_REG|N3": rng_["FULL|N3|0.01"]["half_range_ps"] / rng_["REG"]["half_range_ps"],
            "RTA_over_RTAreg|N2lam": rng_["RTA|N2lam|0.01"]["half_range_ps"] / rng_["RTA_reg"]["half_range_ps"],
            "RTA_over_RTAreg|N3": rng_["RTA|N3|0.01"]["half_range_ps"] / rng_["RTA_reg"]["half_range_ps"],
            "FULL_over_RTA|N2lam": rng_["FULL|N2lam|0.01"]["half_range_ps"] / rng_["RTA|N2lam|0.01"]["half_range_ps"],
            "REF_FULL_over_REG|N3": rng_["FULL|N3|0.01"]["ref_ps"] / rng_["REG"]["ref_ps"],
            "REF_RTA_over_RTAreg|N3": rng_["RTA|N3|0.01"]["ref_ps"] / rng_["RTA_reg"]["ref_ps"]}
        rl["ranges"] = rng_
        # ---- sigma(delta) at the base points
        sig = {}
        for stack, (R1, R2) in fb.STACKS.items():
            bases = {}
            kr = R2B["full"]["REF"][f"base_{L}"]
            bases["REG"] = {"tauR": kr["tauR"] * 1e-12, "taul": kr["taul"] * 1e-12, "lamz": kr["lamz"], "C": kr["C"],
                            "R1": R1, "R2": R2, "d": d}
            for op, tag in (("rta", "RTA"), ("full", "FULL")):
                b = R2B[op]["REF"][f"base_{L}"]
                ft = R2B[op]["REF"][f"fit_{L}|wall|{stack}"]
                for N in ("N2lam", "N3"):
                    f_ = ft[f"{N}|0.01"]
                    bases[f"{tag}@{N}"] = {"tauR": b["tauR"] * 1e-12, "taul": (b["tauR"] - f_["delta_fit_ps"]) * 1e-12,
                                           "lamz": f_["lamz_fit"], "C": f_["C_fit"], "R1": f_["R1_fit"], "R2": f_["R2_fit"],
                                           "d": d}
                f_ = ft["N2lam|0.01"]
                bases[f"{tag}@audit"] = {"tauR": b["tauR"] * 1e-12, "taul": (b["tauR"] - f_["delta_fit_ps"]) * 1e-12,
                                         "lamz": kr["lamz"], "C": kr["C"], "R1": R1, "R2": R2, "d": d}
            for bt, base in bases.items():
                J = a5.jac(base)
                for N, nuis in ALLSETS.items():
                    for sp in SIGS:
                        sig[f"{stack}|{bt}|{N}|{sp}"] = a5.sigma_dir(J, (1, 0), nuis, sp)
        rl["sigma_ps"] = sig
        # ---- margins (note-16 and favourable); range of the same reading; N1/N2 use the N2lam range in slab readings
        marg = {}
        for stack in fb.STACKS:
            for N in ("N1", "N2", "N2lam", "N3"):
                for sp in SIGS:
                    s_reg = sig[f"{stack}|REG|{N}|{sp}"]
                    m = {"REG": s_reg / rng_["REG"]["half_range_ps"]}
                    for tag in ("RTA", "FULL"):
                        Nr = N if N in ("N2lam", "N3") else "N2lam"
                        bt = f"{tag}@{Nr}"
                        hr = rng_[f"{tag}|{Nr}|0.01"]["half_range_ps"]
                        m[tag] = sig[f"{stack}|{bt}|{N}|{sp}"] / hr
                        m[f"{tag}_auditbase"] = sig[f"{stack}|{tag}@audit|{N}|{sp}"] / hr
                    marg[f"{stack}|{N}|{sp}"] = m
        rl["margins"] = marg
        # ---- 9x9x6 projection of the FULL reading (labelled projection)
        t4 = AU4["table"][f"{float(L):.6g}"]["all"]
        hr_c = t4["base_point"]["halfrange_delta_ps_campaign_rule"]
        fmesh = t4["x12_geometric"]["halfrange_delta_over_tau0"] / t4["robust"]["halfrange_delta_over_tau0"]
        proj = {}
        for N in ("N2lam", "N3"):
            T = rng_["transfer"][f"FULL_over_REG|{N}"]
            for sp in SIGS:
                s_ = sig[f"note16|FULL@{N}|{N}|{sp}"]
                proj[f"{N}|{sp}"] = {"halfrange_996_projected_ps": hr_c * T, "m": s_ / (hr_c * T),
                                     "m_mesh_continued": s_ / (hr_c * T * fmesh)}
        rl["projection_996"] = {"campaign_halfrange_ps": hr_c, "f_mesh": fmesh, "rows": proj}
        # ---- direct tests (wall port, note-16), FULL and RTA readings
        direct = {}
        p = 2j * np.pi * a5.FREQ
        for op, tag in (("full", "FULL"), ("rta", "RTA")):
            refb = R2B[op]["REF"][f"base_{L}"]
            ft = R2B[op]["REF"][f"fit_{L}|wall|note16"]
            for N in ("N2lam", "N3", "N2", "N1"):
                Nr = N if N in ("N2lam", "N3") else "N2lam"
                f_ = ft[f"{Nr}|0.01"]
                base = {"tauR": refb["tauR"] * 1e-12, "taul": (refb["tauR"] - f_["delta_fit_ps"]) * 1e-12,
                        "lamz": f_["lamz_fit"], "C": f_["C_fit"], "R1": f_["R1_fit"], "R2": f_["R2_fit"], "d": d}
                for sp in (0.01, 0.1):
                    chis = {}
                    for w in WIT:
                        fz, dz, Y, Th, _ = fb.load("m775", "film2", op, FT, w)
                        k = list(np.round(dz * 1e9).astype(int)).index(int(L))
                        Yf, Tf, _ = fb.at_fdtr(fz, Y, Th, FT)
                        y = np.log(a6.Z_stack_slab(list(Yf[k]), list(Tf[k]), p, 1e-8, 2e-8, port="wall"))
                        chi, _ = a5.gn_refit(y, base, ALLSETS[N], sp)
                        chis[w] = chi
                    wm = max(chis, key=chis.get)
                    direct[f"{tag}|{N}|{sp}"] = {"max_dchi2": chis[wm], "at": wm, "m_direct": 1 / np.sqrt(chis[wm]),
                                                 "per_witness": chis}
        rl["direct"] = direct
        # ---- required prior on lam_z (FULL reading, N2 + ln lam_z with prior, note-16)
        pri = {}
        Sg = np.geomspace(1e-4, 0.3, 61)
        for tag in ("FULL", "RTA"):
            op = tag.lower()
            refb = R2B[op]["REF"][f"base_{L}"]
            f_ = R2B[op]["REF"][f"fit_{L}|wall|note16"]["N2lam|0.01"]
            base = {"tauR": refb["tauR"] * 1e-12, "taul": (refb["tauR"] - f_["delta_fit_ps"]) * 1e-12,
                    "lamz": f_["lamz_fit"], "C": f_["C_fit"], "R1": f_["R1_fit"], "R2": f_["R2_fit"], "d": d}
            J = a5.jac(base)
            hr = rng_[f"{tag}|N2lam|0.01"]["half_range_ps"]
            for sp in (0.01, 0.1):
                s_ = np.array([sigma_prior(J, a5.SETS["N2"] + [2], sp, s) for s in Sg])
                m = s_ / hr
                ok = Sg[m <= 1.0]
                pri[f"{tag}|{sp}"] = {"s_star_rel": float(ok.max()) if ok.size else None,
                                      "m_lam_exact": float(sigma_prior(J, a5.SETS["N2"] + [2], sp, 1e-9) / hr),
                                      "m_at_1pct": float(np.interp(np.log(0.01), np.log(Sg), m)),
                                      "m_at_5pct": float(np.interp(np.log(0.05), np.log(Sg), m))}
        rl["lam_prior"] = pri
        # ---- wall-model sensitivity of the FULL reading (range, REF delta, margins at the wall model's REF fit)
        wm = {}
        walls = fb.walls_of("m775", "film2", "full", FT, "REF")
        for w in walls:
            wt = fb.wall_tag(w)
            sfx = "" if wt == "black" else f"|{wt}"
            ent = {}
            for N in ("N2lam", "N3"):
                key = f"{N}|0.01"
                vals = {nm: R2B["full"][nm][f"fit_{L}|wall|note16" + sfx][key]["delta_fit_ps"] for nm in ["REF"] + WIT}
                hr = halfrange([vals[x] - vals["REF"] for x in WIT])
                f_ = R2B["full"]["REF"][f"fit_{L}|wall|note16" + sfx][key]
                refb = R2B["full"]["REF"][f"base_{L}"]
                base = {"tauR": refb["tauR"] * 1e-12, "taul": (refb["tauR"] - f_["delta_fit_ps"]) * 1e-12,
                        "lamz": f_["lamz_fit"], "C": f_["C_fit"], "R1": f_["R1_fit"], "R2": f_["R2_fit"], "d": d}
                J = a5.jac(base)
                ent[N] = {"ref_ps": vals["REF"], **hr, "transfer_over_REG": hr["half_range_ps"] / rng_["REG"]["half_range_ps"],
                          "lamz_fit": f_["lamz_fit"], "R1_fit": f_["R1_fit"], "R2_fit": f_["R2_fit"], "dchi2_ref": f_["dchi2"],
                          "margins": {f"{NN}|{sp}": a5.sigma_dir(J, (1, 0), ALLSETS[NN], sp) / hr["half_range_ps"]
                                      for NN in ("N1", "N2", "N2lam", "N3") for sp in SIGS}}
            wm[wt] = ent
        rl["walls"] = wm
        print(L, "walls:", {k: {N: (round(v[N]["ref_ps"], 2), round(v[N]["half_range_ps"], 2),
                                    round(v[N]["margins"][f"{N}|0.01"], 2), round(v[N]["margins"][f"{N}|0.1"], 2))
                                for N in v} for k, v in wm.items()}, flush=True)
        res[L] = rl
        print(L, "ranges (ref, half):", {k: (round(v["ref_ps"], 2), round(v["half_range_ps"], 2)) for k, v in rng_.items()
                                         if isinstance(v, dict) and "ref_ps" in v and k.endswith(("0.01", "REG", "RTA_reg"))},
              flush=True)
        print(L, "transfer", {k: round(v, 3) for k, v in rng_["transfer"].items()}, flush=True)
        print(L, "margins note16:", {k.split("|", 1)[1]: {kk: round(vv, 2) for kk, vv in v.items()}
                                     for k, v in marg.items() if k.startswith("note16") and not k.endswith("0.03")}, flush=True)
        print(L, "direct:", {k: (round(v["m_direct"], 2), v["at"]) for k, v in direct.items()}, flush=True)
        print(L, "lam prior:", pri, flush=True)
        (OUT / "rs3_margins.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
