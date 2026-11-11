"""A1: base values of the AlN layer, tau_z(L), kappa_z(L), kappa_x(L), at L = 30, 100, 300, 500 nm (300 K).

Two routes, same heuristic Matthiessen regulator (physical boundary rate |v_z|/L, diffuse grey boundary of a
c-oriented film; fc3-memory E26/E27/E30):
  (a) RTA on the converged phono3py meshes 48x48x26, 56x56x30, 64x64x35 (chunks of fc3-memory E27):
      1/tau_L = 1/tau_3ph + |v_z|/L, tau_mem(L) = sum w tau_L^2 / sum w tau_L, w = mode_kappa / tau_3ph,
      kappa_L = kappa * sum w tau_L / sum w tau_3ph (same formula as rta_chunks.py 'assemble').
      |v_z| is invariant on the star of each irreducible point (point group 6mm and time reversal), so the
      irreducible-point sum is exact for the z regulator.
  (b) conserving event operator, 12x12x8, REF tensor (fc3-memory E30/E37 results_c/film3/*_m12128.json):
      the operator diagonal there is d = |v_z|/(4 pi L) in the THz linewidth units of the operator, i.e. the same
      physical rate |v_z|/L. L = 500 nm is not available: interpolated (log-log) between 300 and 1000 nm; the
      interpolation error is estimated by applying the same interpolation to the RTA curve, where 500 nm is known.
Usage: python a1_base_values.py   (reads D: campaign data, writes out/a1_base_values.json)
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
RTA = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory/results_c/rta")
FILM3 = Path("C:/Users/Koussay/these/theory/aln/fc3_memory/results_c/film3")
LS = [30, 100, 300, 500, 1000, 3000, 10000]


def rta_film(mesh_dir: Path):
    meta = json.loads((mesh_dir / "meta.json").read_text())
    fs = sorted(f for f in mesh_dir.glob("chunk_*.npz") if ".tmp" not in f.name)
    assert len(fs) == meta["n_chunks"], (mesh_dir, len(fs), meta["n_chunks"])
    D = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk")}
    assert int(D["w"].sum()) == meta["sum_weights"]
    ok = (D["gamma"] > 0) & (D["f"] > 1e-3)
    tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, D["gamma"], 1)), 0.0)    # ps
    vz = np.abs(D["gv"][..., 2]) * 100.0                                      # THz*A -> m/s
    res = {"mesh": meta["mesh"], "f_min_THz": float(D["f"][D["f"] > 1e-3].min())}
    for c, nm in ((0, "x"), (2, "z")):
        k = np.where(ok, D["mk"][..., c], 0.0)
        K = k.sum()
        kappa = float(K / meta["sum_weights"])
        wgt = np.where(ok, k / np.where(ok, tau, 1), 0.0)
        res[f"kappa_{nm}_bulk"] = kappa
        res[f"tau_{nm}_bulk"] = float((k * tau).sum() / K)
        res[f"tau_{nm}"], res[f"kappa_{nm}"] = {}, {}
        for L in LS:
            tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vz / (L * 1e-9) * 1e-12), 0.0)
            res[f"tau_{nm}"][str(L)] = float((wgt * tl ** 2).sum() / (wgt * tl).sum())
            res[f"kappa_{nm}"][str(L)] = float(kappa * (wgt * tl).sum() / (wgt * tau).sum())
    # share of the q_z = 0 plane (v_z = 0 exactly by C2z + time reversal; never suppressed by |v_z|/L) in kappa_x(L)
    D2 = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("q",)}
    qz = np.abs(D2["q"][:, 2])
    qz = np.minimum(qz, 1 - qz)
    plane = qz < 1e-9
    k0 = np.where(ok, D["mk"][..., 0], 0.0)
    w0 = np.where(ok, k0 / np.where(ok, tau, 1), 0.0)
    res["qz0_plane_weight"] = float(D["w"][plane].sum() / D["w"].sum())
    res["qz0_plane_max_abs_vz"] = float(vz[plane].max())
    res["kappa_x_qz0_plane"] = {}
    for L in LS:
        tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vz / (L * 1e-9) * 1e-12), 0.0)
        res["kappa_x_qz0_plane"][str(L)] = float(res["kappa_x_bulk"] * (w0 * tl)[plane].sum() / (w0 * tau).sum())
    # regression against the archived assembly (E27)
    ref = json.loads((mesh_dir / "assembled.json").read_text())
    dev = max(abs(res[f"tau_{nm}"][L] / ref[f"film_vz_{nm}"][L] - 1)
              for nm in "xz" for L in ref["film_vz_z"])
    res["max_rel_dev_vs_assembled"] = float(dev)
    res["kappa_rel_dev_vs_assembled"] = float(max(abs(res[f"kappa_{nm}_bulk"] / ref[f"kappa_{nm}"] - 1) for nm in "xz"))
    return res


def loglog(x0, y0, x1, y1, x):
    t = (np.log(x) - np.log(x0)) / (np.log(x1) - np.log(x0))
    return float(np.exp(np.log(y0) + t * (np.log(y1) - np.log(y0))))


def main():
    out = {"regulator": "physical boundary rate |v_z|/L (Matthiessen, heuristic); event operator diag |v_z|/(4 pi L) "
                        "in THz linewidth units = same physical rate", "L_nm": LS, "rta": {}}
    for m in ("m484826", "m565630", "m646435"):
        out["rta"][m] = rta_film(RTA / m)
        r = out["rta"][m]
        print(m, "dev vs assembled", r["max_rel_dev_vs_assembled"],
              {L: (round(r["tau_z"][str(L)], 3), round(r["kappa_z"][str(L)], 2)) for L in LS})
    conv = out["rta"]["m646435"]
    # mesh convergence of the RTA film values (two finest steps)
    out["rta_mesh_change"] = {
        q: {str(L): [out["rta"]["m565630"][q][str(L)] / out["rta"]["m484826"][q][str(L)] - 1,
                     out["rta"]["m646435"][q][str(L)] / out["rta"]["m565630"][q][str(L)] - 1] for L in LS}
        for q in ("tau_z", "kappa_z", "kappa_x")}

    # in-plane film value: linear extrapolation in 1/N_z over the three meshes (the q_z = 0 plane error is O(1/N_z))
    nz = np.array([out["rta"][m]["mesh"][2] for m in ("m484826", "m565630", "m646435")], float)
    out["kappa_x_extrap_1_over_Nz"] = {}
    for L in LS:
        y = np.array([out["rta"][m]["kappa_x"][str(L)] for m in ("m484826", "m565630", "m646435")])
        A = np.vstack([np.ones(3), 1 / nz]).T
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        out["kappa_x_extrap_1_over_Nz"][str(L)] = {"intercept": float(c[0]), "slope": float(c[1]),
                                                   "max_fit_residual": float(np.max(np.abs(A @ c - y)))}
    # event operator, REF, 12x12x8: identical REF blocks in every file (check)
    refs = []
    for fn in sorted(glob.glob(str(FILM3 / "*_m12128.json"))):
        refs.append(json.load(open(fn))["REF"])
    keys = ("tau_z", "kappa_z", "kappa_x", "tau_x", "tauRTA_z")
    spread = max(abs(r[L][q] / refs[0][L][q] - 1) for r in refs for L in refs[0] for q in keys)
    ev = refs[0]
    out["event_m12128_REF"] = {L: {q: ev[L][q] for q in keys} for L in ev}
    out["event_m12128_REF_spread_across_files"] = float(spread)

    # 500 nm by log-log interpolation 300 -> 1000 nm, with error estimate from the RTA curve
    est = {}
    for q in ("tau_z", "kappa_z", "kappa_x"):
        y_ev = loglog(300, ev["300"][q], 1000, ev["1000"][q], 500)
        y_rta_int = loglog(300, conv[q]["300"], 1000, conv[q]["1000"], 500)
        y_rta = conv[q]["500"]
        err_rta = y_rta_int / y_rta - 1
        # ratio method: event/RTA(converged) interpolated, times RTA(converged) at 500 nm
        rat = loglog(300, ev["300"][q] / conv[q]["300"], 1000, ev["1000"][q] / conv[q]["1000"], 500)
        y_ratio = rat * y_rta
        # PCHIP in (ln L, ln y) through all six event-operator thicknesses
        Ls = np.array([30, 100, 300, 1000, 3000, 10000], float)
        ys = np.array([ev[str(int(L))][q] for L in Ls])
        y_pchip = float(np.exp(PchipInterpolator(np.log(Ls), np.log(ys))(np.log(500.0))))
        y_pchip_rta = float(np.exp(PchipInterpolator(np.log(Ls), np.log([conv[q][str(int(L))] for L in Ls]))(np.log(500.))))
        est[q] = {"event_loglog_300_1000": y_ev, "event_ratio_to_RTA": y_ratio, "event_pchip_all_L": y_pchip,
                  "rta_converged_direct": y_rta, "rta_loglog_300_1000": y_rta_int,
                  "rta_pchip_all_L": y_pchip_rta,
                  "interp_error_loglog_from_rta": err_rta,
                  "interp_error_pchip_from_rta": y_pchip_rta / y_rta - 1,
                  "spread_event_estimates": max(y_ev, y_ratio, y_pchip) / min(y_ev, y_ratio, y_pchip) - 1}
    out["L500"] = est

    # base values adopted for the twin (event operator, campaign operator; 500 nm: ratio method, which uses the
    # direct RTA value at 500 nm and interpolates only the smooth event/RTA ratio)
    base = {}
    for L in ("30", "100", "300"):
        base[L] = {"tau_R_ps": ev[L]["tau_z"], "lam_z": ev[L]["kappa_z"], "lam_r": ev[L]["kappa_x"],
                   "source": "event operator REF 12x12x8 (E30 addendum / E37 files)"}
    base["500"] = {"tau_R_ps": est["tau_z"]["event_ratio_to_RTA"], "lam_z": est["kappa_z"]["event_ratio_to_RTA"],
                   "lam_r": est["kappa_x"]["event_ratio_to_RTA"],
                   "source": "event operator REF 12x12x8, event/RTA ratio interpolated log-log 300-1000 nm, times RTA "
                             "64x64x35 at 500 nm"}
    base["500_rta"] = {"tau_R_ps": conv["tau_z"]["500"], "lam_z": conv["kappa_z"]["500"],
                       "lam_r": conv["kappa_x"]["500"], "source": "RTA 64x64x35 direct (sensitivity case)"}
    for L in ("30", "100", "300"):
        base[L + "_rta"] = {"tau_R_ps": conv["tau_z"][L], "lam_z": conv["kappa_z"][L], "lam_r": conv["kappa_x"][L],
                            "source": "RTA 64x64x35 direct (sensitivity case)"}
    for L in ("30", "100", "300", "500"):
        base[L + "_xextrap"] = {"lam_r": out["kappa_x_extrap_1_over_Nz"][L]["intercept"],
                                "source": "RTA kappa_x(L), linear extrapolation in 1/N_z (sensitivity case)"}
    out["base"] = base
    (OUT / "a1_base_values.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(est, indent=1))
    print(json.dumps(base, indent=1))
    print("REF spread across files", spread)
    print(json.dumps(out["rta_mesh_change"], indent=0)[:2000])


if __name__ == "__main__":
    main()
