"""AU5b: FDTR consequence of the k^2 truncation, within the regularised bulk RTA model (64x64x35, V1 completion).

For L = 300, 500 nm: lam_eff(omega)/lam_z from the exact diffusive root of p|e|^2 + k^2 kappa(k, p) = 0 (all orders
in k, au2_rta) at the 80 FDTR frequencies, against the first-order GK film (1 + tau_l p)/(1 + tau_R p) with the same
static times. 1D twin stack of au5 (A1 lam_z, C_AlN), note-16 and favourable interface resistances.
  * Delta chi^2 (noiseless, nuisances N1/N2/N3 refitted by Levenberg-Marquardt) between exact-kernel data and the GK film
    with the static times: the size of the truncation error in noise units;
  * the GK film fitted to the exact-kernel data with tau_l free (tau_R fixed) and with both free: the delta that a
    first-order analysis would report, against the static delta;
  * transfer test of a witness-like change: tau_3ph of all modes scaled by (1 + eta) (eta = +-0.1, which moves the
    static tau_l by about the witness half-range), exact-kernel data of the perturbed operator fitted by the GK film
    with tau_l free: ratio of the fitted delta change to the static delta change.
Output: out/au5b_ktrunc.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au2_rta as a2  # noqa: E402
import au5_fisher1d as a5  # noqa: E402

SETS = a5.SETS
res = {}
md = a2.load_modes("m646435")
c, v, t3 = a2.directed(md)
for L in (300, 500):
    b = a5.A1[str(L)]
    ent = {}
    base_ops = {}
    for eta in (0.0, -0.1, 0.1):
        op, e, ehat = a2.build_variant("V1", c, v, t3 * (1 + eta), L)
        cf = a2.closed_form(op, e, ehat, v)
        droot = a2.diffusive_root(op, e, v, cf["K"], a5.FREQ)
        r_ex = np.array([r["lam_eff_rel"] for r in droot])
        base_ops[eta] = {"tauR": cf["tauR"], "taul": cf["taul"], "ratio": r_ex,
                         "max_newton": int(max(r["newton_iters"] for r in droot))}
    s0 = base_ops[0.0]
    ent["static"] = {str(eta): {"tauR_ps": d["tauR"] * 1e12, "taul_ps": d["taul"] * 1e12,
                                "delta_ps": (d["tauR"] - d["taul"]) * 1e12, "max_newton_iters": d["max_newton"]}
                     for eta, d in base_ops.items()}
    for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
        base = {"tauR": s0["tauR"], "taul": s0["taul"], "lamz": b["lam_z"], "C": a5.C_ALN, "R1": R1, "R2": R2,
                "d": L * 1e-9}
        y_ex = a5.lnH(np.zeros(9), base, s0["ratio"])
        y_gk = a5.lnH(np.zeros(9), base)
        row = {"max_abs_dphase_deg": float(np.degrees(np.max(np.abs((y_ex - y_gk).imag)))),
               "max_abs_dlnamp": float(np.max(np.abs((y_ex - y_gk).real)))}
        for N, nuis in SETS.items():
            for sp in (0.01, 0.1):
                chi0, _ = a5.gn_refit(y_ex, base, nuis, sp)
                chi1, th1 = a5.gn_refit(y_ex, base, nuis, sp, free=(1,))
                chi2, th2 = a5.gn_refit(y_ex, base, nuis, sp, free=(0, 1))
                d_fit1 = (s0["tauR"] - s0["taul"] - th1[1] * 1e-12) * 1e12
                d_fit2 = (s0["tauR"] + th2[0] * 1e-12 - s0["taul"] - th2[1] * 1e-12) * 1e12
                # transfer of a witness-like change (tau_3ph x (1 +- 0.1)): fitted delta change vs static change
                tr = {}
                for eta in (-0.1, 0.1):
                    sw = base_ops[eta]
                    bw = dict(base)
                    bw["tauR"], bw["taul"] = sw["tauR"], sw["taul"]
                    yw = a5.lnH(np.zeros(9), bw, sw["ratio"])
                    chiw, thw = a5.gn_refit(yw, base, nuis, sp, free=(1,))
                    dfit_w = (s0["tauR"] - s0["taul"] - thw[1] * 1e-12) * 1e12
                    d_static_w = (sw["tauR"] - sw["taul"]) * 1e12
                    # the same change seen through the first-order GK film (static times)
                    yw_gk = a5.lnH(np.zeros(9), bw)
                    chig, thg = a5.gn_refit(yw_gk, base, nuis, sp, free=(1,))
                    dfit_g = (s0["tauR"] - s0["taul"] - thg[1] * 1e-12) * 1e12
                    tr[str(eta)] = {"static_delta_change_ps": d_static_w - (s0["tauR"] - s0["taul"]) * 1e12,
                                    "fitted_delta_change_exact_ps": dfit_w - d_fit1,
                                    "fitted_delta_change_gk_ps": dfit_g - (s0["tauR"] - s0["taul"]) * 1e12,
                                    "transfer_exact_over_static": (dfit_w - d_fit1) /
                                    (d_static_w - (s0["tauR"] - s0["taul"]) * 1e12)}
                    # direct Delta chi^2 between the REF and perturbed exact-kernel data (nuisances refitted, film
                    # model = exact kernel of the REF operator, i.e. the distinguishability of the two operators)
                    chi_d, _ = a5.gn_refit(yw, base, nuis, sp, ratio=s0["ratio"])
                    chi_dg, _ = a5.gn_refit(yw_gk, base, nuis, sp)
                    tr[str(eta)]["dchi2_exact"] = chi_d
                    tr[str(eta)]["dchi2_firstorderGK"] = chi_dg
                row[f"{N}|{sp}"] = {"dchi2_gk_static": chi0, "dchi2_gk_taul_free": chi1, "delta_fit_taul_free_ps": d_fit1,
                                    "dchi2_gk_both_free": chi2, "delta_fit_both_free_ps": d_fit2,
                                    "delta_static_ps": (s0["tauR"] - s0["taul"]) * 1e12, "transfer": tr}
        ent[Rtag] = row
        print(L, Rtag, f"max dphase {row['max_abs_dphase_deg']:.4f} deg",
              {k: {"chi_static": round(v_["dchi2_gk_static"], 3), "delta_fit": round(v_["delta_fit_taul_free_ps"], 1),
                   "delta_static": round(v_["delta_static_ps"], 1),
                   "transfer": [round(v_["transfer"][x]["transfer_exact_over_static"], 3) for x in ("-0.1", "0.1")],
                   "dchi2 exact/GK(+0.1)": [round(v_["transfer"]["0.1"]["dchi2_exact"], 2),
                                            round(v_["transfer"]["0.1"]["dchi2_firstorderGK"], 2)]}
               for k, v_ in row.items() if isinstance(v_, dict)}, flush=True)
    res[str(L)] = ent
    (OUT / "au5b_ktrunc.json").write_text(json.dumps(res, indent=1))
