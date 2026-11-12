"""AU9: how well must the film conductivity be known for the delta range to be resolvable?

1D first-order GK film (au5 model), nuisances N2 + ln lam_z (C known) with a Gaussian prior of relative width s on
lam_z; sigma(delta)(s) at the regularised base point (tau_l mesh-corrected) and at the slab-fitted base point
(tau_l = tau_R - delta_fit, au6 nb=30). Margins against (a) the campaign half-range of delta, (b) the slab-visible
half-range (campaign half-range x witness-specific range transfer of au8, 9x9x6). s* = largest prior width with m <= 1.
Output: out/au9_lambda_prior.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au5_fisher1d as a5  # noqa: E402

AU4 = json.loads((OUT / "au4_delta_range.json").read_text())
AU6 = json.loads((OUT / "au6_slab_30.json").read_text())
AU8 = json.loads((OUT / "au8_witness_slab.json").read_text())
S_GRID = np.geomspace(1e-4, 0.3, 61)


def sigma_with_prior(J, nuis, sp, s_rel):
    Jw = a5.whiten(J, sp)
    v = Jw[:, 0]                                   # d/d tau_R at fixed tau_l = d/d delta (ps)
    N = Jw[:, nuis]
    # prior row on ln lam_z (column 2) with width s_rel: append a pseudo-observation
    col = nuis.index(2)
    prow = np.zeros((1, len(nuis)))
    prow[0, col] = 1.0 / s_rel
    A = np.vstack([np.column_stack([v, N]), np.concatenate([[0.0], prow[0]])[None, :]])
    F = A.T @ A
    cov = np.linalg.inv(F)
    return float(np.sqrt(cov[0, 0]))


res = {}
for L in (300, 500):
    b = a5.A1[str(L)]
    t4 = AU4["table"][f"{float(L):.6g}"]["all"]
    hr_c = t4["base_point"]["halfrange_delta_ps_campaign_rule"]
    Tw = AU8["m996"][f"{L}|note16"]["summary"]["range_transfer_N2lam"]
    for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
        for btag in ("regularised", "slab"):
            if btag == "regularised":
                taul = t4["base_point"]["taul0_ps"]
            else:
                taul = b["tau_R_ps"] - AU6[str(L)][f"wall|{Rtag}"]["N2lam|0.01"]["delta_fit_ps"]
            base = {"tauR": b["tau_R_ps"] * 1e-12, "taul": taul * 1e-12, "lamz": b["lam_z"], "C": a5.C_ALN,
                    "R1": R1, "R2": R2, "d": L * 1e-9}
            J = a5.jac(base)
            nuis = a5.SETS["N2"] + [2]
            for sp in (0.01, 0.1):
                sig = np.array([sigma_with_prior(J, nuis, sp, s) for s in S_GRID])
                hr = hr_c if btag == "regularised" else hr_c * Tw
                m = sig / hr
                ok = S_GRID[m <= 1.0]
                res[f"{L}|{Rtag}|{btag}|{sp}"] = {"halfrange_ps": hr, "sigma_N2_ps": float(a5.sigma_dir(J, (1, 0), a5.SETS["N2"], sp)),
                                                   "sigma_free_lam_ps": float(a5.sigma_dir(J, (1, 0), nuis, sp)),
                                                   "s_star_rel": float(ok.max()) if ok.size else None,
                                                   "m_at_s_1pct": float(np.interp(np.log(0.01), np.log(S_GRID), m)),
                                                   "m_at_s_5pct": float(np.interp(np.log(0.05), np.log(S_GRID), m))}
                print(L, Rtag, btag, sp, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[f"{L}|{Rtag}|{btag}|{sp}"].items()},
                      flush=True)
(OUT / "au9_lambda_prior.json").write_text(json.dumps(res, indent=1))
