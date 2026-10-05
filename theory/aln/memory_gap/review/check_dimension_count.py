"""Dimension count of the feasible set F(r, K0) for the AlN 5x5x3 event geometry (13-final-report.md,
"Dimension count").  Claimed: dim F = m_var - n_rows - n_K "when the constraints are independent":
2508 (symmetric rates: 2627 - 117 - 2) and 43,503 (unconstrained: 44,406 - 897 - 6).

Counting principle checked here: at a point g of the OPEN orthant (all rates > 0) where the constraint map
g -> (W g, K_cd(g)) has Jacobian J of rank rho, the feasible set is locally a manifold of dimension
m_var - rho (implicit function theorem).  So the claimed numbers are exact iff rho = n_rows + n_K,
and otherwise the true local dimension is LARGER.  We compute rho at the reference rates.

Read-only use of the campaign modules (no bytecode written).
Run: python -B check_dimension_count.py   (writes check_dimension_count.json)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))

import numpy as np  # noqa: E402

from aln_geometry import load  # noqa: E402
from symmetry import orbits  # noqa: E402

EVENTS = HERE.parent / "results" / "aln" / "events_m553_s0.1.npz"


def numerical_rank(Mx, label):
    s = np.linalg.svd(Mx, compute_uv=False)
    s = s / s[0]
    gaps = s[:-1] / np.maximum(s[1:], 1e-300)
    k = int(np.argmax(gaps))
    return {"label": label, "shape": list(Mx.shape), "sigma_rel_last5": s[-5:].tolist(),
            "sigma_rel_min": float(s[-1]), "largest_gap_after_index": k + 1,
            "largest_gap": float(gaps[k]),
            "rank_tol_1e-10": int(np.sum(s > 1e-10)), "rank_tol_1e-13": int(np.sum(s > 1e-13))}


def main():
    t0 = time.time()
    geom = load(str(EVENTS))
    out = {"n": geom.n, "m": geom.m}
    tie = orbits(geom, geom.labels["maps"])
    sizes = np.bincount(tie["ev_orb"])
    gref = (np.bincount(tie["ev_orb"], weights=geom.gphys) / sizes)[tie["ev_orb"]]
    out["min_gref"] = float(gref.min())
    out["n_event_orbits"] = int(tie["n_ev_orb"])
    out["n_mode_orbits"] = int(tie["n_mode_orb"])

    A = geom.A.tocsc()
    W = geom.W.toarray()                                   # n x m, entries a_{alpha,mu}^2
    r = W @ gref
    Wn = W / r[:, None]                                    # row scaling does not change the rank
    # DC gradients at the reference: dK_cd/dg_alpha = -(a.x_c)(a.x_d)
    X = {c: geom.response(gref, c)["x"] for c in (0, 1, 2)}
    P = {c: A.T @ X[c] for c in (0, 1, 2)}
    pairs_all = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    G = np.array([-(P[c] * P[d]) for (c, d) in pairs_all])
    K = {f"{c}{d}": float(geom.b[:, d] @ X[c]) for (c, d) in pairs_all}
    out["K_ref"] = K
    Gn = G / np.linalg.norm(G, axis=1)[:, None]

    # unconstrained rates
    out["W_full"] = numerical_rank(Wn, "W (897 x 44406)")
    out["W_plus_6K"] = numerical_rank(np.vstack([Wn, Gn]), "[W; grad K_ab] (903 x 44406)")

    # symmetric rates (orbit variables), one row per mode orbit, K_xx and K_zz
    Pm = tie["P"]                                          # (n_orb x m), 0/1
    reps = tie["reps"]
    Wg = (Pm @ W[reps].T).T                                # (117 x 2627)
    Wgn = Wg / (Wg @ (Pm @ gref / sizes))[:, None]
    Gg = np.array([Pm @ G[0], Pm @ G[2]])
    Ggn = Gg / np.linalg.norm(Gg, axis=1)[:, None]
    out["W_sym"] = numerical_rank(Wgn, "W_gamma (117 x 2627)")
    out["W_sym_plus_2K"] = numerical_rank(np.vstack([Wgn, Ggn]), "[W_gamma; grad K_xx, K_zz] (119 x 2627)")
    # consistency of the symmetric reduction: symmetric g gives identical rows within each mode orbit
    rsym = W @ gref
    spread = 0.0
    for o in range(tie["n_mode_orb"]):
        idx = np.nonzero(tie["mode_orb"] == o)[0]
        spread = max(spread, float(np.ptp(rsym[idx]) / rsym[idx].mean()))
    out["max_rel_spread_of_r_within_mode_orbits"] = spread
    rho_full = out["W_plus_6K"]["rank_tol_1e-10"]
    rho_sym = out["W_sym_plus_2K"]["rank_tol_1e-10"]
    out["dim_F_unconstrained"] = geom.m - rho_full
    out["dim_F_symmetric"] = int(tie["n_ev_orb"]) - rho_sym
    out["seconds"] = time.time() - t0
    (HERE / "check_dimension_count.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
