"""Spot-check of the E2 lower bounds on M (07-experiments.md E2: "each value is an explicit basic solution
(certified lower bound)").  results/bt_local_search.json stores the values but not the subsets I or the
conditioning of [a_I^T; b^T], so the word "certified" is not substantiated by the stored data.
Here: re-run memgap.bt_local_search on small geometries, record cond(B_I) of the best basis and re-solve the
basis system at 50 digits (mpmath) from the same float data.
Run: python -B check_pivot_bases.py  (writes check_pivot_bases.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

import mpmath  # noqa: E402
import numpy as np  # noqa: E402

from debye_events import build  # noqa: E402
from memgap import bt_local_search, reduced_basis  # noqa: E402


def verify(geom, I):
    Q = reduced_basis(geom)
    Ared = Q.T @ geom.A.toarray()
    bred = Q.T @ geom.b[:, 0]
    B = np.vstack([Ared[:, I].T, bred[None, :]])
    cond = float(np.linalg.cond(B))
    y = np.linalg.solve(B, np.eye(B.shape[0])[:, -1])
    mpmath.mp.dps = 50
    Bm = mpmath.matrix(B.tolist())
    rhs = mpmath.matrix([0] * (B.shape[0] - 1) + [1])
    ym = mpmath.lu_solve(Bm, rhs)
    nm = mpmath.sqrt(sum(v ** 2 for v in ym))
    return {"cond_B": cond, "norm_float": float(np.linalg.norm(y)), "norm_50digits": mpmath.nstr(nm, 20),
            "rel_diff": float(abs(nm - np.linalg.norm(y)) / nm)}


def main():
    out = {}
    for (d, N, starts) in ((1, 9, 40), (2, 3, 40), (3, 3, 12)):
        geom = build(d, N)
        res = bt_local_search(geom, n_starts=starts, rng=0)
        best, I = res[0]
        out[f"d{d}_N{N}"] = {"M_lower_rerun": best, "subset": [int(v) for v in I], **verify(geom, I),
                             "top5": [r[0] for r in res[:5]]}
    (HERE / "check_pivot_bases.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
