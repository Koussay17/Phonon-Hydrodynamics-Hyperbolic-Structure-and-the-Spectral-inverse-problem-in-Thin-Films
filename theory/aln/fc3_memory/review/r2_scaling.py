"""Review R2: eps-weighted well-conditioned coordinates for a MODIFIED per-orbit bound vector (same construction as
orbit_scaling.py --eps, written to the review folder; campaign files untouched).

Variants of the calibrated k = 1 class:
  cal1nz   : eps_o = max(delta_o, 0.005) on orbits where the folded ALAMODE tensor is nonzero; eps_o = 0.005 (floor)
             on the 266 orbits where ALAMODE is identically zero (delta_o = 1 by its cubic cutoff, not a measured
             difference).
  cal1le5  : eps_o = max(delta_o, 0.005) for orbits with d_max < 5 A, floor beyond (no freedom in the 5-6 and > 6 A
             shells).
usage: python -B r2_scaling.py <rdir> <alamode_orbits.npz> <variant> <outdir>
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

rdir, orbf, variant, outdir = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
outdir.mkdir(parents=True, exist_ok=True)
t0 = time.perf_counter()
A = np.load(orbf)
rel, an, refn, dmax = A["rel"], A["alanorm"], A["refnorm"], A["dmax"]
FLOOR = 0.005
eps = np.maximum(rel, FLOOR)
if variant == "cal1nz":
    zero = an <= 1e-9 * refn
    eps[zero] = FLOOR
elif variant == "cal1le5":
    eps[dmax >= 5.0] = FLOOR
elif variant == "cal1":
    pass
else:
    raise ValueError(variant)
B = np.load(rdir / "basis_B.npy")
CO = np.load(rdir / "cluster_orbits.npz")
refb = CO["block_refnorm"] * eps
Dcol = refb[CO["col_block"]]
BD = B / Dcol[:, None]
del B
M = BD.T @ BD
lam, Q = np.linalg.eigh(M)
Mmh = (Q * lam ** -0.5) @ Q.T
del Q
W = BD @ Mmh
del BD
Z = np.load(rdir / "linearization.npz")
Gu = Z["G"] @ Mmh
pre = f"scaled_{variant}_"
np.save(outdir / f"{pre}Gr.npy", np.ascontiguousarray(Gu.real))
np.save(outdir / f"{pre}Gi.npy", np.ascontiguousarray(Gu.imag))
np.save(outdir / f"{pre}W.npy", W)
np.save(outdir / f"{pre}Mmh.npy", Mmh)
np.save(outdir / f"{pre}epsvec.npy", eps)
rng = np.random.default_rng(0)
v = rng.standard_normal(W.shape[1])
out = {"variant": variant, "n_floor": int((eps == FLOOR).sum()), "eps_quantiles": np.quantile(eps, [0, .5, .9, 1]).tolist(),
       "M_eig_minmax": [float(lam.min()), float(lam.max())], "cond_M": float(lam.max() / lam.min()),
       "isometry_check": float(abs(np.linalg.norm(W @ v) / np.linalg.norm(v) - 1)), "seconds": time.perf_counter() - t0}
(outdir / f"{pre}orbit_scaling.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
