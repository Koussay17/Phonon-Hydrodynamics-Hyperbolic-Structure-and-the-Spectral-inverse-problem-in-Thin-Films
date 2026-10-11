"""Review R2: re-runs of the campaign optimiser (optimize_fc3.run, unchanged) with other starts, larger inner
budgets and modified per-orbit bound vectors. Read-only use of the campaign scripts and matrices.

Coordinates: 'scaled_<tag>' matrices (orbit_scaling.py layout). tag '' = unweighted per-orbit relative bound
(eps scalar); tag 'cal1' etc. = eps-weighted (bound 1.0). Custom weighted scalings built by r2_scaling.py live in
<sdir> with prefix scaled_<tag>_.
Start specs (u-coordinates of the target scaling):
  zero
  scale|<u.npy>|<factor>           u0 = factor * u (same scaling)
  map|<yunscaled.npy>|<clip>       y (unscaled, c = |c0| y) -> w = D^-1 B c (per-orbit weighted), blocks clipped to
                                   <clip> (relative units of the target bound), u0 = W^T w, rescaled so that
                                   max_o ||(W u0)_o|| = clip * bound
usage: python -B r2_opt.py <events> <rdir> <sdir> <tag> <col> <sense> <bound> <start> <maxouter> <maxiter> <outname> [nthreads]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

nth = sys.argv[12] if len(sys.argv) > 12 else "1"
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = nth
sys.dont_write_bytecode = True
CAMP = Path(r"C:\Users\Koussay\ResearchLab\orchestration\campaigns\20261005-fc3-memory\scripts")
sys.path.insert(0, str(CAMP))

import numpy as np  # noqa: E402

from optimize_fc3 import FC3Model, run  # noqa: E402
from fc3class import FOURPI  # noqa: E402


class ScaledModel(FC3Model):
    def __init__(self, ev_path, rdir, sdir, tag):
        super().__init__(ev_path, rdir)
        sdir = Path(sdir)
        pre = f"scaled_{tag}_" if tag else "scaled_"
        self.Gr = np.load(sdir / f"{pre}Gr.npy")
        self.Gi = np.load(sdir / f"{pre}Gi.npy")
        self.n_var = self.Gr.shape[1]
        CO = np.load(Path(rdir) / "cluster_orbits.npz")
        self.col_block = CO["col_block"]
        self.n_blocks = int(self.col_block.max()) + 1
        self.block_refnorm2 = np.ones(self.n_blocks)
        self.Bmat = np.load(sdir / f"{pre}W.npy")
        self.c0n_true = self.c0n
        self.c0n = 1.0
        self.sdir, self.pre, self.rdir = sdir, pre, Path(rdir)

    def weights_col(self, bound_vec):
        CO = np.load(self.rdir / "cluster_orbits.npz")
        return (CO["block_refnorm"] * bound_vec)[CO["col_block"]]


def start_point(model, spec, bound, epsvec):
    if spec == "zero":
        return np.zeros(model.n_var)
    kind, path, fac = spec.split("|")
    fac = float(fac)
    if kind == "scale":
        return fac * np.load(path)
    if kind == "map":
        y = np.load(path)
        c0 = np.load(model.rdir / "basis_c0.npy")
        B = np.load(model.rdir / "basis_B.npy", mmap_mode="r")
        Dcol = model.weights_col(epsvec)
        w = (B @ (np.linalg.norm(c0) * y)) / Dcol
        cb, nb = model.col_block, model.n_blocks
        bn = np.sqrt(np.bincount(cb, weights=w * w, minlength=nb))
        w *= np.minimum(1.0, fac * bound / np.maximum(bn, 1e-300))[cb]
        u = model.Bmat.T @ w
        rel = np.sqrt(model.orbit_rel2(u)[0])
        return u * (fac * bound / rel.max())
    raise ValueError(spec)


def main():
    ev, rdir, sdir, tag, col, sense, bound, spec, maxouter, maxiter, outname = sys.argv[1:12]
    col, bound, maxouter, maxiter = int(col), float(bound), int(maxouter), int(maxiter)
    outname = Path(outname)
    outname.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    model = ScaledModel(ev, rdir, sdir, tag)
    if tag == "":
        epsvec = np.ones(model.n_blocks)          # unweighted: D = refnorm
    else:
        epsvec = np.load(Path(sdir) / f"scaled_{tag}_epsvec.npy")
    u0 = start_point(model, spec, bound, epsvec)
    ev0 = model.evaluate(u0)
    t_init = time.time() - t0
    best, hist = run(model, col, sense, np.inf, 0.10, 0.01, orbit_eps=bound, maxouter=maxouter, maxiter=maxiter,
                     logfile=str(outname.with_suffix(".log")), y_init=u0)
    out = {"args": sys.argv[1:], "tau_base_ps": model.base["tau"][col] / FOURPI,
           "tau_start_ps": ev0["tau"][col] / FOURPI,
           "start_orbit_rel_max": float(np.sqrt(np.max(model.orbit_rel2(u0)[0]))),
           "start_max_lifetime_dev": float(np.max(np.abs(ev0["r"] / model.r0 - 1))),
           "start_kappa_ratio": [ev0["K"][0] / model.K0[0], ev0["K"][2] / model.K0[2]],
           "t_init": t_init, "history": hist, "seconds": time.time() - t0}
    if best is not None:
        val, y, rec = best
        e = model.evaluate(y)
        g, gb, al = e["g_orb"], model.base["g_orb"], model.allowed
        ratio = g[al] / gb[al]
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": e["tau"][0] / FOURPI, "tau_z_ps": e["tau"][2] / FOURPI,
                       "kappa_x_ratio": e["K"][0] / model.K0[0], "kappa_z_ratio": e["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(e["r"] / model.r0 - 1))),
                       "orbit_rel_over_bound_max": float(np.sqrt(np.max(model.orbit_rel2(y)[0])) / bound),
                       "orbit_rate_ratio_min": float(ratio.min()), "orbit_rate_ratio_max": float(ratio.max()),
                       "n_orbits_suppressed_10x": int((ratio < 0.1).sum()), "record": rec}
        np.save(outname.with_name(outname.stem + "_y.npy"), y)
    outname.with_suffix(".json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != "history"}, indent=1), flush=True)


if __name__ == "__main__":
    main()
