"""Task 3 adversarial check: multistart for the per-cluster-orbit bounded maxima (eps = 0.1), lifetimes +-10 %,
kappa +-1 %, no Frobenius bound; same optimiser and settings as t3_batch.py, different initial points u0.

Starts (scaled coordinates u, per-orbit relative change ||(W u)_o||, W isometric):
  rand<s>  : w with an independent Gaussian direction in every cluster orbit block and relative size 0.9 eps U(0,1),
             u0 = W^T w, rescaled so that max_o ||(W u0)_o|| = 0.9 eps;
  aim      : the Frobenius-class witness of the same column (results of the D: campaign, 481.45 ps for x,
             17,121 ps for z), w = D^-1 B c clipped block by block to relative size 0.9 eps, u0 = W^T w, rescaled
             as above (keeps the witness's direction where it is small, caps it where it rewrites a cluster).
usage: python -B t3_multistart.py <events.npz> <resultdir> <outdir> <nworkers> <witness_x_y.npy> <witness_z_y.npy>
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

EPS = 0.1
TASKS = [(0, "aim"), (2, "aim"), (0, "rand1"), (2, "rand1"), (0, "rand2"), (2, "rand2")]
_M = {}


def _cap(model, u, eps):
    rel = np.sqrt(model.orbit_rel2(u)[0])
    return u * (0.9 * eps / rel.max())


def start_point(model, rdir, start, col, wit):
    W, cb, nb = model.Bmat, model.col_block, model.n_blocks
    if start.startswith("rand"):
        rng = np.random.default_rng(int(start[4:]) * 10 + col)
        w = rng.standard_normal(W.shape[0])
        bn = np.sqrt(np.bincount(cb, weights=w * w, minlength=nb))
        w *= (0.9 * EPS * rng.uniform(0, 1, nb) / bn)[cb]
    else:
        CO = np.load(rdir / "cluster_orbits.npz")
        Dcol = CO["block_refnorm"][CO["col_block"]]
        B = np.load(rdir / "basis_B.npy", mmap_mode="r")
        c0 = np.load(rdir / "basis_c0.npy")
        y = np.load(wit[col])
        w = (B @ (np.linalg.norm(c0) * y)) / Dcol
        bn = np.sqrt(np.bincount(cb, weights=w * w, minlength=nb))
        w *= np.minimum(1.0, 0.9 * EPS / np.maximum(bn, 1e-300))[cb]
    return _cap(model, W.T @ w, EPS)


def _init(ev_path, rdir, outdir, wx, wz):
    from optimize_fc3 import FC3ModelScaled
    _M["model"] = FC3ModelScaled(ev_path, rdir)
    _M["rdir"], _M["outdir"], _M["wit"] = Path(rdir), Path(outdir), {0: wx, 2: wz}


def _task(t):
    from optimize_fc3 import run
    from fc3class import FOURPI
    col, start = t
    model, outdir = _M["model"], _M["outdir"]
    name = f"ms_orb{EPS:g}_c{col}_max_{start}"
    if (outdir / f"{name}.json").exists():
        return name, "exists"
    t0 = time.time()
    u0 = start_point(model, _M["rdir"], start, col, _M["wit"])
    ev0 = model.evaluate(u0)
    best, hist = run(model, col, "max", np.inf, 0.10, 0.01, orbit_eps=EPS, maxouter=15, maxiter=200,
                     logfile=str(outdir / f"{name}.log"), y_init=u0)
    out = {"name": name, "col": col, "start": start, "orbit_eps": EPS, "dr": 0.10, "dK": 0.01,
           "tau_base_ps": model.base["tau"][col] / FOURPI, "tau_start_ps": ev0["tau"][col] / FOURPI,
           "start_orbit_rel_max": float(np.sqrt(np.max(model.orbit_rel2(u0)[0]))),
           "history": hist, "seconds": time.time() - t0}
    if best is not None:
        val, y, rec = best
        ev = model.evaluate(y)
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": ev["tau"][0] / FOURPI, "tau_z_ps": ev["tau"][2] / FOURPI,
                       "kappa_x_ratio": ev["K"][0] / model.K0[0], "kappa_z_ratio": ev["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(ev["r"] / model.r0 - 1))),
                       "orbit_rel_max": float(np.sqrt(np.max(model.orbit_rel2(y)[0]))), "record": rec}
        np.save(outdir / f"{name}_y.npy", y)
    (outdir / f"{name}.json").write_text(json.dumps(out, indent=1))
    return name, out.get("best", {}).get("tau_ps")


def main():
    ev_path, rdir, outdir, nw, wx, wz = sys.argv[1:7]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with Pool(int(nw), initializer=_init, initargs=(ev_path, rdir, outdir, wx, wz)) as pool:
        for name, val in pool.imap_unordered(_task, TASKS, chunksize=1):
            print(time.strftime("%H:%M:%S"), name, val, flush=True)


if __name__ == "__main__":
    main()
