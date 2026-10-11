"""Task 3 batch driver (C: mirror): per-cluster-orbit bounded runs, executed in priority order.

Tasks (col, sense, eps): lifetimes +-10 %, kappa +-1 %, no Frobenius bound, every cluster orbit changed by at most
eps relative to its reference norm. Outputs <outdir>/<name>.json, _y.npy, .log (one line per outer iteration).
usage: python -B t3_batch.py <events.npz> <resultdir> <outdir> <nworkers> [maxouter maxiter]
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

_M = {}
TASKS = [(0, "max", 0.1), (2, "max", 0.1), (0, "max", 0.05), (2, "max", 0.05),
         (0, "max", 0.2), (2, "max", 0.2), (0, "min", 0.1), (2, "min", 0.1)]


def _init(ev_path, rdir, outdir, maxouter, maxiter):
    from optimize_fc3 import FC3ModelScaled
    _M["model"] = FC3ModelScaled(ev_path, rdir)
    _M["outdir"], _M["maxouter"], _M["maxiter"] = Path(outdir), maxouter, maxiter


def _task(t):
    from optimize_fc3 import run
    from fc3class import FOURPI
    col, sense, eps = t
    model, outdir = _M["model"], _M["outdir"]
    name = f"orb{eps:g}_c{col}_{sense}"
    if (outdir / f"{name}.json").exists():
        return name, "exists"
    t0 = time.time()
    best, hist = run(model, col, sense, np.inf, 0.10, 0.01, orbit_eps=eps, maxouter=_M["maxouter"],
                     maxiter=_M["maxiter"], logfile=str(outdir / f"{name}.log"))
    out = {"name": name, "col": col, "sense": sense, "orbit_eps": eps, "dr": 0.10, "dK": 0.01,
           "tau_base_ps": model.base["tau"][col] / FOURPI, "history": hist, "seconds": time.time() - t0}
    if best is not None:
        val, y, rec = best
        ev = model.evaluate(y)
        g, gb, al = ev["g_orb"], model.base["g_orb"], model.allowed
        ratio = g[al] / gb[al]
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": ev["tau"][0] / FOURPI, "tau_z_ps": ev["tau"][2] / FOURPI,
                       "kappa_x_ratio": ev["K"][0] / model.K0[0], "kappa_z_ratio": ev["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(ev["r"] / model.r0 - 1))), "u_norm": float(np.linalg.norm(y)),
                       "orbit_rel_max": float(np.sqrt(np.max(model.orbit_rel2(y)[0]))),
                       "orbit_rate_ratio_min": float(ratio.min()), "orbit_rate_ratio_max": float(ratio.max()),
                       "n_orbits_suppressed_10x": int((ratio < 0.1).sum()), "record": rec}
        np.save(outdir / f"{name}_y.npy", y)
    (outdir / f"{name}.json").write_text(json.dumps(out, indent=1))
    return name, out.get("best", {}).get("tau_ps")


def main():
    ev_path, rdir, outdir, nw = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    maxouter = int(sys.argv[5]) if len(sys.argv) > 5 else 15
    maxiter = int(sys.argv[6]) if len(sys.argv) > 6 else 200
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with Pool(nw, initializer=_init, initargs=(ev_path, rdir, outdir, maxouter, maxiter)) as pool:
        for name, val in pool.imap_unordered(_task, TASKS, chunksize=1):
            print(time.strftime("%H:%M:%S"), name, val, flush=True)


if __name__ == "__main__":
    main()
