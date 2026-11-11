"""Task 3, calibrated class: per-cluster-orbit bounds set by the spread between two independent workflows.
eps_o = max(k * delta_o, floor), delta_o = ||(eta_ALAMODE - eta_0)_o|| / ||eta_0,o|| (alamode_compare.py),
k in {1, 2}, floor = 0.005; lifetimes +-10 %, kappa +-1 %, no Frobenius bound; otherwise as t3_batch.py.
Coordinates: ε-weighted scaling scaled_cal<k>_* (orbit_scaling.py --eps), in which every per-orbit constraint is
||(W u)_o|| <= 1 (a first attempt with the eps vector in the unweighted coordinates stalled: 09, item 13).
usage: python -B t3_calibrated.py <events.npz> <resultdir> <outdir> <nworkers> <alamode_orbits.npz> <k>
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

FLOOR = 0.005
TASKS = [(0, "max"), (2, "max"), (0, "min"), (2, "min")]
_M = {}


def _init(ev_path, rdir, outdir, orbf, kf):
    from optimize_fc3 import FC3ModelScaled
    _M["model"] = FC3ModelScaled(ev_path, rdir, tag=f"cal{kf}")
    _M["outdir"], _M["delta"], _M["k"] = Path(outdir), np.load(orbf)["rel"], kf


def _task(t):
    from optimize_fc3 import run
    from fc3class import FOURPI
    col, sense = t
    kf = _M["k"]
    model, outdir = _M["model"], _M["outdir"]
    eps = np.maximum(kf * _M["delta"], FLOOR)
    name = f"cal{kf}_c{col}_{sense}"
    if (outdir / f"{name}.json").exists():
        return name, "exists"
    t0 = time.time()
    best, hist = run(model, col, sense, np.inf, 0.10, 0.01, orbit_eps=1.0, maxouter=15, maxiter=200,
                     logfile=str(outdir / f"{name}.log"))
    out = {"name": name, "col": col, "sense": sense, "k": kf, "floor": FLOOR, "dr": 0.10, "dK": 0.01,
           "eps_quantiles(0,50,90,100%)": np.quantile(eps, [0, .5, .9, 1]).tolist(),
           "tau_base_ps": model.base["tau"][col] / FOURPI, "history": hist, "seconds": time.time() - t0}
    if best is not None:
        val, y, rec = best
        ev = model.evaluate(y)
        g, gb, al = ev["g_orb"], model.base["g_orb"], model.allowed
        ratio = g[al] / gb[al]
        rel = np.sqrt(model.orbit_rel2(y)[0]) * eps       # weighted coordinates: ||(W u)_o|| = rel_o / eps_o
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": ev["tau"][0] / FOURPI, "tau_z_ps": ev["tau"][2] / FOURPI,
                       "kappa_x_ratio": ev["K"][0] / model.K0[0], "kappa_z_ratio": ev["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(ev["r"] / model.r0 - 1))),
                       "orbit_rel_over_eps_max": float(np.max(rel / eps)), "orbit_rel_max": float(rel.max()),
                       "orbit_rate_ratio_min": float(ratio.min()), "orbit_rate_ratio_max": float(ratio.max()),
                       "n_orbits_suppressed_10x": int((ratio < 0.1).sum()), "record": rec}
        np.save(outdir / f"{name}_y.npy", y)
    (outdir / f"{name}.json").write_text(json.dumps(out, indent=1))
    return name, out.get("best", {}).get("tau_ps")


def main():
    ev_path, rdir, outdir, nw, orbf, kf = sys.argv[1:7]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    with Pool(int(nw), initializer=_init, initargs=(ev_path, rdir, outdir, orbf, int(kf))) as pool:
        for name, val in pool.imap_unordered(_task, TASKS, chunksize=1):
            print(time.strftime("%H:%M:%S"), name, val, flush=True)


if __name__ == "__main__":
    main()
