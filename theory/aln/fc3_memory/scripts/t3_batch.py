"""Task 3 batch driver: constrained tau_mem extremisation over symmetric fc3 perturbations (optimize_fc3.run).

Each worker loads the model once. Tasks: (col, sense, rho, dr, dK, seed, init_scale).
usage: python -B t3_batch.py <events.npz> <resultdir> <outdir> <nworkers> <plan>
plans: test | main | dr05 | restarts
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

_MODEL = None
_ARGS = None


def _init(ev_path, rdir, outdir, subspace, orbit_basisdir=None):
    global _MODEL, _ARGS
    from optimize_fc3 import FC3Model
    _MODEL = FC3Model(ev_path, rdir, subspace=subspace)
    if orbit_basisdir:
        _MODEL.load_orbit_structure(orbit_basisdir)
    _ARGS = (ev_path, rdir, outdir)


def _task(t):
    from optimize_fc3 import run
    from fc3class import FOURPI
    col, sense, rho, dr, dK, seed, init_scale = t[:7]
    oeps = t[7] if len(t) > 7 else None
    outdir = Path(_ARGS[2])
    name = f"c{col}_{sense}_rho{rho:g}_dr{dr:g}_dK{dK:g}_s{seed}_i{init_scale:g}" + (f"_orb{oeps:g}" if oeps else "")
    fn = outdir / f"{name}.json"
    if fn.exists():
        return name, "exists"
    t0 = time.time()
    model = _MODEL
    best, hist = run(model, col, sense, rho, dr, dK, seed, init_scale, log=False, orbit_eps=oeps,
                     maxouter=int(os.environ.get("FC3_MAXOUTER", "25")),
                     maxiter=int(os.environ.get("FC3_MAXITER", "300")), logfile=str(outdir / f"{name}.log"))
    out = {"name": name, "col": col, "sense": sense, "rho": rho, "dr": dr, "dK": dK, "seed": seed,
           "init_scale": init_scale, "orbit_eps": oeps, "tau_base_ps": model.base["tau"][col] / FOURPI,
           "history": hist, "seconds": time.time() - t0}
    if best is not None:
        val, y, rec = best
        ev = model.evaluate(y, cols=(0, 2))
        g, gb, al = ev["g_orb"], model.base["g_orb"], model.allowed
        ratio = g[al] / gb[al]
        out["best"] = {"tau_ps": val / FOURPI, "tau_x_ps": ev["tau"][0] / FOURPI, "tau_z_ps": ev["tau"][2] / FOURPI,
                       "kappa_x_ratio": ev["K"][0] / model.K0[0], "kappa_z_ratio": ev["K"][2] / model.K0[2],
                       "max_lifetime_dev": float(np.max(np.abs(ev["r"] / model.r0 - 1))),
                       "rho": float(np.linalg.norm(y)),
                       "orbit_rate_ratio_min": float(ratio.min()), "orbit_rate_ratio_max": float(ratio.max()),
                       "n_orbits_suppressed_10x": int((ratio < 0.1).sum()),
                       "n_orbits_suppressed_100x": int((ratio < 0.01).sum()), "record": rec}
        np.save(outdir / f"{name}_y.npy", y)
    fn.write_text(json.dumps(out, indent=1))
    return name, out.get("best", {}).get("tau_ps")


def plan(name):
    T = []
    inf = float("inf")
    if name.startswith("cut"):          # cut35, cut50, ...: restricted subspace, dr = 10 %
        for rho in (0.05, 0.1, 0.2, 0.5, inf):
            for col in (0, 2):
                for sense in ("max", "min"):
                    T.append((col, sense, rho, 0.10, 0.01, 0, 0.0))
        return T
    if name == "m775":                  # key runs on the 7x7x5 mesh
        return [(0, "max", inf, 0.10, 0.01, 0, 0.0, 0.1), (2, "max", inf, 0.10, 0.01, 0, 0.0, 0.1),
                (0, "max", 0.1, 0.10, 0.01, 0, 0.0)]
    if name == "main_red":              # reduced Frobenius plan: small rho, maxima
        return [(c, "max", r, 0.10, 0.01, 0, 0.0) for r in (0.02, 0.05) for c in (0, 2)] + [(2, "max", 0.1, 0.10, 0.01, 0, 0.0)]
    if name == "orb2":                  # prioritised per-cluster plan (executed in this order)
        return [(0, "max", inf, 0.10, 0.01, 0, 0.0, 0.1), (2, "max", inf, 0.10, 0.01, 0, 0.0, 0.1),
                (0, "max", inf, 0.10, 0.01, 0, 0.0, 0.05), (2, "max", inf, 0.10, 0.01, 0, 0.0, 0.05),
                (0, "max", inf, 0.10, 0.01, 0, 0.0, 0.2), (2, "max", inf, 0.10, 0.01, 0, 0.0, 0.2),
                (0, "min", inf, 0.10, 0.01, 0, 0.0, 0.1), (2, "min", inf, 0.10, 0.01, 0, 0.0, 0.1)]
    if name.startswith("orb"):          # per-cluster-orbit relative bounds, full space, no Frobenius bound
        eps_list = {"orb": (0.05, 0.1, 0.2, 0.5), "orb_small": (0.01, 0.02), "orb_dr05": (0.1, 0.2)}[name]
        drv = 0.05 if name == "orb_dr05" else 0.10
        for eps in eps_list:
            for col in (0, 2):
                for sense in ("max", "min"):
                    T.append((col, sense, inf, drv, 0.01, 0, 0.0, eps))
        return T
    if name == "test":
        T = [(0, "max", 0.1, 0.1, 0.01, 0, 0.0)]
    elif name == "main":
        for rho in (0.02, 0.05, 0.1, 0.2, 0.5, inf):
            for col in (0, 2):
                for sense in ("max", "min"):
                    T.append((col, sense, rho, 0.10, 0.01, 0, 0.0))
    elif name == "main_small":
        for rho in (0.005, 0.01):
            for col in (0, 2):
                for sense in ("max", "min"):
                    T.append((col, sense, rho, 0.10, 0.01, 0, 0.0))
    elif name == "dr05":
        for rho in (0.05, 0.2, inf):
            for col in (0, 2):
                for sense in ("max", "min"):
                    T.append((col, sense, rho, 0.05, 0.01, 0, 0.0))
    elif name == "restarts":
        for rho in (0.1, 0.2, inf):
            for col in (0, 2):
                for seed in (1, 2):
                    T.append((col, "max", rho, 0.10, 0.01, seed, 0.5))
    return T


def main():
    ev_path, rdir, outdir, nw, pl = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
    subspace = {"cut35": 3.5, "cut40": 4.0, "cut45": 4.5, "cut50": 5.0}.get(pl)
    orbit_basisdir = os.environ.get("FC3_BASISDIR", rdir) if (pl.startswith("orb") or pl == "m775") else None
    keep_order = pl in ("orb2", "m775")
    Path(outdir).mkdir(parents=True, exist_ok=True)
    tasks = plan(pl)
    # longest first: unbounded and large rho
    if not keep_order:
        tasks.sort(key=lambda t: -((t[2] if np.isfinite(t[2]) else 1e9) + (t[7] if len(t) > 7 else 0)))
    with Pool(nw, initializer=_init, initargs=(ev_path, rdir, outdir, subspace, orbit_basisdir)) as pool:
        for name, val in pool.imap_unordered(_task, tasks, chunksize=1):
            print(time.strftime("%H:%M:%S"), name, val, flush=True)


if __name__ == "__main__":
    main()
