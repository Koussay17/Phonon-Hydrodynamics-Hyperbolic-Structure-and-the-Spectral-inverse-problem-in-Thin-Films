"""Parallel AlN interval runs (5x5x3 event geometry at 300 K, or another export).

Task: (tie in {full, none}, objective column o in {0 (x, basal), 2 (z, c axis)}, sign, start spec).
Fixed data: diag C = r (all mode lifetimes), and the DC tensor components
  full tie : K_xx, K_zz (K_yy = K_xx and off-diagonals vanish by 6mm x time-reversal symmetry),
  none     : all six components K_ab.
Reference rates: orbit average (6mm x time reversal, 24 operations) of the phono3py-derived rates.
Outputs: results/aln/scan_<label>/<task>.json and .npz (witness gamma).
usage: python aln_scan.py <events.npz> <label> <nproc> [maxit]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "1")
sys.dont_write_bytecode = True

import numpy as np

HERE = Path(__file__).resolve().parent

_CACHE = {}


def setup(events_path, tie_kind):
    key = (events_path, tie_kind)
    if key in _CACHE:
        return _CACHE[key]
    from aln_geometry import load
    from symmetry import orbits
    from optim_multi import MultiSlice
    geom = load(events_path)
    tie_full = orbits(geom, geom.labels["maps"])
    sizes = np.bincount(tie_full["ev_orb"])
    gref = (np.bincount(tie_full["ev_orb"], weights=geom.gphys) / sizes)[tie_full["ev_orb"]]
    r = geom.W @ gref
    if tie_kind == "full":
        kp = [(0, 0), (2, 2)]
        tie = tie_full
    else:
        kp = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        tie = None
    Kfull = {}
    for (c, d) in kp:
        xc = geom.response(gref, c)["x"]
        Kfull[(c, d)] = float(geom.b[:, d] @ xc)
    slices = {o: MultiSlice(geom, r, Kfull, kp, col=o, tie=tie) for o in (0, 2)}
    _CACHE[key] = (geom, gref, r, Kfull, slices)
    return _CACHE[key]


def run_task(task):
    out = Path(task["outdir"]) / f"{task['name']}.json"
    if out.exists():
        return str(out)
    from aln_geometry import tau_ps, kappa_factor
    geom, gref, r, K0, slices = setup(task["events"], task["tie"])
    S = slices[task["col"]]
    gam_ref = S.gamma_of(gref)
    rng = np.random.default_rng(task["seed"])
    t0 = time.time()
    kind = task["start"]
    start_ok = True
    feas_start = False
    if kind == "ref":
        gs = gam_ref
    elif kind.startswith("logn"):
        gs = gam_ref * np.exp(rng.normal(0, float(kind[4:]), len(gam_ref)))
    elif kind == "vertex":
        v = S.lp_vertex(rng, sparse_cost=(task["seed"] % 2 == 0))
        gs = S.mix_start(v, gam_ref) if v is not None else None
        feas_start = True
        start_ok = gs is not None
    if not start_ok:
        res = None
    else:
        res = S.run(gs, task["sign"], maxit=task["maxit"], already_feasible=feas_start)
    rec = {"task": task, "time": time.time() - t0}
    if res is None:
        rec["result"] = None
    else:
        kf = kappa_factor(geom)
        rec["result"] = {"tau_thz_units": res["tau"], "tau_ps": tau_ps(res["tau"]), "feas": res["feas"],
                         "iters": res["iters"], "K": {f"{c}{d}": v for (c, d), v in res["K"].items()},
                         "kappa_obj_WmK": kf * float(geom.b[:, task['col']] @ (S.evaluate(res['g'], grad=False)['X'][:, S.evaluate(res['g'], grad=False)['ci'][task['col']]]))}
        np.savez_compressed(Path(task["outdir"]) / f"{task['name']}.npz", gamma=res["g"])
    out.write_text(json.dumps(rec, indent=1))
    return str(out)


def main():
    events, label, nproc = sys.argv[1], sys.argv[2], int(sys.argv[3])
    maxit = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    outdir = HERE.parent / "results" / "aln" / f"scan_{label}"
    outdir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for tie, starts in (("full", ["ref", "logn0.5", "logn1.0", "logn2.0", "logn3.0"] + ["vertex"] * 5),
                        ("none", ["ref", "logn1.0", "logn2.0", "vertex", "vertex"])):
        for col in (0, 2):
            for sign in (+1, -1):
                for i, st in enumerate(starts):
                    name = f"{tie}_c{col}_{'max' if sign > 0 else 'min'}_{i:02d}_{st}"
                    tasks.append({"events": str(events), "tie": tie, "col": col, "sign": sign, "start": st,
                                  "seed": 1000 * col + 10 * i + (sign > 0), "maxit": maxit, "outdir": str(outdir),
                                  "name": name})
    # interleave so that every (tie, col, sign) gets early coverage
    tasks.sort(key=lambda t: (t["name"].split("_")[3], t["tie"] != "full"))
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, p in enumerate(pool.imap_unordered(run_task, tasks)):
            print(f"[{time.time()-t0:8.0f}s] {i+1}/{len(tasks)} {p}", flush=True)


if __name__ == "__main__":
    main()
