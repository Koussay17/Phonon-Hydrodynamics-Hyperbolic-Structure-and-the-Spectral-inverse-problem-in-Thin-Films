"""Parallel scan of the inner tau_mem interval over (d, N, tie, alpha).

Task = (config, seed). Each task runs a block of starts (lognormal + LP-vertex) for both signs
and saves the best witnesses. Results: results/scan/<tag>__s<seed>.npz/.json.
Usage: python scan.py <tasks.json> [nproc]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "1")

import numpy as np

OUT = Path(__file__).resolve().parent.parent / "results" / "scan"


def setup(cfg):
    from debye_events import build
    from memgap import rta
    from optim import Slice
    from symmetry import debye_mode_maps, orbits, time_reversal_maps

    geom = build(cfg["d"], cfg["N"], speeds=tuple(cfg.get("speeds", (1.0, 2.0))),
                 kT=cfg.get("kT"), tol=cfg.get("tol"), amp_power=cfg.get("amp_power", 1.0))
    tie = None
    if cfg["tie"] == "tr":
        tie = orbits(geom, time_reversal_maps(geom))
    elif cfg["tie"] == "full":
        tie = orbits(geom, debye_mode_maps(geom))
    g_phys = geom.gphys
    r_phys = geom.W @ g_phys
    alpha = cfg.get("alpha")
    if alpha is None:
        g_ref, r = g_phys, r_phys
    else:
        r = geom.eps ** alpha
        r = r * (r_phys.sum() / r.sum())
        S0 = Slice(geom, r, 1.0, tie=tie)
        gam, ok = S0.project_diag(S0.gamma_of(g_phys))
        if not ok:
            return None
        g_ref = S0.g_of(gam)
    R = geom.response(g_ref)
    K0 = R["K"]
    S = Slice(geom, r, K0, tie=tie)
    KR, NR, tR = rta(geom, r)
    b = geom.b[:, 0]
    info = {"n": geom.n, "m": geom.m, "m_var": S.mvar, "n_rows": S.nr, "n_umklapp": int(geom.kind.sum()),
            "K0": K0, "tau_ref": R["tau"], "tau_RTA": tR, "K_RTA": KR, "CS_lower": K0 / float(b @ b),
            "b2": float(b @ b)}
    # emergent infrared exponent of r on the softest 30 % of each branch
    mb = geom.labels["modes_b"]
    al = []
    for br in np.unique(mb):
        sel = (mb == br) & (geom.eps <= np.quantile(geom.eps[mb == br], 0.3))
        if sel.sum() >= 3 and np.ptp(np.log(geom.eps[sel])) > 0:
            al.append(float(np.polyfit(np.log(geom.eps[sel]), np.log(r[sel]), 1)[0]))
    info["alpha_fit"] = al
    return geom, S, g_ref, r, info


def run_task(task):
    cfg, seed = task["cfg"], task["seed"]
    tag = task["tag"]
    out_json = OUT / f"{tag}__s{seed}.json"
    if out_json.exists():
        return str(out_json)
    from interval import run_interval
    t0 = time.time()
    st = setup(cfg)
    if st is None:
        out_json.write_text(json.dumps({"cfg": cfg, "seed": seed, "error": "infeasible r"}))
        return str(out_json)
    geom, S, g_ref, r, info = st
    res = run_interval(S, g_ref, seed=seed, n_logn=task.get("n_logn", 6), n_vert=task.get("n_vert", 6),
                       maxit=task.get("maxit", 1500))
    rec = {"cfg": cfg, "seed": seed, "info": info, "time": time.time() - t0, "n_starts": res["n_starts"]}
    for key in ("min", "max"):
        b = res[key]
        if b is None:
            rec[key] = None
            continue
        rec[key] = {"tau": b["tau"], "feas": b["feas"], "start": b["start"], "iters": b["iters"]}
    rec["log"] = [(l[0], l[1], l[2] if isinstance(l[2], float) else str(l[2])) for l in res["log"]]
    np.savez_compressed(OUT / f"{tag}__s{seed}.npz",
                        gmin=res["min"]["g_full"] if res["min"] else np.zeros(0),
                        gmax=res["max"]["g_full"] if res["max"] else np.zeros(0),
                        g_ref=g_ref, r=r)
    out_json.write_text(json.dumps(rec, indent=1))
    return str(out_json)


def main():
    tasks = json.loads(Path(sys.argv[1]).read_text())
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    OUT.mkdir(parents=True, exist_ok=True)
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, path in enumerate(pool.imap_unordered(run_task, tasks)):
            print(f"[{time.time()-t0:8.0f}s] {i+1}/{len(tasks)} {path}", flush=True)


if __name__ == "__main__":
    main()
