"""Multi-start inner approximation of [tau_min, tau_max] on F(r, K0) for a given geometry."""
from __future__ import annotations

import os
import time

import numpy as np

from optim import Slice


def _starts(S, g_ref, rng, n_logn, n_vert, spread=(0.7, 1.5, 3.0)):
    starts = [("ref", g_ref.copy())]
    for i in range(n_logn):
        sd = spread[i % len(spread)]
        starts.append((f"logn{sd}", g_ref * np.exp(rng.normal(0, sd, len(g_ref)))))
    vstarts = []
    if n_vert > 0:
        g_lo, K_lo = S.k_minimizer(g_ref)
        for i in range(4 * n_vert):
            if len(vstarts) >= n_vert:
                break
            v = S.lp_vertex(rng, sparse_cost=(i % 2 == 0))
            if v is None:
                continue
            g = S.segment_start(v, g_lo)
            if g is not None:
                vstarts.append(("vertex", g))
    return starts + vstarts


def run_interval(S, g_ref, seed=0, n_logn=10, n_vert=10, maxit=600, verbose=False):
    """g_ref is given in full event space; it is converted to the slice variables."""
    rng = np.random.default_rng(seed)
    t0 = time.time()
    gam_ref = S.gamma_of(g_ref)
    starts = _starts(S, gam_ref, rng, n_logn, n_vert)
    best = {+1: None, -1: None}
    log = []
    for name, gs in starts:
        for sign in (+1, -1):
            feas_start = name == "vertex"
            res = S.run(gs, sign, maxit=maxit, already_feasible=feas_start)
            if res is None:
                log.append((name, sign, None))
                continue
            if res["feas"] > 1e-9:
                log.append((name, sign, "infeasible"))
                continue
            log.append((name, sign, res["tau"]))
            cur = best[sign]
            if cur is None or sign * (res["tau"] - cur["tau"]) > 0:
                best[sign] = dict(res, start=name)
            if verbose:
                print(name, sign, res["tau"], res["iters"], flush=True)
    for key in (+1, -1):
        if best[key] is not None:
            best[key]["g_full"] = S.g_of(best[key]["g"])
    return {"max": best[+1], "min": best[-1], "log": log, "time": time.time() - t0,
            "n_starts": len(starts)}
