"""Exploration: tau_mem extremes for the exact 1D two-branch Debye event set (small N).

SLSQP on g >= 0 with constraints W g = r, log K(g) = log K0 (r, K0 from the reference g).
"""
import sys
import time

import numpy as np
from scipy.optimize import minimize

from debye_events import build
from memgap import iproject, rta

N = int(sys.argv[1]) if len(sys.argv) > 1 else 9
nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 20
geom = build(1, N)
g0 = geom.gphys
ref = geom.response(g0)
r = geom.W @ g0
K0 = ref["K"]
KR, NR, tR = rta(geom, r)
print(f"N={N} n={geom.n} m={geom.m} K0={K0:.6g} tau_ref={ref['tau']:.6g} tau_RTA={tR:.6g} "
      f"K_RTA/K0={KR/K0:.4f} |b|^2={np.sum(geom.b[:,0]**2):.6g} CS lower={K0/np.sum(geom.b[:,0]**2):.6g}")

W = geom.W.toarray()
rs = r.copy()


def solve(sign, gstart):
    def f(g):
        res = geom.response(g, want_grad=True)
        return sign * np.log(res["N"]), sign * res["dN"] / res["N"]

    cons = [
        {"type": "eq", "fun": lambda g: (W @ g - r) / rs, "jac": lambda g: W / rs[:, None]},
        {"type": "eq", "fun": lambda g: np.array([np.log(geom.response(g)["K"] / K0)]),
         "jac": lambda g: (lambda R: (R["dK"] / R["K"])[None, :])(geom.response(g, want_grad=True))},
    ]
    try:
        out = minimize(f, gstart, jac=True, method="SLSQP", bounds=[(0, None)] * geom.m,
                       constraints=cons, options={"maxiter": 2000, "ftol": 1e-14})
    except Exception as exc:  # singular C
        return None, str(exc)
    g = np.maximum(out.x, 0)
    try:
        R = geom.response(g)
    except Exception as exc:
        return None, str(exc)
    feas = max(np.max(np.abs(W @ g - r) / rs), abs(R["K"] / K0 - 1))
    return (g, R["tau"], feas, out.success), out.message


rng = np.random.default_rng(1)
best = {+1: (np.inf, None), -1: (-np.inf, None)}
t0 = time.time()
for s in range(nstart):
    gs = g0 * np.exp(rng.normal(0, 1.5, geom.m))
    gs, info = iproject(geom, gs, r, K0=K0)
    if info["resid_inf"] > 1e-9:
        continue
    for sign in (+1, -1):
        res, msg = solve(sign, gs)
        if res is None:
            continue
        g, tau, feas, ok = res
        if feas > 1e-7:
            continue
        if sign == +1 and tau < best[+1][0]:
            best[+1] = (tau, g)
        if sign == -1 and tau > best[-1][0]:
            best[-1] = (tau, g)
print(f"time {time.time()-t0:.1f}s")
tmin, gmin = best[+1]
tmax, gmax = best[-1]
print(f"tau_min(found)={tmin:.6g}  tau_max(found)={tmax:.6g}  ratio={tmax/tmin:.4g}  tau_ref={ref['tau']:.6g}")
for name, g in (("min", gmin), ("max", gmax)):
    act = g > 1e-10 * g.max()
    print(name, "active events", int(act.sum()), "of", geom.m, " umklapp active", int((act & (geom.kind == 1)).sum()))
np.savez(f"../results/explore_1d_N{N}.npz", gmin=gmin, gmax=gmax, g0=g0, r=r, K0=K0)
