import sys
import time

import numpy as np

from debye_events import build
from optim import Slice

d = int(sys.argv[1]); N = int(sys.argv[2]); nstart = int(sys.argv[3]) if len(sys.argv) > 3 else 10
spread = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
geom = build(d, N)
g0 = geom.gphys
ref = geom.response(g0)
r = geom.W @ g0
K0 = ref["K"]
S = Slice(geom, r, K0)
rng = np.random.default_rng(0)
t0 = time.time()
out = {+1: [], -1: []}
for s in range(nstart):
    gs = g0 * np.exp(rng.normal(0, spread, geom.m)) if s > 0 else g0.copy()
    for sign in (+1, -1):
        res = S.run(gs, sign, maxit=3000)
        if res is None:
            continue
        out[sign].append((res["tau"], res["feas"], res["iters"]))
print(f"d={d} N={N} n={geom.n} m={geom.m} tau_ref={ref['tau']:.6g} time={time.time()-t0:.1f}s")
for sign in (+1, -1):
    vals = sorted(out[sign], key=lambda t: t[0], reverse=(sign > 0))
    print("max" if sign > 0 else "min", [f"{v[0]:.5g}({v[2]})" for v in vals[:10]], "feas", max(v[1] for v in vals))
