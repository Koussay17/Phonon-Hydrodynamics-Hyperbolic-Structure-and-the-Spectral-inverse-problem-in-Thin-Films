"""Ball-arithmetic certification (python-flint arb, 256 bits) of AlN witnesses: tau, K_xx, K_zz, diag residual."""
import json, sys, time
sys.dont_write_bytecode = True
import numpy as np
from flint import arb, arb_mat, ctx
from aln_scan import setup
from aln_geometry import tau_ps
from certify import model_balls

def cert(geom, g, cols, prec=256):
    ctx.prec = prec
    eps, dsc = model_balls(geom, prec)
    n = geom.n
    rows = [[arb(0)] * n for _ in range(n)]
    rows = [list(r) for r in rows]
    for a in range(geom.m):
        if float(g[a]) == 0.0:
            continue
        ga = arb(float(g[a]))
        ids, vals = [], []
        for s in range(3):
            mu = int(geom.events[a, s])
            if mu >= 0:
                ids.append(mu); vals.append(arb(float(geom.coef[a, s])) / dsc[mu])
        for i, vi in zip(ids, vals):
            for j, vj in zip(ids, vals):
                rows[i][j] = rows[i][j] + ga * vi * vj
    e = [dsc[i] * eps[i] for i in range(n)]
    en = sum((x * x for x in e), arb(0)).sqrt()
    eh = [x / en for x in e]
    gam = sum((rows[i][i] for i in range(n)), arb(0)) / n
    C = arb_mat([[rows[i][j] + gam * eh[i] * eh[j] for j in range(n)] for i in range(n)])
    B = arb_mat([[dsc[i] * eps[i] * arb(float(geom.vel[i, c])) for c in cols] for i in range(n)])
    X = C.solve(B)
    out = {}
    for k, c in enumerate(cols):
        K = sum((B[i, k] * X[i, k] for i in range(n)), arb(0))
        N = sum((X[i, k] * X[i, k] for i in range(n)), arb(0))
        t = N / K
        out[c] = {"K_mid": float(K.mid()), "K_rad": float(K.rad()), "tau_mid": float(t.mid()), "tau_rad": float(t.rad()),
                  "tau_ps_mid": float(t.mid()) / (4 * np.pi), "tau_ps_rad": float(t.rad()) / (4 * np.pi)}
    return out

if __name__ == "__main__":
    files = sys.argv[1:]
    geom, gref, r, K0, slices = setup('../results/aln/events_m553_s0.1.npz', 'full')
    S = slices[0]
    res = {}
    for f in files:
        t0 = time.time()
        g = S.g_of(np.load(f)['gamma'])
        c = cert(geom, g, [0, 2])
        c["K0"] = {str(k): v for k, v in K0.items()}
        c["seconds"] = time.time() - t0
        res[f] = c
        print(f, json.dumps(c), flush=True)
        json.dump(res, open('../results/aln/certify_arb_m553.json', 'w'), indent=1)
