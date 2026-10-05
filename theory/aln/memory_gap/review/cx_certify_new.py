"""(copied certify() from cx_attack3_arb.py)
Attack 3: rigorous (arb ball) evaluation of witnesses with EXACT energy conservation, and
sensitivity of the near-singular witness to tiny perturbations of the harmonic data.

Differences from the campaign's certify_aln.py: the conserving coefficients l_p = sqrt((f_a+f_b)/f_p),
l_d = 1/l_p are computed in ball arithmetic from the (exact dyadic) float frequencies, so s^T eps = 0
holds exactly for the geometry defined by the float frequencies; the regularised solve
(C + gam e e^T) x = b therefore returns C^+ b exactly (up to the b.e ~ 1e-17 rounding of the
time-reversal symmetry of the velocities, whose effect is reported).
Perturbation test: frequencies f -> f (1 + 1e-10 xi), xi ~ U(-1,1), applied consistently to D, l_p, b, e
(the rates g are kept); tau_x recomputed in arb.
usage: python cx_attack3_arb.py [prec]
Output: review/cx_attack3_arb.json
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, their_ev_orb, load_witness  # noqa: E402

prec = int(sys.argv[1]) if len(sys.argv) > 1 else 256
ctx.prec = prec
M = AlN()
gref = M.g_ref()
ev_orb_theirs, _ = their_ev_orb()


def certify(g, eps_f, cols=(0, 2), exact_lp=True):
    n = M.n
    kT = arb(float(M.KB)) * arb(float(M.T)) / arb(float(M.THzToEv))
    eps = [arb(float(x)) for x in eps_f]
    dsc = []
    for e_ in eps:
        nocc = 1 / arb.expm1(e_ / kT)
        dsc.append((nocc * (1 + nocc)).sqrt())
    rows = [[arb(0)] * n for _ in range(n)]
    rows = [list(r) for r in rows]
    P, A, B, rep = M.P, M.A, M.B, M.rep
    for k in range(M.m):
        gk = float(g[k])
        if gk == 0.0:
            continue
        ga = arb(gk)
        p, a, b = int(P[k]), int(A[k]), int(B[k])
        if exact_lp:
            lp = ((eps[a] + eps[b]) / eps[p]).sqrt()
            ld = 1 / lp
        else:
            lp = arb(float(M.lp[k])); ld = arb(float(M.ld[k]))
        if rep[k]:
            ids = [p, a]; vals = [-lp / dsc[p], 2 * ld / dsc[a]]
        else:
            ids = [p, a, b]; vals = [-lp / dsc[p], ld / dsc[a], ld / dsc[b]]
        for i, vi in zip(ids, vals):
            for j, vj in zip(ids, vals):
                rows[i][j] = rows[i][j] + ga * vi * vj
    e = [dsc[i] * eps[i] for i in range(n)]
    en = sum((x * x for x in e), arb(0)).sqrt()
    eh = [x / en for x in e]
    gam = sum((rows[i][i] for i in range(n)), arb(0)) / n
    # energy residual of C (should contain 0 when exact_lp)
    Ce = [sum((rows[i][j] * e[j] for j in range(n)), arb(0)) for i in range(0, n, 97)]
    Creg = arb_mat([[rows[i][j] + gam * eh[i] * eh[j] for j in range(n)] for i in range(n)])
    Bm = arb_mat([[dsc[i] * eps[i] * arb(float(M.vel[i, c])) for c in cols] for i in range(n)])
    X = Creg.solve(Bm)
    out = {"Ce_sample_abs_mid_max": max(abs(float(v.mid())) for v in Ce),
           "Ce_sample_rad_max": max(float(v.rad()) for v in Ce),
           "Ce_sample_contains_zero": all(v.contains(0) for v in Ce)}
    for kk, c in enumerate(cols):
        K = sum((Bm[i, kk] * X[i, kk] for i in range(n)), arb(0))
        N = sum((X[i, kk] * X[i, kk] for i in range(n)), arb(0))
        be = sum((Bm[i, kk] * e[i] for i in range(n)), arb(0))
        tau = N / K
        out[c] = {"K_mid": float(K.mid()), "K_rad": float(K.rad()), "tau_ps_mid": float(tau.mid()) / (4 * np.pi),
                  "tau_ps_rad": float(tau.rad()) / (4 * np.pi), "b_dot_e_mid": float(be.mid())}
    return out



# ---- certification of the review's physical-class witnesses (forbidden events exactly zero)
import json as _json
o = M.orbits()
forb = gref < 1e-12 * gref.max()
RUNS = CAMP / "review" / "cx_runs"
res = {"prec": prec}
for jid in sys.argv[2:]:
    t0 = time.time()
    u = np.load(RUNS / f"{jid}_u.npy")
    g = gref * (~forb) * np.exp(u[o["ev_orb"]])
    r = certify(g, M.eps, cols=(0, 2), exact_lp=True)
    r["seconds"] = time.time() - t0
    mo = M.moments(gref, (0, 2))
    for c in (0, 2):
        r[c]["K_rel_to_K0"] = r[c]["K_mid"] / mo[c]["K"] - 1
    res[jid] = r
    print(jid, _json.dumps(r), flush=True)
    (CAMP / "review" / "cx_certify_new.json").write_text(_json.dumps(res, indent=1, default=str))
