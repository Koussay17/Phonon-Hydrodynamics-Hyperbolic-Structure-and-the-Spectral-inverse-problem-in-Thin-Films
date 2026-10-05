"""Attack 3: rigorous (arb ball) evaluation of witnesses with EXACT energy conservation, and
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


res = {"prec": prec}
mo = M.moments(gref, (0, 2))
res["K0"] = {c: mo[c]["K"] for c in (0, 2)}
jobs = [("full_c0_max_04_logn3.0", True), ("full_c0_max_04_logn3.0", False), ("full_c0_max_08_vertex", True)]
for name, exact in jobs:
    t0 = time.time()
    g = load_witness(name, ev_orb_theirs)
    r = certify(g, M.eps, exact_lp=exact)
    r["seconds"] = time.time() - t0
    res[f"{name}|exact_lp={exact}"] = r
    print(name, exact, json.dumps(r), flush=True)
    (CAMP / "review" / "cx_attack3_arb.json").write_text(json.dumps(res, indent=1, default=str))

# perturbation test for the near-singular witness
rng = np.random.default_rng(11)
g = load_witness("full_c0_max_04_logn3.0", ev_orb_theirs)
for trial in range(3):
    xi = rng.uniform(-1, 1, M.n)
    eps_p = M.eps * (1 + 1e-10 * xi)
    t0 = time.time()
    r = certify(g, eps_p, cols=(0,), exact_lp=True)
    r["seconds"] = time.time() - t0
    res[f"near_singular_perturbed_1e-10_trial{trial}"] = r
    print("perturbed", trial, json.dumps(r), flush=True)
    (CAMP / "review" / "cx_attack3_arb.json").write_text(json.dumps(res, indent=1, default=str))
g8 = load_witness("full_c0_max_08_vertex", ev_orb_theirs)
xi = rng.uniform(-1, 1, M.n)
r = certify(g8, M.eps * (1 + 1e-10 * xi), cols=(0,), exact_lp=True)
res["robust_witness_perturbed_1e-10"] = r
print("robust perturbed", json.dumps(r), flush=True)
(CAMP / "review" / "cx_attack3_arb.json").write_text(json.dumps(res, indent=1, default=str))
print("done")
