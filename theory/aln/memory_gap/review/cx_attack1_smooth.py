"""Attack 1e: smooth-amplitude prior (the class named in 13-final-report section 4 as missing).

Rates g = g_ref * exp(P(w_p, w_a, w_b)), P a polynomial of total degree <= k in the scaled frequencies
x = 2 w / w_max - 1, symmetric under a <-> b (so automatically constant on symmetry orbits).
Data fixed as in the campaign: all lifetimes (117 orbit rows) and K_xx, K_zz. Maximise/minimise tau.
Optional box |P| <= ln F on every event (linear inequalities). SLSQP in coefficient space.
usage: python cx_attack1_smooth.py k [F or 0] [col] [sign]
Output: review/cx_runs/smooth_k{k}_F{F}_c{col}_{max|min}.json
"""
import itertools
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP  # noqa: E402
from cx_opt import Block, Problem  # noqa: E402

k = int(sys.argv[1])
F = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
col = int(sys.argv[3]) if len(sys.argv) > 3 else 0
sign = int(sys.argv[4]) if len(sys.argv) > 4 else 1

M = AlN()
o = M.orbits()
gref = M.g_ref()
r = M.diag(gref)
mo = M.moments(gref, (0, 2))
K0 = {c: mo[c]["K"] for c in (0, 2)}
wmax = M.eps.max()
# orbit representative event -> frequencies
rep_ev = np.array([np.nonzero(o["ev_orb"] == q)[0][0] for q in range(o["n_ev_orb"])])
xp = 2 * M.eps[M.P[rep_ev]] / wmax - 1
xa = 2 * M.eps[M.A[rep_ev]] / wmax - 1
xb = 2 * M.eps[M.B[rep_ev]] / wmax - 1
cols = []
names = []
for i in range(k + 1):
    for j in range(k + 1 - i):
        for l in range(j, k + 1 - i - j):
            v = xp ** i * (xa ** j * xb ** l + xa ** l * xb ** j)
            cols.append(v)
            names.append((i, j, l))
Bm = np.array(cols).T                      # (n_orb, n_coef)
# orthonormalise columns for conditioning (keeps the span)
Qb, Rb = np.linalg.qr(Bm)
rank = int(np.sum(np.abs(np.diag(Rb)) > 1e-10 * np.abs(Rb[0, 0])))
Bq = Qb[:, :rank] * np.sqrt(len(xp))       # O(1) entries
print("degree", k, "monomials", len(names), "rank", rank, flush=True)

blk = Block(M, gref, r[o["mode_reps"]], K0, kcols=(0, 2))
P = Problem([blk], o["ev_orb"], o["n_ev_orb"], o["mode_reps"], -np.inf, np.inf, obj_col=col, sign=sign)
cache = {}


def ev(c):
    key = c.tobytes()
    if key not in cache:
        cache.clear()
        try:
            cache[key] = P.evaluate(Bq @ c)
        except Exception:  # failed factorisation: large residual so that SLSQP backtracks
            q = len(P.rows) + 2
            cache[key] = {"c": np.full(q, 1e3), "J": np.zeros((q, P.p)), "f": -1e3, "grad": np.zeros(P.p),
                          "tau_ps": float("nan")}
    return cache[key]


cons = [{"type": "eq", "fun": lambda c: ev(c)["c"], "jac": lambda c: ev(c)["J"] @ Bq}]
L = np.log(F) if F > 0 else 30.0          # 'unbounded' runs keep |P| <= 30 so that C stays positive definite
if True:
    cons.append({"type": "ineq", "fun": lambda c: np.concatenate([L - Bq @ c, L + Bq @ c]),
                 "jac": lambda c: np.vstack([-Bq, Bq])})
t0 = time.time()
hist = []


def cb(c):
    hist.append(float(np.exp(sign * ev(c)["f"]) / (4 * np.pi)))


res = minimize(lambda c: -ev(c)["f"], np.zeros(rank), jac=lambda c: -(Bq.T @ ev(c)["grad"]), constraints=cons,
               method="SLSQP", options={"maxiter": 400, "ftol": 1e-12}, callback=cb)
c = res.x
I = P.evaluate(Bq @ c)
u = Bq @ c
mo2 = M.moments(P.rates(u, 0), (0, 2))
out = {"degree": k, "n_coef": rank, "F": F, "col": col, "sign": sign, "message": res.message, "nit": int(res.nit),
       "tau_ps": I["tau_ps"], "tau_x_ps": mo2[0]["tau_ps"], "tau_z_ps": mo2[2]["tau_ps"],
       "c_max": float(np.max(np.abs(I["c"]))), "u_min": float(u.min()), "u_max": float(u.max()),
       "frac_orbits_suppressed_below_1e-2": float(np.mean(u < np.log(1e-2))), "seconds": time.time() - t0,
       "hist_tail": hist[-5:]}
print(json.dumps(out), flush=True)
(CAMP / "review" / "cx_runs").mkdir(exist_ok=True)
tag = f"smooth_k{k}_F{int(F)}_c{col}_{'max' if sign > 0 else 'min'}"
(CAMP / "review" / "cx_runs" / f"{tag}.json").write_text(json.dumps(out, indent=1))
np.save(CAMP / "review" / "cx_runs" / f"{tag}_u.npy", u)
