"""Attack 1e (cross-check): smooth-amplitude prior optimised by projected gradient in coefficient space.

Same class as cx_attack1_smooth.py (g = g_ref * exp(P), P a polynomial of degree <= k in the three
scaled frequencies, symmetric in a <-> b), but optimised with the review's gradient-projection method:
variables c, u = B c; constraints c(u) = 0 (117 lifetimes + K_xx + K_zz); optional box |u| <= L handled
by halving the step if violated. Reports tau along the path and the final witness.
usage: python cx_attack1_smooth_pg.py k col [maxit] [L]
Output: review/cx_runs/smoothpg_k{k}_c{col}.json
"""
import json
import sys
import time

import numpy as np
import scipy.linalg as sla

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP  # noqa: E402
from cx_opt import Block, Problem  # noqa: E402

k = int(sys.argv[1]); col = int(sys.argv[2])
maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 300
Lbox = float(sys.argv[4]) if len(sys.argv) > 4 else 30.0
M = AlN(); o = M.orbits(); gref = M.g_ref(); r = M.diag(gref)
mo = M.moments(gref, (0, 2)); K0 = {c: mo[c]["K"] for c in (0, 2)}
wmax = M.eps.max()
rep_ev = np.array([np.nonzero(o["ev_orb"] == q)[0][0] for q in range(o["n_ev_orb"])])
xp = 2 * M.eps[M.P[rep_ev]] / wmax - 1
xa = 2 * M.eps[M.A[rep_ev]] / wmax - 1
xb = 2 * M.eps[M.B[rep_ev]] / wmax - 1
cols = []
for i in range(k + 1):
    for j in range(k + 1 - i):
        for l in range(j, k + 1 - i - j):
            cols.append(xp ** i * (xa ** j * xb ** l + xa ** l * xb ** j))
Bm = np.array(cols).T
Qb, Rb = np.linalg.qr(Bm)
rank = int(np.sum(np.abs(np.diag(Rb)) > 1e-10 * np.abs(Rb[0, 0])))
B = Qb[:, :rank] * np.sqrt(len(xp))
blk = Block(M, gref, r[o["mode_reps"]], K0, kcols=(0, 2))
P = Problem([blk], o["ev_orb"], o["n_ev_orb"], o["mode_reps"], -np.inf, np.inf, obj_col=col, sign=1)


def ev(c):
    I = P.evaluate(B @ c)
    I["Jc"] = I["J"] @ B
    I["gc"] = B.T @ I["grad"]
    return I


def restore(c, I=None, tol=1e-11):
    if I is None:
        I = ev(c)
    for _ in range(30):
        nc = np.max(np.abs(I["c"]))
        if nc < tol:
            return c, I, True
        Jc = I["Jc"]
        dc = -np.linalg.lstsq(Jc, I["c"], rcond=None)[0]
        t = 1.0
        while t > 1e-4:
            ct = c + t * dc
            if np.max(np.abs(B @ ct)) <= Lbox:
                try:
                    It = ev(ct)
                    if np.max(np.abs(It["c"])) < (1 - 1e-4 * t) * nc:
                        break
                except Exception:
                    pass
            t *= 0.5
        else:
            return c, I, False
        c, I = ct, It
    return c, I, np.max(np.abs(I["c"])) < 10 * tol


t0 = time.time()
c = np.zeros(rank)
I = ev(c)
hist = [I["tau_ps"]]
step = 0.05
tolb = 1e-9


def restore_box(c, I=None, tol=1e-11):
    """Gauss-Newton on c(u)=0 keeping the active box rows fixed (min-norm in the free subspace)."""
    if I is None:
        I = ev(c)
    for _ in range(30):
        nc = np.max(np.abs(I["c"]))
        if nc < tol:
            return c, I, True
        u = B @ c
        act = np.abs(np.abs(u) - Lbox) < 1e-7
        Aact = B[act]
        # null space of the active rows
        if act.any():
            Z = sla.null_space(Aact)
        else:
            Z = np.eye(rank)
        Jz = I["Jc"] @ Z
        dz = -np.linalg.lstsq(Jz, I["c"], rcond=None)[0]
        dc = Z @ dz
        t = 1.0
        while t > 1e-4:
            ct = c + t * dc
            ut = B @ ct
            if np.max(np.abs(ut)) <= Lbox + 1e-6:
                try:
                    It = ev(ct)
                    if np.max(np.abs(It["c"])) < (1 - 1e-4 * t) * nc:
                        break
                except Exception:
                    pass
            t *= 0.5
        else:
            return c, I, False
        c, I = ct, It
    return c, I, np.max(np.abs(I["c"])) < 10 * tol


for it in range(maxit):
    u = B @ c
    sgn = np.sign(u)
    act = np.abs(np.abs(u) - Lbox) < 1e-7
    Jc, gc = I["Jc"], I["gc"]
    released = True
    while released:
        Aact = B[act] * sgn[act][:, None]           # outward normals of active rows
        G = np.vstack([Jc, Aact])
        lam = np.linalg.lstsq(G.T, gc, rcond=None)[0]
        d = gc - G.T @ lam
        mu = lam[Jc.shape[0]:]
        released = False
        if act.any() and np.any(mu < -1e-12 * max(np.abs(mu).max(), 1e-300)):
            idx = np.nonzero(act)[0]
            act[idx[np.argmin(mu)]] = False           # release the most negative multiplier
            released = True
    nd = np.linalg.norm(d)
    if nd < 1e-12:
        break
    ok = False
    while step > 1e-8:
        ct = c + step * d / nd
        ut = B @ ct
        over = np.abs(ut) > Lbox
        if over.any():
            # scale back to the first boundary hit
            uc = B @ c
            du = ut - uc
            with np.errstate(divide="ignore", invalid="ignore"):
                tt = np.where(du > 0, (Lbox - uc) / du, np.where(du < 0, (-Lbox - uc) / du, np.inf))
            frac = float(np.clip(np.min(tt[over | (tt > 0)]), 0, 1))
            ct = c + frac * step * d / nd
        try:
            cn, In, okr = restore_box(ct)
        except Exception:
            okr = False
        if okr and In["f"] > I["f"] + 1e-12:
            ok = True
            break
        step *= 0.5
    if not ok:
        break
    c, I = cn, In
    hist.append(I["tau_ps"])
    step = min(step * 1.5, 5.0)
    if it % 20 == 0:
        print(it, I["tau_ps"], step, float(np.max(np.abs(B @ c))), int(np.sum(np.abs(np.abs(B @ c) - Lbox) < 1e-7)), flush=True)
u = B @ c
mo2 = M.moments(P.rates(u, 0), (0, 2))
out = {"degree": k, "n_coef": rank, "col": col, "Lbox": Lbox, "iters": it, "tau_ps": I["tau_ps"],
       "tau_x_ps": mo2[0]["tau_ps"], "tau_z_ps": mo2[2]["tau_ps"], "c_max": float(np.max(np.abs(I["c"]))),
       "u_min": float(u.min()), "u_max": float(u.max()), "seconds": time.time() - t0, "hist_tail": hist[-5:],
       "hist_every20": hist[::20]}
print(json.dumps(out), flush=True)
(CAMP / "review" / "cx_runs" / f"smoothpg_k{k}_c{col}_L{int(Lbox)}.json").write_text(json.dumps(out, indent=1))
np.save(CAMP / "review" / "cx_runs" / f"smoothpg_k{k}_c{col}_L{int(Lbox)}_u.npy", u)
