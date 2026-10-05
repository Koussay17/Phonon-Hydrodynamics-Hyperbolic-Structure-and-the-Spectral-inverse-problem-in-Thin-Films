"""General constrained optimiser for the review attacks (independent of the campaign's optim*.py).

Variables: u (one per event orbit, or per event), bounds lo <= u <= hi; rates at temperature block k:
    g_k = g0_k * exp(u[orb])          (g0_k = w0 * Bf(T_k); w0 temperature independent)
Constraints (all scaled to relative residuals):
    lifetimes   (W_k g_k)[rows] / r_k[rows] - 1 = 0          for each block k
    DC response K_k[c] / K0_k[c] - 1 = 0                      for each block k and column c in kcols
    optional extras on block 0:
      ('Kz', c, z)     : K(z) = b^T (z + C)^{-1} b fixed      (z in THz, model units)
      ('S2', c)        : second moment b^T C b fixed           (linear in g)
      ('acc', c, fcut) : accumulated DC response sum_{f_mu < fcut} b_mu x_mu fixed
Objective: sign * log tau_c(block 0).
Method: gradient projection on the free variables with the linearised constraints, line search,
exact Gauss-Newton restoration (min-norm on free variables), bound clipping.
"""
from __future__ import annotations

import math
import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp

ERR = (np.linalg.LinAlgError, sla.LinAlgError, ValueError, FloatingPointError)


def gsolve(G, rhs, rcond=None):
    """Solve G y = rhs for symmetric PSD G; with rcond, use the truncated pseudo-inverse."""
    if rcond is None:
        G = G.copy()
        G[np.diag_indices_from(G)] += 1e-14 * max(np.max(np.diag(G)), 1e-300)
        try:
            return sla.solve(G, rhs, assume_a="pos")
        except ERR:
            return np.linalg.lstsq(G, rhs, rcond=None)[0]
    w, V = np.linalg.eigh(G)
    keep = w > rcond * max(w.max(), 1e-300)
    return V[:, keep] @ ((V[:, keep].T @ rhs) / w[keep])


class Block:
    def __init__(self, M, g0, r_target, K0, kcols=(0, 2), fix_lifetimes=True):
        self.M = M
        self.fix_lifetimes = bool(fix_lifetimes)
        self.g0 = np.asarray(g0, float)
        self.kcols = tuple(kcols)
        self.K0 = {c: float(K0[c]) for c in self.kcols}
        self.r = np.asarray(r_target, float)


class Problem:
    def __init__(self, blocks, orb, n_var, rows, lo, hi, obj_col=0, sign=+1, extras=(), symmetric=True):
        self.blocks = blocks
        self.orb = np.asarray(orb)
        self.p = int(n_var)
        m = len(self.orb)
        self.Porb = sp.csr_matrix((np.ones(m), (self.orb, np.arange(m))), shape=(self.p, m))
        self.rows = np.asarray(rows)
        self.lo = np.broadcast_to(np.asarray(lo, float), (self.p,)).copy()
        self.hi = np.broadcast_to(np.asarray(hi, float), (self.p,)).copy()
        self.obj_col = obj_col
        self.sign = sign
        self.extras = list(extras)
        self.Wrows = [blk.M.W[self.rows, :].tocsr() for blk in blocks]
        self.nfev = 0
        # extras targets are set by set_extra_targets(u_ref)
        self.extra_t = None

    # ---------------------------------------------------------------- evaluation
    def rates(self, u, k):
        return self.blocks[k].g0 * np.exp(u[self.orb])

    def _block_eval(self, k, g, need_grad, cols_needed):
        blk = self.blocks[k]
        M = blk.M
        C = M.C(g)
        gam = np.trace(C) / M.n
        cf = sla.cho_factor(C + gam * np.outer(M.ehat, M.ehat), lower=True, check_finite=False)
        cols = sorted(set(cols_needed))
        X = sla.cho_solve(cf, M.b[:, cols], check_finite=False)
        ci = {c: i for i, c in enumerate(cols)}
        out = {"C": C, "cf": cf, "X": X, "ci": ci}
        out["K"] = {c: float(M.b[:, c] @ X[:, ci[c]]) for c in cols}
        if need_grad:
            out["AX"] = {c: M.AT @ X[:, ci[c]] for c in cols}
        return out

    def evaluate(self, u, need_grad=True):
        self.nfev += 1
        res_list, jac_list = [], []
        info = {}
        for k, blk in enumerate(self.blocks):
            g = self.rates(u, k)
            cols_needed = list(blk.kcols) + ([self.obj_col] if k == 0 else [])
            E = self._block_eval(k, g, need_grad, cols_needed)
            # lifetimes
            if blk.fix_lifetimes:
                d = self.Wrows[k] @ g
                res_list.append(d / blk.r - 1.0)
                if need_grad:
                    Jl = (self.Wrows[k].multiply(g[None, :]) @ self.Porb.T).toarray() / blk.r[:, None]
                    jac_list.append(Jl)
            for c in blk.kcols:
                res_list.append(np.array([E["K"][c] / blk.K0[c] - 1.0]))
                if need_grad:
                    dK = -(E["AX"][c] ** 2) * g / blk.K0[c]
                    jac_list.append((self.Porb @ dK)[None, :])
            if k == 0:
                c = self.obj_col
                M = blk.M
                x = E["X"][:, E["ci"][c]]
                K = E["K"][c]
                N = float(x @ x)
                info.update(K=K, N=N, tau=N / K, tau_ps=N / K / (4 * np.pi), E0=E, g0=g)
                if need_grad:
                    y = sla.cho_solve(E["cf"], x, check_finite=False)
                    p_ = E["AX"][c]
                    q_ = M.AT @ y
                    dN = -2.0 * p_ * q_
                    dK = -(p_ ** 2)
                    info["grad"] = self.sign * (self.Porb @ ((dN / N - dK / K) * g))
                # extras on block 0
                for j, ex in enumerate(self.extras):
                    val, grad = self._extra(ex, E, g, need_grad)
                    t = self.extra_t[j] if self.extra_t is not None else val
                    res_list.append(np.array([val / t - 1.0]))
                    if need_grad:
                        jac_list.append((self.Porb @ (grad * g / t))[None, :])
        info["c"] = np.concatenate(res_list)
        if need_grad:
            info["J"] = np.vstack(jac_list)
        info["f"] = self.sign * math.log(info["tau"])
        return info

    def _extra(self, ex, E, g, need_grad):
        M = self.blocks[0].M
        kind = ex[0]
        if kind == "Kz":
            _, c, z = ex
            A = E["C"] + z * np.eye(M.n)
            cf = sla.cho_factor(A, lower=True, check_finite=False)
            xz = sla.cho_solve(cf, M.b[:, c], check_finite=False)
            val = float(M.b[:, c] @ xz)
            grad = -((M.AT @ xz) ** 2) if need_grad else None
            return val, grad
        if kind == "S2":
            _, c = ex
            ab = M.AT @ M.b[:, c]
            return float(g @ ab ** 2), (ab ** 2 if need_grad else None)
        if kind == "acc":
            _, c, fcut = ex
            mask = (M.eps < fcut).astype(float)
            pb = mask * M.b[:, c]
            # accumulated response = (Pi b)^T C^+ b  (symmetric in the two vectors)
            v = sla.cho_solve(E["cf"], pb - M.ehat * (M.ehat @ pb), check_finite=False)
            x = sla.cho_solve(E["cf"], M.b[:, c], check_finite=False)
            val = float(pb @ x)
            grad = -((M.AT @ v) * (M.AT @ x)) if need_grad else None
            return val, grad
        raise ValueError(kind)

    def set_extra_targets(self, u_ref):
        self.extra_t = None
        I = self.evaluate(u_ref, need_grad=False)
        # extras residuals were computed relative to themselves (=0); recompute values
        vals = []
        E = I["E0"]
        g = I["g0"]
        for ex in self.extras:
            vals.append(self._extra(ex, E, g, False)[0])
        self.extra_t = vals
        return vals

    # ---------------------------------------------------------------- restoration
    def restore(self, u, free=None, tol=None, maxit=40, I=None):
        """Gauss-Newton (min-norm on variables not at a bound) onto c(u) = 0. Reuses evaluations."""
        if tol is None:
            tol = getattr(self, "tol", 1e-11)
        u = np.clip(u, self.lo, self.hi)
        if I is None:
            try:
                I = self.evaluate(u, need_grad=True)
            except ERR:
                return u, False, None
        for it in range(maxit):
            c = I["c"]
            nc = float(np.max(np.abs(c)))
            if not np.isfinite(nc):
                return u, False, None
            if nc < tol:
                return u, True, I
            fr = ~((u <= self.lo + 1e-12) | (u >= self.hi - 1e-12))
            JF = I["J"][:, fr]
            G = JF @ JF.T
            lam = gsolve(G, c, getattr(self, "rcond", None))
            du = np.zeros(self.p)
            du[fr] = -JF.T @ lam
            t = 1.0
            ok = False
            while t > 1e-4:
                ut = np.clip(u + t * du, self.lo, self.hi)
                try:
                    It = self.evaluate(ut, need_grad=True)
                except ERR:
                    t *= 0.5
                    continue
                if np.all(np.isfinite(It["c"])) and np.max(np.abs(It["c"])) < (1 - 1e-4 * t) * nc:
                    ok = True
                    break
                t *= 0.5
            if not ok:
                return u, False, I
            u, I = ut, It
        return u, bool(np.max(np.abs(I["c"])) < 10 * tol), I


    # ---------------------------------------------------------------- homotopy start
    def homotopy(self, u_target, ds0=0.05, ds_min=1e-4, time_limit=None):
        """Walk from u = 0 (feasible reference) toward u_target, restoring at each step.
        Returns the last feasible point reached and the fraction s of the path covered."""
        t_start = time.time()
        u, ok, I = self.restore(np.zeros(self.p))
        if not ok:
            return None, 0.0
        s, ds = 0.0, ds0
        while s < 1.0 and ds >= ds_min:
            if time_limit is not None and time.time() - t_start > time_limit:
                break
            s_try = min(1.0, s + ds)
            ut = np.clip(u + (s_try - s) * (u_target - u) / max(1.0 - s, 1e-12), self.lo, self.hi)
            un, okr, In = self.restore(ut)
            if okr:
                u, I, s = un, In, s_try
                ds = min(2 * ds, 0.25)
            else:
                ds *= 0.5
        return u, s

    # ---------------------------------------------------------------- ascent
    def direction(self, u, I):
        grad = I["grad"]
        J = I["J"]
        at_lo = u <= self.lo + 1e-12
        at_hi = u >= self.hi - 1e-12
        free = ~(at_lo | at_hi)
        for _ in range(50):
            JF = J[:, free]
            if getattr(self, "rcond", None) is not None:
                # exact tangent projection (no truncation) when constraints are nearly dependent
                lam = np.linalg.lstsq(JF.T, grad[free], rcond=None)[0]
            else:
                G = JF @ JF.T
                lam = gsolve(G, JF @ grad[free], None)
            red = grad - J.T @ lam
            scale = max(np.max(np.abs(red)), 1e-300)
            release = (at_lo & ~free & (red > 1e-10 * scale)) | (at_hi & ~free & (red < -1e-10 * scale))
            if not release.any():
                break
            free = free | release
        d = np.zeros(self.p)
        d[free] = red[free]
        bad = (at_lo & (d < 0)) | (at_hi & (d > 0))
        d[bad] = 0.0
        return d, free

    def run(self, u0, maxit=500, step0=0.5, window=25, rel_stall=1e-7, verbose=False, time_limit=None,
            max_du=3.0):
        """Projected-gradient ascent with Barzilai-Borwein trial steps (max |du| capped by max_du)."""
        t_start = time.time()
        u, ok, I = self.restore(u0)
        if not ok:
            return None
        hist = [I["f"]]
        d_prev = u_prev = None
        tmax = step0
        it = 0
        for it in range(maxit):
            if time_limit is not None and time.time() - t_start > time_limit:
                break
            if len(hist) > window and (hist[-1] - hist[-1 - window]) < rel_stall:
                break
            d, free = self.direction(u, I)
            sc = np.max(np.abs(d))
            if not np.isfinite(sc) or sc == 0:
                break
            # trial step length (in units of max |du|)
            if d_prev is not None:
                s = u - u_prev
                y = d - d_prev
                sy = float(s @ y)
                t_bb = abs(float(s @ s) / sy) if sy != 0 else tmax
                t = min(max(t_bb * sc, 1e-6), max_du)
            else:
                t = min(tmax, max_du)
            improved = False
            while t > 1e-10:
                ut = np.clip(u + t * d / sc, self.lo, self.hi)
                un, okr, In = self.restore(ut)
                if okr and In["f"] > I["f"]:
                    improved = True
                    break
                t *= 0.3
            if not improved:
                break
            u_prev, d_prev = u, d
            u, I = un, In
            tmax = t
            hist.append(I["f"])
            if verbose and it % 10 == 0:
                print(f"  it {it} tau_ps {I['tau_ps']:.6g} step {t:.3g} nfev {self.nfev}", flush=True)
        return {"u": u, "tau_ps": I["tau_ps"], "f": I["f"], "c_max": float(np.max(np.abs(I["c"]))),
                "iters": it, "nfev": self.nfev, "hist": hist}
