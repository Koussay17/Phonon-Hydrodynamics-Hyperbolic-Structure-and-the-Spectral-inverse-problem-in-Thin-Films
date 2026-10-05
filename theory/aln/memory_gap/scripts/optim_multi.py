"""Multi-constraint version of optim.Slice: several DC components fixed simultaneously.

Constraints: diag C(g) = r (representative rows), and K_cd(g) = b_c^T C^+ b_d = K0_cd for every
pair (c, d) in `kpairs`. Objective: N_o(g) = |C^+ b_o|^2 for the column o = `col`.
Same algorithm as optim.Slice (multiplicative projected gradient + exact Newton restoration).
"""
from __future__ import annotations

import math

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import linprog

from fastops import FastOps

_ERR = (np.linalg.LinAlgError, sla.LinAlgError, ValueError, FloatingPointError)


class MultiSlice:
    def __init__(self, geom, r, K0, kpairs, col=0, tie=None):
        """K0: dict {(c, d): value}; kpairs: list of (c, d)."""
        self.geom, self.col = geom, col
        self.ops = FastOps(geom, tie)
        self.r_full = np.asarray(r, float)
        self.r = self.r_full[self.ops.rows]
        self.rs = np.maximum(self.r, 1e-300)
        self.nr = self.ops.n_rows
        self.mvar = self.ops.m_var
        self.kpairs = list(kpairs)
        self.K0 = {k: float(K0[k]) for k in self.kpairs}
        self.cols = sorted(set([c for p in self.kpairs for c in p] + [col]))
        self.kscale = {(c, d): math.sqrt(abs(self.K0[(c, c)] * self.K0[(d, d)])) if (c, c) in self.K0 and (d, d) in self.K0
                       else abs(self.K0[(c, d)]) for (c, d) in self.kpairs}
        self._Wsp = None

    def g_of(self, gam):
        return self.ops.g_of(gam)

    def gamma_of(self, g):
        return self.ops.gamma_of(g)

    # ------------------------------------------------------------------ evaluation
    def evaluate(self, gam, grad=True):
        ops = self.ops
        g = ops.g_of(gam)
        C = ops.C(g)
        greg = max(np.trace(C) / ops.n, 1e-300)
        cf = sla.cho_factor(C + greg * np.outer(ops.ehat, ops.ehat), lower=True, check_finite=False)
        B = self.geom.b[:, self.cols]
        X = sla.cho_solve(cf, B, check_finite=False)
        ci = {c: i for i, c in enumerate(self.cols)}
        R = {"X": X, "ci": ci}
        R["K"] = {(c, d): float(B[:, ci[c]] @ X[:, ci[d]]) for (c, d) in self.kpairs}
        o = ci[self.col]
        xo = X[:, o]
        R["Kobj"] = float(B[:, o] @ xo)
        R["N"] = float(xo @ xo)
        R["tau"] = R["N"] / R["Kobj"]
        if grad:
            P = np.column_stack([ops.ATx(X[:, i]) for i in range(len(self.cols))])  # (m, ncols)
            u = sla.cho_solve(cf, xo, check_finite=False)
            q = ops.ATx(u)
            R["dN"] = ops.P_dot(-2.0 * P[:, o] * q)
            R["dK"] = {(c, d): ops.P_dot(-P[:, ci[c]] * P[:, ci[d]]) for (c, d) in self.kpairs}
            R["h"] = {(c, d): ops.P_dot(P[:, ci[c]] * P[:, ci[d]]) for (c, d) in self.kpairs}
        return R

    def kres(self, R):
        return np.array([(R["K"][k] - self.K0[k]) / self.kscale[k] for k in self.kpairs])

    def feas(self, gam, R=None):
        if R is None:
            R = self.evaluate(gam, grad=False)
        g = self.g_of(gam)
        rf = np.maximum(self.r_full, 1e-300)
        return max(float(np.max(np.abs(self.ops.diag_full(g) - self.r_full) / rf)),
                   float(np.max(np.abs(self.kres(R)))))

    def safe_eval(self, gam):
        try:
            return self.evaluate(gam, grad=False)
        except _ERR:
            return None

    # ------------------------------------------------------------------ diag projection
    def project_diag(self, g0, tol=1e-13, maxit=200):
        ops, r, rs = self.ops, self.r, self.rs
        n = self.nr
        lam = np.zeros(n)

        def gl(lam):
            return g0 * np.exp(np.clip(ops.WTv(lam), -600, 600))

        g = gl(lam)
        f = float(np.sum(g))
        for it in range(maxit):
            F = ops.Wv(g) - r
            if np.max(np.abs(F) / rs) < tol:
                return g, True
            H = ops.WGW(g)
            H[np.diag_indices(n)] += 1e-15 * np.max(np.diag(H))
            try:
                step = -sla.solve(H, F, assume_a="pos")
            except _ERR:
                step = -np.linalg.lstsq(H, F, rcond=None)[0]
            t = 1.0
            dec = F @ step
            while t > 1e-12:
                lt = lam + t * step
                gt = gl(lt)
                ft = float(np.sum(gt) - r @ lt)
                if np.isfinite(ft) and ft <= f + 1e-4 * t * dec:
                    break
                t *= 0.5
            lam, g, f = lt, gt, ft
        return g, bool(np.max(np.abs(ops.Wv(g) - r) / rs) < 1e3 * tol)

    # ------------------------------------------------------------------ restoration
    def restore(self, g0, tol=1e-12, maxit=60):
        ops, r, rs = self.ops, self.r, self.rs
        n, nk = self.nr, len(self.kpairs)
        try:
            R0 = self.evaluate(g0)
        except _ERR:
            return g0, False, None
        H = np.column_stack([R0["h"][k] / max(np.max(np.abs(R0["h"][k])), 1e-300) for k in self.kpairs])
        lam = np.zeros(n)
        eta = np.zeros(nk)

        def make(lam, eta):
            return g0 * np.exp(np.clip(ops.WTv(lam) + H @ eta, -600, 600))

        def resid(g):
            R = self.evaluate(g, grad=True)
            return np.concatenate([(ops.Wv(g) - r) / rs, self.kres(R)]), R

        g = make(lam, eta)
        try:
            F, R = resid(g)
        except _ERR:
            return g0, False, None
        for it in range(maxit):
            if np.max(np.abs(F)) < tol:
                return g, True, R
            dK = np.column_stack([R["dK"][k] / self.kscale[k] for k in self.kpairs])   # (mvar, nk)
            J = np.empty((n + nk, n + nk))
            J[:n, :n] = ops.WGW(g) / rs[:, None]
            for j in range(nk):
                J[:n, n + j] = ops.Wv(g * H[:, j]) / rs
                J[n + j, :n] = ops.Wv(g * dK[:, j])
            J[n:, n:] = (dK * g[:, None]).T @ H
            try:
                step = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -F, rcond=None)[0]
            t = 1.0
            nF = np.linalg.norm(F)
            accepted = False
            while t > 1e-6:
                lam_t = lam + t * step[:n]
                eta_t = eta + t * step[n:]
                g_t = make(lam_t, eta_t)
                try:
                    F_t, R_t = resid(g_t)
                except _ERR:
                    t *= 0.5
                    continue
                if np.all(np.isfinite(F_t)) and np.linalg.norm(F_t) < (1 - 1e-4 * t) * nF:
                    accepted = True
                    break
                t *= 0.5
            if not accepted:
                return g, False, R
            lam, eta, g, F, R = lam_t, eta_t, g_t, F_t, R_t
        return g, bool(np.max(np.abs(F)) < tol * 10), R

    def restore_continuation(self, g0, nsteps_max=40):
        g1, ok = self.project_diag(g0)
        if not ok:
            return None
        R1 = self.safe_eval(g1)
        if R1 is None:
            return None
        K1 = dict(R1["K"])
        K0 = dict(self.K0)
        s, ds = 0.0, 1.0
        g = g1
        for _ in range(nsteps_max):
            s_try = min(1.0, s + ds)
            self.K0 = {k: (1 - s_try) * K1[k] + s_try * K0[k] for k in self.kpairs}
            try:
                g_new, ok, R = self.restore(g)
            finally:
                self.K0 = K0
            if ok:
                g, s = g_new, s_try
                if s >= 1.0:
                    return g
                ds = min(1.0, 2 * ds)
            else:
                ds *= 0.5
                if ds < 1e-4:
                    return None
        return None

    # ------------------------------------------------------------------ starts
    def W_sparse(self):
        if self._Wsp is None:
            ops = self.ops
            if ops.Wg_dense is not None:
                self._Wsp = sp.csr_matrix(ops.Wg_dense)
            else:
                rows = ops.rp.ravel()
                cols = np.repeat(np.arange(ops.m), 3)
                vals = ops.vals2.ravel()
                keep = rows < ops.n_rows
                self._Wsp = sp.coo_matrix((vals[keep], (rows[keep], cols[keep])), shape=(ops.n_rows, ops.m)).tocsr()
        return self._Wsp

    def lp_vertex(self, rng, cost=None, sparse_cost=True):
        m = self.mvar
        if cost is None:
            cost = rng.exponential(size=m) if sparse_cost else rng.normal(size=m)
        res = linprog(cost, A_eq=self.W_sparse(), b_eq=self.r, bounds=(0, None), method="highs")
        if res.status != 0:
            return None
        return res.x

    def mix_start(self, v, g_ref):
        """Feasible start on the path from vertex v toward g_ref: restore from (1-t) v + t g_ref."""
        for t in (0.02, 0.05, 0.1, 0.2, 0.4, 0.7):
            g = (1 - t) * v + t * g_ref
            g2, ok, R = self.restore(np.maximum(g, 1e-300))
            if ok:
                return g2
            g3 = self.restore_continuation(np.maximum(g, 1e-300))
            if g3 is not None:
                return g3
        return None

    # ------------------------------------------------------------------ ascent
    def tangent_direction(self, g, R, sign):
        ops = self.ops
        n, nk = self.nr, len(self.kpairs)
        grad = sign * R["dN"] / R["N"]
        dK = np.column_stack([R["dK"][k] / self.kscale[k] for k in self.kpairs])
        J = np.empty((n + nk, n + nk))
        J[:n, :n] = ops.WGW(g)
        for j in range(nk):
            J[:n, n + j] = ops.Wv(g * dK[:, j])
            J[n + j, :n] = J[:n, n + j]
        J[n:, n:] = (dK * g[:, None]).T @ dK
        rhs = np.concatenate([ops.Wv(g * grad), (dK * g[:, None]).T @ grad])
        J[np.diag_indices(n + nk)] += 1e-14 * np.max(np.abs(np.diag(J)))
        try:
            sol = sla.solve(J, rhs, assume_a="sym")
        except _ERR:
            sol = np.linalg.lstsq(J, rhs, rcond=None)[0]
        return grad - ops.WTv(sol[:n]) - dK @ sol[n:]

    def run(self, g_init, sign, maxit=1500, t0=0.5, already_feasible=False, window=30, rel_stall=1e-6,
            verbose=False):
        if already_feasible:
            g = g_init
            R = self.evaluate(g)
        else:
            g, ok, R = self.restore(g_init)
            if not ok:
                g = self.restore_continuation(g_init)
                if g is None:
                    return None
                R = self.evaluate(g)
        t = t0
        val = math.log(R["N"])
        vals = [val]
        it = 0
        for it in range(maxit):
            if len(vals) > window and sign * (vals[-1] - vals[-1 - window]) < rel_stall:
                break
            red = self.tangent_direction(g, R, sign)
            scale = np.max(np.abs(red))
            if not np.isfinite(scale) or scale == 0:
                break
            improved = False
            while t > 1e-12:
                g_try = g * np.exp(np.clip(t * red / scale, -60, 60))
                g_new, ok, R_new = self.restore(g_try)
                if ok:
                    v_new = math.log(R_new["N"])
                    if sign * (v_new - val) > 0:
                        improved = True
                        break
                t *= 0.5
            if not improved:
                break
            g, R, val = g_new, R_new, v_new
            vals.append(val)
            t = min(2 * t, 60.0)
            if verbose and it % 25 == 0:
                print(it, R["tau"], t, flush=True)
        return {"g": g, "tau": R["tau"], "N": R["N"], "K": R["K"], "feas": self.feas(g, R), "iters": it}
