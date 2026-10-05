"""Extremization of tau_mem on F = {gamma >= 0 : W g = r, K(g) = K0}, g = P^T gamma.

Method: Riemannian (multiplicative-metric) projected gradient, exponentiated step, and exact
restoration by a Bregman (KL-type) projection onto F (Newton in the n_rows+1 dual variables).
Every accepted iterate is feasible to the restoration tolerance (relative 1e-12), so every
reported value is attained by an explicit witness (inner approximation of [tau_min, tau_max]).

Start generation: (i) lognormal perturbations of the reference; (ii) points on segments between
random LP vertices of P = {gamma >= 0 : W g = r} and the K-minimiser over P, at the crossing K = K0.
Numerical core: fastops.FastOps (bincount assembly, no sparse temporaries).
"""
from __future__ import annotations

import math

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from scipy.optimize import linprog

from fastops import FastOps

_LINALG_ERR = (np.linalg.LinAlgError, sla.LinAlgError, ValueError, FloatingPointError)


class Slice:
    def __init__(self, geom, r, K0, col=0, tie=None):
        self.geom, self.col = geom, col
        self.ops = FastOps(geom, tie)
        self.r_full = np.asarray(r, float)
        self.K0 = float(K0)
        self.r = self.r_full[self.ops.rows]
        self.rs = np.maximum(self.r, 1e-300)
        self.nr = self.ops.n_rows
        self.mvar = self.ops.m_var
        self._Wsp = None

    def g_of(self, gam):
        return self.ops.g_of(gam)

    def gamma_of(self, g):
        return self.ops.gamma_of(g)

    # ---------------------------------------------------------------- evaluation
    def evaluate(self, gam, grad=True):
        return self.ops.response(gam, self.col, want_grad=grad)

    def feas(self, gam, R=None):
        if R is None:
            R = self.evaluate(gam, grad=False)
        g = self.g_of(gam)
        rf = np.maximum(self.r_full, 1e-300)
        return max(float(np.max(np.abs(self.ops.diag_full(g) - self.r_full) / rf)), abs(R["K"] / self.K0 - 1.0))

    def safe_K(self, gam):
        try:
            return self.evaluate(gam, grad=False)["K"]
        except _LINALG_ERR:
            return np.inf

    # ---------------------------------------------------------------- diagonal projection
    def project_diag(self, g0, tol=1e-13, maxit=200):
        """KL projection of g0 (>0 on its support) onto {W g = r}: dual Newton (convex)."""
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
            except _LINALG_ERR:
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

    # ---------------------------------------------------------------- full restoration
    def restore(self, g0, tol=1e-12, maxit=60, h=None):
        """Projection of g0 onto F via g = g0 exp(W^T lam + eta h). Returns (g, ok, R)."""
        ops, r, rs, K0 = self.ops, self.r, self.rs, self.K0
        n = self.nr
        try:
            if h is None:
                R0 = self.evaluate(g0)
                h = R0["p2"]
                h = h / max(h.max(), 1e-300)
        except _LINALG_ERR:
            return g0, False, None
        lam = np.zeros(n)
        eta = 0.0

        def make(lam, eta):
            return g0 * np.exp(np.clip(ops.WTv(lam) + eta * h, -600, 600))

        def resid(g):
            R = self.evaluate(g, grad=True)
            F = np.append((ops.Wv(g) - r) / rs, math.log(R["K"] / K0))
            return F, R

        g = make(lam, eta)
        try:
            F, R = resid(g)
        except _LINALG_ERR:
            return g0, False, None
        for it in range(maxit):
            if np.max(np.abs(F)) < tol:
                return g, True, R
            dlogK = R["dK"] / R["K"]
            J = np.empty((n + 1, n + 1))
            J[:n, :n] = ops.WGW(g) / rs[:, None]
            J[:n, n] = ops.Wv(g * h) / rs
            J[n, :n] = ops.Wv(g * dlogK)
            J[n, n] = float(np.sum(g * h * dlogK))
            try:
                step = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(J, -F, rcond=None)[0]
            t = 1.0
            nF = np.linalg.norm(F)
            accepted = False
            while t > 1e-6:
                lam_t = lam + t * step[:n]
                eta_t = eta + t * step[n]
                g_t = make(lam_t, eta_t)
                try:
                    F_t, R_t = resid(g_t)
                except _LINALG_ERR:
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
        """Project onto W g = r, then move K to K0 by continuation of the K target."""
        g1, ok = self.project_diag(g0)
        if not ok:
            return None
        K1 = self.safe_K(g1)
        if not np.isfinite(K1):
            return None
        K0 = self.K0
        s, ds = 0.0, 1.0
        g = g1
        for _ in range(nsteps_max):
            s_try = min(1.0, s + ds)
            Kt = math.exp((1 - s_try) * math.log(K1) + s_try * math.log(K0))
            self.K0 = Kt
            try:
                g_new, ok, R = self.restore(g)
            finally:
                self.K0 = K0
            if ok:
                g, s = g_new, s_try
                if s >= 1.0:
                    return g
                ds = min(1.0, ds * 2)
            else:
                ds *= 0.5
                if ds < 1e-4:
                    return None
        return None

    # ---------------------------------------------------------------- start generation
    def k_minimizer(self, g_init, maxit=300):
        """Minimiser of K over P (convex problem) by mirror descent + diagonal projection."""
        g, ok = self.project_diag(g_init)
        K = self.safe_K(g)
        t = 0.5
        for it in range(maxit):
            R = self.evaluate(g)
            grad = R["dK"] / R["K"]
            H = self.ops.WGW(g)
            H[np.diag_indices(self.nr)] += 1e-14 * np.max(np.diag(H))
            lam = sla.solve(H, self.ops.Wv(g * grad), assume_a="pos")
            red = grad - self.ops.WTv(lam)
            sc = np.max(np.abs(red))
            if sc < 1e-12:
                break
            improved = False
            while t > 1e-10:
                gt, ok = self.project_diag(g * np.exp(-t * red / sc))
                Kt = self.safe_K(gt) if ok else np.inf
                if Kt < K:
                    improved = True
                    break
                t *= 0.5
            if not improved:
                break
            done = (K - Kt) / K < 1e-9
            g, K = gt, Kt
            if done:
                break
            t = min(2 * t, 50.0)
        return g, K

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
                self._Wsp = sp.coo_matrix((vals[keep], (rows[keep], cols[keep])),
                                          shape=(ops.n_rows, ops.m)).tocsr()
        return self._Wsp

    def lp_vertex(self, rng, sparse_cost=True):
        m = self.mvar
        c = rng.exponential(size=m) if sparse_cost else rng.normal(size=m)
        res = linprog(c, A_eq=self.W_sparse(), b_eq=self.r, bounds=(0, None), method="highs")
        if res.status != 0:
            return None
        return res.x

    def segment_start(self, v, g_lo, Kv=None):
        """Point on [v, g_lo] with K = K0 (K(g_lo) < K0 < K(v)); bisection in t."""
        if Kv is None:
            Kv = self.safe_K(v)
        Klo = self.safe_K(g_lo)
        if not (Klo < self.K0) or not (Kv > self.K0):
            return None
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            Km = self.safe_K((1 - mid) * v + mid * g_lo)
            if Km > self.K0:
                lo = mid
            else:
                hi = mid
            if hi - lo < 1e-13:
                break
        g = (1 - hi) * v + hi * g_lo
        g2, ok, R = self.restore(np.maximum(g, 1e-300 * g.max()))
        return g2 if ok else None

    # ---------------------------------------------------------------- ascent
    def tangent_direction(self, g, R, sign):
        """Multiplicative-metric projected gradient of sign*log N on the tangent space of F."""
        ops = self.ops
        n = self.nr
        grad = sign * R["dN"] / R["N"]
        dK = R["dK"] / R["K"]
        J = np.empty((n + 1, n + 1))
        J[:n, :n] = ops.WGW(g)
        J[:n, n] = ops.Wv(g * dK)
        J[n, :n] = J[:n, n]
        J[n, n] = float(np.sum(g * dK * dK))
        rhs = np.append(ops.Wv(g * grad), float(np.sum(g * dK * grad)))
        J[np.diag_indices(n + 1)] += 1e-14 * np.max(np.abs(np.diag(J)))
        try:
            sol = sla.solve(J, rhs, assume_a="sym")
        except _LINALG_ERR:
            sol = np.linalg.lstsq(J, rhs, rcond=None)[0]
        lam, eta = sol[:n], sol[n]
        return grad - ops.WTv(lam) - eta * dK

    def run(self, g_init, sign, maxit=400, t0=0.5, verbose=False, record=False,
            already_feasible=False, window=30, rel_stall=1e-6):
        """sign=+1 maximises tau, sign=-1 minimises."""
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
        hist = []
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
            t = min(t * 2.0, 60.0)
            if record:
                hist.append((it, R["tau"], t))
            if verbose and it % 20 == 0:
                print(it, R["tau"], t, flush=True)
        return {"g": g, "tau": R["tau"], "K": R["K"], "N": R["N"], "feas": self.feas(g, R),
                "iters": it, "hist": hist}
