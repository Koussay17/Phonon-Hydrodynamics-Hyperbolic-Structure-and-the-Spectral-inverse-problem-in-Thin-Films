"""Core library: identifiability gap of the memory moment in a fixed event cone.

Conventions (entropy coordinates, as in notes 18 and 19 of the thesis repository):
  modes mu = 1..n, energies eps_mu > 0, Bose occupations N_mu at temperature kT,
  entropy scale d_mu = sqrt(N_mu (1 + N_mu)),
  event alpha with stoichiometric vector s_alpha (integer for exact events, slightly
  rescaled "conserving" coefficients for tolerance events), event vector
  a_alpha = D^{-1} s_alpha, operator C(g) = sum_alpha g_alpha a_alpha a_alpha^T, g >= 0,
  energy invariant e = D eps (C e = 0 iff s_alpha^T eps = 0 for every event),
  current b = D (v * eps) (one column per direction), b orthogonal to e.
  K(g) = b^T C^+ b,  N(g) = |C^+ b|^2,  tau(g) = N/K = -K'(0)/K(0).

C^+ b is computed as the solution of (C + gamma e_hat e_hat^T) x = b, which equals C^+ b
whenever ker C = span(e) and b is orthogonal to e (admissible g).
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp


# --------------------------------------------------------------------------- geometry
@dataclass
class EventGeometry:
    eps: np.ndarray            # (n,) mode energies (units of the energy scale)
    kT: float                  # temperature in the same units
    vel: np.ndarray            # (n, ndir) group velocities
    events: np.ndarray         # (m, 3) mode indices (third index may repeat the second: -1 means absent)
    coef: np.ndarray           # (m, 3) stoichiometric coefficients (s_alpha entries)
    kind: np.ndarray           # (m,) 0 normal, 1 umklapp
    gphys: np.ndarray          # (m,) reference ("physical-model") event rates
    labels: dict = field(default_factory=dict)

    def __post_init__(self):
        x = self.eps / self.kT
        self.nocc = 1.0 / np.expm1(x)
        self.dscale = np.sqrt(self.nocc * (1.0 + self.nocc))
        n, m = len(self.eps), len(self.events)
        self.n, self.m = n, m
        rows, cols, vals = [], [], []
        for j in range(3):
            idx = self.events[:, j]
            ok = idx >= 0
            rows.append(idx[ok]); cols.append(np.nonzero(ok)[0])
            vals.append(self.coef[ok, j] / self.dscale[idx[ok]])
        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        A = sp.coo_matrix((vals, (rows, cols)), shape=(n, m)).tocsc()
        A.sum_duplicates()
        self.A = A                                   # columns a_alpha
        self.W = A.multiply(A).tocsr()               # diag(C(g)) = W g
        self.e = self.dscale * self.eps              # energy invariant
        self.ehat = self.e / np.linalg.norm(self.e)
        self.b = (self.dscale * self.eps)[:, None] * self.vel   # (n, ndir)

    # ------------------------------------------------------------------ operator
    def C(self, g):
        Ag = self.A.multiply(np.asarray(g)[None, :])
        return (Ag @ self.A.T).toarray()

    def solve_factory(self, g):
        C = self.C(g)
        gam = max(np.mean(np.diag(C)), 1e-300)
        Creg = C + gam * np.outer(self.ehat, self.ehat)
        cf = sla.cho_factor(Creg, lower=True, check_finite=False)
        return C, cf

    def response(self, g, col=0, want_grad=False):
        """Return dict with K, N, tau, x, and optionally gradients w.r.t. g."""
        C, cf = self.solve_factory(g)
        b = self.b[:, col]
        x = sla.cho_solve(cf, b, check_finite=False)
        K = float(b @ x)
        N = float(x @ x)
        out = {"K": K, "N": N, "tau": N / K, "x": x}
        if want_grad:
            u = sla.cho_solve(cf, x, check_finite=False)       # C^+ x  (x in H, so u in H)
            p = self.A.T @ x                                    # a_alpha^T x
            q = self.A.T @ u                                    # a_alpha^T u
            out["dK"] = -(p * p)
            out["dN"] = -2.0 * p * q
            out["p"] = p
            out["u"] = u
        return out

    def check_kernel(self, g, tol=1e-10):
        C = self.C(g)
        w = np.linalg.eigvalsh(C)
        scale = max(w.max(), 1e-300)
        return int(np.sum(w < tol * scale)), w


# --------------------------------------------------------------------------- utilities
def bose(x):
    return 1.0 / np.expm1(x)


def rta(geom: EventGeometry, r, col=0):
    """RTA-type sums with diagonal r (no conservation): K_RTA = sum b^2/r, N_RTA = sum b^2/r^2."""
    b = geom.b[:, col]
    K = float(np.sum(b * b / r))
    N = float(np.sum(b * b / r**2))
    return K, N, N / K


# --------------------------------------------------------------------------- I-projection
def iproject(geom: EventGeometry, g0, r, K0=None, col=0, tol=1e-13, maxit=200, verbose=False):
    """Bregman (KL) projection of g0 onto {W g = r} (and K(g) = K0 if given).

    g = g0 * exp(W^T lam + eta * h), Newton on (lam, eta).  h is the normalised
    K-gradient magnitude at g0 (increasing g where the potential drop is large lowers K).
    Returns g, info.
    """
    W = geom.W
    n = geom.n
    lam = np.zeros(n)
    eta = 0.0
    use_K = K0 is not None
    if use_K:
        res0 = geom.response(g0, col, want_grad=True)
        h = res0["p"] ** 2
        h = h / max(h.max(), 1e-300)
    rs = np.maximum(r, 1e-300)

    def make(lam, eta):
        z = W.T @ lam
        if use_K:
            z = z + eta * h
        z = np.clip(z, -700, 700)
        return g0 * np.exp(z)

    def resid(g):
        F = (W @ g - r) / rs
        if use_K:
            Kg = geom.response(g, col)["K"]
            F = np.append(F, math.log(Kg / K0))
        return F

    g = make(lam, eta)
    F = resid(g)
    it = 0
    for it in range(maxit):
        nF = np.linalg.norm(F, np.inf)
        if verbose:
            print(it, nF)
        if nF < tol:
            break
        # Jacobian
        Dg = sp.diags(g)
        J11 = ((W @ Dg) @ W.T).toarray() / rs[:, None]
        if use_K:
            resg = geom.response(g, col, want_grad=True)
            dlogK = resg["dK"] / resg["K"]
            J12 = (W @ (g * h)) / rs
            J21 = (W @ (g * dlogK))
            J22 = float(np.sum(g * h * dlogK))
            J = np.block([[J11, J12[:, None]], [J21[None, :], np.array([[J22]])]])
        else:
            J = J11
        try:
            step = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(J, -F, rcond=None)[0]
        t = 1.0
        while t > 1e-8:
            lam_t = lam + t * step[:n]
            eta_t = eta + (t * step[n] if use_K else 0.0)
            g_t = make(lam_t, eta_t)
            try:
                F_t = resid(g_t)
            except (np.linalg.LinAlgError, ValueError):
                t *= 0.5
                continue
            if np.all(np.isfinite(F_t)) and np.linalg.norm(F_t) < (1 - 1e-4 * t) * np.linalg.norm(F):
                break
            t *= 0.5
        lam, eta, g, F = lam_t, eta_t, g_t, F_t
    info = {"iters": it, "resid_inf": float(np.linalg.norm(F, np.inf)), "eta": eta}
    return g, info


# --------------------------------------------------------------------------- BT constant
def reduced_basis(geom: EventGeometry):
    """Orthonormal basis Q (n x (n-1)) of H = e^perp."""
    n = geom.n
    # Householder-style: complete ehat to an orthonormal basis
    Qfull, _ = np.linalg.qr(np.column_stack([geom.ehat, np.eye(n)[:, : n - 1]]))
    Q = Qfull[:, 1:n]
    # make sure columns orthogonal to ehat
    Q = Q - np.outer(geom.ehat, geom.ehat @ Q)
    Q, _ = np.linalg.qr(Q)
    return Q


def bt_basic_solution(Ared, bred, I):
    """y with a_alpha^T y = 0 (alpha in I), b^T y = 1 in reduced coordinates; None if singular."""
    B = np.vstack([Ared[:, list(I)].T, bred[None, :]])
    try:
        cond = np.linalg.cond(B)
        if not np.isfinite(cond) or cond > 1e13:
            return None
        rhs = np.zeros(B.shape[0]); rhs[-1] = 1.0
        return np.linalg.solve(B, rhs)
    except np.linalg.LinAlgError:
        return None


def bt_exact_enumeration(geom: EventGeometry, col=0, max_subsets=5_000_000):
    """Exact M = max_I |y_I| by enumerating all (d-1)-subsets (small sets only)."""
    Q = reduced_basis(geom)
    Ared = Q.T @ geom.A.toarray()
    bred = Q.T @ geom.b[:, col]
    d = Ared.shape[0]
    m = Ared.shape[1]
    total = math.comb(m, d - 1)
    if total > max_subsets:
        raise ValueError(f"{total} subsets exceed limit")
    best, bestI, count = 0.0, None, 0
    for I in itertools.combinations(range(m), d - 1):
        y = bt_basic_solution(Ared, bred, I)
        if y is None:
            continue
        count += 1
        nv = float(np.linalg.norm(y))
        if nv > best:
            best, bestI = nv, I
    return {"M": best, "I": bestI, "n_bases": count, "n_subsets": total, "d": d, "m": m}


def bt_local_search(geom: EventGeometry, col=0, n_starts=20, rng=None, max_pivots=10_000, init=None):
    """Lower bound on M by best-improvement pivoting over bases (each y_I is an explicit witness)."""
    rng = np.random.default_rng(rng)
    Q = reduced_basis(geom)
    Ared = Q.T @ geom.A.toarray()          # d x m
    bred = Q.T @ geom.b[:, col]
    d, m = Ared.shape
    results = []
    for s in range(n_starts):
        # random basis: greedy independent selection in random order
        order = rng.permutation(m) if (init is None or s > 0) else np.asarray(init)
        I = []
        Mcur = np.zeros((0, d))
        for a in order:
            cand = np.vstack([Mcur, Ared[:, a]])
            if np.linalg.matrix_rank(np.vstack([cand, bred])) == cand.shape[0] + 1:
                Mcur = cand
                I.append(a)
                if len(I) == d - 1:
                    break
        if len(I) < d - 1:
            continue
        for piv in range(max_pivots):
            B = np.vstack([Ared[:, I].T, bred[None, :]])
            Binv = np.linalg.inv(B)                     # columns: c_gamma (gamma in I), y_I
            y = Binv[:, -1]
            Cb = Binv[:, :-1]                           # d x (d-1)
            ay = Ared.T @ y                             # (m,)
            aC = Ared.T @ Cb                            # (m, d-1)
            ny2 = y @ y
            cc = np.sum(Cb * Cb, axis=0)                # |c_beta|^2
            yc = y @ Cb                                 # y^T c_beta
            with np.errstate(divide="ignore", invalid="ignore"):
                t = -ay[:, None] / aC
                val = ny2 + t * t * cc[None, :] + 2 * t * yc[None, :]
            bad = ~np.isfinite(val) | (np.abs(aC) < 1e-12 * np.abs(aC).max())
            val[bad] = -np.inf
            val[I, :] = -np.inf
            k = np.unravel_index(np.argmax(val), val.shape)
            if val[k] <= ny2 * (1 + 1e-12):
                break
            I[k[1]] = k[0]
        y = np.linalg.solve(np.vstack([Ared[:, I].T, bred[None, :]]), np.eye(d)[:, -1])
        results.append((float(np.linalg.norm(y)), list(I)))
    results.sort(key=lambda t: -t[0])
    return results
