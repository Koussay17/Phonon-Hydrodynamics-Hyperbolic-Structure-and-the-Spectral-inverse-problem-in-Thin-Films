"""Independent 1D FDTR side for the slab audit (own implementation; no import of au5_fisher1d / au6_slab / rs2b).

Stack (top to bottom): Au 80 nm (150 W/mK, 2.49e6 J/m^3K) / R1 / AlN (slab two-port or GK film) / R2 / sapphire
(35 W/mK, 3.03e6 J/m^3K, semi-infinite); heat absorbed at the Au surface; Z = theta/q at the Au surface.
Observable ln H = ln Z + g + i (phi0 - 2 pi f t_d); noise 1 % on ln|H|, sigma_phi on the phase; 80 log-spaced
frequencies 10 kHz - 200 MHz (the declared scenario).
GK film (local interface temperature): lam_eff = lam_z (1 + tau_l p)/(1 + tau_R p) in a homogeneous layer of capacity C.
Slab two-port, 'wall' port: R1 and R2 act on the wall emission temperatures T_A, T_B;
(q(0), q(d)) = Y (T_A, T_B), T_B = (Z_sub + R2) q(d).
Fits: the three instrument nuisances enter linearly in the whitened residual and are projected out exactly
(variable projection with a fixed orthogonal projector); the nonlinear parameters (tau_l and the set's nuisances
among ln lam_z, ln C, ln R1, ln R2) by Levenberg-Marquardt (scipy least_squares) from several starts; profile of the
misfit over delta on a grid as a global-minimum check.
Fisher: derivatives of ln H by Cauchy-integral differentiation on a circle in the complex parameter plane (trapezoidal
rule, M points), independent of finite differences.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

AU = (150.0, 2.49e6, 80e-9)
SAP = (35.0, 3.03e6)
FREQ = np.geomspace(1e4, 2e8, 80)
SIG_A = 0.01
# parameter vector (nonlinear): 0 tau_R [s], 1 tau_l [s], 2 ln lam_z, 3 ln C, 4 ln R1, 5 ln R2
NL_SETS = {"N1": [], "N2": [4, 5], "N2lam": [2, 4, 5], "N3": [2, 3, 4, 5]}


def quad_layer(Zb, lam, C, d, p):
    g = np.sqrt(C * p / lam)
    ch, sh = np.cosh(g * d), np.sinh(g * d)
    Yl = lam * g
    return (ch * Zb + sh / Yl) / (Yl * sh * Zb + ch)


def Z_gk(par, d, p):
    tR, tl, llam, lC, lR1, lR2 = par
    lam = np.exp(llam) * (1 + tl * p) / (1 + tR * p)
    Zs = 1.0 / np.sqrt(SAP[0] * SAP[1] * p)
    Z = quad_layer(Zs + np.exp(lR2), lam, np.exp(lC), d, p) + np.exp(lR1)
    return quad_layer(Z, AU[0], AU[1], AU[2], p)


def Z_slab(Y, R1, R2, p, port="wall", Th=None):
    """Y, Th: arrays (nf, 2, 2)."""
    Zs = 1.0 / np.sqrt(SAP[0] * SAP[1] * p)
    Zb = Zs + R2
    if port == "wall":
        r = Y[:, 1, 0] * Zb / (1 - Y[:, 1, 1] * Zb)          # T_B per T_A
        q0 = Y[:, 0, 0] + Y[:, 0, 1] * r
        Zt = 1.0 / q0
    else:   # R1, R2 on the energy temperatures at the walls
        row0 = Th[:, 1, 0] - Zb * Y[:, 1, 0]
        row1 = Th[:, 1, 1] - Zb * Y[:, 1, 1]
        r = -row0 / row1
        q0 = Y[:, 0, 0] + Y[:, 0, 1] * r
        Zt = (Th[:, 0, 0] + Th[:, 0, 1] * r) / q0
    return quad_layer(Zt + R1, AU[0], AU[1], AU[2], p)


def whiten(z, sp_deg):
    return np.concatenate([z.real / SIG_A, z.imag / np.deg2rad(sp_deg)])


def lin_projector(sp_deg, freq=FREQ):
    """Orthonormal basis of the whitened instrument columns (g, phi0, t_d)."""
    nf = freq.size
    sp = np.deg2rad(sp_deg)
    cols = np.zeros((2 * nf, 3))
    cols[:nf, 0] = 1.0 / SIG_A
    cols[nf:, 1] = 1.0 / sp
    cols[nf:, 2] = -2 * np.pi * freq * 1e-9 / sp
    Q, _ = np.linalg.qr(cols)
    return Q


class Fitter:
    def __init__(self, y, d, base, sp_deg, freq=FREQ, weights=None):
        self.y, self.d, self.base, self.sp = y, d, np.array(base, float), sp_deg
        self.freq = freq
        self.p = 2j * np.pi * freq
        self.Q = lin_projector(sp_deg, freq)
        self.w = None if weights is None else np.concatenate([weights, weights])
        if self.w is not None:     # weighted: rebuild the projector on the weighted columns
            nf = freq.size
            sp = np.deg2rad(sp_deg)
            cols = np.zeros((2 * nf, 3))
            cols[:nf, 0] = 1.0 / SIG_A
            cols[nf:, 1] = 1.0 / sp
            cols[nf:, 2] = -2 * np.pi * freq * 1e-9 / sp
            cols *= self.w[:, None]
            self.Q, _ = np.linalg.qr(cols)

    def resid(self, par):
        z = np.log(Z_gk(par, self.d, self.p))
        r = whiten(self.y - z, self.sp)
        if self.w is not None:
            r = r * self.w
        return r - self.Q @ (self.Q.T @ r)

    def chi2(self, par):
        r = self.resid(par)
        return float(r @ r)

    def fit(self, free, fixed_par, starts):
        """Minimise over the free nonlinear indices from several starting vectors (full parameter vectors);
        returns best chi2, best parameter vector, and all (chi2, par) results."""
        out = []
        for st in starts:
            par0 = np.array(fixed_par, float)
            par0[free] = np.array(st, float)[free]
            sc = np.array([1e-12 if i in (0, 1) else 1.0 for i in free])

            def fun(x):
                pr = par0.copy()
                pr[free] = x * sc
                r = self.resid(pr)
                return r if np.all(np.isfinite(r)) else np.full(r.size, 1e6)
            if len(free) == 0:
                out.append((self.chi2(par0), par0))
                continue
            sol = least_squares(fun, par0[free] / sc, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
            pr = par0.copy()
            pr[free] = sol.x * sc
            out.append((float(sol.fun @ sol.fun), pr))
        best = min(out, key=lambda t: t[0])
        return best[0], best[1], out


def nuisance_starts(base, nl):
    """Starting vectors for the nonlinear nuisances: base and multiplicative excursions."""
    sts = [np.array(base, float)]
    for fac in (0.75, 1.33):
        for k in nl:
            s = np.array(base, float)
            s[k] += np.log(fac)
            sts.append(s)
    s = np.array(base, float)
    for k in nl:
        s[k] += np.log(0.85) if k != 3 else 0.0
    sts.append(s)
    return sts


def fit_delta(y, d, base, sp_deg, Nset, delta_starts_ps=(0.0, -30.0, -60.0, -100.0, -200.0, 50.0), freq=FREQ,
              weights=None):
    """GK fit with tau_R fixed at base[0] and tau_l free; returns delta_fit (ps), chi2, parameters, spread over
    starts that reach the minimum."""
    F = Fitter(y, d, base, sp_deg, freq, weights)
    nl = NL_SETS[Nset]
    free = [1] + nl
    starts = []
    for dl in delta_starts_ps:
        for s in nuisance_starts(base, nl):
            s = s.copy()
            s[1] = base[0] - dl * 1e-12
            starts.append(s)
    chi, par, allr = F.fit(free, base, starts)
    deltas = [(base[0] - pr[1]) * 1e12 for c, pr in allr if c <= chi + 1e-6 * max(1.0, chi)]
    return {"delta_ps": (base[0] - par[1]) * 1e12, "chi2": chi, "par": par.tolist(), "n_starts": len(starts),
            "n_at_min": len(deltas), "spread_ps": float(max(deltas) - min(deltas)),
            "worst_chi2": float(max(c for c, _ in allr))}


def profile_delta(y, d, base, sp_deg, Nset, grid_ps, warm=None, freq=FREQ):
    """Misfit profile chi2(delta) = min over the set's nonlinear nuisances at fixed tau_l = tau_R - delta."""
    F = Fitter(y, d, base, sp_deg, freq)
    nl = NL_SETS[Nset]
    prof = []
    prev = np.array(base if warm is None else warm, float)
    for dl in grid_ps:
        fixed = prev.copy()
        fixed[1] = base[0] - dl * 1e-12
        sts = nuisance_starts(base, nl) + [prev]
        sts2 = []
        for s in sts:
            s = s.copy()
            s[1] = fixed[1]
            sts2.append(s)
        chi, par, _ = F.fit(nl, fixed, sts2)
        prof.append(chi)
        prev = par
    return np.array(prof)


# ------------------------------------------------------------------------------------------------ Fisher
def cauchy_jac(par, d, M=16, rel=0.05, freq=FREQ):
    """Columns d lnH / d par_k (k = 0..5; tau in seconds) by Cauchy-integral differentiation, plus the three
    instrument columns (g, phi0, t_d [ns])."""
    p = 2j * np.pi * freq
    par = np.array(par, float)
    cols = []
    om = np.exp(2j * np.pi * np.arange(M) / M)
    for k in range(6):
        r = rel * (max(abs(par[k]), 1e-12) if k in (0, 1) else 1.0)
        acc = 0.0
        for w in om:
            pk = par.astype(complex)
            pk[k] = par[k] + r * w
            acc = acc + np.log(Z_gk(pk, d, p)) / (r * w)
        cols.append(acc / M)
    cols.append(np.ones(freq.size, complex))
    cols.append(np.full(freq.size, 1j))
    cols.append(-1j * 2 * np.pi * freq * 1e-9)
    return np.stack(cols, 1)


def sigma_direction(J, target_col, nuis_cols, sp_deg):
    """1/|| (I - P_N) J_target || in whitened units (marginal standard deviation of the target parameter)."""
    Jw = np.concatenate([J.real / SIG_A, J.imag / np.deg2rad(sp_deg)], 0)
    v = Jw[:, target_col]
    Nm = Jw[:, nuis_cols]
    U, s, _ = np.linalg.svd(Nm, full_matrices=False)
    U = U[:, s > s[0] * 1e-13]
    r = v - U @ (U.T @ v)
    return 1.0 / np.linalg.norm(r)


SETCOLS = {"N1": [6, 7, 8], "N2": [6, 7, 8, 4, 5], "N2lam": [6, 7, 8, 4, 5, 2], "N3": [6, 7, 8, 4, 5, 2, 3]}
