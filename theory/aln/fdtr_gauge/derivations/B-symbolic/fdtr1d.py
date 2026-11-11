"""Front-face impedance of a Cattaneo / Guyer-Krumhansl film on a semi-infinite
Fourier substrate (note 7 model), in three independent implementations.

Notation (all square roots principal; on p = i*omega, omega > 0, no branch cut
is crossed, see B-symbolic.md section 1):

    rho = sqrt((1 + tau_R p) / (1 + tau_l p))
    A   = xi1 sqrt(p) rho                 (= sigma e)
    B   = b sqrt(p) / rho                 (= lambda_eff sigma)
    S   = bs sqrt(p)                      (substrate admittance)
    Zs  = R + 1/S                         (substrate impedance, optional interface R)
    w   = B Zs,  t = tanh(A)

Quadrupole of the film  M = [[cosh A, sinh A / B], [B sinh A, cosh A]],
front-face impedance    Z = (M00 Zs + M01) / (M10 Zs + M11)
                          = (1/B) (w + t) / (1 + w t).
With R = 0:             Z = (1/B) (B + S t) / (S + B t).

The expression given in the branch brief, (1/B)(B t + S)/(B + S t), is kept as
Z_brief for the discrepancy check only.

theta = (xi1, b, tau_R, tau_l) [, R].  bs is known.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from fractions import Fraction

import mpmath as mp

PARAMS4 = ("xi1", "b", "tau_R", "tau_l")
PARAMS5 = ("xi1", "b", "tau_R", "tau_l", "R")


# ---------------------------------------------------------------------------
# mpmath implementations (arbitrary precision floating point)
# ---------------------------------------------------------------------------

def _mpf(x):
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / x.denominator
    return mp.mpmathify(x)


def parts_mp(p, xi1, b, tR, tl):
    p = mp.mpmathify(p)
    rho = mp.sqrt((1 + tR * p) / (1 + tl * p))
    sp = mp.sqrt(p)
    return xi1 * sp * rho, b * sp / rho, sp


def Z_mp(p, xi1, b, tR, tl, bs, R=0):
    """Closed form (1/B)(w + t)/(1 + w t)."""
    xi1, b, tR, tl, bs, R = map(_mpf, (xi1, b, tR, tl, bs, R))
    A, B, sp = parts_mp(p, xi1, b, tR, tl)
    Zs = R + 1 / (bs * sp)
    w = B * Zs
    t = mp.tanh(A)
    return (w + t) / (B * (1 + w * t))


def Z_quad_mp(p, xi1, b, tR, tl, bs, R=0):
    """Direct quadrupole product closed by the substrate impedance."""
    xi1, b, tR, tl, bs, R = map(_mpf, (xi1, b, tR, tl, bs, R))
    A, B, sp = parts_mp(p, xi1, b, tR, tl)
    Zs = R + 1 / (bs * sp)
    ch, sh = mp.cosh(A), mp.sinh(A)
    M = mp.matrix([[ch, sh / B], [B * sh, ch]])
    if R != 0:
        M = M * mp.matrix([[1, R], [0, 1]])
        Zs = 1 / (bs * sp)
    return (M[0, 0] * Zs + M[0, 1]) / (M[1, 0] * Zs + M[1, 1])


def Z_brief_mp(p, xi1, b, tR, tl, bs):
    """Expression of the brief: (1/B)(B tanh A + bs sqrt p)/(B + bs sqrt p tanh A)."""
    xi1, b, tR, tl, bs = map(_mpf, (xi1, b, tR, tl, bs))
    A, B, sp = parts_mp(p, xi1, b, tR, tl)
    S = bs * sp
    t = mp.tanh(A)
    return (B * t + S) / (B * (B + S * t))


def jac_fd_mp(p, theta, bs, model="quad", R=None, dps_inner=None):
    """Jacobian dZ/dtheta by mpmath numerical differentiation (no symbolic input).

    theta: 4 or 5 values (xi1, b, tau_R, tau_l[, R]).  Uses the quadrupole
    product (model='quad') or the closed form (model='closed').
    Differentiation is carried out at dps_inner (default 2*mp.dps) digits.
    """
    f = Z_quad_mp if model == "quad" else Z_mp
    theta = [_mpf(x) for x in theta]
    n = len(theta)
    out = []
    dps0 = mp.mp.dps
    with mp.workdps(dps_inner or 2 * dps0):
        for k in range(n):
            def g(x, k=k):
                th = list(theta)
                th[k] = x
                if n == 5:
                    return f(p, th[0], th[1], th[2], th[3], bs, th[4])
                return f(p, th[0], th[1], th[2], th[3], bs, R or 0)
            out.append(mp.diff(g, theta[k]))
    return out


# ---------------------------------------------------------------------------
# python-flint ball arithmetic (rigorous enclosures), forward-mode derivatives
# ---------------------------------------------------------------------------

def _flint():
    import flint
    return flint


def to_arb(x):
    flint = _flint()
    if isinstance(x, Fraction):
        return flint.arb(flint.fmpq(x.numerator, x.denominator))
    if isinstance(x, int):
        return flint.arb(x)
    raise TypeError("ball arithmetic inputs must be exact (int or Fraction)")


def Z_ball(p, xi1, b, tR, tl, bs, R):
    """Z with acb / acb_series arguments (generic in the arithmetic type)."""
    flint = _flint()
    acb, acb_series = flint.acb, flint.acb_series

    def sqrt(x):
        return x.sqrt() if isinstance(x, acb) else acb_series.sqrt(x)

    def exp(x):
        return x.exp() if isinstance(x, acb) else acb_series.exp(x)

    sp = p.sqrt()
    rho = sqrt((1 + tR * p) / (1 + tl * p))
    A = xi1 * sp * rho
    B = b * sp / rho
    Zs = R + 1 / (bs * sp)
    w = B * Zs
    e2 = exp(2 * A)
    t = 1 - 2 / (e2 + 1)          # tanh(A), stable for Re A > 0
    return (w + t) / (B * (1 + w * t))


def jac_ball(omega, theta, bs, R=Fraction(0)):
    """Rigorous enclosures of dZ/dtheta_k at p = i*omega (forward-mode AD).

    theta: tuple of 4 (R fixed, given) or 5 exact rationals.  Returns list of acb.
    """
    flint = _flint()
    acb, acb_series = flint.acb, flint.acb_series
    p = acb(0, to_arb(omega))
    th = [acb(to_arb(x)) for x in theta]
    n = len(theta)
    if n == 4:
        th.append(acb(to_arb(R)))
    bsb = acb(to_arb(bs))
    cols = []
    for k in range(n):
        args = list(th)
        args[k] = acb_series([th[k], 1], prec=2)
        z = Z_ball(p, args[0], args[1], args[2], args[3], bsb, args[4])
        c = z.coeffs()
        cols.append(c[1] if len(c) > 1 else acb(0))
    return cols


def Z_ball_value(omega, theta, bs, R=Fraction(0)):
    flint = _flint()
    acb = flint.acb
    p = acb(0, to_arb(omega))
    th = [acb(to_arb(x)) for x in theta]
    if len(th) == 4:
        th.append(acb(to_arb(R)))
    return Z_ball(p, th[0], th[1], th[2], th[3], acb(to_arb(bs)), th[4])


# ---------------------------------------------------------------------------
# SymPy expression (built on demand)
# ---------------------------------------------------------------------------

def sympy_model(with_R=True):
    import sympy as sp
    p = sp.Symbol("p")                      # complex, no assumption
    xi1, b, bs = sp.symbols("xi1 b b_s", positive=True)
    tR, tl, R = sp.symbols("tau_R tau_l R", nonnegative=True)
    rho = sp.sqrt((1 + tR * p) / (1 + tl * p))
    A = xi1 * sp.sqrt(p) * rho
    B = b * sp.sqrt(p) / rho
    S = bs * sp.sqrt(p)
    Zs = (R if with_R else 0) + 1 / S
    t = sp.tanh(A)
    Z = (Zs + t / B) / (1 + B * Zs * t)
    syms = dict(p=p, xi1=xi1, b=b, bs=bs, tau_R=tR, tau_l=tl, R=R)
    return Z, syms, dict(rho=rho, A=A, B=B, S=S, Zs=Zs)


def sympy_jacobian_mp(with_R=True):
    """Lambdified symbolic derivatives (mpmath backend).

    Returns f(p, xi1, b, tR, tl, bs, R) -> list of dZ/dtheta_k.
    """
    import sympy as sp
    Z, s, _ = sympy_model(with_R=True)
    params = [s["xi1"], s["b"], s["tau_R"], s["tau_l"]] + ([s["R"]] if with_R else [])
    derivs = [sp.diff(Z, q) for q in params]
    args = (s["p"], s["xi1"], s["b"], s["tau_R"], s["tau_l"], s["b_s"] if "b_s" in s else s["bs"], s["R"])
    fn = sp.lambdify(args, derivs, modules="mpmath")

    def f(p, xi1, b, tR, tl, bs, R=0):
        return fn(mp.mpmathify(p), *map(_mpf, (xi1, b, tR, tl, bs, R)))
    return f


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def real_jacobian(rows):
    """rows: list (over frequencies) of complex gradient lists -> real 2N x n."""
    out = []
    for r in rows:
        out.append([mp.re(x) for x in r])
        out.append([mp.im(x) for x in r])
    return mp.matrix(out)


def complex_jacobian(rows):
    return mp.matrix([list(r) for r in rows])


def equilibrate(M):
    """Scale columns to unit 2-norm (rank-preserving)."""
    M = M.copy()
    for j in range(M.cols):
        nrm = mp.sqrt(sum(abs(M[i, j]) ** 2 for i in range(M.rows)))
        if nrm != 0:
            for i in range(M.rows):
                M[i, j] /= nrm
    return M


def singular_values(M, complex_=False):
    if complex_:
        s = mp.svd_c(M, compute_uv=False)
    else:
        s = mp.svd_r(M, compute_uv=False)
    return sorted([abs(s[i]) for i in range(len(s))], reverse=True)
