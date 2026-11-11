"""Shared 1D models for branch AB-A (asymptotic regimes and phase-sensitivity scaling).

Conventions (notes/conventions.md, src/forward_model.py, src/twin_fdtr.py): exp(+i omega t), p = i omega,
quantities per unit area, front-face flux, response = front-face temperature per unit absorbed flux (m^2 K/W).

Film on a semi-infinite Fourier substrate (workstream B):
    g = sqrt((1 + tau_R p)/(1 + tau_l p)),  A = xi1 sqrt(p) g,  B = b sqrt(p)/g,  S = b_s sqrt(p),
    Z = (1/B) (B + S tanh A)/(S + B tanh A)        (= zin of forward_model.response with R = 0)
NB: the brief writes (B tanh A + S)/(B + S tanh A); that form has the thin-film limit S/B^2, not 1/S.

Full 1D stack (workstream A): Au (Fourier) / R1 / AlN (Cattaneo or GK, 1D lambda_eff form) / R2 / sapphire.
The AlN layer is a Fourier layer with the complex conductivity lambda_eff(p) = lambda (1 + tau_l p)/(1 + tau_R p)
at fixed C and thickness (exact in the 1D model).
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np

# --------------------------------------------------------------------------------------------
# Workstream B: film on substrate, analytic derivatives
# --------------------------------------------------------------------------------------------

def film_parts(p, xi1, b, tau_r, tau_l, bs):
    p = np.asarray(p, dtype=complex)
    G = (1.0 + tau_r * p) / (1.0 + tau_l * p)
    g = np.sqrt(G)
    sp = np.sqrt(p)
    A = xi1 * sp * g
    B = b * sp / g
    S = bs * sp
    return G, g, A, B, S


def _tanh_sech2(A):
    """tanh A and sech^2 A for Re A >= 0 without overflow (E = exp(-2A), |E| <= 1)."""
    E = np.exp(-2.0 * A)
    return (1.0 - E) / (1.0 + E), 4.0 * E / (1.0 + E) ** 2


def film_Z(p, xi1, b, tau_r, tau_l, bs):
    G, g, A, B, S = film_parts(p, xi1, b, tau_r, tau_l, bs)
    t, _ = _tanh_sech2(A)
    return (B + S * t) / (B * (S + B * t))


def film_reflection(p, xi1, b, tau_r, tau_l, bs):
    """Gamma = (B - S)/(B + S) and E = exp(-2A); Z = (1/B)(1 + Gamma E)/(1 - Gamma E)."""
    G, g, A, B, S = film_parts(p, xi1, b, tau_r, tau_l, bs)
    return (B - S) / (B + S), np.exp(-2.0 * A), A, B


def film_jacobian(p, xi1, b, tau_r, tau_l, bs, tscale=None):
    """Analytic d ln Z / d(ln xi1, ln b, tau_R/ts, tau_l/ts), complex, shape (len(p), 4).

    Z(A, B) = (B + S t)/(B (S + B t)), t = tanh A.
      dZ/dA = (S^2 - B^2) sech^2 A / (B (S + B t)^2)
      dZ/dB = -t (B^2 + S^2 + 2 B S t) / (B (S + B t))^2
    A ~ xi1 sqrt(G), B ~ b / sqrt(G), G = (1 + tau_R p)/(1 + tau_l p).
    """
    p = np.asarray(p, dtype=complex)
    G, g, A, B, S = film_parts(p, xi1, b, tau_r, tau_l, bs)
    t, sech2 = _tanh_sech2(A)
    D = S + B * t
    Z = (B + S * t) / (B * D)
    dZdA = (S * S - B * B) * sech2 / (B * D * D)
    dZdB = -t * (B * B + S * S + 2 * B * S * t) / (B * D) ** 2
    dZ_dlnxi = A * dZdA
    dZ_dlnb = B * dZdB
    dZ_dG = (A * dZdA - B * dZdB) / (2.0 * G)
    dG_dtr = p / (1.0 + tau_l * p)
    dG_dtl = -p * (1.0 + tau_r * p) / (1.0 + tau_l * p) ** 2
    ts = xi1 ** 2 if tscale is None else tscale
    cols = np.stack([dZ_dlnxi, dZ_dlnb, dZ_dG * dG_dtr * ts, dZ_dG * dG_dtl * ts], axis=-1)
    return Z, cols / Z[:, None]


# --------------------------------------------------------------------------------------------
# Workstream A: Au / R1 / AlN / R2 / sapphire, 1D
# --------------------------------------------------------------------------------------------

BASE = dict(
    lam_au=150.0, C_au=2.49e6, e_au=80e-9,
    R1=1e-8,
    lam=200.0, C=2.41e6, e=500e-9, tau_r=27.5e-12, tau_l=0.0,
    R2=2e-8,
    lam_s=35.0, C_s=3.03e6,
)


def _layer(Zb, A, B):
    """Input impedance of a homogeneous layer (A = sigma e, B = lambda sigma) on load Zb (stable tanh form)."""
    t, _ = _tanh_sech2(A)
    return (Zb + t / B) / (1.0 + B * Zb * t)


def _sech(A):
    E = np.exp(-A)
    return 2.0 * E / (1.0 + E * E)


def _layer_transfer(Zb, A, B):
    """q_bottom / q_top of a layer on load Zb: 1/(cosh A + B Zb sinh A), evaluated stably."""
    # 1/(cosh A + B Zb sinh A) = sech A / (1 + B Zb tanh A)
    t, _ = _tanh_sech2(A)
    return _sech(A) / (1.0 + B * Zb * t)


def aln_AB(p, P):
    lam_eff = P["lam"] * (1.0 + P["tau_l"] * p) / (1.0 + P["tau_r"] * p)
    G = P["lam"] / lam_eff                      # = g^2
    g = np.sqrt(G)
    sp = np.sqrt(p)
    xi1 = P["e"] * np.sqrt(P["C"] / P["lam"])
    b = np.sqrt(P["lam"] * P["C"])
    return xi1 * sp * g, b * sp / g


def stack_Z(omega, P=BASE, return_parts=False):
    """Front-face impedance of the 1D stack, m^2 K/W, at p = i omega."""
    p = 1j * np.asarray(omega, dtype=float)
    bs = np.sqrt(P["lam_s"] * P["C_s"])
    Zs = 1.0 / (bs * np.sqrt(p))
    Z2 = Zs + P["R2"]                            # seen at the AlN bottom face
    Af, Bf = aln_AB(p, P)
    Zf = _layer(Z2, Af, Bf)                      # seen at the AlN top face
    Z1 = Zf + P["R1"]                            # seen at the Au bottom face
    Aa = P["e_au"] * np.sqrt(p * P["C_au"] / P["lam_au"])
    Ba = np.sqrt(P["lam_au"] * P["C_au"] * p)
    Z = _layer(Z1, Aa, Ba)
    if not return_parts:
        return Z
    return dict(p=p, Z=Z, Zs=Zs, Z2=Z2, Zf=Zf, Z1=Z1, Af=Af, Bf=Bf, Aa=Aa, Ba=Ba)


def stack_flux_transfer(omega, P=BASE, n_gauss=24):
    """Flux transfer T(z) = q(z)/q(front) through the AlN film (Gauss-Legendre nodes).

    Returns dict with T1 (AlN top), T2 (AlN bottom), mean of T^2 over the film, nodes.
    """
    parts = stack_Z(omega, P, return_parts=True)
    Ta = _layer_transfer(parts["Z1"], parts["Aa"], parts["Ba"])        # through the Au
    x, wq = np.polynomial.legendre.leggauss(n_gauss)
    u = 0.5 * (x + 1.0)                                                   # depth fraction in the film
    wq = 0.5 * wq
    Af, Bf, Z2 = parts["Af"], parts["Bf"], parts["Z2"]
    # impedance at depth u e inside the film: film of thickness (1-u) e on Z2
    Zu = _layer(Z2[:, None], Af[:, None] * (1.0 - u)[None, :], Bf[:, None])
    # q(u)/q(0) = 1/(cosh(A u) + B Z(u) sinh(A u))
    Au = Af[:, None] * u[None, :]
    Tu = _sech(Au) / (1.0 + Bf[:, None] * Zu * _tanh_sech2(Au)[0])
    Ttot = Ta[:, None] * Tu
    T2 = Ta * _layer_transfer(Z2, Af, Bf)
    meanT2 = (Ttot ** 2) @ wq
    return dict(T1=Ta, T2=T2, meanT2=meanT2, parts=parts)


def with_(P, **kw):
    Q = dict(P)
    Q.update(kw)
    return Q


def log_sensitivity(omega, P, name, h=1e-6):
    """d ln Z / d ln P[name] by a 4-point central difference in ln P[name] (positive parameters)."""
    v = P[name]
    def lz(f):
        return np.log(stack_Z(omega, with_(P, **{name: v * np.exp(f)})))
    # 4th-order central difference
    d = (-lz(2 * h) + 8 * lz(h) - 8 * lz(-h) + lz(-2 * h)) / (12 * h)
    return d


def abs_sensitivity(omega, P, name, h):
    """d ln Z / d P[name] (absolute step h), 4th-order central difference."""
    v = P[name]
    def lz(dv):
        return np.log(stack_Z(omega, with_(P, **{name: v + dv})))
    return (-lz(2 * h) + 8 * lz(h) - 8 * lz(-h) + lz(-2 * h)) / (12 * h)


# --------------------------------------------------------------------------------------------
# Forward-mode derivatives (dual numbers over complex arrays): exact to rounding, no cancellation
# --------------------------------------------------------------------------------------------

class Dual:
    """v + eps d with eps^2 = 0; v, d complex arrays (or scalars)."""
    __array_priority__ = 1000

    def __init__(self, v, d=0.0):
        self.v = np.asarray(v, dtype=complex)
        self.d = np.asarray(d, dtype=complex) + 0 * self.v

    @staticmethod
    def lift(x):
        return x if isinstance(x, Dual) else Dual(x, 0.0)

    def __add__(self, o):
        o = Dual.lift(o); return Dual(self.v + o.v, self.d + o.d)
    __radd__ = __add__

    def __sub__(self, o):
        o = Dual.lift(o); return Dual(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        return Dual.lift(o) - self

    def __mul__(self, o):
        o = Dual.lift(o); return Dual(self.v * o.v, self.d * o.v + self.v * o.d)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = Dual.lift(o); return Dual(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))

    def __rtruediv__(self, o):
        return Dual.lift(o) / self

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __getitem__(self, k):
        return Dual(self.v[k], self.d[k])


def dsqrt(x):
    if isinstance(x, Dual):
        r = np.sqrt(x.v); return Dual(r, x.d / (2 * r))
    return np.sqrt(np.asarray(x, dtype=complex))


def dexp(x):
    if isinstance(x, Dual):
        e = np.exp(x.v); return Dual(e, x.d * e)
    return np.exp(np.asarray(x, dtype=complex))


def _tanh_d(A):
    E = dexp(-2.0 * A)
    return (1.0 - E) / (1.0 + E)


def _layer_d(Zb, A, B):
    t = _tanh_d(A)
    return (Zb + t / B) / (1.0 + B * Zb * t)


def stack_Z_d(omega, P):
    """stack_Z with Dual-valued parameters allowed (same formulas, branch-safe for p = i omega)."""
    p = 1j * np.asarray(omega, dtype=float)
    sp_ = np.sqrt(p)
    bs = dsqrt(P["lam_s"] * P["C_s"])
    Z2 = 1.0 / (bs * sp_) + P["R2"]
    G = (1.0 + P["tau_r"] * p) / (1.0 + P["tau_l"] * p)
    g = dsqrt(G)
    xi1 = P["e"] * dsqrt(P["C"] / P["lam"])
    b = dsqrt(P["lam"] * P["C"])
    Af = xi1 * sp_ * g
    Bf = b * sp_ / g
    Zf = _layer_d(Z2, Af, Bf)
    Z1 = Zf + P["R1"]
    Aa = P["e_au"] * dsqrt(P["C_au"] / P["lam_au"]) * sp_
    Ba = dsqrt(P["lam_au"] * P["C_au"]) * sp_
    return _layer_d(Z1, Aa, Ba)


def dlnZ(omega, P, name, log=True):
    """Exact d ln Z / d ln P[name] (log=True) or d ln Z / d P[name] by forward-mode differentiation."""
    v = P[name]
    Q = dict(P)
    Q[name] = Dual(v, v if log else 1.0)
    Z = stack_Z_d(omega, Q)
    if not isinstance(Z, Dual):
        return np.zeros_like(np.asarray(omega, dtype=complex))
    return Z.d / Z.v


def jacobian_lnZ(omega, P, names, log_flags):
    return np.stack([dlnZ(omega, P, n, lg) for n, lg in zip(names, log_flags)], axis=-1)
