"""Front-face impedance of a Cattaneo/GK film on a semi-infinite Fourier substrate (mpmath).

Conventions (proofs/B-analytic.md, section 1):
  g(p) = (1 + tR p)/(1 + tl p),  A = xi sqrt(p g),  B = b sqrt(p/g),  S = bs sqrt(p)
  principal square roots; on Re p >= 0, p != 0 they have Re > 0 and A B = xi b p.
  Film quadrupole M = [[cosh A, sinh A / B], [B sinh A, cosh A]], load L = R + 1/S,
  Z = (cosh A L + sinh A / B)/(B sinh A L + cosh A) = (1/B)(1 + r E)/(1 - r E),
  r = (B L - 1)/(B L + 1),  E = exp(-2A).
Parameters theta = (xi, b, tR, tl) or (xi, b, tR, tl, R); bs is known.
"""
import mpmath as mp


def parts(p, xi, b, tR, tl, bs):
    p = mp.mpc(p)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g)
    B = b*mp.sqrt(p/g)
    S = bs*mp.sqrt(p)
    return g, A, B, S


def Z_quad(p, xi, b, tR, tl, bs, R=0):
    """Impedance from the quadrupole product (independent of the r/E form)."""
    g, A, B, S = parts(p, xi, b, tR, tl, bs)
    L = R + 1/S
    ch, sh = mp.cosh(A), mp.sinh(A)
    return (ch*L + sh/B)/(B*sh*L + ch)


def Z_refl(p, xi, b, tR, tl, bs, R=0):
    g, A, B, S = parts(p, xi, b, tR, tl, bs)
    L = R + 1/S
    r = (B*L - 1)/(B*L + 1)
    E = mp.exp(-2*A)
    return (1/B)*(1 + r*E)/(1 - r*E)


def jac_analytic(p, xi, b, tR, tl, bs, R=None):
    """Analytic Jacobian dZ/d(xi, b, tR, tl[, R]) from the log-derivative formulas.

    A d(log Z)/dA = -4 A r E/(1 - r^2 E^2)
    B d(log Z)/dB = -1 + 2 E (B dr/dB)/(1 - r^2 E^2),  B dr/dB = 2 B L/(B L + 1)^2
    d(log Z)/dR   = 2 E (dr/dR)/(1 - r^2 E^2),          dr/dR   = 2 B/(B L + 1)^2
    dA/dxi = A/xi, dB/db = B/b, dA/dtR = (A/2) h_R, dB/dtR = -(B/2) h_R, h_tau = p/(1 + tau p),
    dA/dtl = -(A/2) h_l, dB/dtl = (B/2) h_l.
    """
    p = mp.mpc(p)
    g, A, B, S = parts(p, xi, b, tR, tl, bs)
    RR = 0 if R is None else R
    L = RR + 1/S
    r = (B*L - 1)/(B*L + 1)
    E = mp.exp(-2*A)
    Z = (1/B)*(1 + r*E)/(1 - r*E)
    den = 1 - r**2*E**2
    u = Z*(-4*A*r*E/den)                      # A dZ/dA
    v = Z*(-1 + 2*E*(2*B*L/(B*L + 1)**2)/den)  # B dZ/dB
    hR = p/(1 + tR*p)
    hl = p/(1 + tl*p)
    W = u - v
    J = [u/xi, v/b, W*hR/2, -W*hl/2]
    if R is not None:
        J.append(Z*2*E*(2*B/(B*L + 1)**2)/den)
    return J


def jac_numeric(p, xi, b, tR, tl, bs, R=None):
    """Jacobian by high-precision numerical differentiation of the quadrupole form."""
    th = [xi, b, tR, tl] + ([] if R is None else [R])

    def f(*a):
        if R is None:
            return Z_quad(p, a[0], a[1], a[2], a[3], bs)
        return Z_quad(p, a[0], a[1], a[2], a[3], bs, a[4])
    J = []
    for i in range(len(th)):
        def fi(x, i=i):
            a = list(th)
            a[i] = x
            return f(*a)
        J.append(mp.diff(fi, th[i]))   # two-sided; the formula is analytic in each parameter for fixed p = i omega
    return J


def real_jacobian(omegas, theta, bs, R=None, scale=True):
    """Stack [Re J; Im J] over frequencies (p = i omega). Columns scaled by parameter value when > 0."""
    rows = []
    th = list(theta) + ([] if R is None else [R])
    for w in omegas:
        J = jac_analytic(mp.mpc(0, w), *theta, bs, R=R)
        Z = Z_quad(mp.mpc(0, w), *theta, bs, 0 if R is None else R)
        J = [Jk/Z for Jk in J]     # relative (log) sensitivities: d log Z / d theta
        if scale:
            J = [Jk*(th[k] if th[k] != 0 else 1) for k, Jk in enumerate(J)]
        rows.append([mp.re(x) for x in J])
        rows.append([mp.im(x) for x in J])
    return mp.matrix(rows)
