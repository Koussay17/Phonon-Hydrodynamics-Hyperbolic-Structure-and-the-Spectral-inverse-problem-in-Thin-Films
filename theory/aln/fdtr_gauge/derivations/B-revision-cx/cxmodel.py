"""Independent k = 0 model of a GK film with interface temperature theta = (1 + tau p) T.

Written from the definitions only (no import of the revision-branch scripts).

Film law (1D, k = 0):  C p T + q' = 0,  (1 + tau_R p) q = -lam T' + (1 + alpha) ell^2 q''.
Elimination gives q = -(lam/g) T',  g = (1 + tau_R p)/(1 + tau_l p), so the film is a quadrupole
    M = [[cosh A, sinh A / B], [B sinh A, cosh A]],  A = xi s sqrt(g),  B = b s / sqrt(g),  p = s^2.
Faces: theta = N T on the GK side, N = 1 + tau p (tau = c_n C).  Top face factor N0 = 1 + tau_o p,
bottom face factor N = 1 + tau_n p.  Jump theta_above - theta_below = R q.  Substrate: T/q = 1/(b_s s).

Bare film, read-out N0 T(0) (N0 = 1: read-out T(0)):      Z = N0 Mob(M; L/N),  L = R + 1/(b_s s).
Stack (known Fourier transducer xi_A, b_A above R1):      Z = Mob(T_A; R1 + N0 Mob(M; L2/N)).
Mob(M; x) = (M11 x + M12)/(M21 x + M22), evaluated as (1/B)(1 + r E)/(1 - r E),
r = (u - 1)/(u + 1), u = B x, E = exp(-2A) (stable for Re A > 0, true on p = i omega).

Units: b_s = 1 unless given; frequencies omega in the time unit of the taus; s = sqrt(i omega).

Backends: numpy (complex128, vectorised), mpmath (mpc), python-flint (acb).  Derivatives with
respect to the parameters are exact forward-mode (Dual numbers), not finite differences.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import mpmath as mp


# ----------------------------------------------------------------------------------------------
# forward-mode dual numbers over an arbitrary scalar/array backend
# ----------------------------------------------------------------------------------------------
class Dual:
    __array_ufunc__ = None          # numpy defers to the reflected operators
    __slots__ = ("v", "d")

    def __init__(self, v, d):
        self.v = v
        self.d = d                  # list of derivatives (same backend as v)

    def _wrap(self, o):
        return o if isinstance(o, Dual) else Dual(o, None)

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v + o.v, _addd(self.d, o.d))
        return Dual(self.v + o, self.d)
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, None if self.d is None else [-x for x in self.d])

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, Dual):
            d = _addd(_scal(self.d, o.v), _scal(o.d, self.v))
            return Dual(self.v * o.v, d)
        return Dual(self.v * o, _scal(self.d, o))
    __rmul__ = __mul__

    def inv(self):
        iv = 1 / self.v
        return Dual(iv, _scal(self.d, -iv * iv))

    def __truediv__(self, o):
        if isinstance(o, Dual):
            return self * o.inv()
        return Dual(self.v / o, _scal(self.d, 1 / o))

    def __rtruediv__(self, o):
        return self.inv() * o


def _addd(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return [x + y for x, y in zip(a, b)]


def _scal(a, c):
    if a is None:
        return None
    return [x * c for x in a]


class Ops:
    def __init__(self, sqrt, exp):
        self._sqrt, self._exp = sqrt, exp

    def sqrt(self, x):
        if isinstance(x, Dual):
            r = self._sqrt(x.v)
            return Dual(r, _scal(x.d, 1 / (2 * r)))
        return self._sqrt(x)

    def exp(self, x):
        if isinstance(x, Dual):
            e = self._exp(x.v)
            return Dual(e, _scal(x.d, e))
        return self._exp(x)


NP_OPS = Ops(np.sqrt, np.exp)
MP_OPS = Ops(mp.sqrt, mp.exp)


def flint_ops():
    return Ops(lambda z: z.sqrt(), lambda z: z.exp())


# ----------------------------------------------------------------------------------------------
# model specification
# ----------------------------------------------------------------------------------------------
class Spec:
    """Which parameters are unknown and how the face times depend on them.

    kind      'stack' or 'bare'
    conv      'tied'  : tau_o = eps0 * tl, tau_n = eps * tl   (twin: eps0 = eps for the stack)
              'free'  : tau_n is a parameter 'tn'; tau_o = tn (stack, same c_n at both faces)
                        or tau_o = eps0 * tl (bare film; eps0 = 0 is the T(0) read-out)
    eps0, eps face coefficients (tied); for 'free' bare film only eps0 is used
    R_known   bare film only: None -> R is a parameter; float -> fixed value
    transducer (xiA, bA) for the stack; bs substrate effusivity
    """

    def __init__(self, kind="stack", conv="tied", eps0=0.0, eps=0.0, R_known=None,
                 xiA=0.2, bA=3.0, bs=1.0):
        self.kind, self.conv = kind, conv
        self.eps0, self.eps = eps0, eps
        self.R_known = R_known
        self.xiA, self.bA, self.bs = xiA, bA, bs
        names = ["xi", "b", "tR", "tl"]
        if conv == "free":
            names.append("tn")
        if kind == "stack":
            names += ["R1", "R2"]
        elif R_known is None:
            names.append("R")
        self.names = names
        self.n = len(names)

    def label(self):
        if self.conv == "free":
            c = "free"
        elif self.eps0 == self.eps:
            c = f"eps={float(self.eps):g}"
        else:
            c = f"eps0={float(self.eps0):g},eps={float(self.eps):g}"
        R = "" if self.kind == "stack" else (" R unknown" if self.R_known is None else f" R={float(self.R_known):g}")
        return f"{self.kind}[{c}]{R}"

    def unpack(self, th):
        """Dict of physical quantities from the parameter vector (entries may be Duals)."""
        P = dict(zip(self.names, th))
        xi, b, tR, tl = P["xi"], P["b"], P["tR"], P["tl"]
        if self.conv == "free":
            tn = P["tn"]
            to = tn if self.kind == "stack" else self.eps0 * tl
        else:
            tn = self.eps * tl
            to = self.eps0 * tl
        if self.kind == "stack":
            R1, R2 = P["R1"], P["R2"]
        else:
            R1 = 0.0
            R2 = P["R"] if self.R_known is None else self.R_known
        return xi, b, tR, tl, to, tn, R1, R2


def response(spec: Spec, th, s, ops):
    """Z at s = sqrt(i omega) (scalar or array) for parameter vector th (plain or Dual entries)."""
    xi, b, tR, tl, to, tn, R1, R2 = spec.unpack(th)
    bs = spec.bs
    p = s * s
    g = (1 + tR * p) / (1 + tl * p)
    sg = ops.sqrt(g)
    A = xi * s * sg
    B = b * s / sg
    N = 1 + tn * p
    N0 = 1 + to * p
    L = R2 + 1 / (bs * s)
    u = B * L / N
    r = (u - 1) / (u + 1)
    E = ops.exp(-2 * A)
    rE = r * E
    Phi = N0 * (1 + rE) / ((1 - rE) * B)
    if spec.kind == "bare":
        return Phi
    W = R1 + Phi
    sA = spec.bA * s
    uA = sA * W
    rA = (uA - 1) / (uA + 1)
    EA = ops.exp(-2 * spec.xiA * s)
    return (1 + rA * EA) / ((1 - rA * EA) * sA)


def seed(th, backend_one, backend_zero):
    n = len(th)
    out = []
    for k, v in enumerate(th):
        d = [backend_zero] * n
        d[k] = backend_one
        out.append(Dual(v, d))
    return out


# ----------------------------------------------------------------------------------------------
# numpy front end (vectorised over frequency; optionally over a batch of parameter points)
# ----------------------------------------------------------------------------------------------
def s_of(omega):
    omega = np.asarray(omega, dtype=float)
    return np.sqrt(1j * omega)


def Z_np(spec, th, omega):
    return response(spec, [float(x) for x in th], s_of(omega), NP_OPS)


def ZJ_np(spec, th, omega):
    """(Z, dZ/dth) with dZ of shape (n, nf); exact derivatives."""
    s = s_of(omega)
    one = np.ones_like(s)
    zero = np.zeros_like(s)
    th_d = seed([float(x) for x in th], one, zero)
    z = response(spec, th_d, s, NP_OPS)
    return z.v, np.array(z.d)


def ZJ_np_batch(spec, TH, omega):
    """Batch: TH shape (m, n). Returns Z (m, nf) and J (m, n, nf)."""
    TH = np.asarray(TH, dtype=float)
    m, n = TH.shape
    s = s_of(omega)[None, :]
    one = np.ones((m, s.shape[1]), dtype=complex)
    zero = np.zeros((m, s.shape[1]), dtype=complex)
    th_d = []
    for k in range(n):
        d = [zero] * n
        d[k] = one
        th_d.append(Dual(TH[:, k:k + 1] * np.ones((1, s.shape[1])), d))
    z = response(spec, th_d, s, NP_OPS)
    return z.v, np.stack(z.d, axis=1)


def real_jac(Z, J, rel=True):
    """Real Jacobian rows [Re; Im] per frequency of (dZ/Z) (rel) or dZ; J shape (n, nf)."""
    G = J / Z[None, :] if rel else J
    return np.concatenate([G.real.T, G.imag.T], axis=0)


# ----------------------------------------------------------------------------------------------
# mpmath front end
# ----------------------------------------------------------------------------------------------
def Z_mp(spec, th, omega):
    th = [mp.mpf(x) for x in th]
    s = mp.sqrt(mp.mpc(0, omega))
    return response(spec, th, s, MP_OPS)


def ZJ_mp(spec, th, omega):
    th = [mp.mpf(x) for x in th]
    s = mp.sqrt(mp.mpc(0, omega))
    z = response(spec, seed(th, mp.mpc(1), mp.mpc(0)), s, MP_OPS)
    return z.v, list(z.d)


def realjac_mp(spec, th, omegas, rel=True):
    """mpmath real Jacobian (2 nf x n) of Z (or log Z if rel)."""
    rows_re, rows_im = [], []
    for w in omegas:
        z, d = ZJ_mp(spec, th, w)
        if rel:
            d = [x / z for x in d]
        rows_re.append([x.real for x in d])
        rows_im.append([x.imag for x in d])
    return mp.matrix(rows_re + rows_im)


def svd_mp(Jm, normalise=True):
    """Singular values (descending) of an mpmath matrix, columns normalised to unit norm."""
    J = Jm.copy()
    m, n = J.rows, J.cols
    scales = []
    if normalise:
        for k in range(n):
            nk = mp.sqrt(mp.fsum(J[i, k] ** 2 for i in range(m)))
            scales.append(nk)
            if nk != 0:
                for i in range(m):
                    J[i, k] /= nk
    S = mp.svd_r(J, compute_uv=False)
    sv = sorted([S[i] for i in range(len(S))], reverse=True)
    return sv, scales


def kernel_residual_mp(Jm, c):
    """|J c| / (|J| |c|) for a predicted kernel vector c (unnormalised columns)."""
    m, n = Jm.rows, Jm.cols
    c = [mp.mpf(x) for x in c]
    Jc = [mp.fsum(Jm[i, k] * c[k] for k in range(n)) for i in range(m)]
    nrm = mp.sqrt(mp.fsum(x ** 2 for x in Jc))
    nJ = mp.sqrt(mp.fsum(Jm[i, k] ** 2 for i in range(m) for k in range(n)))
    nc = mp.sqrt(mp.fsum(x ** 2 for x in c))
    return nrm / (nJ * nc)


# ----------------------------------------------------------------------------------------------
# physical (note-16) values for the Au/AlN/sapphire stack, nondimensionalised
# ----------------------------------------------------------------------------------------------
def aln_stack_nd(tau_R=25e-12, tau_l=10e-12, R1=1e-8, R2=2e-8, e=500e-9, lam=188.6, C=2.41e6,
                 eA=80e-9, lamA=150.0, CA=2.49e6, lams=35.0, Cs=3.03e6):
    """Nondimensional (xi, b, tR, tl, R1, R2), transducer (xiA, bA), time unit t0 = e^2/a."""
    a = lam / C
    t0 = e * e / a
    bs = np.sqrt(lams * Cs)
    b = np.sqrt(lam * C) / bs
    bA = np.sqrt(lamA * CA) / bs
    xiA = (eA / np.sqrt(lamA / CA)) / np.sqrt(t0)
    Rs = bs / np.sqrt(t0)
    th = (1.0, b, tau_R / t0, tau_l / t0, R1 * Rs, R2 * Rs)
    return th, (xiA, bA), t0
