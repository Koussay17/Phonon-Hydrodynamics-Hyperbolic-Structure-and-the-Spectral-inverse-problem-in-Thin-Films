"""Front-face response of the 1D film problem under general interface-temperature conventions (mpmath).

Nondimensional units: b_s = 1, time unit fixed by the frequencies (omega xi^2 in the examples).
    g = (1 + tR p)/(1 + tl p),  A = xi sqrt(p g),  B = b sqrt(p/g),  S = sqrt(p)  (principal roots)
    M = [[cosh A, sinh A/B], [B sinh A, cosh A]]                         (film quadrupole, det 1)
    N = 1 + tau_n p  (theta-jump at the film/substrate face),  N0 = 1 + tau_o p (read-out / top-face coefficient)
    bare film:      Z = N0 * Mobius(M; (R + 1/S)/N)
    transducer:     Z_tot = Mobius(T_A; R1 + N0 * Mobius(M; (R2 + 1/S)/N)),  T_A = Fourier layer (xi_A, b_A) known
Conventions: tau_n = eps tau_l (local eps = 0, entropic eps = 1, kinetic eps real) or free;
             tau_o = 0 (T read at the heated face), tau_l (theta_e), tau_n (theta-quadrupole), or free.
The function `response(par, w, spec)` takes a parameter dict and a frequency omega (p = i omega).
"""
import mpmath as mp


def mobius(M, x):
    (a, b_), (c, d) = M
    return (a * x + b_) / (c * x + d)


def film_quad(xi, b, tR, tl, p):
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * mp.sqrt(p * g)
    B = b * mp.sqrt(p / g)
    return ((mp.cosh(A), mp.sinh(A) / B), (B * mp.sinh(A), mp.cosh(A)))


def taus(par, spec):
    tl = par["tl"]
    iface = spec.get("iface", "local")
    if iface == "local":
        tn = 0
    elif iface == "eps":
        tn = spec["eps"] * tl
    elif iface == "free":
        tn = par["tn"]
    else:
        raise ValueError(iface)
    obs = spec.get("obs", "T")
    if obs == "T":
        to = 0
    elif obs == "theta_e":
        to = tl
    elif obs == "theta_q":
        to = tn
    elif obs == "eps0":
        to = spec["eps0"] * tl
    elif obs == "free_o":
        to = par["to"]
    else:
        raise ValueError(obs)
    return tn, to


def bare(par, p, spec):
    tn, to = taus(par, spec)
    S = mp.sqrt(p)
    L = par.get("R", 0) + 1 / S
    Mf = film_quad(par["xi"], par["b"], par["tR"], par["tl"], p)
    return (1 + to * p) * mobius(Mf, L / (1 + tn * p))


def response(par, w, spec):
    p = mp.mpc(0, w)
    if spec.get("transducer"):
        xa, ba = spec["xiA"], spec["bA"]
        par2 = dict(par)
        par2["R"] = par["R2"]
        Win = par["R1"] + bare(par2, p, spec)
        A = xa * mp.sqrt(p)
        B = ba * mp.sqrt(p)
        return mobius(((mp.cosh(A), mp.sinh(A) / B), (B * mp.sinh(A), mp.cosh(A))), Win)
    return bare(par, p, spec)


def jacobian(par, names, ws, spec, h=None, log_response=True):
    """Real Jacobian [Re; Im] of log Z (or Z) w.r.t. the parameters `names`, central differences.
    Each column is the derivative w.r.t. the raw parameter (no log scaling)."""
    if h is None:
        h = mp.mpf(10) ** (-(mp.mp.dps // 3))
    f = (lambda q, w: mp.log(response(q, w, spec))) if log_response else (lambda q, w: response(q, w, spec))
    cols = []
    for nm in names:
        qp, qm = dict(par), dict(par)
        qp[nm] = par[nm] + h
        qm[nm] = par[nm] - h
        cols.append([(f(qp, w) - f(qm, w)) / (2 * h) for w in ws])
    J = mp.matrix(2 * len(ws), len(names))
    for j, col in enumerate(cols):
        for i, d in enumerate(col):
            J[2 * i, j] = mp.re(d)
            J[2 * i + 1, j] = mp.im(d)
    return J


def normalized_svals(J, zero_tol=None):
    """Singular values (relative to the largest) of the column-normalised matrix.
    A column whose norm is below zero_tol times the largest column norm is identically zero to working accuracy
    (e.g. d/dxi where the response does not depend on xi); it is set to zero instead of being normalised, so that it
    shows as a zero singular value rather than as amplified rounding noise."""
    if zero_tol is None:
        zero_tol = mp.mpf(10) ** (-(mp.mp.dps // 2))
    Jn = J.copy()
    norms = [mp.sqrt(sum(Jn[i, j] ** 2 for i in range(Jn.rows))) for j in range(Jn.cols)]
    big = max(norms)
    for j in range(Jn.cols):
        nrm = norms[j]
        for i in range(Jn.rows):
            Jn[i, j] = Jn[i, j] / nrm if nrm > zero_tol * big else mp.mpf(0)
    s = sorted([abs(x) for x in mp.svd_r(Jn, compute_uv=False)], reverse=True)
    return [x / s[0] for x in s]


def kernel_residual(J, c):
    """|J c| / (|J|_F |c|) for a predicted kernel vector c."""
    v = mp.matrix(c)
    r = J * v
    nr = mp.sqrt(sum(r[i] ** 2 for i in range(r.rows)))
    nJ = mp.sqrt(sum(J[i, j] ** 2 for i in range(J.rows) for j in range(J.cols)))
    nc = mp.sqrt(sum(x ** 2 for x in c))
    return nr / (nJ * nc)
