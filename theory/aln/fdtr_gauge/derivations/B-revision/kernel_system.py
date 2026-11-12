"""Polynomial kernel system of the 1D front-face map under general interface-temperature conventions.

Model family (k = 0, film Cattaneo/GK on a Fourier half-space, p Laplace variable, q = sqrt(p)):
    Z = N0 * Mobius(M; L/N),   M = [[cosh A, sinh A/B], [B sinh A, cosh A]],
    A^2 = xi^2 p g,  B^2 = b^2 p/g,  g = (1 + tR p)/(1 + tl p),  L = R + 1/(bs q),
    N = 1 + tau_n p   (theta-jump at the film/substrate face, tau_n = c_n C),
    N0 = 1 + tau_o p  (temperature read at the film top: T -> tau_o = 0; theta -> tau_o = tau of that face).
Equivalently Z = (1/Bt)(1 + r E)/(1 - r E), Bt = B/N0, u = B L/N, r = (u - 1)/(u + 1), E = exp(-2A).

LEMMA K (proofs/B-revision.md, Sec. 3): for a real direction c of the parameters, D_c Z == 0 as a function iff
   (i)  c_b = 0  and  D_c log g = -2 dtau_o p/N0,
   (O)  N dtau_o - N0 dtau_n = 2 b R kappa,           kappa = c_xi N0 - xi dtau_o p,
   (E)  c_R N N0 = kappa (b/bs^2 - b R^2 p - g N^2/b),
where dtau_o = D_c tau_o, dtau_n = D_c tau_n, c_R = D_c R. All three are polynomial identities in p after
clearing the denominators (1 + tR p)(1 + tl p)(1 + tau_o p) and b (1 + tl p); their p-coefficients give a linear
system Mmat(theta) c = 0, and rank D Z(theta) = n - dim ker Mmat(theta).

This module builds Mmat for a model specification. Pure SymPy, exact.
"""
import itertools
import sympy as sp

p = sp.Symbol("p")
xi, b, bs = sp.symbols("xi b b_s", positive=True)
tR, tl = sp.symbols("tau_R tau_l", nonnegative=True)
R = sp.Symbol("R", nonnegative=True)
eps, eps0 = sp.symbols("epsilon epsilon_0", real=True)
tn_s, to_s = sp.symbols("tau_n tau_o", real=True)

C_XI, C_B, C_TR, C_TL, C_TN, C_TO, C_R = sp.symbols("c_xi c_b c_tR c_tl c_tn c_to c_R")


def model(obs="T", iface="local", R_unknown=False):
    """Return dict with tau_o, tau_n (values), dtau_o, dtau_n (directional derivatives), the list of
    direction components, and labels.

    obs   : 'T' (tau_o = 0), 'theta_e' (tau_o = tau_l), 'theta_q' (tau_o = tau_n: theta-jumps at both faces),
            'eps0' (tau_o = eps0 tau_l), 'free_o' (tau_o a free parameter).
    iface : 'local' (tau_n = 0), 'eps' (tau_n = eps tau_l), 'free' (tau_n a free parameter).
    """
    comps = [C_XI, C_B, C_TR, C_TL]
    names = ["xi", "b", "tau_R", "tau_l"]
    if iface == "local":
        tau_n, d_n = sp.Integer(0), sp.Integer(0)
    elif iface == "eps":
        tau_n, d_n = eps * tl, eps * C_TL
    elif iface == "free":
        tau_n, d_n = tn_s, C_TN
        comps.append(C_TN); names.append("tau_n")
    else:
        raise ValueError(iface)
    if obs == "T":
        tau_o, d_o = sp.Integer(0), sp.Integer(0)
    elif obs == "theta_e":
        tau_o, d_o = tl, C_TL
    elif obs == "theta_q":
        tau_o, d_o = tau_n, d_n
    elif obs == "eps0":
        tau_o, d_o = eps0 * tl, eps0 * C_TL
    elif obs == "free_o":
        tau_o, d_o = to_s, C_TO
        comps.append(C_TO); names.append("tau_o")
    else:
        raise ValueError(obs)
    cR = sp.Integer(0)
    if R_unknown:
        comps.append(C_R); names.append("R")
        cR = C_R
    return dict(obs=obs, iface=iface, R_unknown=R_unknown, tau_o=tau_o, tau_n=tau_n, d_o=d_o, d_n=d_n,
                cR=cR, comps=comps, names=names)


def identities(m):
    """The three polynomial identities (in p) of Lemma K, denominators cleared."""
    N = 1 + m["tau_n"] * p
    N0 = 1 + m["tau_o"] * p
    g_num, g_den = 1 + tR * p, 1 + tl * p
    # (i): D_c log(Bt^2) = 2 c_b/b - [c_tR p/(1+tR p) - c_tl p/(1+tl p)] - 2 d_o p/N0 == 0
    I1 = sp.expand(2 * C_B / b * g_num * g_den * N0 - C_TR * p * g_den * N0 + C_TL * p * g_num * N0
                   - 2 * m["d_o"] * p * g_num * g_den)
    kappa = C_XI * N0 - xi * m["d_o"] * p
    I2 = sp.expand(N * m["d_o"] - N0 * m["d_n"] - 2 * b * R * kappa)
    # (E) times b (1 + tl p):  c_R N N0 b g_den - kappa (b^2 g_den/bs^2 - b^2 R^2 p g_den - g_num N^2)
    I3 = sp.expand(m["cR"] * N * N0 * b * g_den - kappa * (b ** 2 * g_den / bs ** 2 - b ** 2 * R ** 2 * p * g_den
                                                         - g_num * N ** 2))
    return I1, I2, I3


def matrix(m):
    """Coefficient matrix Mmat (rows: p-coefficients of the identities; columns: direction components)."""
    rows = []
    for I in identities(m):
        P = sp.Poly(I, p)
        for coeff in P.all_coeffs():
            row = [sp.expand(sp.diff(coeff, c)) for c in m["comps"]]
            # the identities are linear and homogeneous in the components
            assert sp.simplify(coeff - sum(r * c for r, c in zip(row, m["comps"]))) == 0
            if any(r != 0 for r in row):
                rows.append(row)
    return sp.Matrix(rows)


def maximal_minors(Mm):
    n = Mm.shape[1]
    out = []
    for idx in itertools.combinations(range(Mm.shape[0]), n):
        d = sp.factor(Mm.extract(list(idx), list(range(n))).det(method="berkowitz"))
        if d != 0:
            out.append(d)
    return out


def kernel_at(m, subs):
    """Exact kernel of Mmat at a (possibly partially symbolic) substitution."""
    Mm = matrix(m).subs(subs)
    Mm = Mm.applyfunc(sp.simplify)
    return Mm.rank(simplify=True), [sp.simplify(v) for v in Mm.nullspace(simplify=True)]
