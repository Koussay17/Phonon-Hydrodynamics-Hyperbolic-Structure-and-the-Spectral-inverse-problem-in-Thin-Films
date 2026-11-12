"""Krawczyk certification of twin roots with python-flint ball arithmetic (independent implementation).

F(th) = [Re; Im](Z(omega_j; th)/Z(omega_j; th*) - 1), j = 1..N (square: N = 3 for the 6-parameter stack;
for joint pairs the unknowns are (th, th') with some coordinates fixed).
K(X) = m - Y F(m) + (I - Y J(X)) (X - m);  K(X) inside int(X)  =>  exactly one root in X, and it is regular.
"""
import flint
from flint import arb, acb, arb_mat
import mpmath as mp
import cxmodel as cx

OPS = cx.flint_ops()


def _acb_par(x, rad=None):
    a = arb(x) if rad is None else arb(x, rad)
    return acb(a)


def Z_flint(spec, th, omega, dual=False):
    s = acb(0, arb(omega)).sqrt()
    if dual:
        th_d = cx.seed(th, acb(1), acb(0))
        return cx.response(spec, th_d, s, OPS)
    return cx.response(spec, th, s, OPS)


def krawczyk_target(spec, target_str, omega_str, mid_str, rad="1e-40", prec=400):
    """Certify a root of the target-fixed system. Returns dict with containment margin."""
    flint.ctx.prec = prec
    n = spec.n
    tg = [_acb_par(x) for x in target_str]
    Zt = [Z_flint(spec, tg, w) for w in omega_str]
    m = [arb(x) for x in mid_str]
    r = arb(rad)
    # F at the midpoint
    Fm = []
    for w, zt in zip(omega_str, Zt):
        z = Z_flint(spec, [acb(x) for x in m], w)
        q = z / zt - 1
        Fm.append(q)
    Fm = [q.real for q in Fm] + [q.imag for q in Fm]
    # Jacobian over the box
    X = [acb(arb(mid_str[k], rad)) for k in range(n)]
    rows_re, rows_im = [], []
    for w, zt in zip(omega_str, Zt):
        zd = Z_flint(spec, X, w, dual=True)
        rows_re.append([(d / zt).real for d in zd.d])
        rows_im.append([(d / zt).imag for d in zd.d])
    JX = arb_mat(rows_re + rows_im)
    # Y: inverse of the midpoint Jacobian (floating approximation, then exact as arb)
    Jmid = arb_mat([[arb(JX[i, k].mid()) for k in range(n)] for i in range(n)])
    Y = Jmid.inv()
    Y = arb_mat([[arb(Y[i, k].mid()) for k in range(n)] for i in range(n)])
    I = arb_mat([[arb(1) if i == k else arb(0) for k in range(n)] for i in range(n)])
    Fv = arb_mat([[x] for x in Fm])
    Xm = arb_mat([[arb(0, rad)] for k in range(n)])
    K = arb_mat([[x] for x in m]) - Y * Fv + (I - Y * JX) * Xm
    inside = True
    margins = []
    for k in range(n):
        dev = abs(K[k, 0] - m[k])
        up = dev.upper()          # arb upper bound of |K_k - m_k| (including radius)
        margins.append(float(up / r.mid()))
        if not (up < r):
            inside = False
    return dict(certified=bool(inside), max_K_dev_over_radius=max(margins),
                F_mid_abs_max=float(max(abs(x).upper() for x in Fm)))
