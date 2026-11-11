"""Thin-film / low-frequency hierarchy of identifiable combinations (symbolic).

Scaling: sqrt(p) -> eps sqrt(p) (so p -> eps^2 p, tau p -> eps^2 tau p).  W = bs sqrt(p) Z expanded in eps.
The k-th order coefficient c_k(p) is a polynomial in sqrt(p); its monomials give real functions of
theta.  Rank of the Jacobian of {all monomial coefficients up to order k} w.r.t. (xi1, beta, tR, tl):
predicts singular values s_j ~ (omega xi1^2)^{k_j/2} where k_j = first order at which the rank reaches j.
Numerical fits (practical_bands.py) gave exponents 1/2, 1, 3/2, 5/2 for s1..s4, i.e. k = 1, 2, 3, 5.
Also the matched-effusivity case beta = 1.
"""
import sympy as sp

eps, q, xi, beta, tR, tl = sp.symbols("epsilon q xi beta tau_R tau_l", positive=True)   # q = sqrt(p)
N = 6
r2 = (1 + eps ** 2 * tR * q ** 2) / (1 + eps ** 2 * tl * q ** 2)
r = sp.sqrt(r2)
A = xi * eps * q * r
y = beta / r
t = sp.tanh(A)
W = (1 / y) * (y + t) / (1 + y * t)


def hierarchy(Wexpr, params, label, point):
    ser = sp.series(Wexpr, eps, 0, N + 1).removeO()
    ser = sp.expand(ser)
    coeffs = []
    ranks = []
    for k in range(1, N + 1):
        ck = sp.expand(ser.coeff(eps, k))
        poly = sp.Poly(ck, q)
        for mon, cf in zip(poly.monoms(), poly.coeffs()):
            coeffs.append(cf)
        J = sp.Matrix([[sp.diff(cf, pp) for pp in params] for cf in coeffs])
        ranks.append(J.subs(point).rank())
        print(f"{label}: order eps^{k}: c_{k} = {sp.factor(ck.subs(beta, point[beta]) if label.startswith('beta = 1') else ck)}"
              f"   cumulative rank = {ranks[-1]}")
    return ranks


generic = {xi: sp.Rational(7, 5), beta: sp.Rational(9, 4), tR: sp.Rational(3, 7), tl: sp.Rational(2, 9)}
print("generic beta:")
hierarchy(W, [xi, beta, tR, tl], "beta free", generic)
print()
print("matched effusivity beta = 1 (derivatives in all four parameters, evaluated at beta = 1):")
matched = dict(generic); matched[beta] = 1
hierarchy(W, [xi, beta, tR, tl], "beta = 1", matched)
