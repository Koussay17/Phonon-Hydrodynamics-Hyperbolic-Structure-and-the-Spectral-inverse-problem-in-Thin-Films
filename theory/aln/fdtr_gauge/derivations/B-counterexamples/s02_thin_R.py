"""Thin-film expansion with interface resistance R (units xi1 = bs = 1 at the point; R in units of xi1/bs).
W = bs sqrt(p) Z, sqrt(p) -> eps sqrt(p), tau p -> eps^2 tau p, R fixed.  First-order coefficient and the
cumulative rank of the coefficient map in (xi1, beta, tR, tl, R) order by order."""
import sympy as sp
eps, q, xi, beta, tR, tl, R = sp.symbols("epsilon q xi beta tau_R tau_l R", positive=True)
r = sp.sqrt((1 + eps**2 * tR * q**2) / (1 + eps**2 * tl * q**2))
A = xi * eps * q * r
B = beta * eps * q / r                 # b sqrt(p) (1+tl p)^(1/2)/(1+tR p)^(1/2), bs = 1
Zs = R + 1 / (eps * q)
t = sp.tanh(A)
Z = (Zs + t / B) / (1 + B * Zs * t)
W = sp.expand(sp.series(eps * q * Z, eps, 0, 6).removeO())
params = [xi, beta, tR, tl, R]
pt = {xi: sp.Rational(7, 5), beta: sp.Rational(9, 4), tR: sp.Rational(3, 7), tl: sp.Rational(2, 9), R: sp.Rational(5, 3)}
coeffs = []
for k in range(1, 6):
    ck = sp.expand(W.coeff(eps, k))
    if k == 1:
        print("c1 =", sp.factor(ck), "  [expected q*(R + xi/beta - xi*beta)]")
    for cf in sp.Poly(ck, q).coeffs():
        coeffs.append(cf)
    J = sp.Matrix([[sp.diff(c, pp) for pp in params] for c in coeffs])
    print(f"order eps^{k}: cumulative rank {J.subs(pt).rank()}")
