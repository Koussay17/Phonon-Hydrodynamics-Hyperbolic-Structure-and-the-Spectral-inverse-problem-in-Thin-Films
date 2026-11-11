"""Audit of proofs/B-analytic.md: remaining algebraic steps.

(1) Theorem 1 Step 2: 2x2 minors of the two-load matrix and their common zeros.
(2) Theorem 1 Step 5: (1+tR p)(1+tl' p) == (1+tR' p)(1+tl p)  <=>  tR - tl = tR' - tl' and delta (tl' - tl) = 0.
(3) Theorem 2 Proof A, D2 case: D_c Z at b = b_s, tau_R = tau_l = t, by numerical differentiation (mpmath),
    compared with -(1 - e^{-2A})/(2B) (2 c_b/b - (c_R - c_l) p/(1 + t p)).
(4) Theorem G: invariance of Phi under g_mu, kernel of D Phi, rank of D Phi, transversality of the constraints.
(5) Theorem 2 Proof A: d/dw cosh(sqrt w) = S(w)/2 (series identity).
Run:  python 03_gauge_misc.py > 03_gauge_misc.out
"""
import sympy as sp
import mpmath as mp

print("=== (1) two-load minors ===")
F, Ft, S = sp.symbols('F Ft S')
L = sp.Matrix([[F*S - 1, -S, F], [-Ft*S - 1, S, Ft]])
ms = [sp.factor(L.extract([0, 1], [i, j]).det()) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
print("minors:", ms)
print("common zeros (S != 0):", sp.solve(ms, [F, Ft], dict=True))

print("\n=== (2) rational identity of Step 5 ===")
p, a1, b1, a2, b2 = sp.symbols('p tR tl tRp tlp')
poly = sp.Poly(sp.expand((1 + a1*p)*(1 + b2*p) - (1 + a2*p)*(1 + b1*p)), p)
eqs = poly.coeffs()
print("coefficient equations:", eqs)
sol = sp.solve(eqs, [a2, b2], dict=True)
print("solutions for (tau_R', tau_l'):", sol)

print("\n=== (3) D2 derivative formula (numerical differentiation) ===")
mp.mp.dps = 40


def Z(p, xi, b, tR, tl, bs=1):
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*mp.sqrt(p*g); B = b*mp.sqrt(p/g); Sb = bs*mp.sqrt(p)
    return (B*mp.cosh(A) + Sb*mp.sinh(A))/(B*(B*mp.sinh(A) + Sb*mp.cosh(A)))


worst = 0
for (xi, t, w) in [(mp.mpf('0.7'), mp.mpf('0.4'), mp.mpf('0.3')), (mp.mpf('1.9'), mp.mpf(0), mp.mpf(5)),
                   (mp.mpf('1.2'), mp.mpf('2.5'), mp.mpf('1.1'))]:
    p = mp.mpc(0, w)
    th = [xi, mp.mpf(1), t, t]
    grads = []
    for k in range(4):
        def f(x, k=k):
            a = list(th); a[k] = x
            return Z(p, *a)
        grads.append(mp.diff(f, th[k]))
    for c in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (mp.mpf('0.3'), -2, 1, mp.mpf('0.5'))]:
        lhs = sum(ci*gi for ci, gi in zip(c, grads))
        A = xi*mp.sqrt(p); B = mp.sqrt(p)
        rhs = -(1 - mp.exp(-2*A))/(2*B)*(2*c[1]/1 - (c[2] - c[3])*p/(1 + t*p))
        worst = max(worst, abs(lhs - rhs)/max(abs(rhs), mp.mpf(10)**-30))
print("max relative discrepancy of the D2 formula:", mp.nstr(worst, 3))

print("\n=== (4) Theorem G ===")
e, lamb, C, Lam, tR, mu = sp.symbols('e lambda C Lambda tau_R mu', positive=True)
Phi = sp.Matrix([e*sp.sqrt(C/lamb), sp.sqrt(lamb*C), tR, Lam*C/lamb])
gmu = {e: mu*e, lamb: mu*lamb, C: C/mu, Lam: mu**2*Lam}
print("Phi(g_mu P) - Phi(P) =", sp.simplify(Phi.subs(gmu, simultaneous=True) - Phi).T)
vars_ = [e, lamb, C, Lam, tR]
DPhi = Phi.jacobian(vars_)
v = sp.Matrix([e, lamb, -C, 2*Lam, 0])
print("D Phi . (e, lambda, -C, 2 Lambda, 0) =", sp.simplify(DPhi*v).T)
print("rank D Phi (generic) =", DPhi.rank(), "; at Lambda = 0:", DPhi.subs(Lam, 0).rank())
# transversality: derivative of each constraint along the orbit generator
a = lamb/C
v0 = sp.symbols('v0', positive=True)
for name, cons in [("e", e), ("lambda", lamb), ("C", C), ("a", a), ("Lambda", Lam), ("a - v0^2 tau_R/3", a - v0**2*tR/3)]:
    d = sp.simplify(sum(sp.diff(cons, x)*vx for x, vx in zip(vars_, v)))
    print(f"   derivative of constraint {name} along the orbit: {d}")

print("\n=== (5) d/dw cosh(sqrt w) = S(w)/2 as power series ===")
w = sp.symbols('w')
Cs = sum(w**n/sp.factorial(2*n) for n in range(12))
Ss = sum(w**n/sp.factorial(2*n + 1) for n in range(12))
print("difference up to w^10:", sp.expand(sp.diff(Cs, w) - Ss/2 + sp.Rational(1, 2)*w**11/sp.factorial(23)))
