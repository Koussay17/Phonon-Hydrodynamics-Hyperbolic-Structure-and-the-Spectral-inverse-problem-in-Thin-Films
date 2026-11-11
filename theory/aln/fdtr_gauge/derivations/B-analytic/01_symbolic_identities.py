"""Workstream B, analytic branch: symbolic identities used in proofs/B-analytic.md.

Each check is done symbolically (SymPy simplification) and, independently,
by evaluation at random rational points with 50-digit arithmetic.
Run:  python 01_symbolic_identities.py  > 01_symbolic_identities.out
"""
import random
import sympy as sp

random.seed(20261111)
xi, b, bs, tR, tl, p, R, t = sp.symbols('xi b b_s tau_R tau_l p R t', positive=True)
s = sp.symbols('s', positive=True)

g = (1 + tR*p)/(1 + tl*p)
A = xi*sp.sqrt(p*g)
B = b*sp.sqrt(p/g)
S = bs*sp.sqrt(p)
E = sp.exp(-2*A)

results = []


def check(name, expr, subs_extra=None, tol=sp.Rational(1, 10**35), symbolic=True, n=6):
    """expr must vanish identically. Test symbolically (if cheap) and at random points."""
    ok_sym = None
    if symbolic:
        try:
            ok_sym = sp.simplify(sp.expand(expr)) == 0
        except Exception:
            ok_sym = None
    worst = 0
    for _ in range(n):
        vals = {xi: sp.Rational(random.randint(1, 40), 10), b: sp.Rational(random.randint(1, 40), 10),
                bs: sp.Rational(random.randint(1, 40), 10), tR: sp.Rational(random.randint(0, 40), 10),
                tl: sp.Rational(random.randint(0, 40), 10), p: sp.Rational(random.randint(1, 60), 10),
                R: sp.Rational(random.randint(0, 40), 10), t: sp.Rational(random.randint(0, 40), 10),
                s: sp.Rational(random.randint(1, 30), 10)}
        for extra in ('q', 'k', 'cxi', 'cb', 'cR', 'cl'):
            if extra in globals():
                vals[globals()[extra]] = sp.Rational(random.randint(1, 30), 10)
        if subs_extra:
            vals.update({k: v.subs(vals) if hasattr(v, 'subs') else v for k, v in subs_extra.items()})
        v = sp.N(expr.subs(vals), 60)
        worst = max(worst, abs(complex(v)))
    ok_num = worst < 1e-35
    results.append((name, ok_sym, ok_num, worst))
    print(f"[{'OK' if ok_num else 'FAIL'}] {name}: symbolic={ok_sym}, max|residual| at random points = {worst:.3e}")
    return ok_num


print("=== 1. Front-face impedance from the quadrupole ===")
ch, sh = sp.cosh(A), sp.sinh(A)
M = sp.Matrix([[ch, sh/B], [B*sh, ch]])
Zs = 1/S
Z_quad = (M[0, 0]*Zs + M[0, 1])/(M[1, 0]*Zs + M[1, 1])
Z_true = (B*ch + S*sh)/(B*(B*sh + S*ch))
Z_task = (B*sh + S*ch)/(B*(B*ch + S*sh))      # formula as stated in the brief
check("Z_quad == (B cosh A + S sinh A)/(B (B sinh A + S cosh A))", Z_quad - Z_true)
check("Z_task * Z_true == 1/B^2  (the stated formula is the reciprocal form)", Z_task*Z_true - 1/B**2)
d = sp.N((Z_task - Z_true).subs({xi: 1, b: 2, bs: 1, tR: 1, tl: 2, p: 3}), 20)
print("   Z_task - Z_true at a sample point =", d, "(nonzero: the stated formula is not the impedance)")
# thin-film limit xi -> 0
print("   lim_{xi->0} Z_true =", sp.simplify(sp.limit(Z_true, xi, 0)), "   lim_{xi->0} Z_task =", sp.simplify(sp.limit(Z_task, xi, 0)))
Zr = (1/B)*(1 + (B - S)/(B + S)*E)/(1 - (B - S)/(B + S)*E)
check("Z_true == (1/B)(1 + r E)/(1 - r E), r=(B-S)/(B+S), E=exp(-2A)", Z_true - Zr, symbolic=False)

print("\n=== 2. Film quadrupole as a single-valued function of z = A^2 ===")
z = xi**2*p*g
check("A*B == xi*b*p", A*B - xi*b*p)
check("A/B == (xi/b) g", A/B - xi/b*g)
check("A^2 == z", A**2 - z)
Sc = sp.sinh(A)/A
check("sinh(A)/B == (xi/b) g sinh(A)/A", sh/B - xi/b*g*Sc)
check("B sinh(A) == xi b p sinh(A)/A", B*sh - xi*b*p*Sc)
C_, X_, Y_ = ch, xi/b*g*Sc, xi*b*p*Sc
check("det M = C^2 - X Y == 1", C_**2 - X_*Y_ - 1, symbolic=False)
# check that M entries are even in A (single-valued): replace A by -A in series sense
zz = sp.symbols('zz')
Cz = sp.cosh(sp.sqrt(zz)); Scz = sp.sinh(sp.sqrt(zz))/sp.sqrt(zz)
print("   cosh(sqrt z) series:", sp.series(Cz, zz, 0, 4))
print("   sinh(sqrt z)/sqrt z series:", sp.series(Scz, zz, 0, 4))
print("   (both are power series in z: entire functions of z, so M is single-valued in p)")

print("\n=== 3. Two-load identities on the s-plane (p = s^2, S = b_s s) ===")
Fp = (C_ + X_*S)/(Y_ + C_*S)
check("F - 1/S == X (S^2 - B^2)/(S (Y + C S))", Fp - 1/S - X_*(S**2 - B**2)/(S*(Y_ + C_*S)), symbolic=False)
# linear system rows for unknowns (C, X, Y):  row1 = (F S - 1, -S, F), row2 = (-Ft S - 1, S, Ft)
F, Ft, SS = sp.symbols('F Ft SS')
L = sp.Matrix([[F*SS - 1, -SS, F], [-Ft*SS - 1, SS, Ft]])
minors = [sp.factor(L[:, [i, j]].det()) for (i, j) in [(0, 1), (0, 2), (1, 2)]]
print("   2x2 minors of the system matrix:", minors)
sol = sp.solve([m for m in minors], [F, Ft], dict=True)
print("   common zeros of all minors (rank < 2):", sol)
# the true (C, X, Y) solves both rows when F = F(s), Ft = F(-s)
Fm = (C_ - X_*S)/(Y_ - C_*S)
row1 = (Fp*S - 1)*C_ - S*X_ + Fp*Y_
row2 = (-Fm*S - 1)*C_ + S*X_ + Fm*Y_
check("true (C,X,Y) satisfies row 1", row1, symbolic=False)
check("true (C,X,Y) satisfies row 2", row2, symbolic=False)
check("F(s)+F(-s) == 2 C X (S^2 - B^2)/(C^2 S^2 - Y^2)", Fp + Fm - 2*C_*X_*(S**2 - B**2)/(C_**2*S**2 - Y_**2), symbolic=False)

print("\n=== 4. Log-derivatives and the ratio u/v ===")
Aa, Bb, Ss = sp.symbols('Aa Bb Ss', positive=True)
Zab = (Bb*sp.cosh(Aa) + Ss*sp.sinh(Aa))/(Bb*(Bb*sp.sinh(Aa) + Ss*sp.cosh(Aa)))
Ee = sp.exp(-2*Aa); rr = (Bb - Ss)/(Bb + Ss)
u = Aa*sp.diff(Zab, Aa); v = Bb*sp.diff(Zab, Bb)
uform = Zab*(-4*Aa*rr*Ee/(1 - rr**2*Ee**2))
vform = Zab*(-1 + 4*Bb*Ss*Ee/((Bb + Ss)**2*(1 - rr**2*Ee**2)))


def num_check_ABS(name, expr, n=8):
    worst = 0
    for _ in range(n):
        vals = {Aa: sp.Rational(random.randint(1, 50), 10), Bb: sp.Rational(random.randint(1, 50), 10),
                Ss: sp.Rational(random.randint(1, 50), 10)}
        worst = max(worst, abs(complex(sp.N(expr.subs(vals), 60))))
    # complex test point too
    vals = {Aa: sp.Rational(7, 5) + sp.I*sp.Rational(9, 4), Bb: sp.Rational(3, 2) + sp.I*sp.Rational(1, 3), Ss: sp.Rational(1, 2) + sp.I*sp.Rational(1, 2)}
    worst = max(worst, abs(complex(sp.N(expr.subs(vals), 60))))
    ok = worst < 1e-35
    results.append((name, None, ok, worst))
    print(f"[{'OK' if ok else 'FAIL'}] {name}: max|residual| = {worst:.3e}")


num_check_ABS("A dZ/dA == Z * (-4 A r E/(1 - r^2 E^2))", u - uform)
num_check_ABS("B dZ/dB == Z * (-1 + 4 B S E/((B+S)^2 (1 - r^2 E^2)))", v - vform)
num_check_ABS("u/v == 4 A r E/((1 - E)(1 + r^2 E))", u/v - 4*Aa*rr*Ee/((1 - Ee)*(1 + rr**2*Ee)))
num_check_ABS("1 - r^2 == 4 B S/(B+S)^2", 1 - rr**2 - 4*Bb*Ss/(Bb + Ss)**2)

print("\n=== 5. Parameter Jacobian structure ===")
# J = (u/xi, v/b, W h_R/2, -W h_l/2), W = u - v, h_tau = p/(1+tau p)
Zfull = Z_true
JR = sp.diff(Zfull, tR); Jl = sp.diff(Zfull, tl)
check("dZ/dtau_R * (1+tau_R p) == - dZ/dtau_l * (1+tau_l p)", JR*(1 + tR*p) + Jl*(1 + tl*p), symbolic=False)
check("on tau_R = tau_l = t: dZ/dtau_R + dZ/dtau_l == 0", (JR + Jl).subs({tR: t, tl: t}), symbolic=False)

print("\n=== 6. Small-s (low-frequency) Laurent expansion, p = s^2 ===")
Zs_expr = Z_true.subs(p, s**2)
ser = sp.series(Zs_expr, s, 0, 3).removeO()
ser = sp.expand(ser)
cm1 = sp.simplify(ser.coeff(s, -1)); c0 = sp.simplify(ser.coeff(s, 0)); c1 = sp.simplify(ser.coeff(s, 1)); c2 = sp.simplify(ser.coeff(s, 2))
print("   c_-1 =", cm1)
print("   c_0  =", sp.factor(c0))
print("   c_1  =", sp.factor(c1))
print("   c_2  =", sp.factor(c2))
check("c_0 == xi (1/b - b/b_s^2)", c0 - xi*(1/b - b/bs**2))
check("c_1 == xi^2 (b^2 - b_s^2)/b_s^3", c1 - xi**2*(b**2 - bs**2)/bs**3)

print("\n=== 7. Doubly degenerate stratum D2 (tau_R = tau_l = t, b = b_s) ===")
cxi, cb, cR, cl = sp.symbols('c_xi c_b c_R c_l')
DZ = cxi*sp.diff(Zfull, xi) + cb*sp.diff(Zfull, b) + cR*sp.diff(Zfull, tR) + cl*sp.diff(Zfull, tl)
DZ_D2 = DZ.subs({tR: t, tl: t}).subs({b: bs})
pred = -(1 - sp.exp(-2*xi*sp.sqrt(p)))/(2*bs*sp.sqrt(p))*(2*cb/bs - (cR - cl)*p/(1 + t*p))
check("D_c Z on D2 == -(1-e^{-2A})/(2B) (2 c_b/b - (c_R - c_l) p/(1+t p))", DZ_D2 - pred, symbolic=False)

print("\n=== 8. Interface resistance R ===")
Lr = R + 1/S
Z_R = (C_*Lr + X_)/(Y_*Lr + C_)
rR = (B*Lr - 1)/(B*Lr + 1)
check("Z_R == (1/B)(1 + r_R E)/(1 - r_R E), r_R = (B L - 1)/(B L + 1), L = R + 1/S",
      Z_R - (1/B)*(1 + rR*E)/(1 - rR*E), symbolic=False)
check("dr_R/dR / r_R == 2B/(B^2 L^2 - 1)", sp.diff(rR, R)/rR - 2*B/(B**2*Lr**2 - 1), symbolic=False)
# Fourier film, R = 0: dZ/dR proportional to dZ/dxi with constant ratio
ZF = Z_R.subs({tR: 0, tl: 0})
ratio = (sp.diff(ZF, R)/sp.diff(ZF, xi)).subs(R, 0)
vals1 = {xi: sp.Rational(3, 2), b: sp.Rational(5, 2), bs: 1}
rv = [sp.N(ratio.subs(vals1).subs(p, pv), 30) for pv in [sp.Rational(1, 3), 2, 7, sp.Rational(31, 2)]]
print("   Fourier, R=0: (dZ/dR)/(dZ/dxi) at p = 1/3, 2, 7, 15.5:", [sp.nsimplify(x, tolerance=1e-25) for x in rv])
print("   predicted constant b b_s^2/(b_s^2 - b^2) at the same point:",
      sp.N((b*bs**2/(bs**2 - b**2)).subs(vals1), 30))
check("Fourier, R=0: dZ/dR - b b_s^2/(b_s^2 - b^2) dZ/dxi == 0  (thin-layer equivalence)",
      (sp.diff(ZF, R) - b*bs**2/(bs**2 - b**2)*sp.diff(ZF, xi)).subs(R, 0), symbolic=False)
# non-Fourier R = 0: ratio is not constant
ZNF = Z_R
ratio2 = (sp.diff(ZNF, R)/sp.diff(ZNF, xi)).subs(R, 0)
vals2 = {xi: sp.Rational(3, 2), b: sp.Rational(5, 2), bs: 1, tR: 1, tl: sp.Rational(1, 2)}
print("   non-Fourier (tau_R=1, tau_l=1/2), R=0: (dZ/dR)/(dZ/dxi) at p = 1/3, 2, 7:",
      [sp.N(ratio2.subs(vals2).subs(p, pv), 15) for pv in [sp.Rational(1, 3), 2, 7]])
# algebraic step: c_R ∂_R r = 2 c_xi r sqrt(p g)  <=>  b + k g = k b^2 (R q + 1/b_s)^2, q = sqrt(p), k = c_xi/c_R
q, k = sp.symbols('q k', positive=True)
gq = (1 + tR*q**2)/(1 + tl*q**2)
Bq = b*q/sp.sqrt(gq); Lq = R + 1/(bs*q)
lhs = 2*Bq/(Bq**2*Lq**2 - 1) - 2*k*q*sp.sqrt(gq)        # = 0  <=> relation
rhs_form = (b - k*(b**2*(R*q + 1/bs)**2 - gq))         # claimed equivalent numerator
check("2B/(B^2L^2-1) - 2k q sqrt(g) == 2 q sqrt(g) * [b - k(b^2(Rq+1/b_s)^2 - g)] / (b^2(Rq+1/b_s)^2 - g)",
      lhs - 2*q*sp.sqrt(gq)*rhs_form/(b**2*(R*q + 1/bs)**2 - gq), subs_extra=None, symbolic=False)

print("\nSUMMARY:", sum(1 for r_ in results if r_[2]), "of", len(results), "numerical checks passed")
