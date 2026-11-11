"""Audit of proofs/B-analytic.md: independent derivation of the low-frequency Laurent coefficients
and of the moment minors used in Theorem 2 (Proof C), Theorem N2 and section 7.

Route (independent of the closed-form quadrupole and of laurent.py / 07_moment_jacobian.py):
the impedance W(xi) = T/phi seen at Liouville depth xi = z/sqrt(a) obeys the Riccati equation
    dW/dxi = b p W^2 - g/b,      W(xi_1) = L = R + 1/(b_s s),      Z = W(0),
which follows from  d_xi phi = -b p T,  d_xi T = -(g/b) phi  (film law of note 7).
With p = s^2 and W = 1/(b_s s) + V:
    V' = b/b_s^2 - g/b + (2 b s/b_s) V + b s^2 V^2,   V(xi_1) = R,
so V = sum_k v_k(xi) s^k with polynomial v_k obtained order by order, and c_k = v_k(0).
Run:  python 01_riccati_moments.py > 01_riccati_moments.out
"""
from itertools import combinations
import sympy as sp

xi, b, bs, tR, tl, R, eta, lam, s = sp.symbols('xi b b_s tau_R tau_l R eta lam s', positive=True)
NK = 8                                   # c_0 .. c_7


def g_coeffs(n):
    """Coefficients g_k of g(s^2) = (1 + tR s^2)/(1 + tl s^2) in powers of s (k < n)."""
    ser = sp.series((1 + tR*s**2)/(1 + tl*s**2), s, 0, n).removeO()
    return [sp.expand(ser.coeff(s, k)) for k in range(n)]


def riccati_moments(n, with_R):
    gk = g_coeffs(n)
    v = []
    for k in range(n):
        rhs = -gk[k]/b
        if k == 0:
            rhs += b/bs**2
        if k >= 1:
            rhs += 2*b/bs*v[k - 1]
        if k >= 2:
            rhs += b*sum(v[i]*v[k - 2 - i] for i in range(k - 1))
        P = sp.integrate(sp.expand(rhs), eta)           # antiderivative in eta
        vk = (R if (k == 0 and with_R) else 0) - (P.subs(eta, xi) - P)
        v.append(sp.expand(vk))
    return [sp.factor(sp.expand(vk.subs(eta, 0))) for vk in v]


# ---------------------------------------------------------------------------------------------
print("=== (1) Laurent coefficients c_0..c_7 by the Riccati route (R = 0) ===")
c = riccati_moments(NK, with_R=False)
for k, ck in enumerate(c):
    print(f"c_{k} =", ck)

# document values (07_moment_jacobian.out, copied verbatim)
doc = {
    0: "-xi*(b - b_s)*(b + b_s)/(b*b_s**2)",
    1: "xi**2*(b - b_s)*(b + b_s)/b_s**3",
    2: "-xi*(3*b**4*xi**2 - 4*b**2*b_s**2*xi**2 - 3*b_s**4*tau_R + 3*b_s**4*tau_l + b_s**4*xi**2)/(3*b*b_s**4)",
    3: "xi**2*(3*b**4*xi**2 - 5*b**2*b_s**2*xi**2 - 3*b_s**4*tau_R + 3*b_s**4*tau_l + 2*b_s**4*xi**2)/(3*b_s**5)",
    4: "-xi*(15*b**6*xi**4 - 30*b**4*b_s**2*xi**4 - 20*b**2*b_s**4*tau_R*xi**2 + 20*b**2*b_s**4*tau_l*xi**2 + 17*b**2*b_s**4*xi**4 + 15*b_s**6*tau_R*tau_l + 10*b_s**6*tau_R*xi**2 - 15*b_s**6*tau_l**2 - 10*b_s**6*tau_l*xi**2 - 2*b_s**6*xi**4)/(15*b*b_s**6)",
    5: "xi**2*(45*b**6*xi**4 - 105*b**4*b_s**2*xi**4 - 75*b**2*b_s**4*tau_R*xi**2 + 75*b**2*b_s**4*tau_l*xi**2 + 77*b**2*b_s**4*xi**4 + 45*b_s**6*tau_R*tau_l + 60*b_s**6*tau_R*xi**2 - 45*b_s**6*tau_l**2 - 60*b_s**6*tau_l*xi**2 - 17*b_s**6*xi**4)/(45*b_s**7)",
    6: "-xi*(315*b**8*xi**6 - 840*b**6*b_s**2*xi**6 - 630*b**4*b_s**4*tau_R*xi**4 + 630*b**4*b_s**4*tau_l*xi**4 + 756*b**4*b_s**4*xi**6 + 420*b**2*b_s**6*tau_R*tau_l*xi**2 + 714*b**2*b_s**6*tau_R*xi**4 - 420*b**2*b_s**6*tau_l**2*xi**2 - 714*b**2*b_s**6*tau_l*xi**4 - 248*b**2*b_s**6*xi**6 + 105*b_s**8*tau_R**2*xi**2 - 315*b_s**8*tau_R*tau_l**2 - 420*b_s**8*tau_R*tau_l*xi**2 - 126*b_s**8*tau_R*xi**4 + 315*b_s**8*tau_l**3 + 315*b_s**8*tau_l**2*xi**2 + 126*b_s**8*tau_l*xi**4 + 17*b_s**8*xi**6)/(315*b*b_s**8)",
}
loc = {'xi': xi, 'b': b, 'b_s': bs, 'tau_R': tR, 'tau_l': tl}
print("\n  agreement with the document's c_0..c_6:")
for k, ex in doc.items():
    diff = sp.simplify(c[k] - sp.sympify(ex, locals=loc))
    print(f"   c_{k}: difference = {diff}")

# ---------------------------------------------------------------------------------------------
print("\n=== (2) moment minors G_K (4 parameters) ===")
params = [xi, b, tR, tl]
G = sp.Matrix([[sp.diff(ck, q) for q in params] for ck in c])


def minor(K, M=G, ncol=4):
    return sp.factor(sp.cancel(M.extract(list(K), list(range(ncol))).det()))


G0123 = minor((0, 1, 2, 3))
G0124 = minor((0, 1, 2, 4))
G0234 = minor((0, 2, 3, 4))
print("G_0123 =", G0123)
print("G_0124 =", G0124)
print("  document: 2 xi^4 (b^2-b_s^2)(tau_l-tau_R)/(b^4 b_s^3); difference =",
      sp.simplify(G0124 - 2*xi**4*(b**2 - bs**2)*(tl - tR)/(b**4*bs**3)))
print("G_0234 =", G0234)
G0234bs = sp.factor(sp.cancel(G0234.subs(b, bs)))
print("G_0234 at b = b_s =", G0234bs, "; document: -2 xi^4 (tau_R-tau_l)^2/b_s^5; difference =",
      sp.simplify(G0234bs + 2*xi**4*(tR - tl)**2/bs**5))
print("grad c_0 at b = b_s:", [sp.simplify(x.subs(b, bs)) for x in G.row(0)])
print("grad c_1 at b = b_s:", [sp.simplify(x.subs(b, bs)) for x in G.row(1)])
print("all 4-subsets K with |K| <= 9 and their minors at b = b_s:")
for K in combinations(range(NK), 4):
    if sum(K) <= 9:
        print(f"   K = {K}, |K| = {sum(K)}: G_K(b=b_s) = {sp.factor(sp.cancel(minor(K).subs(b, bs)))}")
print("all 4-subsets K with |K| <= 7 and their generic minors:")
for K in combinations(range(NK), 4):
    if sum(K) <= 7:
        print(f"   K = {K}, |K| = {sum(K)}: G_K = {minor(K)}")
print("every 4x4 minor with |K| <= 12 vanishes on tau_R = tau_l:",
      all(sp.cancel(minor(K).subs(tR, tl)) == 0 for K in combinations(range(NK), 4) if sum(K) <= 12))
# rank of the moment matrix on the resonance (expected 3 off D2, 2 on D2)
for lab, sub in [("Delta, b != b_s", {tR: tl}), ("D2", {tR: tl, b: bs})]:
    Mn = G.subs(sub).subs({xi: sp.Rational(3, 2), b: sp.Rational(5, 2), bs: 1, tl: sp.Rational(2, 7)}) \
        if lab != "D2" else G.subs(sub).subs({xi: sp.Rational(3, 2), bs: 1, tl: sp.Rational(2, 7)})
    print(f"   rank of [grad c_0..c_7] on {lab} at a rational point: {Mn.rank()}")

# ---------------------------------------------------------------------------------------------
print("\n=== (3) Cauchy-Binet factors t_K(lam) ===")


def tK(K):
    rows = [[sp.cos(k*sp.pi/4) for k in K], [sp.sin(k*sp.pi/4) for k in K],
            [lam**k*sp.cos(k*sp.pi/4) for k in K], [lam**k*sp.sin(k*sp.pi/4) for k in K]]
    return sp.factor(sp.nsimplify(sp.simplify(sp.Matrix(rows).det())))


for K, claim in [((0, 1, 2, 3), None),
                 ((0, 1, 2, 4), lam*(lam - 1)**2*(lam + 1)*(lam**2 + 1)/sp.sqrt(2)),
                 ((0, 2, 3, 4), lam**2*(lam - 1)**2*(lam + 1)*(lam**2 + 1)/sp.sqrt(2))]:
    t = tK(K)
    msg = "" if claim is None else f"; difference from document = {sp.simplify(t - claim)}"
    print(f"t_{K}(lam) = {t}{msg}")
# number of 4-subsets of N with given index sum m versus the bound (m+1)^3
cnt = {}
for K in combinations(range(40), 4):
    cnt[sum(K)] = cnt.get(sum(K), 0) + 1
print("#{K: |K| = m} <= (m+1)^3 for m <= 39:", all(cnt.get(m, 0) <= (m + 1)**3 for m in range(40)))

# ---------------------------------------------------------------------------------------------
print("\n=== (4) interface resistance: 5-parameter moments ===")
cR = riccati_moments(NK, with_R=True)
print("c_0 (with R) =", cR[0], "; document: R + xi/b - xi b/b_s^2; difference =",
      sp.simplify(cR[0] - (R + xi/b - xi*b/bs**2)))
paramsR = [xi, b, tR, tl, R]
GR = sp.Matrix([[sp.diff(ck, q) for q in paramsR] for ck in cR])
GR01234 = sp.factor(sp.cancel(GR.extract([0, 1, 2, 3, 4], [0, 1, 2, 3, 4]).det()))
Q = 3*bs**2*(tR - tl) + 2*R*b**3*xi - R**2*b**2*bs**2
print("G^R_01234 =", GR01234)
print("  document: -2 xi^5 (tau_R - tau_l) Q/(3 b^3 b_s^4); difference =",
      sp.simplify(GR01234 + 2*xi**5*(tR - tl)*Q/(3*b**3*bs**4)))
print("  at R = 0:", sp.factor(GR01234.subs(R, 0)))
kvec = sp.Matrix([b*bs**2/(b**2 - bs**2), 0, 0, 0, 1])
res = [sp.simplify((GR.row(k)*kvec)[0].subs({tR: tl, R: 0})) for k in range(NK)]
print("grad c_k . (b b_s^2/(b^2-b_s^2), 0,0,0,1) on Delta at R = 0, k = 0..7:", res)
res2 = [sp.simplify((GR.row(k)*sp.Matrix([1, 0, 0, 0, 0]))[0].subs({tR: tl, R: 0, b: bs})) for k in range(NK)]
print("grad c_k . (1,0,0,0,0) on D2 at R = 0, k = 0..7:", res2)
res3 = [sp.simplify((GR.row(k)*sp.Matrix([0, 0, 1, 1, 0]))[0].subs({tR: tl})) for k in range(NK)]
print("grad c_k . (0,0,1,1,0) on Delta, k = 0..7:", res3)
pt = {xi: sp.Rational(3, 2), b: sp.Rational(5, 2), bs: 1, tl: sp.Rational(2, 7)}
for lab, sub in [("Delta, R = 1/3", {tR: tl, R: sp.Rational(1, 3)}), ("Delta, R = 0, b != b_s", {tR: tl, R: 0}),
                 ("Delta, R = 0, b = b_s", {tR: tl, R: 0, b: bs}), ("off Delta, R = 1/3", {tR: sp.Rational(9, 7), R: sp.Rational(1, 3)}),
                 ("off Delta, R = 0, b = b_s", {tR: sp.Rational(9, 7), R: 0, b: bs})]:
    Mn = GR.subs(sub).subs(pt)
    print(f"   rank of [grad c_0..c_7] (5 params) on {lab}: {Mn.rank()}")

# ---------------------------------------------------------------------------------------------
print("\n=== (5) Theorem 2' step (star): algebraic equivalence ===")
q, k = sp.symbols('q k', positive=True)
gq = (1 + tR*q**2)/(1 + tl*q**2)
Bq = b*q/sp.sqrt(gq)
Lq = R + 1/(bs*q)
lhs = 2*Bq/(Bq**2*Lq**2 - 1) - 2*k*q*sp.sqrt(gq)
target = b + k*gq - k*b**2*(R*q + 1/bs)**2
ratio = sp.simplify(lhs*(b**2*q**2*Lq**2 - gq)/(2*q*sp.sqrt(gq)) - target)
print("lhs * (b^2 q^2 L^2 - g)/(2 q sqrt g) - [b + k g - k b^2 (R q + 1/b_s)^2] =", ratio)
print("odd part in q of b + k g(q^2) - k b^2 (R q + 1/b_s)^2:",
      sp.factor(sp.expand((target - target.subs(q, -q))/2)))
# d r / dR identity
Bv, Lv = sp.symbols('B L')
r = (Bv*Lv - 1)/(Bv*Lv + 1)
print("dr/dR - 2B/(BL+1)^2 =", sp.simplify(sp.diff(r.subs(Lv, R + Lv), R).subs(R, 0) - 2*Bv/(Bv*Lv + 1)**2))
