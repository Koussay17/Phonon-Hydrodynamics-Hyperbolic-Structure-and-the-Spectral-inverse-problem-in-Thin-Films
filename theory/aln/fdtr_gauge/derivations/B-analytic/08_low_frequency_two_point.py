"""Workstream B: two low frequencies suffice at every theta off the resonance (Theorem N2).

With s = e^{i pi/4} sigma (p = i sigma^2) and J = dZ/dtheta = sum_k grad c_k s^k (the principal part
1/(b_s s) does not depend on theta), the real 4x4 matrix at sigma_1 = sigma, sigma_2 = lam*sigma factorises
by Cauchy-Binet: D = sum_K t_K(lam) sigma^{|K|} G_K, with t_K the 4x4 determinant of rows
(cos k pi/4), (sin k pi/4), (lam^k cos k pi/4), (lam^k sin k pi/4), k in K, and G_K the minor of [grad c_k].
Leading terms: |K| = 7, K = (0,1,2,4) if b != b_s;  |K| = 9, K = (0,2,3,4) if b = b_s.
(a) SymPy: t_K(lam) and the G-minors with |K| <= 9 (generic b and b = b_s).
(b) mpmath (80 digits): D(sigma)/(sigma^{|K|} t_K G_K) -> 1 on the interior, b = b_s, and both boundary strata.
(c) Interface resistance: 5x5 minors of the moment Jacobian in (xi, b, tau_R, tau_l, R).
Laurent coefficients from laurent.py (truncated power-series division; cross-checked against 07).
Run:  python 08_low_frequency_two_point.py > 08_low_frequency_two_point.out
"""
from itertools import combinations
import sympy as sp
import mpmath as mp
from model import jac_analytic
from laurent import laurent_coefficients

xi, b, bs, tR, tl, s, R, lam = sp.symbols('xi b b_s tau_R tau_l s R lam', positive=True)
NK = 9                                                   # q_0..q_8  ->  c_{-1}..c_7

print("=== (a) T-determinants ===")


def Tdet(K):
    rows = [[sp.cos(k*sp.pi/4) for k in K], [sp.sin(k*sp.pi/4) for k in K],
            [lam**k*sp.cos(k*sp.pi/4) for k in K], [lam**k*sp.sin(k*sp.pi/4) for k in K]]
    return sp.factor(sp.simplify(sp.Matrix(rows).det()))


for K in [(0, 1, 2, 3), (0, 1, 2, 4), (0, 2, 3, 4), (0, 1, 2, 5), (0, 1, 3, 4)]:
    print(f"  t_{K}(lam) =", Tdet(K))

q = laurent_coefficients(xi, b, bs, tR, tl, 0, s, NK)
cs = q[1:]                                              # c_0 .. c_7
params = [xi, b, tR, tl]
G = sp.Matrix([[sp.diff(c, v) for v in params] for c in cs])
print("\n  G-minors with |K| <= 9 (generic b; and at b = b_s):")
for K in combinations(range(len(cs)), 4):
    if sum(K) > 9:
        continue
    m = sp.factor(sp.cancel(G.extract(list(K), [0, 1, 2, 3]).det()))
    mb = sp.factor(sp.cancel(m.subs(b, bs)))
    print(f"   K={K} |K|={sum(K)}:  G_K = {m}\n        at b = b_s: {mb}")

print("\n=== (b) high-precision check of the leading term ===")
mp.mp.dps = 80
bsn = mp.mpf(1)
Gf = sp.lambdify((xi, b, bs, tR, tl), G, modules='mpmath')


def Dsig(th, sig, lamv):
    rows = []
    for sg in (sig, lamv*sig):
        J = jac_analytic(mp.mpc(0, sg**2), *th, bsn)
        rows.append([mp.re(x) for x in J])
        rows.append([mp.im(x) for x in J])
    return mp.det(mp.matrix(rows))


for label, th in [("interior", (1, 2, '0.7', '0.3')), ("b = b_s", (1, 1, '0.7', '0.3')),
                  ("tau_l = 0", (1, 2, '0.7', 0)), ("tau_R = 0", (1, 2, 0, '0.5')),
                  ("b = b_s, tau_l = 0", (1, 1, '0.7', 0)), ("b = b_s, tau_R = 0", (1, 1, 0, '0.5')),
                  ("AlN-like", (1, '2.70', '0.532', '0.266'))]:
    th = [mp.mpf(x) for x in th]
    Gn = Gf(th[0], th[1], bsn, th[2], th[3])
    K = (0, 1, 2, 4) if th[1] != bsn else (0, 2, 3, 4)
    Gl = Gn.tolist()
    GK = mp.det(mp.matrix([[mp.mpf(Gl[k][j]) for j in range(4)] for k in K]))
    lamv = mp.mpf(2)
    TK = mp.mpf(str(sp.N(Tdet(K).subs(lam, 2), 90)))
    out = []
    for sig in [mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4')]:
        out.append(Dsig(th, sig, lamv)/(sig**sum(K)*TK*GK))
    print(f"  {label}: K={K}, G_K={mp.nstr(GK, 6)};  D/(sigma^|K| t_K G_K) at sigma = 1e-2, 1e-3, 1e-4:",
          [mp.nstr(x, 10) for x in out])

print("\n=== (c) interface resistance: 5x5 minors of the moment Jacobian ===")
qR = laurent_coefficients(xi, b, bs, tR, tl, R, s, NK)
csR = qR[1:]
paramsR = [xi, b, tR, tl, R]
GR = sp.Matrix([[sp.diff(c, v) for v in paramsR] for c in csR])
print("  c_0 =", sp.factor(csR[0]))
shown = 0
for K in combinations(range(len(csR)), 5):
    if sum(K) > 13:
        continue
    m = sp.factor(sp.cancel(GR.extract(list(K), [0, 1, 2, 3, 4]).det()))
    if m == 0:
        print(f"   K={K} |K|={sum(K)}: 0")
        continue
    mb = sp.factor(sp.cancel(m.subs(b, bs)))
    m0 = sp.factor(sp.cancel(m.subs(R, 0)))
    print(f"   K={K} |K|={sum(K)}: {m}\n        at b=b_s: {mb}\n        at R=0: {m0}")
    shown += 1
    if shown >= 8:
        break
# rank on the resonance strata at a sample point
for label, sub in [("resonance, R > 0", {tR: tl}), ("resonance, R = 0, b != b_s", {tR: tl, R: 0}),
                   ("resonance, R = 0, b = b_s", {tR: tl, R: 0, b: bs}), ("off resonance, R > 0", {}),
                   ("off resonance, R = 0, b = b_s", {R: 0, b: bs})]:
    vals = {xi: sp.Rational(1), b: sp.Rational(2), bs: 1, tR: sp.Rational(7, 10), tl: sp.Rational(1, 2), R: sp.Rational(3, 10)}
    Mr = GR.subs(sub).subs(vals)
    print(f"   {label}: rank of [grad c_0..grad c_7] at a sample point = {Mr.rank()}")
