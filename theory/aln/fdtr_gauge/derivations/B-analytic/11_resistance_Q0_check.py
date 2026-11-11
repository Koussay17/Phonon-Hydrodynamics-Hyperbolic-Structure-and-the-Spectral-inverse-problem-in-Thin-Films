"""Workstream B: moment-Jacobian rank of the resistance model on the hypersurface Q = 0.

G^R_(0,1,2,3,4) = -2 xi^5 (tau_R - tau_l) Q / (3 b^3 b_s^4),  Q = 3 b_s^2 (tau_R - tau_l) + 2 R b^3 xi - R^2 b^2 b_s^2,
and every minor containing rows (0,1,2,3) carries the factor Q.
(1) Exact rank of [grad c_0; ...; grad c_7] in (xi, b, tau_R, tau_l, R) at rational points with Q = 0, tau_R != tau_l.
(2) Elimination: on Q = 0 put tau_R = tau_l + d(R), d = (R^2 b^2 b_s^2 - 2 R b^3 xi)/(3 b_s^2); factor the minors
    (0,1,2,4,5), (0,1,3,4,5), (0,1,2,4,6) there and report their common factors (gcd).
Run:  python 11_resistance_Q0_check.py > 11_resistance_Q0_check.out
"""
import sympy as sp
from laurent import laurent_coefficients

xi, b, bs, tR, tl, s, R = sp.symbols('xi b b_s tau_R tau_l s R', positive=True)
qR = laurent_coefficients(xi, b, bs, tR, tl, R, s, 9)
csR = qR[1:]
GR = sp.Matrix([[sp.diff(c, v) for v in (xi, b, tR, tl, R)] for c in csR])
Q = 3*bs**2*(tR - tl) + 2*R*b**3*xi - R**2*b**2*bs**2

print("=== (1) exact ranks at rational points with Q = 0 ===")
for pt in [(1, 2, 1, 5, sp.Rational(1, 2)), (1, 1, 1, 3, sp.Rational(1, 3)),
           (sp.Rational(3, 2), sp.Rational(1, 2), 1, 7, 1), (1, 2, 1, sp.Rational(9, 2), 0)]:
    xv, bv, bsv, Rv, tlv = [sp.nsimplify(v) for v in pt]          # exact rationals throughout
    tRv = tlv + (Rv**2*bv**2*bsv**2 - 2*Rv*bv**3*xv)/(3*bsv**2)
    if tRv < 0:
        print("  skip (tau_R < 0)", (xv, bv, bsv, Rv, tlv)); continue
    vals = {xi: xv, b: bv, bs: bsv, R: Rv, tl: tlv, tR: tRv}
    assert sp.simplify(Q.subs(vals)) == 0
    M = GR.subs(vals)
    print(f"  xi={xv}, b={bv}, b_s={bsv}, R={Rv}, tau_l={tlv}, tau_R={tRv}: Q = 0, "
          f"rank[grad c_0..c_7] = {M.rank()}, rank[grad c_0..c_4] = {M[:5, :].rank()}")

print("\n=== (2) minors restricted to Q = 0 ===")
dR = (R**2*b**2*bs**2 - 2*R*b**3*xi)/(3*bs**2)
restricted = []
for K in [(0, 1, 2, 4, 5), (0, 1, 3, 4, 5), (0, 1, 2, 4, 6)]:
    m = sp.cancel(GR.extract(list(K), [0, 1, 2, 3, 4]).det())
    mq = sp.factor(sp.cancel(m.subs(tR, tl + dR)))
    restricted.append(sp.numer(sp.together(mq)))
    print(f"  K={K} on Q=0: {mq}")
g = sp.gcd_list([sp.expand(x) for x in restricted])
print("  gcd of the numerators on Q = 0:", sp.factor(g))
