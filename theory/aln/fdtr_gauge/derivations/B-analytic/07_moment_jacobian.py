"""Workstream B: low-frequency (Laurent) coefficients and their parameter Jacobian.

F(s) = Z(p = s^2) = 1/(b_s s) + sum_{k>=0} c_k s^k near s = 0 (the expansion is analytic in theta,
including the boundary tau = 0). If the matrix [grad c_0; ...; grad c_K] has rank 4 at theta, the
function-level Jacobian has rank >= 4 there (Taylor coefficients are continuous linear functionals of
the response, determined by the data through analytic continuation). The 4x4 minors are rational
functions of theta; their common zero set in Theta is read off the minors printed here.
Run:  python 07_moment_jacobian.py > 07_moment_jacobian.out
"""
import sympy as sp

xi, b, bs, tR, tl, s = sp.symbols('xi b b_s tau_R tau_l s', positive=True)
d = sp.symbols('delta')            # delta = tau_R - tau_l
K = 7                              # number of coefficients c_0..c_{K-1}

p = s**2
g = (1 + tR*p)/(1 + tl*p)
z = xi**2*p*g
w = sp.symbols('w')
Cw = sum(w**n/sp.factorial(2*n) for n in range(K + 2))          # cosh(sqrt w), truncated
Scw = sum(w**n/sp.factorial(2*n + 1) for n in range(K + 2))     # sinh(sqrt w)/sqrt w, truncated
C = Cw.subs(w, z); Sc = Scw.subs(w, z)
X = xi/b*g*Sc
Y = xi*b*p*Sc
S = bs*s
F = (C + X*S)/(Y + C*S)
ser = sp.series(F, s, 0, K).removeO()
ser = sp.expand(ser)
coeffs = [sp.factor(sp.simplify(ser.coeff(s, k))) for k in range(K)]
print("c_-1 =", sp.simplify(ser.coeff(s, -1)))
for k, c in enumerate(coeffs):
    print(f"c_{k} =", c)

params = [xi, b, tR, tl]
Jm = sp.Matrix([[sp.diff(c, q) for q in params] for c in coeffs])

print("\n--- dependence on the times: substitute tau_R = tau_l + delta ---")
for k, c in enumerate(coeffs):
    cs = sp.expand(sp.simplify(c.subs(tR, tl + d)))
    print(f"c_{k}: depends on tau_l at fixed delta? {sp.simplify(sp.diff(cs, tl)) != 0}")

print("\n--- block structure: (xi, b) block of c_0, c_1 ---")
blk = sp.factor(sp.Matrix([[Jm[0, 0], Jm[0, 1]], [Jm[1, 0], Jm[1, 1]]]).det())
print("det =", blk)

print("\n--- 4x4 minors of the moment Jacobian (rows among c_0..c_6 with index sum <= 10) ---")
from itertools import combinations
minors = {}
for rows in combinations(range(K), 4):
    if sum(rows) > 10:          # only the low-order minors are used (Theorem N2, Proof C)
        continue
    m = sp.factor(sp.cancel(Jm.extract(list(rows), [0, 1, 2, 3]).det()))
    minors[rows] = m
    if rows in [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4), (0, 1, 2, 5), (0, 2, 3, 4), (1, 2, 3, 4)]:
        print(f"rows {rows}: {m}")

print("\n--- common zero set (read off the displayed minors; the brute-force gcd over all 35 minors is not needed) ---")
print("off the resonance: G_(0,1,2,4) != 0 iff b != b_s;  at b = b_s the minor G_(0,2,3,4) is:")
print("   G_(0,2,3,4)|_{b=b_s} =", sp.factor(sp.simplify(minors[(0, 2, 3, 4)].subs(b, bs))))
print("every computed 4x4 minor vanishes on the resonance tau_R = tau_l:",
      all(sp.cancel(m.subs(tR, tl)) == 0 for m in minors.values()))
