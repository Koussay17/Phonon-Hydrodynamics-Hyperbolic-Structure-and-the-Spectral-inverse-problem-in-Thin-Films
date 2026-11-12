"""Workstream B, formal check: provenance of the Laurent coefficients used in the Lean files.

(1) Recompute c_{-1}..c_7 with derivations/B-analytic/laurent.py (power-series division), R = 0.
(2) Compare with the closed forms printed in derivations/B-analytic/07_moment_jacobian.out (c_0..c_6)
    and with the forms displayed in proofs/B-analytic.md, section 5, Proof C (c_0..c_3).
(3) Independent route: the front-face impedance obeys the Riccati equation in the reduced thickness,
        d/dxi F(s) = g(s^2)/b - b s^2 F(s)^2,      F|_{xi=0} = R + 1/(b_s s),
    (d/dA of the Moebius form of Lemma 1 equals (1 - B^2 Z^2)/B; A/B = (xi/b) g, A B = xi b p).
    With s F = sum_k q_k s^k (q_0 = 1/b_s, q_{k+1} = c_k) this is the triangular hierarchy
        d q_{n+1}/dxi = [n even] gamma_{n/2}/b - b sum_{j=0}^{n} q_j q_{n-j},
        q_0 = 1/b_s, q_1(0) = R, q_k(0) = 0 (k >= 2),
    gamma_0 = 1, gamma_m = (tau_R - tau_l)(-tau_l)^(m-1) (m >= 1), g(s^2) = sum_m gamma_m s^(2m).
    The hierarchy is integrated here directly (no series division) and compared with (1).
(4) Recompute the minors G_0123, G_0124, G_0234 and G_0234|_{b=b_s} from the coefficients of (1).
(5) Print the coefficients and their partial derivatives in the variable names of the Lean files.
Run:  python 01_coefficients.py > 01_coefficients.out
"""
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'B-analytic'))
from laurent import laurent_coefficients  # noqa: E402

xi, b, bs, tR, tl, s, R = sp.symbols('xi b b_s tau_R tau_l s R', positive=True)
NK = 9                                    # q_0..q_8  ->  c_{-1}..c_7

# ---------------------------------------------------------------- (1) series division, R = 0
q = laurent_coefficients(xi, b, bs, tR, tl, 0, s, NK)
c_div = q[1:]                             # c_0..c_7
print("=== (1) laurent.py, R = 0 ===")
print("c_-1 =", sp.factor(q[0]))
for k, ck in enumerate(c_div):
    print(f"c_{k} =", sp.factor(ck))

# ---------------------------------------------------------------- (2) comparison with 07 and the proof
loc = {'xi': xi, 'b': b, 'b_s': bs, 'tau_R': tR, 'tau_l': tl}
c07_txt = [
    "-xi*(b - b_s)*(b + b_s)/(b*b_s**2)",
    "xi**2*(b - b_s)*(b + b_s)/b_s**3",
    "-xi*(3*b**4*xi**2 - 4*b**2*b_s**2*xi**2 - 3*b_s**4*tau_R + 3*b_s**4*tau_l + b_s**4*xi**2)/(3*b*b_s**4)",
    "xi**2*(3*b**4*xi**2 - 5*b**2*b_s**2*xi**2 - 3*b_s**4*tau_R + 3*b_s**4*tau_l + 2*b_s**4*xi**2)/(3*b_s**5)",
    "-xi*(15*b**6*xi**4 - 30*b**4*b_s**2*xi**4 - 20*b**2*b_s**4*tau_R*xi**2 + 20*b**2*b_s**4*tau_l*xi**2"
    " + 17*b**2*b_s**4*xi**4 + 15*b_s**6*tau_R*tau_l + 10*b_s**6*tau_R*xi**2 - 15*b_s**6*tau_l**2"
    " - 10*b_s**6*tau_l*xi**2 - 2*b_s**6*xi**4)/(15*b*b_s**6)",
    "xi**2*(45*b**6*xi**4 - 105*b**4*b_s**2*xi**4 - 75*b**2*b_s**4*tau_R*xi**2 + 75*b**2*b_s**4*tau_l*xi**2"
    " + 77*b**2*b_s**4*xi**4 + 45*b_s**6*tau_R*tau_l + 60*b_s**6*tau_R*xi**2 - 45*b_s**6*tau_l**2"
    " - 60*b_s**6*tau_l*xi**2 - 17*b_s**6*xi**4)/(45*b_s**7)",
    "-xi*(315*b**8*xi**6 - 840*b**6*b_s**2*xi**6 - 630*b**4*b_s**4*tau_R*xi**4 + 630*b**4*b_s**4*tau_l*xi**4"
    " + 756*b**4*b_s**4*xi**6 + 420*b**2*b_s**6*tau_R*tau_l*xi**2 + 714*b**2*b_s**6*tau_R*xi**4"
    " - 420*b**2*b_s**6*tau_l**2*xi**2 - 714*b**2*b_s**6*tau_l*xi**4 - 248*b**2*b_s**6*xi**6"
    " + 105*b_s**8*tau_R**2*xi**2 - 315*b_s**8*tau_R*tau_l**2 - 420*b_s**8*tau_R*tau_l*xi**2"
    " - 126*b_s**8*tau_R*xi**4 + 315*b_s**8*tau_l**3 + 315*b_s**8*tau_l**2*xi**2 + 126*b_s**8*tau_l*xi**4"
    " + 17*b_s**8*xi**6)/(315*b*b_s**8)",
]
c07 = [sp.sympify(t, locals=loc) for t in c07_txt]
proof_txt = [
    "xi*(b_s**2 - b**2)/(b*b_s**2)",
    "xi**2*(b**2 - b_s**2)/b_s**3",
    "xi*(tau_R - tau_l)/b - xi**3*(3*b**4 - 4*b**2*b_s**2 + b_s**4)/(3*b*b_s**4)",
    "xi**2*(3*b_s**4*(tau_l - tau_R) + (3*b**4 - 5*b**2*b_s**2 + 2*b_s**4)*xi**2)/(3*b_s**5)",
]
cproof = [sp.sympify(t, locals=loc) for t in proof_txt]
print("\n=== (2) comparison ===")
for k in range(7):
    print(f"c_{k}: laurent.py - 07 output == 0 :", sp.simplify(c_div[k] - c07[k]) == 0)
for k in range(4):
    print(f"c_{k}: laurent.py - proof text == 0 :", sp.simplify(c_div[k] - cproof[k]) == 0)

# ---------------------------------------------------------------- (3) Riccati hierarchy, general R
print("\n=== (3) Riccati hierarchy (independent of the series division) ===")
x = sp.symbols('x')


def gamma(m):
    return sp.Integer(1) if m == 0 else (tR - tl)*(-tl)**(m - 1)


def riccati_coefficients(Rval, N):
    """q_0..q_{N-1} by integrating the hierarchy in xi from the initial data at xi = 0."""
    qs = [1/bs]
    for n in range(N - 1):
        rhs = (gamma(n//2)/b if n % 2 == 0 else 0) - b*sum(qs[j]*qs[n - j] for j in range(n + 1))
        init = Rval if n == 0 else 0
        qn1 = init + sp.integrate(sp.expand(rhs).subs(xi, x), (x, 0, xi))
        qs.append(sp.expand(qn1))
    return qs


q_ric = riccati_coefficients(0, NK)
for k in range(NK):
    print(f"q_{k}: Riccati - laurent.py == 0 :", sp.simplify(q_ric[k] - q[k]) == 0)
qR = laurent_coefficients(xi, b, bs, tR, tl, R, s, NK)
q_ricR = riccati_coefficients(R, NK)
for k in range(NK):
    print(f"q_{k} (with R): Riccati - laurent.py == 0 :", sp.simplify(q_ricR[k] - qR[k]) == 0)
# the hierarchy residual of the laurent.py coefficients themselves
for n in range(NK - 1):
    rhs = (gamma(n//2)/b if n % 2 == 0 else 0) - b*sum(q[j]*q[n - j] for j in range(n + 1))
    res = sp.simplify(sp.diff(q[n + 1], xi) - rhs)
    print(f"hierarchy residual n={n} (R = 0):", res, "; initial value q_{n+1}(0) =", sp.simplify(q[n + 1].subs(xi, 0)))

# ---------------------------------------------------------------- (4) minors
print("\n=== (4) minors from the laurent.py coefficients ===")
params = [xi, b, tR, tl]
Jm = sp.Matrix([[sp.diff(ck, v) for v in params] for ck in c_div])


def G(rows):
    return sp.factor(sp.cancel(Jm.extract(list(rows), [0, 1, 2, 3]).det()))


print("G_0123 =", G((0, 1, 2, 3)))
print("G_0124 =", G((0, 1, 2, 4)))
print("G_0234 =", G((0, 2, 3, 4)))
print("G_0234|_{b=b_s} =", sp.factor(sp.cancel(G((0, 2, 3, 4)).subs(b, bs))))
print("G_0124 - 2 xi^4 (b^2-b_s^2)(tau_l-tau_R)/(b^4 b_s^3) == 0 :",
      sp.simplify(G((0, 1, 2, 4)) - 2*xi**4*(b**2 - bs**2)*(tl - tR)/(b**4*bs**3)) == 0)
print("G_0234|_{b=b_s} + 2 xi^4 (tau_R-tau_l)^2/b_s^5 == 0 :",
      sp.simplify(G((0, 2, 3, 4)).subs(b, bs) + 2*xi**4*(tR - tl)**2/bs**5) == 0)

# ---------------------------------------------------------------- (5) Lean-syntax output
print("\n=== (5) Lean syntax (xi b bs tR tl) ===")
ren = {xi: sp.Symbol('xi'), b: sp.Symbol('b'), bs: sp.Symbol('bs'), tR: sp.Symbol('tR'), tl: sp.Symbol('tl')}


def lean(e):
    return sp.sstr(e.subs(ren)).replace('**', '^')


for k in range(8):
    print(f"c{k} := {lean(sp.factor(c_div[k]))}")
for k in range(5):
    for v, nm in zip(params, ['xi', 'b', 'tR', 'tl']):
        print(f"d c{k} / d {nm} = {lean(sp.factor(sp.diff(c_div[k], v)))}")
