"""Exact lower bounds on the strata of Sigma from the small-p Puiseux coefficients (Proof B).

On Sigma = {tau_R = tau_l}:            3 x 3 minor d(z_0, z_1, z_2)/d(xi1, b, tau_R)
On Sigma' = Sigma and {b = bs}:        2 x 2 minor d(z_0, z_2)/d(b, tau_R)
Off Sigma (recorded for completeness): 4 x 4 minors of d(z_0..z_4)/d theta, at b != bs and b = bs.
Output: results/06b_sigma_minors.txt
"""
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
q = sp.Symbol("q")
xi1, b, bs = sp.symbols("xi1 b b_s", positive=True)
tR, tl, tau = sp.symbols("tau_R tau_l tau", nonnegative=True)
rho = sp.sqrt((1 + tR * q**2) / (1 + tl * q**2))
A, B = xi1 * q * rho, b * q / rho
w = B / (bs * q)
Z = (w + sp.tanh(A)) / (B * (1 + w * sp.tanh(A)))
ser = sp.expand(sp.series(Z * q, q, 0, 6).removeO())
z = {n - 1: sp.expand(sp.simplify(ser.coeff(q, n))) for n in range(6)}
out = []
J3 = sp.Matrix([[sp.diff(z[k], v) for v in (xi1, b, tR)] for k in (0, 1, 2)])
d3 = sp.factor(sp.simplify(J3.det().subs(tl, tR)))
out.append(f"Sigma: det d(z_0, z_1, z_2)/d(xi1, b, tau_R) at tau_l = tau_R = {d3}")
J2 = sp.Matrix([[sp.diff(z[k], v) for v in (b, tR)] for k in (0, 2)])
d2 = sp.factor(sp.simplify(J2.det().subs({tl: tR}).subs(b, bs)))
out.append(f"Sigma': det d(z_0, z_2)/d(b, tau_R) at tau_l = tau_R, b = bs = {d2}")
# upper bounds on Sigma for the coefficient map (consistency): all 4x4 minors vanish on tau_R = tau_l
J4 = sp.Matrix([[sp.diff(z[k], v) for v in (xi1, b, tR, tl)] for k in range(0, 5)])
mins = [sp.simplify(J4.extract(list(rows), [0, 1, 2, 3]).det().subs(tl, tR)) for rows in
        [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4), (0, 2, 3, 4), (1, 2, 3, 4)]]
out.append(f"all 4x4 minors of d(z_0..z_4)/d theta on tau_l = tau_R: {mins}")
J3b = J4.subs(tl, tR).subs(b, bs)
m3 = []
import itertools
for rows in itertools.combinations(range(5), 3):
    for cols in itertools.combinations(range(4), 3):
        m3.append(sp.simplify(J3b.extract(list(rows), list(cols)).det()))
out.append(f"all 3x3 minors of d(z_0..z_4)/d theta on tau_l = tau_R, b = bs vanish: {all(m == 0 for m in m3)}")
(HERE / "results" / "06b_sigma_minors.txt").write_text("\n".join(out) + "\n", encoding="utf-8")
print("\n".join(out))
