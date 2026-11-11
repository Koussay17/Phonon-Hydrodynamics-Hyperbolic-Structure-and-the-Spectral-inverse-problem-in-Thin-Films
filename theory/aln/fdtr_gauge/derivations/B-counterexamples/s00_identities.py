"""Symbolic checks (SymPy) of the identities used by the counterexample branch.

(1) transmission-line form = (1/B)(1 + rho E)/(1 - rho E), E = exp(-2A), rho = (y-1)/(y+1), y = B Zs;
    the formula as typed in the brief equals the same expression with y -> 1/y.
(2) u = A dZ/dA and v = B dZ/dB (A fixed, Zs fixed): closed forms used in zmodel.py.
(3) Q := -v/u = sinh A [cosh A (1 + y^2) + 2 y sinh A] / (A (1 - y^2)).
(4) tau-columns: dZ/dtR = gR (u - v), dZ/dtl = -gl (u - v), g = p/(2(1 + t p)).
(5) two-frequency (case 2) algebra: with zeta = 2c(1 + tR p)(1 + tl p)/p on p = i w:
    Re zeta = 2c(tR + tl), Re(zeta gR) = Re(zeta gl) = c.
"""
import sympy as sp

A, B, y, E, Zs, s, bs = sp.symbols("A B y E Z_s s b_s", positive=True)
t = sp.tanh(A)
# (1)
Z_tl = (Zs + t / B) / (1 + B * Zs * t)
rho = (B * Zs - 1) / (B * Zs + 1)
Ex = sp.exp(-2 * A)
Z_rho = (1 + rho * Ex) / (B * (1 - rho * Ex))
d1 = sp.simplify((Z_tl - Z_rho).rewrite(sp.exp))
print("(1) transmission-line minus rho-form:", d1)
Z_brief = (1 / B) * (B * t + bs * s) / (B + bs * s * t)
Zs0 = 1 / (bs * s)
rho_b = (1 / (B * Zs0) - 1) / (1 / (B * Zs0) + 1)
Z_brief_rho = (1 + rho_b * Ex) / (B * (1 - rho_b * Ex))
print("(1) brief formula minus (y -> 1/y) rho-form:", sp.simplify((Z_brief - Z_brief_rho).rewrite(sp.exp)))
print("(1) brief formula at A -> 0:", sp.simplify(Z_brief.subs(A, 0)), " (physical value 1/(bs s))")
# (2)
Zf = lambda AA, BB: (Zs + sp.tanh(AA) / BB) / (1 + BB * Zs * sp.tanh(AA))
u = sp.simplify(A * sp.diff(Zf(A, B), A))
v = sp.simplify(B * sp.diff(Zf(A, B), B))
yy = B * Zs
rr = (yy - 1) / (yy + 1)
u_cf = -4 * A * rr * Ex / (B * (1 - rr * Ex) ** 2)
v_cf = -(1 - Ex) * ((1 + Ex) * (1 + yy ** 2) + 2 * yy * (1 - Ex)) / (B * (1 - rr * Ex) ** 2 * (yy + 1) ** 2)
print("(2) u - closed form:", sp.simplify((u - u_cf).rewrite(sp.exp)))
print("(2) v - closed form:", sp.simplify((v - v_cf).rewrite(sp.exp)))
# (3)
Q = -v / u
Q_cf = sp.sinh(A) * (sp.cosh(A) * (1 + yy ** 2) + 2 * yy * sp.sinh(A)) / (A * (1 - yy ** 2))
print("(3) -v/u - closed form:", sp.simplify((Q - Q_cf).rewrite(sp.exp)))
# (4)
p, xi, b, tR, tl, w = sp.symbols("p xi b tau_R tau_l omega", positive=True)
Ap = xi * sp.sqrt(p * (1 + tR * p) / (1 + tl * p))
Bp = b * sp.sqrt(p * (1 + tl * p) / (1 + tR * p))
dA_R = sp.simplify(sp.diff(Ap, tR) / Ap); dB_R = sp.simplify(sp.diff(Bp, tR) / Bp)
dA_l = sp.simplify(sp.diff(Ap, tl) / Ap); dB_l = sp.simplify(sp.diff(Bp, tl) / Bp)
print("(4) dlnA/dtR =", dA_R, "; dlnB/dtR =", dB_R, "; dlnA/dtl =", dA_l, "; dlnB/dtl =", dB_l)
print("(4) A B =", sp.simplify(Ap * Bp), "; A/B =", sp.simplify(Ap / Bp))
# (5)
c = sp.symbols("c", real=True)
pp = sp.I * w
zeta = 2 * c * (1 + tR * pp) * (1 + tl * pp) / pp
gR = pp / (2 * (1 + tR * pp)); gl = pp / (2 * (1 + tl * pp))
print("(5) Re zeta =", sp.simplify(sp.re(sp.expand(zeta))), "; Re(zeta gR) =", sp.simplify(sp.re(sp.expand(zeta * gR))),
      "; Re(zeta gl) =", sp.simplify(sp.re(sp.expand(zeta * gl))))
print("(5) Re[(1+tR p)(1+tl p)/p * (-1/2)] =", sp.simplify(sp.re(sp.expand(-(1 + tR * pp) * (1 + tl * pp) / (2 * pp)))))
