"""Audit of Lemma K (proofs/B-revision.md, Sec. 2): the algebra between D_c Phi == 0 and (i), (O), (E).

Independent of kernel_system.py.  A derivation D (directional derivative in parameter space) is built by the chain
rule on base quantities:  xi -> c_xi, b -> c_b, R -> c_R, G = sqrt(g) -> dG, N = 1 + tau_n p -> dtn p,
N0 = 1 + tau_o p -> dto p;  s is not differentiated.  Checks:
 (1) (K.1):  B~^2 (1 - Q)^2 D Phi  ==  -D B~ + 2 B~ phi E + r^2 D B~ E^2,  Phi = (1/B~)(1+Q)/(1-Q), Q = r E,
     E = exp(-2A), phi = D r - 2 r D A.
 (2) Invariance of Phi under (A, B, r, E) -> (-A, -B, 1/r, 1/E) and the factor relating the two forms of (K.1).
 (3) Stack form: B~ (1-Q)^2 D W  with W = R1 + Phi equals  c_R1 B~ (1-Q)^2 - (D B~/B~)(1-Q^2) + 2 D Q.
 (4) With (i) imposed (c_b = 0, dG/G = -dto p/N0):  D B~ == 0, and
       N0 N L (D u - (u^2 - 1) D A) / u  ==  pL (N dto - N0 dtn) + c_R N N0 - (b p L^2 - g N^2/b) kappa,
     kappa = c_xi N0 - xi dto p   (so the identity of the text is exactly D u = (u^2-1) D A).
 (5) Odd/even split in s with L = R + 1/(bs s): the s-odd part is (O)/bs * s-factor and the s-even part minus
     R s^2 (O) is (E).
 (6) Special path for tau_R = 0 < tau_l: A^2 = xi^2 p/(1 + tl p) on s = i tl^(-1/2) + rho e^{i pi/4}
     behaves as c rho^-1 e^{i pi/4}, c = xi^2/(2 tl^{3/2}) (series in rho).
Output: a02_lemmaK_symbolic.out
"""
import os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.append(s)


# ------------------------------------------------------------------ (1)-(3) abstract form
Bt, r, A, dBt, dr, dA, cR1 = sp.symbols("Bt r A dBt dr dA c_R1")
E = sp.exp(-2 * A)
Q = r * E
dQ = dr * E + r * (-2 * dA) * E
Phi = (1 + Q) / (Bt * (1 - Q))
DPhi = sp.diff(Phi, Bt) * dBt + sp.diff(Phi, r) * dr + sp.diff(Phi, A) * dA
phi = dr - 2 * r * dA
K1 = -dBt + 2 * Bt * phi * E + r ** 2 * dBt * E ** 2
log("(1) B~^2(1-Q)^2 DPhi - K.1 =", sp.simplify(sp.expand(Bt ** 2 * (1 - Q) ** 2 * DPhi - K1)))
# (2) invariance
Phi_sub = Phi.subs({Bt: -Bt, r: 1 / r, A: -A}, simultaneous=True)
log("(2) Phi(-B,1/r,-A) - Phi =", sp.simplify(Phi_sub - Phi))
phi_sub = (sp.diff(1 / r, r) * dr) - 2 * (1 / r) * (-dA)
K1_sub = -(-dBt) + 2 * (-Bt) * phi_sub * sp.exp(2 * A) + r ** -2 * (-dBt) * sp.exp(4 * A)
log("(2) K.1(transformed) * (r^2 E^2) - K.1 =", sp.simplify(sp.expand(K1_sub * (r ** 2 * E ** 2) - K1)))
# (3) stack
W = cR1 * 0 + Phi
DW = cR1 + DPhi
stack = cR1 * Bt * (1 - Q) ** 2 - (dBt / Bt) * (1 - Q ** 2) + 2 * dQ
log("(3) B~(1-Q)^2 DW - stack form =", sp.simplify(sp.expand(Bt * (1 - Q) ** 2 * DW - stack)))
# coefficients in E of the stack form
st = sp.expand(stack.subs(E, sp.Symbol("EE")).subs(sp.exp(-2 * A), sp.Symbol("EE")))
EE = sp.Symbol("EE")
st = sp.expand((cR1 * Bt * (1 - r * EE) ** 2 - (dBt / Bt) * (1 - r ** 2 * EE ** 2) + 2 * EE * phi))
log("(3) stack E^0, E^1, E^2 coefficients:", [sp.factor(st.coeff(EE, k)) for k in range(3)])

# ------------------------------------------------------------------ (4)-(5) explicit translation
s, xi, b, bs, R, G, tn, to, tR, tl = sp.symbols("s xi b b_s R G tau_n tau_o tau_R tau_l", positive=True)
cxi, cb, cR, dG, dtn, dto = sp.symbols("c_xi c_b c_R dG dtn dto")
p = s ** 2
N = 1 + tn * p
N0 = 1 + to * p
L = R + 1 / (bs * s)
base = {xi: cxi, b: cb, R: cR, G: dG, tn: dtn, to: dto}


def Dc(expr):
    return sp.expand(sum(sp.diff(expr, v) * dv for v, dv in base.items()))


Aex = xi * s * G
Bex = b * s / G
Btex = Bex / N0
uex = Bex * L / N
imp_i = {cb: 0, dG: -G * dto * p / N0}
DBt_i = sp.simplify(Dc(Btex).subs(imp_i))
log("(4) D B~ with (i) imposed =", DBt_i)
lhs = sp.simplify(N0 * N * L * (Dc(uex) - (uex ** 2 - 1) * Dc(Aex)) / uex)
lhs = sp.simplify(lhs.subs(imp_i))
g = G ** 2
kappa = cxi * N0 - xi * dto * p
rhs = p * L * (N * dto - N0 * dtn) + cR * N * N0 - (b * p * L ** 2 - g * N ** 2 / b) * kappa
log("(4) N0 N L (Du - (u^2-1) DA)/u - [pL(N dto - N0 dtn) + c_R N N0 - (b p L^2 - g N^2/b) kappa] =",
    sp.simplify(sp.expand(lhs - rhs)))
# (5) odd / even split of the combined identity in s (g, N, N0, kappa are even in s)
comb = sp.expand(rhs * s ** 0)
odd = sp.expand((comb - comb.subs(s, -s)) / 2)
even = sp.expand((comb + comb.subs(s, -s)) / 2)
O_id = N * dto - N0 * dtn - 2 * b * R * kappa
E_id = cR * N * N0 - kappa * (b / bs ** 2 - b * R ** 2 * p - g * N ** 2 / b)
log("(5) odd part - (s/bs)*(O) =", sp.simplify(odd - (s / bs) * O_id),
    "  [the odd part of pL is s/bs]")
log("(5) even part - R s^2 (O) - (E) =", sp.simplify(even - R * s ** 2 * O_id - E_id))

# ------------------------------------------------------------------ (6) special path
rho = sp.Symbol("rho", positive=True)
sp_path = sp.I / sp.sqrt(tl) + rho * sp.exp(sp.I * sp.pi / 4)
A2 = xi ** 2 * sp_path ** 2 / (1 + tl * sp_path ** 2)
lead = sp.series(sp.simplify(A2 * rho), rho, 0, 2).removeO()
log("(6) rho * A^2 on the path, series to O(rho^2):", sp.simplify(lead))
log("    expected leading term xi^2/(2 tl^(3/2)) e^{i pi/4} =", sp.simplify(xi ** 2 / (2 * tl ** sp.Rational(3, 2)) * sp.exp(sp.I * sp.pi / 4)))
log("    difference of the rho^0 term:",
    sp.simplify(sp.expand_complex(lead.subs(rho, 0) - xi ** 2 / (2 * tl ** sp.Rational(3, 2)) * sp.exp(sp.I * sp.pi / 4))))

# ------------------------------------------------------------------ r == 0 is impossible on {tau_R = 0 < tau_l}
# u == 1  <=>  b (R s + 1/bs) sqrt(1 + tl s^2) == 1 + tn s^2 ; square and compare
tnn = sp.Symbol("tn_real", real=True)
expr = sp.expand(b ** 2 * (R * s + 1 / bs) ** 2 * (1 + tl * s ** 2) - (1 + tnn * s ** 2) ** 2)
log("(7) tau_R = 0: b^2 (R s + 1/bs)^2 (1 + tl s^2) - (1 + tn s^2)^2 coefficients in s:",
    sp.Poly(expr, s).all_coeffs())
log("    (s^3 and s^1 coefficients force R = 0; then s^4: -tn^2 = 0, s^2: b^2 tl/bs^2 - 2 tn = 0 forces tl = 0)")

with open(os.path.join(HERE, "a02_lemmaK_symbolic.out"), "w") as fh:
    fh.write("\n".join(out) + "\n")
