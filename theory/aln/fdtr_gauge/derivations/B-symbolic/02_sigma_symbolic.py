"""Task 2, symbolic part.  Structure of the Jacobian columns and the degenerate locus.

Z = Phi(A, B, Zs) = (1/B)(w + t)/(1 + w t),  w = B Zs,  t = tanh A,  Zs = R + 1/S.

(A1) Phi_A, Phi_B, Phi_Zs in closed form:
       U := A Phi_A = A (1-t^2)(1-w^2) / (B (1+wt)^2)
       V := B Phi_B = -t (1 + 2wt + w^2) / (B (1+wt)^2)
       Y := Phi_Zs  = (1-t^2)/(1+wt)^2              (= dZ/dR)
(A2) log-derivatives of A and B with respect to the four parameters (full expressions):
       d ln A/d xi1 = 1/xi1, d ln B/d b = 1/b,
       d ln A/d tau_R = -d ln B/d tau_R = p/(2(1+tau_R p)) =: g_R,
       d ln A/d tau_l = -d ln B/d tau_l = -p/(2(1+tau_l p)) =: -g_l,
     hence the columns  [U/xi1, V/b, g_R (U-V), -g_l (U-V), Y].
(A3) the same columns obtained by differentiating the full expression of Z directly (SymPy diff),
     compared with (A1)-(A2) at random points (50 digits).
(A4) B (1+wt)^2 * sum_k c_k col_k is a quadratic polynomial in t; its coefficients.
(A5) strata identities: w^2 = 1 on {b = bs, tau_R = tau_l, R = 0}; U/xi1 = ((1-w^2)/b) Y on
     {tau_R = tau_l, R = 0}; d_tauR Z + d_taul Z = 0 on {tau_R = tau_l} (direct, full expression).
(B)  factorisation of the 4 x 4 determinant det[col_k(p_m)]_{m,k} (four arbitrary points):
       det = -(tau_R - tau_l)/(xi1 b) * prod_m W(p_m) * det[1, V/W, g_R, kappa](p_m),
       W = U - V,  kappa(p) = p^2 / (2 (1+tau_R p)(1+tau_l p)),
     and det[Re r1; Im r1; Re r2; Im r2] = -(1/4) det[r1; conj r1; r2; conj r2].
(C)  low-order Taylor test of the transcendence step: the t-coefficients are not identically zero
     unless w^2 = 1 (checked by series in p).
Output: results/02_sigma_symbolic.txt
"""
import random
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

lines = []
import time as _time
_T0 = _time.time()


def log(*a):
    s = " ".join(str(x) for x in a) + f"   [{_time.time() - _T0:.1f} s]"
    print(s, flush=True)
    lines.append(s)


# ---------------------------------------------------------------- (A1)
Asym, Bsym, Zs, T = sp.symbols("A B Z_s t")
Phi = (Bsym * Zs + sp.tanh(Asym)) / (Bsym * (1 + Bsym * Zs * sp.tanh(Asym)))
w = Bsym * Zs
U_f = Asym * (1 - T**2) * (1 - w**2) / (Bsym * (1 + w * T) ** 2)
V_f = -T * (1 + 2 * w * T + w**2) / (Bsym * (1 + w * T) ** 2)
Y_f = (1 - T**2) / (1 + w * T) ** 2
sub = {sp.tanh(Asym): T}
r1 = sp.cancel((Asym * sp.diff(Phi, Asym)).subs(sub) - U_f)
r2 = sp.cancel((Bsym * sp.diff(Phi, Bsym)).subs(sub) - V_f)
r3 = sp.cancel(sp.diff(Phi, Zs).subs(sub) - Y_f)
log("(A1) A Phi_A - U =", r1, ";  B Phi_B - V =", r2, ";  Phi_Zs - Y =", r3)

# ---------------------------------------------------------------- (A2)
Z, s, parts = M.sympy_model(with_R=True)
p, xi1, b, bs, tR, tl, R = (s[k] for k in ("p", "xi1", "b", "bs", "tau_R", "tau_l", "R"))
A_e, B_e = parts["A"], parts["B"]
gR = p / (2 * (1 + tR * p))
gl = p / (2 * (1 + tl * p))
X = (1 + tR * p) / (1 + tl * p)
lnA = sp.log(xi1) + sp.log(p) / 2 + sp.log(X) / 2      # ln A on p = i omega (principal branches add)
lnB = sp.log(b) + sp.log(p) / 2 - sp.log(X) / 2
log("(A2) A/(xi1 sqrt(p) X^(1/2)) - 1 =", sp.cancel(A_e / (xi1 * sp.sqrt(p) * sp.sqrt(X)) - 1),
    "; B/(b sqrt(p) X^(-1/2)) - 1 =", sp.cancel(B_e * sp.sqrt(X) / (b * sp.sqrt(p)) - 1))
checks = {
    "dlnA/dxi1 - 1/xi1": sp.diff(lnA, xi1) - 1 / xi1,
    "dlnB/db - 1/b": sp.diff(lnB, b) - 1 / b,
    "dlnA/dtauR - gR": sp.diff(lnA, tR) - gR,
    "dlnB/dtauR + gR": sp.diff(lnB, tR) + gR,
    "dlnA/dtaul + gl": sp.diff(lnA, tl) + gl,
    "dlnB/dtaul - gl": sp.diff(lnB, tl) - gl,
}
for k, v in checks.items():
    log(f"(A2) {k} =", sp.cancel(v))
log("     (A does not depend on b, B does not depend on xi1.)")

# ---------------------------------------------------------------- (A3)
mp.mp.dps = 50
params = [xi1, b, tR, tl, R]
dZ = [sp.diff(Z, q) for q in params]
fZ = sp.lambdify((p, xi1, b, tR, tl, bs, R), dZ, modules="mpmath")
rng = random.Random(7)
worst = mp.mpf(0)
for _ in range(30):
    vals = [mp.mpf(rng.uniform(0.1, 3)) for _ in range(6)]
    vx1, vb, vtR, vtl, vbs, vR = vals
    if _ % 3 == 0:
        vtl = mp.mpf(0)
    if _ % 4 == 0:
        vR = mp.mpf(0)
    pv = mp.mpc(0, rng.uniform(0.01, 30))
    Av, Bv, sq = M.parts_mp(pv, vx1, vb, vtR, vtl)
    Zsv = vR + 1 / (vbs * sq)
    wv = Bv * Zsv
    tv = mp.tanh(Av)
    Uv = Av * (1 - tv**2) * (1 - wv**2) / (Bv * (1 + wv * tv) ** 2)
    Vv = -tv * (1 + 2 * wv * tv + wv**2) / (Bv * (1 + wv * tv) ** 2)
    Yv = (1 - tv**2) / (1 + wv * tv) ** 2
    gRv, glv = pv / (2 * (1 + vtR * pv)), pv / (2 * (1 + vtl * pv))
    cols = [Uv / vx1, Vv / vb, gRv * (Uv - Vv), -glv * (Uv - Vv), Yv]
    ref = fZ(pv, vx1, vb, vtR, vtl, vbs, vR)
    sc = max(abs(x) for x in ref)
    worst = max(worst, max(abs(x - y) for x, y in zip(cols, ref)) / sc)
log(f"(A3) columns [U/xi1, V/b, gR(U-V), -gl(U-V), Y] vs direct SymPy derivatives, 30 random points: "
    f"max rel. difference {mp.nstr(worst, 3)}")

# ---------------------------------------------------------------- (A4)
c1, c2, c3, c4, c5, h, Aa, Bb, ww, tt, xx, bb = sp.symbols("c1 c2 c3 c4 c5 h A B w t xi1 b")
U_g = Aa * (1 - tt**2) * (1 - ww**2) / (Bb * (1 + ww * tt) ** 2)
V_g = -tt * (1 + 2 * ww * tt + ww**2) / (Bb * (1 + ww * tt) ** 2)
Y_g = (1 - tt**2) / (1 + ww * tt) ** 2
# h stands for c3 gR - c4 gl
Lc = c1 * U_g / xx + c2 * V_g / bb + h * (U_g - V_g) + c5 * Y_g
poly = sp.expand(sp.cancel(Lc * Bb * (1 + ww * tt) ** 2))
coeffs = sp.Poly(poly, tt).all_coeffs()[::-1]
log("(A4) B(1+wt)^2 * sum c_k col_k = sum_j e_j t^j with")
for j, cj in enumerate(coeffs):
    log(f"      e_{j} =", sp.factor(cj))
claim = {0: (c1 / xx + h) * Aa * (1 - ww**2) + c5 * Bb,
         1: -(c2 / bb - h) * (1 + ww**2),
         2: -(c1 / xx + h) * Aa * (1 - ww**2) - 2 * (c2 / bb - h) * ww - c5 * Bb}
log("      check against the stated forms:", [sp.cancel(sp.expand(coeffs[j] - claim[j])) for j in range(3)])

# ---------------------------------------------------------------- (A5)
w_diag = (B_e * (1 / (bs * sp.sqrt(p)))).subs({tl: tR, b: bs})
log("(A5a) w on {b = bs, tau_R = tau_l, R = 0}:", sp.simplify(w_diag))
# U/xi1 = ((1 - w^2)/b) Y on {tau_R = tau_l, R = 0}: here w = b/bs, A = xi1 sqrt p, B = b sqrt p
qq = sp.Symbol("q", positive=True)
wF = b / bs
U_F = (xi1 * qq) * (1 - T**2) * (1 - wF**2) / (b * qq * (1 + wF * T) ** 2)
Y_F = (1 - T**2) / (1 + wF * T) ** 2
log("(A5b) U/xi1 - ((1-w^2)/b) Y on {tau_R = tau_l, R = 0}:", sp.cancel(U_F / xi1 - (1 - wF**2) / b * Y_F))
log("      so dZ/dxi1 + c5 dZ/dR = 0 with c5 = -(1 - b^2/bs^2)/b there.")
# direct check of d_tauR Z + d_taul Z = 0 on the diagonal with the full expression (numerically, 50 digits)
fsum = sp.lambdify((p, xi1, b, tR, bs, R), (dZ[2] + dZ[3]).subs(tl, tR), modules="mpmath")
wd = mp.mpf(0)
for _ in range(20):
    v = [mp.mpf(rng.uniform(0.1, 3)) for _ in range(5)]
    pv = mp.mpc(0, rng.uniform(0.01, 30))
    ref = abs(fZ(pv, v[0], v[1], v[2], v[2], v[3], v[4])[2])
    wd = max(wd, abs(fsum(pv, *v)) / ref)
log(f"(A5c) |d_tauR Z + d_taul Z| / |d_tauR Z| on tau_R = tau_l (20 random points, 50 digits): {mp.nstr(wd, 3)}")
# symbolic: on the diagonal rho = 1 identically, so A, B do not depend on the common time
tau = sp.Symbol("tau", nonnegative=True)
log("(A5d) rho on the diagonal:", sp.simplify(parts["rho"].subs({tR: tau, tl: tau})),
    "; A, B there:", sp.simplify(A_e.subs({tR: tau, tl: tau})), ",", sp.simplify(B_e.subs({tR: tau, tl: tau})))

# ---------------------------------------------------------------- (B)
from fractions import Fraction as Fr
import random as _r


def det_exact(Mx):
    """Exact determinant of a square list-of-lists of Fractions (Gaussian elimination)."""
    Mx = [row[:] for row in Mx]
    n, d = len(Mx), Fr(1)
    for k in range(n):
        piv = next((i for i in range(k, n) if Mx[i][k] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != k:
            Mx[k], Mx[piv] = Mx[piv], Mx[k]
            d = -d
        d *= Mx[k][k]
        for i in range(k + 1, n):
            f = Mx[i][k] / Mx[k][k]
            for j in range(k, n):
                Mx[i][j] -= f * Mx[k][j]
    return d


_rr = _r.Random(99)
nz, ntest = 0, 200
for _ in range(ntest):
    rq = lambda: Fr(_rr.randint(-10**6, 10**6), _rr.randint(1, 10**4))
    vx1, vb, vtR, vtl = rq(), rq(), rq(), rq()
    Pm = [rq() for _ in range(4)]
    Um = [rq() for _ in range(4)]
    Vm = [rq() for _ in range(4)]
    gRv = lambda z: z / (2 * (1 + vtR * z))
    glv = lambda z: z / (2 * (1 + vtl * z))
    kav = lambda z: z * z / (2 * (1 + vtR * z) * (1 + vtl * z))
    lhs = det_exact([[Um[m] / vx1, Vm[m] / vb, gRv(Pm[m]) * (Um[m] - Vm[m]), -glv(Pm[m]) * (Um[m] - Vm[m])]
                     for m in range(4)])
    Wv = [Um[m] - Vm[m] for m in range(4)]
    Dt = det_exact([[Fr(1), Vm[m] / Wv[m], gRv(Pm[m]), kav(Pm[m])] for m in range(4)])
    rhs = -(vtR - vtl) / (vx1 * vb) * Wv[0] * Wv[1] * Wv[2] * Wv[3] * Dt
    nz += (lhs != rhs)
log(f"(B1) det[col_k(p_m)] = -(tauR-taul)/(xi1 b) prod_m W(p_m) det[1, V/W, gR, kappa](p_m): exact rational "
    f"evaluation at {ntest} random points (16 free quantities): violated in {nz} cases "
    f"(the symbolic identity is also checked in 02_sigma.wl)")
# real / complex determinant relation, generic complex rows
a = sp.symbols("a1:9", real=True)
c = sp.symbols("c1:9", real=True)
r1v = [a[k] + sp.I * c[k] for k in range(4)]
r2v = [a[4 + k] + sp.I * c[4 + k] for k in range(4)]
Mr = sp.Matrix([[sp.re(x) for x in r1v], [sp.im(x) for x in r1v], [sp.re(x) for x in r2v], [sp.im(x) for x in r2v]])
Mc = sp.Matrix([r1v, [sp.conjugate(x) for x in r1v], r2v, [sp.conjugate(x) for x in r2v]])
log("(B2) det[Re r1; Im r1; Re r2; Im r2] + det[r1; conj r1; r2; conj r2]/4 =",
    sp.expand(Mr.det() + Mc.det() / 4))

# ---------------------------------------------------------------- (C)
# The t-coefficient e_1 = -(c2/b - h)(1 + w^2): 1 + w^2 -> 1 + b^2/bs^2 at p -> 0, and (with R)
w_full = B_e * (R + 1 / (bs * sp.sqrt(p)))
lim = sp.limit((w_full**2).subs({p: qq**2}), qq, 0)
log("(C) lim_{p->0} w^2 =", sp.simplify(lim), " (so 1 + w^2 is not identically zero)")
hcoef = sp.Poly(sp.expand(c3 * (1 + tl * p) - c4 * (1 + tR * p)), p).all_coeffs()
hsys = sp.Matrix([[sp.diff(e, c3), sp.diff(e, c4)] for e in hcoef])
log("    h(p) = c3 gR - c4 gl vanishes at p = 0; h = 0 identically iff c3 (1+tau_l p) = c4 (1+tau_R p),")
log("    i.e. the linear system", hcoef, "= 0, with determinant", sp.factor(hsys.det()),
    ": only c3 = c4 = 0 unless tau_R = tau_l, where (c3, c4) = (1, 1) spans the solutions;",
    "null space at tau_l = tau_R:", hsys.subs(tl, tR).nullspace())
# rho^2 (1 - w^2) contains an odd power of sqrt(p) when R > 0:
expr = sp.expand(((1 + tR * qq**2) / (1 + tl * qq**2)) - b**2 * (R * qq + 1 / bs) ** 2)
log("(C) rho^2 - b^2 (R sqrt p + 1/bs)^2, odd part in q = sqrt p:",
    sp.simplify((expr - expr.subs(qq, -qq)) / 2))

(HERE / "results" / "02_sigma_symbolic.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
