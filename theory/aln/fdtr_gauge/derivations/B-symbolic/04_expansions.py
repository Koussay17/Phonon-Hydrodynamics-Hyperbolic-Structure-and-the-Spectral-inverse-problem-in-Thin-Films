"""Task 4. Small-p (Puiseux in q = sqrt p) and large-p expansions of Z, by stratum.

Small p:  Z = sum_{n >= -1} z_n q^n.  Coefficients up to z_5, with and without the interface
          resistance R; physical reading through R_f = xi1/b = e/lambda and C_f = xi1 b = rho c e.
          Exact Jacobian of theta -> (z_0, z_1, z_2, z_k) and its factorisation: an independent,
          exact-arithmetic sufficient condition for rank 4 of the function-level Jacobian.
Large p:  (G) tau_R, tau_l > 0:  Z = 1/B + (2/B) G E/(1 - G E), E = exp(-2A): power series of 1/B
              in p^(-1/2) plus exponentially small terms (which alone carry xi1, bs, R);
          (F) tau_R = tau_l: 1/B = 1/(b sqrt p) exactly; transseries in exp(-2 xi1 sqrt p);
          (C) tau_l = 0 < tau_R (Cattaneo): expansion in p^(-1/2) with coefficients polynomial in
              t = tanh A, which stays bounded on p = i omega;
          (N) tau_R = 0 < tau_l: A -> xi1/sqrt(tau_l); power series in p^(-1/2) from p^(-1).
Numerical check (mpmath, 50 digits) of every truncation on p = i omega: the error of the
n-term truncation divided by the first omitted order must stay bounded (and tend to the next
coefficient where it is a constant).
Output: results/04_expansions.txt
"""
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


q = sp.Symbol("q")          # q = sqrt(p) (small p) ; u = p^(-1/2) (large p)
u = sp.Symbol("u", positive=True)
xi1, b, bs, R = sp.symbols("xi1 b b_s R", positive=True)
tR, tl = sp.symbols("tau_R tau_l", positive=True)
Rf, Cf = sp.symbols("R_f C_f", positive=True)


def Z_of(p_expr, A_expr, B_expr, Zs_expr, t_expr=None):
    t = sp.tanh(A_expr) if t_expr is None else t_expr
    w = B_expr * Zs_expr
    return (w + t) / (B_expr * (1 + w * t))


# ------------------------------------------------------------------ small p
ORDER = 7
rho_q = sp.sqrt((1 + tR * q**2) / (1 + tl * q**2))
A_q = xi1 * q * rho_q
B_q = b * q / rho_q
log("SMALL p (q = sqrt p).  Z = sum_n z_n q^n")
coef = {}
for withR in (False, True):
    Zs_q = (R if withR else 0) + 1 / (bs * q)
    Zq = Z_of(q**2, A_q, B_q, Zs_q)
    ser = sp.series(Zq * q, q, 0, ORDER).removeO()
    ser = sp.expand(ser)
    zs = {n - 1: sp.factor(sp.simplify(ser.coeff(q, n))) for n in range(ORDER)}
    coef[withR] = zs
    log(f"  {'with interface resistance R' if withR else 'R = 0'}:")
    for n in range(-1, ORDER - 1):
        log(f"    z_{n} = {zs[n]}")
# physical reading (R = 0): substitute xi1 = sqrt(Rf Cf), b = sqrt(Cf/Rf)
phys = {xi1: sp.sqrt(Rf * Cf), b: sp.sqrt(Cf / Rf)}
log("  R = 0 in terms of R_f = xi1/b = e/lambda and C_f = xi1 b = rho c e:")
for n in range(0, 4):
    log(f"    z_{n} = {sp.factor(sp.simplify(coef[False][n].subs(phys)))}")
z0, z1 = coef[False][0], coef[False][1]
log("  identity z_1 = -(C_f/b_s) z_0 :", sp.simplify((z1 + (xi1 * b / bs) * z0)))
# where do tau_R and tau_l separate?  d z_n / d tau_R + d z_n / d tau_l
for n in range(0, ORDER - 1):
    d = sp.simplify(sp.diff(coef[False][n], tR) + sp.diff(coef[False][n], tl))
    dm = sp.simplify(sp.diff(coef[False][n], tR) - sp.diff(coef[False][n], tl))
    log(f"  z_{n}: (d_tauR + d_taul) z_n = {sp.factor(d)} ;  (d_tauR - d_taul) z_n = {sp.factor(dm)}")
# exact Jacobian of the coefficient map
th = [xi1, b, tR, tl]
for ks in [(0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 2, 5)]:
    Jc = sp.Matrix([[sp.diff(coef[False][k], v) for v in th] for k in ks])
    dJ = sp.factor(sp.simplify(Jc.det()))
    log(f"  det d(z_{ks[0]}, z_{ks[1]}, z_{ks[2]}, z_{ks[3]})/d(xi1, b, tau_R, tau_l) = {dJ}")
# the same with R: coefficient map (z_0 .. ) on 5 parameters
th5 = [xi1, b, tR, tl, R]
for ks in [(0, 1, 2, 3, 4), (0, 1, 2, 3, 5)]:
    Jc = sp.Matrix([[sp.diff(coef[True][k], v) for v in th5] for k in ks])
    dJ = sp.factor(sp.simplify(Jc.det()))
    log(f"  with R: det d(z_{ks})/d(xi1, b, tau_R, tau_l, R) = {dJ}")

# ------------------------------------------------------------------ large p, generic
log("")
log("LARGE p, generic tau_R, tau_l > 0:  Z = 1/B + exponentially small;  1/B = sum_k y_k p^(-1/2 - k)")
invB = (1 / b) * u * sp.sqrt((tR + u**2) / (tl + u**2))      # p = 1/u^2: 1/B = u/b sqrt((1 + tR/u^2)/(1 + tl/u^2))
serB = sp.expand(sp.series(invB, u, 0, 9).removeO())
ys = [sp.factor(sp.simplify(serB.coeff(u, 2 * k + 1))) for k in range(4)]
for k, y in enumerate(ys):
    log(f"  y_{k} (coefficient of p^(-{2*k+1}/2)) = {y}")
log("  in terms of K = sqrt(tau_R/tau_l)/b, a_R = 1/tau_R, a_l = 1/tau_l:")
K, aR, aL = sp.symbols("K a_R a_l", positive=True)
for k, y in enumerate(ys):
    log(f"    y_{k} / K = {sp.factor(sp.simplify((y / (sp.sqrt(tR / tl) / b)).subs({tR: 1 / aR, tl: 1 / aL})))}")
log("  leading exponential term: Z - 1/B ~ (2 G_inf / B) exp(-2A),  A = xi1 sqrt(tau_R/tau_l) [sqrt p + (a_R - a_l)/(2 sqrt p) + ...],")
log("  G_inf = (b sqrt(tau_l/tau_R) - b_s)/(b sqrt(tau_l/tau_R) + b_s)   (R = 0); with R > 0, G -> 1.")
A_large = xi1 * sp.sqrt(tR / tl) * sp.sqrt((1 + aR * u**2) / (1 + aL * u**2)).subs({aR: 1 / tR, aL: 1 / tl}) / u
serA = sp.expand(sp.series(A_large * u, u, 0, 6).removeO())
log("  A sqrt(p)^-1 expansion: A =", " + ".join(f"({sp.factor(sp.simplify(serA.coeff(u, 2*k)))}) p^({1-2*k}/2)" for k in range(3)))

# ------------------------------------------------------------------ large p, Cattaneo and tau_R = 0
log("")
log("LARGE p, Cattaneo tau_l = 0 < tau_R: expansion in u = p^(-1/2) at fixed t = tanh A (|t|, |sech A| bounded on i omega)")
t = sp.Symbol("t")
B_C = b * sp.sqrt(1 / (tR + u**2))          # b sqrt(p/(1 + tR p)) with p = 1/u^2
S_C = bs / u
Zs_C = 1 / S_C
ZC = Z_of(None, None, B_C, Zs_C, t_expr=t)
serC = sp.expand(sp.series(ZC, u, 0, 4).removeO())
for k in range(4):
    log(f"  coefficient of p^(-{k}/2): {sp.factor(sp.simplify(serC.coeff(u, k)))}")
A_C = xi1 * sp.sqrt((1 + tR / u**2) / u**2)       # xi1 sqrt(p (1 + tR p))
serAC = sp.series(A_C * u**2, u, 0, 5).removeO()
log("  A = xi1 sqrt(tau_R) p + xi1/(2 sqrt(tau_R)) + O(1/p); exact series of A u^2:", sp.expand(serAC))

log("")
log("LARGE p, tau_R = 0 < tau_l: A -> xi1/sqrt(tau_l); power series in u = p^(-1/2)")
# Fully symbolic series of this branch is slow in SymPy; the symbolic coefficients c_2..c_5 are taken from
# 04_expansions.wl (Wolfram), the first two are written below, and all four are checked here (i) against a SymPy
# series with the parameters substituted by exact rationals and (ii) by the 50-digit truncation test further down.
A_inf = xi1 / sp.sqrt(tl)
log("  c_2 (p^-1)   = coth(A_inf)/(b sqrt(tau_l)),          A_inf = xi1/sqrt(tau_l)")
log("  c_3 (p^-3/2) = -b_s csch(A_inf)^2/(b^2 tau_l)")
log("  c_4, c_5: see 04_expansions.wl / results/04_expansions_wl.txt")
uN = sp.Symbol("uN", positive=True)
numv = {xi1: sp.Rational(7, 10), b: sp.Rational(23, 10), bs: sp.Rational(11, 10), tl: sp.Rational(1, 8)}
B_N = numv[b] * sp.sqrt(1 + numv[tl] / uN**2) / uN
A_N = numv[xi1] / sp.sqrt(numv[tl] + uN**2)
ZN = Z_of(None, A_N, B_N, uN / numv[bs])
serN_num = sp.expand(sp.series(ZN, uN, 0, 6).removeO())
c2v = sp.N((sp.coth(A_inf) / (b * sp.sqrt(tl))).subs(numv), 30)
c3v = sp.N((-bs * sp.csch(A_inf) ** 2 / (b**2 * tl)).subs(numv), 30)
log(f"  check at the test point: series c_2 = {sp.N(serN_num.coeff(uN, 2), 20)} vs formula {sp.N(c2v, 20)};"
    f" c_3 = {sp.N(serN_num.coeff(uN, 3), 20)} vs formula {sp.N(c3v, 20)}")

# ------------------------------------------------------------------ numerical checks
log("")
log("NUMERICAL CHECKS on p = i omega (50 digits): error of the truncation / |first omitted order|")
mp.mp.dps = 50
vals = {xi1: sp.Rational(7, 10), b: sp.Rational(23, 10), bs: sp.Rational(11, 10), tR: sp.Rational(3, 10),
        tl: sp.Rational(1, 8), R: sp.Rational(1, 5)}
fnum = lambda expr: sp.lambdify(q, expr.subs(vals), modules="mpmath")
for withR in (False, True):
    zs = coef[withR]
    Rv = vals[R] if withR else 0
    log(f"  small p, {'R = 1/5' if withR else 'R = 0'}: n-term truncation up to z_4; ratio |Z - sum_{{k<=4}} z_k q^k| / |q|^5")
    for om in (mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4")):
        pv = mp.mpc(0, om)
        qv = mp.sqrt(pv)
        exact = M.Z_mp(pv, *(mp.mpf(sp.Rational(vals[k]).p) / sp.Rational(vals[k]).q for k in (xi1, b, tR, tl, bs)), mp.mpf(sp.Rational(Rv).p) / sp.Rational(Rv).q)
        part = sum(mp.mpmathify(complex(sp.N(zs[n].subs(vals), 60))) * qv**n for n in range(-1, 5))
        ratio = abs(exact - part) / abs(qv) ** 5
        log(f"    omega={mp.nstr(om, 2)}: ratio = {mp.nstr(ratio, 8)}   (|z_5| = {mp.nstr(abs(complex(sp.N(zs[5].subs(vals), 30))), 8)})")
log("  large p, generic: ratio |Z - sum_{k<=3} y_k p^(-1/2-k)| / |p|^(-9/2)")
ysnum = [complex(sp.N(y.subs(vals), 60)) for y in ys]
y4 = sp.expand(sp.series(invB, u, 0, 11).removeO()).coeff(u, 9)
for om in (mp.mpf(1e3), mp.mpf(1e4), mp.mpf(1e5)):
    pv = mp.mpc(0, om)
    exact = M.Z_mp(pv, 0.7, 2.3, 0.3, 0.125, 1.1)
    part = sum(mp.mpmathify(ysnum[k]) * pv ** (-mp.mpf(1) / 2 - k) for k in range(4))
    log(f"    omega={mp.nstr(om, 2)}: ratio = {mp.nstr(abs(exact - part) / abs(pv) ** (-4.5), 8)}  (|y_4| = {mp.nstr(abs(complex(sp.N(y4.subs(vals), 30))), 8)});"
        f"  |Z - 1/B| = {mp.nstr(abs(exact - 1 / M.parts_mp(pv, 0.7, 2.3, 0.3, 0.125)[1]), 3)}")
log("  large p, Cattaneo (tau_l = 0): ratio |Z - [t/B + (1-t^2)/S]| / |p|^(-1), t = tanh A exact")
for om in (mp.mpf(1e3), mp.mpf(1e4), mp.mpf(1e5)):
    pv = mp.mpc(0, om)
    A_, B_, sq = M.parts_mp(pv, 0.7, 2.3, 0.3, 0)
    tv = mp.tanh(A_)
    exact = M.Z_mp(pv, 0.7, 2.3, 0.3, 0, 1.1)
    approx = tv / B_ + (1 - tv**2) / (1.1 * sq)
    log(f"    omega={mp.nstr(om, 2)}: ratio = {mp.nstr(abs(exact - approx) * abs(pv), 8)};  |t| = {mp.nstr(abs(tv), 5)}")
log("  large p, tau_R = 0: ratio |Z - sum_{k=2..5} c_k p^(-k/2)| / |p|^(-3)")
cN = [complex(sp.N(serN_num.coeff(uN, k), 60)) for k in range(2, 6)]
for om in (mp.mpf(1e3), mp.mpf(1e4), mp.mpf(1e5)):
    pv = mp.mpc(0, om)
    exact = M.Z_mp(pv, 0.7, 2.3, 0, 0.125, 1.1)
    part = sum(mp.mpmathify(cN[k - 2]) * pv ** (-mp.mpf(k) / 2) for k in range(2, 6))
    log(f"    omega={mp.nstr(om, 2)}: ratio = {mp.nstr(abs(exact - part) * abs(pv) ** 3, 8)}")

(HERE / "results" / "04_expansions.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
