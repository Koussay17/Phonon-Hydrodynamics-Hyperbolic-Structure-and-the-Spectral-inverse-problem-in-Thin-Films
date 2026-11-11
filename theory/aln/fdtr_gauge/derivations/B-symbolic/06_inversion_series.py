"""By-product of Task 4: constructive inversion from the small-p Puiseux coefficients (claim B2).

Z(i omega) determines its expansion Z = sum_{n >= -1} z_n q^n, q = sqrt(p) (convergent near q = 0),
hence z_0 ... z_4.  With delta = tau_R - tau_l:
  z_0 = xi1 (1/b - b/bs^2)                         (no tau)
  z_1 = -(xi1 b / bs) z_0                          (no tau)
  z_2 = xi1 delta / b + F2(xi1, b)
  z_3 = -xi1^2 delta / bs + F3(xi1, b)
  z_4 = -xi1 tau_l delta / b + G4(xi1, b, delta)
Checks: (I1) these structural statements (SymPy, exact); (I2) the inversion formulas applied to the
symbolic coefficients return theta (exact simplification) in the two cases b != bs and b = bs;
(I3) exact Jacobian determinants of coefficient subsets, factored, including at b = bs;
(I4) end-to-end numerical inversion: z_n extracted from Z(i omega) values alone (Cauchy integral on a
small circle in q, 50 digits) and inverted, at random points.
Output: results/06_inversion_series.txt
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


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


q = sp.Symbol("q")
xi1, b, bs = sp.symbols("xi1 b b_s", positive=True)
tR, tl, dl = sp.symbols("tau_R tau_l delta")
rho = sp.sqrt((1 + tR * q**2) / (1 + tl * q**2))
A = xi1 * q * rho
B = b * q / rho
w = B / (bs * q)
Z = (w + sp.tanh(A)) / (B * (1 + w * sp.tanh(A)))
ser = sp.expand(sp.series(Z * q, q, 0, 7).removeO())
z = {n - 1: sp.expand(sp.simplify(ser.coeff(q, n))) for n in range(7)}
zd = {n: sp.expand(z[n].subs(tR, tl + dl)) for n in z}
log("(I1) structure with tau_R = tau_l + delta:")
log("   d z_0/d(tau_l, delta) =", [sp.simplify(sp.diff(zd[0], v)) for v in (tl, dl)])
log("   d z_1/d(tau_l, delta) =", [sp.simplify(sp.diff(zd[1], v)) for v in (tl, dl)], ";  z_1 + (xi1 b/bs) z_0 =",
    sp.simplify(zd[1] + xi1 * b / bs * zd[0]))
log("   d z_2/d tau_l =", sp.simplify(sp.diff(zd[2], tl)), ";  d z_2/d delta =", sp.simplify(sp.diff(zd[2], dl)))
log("   d z_3/d tau_l =", sp.simplify(sp.diff(zd[3], tl)), ";  d z_3/d delta =", sp.simplify(sp.diff(zd[3], dl)))
log("   d z_4/d tau_l =", sp.factor(sp.diff(zd[4], tl)), ";  d^2 z_4/d tau_l^2 =", sp.simplify(sp.diff(zd[4], tl, 2)))
F2 = sp.simplify(zd[2] - xi1 * dl / b)
F3 = sp.simplify(zd[3] + xi1**2 * dl / bs)
G4 = sp.simplify(zd[4] + xi1 * tl * dl / b)
log("   F2 =", sp.factor(F2), "  (free of tau_l, delta:", sp.simplify(sp.diff(F2, tl)) == 0 and sp.simplify(sp.diff(F2, dl)) == 0, ")")
log("   F3 =", sp.factor(F3), "  (free of tau_l, delta:", sp.simplify(sp.diff(F3, tl)) == 0 and sp.simplify(sp.diff(F3, dl)) == 0, ")")
log("   G4 =", sp.factor(G4), "  (free of tau_l:", sp.simplify(sp.diff(G4, tl)) == 0, ")")

# (I2) inversion formulas
Z0, Z1, Z2, Z3, Z4 = sp.symbols("Z0 Z1 Z2 Z3 Z4")
Cf_inv = -bs * Z1 / Z0
Rf_inv = Z0 + Cf_inv / bs**2
xi1_inv = sp.sqrt(Cf_inv * Rf_inv)
b_inv = sp.sqrt(Cf_inv / Rf_inv)
delta_inv = (Z2 - F2.subs({xi1: xi1_inv, b: b_inv})) * b_inv / xi1_inv
taul_inv = -(Z4 - G4.subs({xi1: xi1_inv, b: b_inv, dl: delta_inv})) * b_inv / (xi1_inv * delta_inv)
subsZ = {Z0: zd[0], Z1: zd[1], Z2: zd[2], Z3: zd[3], Z4: zd[4]}
# exact check at random rational points (b != bs), and symbolic simplification of the first steps
rng = random.Random(8)
bad = 0
for _ in range(20):
    vals = {xi1: sp.Rational(rng.randint(1, 300), 100), b: sp.Rational(rng.randint(1, 300), 100),
            bs: sp.Rational(rng.randint(1, 300), 100), tl: sp.Rational(rng.randint(0, 300), 100),
            dl: sp.Rational(rng.randint(-300, 300), 100)}
    if vals[b] == vals[bs] or vals[dl] == 0:
        continue
    zv = {k: v.subs(vals) for k, v in subsZ.items()}
    rec = [sp.nsimplify(sp.simplify(e.subs(zv).subs(bs, vals[bs]))) for e in (xi1_inv, b_inv, delta_inv, taul_inv)]
    bad += rec != [vals[xi1], vals[b], vals[dl], vals[tl]]
log(f"(I2) case b != bs: inversion (C_f, R_f) -> (xi1, b) -> delta -> tau_l reproduces theta exactly at random rational points: failures {bad}")
log("     symbolic: C_f(z) - xi1 b =", sp.simplify(Cf_inv.subs(subsZ) - xi1 * b), ";  R_f(z) - xi1/b =",
    sp.simplify(Rf_inv.subs(subsZ) - xi1 / b))
# case b = bs
zb = {n: sp.expand(zd[n].subs(b, bs)) for n in zd}
log("(I2) case b = bs:  z_0 =", sp.simplify(zb[0]), ", z_1 =", sp.simplify(zb[1]), ", z_2 =", sp.factor(zb[2]),
    ", z_3 =", sp.factor(zb[3]), ", z_4 =", sp.factor(zb[4]))
log("     so xi1 = -z_3/z_2, delta = bs z_2/xi1, tau_l from z_4 (when delta != 0): check",
    sp.simplify(-zb[3] / zb[2] - xi1), sp.simplify(bs * zb[2] / xi1 - dl))

# (I3) exact Jacobian determinants of coefficient subsets
th = [xi1, b, tR, tl]
for ks in [(0, 1, 2, 4), (0, 2, 3, 4), (1, 2, 3, 4), (0, 1, 3, 4), (2, 3, 4, 5)]:
    Jc = sp.Matrix([[sp.diff(z[k], v) for v in th] for k in ks])
    d = sp.factor(sp.simplify(Jc.det()))
    log(f"(I3) det d(z_{ks})/d(xi1, b, tau_R, tau_l) = {d}")
    log(f"       at b = bs: {sp.factor(sp.simplify(d.subs(b, bs)))}")

# (I4) end-to-end: coefficients from Z values on a circle in q (Cauchy), then inversion
mp.mp.dps = 50
log("(I4) end-to-end numerical inversion (z_n from Z on |q| = r by trapezoidal Cauchy integrals, 50 digits):")
for _ in range(6):
    tv = [mp.mpf(rng.uniform(0.2, 2)) for _ in range(3)]
    x1v, bv, bsv = tv
    if _ == 5:
        bsv = bv
    tlv, tRv = mp.mpf(rng.uniform(0, 1)), mp.mpf(rng.uniform(0, 1))
    r, K = mp.mpf("0.05") / x1v, 64

    def Zq(qq):
        # branch sqrt(p) = q (the analytic continuation in q through q = 0), rho principal near 1
        rh = mp.sqrt((1 + tRv * qq**2) / (1 + tlv * qq**2))
        A_, B_ = x1v * qq * rh, bv * qq / rh
        w_ = B_ / (bsv * qq)
        t_ = mp.tanh(A_)
        return (w_ + t_) / (B_ * (1 + w_ * t_))

    # consistency with the principal-branch implementation on the right half q-plane
    q0 = mp.mpf("0.03") * mp.expj(mp.pi / 4)
    assert abs(Zq(q0) - M.Z_mp(q0**2, x1v, bv, tRv, tlv, bsv)) < mp.mpf("1e-40") * abs(Zq(q0))
    coeffs = {}
    for n in range(-1, 6):
        acc = 0
        for k in range(K):
            qk = r * mp.expj(2 * mp.pi * k / K)
            acc += Zq(qk) * qk ** (-n)
        coeffs[n] = acc / K
    c = {n: mp.re(coeffs[n]) for n in coeffs}
    if abs(c[0]) > mp.mpf("1e-30"):
        Cf_ = -bsv * c[1] / c[0]
        Rf_ = c[0] + Cf_ / bsv**2
        x1r, br = mp.sqrt(Cf_ * Rf_), mp.sqrt(Cf_ / Rf_)
    else:
        x1r, br = -c[3] / c[2], bsv
    F2n = sp.lambdify((xi1, b, bs), F2, modules="mpmath")
    dr = (c[2] - F2n(x1r, br, bsv)) * br / x1r
    G4n = sp.lambdify((xi1, b, bs, dl), G4, modules="mpmath")
    tlr = -(c[4] - G4n(x1r, br, bsv, dr)) * br / (x1r * dr)
    err = max(abs(x1r - x1v), abs(br - bv), abs(dr - (tRv - tlv)), abs(tlr - tlv))
    log(f"   true (xi1, b, tau_R, tau_l, bs) = ({mp.nstr(x1v, 6)}, {mp.nstr(bv, 6)}, {mp.nstr(tRv, 6)}, {mp.nstr(tlv, 6)}, {mp.nstr(bsv, 6)})"
        f" -> max reconstruction error {mp.nstr(err, 3)}")
(HERE / "results" / "06_inversion_series.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
