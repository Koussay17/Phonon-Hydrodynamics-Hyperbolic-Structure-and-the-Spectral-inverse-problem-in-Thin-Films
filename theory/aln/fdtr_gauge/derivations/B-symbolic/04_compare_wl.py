"""Exact comparison of the SymPy expansion coefficients (04_expansions.py) with the Wolfram ones
(results/04_expansions_wl.txt) at the test point xi1 = 7/10, b = 23/10, bs = 11/10, tau_R = 3/10,
tau_l = 1/8, R = 1/5.  Small-p coefficients are rational there and are compared exactly; the 1/B and
tau_R = 0 coefficients are compared at 30 digits.  Output: results/04_compare_wl.txt"""
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
q, u = sp.symbols("q u", positive=True)
xi1, b, bs, R, tR, tl, T0 = sp.symbols("xi1 b b_s R tau_R tau_l T0", positive=True)
pt = {xi1: sp.Rational(7, 10), b: sp.Rational(23, 10), bs: sp.Rational(11, 10), tR: sp.Rational(3, 10),
      tl: sp.Rational(1, 8), R: sp.Rational(1, 5)}
wl0 = [sp.Rational(10, 11), sp.Rational(-2856, 2783), sp.Rational(9996, 6655), sp.Rational(-133184891, 67348600),
       sp.Rational(34130901, 12884080), sp.Rational(-86473015525463, 24447541800000),
       sp.Rational(183545921959029, 38974342000000)]
wlR = [sp.Rational(10, 11), sp.Rational(-11497, 13915), sp.Rational(30499, 33275), sp.Rational(-1438913301, 1683715000),
       sp.Rational(19470795673, 24157650000), sp.Rational(-456989146819289, 611188545000000),
       sp.Rational(664403087684609, 974358550000000)]
wlN = [sp.Float("1.27755769170825977337981950167881818950", 35), sp.Float("-0.13185295315430055849709208082329289643", 35),
       sp.Float("-4.15299692935073592154135972082225511917", 35), sp.Float("-1.37520639061853173738441604842824718841", 35)]
out = []
# small p: substitute the numbers first (exact rationals), then series
for withR, ref in ((False, wl0), (True, wlR)):
    rho = sp.sqrt((1 + pt[tR] * q**2) / (1 + pt[tl] * q**2))
    A = pt[xi1] * q * rho
    B = pt[b] * q / rho
    Zs = (pt[R] if withR else 0) + 1 / (pt[bs] * q)
    w = B * Zs
    Z = (w + sp.tanh(A)) / (B * (1 + w * sp.tanh(A)))
    ser = sp.expand(sp.series(Z * q, q, 0, 7).removeO())
    mine = [sp.nsimplify(ser.coeff(q, n)) for n in range(7)]
    out.append(f"small p {'with R' if withR else 'R = 0'}: SymPy z_-1..z_5 at the point = {mine}")
    out.append(f"   exact agreement with Wolfram: {[m - r == 0 for m, r in zip(mine, ref)]}")
# tau_R = 0 branch, numerically at 30 digits
uN = sp.Symbol("uN", positive=True)
B_N = pt[b] * sp.sqrt(1 + pt[tl] / uN**2) / uN
A_N = pt[xi1] / sp.sqrt(pt[tl] + uN**2)
Zn = (B_N * uN / pt[bs] + sp.tanh(A_N)) / (B_N * (1 + B_N * uN / pt[bs] * sp.tanh(A_N)))
serN = sp.series(Zn, uN, 0, 6).removeO()
mineN = [sp.N(sp.expand(serN).coeff(uN, k), 32) for k in range(2, 6)]
out.append(f"tau_R = 0 branch, coefficients of p^(-k/2), k=2..5, SymPy: {mineN}")
out.append(f"   max |SymPy - Wolfram| = {max(abs(m - r) for m, r in zip(mineN, wlN))}")
(HERE / "results" / "04_compare_wl.txt").write_text("\n".join(map(str, out)) + "\n", encoding="utf-8")
print("\n".join(map(str, out)))
