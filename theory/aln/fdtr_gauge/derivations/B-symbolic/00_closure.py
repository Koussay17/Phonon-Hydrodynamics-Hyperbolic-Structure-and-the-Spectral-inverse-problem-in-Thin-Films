"""Task 0. Front-face impedance from the 1D quadrupole closed by the substrate.

Checks
  (S1) SymPy: (M00 Zs + M01)/(M10 Zs + M11) with M = [[cosh A, sinh A/B],[B sinh A, cosh A]]
       equals (1/B)(w + t)/(1 + w t), w = B Zs, t = tanh A; with Zs = 1/S equals
       (1/B)(B + S t)/(S + B t).
  (S2) SymPy: the brief's expression Zb = (1/B)(B t + S)/(B + S t) satisfies Zb = 1/(B^2 Z).
  (S3) limits A -> 0 (thin film: must give the bare substrate 1/S) and t -> 1 (thick film: 1/B).
  (S4) multiple-reflection form Z = (1/B)(1 + G E)/(1 - G E), G = (w-1)/(w+1), E = exp(-2A).
  (N1) mpmath, 50 digits, random rational points on p = i omega: quadrupole product vs closed form
       vs the brief's expression.
  (N2) the project's own forward model (src/forward_model.py, double precision) vs the closed form
       and vs the brief's expression.
Output: results/00_closure.txt
"""
import random
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

out_lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out_lines.append(s)


# ----------------------------------------------------------------------------- S1-S4
A, B, S, Zs = sp.symbols("A B S Z_s")          # complex symbols, no assumptions
ch, sh, t = sp.cosh(A), sp.sinh(A), sp.tanh(A)
quad = (ch * Zs + sh / B) / (B * sh * Zs + ch)
w = B * Zs
closed = (w + t) / (B * (1 + w * t))
d1 = sp.simplify((quad - closed).rewrite(sp.exp))
log("S1a  quad - (1/B)(w+t)/(1+wt)  ->", d1)
quadS = quad.subs(Zs, 1 / S)
closedS = (B + S * t) / (B * (S + B * t))
d2 = sp.simplify((quadS - closedS).rewrite(sp.exp))
log("S1b  quad(Zs=1/S) - (1/B)(B+St)/(S+Bt)  ->", d2)

brief = (B * t + S) / (B * (B + S * t))
d3 = sp.simplify((brief * closedS * B ** 2 - 1).rewrite(sp.exp))
log("S2   Zbrief * Z * B^2 - 1  ->", d3)
d3b = sp.simplify((brief - closedS).rewrite(sp.exp))
log("S2b  Zbrief - Z  ->", sp.factor(sp.simplify(sp.together(brief - closedS))))

# limits
tt = sp.Symbol("t")
Zt = (B + S * tt) / (B * (S + B * tt))
Zbt = (B * tt + S) / (B * (B + S * tt))
log("S3   t->0 (thin film):   Z ->", sp.simplify(Zt.subs(tt, 0)), "   Zbrief ->", sp.simplify(Zbt.subs(tt, 0)))
log("S3   t->1 (thick film):  Z ->", sp.simplify(Zt.subs(tt, 1)), "   Zbrief ->", sp.simplify(Zbt.subs(tt, 1)))
log("S3   the bare-substrate value 1/S is reproduced by Z only.")

E = sp.Symbol("E")
G = (w - 1) / (w + 1)
refl = (1 + G * E) / (B * (1 - G * E))
d4 = sp.simplify(sp.together(refl.subs(E, sp.exp(-2 * A)) - closed).rewrite(sp.exp))
log("S4   (1/B)(1+G e^{-2A})/(1-G e^{-2A}) - Z  ->", d4)

# ----------------------------------------------------------------------------- N1
mp.mp.dps = 50
rng = random.Random(20261111)


def rnd_log(lo, hi, den=1000):
    u = rng.uniform(lo, hi)
    return Fraction(max(1, round(10 ** u * den)), den)


log("")
log("N1  50-digit comparison on p = i omega (random rational points)")
log("    columns: xi1, b, tau_R, tau_l, bs, omega | |Zquad-Zclosed|/|Z| | |Zbrief-Z|/|Z|")
maxd_q, mind_b = mp.mpf(0), mp.mpf(10)
for k in range(40):
    xi1 = rnd_log(-1, 0.7)
    b, bs = rnd_log(-1, 1), rnd_log(-1, 1)
    T = xi1 ** 2
    tR = rnd_log(-2, 1) * T if k % 5 else Fraction(0)
    tl = rnd_log(-2, 1) * T if k % 7 else Fraction(0)
    om = rnd_log(-2, 2.5) / T
    p = mp.mpc(0, M._mpf(om))
    zq = M.Z_quad_mp(p, xi1, b, tR, tl, bs)
    zc = M.Z_mp(p, xi1, b, tR, tl, bs)
    zb = M.Z_brief_mp(p, xi1, b, tR, tl, bs)
    dq = abs(zq - zc) / abs(zc)
    db = abs(zb - zc) / abs(zc)
    maxd_q = max(maxd_q, dq)
    mind_b = min(mind_b, db)
    if k < 10:
        log(f"    {float(xi1):.4g} {float(b):.4g} {float(tR):.4g} {float(tl):.4g} {float(bs):.4g} {float(om):.4g} |"
            f" {mp.nstr(dq, 3)} | {mp.nstr(db, 3)}")
log(f"    max over 40 points: quad vs closed {mp.nstr(maxd_q, 3)};  min over 40 points: brief vs closed {mp.nstr(mind_b, 3)}")

# with interface resistance: quadrupole product with [[1,R],[0,1]] vs closed form with Zs = R + 1/S
maxd_r = mp.mpf(0)
for k in range(20):
    xi1 = rnd_log(-1, 0.7)
    b, bs = rnd_log(-1, 1), rnd_log(-1, 1)
    T = xi1 ** 2
    tR, tl = rnd_log(-2, 1) * T, rnd_log(-2, 1) * T
    R = rnd_log(-2, 1) * xi1 / b
    om = rnd_log(-2, 2.5) / T
    p = mp.mpc(0, M._mpf(om))
    zq = M.Z_quad_mp(p, xi1, b, tR, tl, bs, R)
    zc = M.Z_mp(p, xi1, b, tR, tl, bs, R)
    maxd_r = max(maxd_r, abs(zq - zc) / abs(zc))
log(f"    with interface resistance (20 points): max |Zquad - Zclosed|/|Z| = {mp.nstr(maxd_r, 3)}")

# ----------------------------------------------------------------------------- N2
log("")
log("N2  project forward model (src/forward_model.py, double precision), AlN/sapphire-like values")
sys.path.insert(0, r"C:\Users\Koussay\these\src")
try:
    import numpy as np
    import forward_model as fm
    worst_c, best_b = 0.0, 1e9
    for lam, rc, e, tR, tl, lams, rcs, Rc in [
        (321.0, 2.41e6, 500e-9, 1e-9, 0.0, 35.0, 3.03e6, 0.0),
        (60.0, 2.41e6, 500e-9, 0.5e-9, 0.2e-9, 35.0, 3.03e6, 0.0),
        (150.0, 2.41e6, 100e-9, 2e-9, 5e-9, 35.0, 3.03e6, 1e-8),
        (321.0, 2.41e6, 300e-9, 0.0, 1e-9, 1.4, 1.6e6, 2e-8),
        (321.0, 2.41e6, 500e-9, 0.0, 0.0, 35.0, 3.03e6, 0.0),
    ]:
        s = fm.Sample(film_lam=lam, film_rho_c=rc, thickness=e, sub_lam=lams, sub_rho_c=rcs,
                      relaxation_time=tR, nonlocal_time=tl, contact_resistance=Rc)
        xi1 = e / np.sqrt(lam / rc)
        b = np.sqrt(lam * rc)
        bs = np.sqrt(lams * rcs)
        for f in [1e4, 1e5, 1e6, 1e7, 1e8, 2e8]:
            p = 2j * np.pi * f
            zr = complex(fm.response(np.array([p]), s)[0])
            zc = complex(M.Z_mp(mp.mpc(0, 2 * mp.pi * f), xi1, b, tR, tl, bs, Rc))
            worst_c = max(worst_c, abs(zr - zc) / abs(zc))
            if Rc == 0:
                zb = complex(M.Z_brief_mp(mp.mpc(0, 2 * mp.pi * f), xi1, b, tR, tl, bs))
                best_b = min(best_b, abs(zr - zb) / abs(zr))
        log(f"    lam={lam}, e={e}, tR={tR}, tl={tl}, R={Rc}: running max rel diff (repo vs closed) = {worst_c:.2e}")
    log(f"    max rel diff repo vs closed form: {worst_c:.2e};  min rel diff repo vs brief (R=0 cases): {best_b:.2e}")
except Exception as exc:  # pragma: no cover
    log("    project model not importable:", repr(exc))

# thin-film physics check: e -> 0 must give the bare substrate
log("")
log("N3  thin-film limit (xi1 = 1e-8, Fourier film, b = 2, bs = 1, omega = 1):")
p = mp.mpc(0, 1)
for name, f in [("closed", M.Z_mp), ("brief", M.Z_brief_mp)]:
    z = f(p, mp.mpf("1e-8"), 2, 0, 0, 1)
    log(f"    {name:7s} Z = {mp.nstr(z, 12)}   bare substrate 1/(bs sqrt p) = {mp.nstr(1/mp.sqrt(p), 12)}"
        f"   S/B^2 = {mp.nstr(mp.sqrt(p)/(4*p), 12)}")

(HERE / "results").mkdir(exist_ok=True)
(HERE / "results" / "00_closure.txt").write_text("\n".join(out_lines) + "\n", encoding="utf-8")
