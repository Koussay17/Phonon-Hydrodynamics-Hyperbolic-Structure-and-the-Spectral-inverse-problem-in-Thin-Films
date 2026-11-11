"""Task 1. Rank of the Jacobian of theta -> (Z(i omega_j))_j, theta = (xi1, b, tau_R, tau_l).

Two independent derivative computations at 50 significant digits:
  (a) SymPy symbolic derivatives of the closed form, evaluated with mpmath;
  (b) mpmath numerical differentiation (100-digit working precision) of the direct
      quadrupole matrix product (no symbolic input).
Then singular values (mpmath SVD) of the column-equilibrated real Jacobian
[Re; Im] (2N x 4) and of the complex Jacobian (N x 4), N = 1, 2, 3, 4, 8.
The rational points are written to results/01_points.json for the ball-arithmetic
certification (01_rank_ball.py) and for the Wolfram check (01_rank.wl).
Output: results/01_rank_mp.txt
"""
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

mp.mp.dps = 50
rng = random.Random(4242)
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def rnd_log(lo, hi, den=1000):
    u = rng.uniform(lo, hi)
    return Fraction(max(1, round(10 ** u * den)), den)


def frac_str(x):
    return f"{x.numerator}/{x.denominator}"


# ------------------------------------------------------------------ parameter points
points = []
# generic random points (dimensionless family; bs random as well)
for k in range(24):
    xi1 = rnd_log(-1, 0.7)
    T = xi1 ** 2
    pt = dict(kind="random", xi1=xi1, b=rnd_log(-1, 1), bs=rnd_log(-1, 1),
              tau_R=rnd_log(-2, 1.3) * T, tau_l=rnd_log(-2, 1.3) * T)
    pt["omegas"] = sorted(rnd_log(-2, 2) / T for _ in range(8))
    points.append(pt)
# boundary strata of the domain, not on tau_R = tau_l: Cattaneo (tau_l = 0) and tau_R = 0
for k in range(4):
    xi1 = rnd_log(-1, 0.7)
    T = xi1 ** 2
    pt = dict(kind="cattaneo tau_l=0", xi1=xi1, b=rnd_log(-1, 1), bs=rnd_log(-1, 1),
              tau_R=rnd_log(-2, 1.3) * T, tau_l=Fraction(0))
    pt["omegas"] = sorted(rnd_log(-2, 2) / T for _ in range(8))
    points.append(pt)
for k in range(4):
    xi1 = rnd_log(-1, 0.7)
    T = xi1 ** 2
    pt = dict(kind="tau_R=0", xi1=xi1, b=rnd_log(-1, 1), bs=rnd_log(-1, 1),
              tau_R=Fraction(0), tau_l=rnd_log(-2, 1.3) * T)
    pt["omegas"] = sorted(rnd_log(-2, 2) / T for _ in range(8))
    points.append(pt)
# matched effusivity b = bs, off the diagonal
for k in range(4):
    xi1 = rnd_log(-1, 0.7)
    T = xi1 ** 2
    bb = rnd_log(-1, 1)
    pt = dict(kind="b=bs", xi1=xi1, b=bb, bs=bb,
              tau_R=rnd_log(-2, 1.3) * T, tau_l=rnd_log(-2, 1.3) * T)
    pt["omegas"] = sorted(rnd_log(-2, 2) / T for _ in range(8))
    points.append(pt)
# AlN / sapphire-like SI values (rational approximations): e = 500 nm,
# lambda = 321 or 60 W/mK, C = 2.41e6 J/m3K; sapphire 35 W/mK, 3.03e6 J/m3K.
for lam, tR_ns, tl_ns in [(321, 1, 0), (321, 1, Fraction(1, 2)), (60, Fraction(1, 10), 2), (60, 10, 3)]:
    a = Fraction(lam) / Fraction(241 * 10 ** 4)
    xi1 = Fraction(round(5e-7 / float(a) ** 0.5 * 1e12), 10 ** 12)
    b = Fraction(round((lam * 2.41e6) ** 0.5 * 1000), 1000)
    bs = Fraction(round((35 * 3.03e6) ** 0.5 * 1000), 1000)
    pt = dict(kind=f"AlN/sapphire lam={lam} tauR={float(tR_ns)}ns taul={float(tl_ns)}ns",
              xi1=xi1, b=b, bs=bs, tau_R=Fraction(tR_ns) * Fraction(1, 10 ** 9),
              tau_l=Fraction(tl_ns) * Fraction(1, 10 ** 9))
    # 8 log-spaced frequencies 10 kHz .. 200 MHz (omega = 2 pi f, rationalised)
    fs = [10 ** (4 + k * (8.30103 - 4) / 7) for k in range(8)]
    pt["omegas"] = [Fraction(round(2 * 3.141592653589793 * f)) for f in fs]
    points.append(pt)

json.dump([{k: (frac_str(v) if isinstance(v, Fraction) else
                ([frac_str(x) for x in v] if isinstance(v, list) else v))
            for k, v in pt.items()} for pt in points],
          open(HERE / "results" / "01_points.json", "w"), indent=1)

# ------------------------------------------------------------------ Jacobians
jsym = M.sympy_jacobian_mp(with_R=False)
NS = [1, 2, 3, 4, 8]
SUBSETS = {1: [3], 2: [1, 5], 3: [0, 3, 6], 4: [0, 2, 5, 7], 8: list(range(8))}

log("Rank of the front-face Jacobian, theta = (xi1, b, tau_R, tau_l); 50-digit arithmetic.")
log("Singular values of the column-equilibrated Jacobians, normalised by the largest (s_i/s_1).")
log("Numerical rank: number of s_i/s_1 > 1e-30 (an exact zero shows up below 1e-45 at this precision).")
log("")
t0 = time.time()
summary = {}
worst_ab = mp.mpf(0)
for ip, pt in enumerate(points):
    xi1, b, bs, tR, tl = pt["xi1"], pt["b"], pt["bs"], pt["tau_R"], pt["tau_l"]
    rows_a, rows_b = [], []
    for om in pt["omegas"]:
        p = mp.mpc(0, M._mpf(om))
        ra = jsym(p, xi1, b, tR, tl, bs, 0)
        rb = M.jac_fd_mp(p, (xi1, b, tR, tl), bs, model="quad", dps_inner=100)
        rows_a.append(ra)
        rows_b.append(rb)
        for x, y in zip(ra, rb):
            scale = max(abs(z) for z in ra)
            worst_ab = max(worst_ab, abs(x - y) / scale)
    log(f"[{ip:02d}] {pt['kind']}: xi1={float(xi1):.4g} b={float(b):.4g} bs={float(bs):.4g} "
        f"tau_R={float(tR):.4g} tau_l={float(tl):.4g}")
    for N in NS:
        sel = SUBSETS[N]
        Jr = M.equilibrate(M.real_jacobian([rows_a[i] for i in sel]))
        Jc = M.equilibrate(M.complex_jacobian([rows_a[i] for i in sel]))
        sr = M.singular_values(Jr)
        sc = M.singular_values(Jc, complex_=True)
        nr = [x / sr[0] for x in sr]
        nc = [x / sc[0] for x in sc]
        rk_r = sum(1 for x in nr if x > mp.mpf("1e-30"))
        rk_c = sum(1 for x in nc if x > mp.mpf("1e-30"))
        # the same with derivative set (b)
        Jr_b = M.equilibrate(M.real_jacobian([rows_b[i] for i in sel]))
        sr_b = M.singular_values(Jr_b)
        dsv = max(abs(x - y) / sr[0] for x, y in zip(sr, sr_b))
        summary.setdefault((pt["kind"].split()[0], N), []).append((rk_r, rk_c, nr[-1], nc[-1]))
        log(f"     N={N}: real rank {rk_r}  s/s1 = [{', '.join(mp.nstr(x, 3) for x in nr)}]"
            f" | complex rank {rk_c}  s/s1 = [{', '.join(mp.nstr(x, 3) for x in nc)}]"
            f" | (a)-(b) singular values {mp.nstr(dsv, 2)}")
log("")
log(f"max relative difference of Jacobian entries, symbolic (a) vs numerical differentiation (b): {mp.nstr(worst_ab, 3)}")
log(f"elapsed {time.time() - t0:.1f} s")
log("")
log("Summary by family and N: (real ranks found, complex ranks found, min s_min/s_1 real, min s_min/s_1 complex)")
for key in sorted(summary):
    vals = summary[key]
    rr = sorted(set(v[0] for v in vals))
    cc = sorted(set(v[1] for v in vals))
    mr = min(v[2] for v in vals)
    mc = min(v[3] for v in vals)
    log(f"  {key[0]:12s} N={key[1]}: real ranks {rr}, complex ranks {cc}, min s_min/s1 real {mp.nstr(mr, 3)}, complex {mp.nstr(mc, 3)}")

(HERE / "results" / "01_rank_mp.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
