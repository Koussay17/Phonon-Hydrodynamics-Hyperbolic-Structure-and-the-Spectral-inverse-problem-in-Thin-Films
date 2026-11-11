"""Task 2, rigorous rank lower bounds on the strata (Arb ball arithmetic, forward-mode AD).

Upper bounds on the strata are the exact identities (02_sigma_symbolic.py, 02_sigma.wl):
  tau_R = tau_l           : col3 + col4 = 0               -> rank <= 3
  tau_R = tau_l, b = bs   : col1 = 0 and col3 + col4 = 0   -> rank <= 2
Lower bounds: certified nonzero minors at exact rational points.  For each stratum and each
frequency count N the best minor of the expected size is enclosed; for the real Jacobian also
the next size up is reported (its ball must contain 0 if the upper bound is right).
Output: results/02_sigma_ball.txt
"""
import itertools
import random
import sys
from fractions import Fraction as Fr
from pathlib import Path

from flint import acb, arb_mat, acb_mat, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

ctx.prec = 320
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def best_minor(rows, k, complex_=False):
    best = None
    for ri in itertools.combinations(range(len(rows)), k):
        for ci in itertools.combinations(range(len(rows[0])), k):
            sub = [[rows[i][j] for j in ci] for i in ri]
            d = (acb_mat(sub) if complex_ else arb_mat(sub)).det()
            mag = abs(d).mid() if complex_ else abs(d.mid())
            if best is None or mag > best[0]:
                best = (mag, d, ci)
    return best[1], best[2]


def fmt(d):
    if isinstance(d, acb):
        return f"|.|={abs(d).mid().str(4, radius=False)} rad({d.real.rad().str(2, radius=False)},{d.imag.rad().str(2, radius=False)})"
    return f"{d.mid().str(5, radius=False)} +/- {d.rad().str(2, radius=False)}"


rng = random.Random(31)


def rq(lo, hi, den=100):
    return Fr(rng.randint(int(lo * den), int(hi * den)), den)


strata = []
for k in range(4):
    tau = rq(0.01, 3)
    strata.append(("tau_R = tau_l > 0, b != bs", (rq(0.2, 2), rq(0.2, 3), tau, tau), rq(3.1, 5), 3))
for k in range(2):
    strata.append(("Fourier tau_R = tau_l = 0, b != bs", (rq(0.2, 2), rq(0.2, 3), Fr(0), Fr(0)), rq(3.1, 5), 3))
for k in range(3):
    tau, bb = rq(0.01, 3), rq(0.2, 3)
    strata.append(("tau_R = tau_l > 0, b = bs", (rq(0.2, 2), bb, tau, tau), bb, 2))
for k in range(2):
    bb = rq(0.2, 3)
    strata.append(("Fourier tau = 0, b = bs", (rq(0.2, 2), bb, Fr(0), Fr(0)), bb, 2))

all_ok = True
log("Certified lower bounds on the strata (320-bit Arb balls).")
for name, th, bs, r_exp in strata:
    xi1 = th[0]
    oms = [rq(0.05, 30) / (xi1 * xi1) for _ in range(4)]
    rows = [M.jac_ball(om, th, bs) for om in oms]
    log(f"{name}: xi1={th[0]}, b={th[1]}, tau={th[2]}, bs={bs}; expected function-level rank {r_exp}")
    # real Jacobian, N = 1, 2, 4
    for N in (1, 2, 4):
        R = []
        for r in rows[:N]:
            R.append([c.real for c in r])
            R.append([c.imag for c in r])
        k = min(r_exp, 2 * N)
        d, ci = best_minor(R, k)
        ok = not d.contains(0)
        msg = f"   real N={N}: best {k}x{k} minor (cols {ci}) {fmt(d)} {'ok' if ok else 'NOT CERTIFIED'}"
        if k + 1 <= min(2 * N, 4):
            d2, _ = best_minor(R, k + 1)
            msg += f";  largest {k+1}x{k+1} minor {fmt(d2)} (contains 0: {d2.contains(0)})"
            ok &= d2.contains(0)
        log(msg)
        all_ok &= ok
    # complex Jacobian, N = r_exp
    C = [list(r) for r in rows[:r_exp]]
    d, ci = best_minor(C, r_exp, complex_=True)
    ok = not d.contains(0)
    log(f"   complex N={r_exp}: best {r_exp}x{r_exp} minor (cols {ci}) {fmt(d)} {'ok' if ok else 'NOT CERTIFIED'}")
    all_ok &= ok
    # identities in ball arithmetic
    s34 = max(abs(r[2] + r[3]).mid() for r in rows)
    log(f"   max |col3 + col4| over the 4 frequencies (midpoint): {s34.str(3, radius=False)}"
        + (f";  max |col1|: {max(abs(r[0]).mid() for r in rows).str(3, radius=False)}" if r_exp == 2 else ""))
log("")
log(f"all lower-bound certificates valid and all next-size minors contain 0: {all_ok}")
(HERE / "results" / "02_sigma_ball.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
