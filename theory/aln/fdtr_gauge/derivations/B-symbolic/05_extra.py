"""Task 5, additional certificates.

(E1) A point on the hypersurface Q := 3 bs^2 (tau_R - tau_l) + 2 R b^3 xi1 - R^2 b^2 bs^2 = 0, where the
     small-p coefficient map (z_0..z_4) and (z_0..z_3, z_5) of 04_expansions.py is singular: the
     function-level Jacobian still has rank 5 (certified 5 x 5 minor).
(E2) R = 0, tau_R != tau_l, beta = b/bs = 100 (smallest s5/s1 ~ 1e-8 in 05_scaling): rank 5 certified.
(E3) b = bs, tau_R != tau_l, R = 0, where the small-p coefficient map of the 4-parameter model is
     singular (factor b - bs): rank 4 certified (4-parameter model).
Output: results/05_extra.txt
"""
import itertools
import sys
from fractions import Fraction as Fr
from pathlib import Path

from flint import arb_mat, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

ctx.prec = 400
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def best_minor(rows_, k):
    best = None
    for ri in itertools.combinations(range(len(rows_)), k):
        for ci in itertools.combinations(range(len(rows_[0])), k):
            d = arb_mat([[rows_[i][j] for j in ci] for i in ri]).det()
            if best is None or abs(d.mid()) > best[0]:
                best = (abs(d.mid()), d)
    return best[1]


def realrows(cr):
    out = []
    for r in cr:
        out.append([c.real for c in r])
        out.append([c.imag for c in r])
    return out


def certify(name, th, bs, k, oms):
    rows = realrows([M.jac_ball(o, th, bs) for o in oms])
    d = best_minor(rows, k)
    ok = not d.contains(0)
    log(f"{name}: theta={tuple(str(x) for x in th)}, bs={bs}; best {k}x{k} minor of the real Jacobian at "
        f"{len(oms)} frequencies: {d.mid().str(6, radius=False)} +/- {d.rad().str(2, radius=False)}  {'ok' if ok else 'NOT CERTIFIED'}")
    return ok


ok = True
# (E1) xi1 = b = bs = 1, R = 3:  R^2 b^2 bs^2 - 2 R b^3 xi1 = 9 - 6 = 3 = 3 bs^2 (tau_R - tau_l) -> tau_R - tau_l = 1
th = (Fr(1), Fr(1), Fr(3, 2), Fr(1, 2), Fr(3))
Q = 3 * 1 * (th[2] - th[3]) + 2 * th[4] * th[0] - th[4] ** 2
log(f"(E1) Q at the point = {Q}")
ok &= certify("(E1) on Q = 0", th, Fr(1), 5, [Fr(1, 10), Fr(1), Fr(10)])
th = (Fr(7, 10), Fr(2), Fr(2, 5), Fr(1, 10), None)
# second Q = 0 point: bs = 1, choose R from the quadratic R^2 b^2 - 2 R b^3 xi1 - 3 (tau_R - tau_l) = 0 -> R = (2 b^3 xi1 +- sqrt(...))/(2 b^2)
import math
bq, x1, dt = 2.0, 0.7, 0.3
Rpos = (2 * bq ** 3 * x1 + math.sqrt((2 * bq ** 3 * x1) ** 2 + 12 * bq ** 2 * dt)) / (2 * bq ** 2)
Rr = Fr(Rpos).limit_denominator(10 ** 12)
th = (Fr(7, 10), Fr(2), Fr(2, 5), Fr(1, 10), Rr)
Q = 3 * (th[2] - th[3]) + 2 * th[4] * th[1] ** 3 * th[0] - th[4] ** 2 * th[1] ** 2
log(f"(E1) second point: R = {Rr} (rational approximation of the root), Q = {float(Q):.3e}")
ok &= certify("(E1) on Q ~ 0", th, Fr(1), 5, [Fr(1, 10), Fr(1), Fr(10)])
# (E2)
ok &= certify("(E2) R = 0, beta = 100, T_R = 1, T_l = 3/10", (Fr(1), Fr(100), Fr(1), Fr(3, 10), Fr(0)), Fr(1), 5,
              [Fr(1, 100), Fr(1, 3), Fr(3), Fr(30)])
# (E3) 4-parameter model at b = bs off the diagonal
for th, bs in [((Fr(1), Fr(2), Fr(1), Fr(3, 10)), Fr(2)), ((Fr(1, 3), Fr(5), Fr(1, 50), Fr(1, 7)), Fr(5))]:
    ok &= certify("(E3) b = bs, tau_R != tau_l, R = 0 (4-parameter model)", th, bs, 4,
                  [Fr(1, 10) / th[0] ** 2, Fr(3) / th[0] ** 2])
log(f"all certificates valid: {ok}")
(HERE / "results" / "05_extra.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
