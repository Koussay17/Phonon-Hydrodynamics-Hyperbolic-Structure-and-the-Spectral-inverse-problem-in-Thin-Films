"""Task 1, rigorous part. Certified rank lower bounds by Arb ball arithmetic (python-flint).

At the exact rational points of results/01_points.json the Jacobian entries are
enclosed rigorously by forward-mode differentiation in acb_series arithmetic
(no symbolic derivative is used).  A minor whose enclosing ball excludes zero
certifies the corresponding rank lower bound; the upper bounds min(rows, 4) are
trivial.  Also checked: det[Re r1; Im r1; Re r2; Im r2] = -(1/4) det[r(i w1); r(-i w1); r(i w2); r(-i w2)].
Output: results/01_rank_ball.txt
"""
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import flint
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


def F(s):
    n, d = s.split("/")
    return Fraction(int(n), int(d))


def best_minor(rows, k, complex_=False):
    """Largest-|mid| k x k minor over all row/column choices; returns (ball, rows, cols)."""
    best = None
    nr, nc = len(rows), len(rows[0])
    for ri in itertools.combinations(range(nr), k):
        for ci in itertools.combinations(range(nc), k):
            sub = [[rows[i][j] for j in ci] for i in ri]
            d = (acb_mat(sub) if complex_ else arb_mat(sub)).det()
            mag = abs(d).mid() if complex_ else abs(d.mid())
            if best is None or mag > best[0]:
                best = (mag, d, ri, ci)
    return best[1], best[2], best[3]


def certified(d):
    return not d.contains(0)


def fmt(d):
    """midpoint (6 digits) and the true radius of the enclosing ball"""
    if isinstance(d, acb):
        return (f"({d.real.mid().str(6, radius=False)} + {d.imag.mid().str(6, radius=False)}j) +/- "
                f"({d.real.rad().str(2, radius=False)}, {d.imag.rad().str(2, radius=False)})")
    return f"{d.mid().str(6, radius=False)} +/- {d.rad().str(2, radius=False)}"


pts = json.load(open(HERE / "results" / "01_points.json"))
SUB = {1: [3], 2: [1, 5], 3: [0, 3, 6], 4: [0, 2, 5, 7]}
all_ok = True
log("Certified minors (Arb, 320-bit balls).  'ok' = enclosing ball excludes 0.")
for ip, pt in enumerate(pts):
    th = tuple(F(pt[k]) for k in ("xi1", "b", "tau_R", "tau_l"))
    bs = F(pt["bs"])
    oms = [F(x) for x in pt["omegas"]]
    rows = {}
    for idx in sorted(set(sum(SUB.values(), []))):
        rows[idx] = M.jac_ball(oms[idx], th, bs)
    log(f"[{ip:02d}] {pt['kind']}")
    # real N = 1: 2 x 4, rank 2
    r1 = rows[SUB[1][0]]
    R1 = [[c.real for c in r1], [c.imag for c in r1]]
    d, ri, ci = best_minor(R1, 2)
    ok1 = certified(d)
    # real N = 2: 4 x 4 determinant
    R2 = []
    for idx in SUB[2]:
        R2.append([c.real for c in rows[idx]])
        R2.append([c.imag for c in rows[idx]])
    d2 = arb_mat(R2).det()
    ok2 = certified(d2)
    # identity with the complex determinant at +-i omega (conjugate rows)
    C4 = []
    for idx in SUB[2]:
        C4.append(list(rows[idx]))
        C4.append([c.conjugate() for c in rows[idx]])
    dc = acb_mat(C4).det()
    rel = (acb(d2) + dc / 4)
    okid = rel.contains(0)
    # complex N = 1..4
    okc = []
    for N in (1, 2, 3, 4):
        Cr = [list(rows[i]) for i in SUB[N]]
        dN, _, _ = best_minor(Cr, N, complex_=True)
        okc.append((N, certified(dN), dN))
    log(f"     real N=1 2x2 minor {fmt(d)} {'ok' if ok1 else 'NOT CERTIFIED'}")
    log(f"     real N=2 det       {fmt(d2)} {'ok' if ok2 else 'NOT CERTIFIED'};"
        f"  det_real + det_complex(+-i w)/4 contains 0: {okid}")
    for N, okN, dN in okc:
        log(f"     complex N={N} {N}x{N} minor {fmt(dN)} {'ok' if okN else 'NOT CERTIFIED'}")
    all_ok &= ok1 and ok2 and okid and all(o for _, o, _ in okc)
log("")
log(f"all certificates valid: {all_ok}")
log("Consequence: at every listed point the real Jacobian has rank 2 (N=1) and 4 (N=2, hence N>=2),")
log("the complex Jacobian has rank N for N<=4 and 4 for N>=4.  By real-analyticity in (theta, omega)")
log("these ranks are the generic ones, and the exceptional set is a proper analytic subset.")
(HERE / "results" / "01_rank_ball.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
