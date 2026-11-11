"""Task 5 follow-up: scaling of the small singular values near the degenerate strata (40 digits).

(i)   R = 0, T_R = T, T_l = T + eps: two singular values vanish on the diagonal (rank 3); their
      orders in eps.
(ii)  T_R = T_l = T, Rhat = delta -> 0: rank 4 for delta > 0, rank 3 at delta = 0; order in delta.
(iii) off the diagonal, Rhat = delta -> 0: rank 5 throughout (no degeneracy at R = 0).
Equilibrated columns (no identically-zero column occurs at these points), 65 frequencies in [1e-4, 1e4].
Output: results/05_scaling.txt
"""
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

mp.mp.dps = 40
jsym = M.sympy_jacobian_mp(with_R=True)
Ps = [mp.mpf(10) ** (-4 + 8 * k / 64) for k in range(65)]
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def sv(beta, TR, TL, Rh):
    rows = [jsym(mp.mpc(0, P), 1, beta, TR, TL, 1, Rh) for P in Ps]
    s = M.singular_values(M.equilibrate(M.real_jacobian(rows)))
    return [x / s[0] for x in s]


def slope(xs, ys):
    import math
    lx = [math.log10(float(x)) for x in xs]
    ly = [math.log10(float(y)) for y in ys]
    n = len(lx)
    mx, my = sum(lx) / n, sum(ly) / n
    return sum((a - mx) * (b - my) for a, b in zip(lx, ly)) / sum((a - mx) ** 2 for a in lx)


eps = [mp.mpf(10) ** (-k) for k in range(1, 6)]
log("(i) R = 0, T_R = T, T_l = T + eps: s4/s1, s5/s1 and fitted orders in eps")
for beta in (mp.mpf("0.1"), mp.mpf(1), mp.mpf(100)):
    for T in (mp.mpf(0), mp.mpf(1)):
        s4, s5 = [], []
        for e in eps:
            s = sv(beta, T, T + e, 0)
            s4.append(s[3])
            s5.append(s[4])
        log(f"   beta={mp.nstr(beta, 3)} T={mp.nstr(T, 3)}: s4 = [{', '.join(mp.nstr(x, 3) for x in s4)}] (order {slope(eps, s4):.2f});"
            f" s5 = [{', '.join(mp.nstr(x, 3) for x in s5)}] (order {slope(eps, s5):.2f})")
log("(ii) T_R = T_l = T, Rhat = delta: s4/s1 (s5 is the exact zero of the diagonal)")
for beta in (mp.mpf("0.1"), mp.mpf(1), mp.mpf(100)):
    for T in (mp.mpf(0), mp.mpf(1)):
        s4, s5 = [], []
        for d in eps:
            s = sv(beta, T, T, d)
            s4.append(s[3])
            s5.append(s[4])
        log(f"   beta={mp.nstr(beta, 3)} T={mp.nstr(T, 3)}: s4 = [{', '.join(mp.nstr(x, 3) for x in s4)}] (order {slope(eps, s4):.2f});"
            f" max s5 = {mp.nstr(max(s5), 3)}")
log("(iii) off the diagonal (T_R = 1, T_l = 0.3), Rhat = delta: s5/s1 stays bounded away from 0")
for beta in (mp.mpf("0.1"), mp.mpf(1), mp.mpf(100)):
    s5 = [sv(beta, mp.mpf(1), mp.mpf("0.3"), d)[4] for d in eps + [mp.mpf(0)]]
    log(f"   beta={mp.nstr(beta, 3)}: s5 (delta = 1e-1 .. 1e-5, 0) = [{', '.join(mp.nstr(x, 3) for x in s5)}]")
(HERE / "results" / "05_scaling.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
