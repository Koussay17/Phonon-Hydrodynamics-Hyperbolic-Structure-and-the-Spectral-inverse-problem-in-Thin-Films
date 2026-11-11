"""Task 2, follow-up of the grid search.

(a) Equilibrated columns amplify an identically zero column (round-off) to unit norm; on
    {beta = 1, T_R = T_l} the xi1 column is exactly zero.  Here the natural scaling
    [xi1 dZ/dxi1, b dZ/db, dZ/dT_R, dZ/dT_l] (xi1 = bs = 1) is used instead, at 40 digits.
(b) The smallest off-diagonal s4/s1 of the grid (~1e-7) occurs at T_R = 1e3, T_l = 0, beta ~ 30.
    Is that a band effect or an asymptotic degeneracy as T_R -> infinity?  s4/s1 is computed at
    40 digits for T_R = 1e1 ... 1e5 with (i) the grid band P in [1e-4, 1e4] and (ii) a wide band
    P in [1e-10, 1e6], and the null vector is reported.
Output: results/02_sigma_corner.txt
"""
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

mp.mp.dps = 40
jsym = M.sympy_jacobian_mp(with_R=False)
lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def natural_rows(beta, TR, TL, Ps):
    rows = []
    for P in Ps:
        r = jsym(mp.mpc(0, P), 1, beta, TR, TL, 1, 0)
        rows.append([r[0], beta * r[1], r[2], r[3]])
    return rows


def logspace(a, b, n):
    return [mp.mpf(10) ** (a + (b - a) * k / (n - 1)) for k in range(n)]


log("(a) natural scaling on the diagonal, 65 frequencies in [1e-4, 1e4], 40 digits")
for beta in (mp.mpf(1), mp.mpf("0.3"), mp.mpf(3)):
    for T in (mp.mpf(0), mp.mpf("0.1"), mp.mpf(10)):
        rows = natural_rows(beta, T, T, logspace(-4, 4, 65))
        s = M.singular_values(M.real_jacobian(rows))
        log(f"   beta={mp.nstr(beta, 3)} T={mp.nstr(T, 3)}: s/s1 = [{', '.join(mp.nstr(x / s[0], 3) for x in s)}]")

log("")
log("(b) Cattaneo corner T_l = 0: s4/s1 of the column-equilibrated real Jacobian versus T_R and band;")
log("    weakest direction in the coordinates (ln xi1, ln b, ln T_R, T_l/T_R), normalised to max 1")


def logrows(beta, TR, Ps):
    rows = []
    for P in Ps:
        r = jsym(mp.mpc(0, P), 1, beta, TR, 0, 1, 0)
        rows.append([r[0], beta * r[1], TR * r[2], TR * r[3]])
    return rows


for beta in (mp.mpf("31.6"), mp.mpf(1), mp.mpf("0.1")):
    for TR in (mp.mpf(1), mp.mpf(10), mp.mpf(100), mp.mpf(1000), mp.mpf(10000), mp.mpf(100000)):
        out, vec = [], None
        for (a, b, n) in ((-4, 4, 65), (-10, 6, 161)):
            Jr = M.real_jacobian(logrows(beta, TR, logspace(a, b, n)))
            s = M.singular_values(M.equilibrate(Jr))
            out.append(s[3] / s[0])
            if n == 161:
                U, S, V = mp.svd_r(Jr)
                v = [V[3, j] for j in range(4)]
                vmax = max(abs(x) for x in v)
                vec = [x / vmax for x in v]
        log(f"   beta={mp.nstr(beta, 3)} T_R={mp.nstr(TR, 3)}: s4/s1 band[1e-4,1e4] = {mp.nstr(out[0], 3)},"
            f" band[1e-10,1e6] = {mp.nstr(out[1], 3)};  weakest direction ~ [{', '.join(mp.nstr(x, 3) for x in vec)}]")
(HERE / "results" / "02_sigma_corner.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
