"""cx_08b: Proposition LF (low-frequency coefficients of the film block Phi = N0 Mob(M; L/N)) checked by trapezoidal
Cauchy integrals in s at 50 digits, for random (xi, b, tR, tl, R) and independent face times (tau_o, tau_n), b_s = 1.3
(b_s != 1 to test the b_s dependence)."""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
rng = np.random.default_rng(8)
bs = mp.mpf("1.3")


def Phi(xi, b, tR, tl, to, tn, R, s):
    p = s * s
    g = (1 + tR * p) / (1 + tl * p); sg = mp.sqrt(g)
    A = xi * s * sg; B = b * s / sg; N = 1 + tn * p; N0 = 1 + to * p; L = R + 1 / (bs * s)
    u = B * L / N; r = (u - 1) / (u + 1); E = mp.exp(-2 * A)
    return N0 * (1 + r * E) / ((1 - r * E) * B)


rows = []
worst = 0
for i in range(8):
    xi, b, tR, tl, R = [mp.mpf(float(x)) for x in (rng.uniform(.4, 1.2), rng.uniform(.5, 2.5), rng.uniform(.05, .8),
                                                   rng.uniform(.05, .8), rng.uniform(0, .5))]
    to, tn = [mp.mpf(float(x)) for x in rng.uniform(-0.3, 0.8, 2)]
    if i % 2 == 0:
        R = mp.mpf(0)
    rho = mp.mpf("0.08")
    M = 200
    c = [mp.mpc(0)] * 4
    for j in range(M):
        s = rho * mp.expjpi(mp.mpf(2 * j) / M)
        v = Phi(xi, b, tR, tl, to, tn, R, s)
        for k in range(-1, 3):
            c[k + 1] += v * s ** (-k) / M
    cm1 = 1 / bs
    c0 = R + xi * (bs ** 2 - b ** 2) / (b * bs ** 2)
    c1 = xi ** 2 * (b ** 2 - bs ** 2) / bs ** 3 - 2 * b * R * xi / bs - (tn - to) / bs
    c2 = (xi / b) * ((tR - tl) + 2 * b ** 2 / bs ** 2 * tn - (b ** 2 - bs ** 2) / bs ** 2 * to) \
        - xi ** 3 * (b ** 2 - bs ** 2) * (3 * b ** 2 - bs ** 2) / (3 * b * bs ** 4)
    e = [abs(c[0] - cm1), abs(c[1] - c0), abs(c[2] - c1)] + ([abs(c[3] - c2)] if R == 0 else [])
    rel = [float(x / (abs(y) + 1e-30)) for x, y in zip(e, [cm1, c0, c1, c2])]
    worst = max(worst, max(rel))
    rows.append(dict(point=[float(x) for x in (xi, b, tR, tl, to, tn, R)], rel_err=rel))
    print("R=%.3f" % float(R), "rel. errors c_-1..c_2:", ["%.1e" % x for x in rel], flush=True)
print("worst relative error", worst)
json.dump(dict(bs=float(bs), rows=rows, worst=worst),
          open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_08b_laurent.json", "w"), indent=1)
