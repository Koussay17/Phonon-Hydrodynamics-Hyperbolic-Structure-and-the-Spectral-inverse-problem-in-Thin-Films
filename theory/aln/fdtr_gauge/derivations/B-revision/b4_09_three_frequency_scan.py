"""B-4, step 9: do three frequencies give local rank 6 at every point off the degenerate set? (numerical scan)

For random parameter points (interior, b = b_s, Cattaneo tau_l = 0, tau_R = 0, R1 = 0, R2 = 0, near-resonant, resonant
with R2 > 0 for the entropic convention) the best inverse condition number of the 6 x 6 real Jacobian J_3 (rows Re/Im
of d log Z_tot/d log-parameters at three frequencies; times in absolute units when a boundary stratum sets them to 0)
is searched over a grid of frequency triples in [1e-2, 1e2]^3. An exact rank drop at a point would show as a best
inverse condition number at the rounding level (1e-14); the degenerate strata serve as controls.
Double precision, analytic model of b4_06. Output: b4_09_three_frequency_scan.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import itertools
import numpy as np

XA, BA = 0.2, 3.0
grid = np.logspace(-2, 2, 13)
TRIPLES = [t for t in itertools.combinations(grid, 3)]


def ztot(th, w, conv):
    xi, b, tR, tl, R1, R2 = th
    p = 1j * w
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(p)
    N = 1 + (tl * p if conv == "entropic" else 0)
    L = (R2 + 1 / S) / N
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    W = R1 + N * (1 + r * E) / ((1 - r * E) * B)
    At, Bt = XA * np.sqrt(p), BA * np.sqrt(p)
    return (np.cosh(At) * W + np.sinh(At) / Bt) / (Bt * np.sinh(At) * W + np.cosh(At))


def dlog_columns(th, conv):
    """d log Z / d(log xi, log b, tau_R, tau_l, R1, R2) on the whole grid (times and resistances absolute)."""
    th = np.asarray(th, float)
    cols = []
    for j in range(6):
        if j < 2:
            h = 1e-6
            tp, tm = th.copy(), th.copy()
            tp[j] *= np.exp(h)
            tm[j] *= np.exp(-h)
        else:
            h = 1e-7 * max(1.0, th[j])
            tp, tm = th.copy(), th.copy()
            tp[j] += h
            tm[j] -= h
        cols.append((np.log(ztot(tp, grid, conv)) - np.log(ztot(tm, grid, conv))) / (2 * h))
    return np.array(cols).T          # (len(grid), 6)


def best_icond(th, conv):
    C = dlog_columns(th, conv)
    C = C / np.linalg.norm(np.concatenate([C.real, C.imag]), axis=0)
    best = 0.0
    for t in itertools.combinations(range(len(grid)), 3):
        Jr = np.concatenate([C[list(t)].real, C[list(t)].imag])
        s = np.linalg.svd(Jr, compute_uv=False)
        best = max(best, s[-1] / s[0])
    return best


rng = np.random.default_rng(7)


def rnd():
    xi = 1.0
    b = np.exp(rng.uniform(-1.2, 1.2))
    tR, tl = np.exp(rng.uniform(np.log(1e-2), np.log(3.0), 2))
    R1, R2 = np.exp(rng.uniform(np.log(1e-2), np.log(1.0), 2))
    return [xi, b, tR, tl, R1, R2]


strata = {
    "interior": lambda q: q,
    "b=bs": lambda q: [q[0], 1.0] + q[2:],
    "tl=0": lambda q: q[:3] + [0.0] + q[4:],
    "tR=0": lambda q: q[:2] + [0.0] + q[3:],
    "R1=0": lambda q: q[:4] + [0.0, q[5]],
    "R2=0": lambda q: q[:5] + [0.0],
    "near-resonant (tl = 1.01 tR)": lambda q: q[:3] + [1.01 * q[2]] + q[4:],
    "resonant R2>0 (control for local)": lambda q: q[:3] + [q[2]] + q[4:],
}
OUT = {}
for conv in ("local", "entropic"):
    for name, f in strata.items():
        vals = [best_icond(f(rnd()), conv) for _ in range(40)]
        OUT[f"{conv}:{name}"] = {"min": float(np.min(vals)), "median": float(np.median(vals)),
                                 "n_below_1e-10": int(np.sum(np.array(vals) < 1e-10))}
        print(f"{conv}:{name:36s} min={np.min(vals):.2e} median={np.median(vals):.2e}")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_09_three_frequency_scan.json"), "w") as fh:
    json.dump(dict(grid=list(grid), results=OUT), fh, indent=1)
