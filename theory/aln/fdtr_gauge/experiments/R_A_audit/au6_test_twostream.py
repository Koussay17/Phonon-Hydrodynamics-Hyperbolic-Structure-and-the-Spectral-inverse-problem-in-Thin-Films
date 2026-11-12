"""AU6 test: the modal slab solver of au6_slab.py against the closed-form two-stream (single class) black-wall slab.
q(0)/T_A = 2 c v kr (ch + kr sh) / [(1 + kr^2) sh + 2 kr ch],  kr = sqrt(p/g), lam = sqrt(g p)/v, ch/sh at lam d,
static limit G = 2 c v Lambda / (d + 2 Lambda). Also a two-class case against an independent finite-difference
(Chebyshev collocation) solution of the same boundary-value problem. Output: out/au6_test_twostream.json"""
import json, os, sys
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from au6_slab import slab_twoport

res = {}
c, v, tau, d = 1.3e6, 4000.0, 80e-12, 400e-9
cl = {"c": np.array([c]), "tau": np.array([tau]), "v": np.array([v]), "c0": np.zeros(0), "tau0": np.zeros(0)}
rows = []
for f in (1e3, 1e6, 1e8, 1e9, 1e10):
    p = 2j * np.pi * f
    Y, Th, dg = slab_twoport(cl, d, p)
    g = p + 1 / tau
    kr = np.sqrt(p / g)
    lam = np.sqrt(g * p) / v
    ch, sh = np.cosh(lam * d), np.sinh(lam * d)
    q0 = 2 * c * v * kr * (ch + kr * sh) / ((1 + kr ** 2) * sh + 2 * kr * ch)
    rows.append({"f": f, "rel_err_Y11": float(abs(Y[0, 0] / q0 - 1)), "energy_balance": dg["energy_balance_rel"]})
Lam = v * tau
G_static = 2 * c * v * Lam / (d + 2 * Lam)
Y, _, _ = slab_twoport(cl, d, 2j * np.pi * 1.0)
res["single_class"] = {"rows": rows, "static_rel_err": float(abs(Y[0, 0].real / G_static - 1))}
print(res["single_class"])

# two classes + one local class: Chebyshev collocation of the first-order system for T+ and T- (independent method)
cl2 = {"c": np.array([0.8e6, 0.5e6]), "tau": np.array([20e-12, 400e-12]), "v": np.array([6000.0, 2500.0]),
       "c0": np.array([0.3e6]), "tau0": np.array([150e-12])}
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    cc = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(cc, 1 / cc) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x
def colloc(cl, d, p, N=120):
    D, x = cheb(N)
    z = (1 - x) * d / 2                     # z from 0 (x=1) to d (x=-1)
    Dz = -2 / d * D
    nc = cl["c"].size
    n = N + 1
    W = 2 * np.sum(cl["c"] / cl["tau"]) + np.sum(cl["c0"] * p / (1 + p * cl["tau0"]))
    beta = (cl["c"] / cl["tau"]) / W
    A = np.zeros((2 * nc * n, 2 * nc * n), complex)
    rhs = np.zeros(2 * nc * n, complex)
    def idx(i, s):  # class i, s = 0 (+), 1 (-)
        return slice((2 * i + s) * n, (2 * i + s + 1) * n)
    for i in range(nc):
        g = p + 1 / cl["tau"][i]
        vi = cl["v"][i]
        for s, sgn in ((0, 1.0), (1, -1.0)):
            r = idx(i, s)
            A[r, r] += sgn * vi * Dz + g * np.eye(n)
            for j in range(nc):
                for s2 in (0, 1):
                    A[r, idx(j, s2)] += -(1 / cl["tau"][i]) * beta[j] * np.eye(n)
    # boundary rows: T+ at z = 0 (x = 1, index 0) = T_A ; T- at z = d (x = -1, index N) = T_B
    out = []
    for TA, TB in ((1.0, 0.0), (0.0, 1.0)):
        Ab, rb = A.copy(), rhs.copy()
        for i in range(nc):
            rp = idx(i, 0).start + 0
            Ab[rp, :] = 0; Ab[rp, rp] = 1; rb[rp] = TA
            rm = idx(i, 1).start + N
            Ab[rm, :] = 0; Ab[rm, rm] = 1; rb[rm] = TB
        sol = np.linalg.solve(Ab, rb)
        q0 = sum(cl["c"][i] * cl["v"][i] * (sol[idx(i, 0)][0] - sol[idx(i, 1)][0]) for i in range(nc))
        qd = sum(cl["c"][i] * cl["v"][i] * (sol[idx(i, 0)][N] - sol[idx(i, 1)][N]) for i in range(nc))
        out.append((q0, qd))
    return np.array([[out[0][0], out[1][0]], [out[0][1], out[1][1]]])
rows2 = []
for f in (1e4, 1e7, 2e8, 3e9):
    p = 2j * np.pi * f
    Ym, _, dg = slab_twoport(cl2, d, p)
    Yc = colloc(cl2, d, p)
    Yc2 = colloc(cl2, d, p, N=160)
    rows2.append({"f": f, "max_rel_diff_modal_vs_colloc": float(np.max(np.abs(Ym - Yc)) / np.max(np.abs(Ym))),
                  "colloc_N120_vs_N160": float(np.max(np.abs(Yc - Yc2)) / np.max(np.abs(Yc))),
                  "energy_balance": dg["energy_balance_rel"]})
res["two_class_plus_local"] = rows2
print(rows2)
(HERE / "out" / "au6_test_twostream.json").write_text(json.dumps(res, indent=1))
