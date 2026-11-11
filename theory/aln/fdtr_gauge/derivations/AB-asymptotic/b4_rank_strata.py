"""B: numerical rank of the function-level Jacobian on the candidate strata (independent check of the
asymptotic rank argument).

Units: xi1 = 1, b_s = 1, times in units of xi1^2. Columns d ln Z / d(ln xi1, ln b, tau_R, tau_l), rows
Re and Im of d ln Z (i.e. d ln|Z| and d arg Z) on a log grid Omega = omega xi1^2 in [1e-6, 1e8].
Columns are normalised to unit norm before the SVD (separates weak sensitivity from dependence).
Predictions (asymptotic argument, derivations/AB-asymptotic.md section B4):
  rank 4 for every theta with tau_R != tau_l (including b = b_s, tau_l = 0, tau_R = 0);
  rank 3 on {tau_R = tau_l, b != b_s}; rank 2 on {tau_R = tau_l, b = b_s};
  sigma_min ~ |delta| as delta = tau_R - tau_l -> 0 (exact: (d_tR + d_tl) Z = -dZ/dG p^2 delta/(1+tau_l p)^2).
"""
import json
import os

import numpy as np
import mpmath as mp

from ab_common import film_jacobian, film_Z

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results", "b4_rank_strata.json")
Om = np.logspace(-6, 8, 561)
p = 1j * Om


def svals(b, tr, tl, grid=p):
    Z, J = film_jacobian(grid, 1.0, b, tr, tl, 1.0, tscale=1.0)
    R = np.vstack([J.real, J.imag])
    nrm = np.linalg.norm(R, axis=0)
    s = np.linalg.svd(R / nrm, compute_uv=False)
    return s, nrm


# finite-difference check of the analytic Jacobian
def fd_check():
    th = (1.0, 2.13, 0.01, 0.004)
    Z, J = film_jacobian(p, *th, 1.0, tscale=1.0)
    h = 1e-6
    out = []
    for i in range(4):
        tp, tm = list(th), list(th)
        if i < 2:
            tp[i] *= np.exp(h); tm[i] *= np.exp(-h)
        else:
            tp[i] += h * 0.01; tm[i] -= h * 0.01
        d = (np.log(film_Z(p, *tp, 1.0)) - np.log(film_Z(p, *tm, 1.0))) / (2 * h * (1 if i < 2 else 0.01))
        out.append(float(np.max(np.abs(d - J[:, i])) / np.max(np.abs(J[:, i]))))
    return out


res = {"fd_check_rel_err": fd_check()}
print("analytic vs finite-difference Jacobian, rel. err per column:", ["%.1e" % v for v in res["fd_check_rel_err"]])

cases = {
    "generic GK (b=2.13, tR=0.01, tl=0.004)": (2.13, 0.01, 0.004),
    "matched b=b_s, GK (tR=0.01, tl=0.004)": (1.0, 0.01, 0.004),
    "Cattaneo (b=2.13, tR=0.01, tl=0)": (2.13, 0.01, 0.0),
    "Cattaneo at b=b_s": (1.0, 0.01, 0.0),
    "tau_R=0, tl=0.01 (nonlocal only)": (2.13, 0.0, 0.01),
    "tau_l > tau_R (tR=0.004, tl=0.01)": (2.13, 0.004, 0.01),
    "resonance tR=tl=0.01": (2.13, 0.01, 0.01),
    "resonance and b=b_s": (1.0, 0.01, 0.01),
    "Fourier tR=tl=0": (2.13, 0.0, 0.0),
    "large times (tR=3, tl=0.5), b=0.6": (0.6, 3.0, 0.5),
}
res["cases"] = {}
for name, (b, tr, tl) in cases.items():
    s, nrm = svals(b, tr, tl)
    res["cases"][name] = {"singular_values": s.tolist(), "column_norms": nrm.tolist()}
    print("%-45s s = %s" % (name, " ".join("%.2e" % v for v in s)))

# delta -> 0 scans
scan = {}
for label, b in (("b=2.13", 2.13), ("b=b_s", 1.0)):
    rows = []
    for d in np.logspace(-1, -7, 13):
        s, _ = svals(b, 0.01 + d * 0.01, 0.01)
        rows.append([d] + s.tolist())
    rows = np.array(rows)
    sl3 = np.polyfit(np.log(rows[-6:, 0]), np.log(rows[-6:, 3]), 1)[0]
    sl4 = np.polyfit(np.log(rows[-6:, 0]), np.log(rows[-6:, 4]), 1)[0]
    scan[label] = {"rows_reldelta_s1..s4": rows.tolist(), "slope_s3": sl3, "slope_s4": sl4}
    print("delta scan %s: slope(s3) = %.3f, slope(s4) = %.3f" % (label, sl3, sl4))
res["delta_scan"] = scan

# b -> b_s scan at fixed delta != 0
rows = []
for db in [0.5, 0.1, 1e-2, 1e-3, 1e-4, 0.0, -1e-3, -0.1]:
    s, _ = svals(1.0 + db, 0.01, 0.004)
    rows.append([db] + s.tolist())
res["b_scan_delta_fixed"] = rows
print("b -> b_s scan (delta = 0.006): s4 =", ["%.3e" % r[4] for r in rows])

# along b = b_s: sigma_min vs delta (predicted ~ delta: the xi1 sensitivity there is O(delta omega))
rows = []
for d in np.logspace(-1, -5, 9):
    s, _ = svals(1.0, 0.004 + d, 0.004)
    rows.append([d] + s.tolist())
rows = np.array(rows)
res["b_eq_bs_delta_scan"] = rows.tolist()
print("b=b_s, delta scan: slopes s3, s4 =",
      "%.3f %.3f" % (np.polyfit(np.log(rows[-5:, 0]), np.log(rows[-5:, 3]), 1)[0],
                     np.polyfit(np.log(rows[-5:, 0]), np.log(rows[-5:, 4]), 1)[0]))

# 60-digit cross-check of the smallest singular value at two delicate points (Jacobian by mp.diff)
mp.mp.dps = 40


def Zmp(w, X, Bb, TR, TL):
    pp = mp.mpc(0, w)
    gm = mp.sqrt((1 + TR * pp) / (1 + TL * pp))
    sq = mp.sqrt(pp)
    A = X * sq * gm
    B = Bb * sq / gm
    t = mp.tanh(A)
    return (B + sq * t) / (B * (sq + B * t))


def mp_smin(b, tr, tl, npts=141):
    grid = [mp.mpf(10) ** e for e in np.linspace(-6, 8, npts)]
    rows = []
    for w in grid:
        f = lambda X, Bb, TR, TL: mp.log(Zmp(w, X, Bb, TR, TL))
        th = [mp.mpf(1), mp.mpf(b), mp.mpf(tr), mp.mpf(tl)]
        col = []
        for i in range(4):
            def fi(v, i=i):
                a = list(th)
                a[i] = v
                return f(*a)
            dv = mp.diff(fi, th[i])
            col.append(dv * (th[i] if i < 2 else 1))
        rows.append(col)
    M = np.array([[complex(c) for c in r] for r in rows])
    R = np.vstack([M.real, M.imag])
    R = R / np.linalg.norm(R, axis=0)
    return np.linalg.svd(R, compute_uv=False)


for name, (b, tr, tl) in (("matched b=b_s, GK", (1.0, 0.01, 0.004)), ("resonance", (2.13, 0.01, 0.01))):
    s_mp = mp_smin(b, tr, tl)
    s_np, _ = svals(b, tr, tl, grid=1j * np.logspace(-6, 8, 141))
    res.setdefault("mp_crosscheck", {})[name] = {"mp": s_mp.tolist(), "np": s_np.tolist()}
    print("mp cross-check %-20s mp: %s | np: %s" % (name, " ".join("%.3e" % v for v in s_mp),
                                                     " ".join("%.3e" % v for v in s_np)))

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
print("written", OUT)

# ------------------------------------------------------------------------------------------------
# Column norms (normalisation hides weak columns): un-normalised singular values with a fixed
# reference scaling (column norms of the generic case), and ||d ln Z/d ln xi1|| near b = b_s.
# ------------------------------------------------------------------------------------------------
_, nref = svals(2.13, 0.01, 0.004)


def svals_fixed(b, tr, tl):
    Z, J = film_jacobian(p, 1.0, b, tr, tl, 1.0, tscale=1.0)
    R = np.vstack([J.real, J.imag]) / nref
    return np.linalg.svd(R, compute_uv=False), np.linalg.norm(R, axis=0)


extra = {}
for name, (b, tr, tl) in (("resonance and b=b_s", (1.0, 0.01, 0.01)), ("resonance", (2.13, 0.01, 0.01)),
                          ("matched b=b_s, GK", (1.0, 0.01, 0.004)), ("generic GK", (2.13, 0.01, 0.004))):
    s, cn = svals_fixed(b, tr, tl)
    extra[name] = {"s_fixed_scaling": s.tolist(), "column_norms_rel_generic": cn.tolist()}
    print("%-22s fixed-scale s = %s | col norms = %s" % (name, " ".join("%.2e" % v for v in s),
                                                          " ".join("%.2e" % v for v in cn)))
rows = []
for d in np.logspace(-1, -5, 9):
    s, cn = svals_fixed(1.0, 0.004 + d, 0.004)
    rows.append([d] + s.tolist() + cn.tolist())
rows = np.array(rows)
sl = [np.polyfit(np.log(rows[-5:, 0]), np.log(rows[-5:, j]), 1)[0] for j in (3, 4, 5)]
print("b=b_s delta scan, fixed scaling: slopes s3 %.3f s4 %.3f ||d/dln xi1|| %.3f" % tuple(sl))
extra["b_eq_bs_delta_scan_fixed"] = {"rows_delta_s1..s4_colnorms1..4": rows.tolist(), "slopes_s3_s4_xi1col": sl}
rows = []
for db in [0.5, 0.1, 1e-2, 1e-3, 1e-4, 0.0]:
    s, cn = svals_fixed(1.0 + db, 0.01, 0.004)
    rows.append([db] + s.tolist() + cn.tolist())
extra["b_scan_fixed"] = rows
print("b -> b_s (delta=0.006), ||d/dln xi1|| rel generic:", ["%.3e" % r[5] for r in rows])
res["fixed_scaling"] = extra
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=1)
