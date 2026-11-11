"""B1: rank of the full-band Jacobian, random log-uniform scan.

Nondimensionalisation.  b_s is kept fixed and is used as the effusivity unit (beta = b/b_s, Z b_s);
this is a choice of units, not a parameter reduction (beta stays free).  The time scaling
    (xi1, tR, tl, omega) -> (lam^(1/2) xi1, lam tR, lam tl, omega/lam),  Z -> lam^(1/2) Z
maps the set "all omega > 0" onto itself, so the full-band rank at theta equals the rank at the
scaled point; we put xi1 = 1 at the point examined (t = tau/xi1^2).  The Jacobian is still taken
with respect to all four parameters, xi1 included.

Data: relative response, rows (Re, Im) of dZ/Z = d ln|Z| + i d arg Z, weighted by sqrt(d ln omega),
on a log grid covering every characteristic frequency with 4 decades of padding.  Columns: d/dln xi1,
d/dln b, tau d/dtau (tau > 0) or xi1^2 d/dtau (tau = 0).
Columns are equilibrated to unit norm, so s_4 in [0, 1] is the inverse collinearity index (s_4 = 1:
orthogonal sensitivities; s_4 = 0: exact linear dependence).  Raw (unequilibrated) s_4/s_1 is not
scale-free: at tau = 0 the column d/dtau grows like omega and dominates s_1.
Output: equilibrated singular values s_k/s_1, the smallest one, and its right singular vector.
"""
import sys, os, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


GRID = "wide"      # "wide": log grid, 24 points/decade, +-4 decades (first pass);
                   # "resolved": zm.resolved_grid, +-1 decade, thermal-wave resonances resolved (second pass)


def full_band_svd(beta, tR, tl, variant="code", tscale="log", pad=4.0, ppd=24):
    if GRID == "wide":
        w, wts = zm.wide_grid(1.0, beta, tR, tl, 1.0, pad=pad, ppd=ppd)
        ok = True
    else:
        w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, pad=1.0, nmax=60000, variant=variant)
    M = zm.real_rel_jac_np(w, 1.0, beta, tR, tl, 1.0, 0.0, variant, False, tscale, wts)
    M = M / np.linalg.norm(M, axis=0)          # column equilibration: scale-free collinearity measure
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    if not ok:
        S = S * np.nan
    return S, Vt[-1]


def main(n=20000, seed=1, variant="code"):
    rng = np.random.default_rng(seed)
    rows = []
    t0 = time.time()
    for i in range(n):
        kind = i % 10
        beta = 10 ** rng.uniform(-3, 3)
        tR = 10 ** rng.uniform(-TMAX, TMAX)
        tl = 10 ** rng.uniform(-TMAX, TMAX)
        if kind == 7:
            tl = 0.0                     # Cattaneo stratum
        elif kind == 8:
            tR = 0.0                     # GK without relaxation
        elif kind == 9:
            beta = 1.0                   # matched effusivity
        S, v = full_band_svd(beta, tR, tl, variant)
        rows.append([beta, tR, tl, *(S / S[0]), *v])
    rows = np.array(rows)
    rows = rows[np.isfinite(rows[:, 6])]
    n = len(rows)
    np.save(os.path.join(HERE, "data", f"b1_scan_{variant}{SUFFIX}.npy"), rows)
    s4 = rows[:, 6]
    lr = np.where((rows[:, 1] > 0) & (rows[:, 2] > 0), np.abs(np.log(np.maximum(rows[:, 1], 1e-300) /
                                                                  np.maximum(rows[:, 2], 1e-300))), np.inf)
    summary = {"n": n, "variant": variant, "seconds": time.time() - t0}
    for name, mask in (("all", np.ones(n, bool)),
                       ("interior", (rows[:, 1] > 0) & (rows[:, 2] > 0) & (rows[:, 0] != 1.0)),
                       ("tl=0", rows[:, 2] == 0), ("tR=0", rows[:, 1] == 0),
                       ("beta=1", rows[:, 0] == 1.0)):
        m = mask
        summary[name] = {"count": int(m.sum()), "min_s4": float(s4[m].min()),
                         "median_s4": float(np.median(s4[m])),
                         "argmin": rows[m][np.argmin(s4[m]), :3].tolist()}
    # smallest s4 among interior points at fixed distance from resonance
    for c in (0.01, 0.1, 0.5, 1.0, 3.0):
        m = (lr > c) & (rows[:, 1] > 0) & (rows[:, 2] > 0)
        summary[f"interior_|ln tR/tl|>{c}"] = {"count": int(m.sum()), "min_s4": float(s4[m].min()),
                                               "argmin": rows[m][np.argmin(s4[m]), :3].tolist()}
    print(json.dumps(summary, indent=1))
    with open(os.path.join(HERE, "data", f"b1_scan_{variant}{SUFFIX}.json"), "w") as fh:
        json.dump(summary, fh, indent=1)


TMAX = 5.0
SUFFIX = ""

if __name__ == "__main__":
    variant = sys.argv[1] if len(sys.argv) > 1 else "code"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    if len(sys.argv) > 3 and sys.argv[3] == "resolved":
        GRID, TMAX, SUFFIX = "resolved", 4.0, "_resolved"
    main(n=n, variant=variant)
