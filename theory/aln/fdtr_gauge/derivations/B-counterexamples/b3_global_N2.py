"""B3 (global aspect): with N = 2 frequencies the data map theta -> (Z(w1), Z(w2)) is square (4 real -> 4 real).
Generic local injectivity (rank 4) does not imply global injectivity.  Search for distinct exact solutions
theta' != theta of Z(theta'; wk) = Z(theta; wk), k = 1, 2, by multistart Newton/least squares, then refine
in 40-digit arithmetic (mp.findroot) to separate exact roots from near-misses.  Same for N = 3 (overdetermined).
"""
import sys, os, json, time
import numpy as np
import mpmath as mp
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


def res_jac(x, ws, Zt):
    xi1, b, tR, tl = np.exp(x)
    Z, J = zm.jac_np(ws, xi1, b, tR, tl)
    r = Z / Zt - 1.0
    Jr = J[:, :4] * np.array([1, 1, tR, tl]) / Zt[:, None]
    return np.concatenate([r.real, r.imag]), np.concatenate([Jr.real, Jr.imag], axis=0)


def solve_all(theta, ws, nstart, rng):
    ws = np.asarray(ws, float)
    Zt = zm.Z_np(ws, *theta)
    lo = np.log([0.05, 0.01, 1e-4, 1e-4]); hi = np.log([20, 100, 1e4, 1e4])
    sols = []
    for k in range(nstart):
        x0 = rng.uniform(lo, hi)
        r = least_squares(lambda x: res_jac(x, ws, Zt)[0], x0, jac=lambda x: res_jac(x, ws, Zt)[1],
                          method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000,
                          bounds=(lo - 3, hi + 3))
        if np.sqrt(2 * r.cost) < 1e-11:
            th = np.exp(r.x)
            if not any(np.max(np.abs(np.log(th / s))) < 1e-6 for s in sols):
                sols.append(th)
    return sols


def refine_mp(theta, ws, x0, dps=40):
    mp.mp.dps = dps
    Zt = [zm.Z_mp(mp.mpf(w), *[mp.mpf(v) for v in theta]) for w in ws]
    def F(a, b, c, d):
        th = [mp.exp(a), mp.exp(b), mp.exp(c), mp.exp(d)]
        out = []
        for w, zt in zip(ws, Zt):
            z = zm.Z_mp(mp.mpf(w), *th) / zt - 1
            out += [z.real, z.imag]
        return out
    x = mp.findroot(F, [mp.log(mp.mpf(v)) for v in x0], tol=mp.mpf(10) ** (-dps + 5))
    th = [mp.exp(v) for v in x]
    resid = max(abs(v) for v in F(*x))
    return [float(v) for v in th], float(resid)


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(23)
    out = {"N2": [], "N3": []}
    cases = [((1.0, 2.0, 0.5, 0.1), (0.3, 3.0)), ((1.0, 2.0, 0.5, 0.1), (0.05, 1.0)),
             ((1.0, 0.5, 0.1, 0.5), (0.3, 3.0)), ((1.0, 1.0, 0.5, 0.1), (0.2, 2.0)),
             ((1.0, 10.0, 3.0, 1.0), (0.1, 1.0))]
    for theta, ws in cases:
        sols = solve_all(theta, ws, 300, rng)
        rec = {"theta": theta, "ws": ws, "n_solutions": len(sols), "solutions": []}
        for s in sols:
            try:
                thr, resid = refine_mp(theta, ws, s)
            except Exception as e:
                thr, resid = list(map(float, s)), None
            rec["solutions"].append({"theta_p": thr, "mp_residual": resid,
                                     "dist": float(np.max(np.abs(np.log(np.array(thr) / np.array(theta)))))})
        out["N2"].append(rec)
        print("N=2 theta=%s ws=%s: %d distinct exact solutions" % (theta, ws, len(sols)))
        for s in rec["solutions"]:
            print("    theta'=%s  dist=%.3f  40-digit residual=%s" % (np.array2string(np.array(s["theta_p"]), precision=6),
                                                                 s["dist"], s["mp_residual"]))
    for theta, ws in [((1.0, 2.0, 0.5, 0.1), (0.05, 0.5, 5.0)), ((1.0, 0.5, 0.1, 0.5), (0.1, 1.0, 10.0))]:
        sols = solve_all(theta, ws, 300, rng)
        out["N3"].append({"theta": theta, "ws": ws, "n_solutions": len(sols),
                          "solutions": [list(map(float, s)) for s in sols]})
        print("N=3 theta=%s ws=%s: %d distinct exact solutions" % (theta, ws, len(sols)))
        for s in sols:
            print("    theta'=%s" % np.array2string(s, precision=6))
    out["seconds"] = time.time() - t0
    with open(os.path.join(HERE, "data", "b3_global_N2.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("seconds", out["seconds"])
