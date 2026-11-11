"""B2: global search for distinct parameter sets with (nearly) equal front-face response.

Misfit  F(theta'; theta, band) = [ sum_k w_k |Z(theta'; omega_k)/Z(theta; omega_k) - 1|^2 / sum_k w_k ]^(1/2)
(RMS relative complex misfit, weights d ln omega).  theta = target (xi1 = bs = 1), theta' free.
Variables of theta': ln xi1', ln beta', ln tR', ln tl'  (strata: one time fixed at 0; Fourier: tR' = tl' = 0).
Multistart trust-region least squares with analytic Jacobian from log-uniform starts in a wide box
(xi1' in [0.1, 10], beta' in [0.01, 100], times in [1e-4, 1e4]); the solver is confined to the box
enlarged by e^2 on each side, and minima on that outer face are flagged (asymptotic directions).
Each local minimum is recorded with its misfit and its distance from theta,
    dist = max_k |ln theta'_k - ln theta_k|  over (xi1, beta) and over the times when both are > 0.
A counterexample to B2 would be a minimum with dist >> 0 and misfit at round-off (then refined at 40 digits).
Bands: "full" (target characteristic frequencies +- 1 decade, thermal-wave resonances resolved) and
three FDTR-like windows of 4.3 decades (the ratio 200 MHz / 10 kHz) placed relative to 1/xi1^2.
"""
import sys, os, json, time
import numpy as np
import mpmath as mp
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm

BOX = {"lxi": (np.log(0.1), np.log(10.0)), "lb": (np.log(1e-2), np.log(1e2)),
       "lt": (np.log(1e-4), np.log(1e4))}


def unpack(x, stratum):
    if stratum == "gk":
        lxi, lb, ltR, ltl = x
        return np.exp(lxi), np.exp(lb), np.exp(ltR), np.exp(ltl)
    if stratum == "cat":            # tl' = 0
        lxi, lb, ltR = x
        return np.exp(lxi), np.exp(lb), np.exp(ltR), 0.0
    if stratum == "gk0":            # tR' = 0
        lxi, lb, ltl = x
        return np.exp(lxi), np.exp(lb), 0.0, np.exp(ltl)
    if stratum == "fourier":
        lxi, lb = x
        return np.exp(lxi), np.exp(lb), 0.0, 0.0
    raise ValueError(stratum)


def residual_and_jac(x, stratum, w, sw, Zt):
    xi1, b, tR, tl = unpack(x, stratum)
    Z, J = zm.jac_np(w, xi1, b, tR, tl, 1.0, 0.0)
    r = (Z / Zt - 1.0) * sw
    cols = [J[:, 0], J[:, 1]]
    if stratum in ("gk", "cat"):
        cols.append(J[:, 2] * tR)
    if stratum in ("gk", "gk0"):
        cols.append(J[:, 3] * tl)
    Jr = np.stack(cols, axis=1) / Zt[:, None] * sw[:, None]
    res = np.concatenate([r.real, r.imag])
    jac = np.concatenate([Jr.real, Jr.imag], axis=0)
    return res, jac


def misfit(theta_p, w, wts, Zt):
    Z = zm.Z_np(w, *theta_p)
    return float(np.sqrt(np.sum(wts * np.abs(Z / Zt - 1) ** 2) / np.sum(wts)))


def distance(tp, t):
    d = [abs(np.log(tp[0] / t[0])), abs(np.log(tp[1] / t[1]))]
    for k in (2, 3):
        if tp[k] > 0 and t[k] > 0:
            d.append(abs(np.log(tp[k] / t[k])))
        elif (tp[k] > 0) != (t[k] > 0):
            d.append(np.inf if max(tp[k], t[k]) > 1e-3 else 0.0)
    return float(max(d))


def search(theta, w, wts, stratum, nstart, rng):
    Zt = zm.Z_np(w, *theta)
    sw = np.sqrt(wts / np.sum(wts))
    dim = {"gk": 4, "cat": 3, "gk0": 3, "fourier": 2}[stratum]
    lo = np.array([BOX["lxi"][0], BOX["lb"][0]] + [BOX["lt"][0]] * (dim - 2))
    hi = np.array([BOX["lxi"][1], BOX["lb"][1]] + [BOX["lt"][1]] * (dim - 2))
    sols = []
    for k in range(nstart):
        x0 = rng.uniform(lo, hi)
        f = lambda x: residual_and_jac(x, stratum, w, sw, Zt)[0]
        g = lambda x: residual_and_jac(x, stratum, w, sw, Zt)[1]
        try:
            r = least_squares(f, x0, jac=g, method="trf", x_scale="jac", xtol=1e-15, ftol=1e-15,
                              gtol=1e-15, max_nfev=3000, bounds=(lo - 2.0, hi + 2.0))
        except Exception as e:      # overflow at extreme points
            continue
        tp = unpack(r.x, stratum)
        on_bound = bool(np.any(r.x <= lo - 2.0 + 1e-6) or np.any(r.x >= hi + 2.0 - 1e-6))
        sols.append({"theta_p": [float(v) for v in tp], "misfit": float(np.sqrt(2 * r.cost)),
                     "dist": distance(tp, theta), "on_bound": on_bound})
    return sols


def cluster(sols, tol=1e-3):
    out = []
    for s in sorted(sols, key=lambda z: z["misfit"]):
        lp = np.log(np.maximum(np.array(s["theta_p"]), 1e-300))
        if any(np.max(np.abs(lp - np.log(np.maximum(np.array(o["theta_p"]), 1e-300)))) < tol for o in out):
            continue
        out.append(s)
    return out


def bands(theta):
    w, wts, ok = zm.resolved_grid(*theta, 1.0, pad=1.0, nmax=60000)
    B = {"full": (w, wts)}
    for name, top in (("fdtr_thin", 1e-2), ("fdtr_mid", 2.0), ("fdtr_thick", 1e2)):
        lo = top / 2e4
        w2, wts2, ok2 = zm.resolved_grid(*theta, 1.0, wmin=lo, wmax=top, nmax=60000)
        B[name] = (w2, wts2)
    return B


TARGETS = {
    "T1_gk": (1.0, 2.0, 0.5, 0.1),
    "T2_gk_tl>tR": (1.0, 0.5, 0.1, 0.5),
    "T3_cattaneo": (1.0, 2.0, 0.5, 0.0),
    "T4_beta1": (1.0, 1.0, 0.5, 0.1),
    "T5_near_res": (1.0, 2.0, 0.5 * np.exp(0.05), 0.5 * np.exp(-0.05)),
    "T6_ballistic": (1.0, 10.0, 100.0, 10.0),
    "T7_lowbeta": (1.0, 0.1, 10.0, 1.0),
    "T8_gk0": (1.0, 2.0, 0.0, 0.3),
}

if __name__ == "__main__":
    t0 = time.time()
    nstart = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rng = np.random.default_rng(5)
    results = {}
    for tname, theta in TARGETS.items():
        results[tname] = {}
        for bname, (w, wts) in bands(theta).items():
            Zt = zm.Z_np(w, *theta)
            rec = {}
            for stratum in ("gk", "cat", "gk0", "fourier"):
                sols = cluster(search(theta, w, wts, stratum, nstart if stratum == "gk" else nstart // 2, rng))
                far = [s for s in sols if s["dist"] > 0.05]
                rec[stratum] = {"n_minima": len(sols), "best_far": far[:4],
                                "best_any": sols[:2]}
            # explicit swap tR <-> tl, refit (xi1, beta) then all four
            sw = np.sqrt(wts / np.sum(wts))
            if theta[2] > 0 and theta[3] > 0:
                x0 = np.log([theta[0], theta[1], theta[3], theta[2]])
                f = lambda x: residual_and_jac(np.concatenate([x, x0[2:]]), "gk", w, sw, Zt)[0]
                g = lambda x: residual_and_jac(np.concatenate([x, x0[2:]]), "gk", w, sw, Zt)[1][:, :2]
                r2 = least_squares(f, x0[:2], jac=g, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15)
                tp2 = unpack(np.concatenate([r2.x, x0[2:]]), "gk")
                f4 = lambda x: residual_and_jac(x, "gk", w, sw, Zt)[0]
                g4 = lambda x: residual_and_jac(x, "gk", w, sw, Zt)[1]
                r4 = least_squares(f4, np.concatenate([r2.x, x0[2:]]), jac=g4, method="trf",
                                   xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
                tp4 = unpack(r4.x, "gk")
                rec["swap"] = {"refit_xi_beta": {"theta_p": list(map(float, tp2)),
                                                 "misfit": float(np.sqrt(2 * r2.cost))},
                               "refit_all_from_swap": {"theta_p": list(map(float, tp4)),
                                                       "misfit": float(np.sqrt(2 * r4.cost)),
                                                       "dist": distance(tp4, theta)}}
            results[tname][bname] = rec
            bf = rec["gk"]["best_far"][0] if rec["gk"]["best_far"] else None
            print(f"{tname:14s} {bname:10s} n_omega={len(w):5d}  best far GK minimum:",
                  None if bf is None else "misfit=%.3e dist=%.3f theta'=%s" % (
                      bf["misfit"], bf["dist"], np.round(bf["theta_p"], 4)),
                  "| Fourier best: %.3e" % rec["fourier"]["best_any"][0]["misfit"],
                  "| swap refit: %s" % ("%.3e" % rec["swap"]["refit_all_from_swap"]["misfit"]
                                        if "swap" in rec else "-"))
            sys.stdout.flush()
    results["seconds"] = time.time() - t0
    with open(os.path.join(HERE, "data", "b2_injectivity.json"), "w") as fh:
        json.dump(results, fh, indent=1)
    print("seconds", results["seconds"])
