"""B1: adversarial minimisation of the smallest singular value of the full-band Jacobian.

Search variables (xi1 = bs = 1): (log10 beta, log10 tl, d) with tR = tl e^d, |d| >= c (distance from
the Fourier resonance), log10 beta in [-2, 2], log10 tl, log10 tR in [-3, 3] (penalised, never clipped:
clipping both times to the same face lands exactly on the resonance and is a spurious minimum).
Objective: log10 s4 on the resolved full-band grid (characteristic frequencies +- 1 decade, thermal-
wave resonances resolved), metric "eq" (column-equilibrated) or "raw" (log parameters, weights
d ln omega).  Local minima are classified
  interior   (no constraint active)   -> candidate hidden degeneracy, re-evaluated at 50 digits
  box        (a box face active)      -> asymptotic (practical) degeneracy
  resonance  (|d| = c active)         -> the known stratum
Strata tl = 0 and tR = 0 (variables log10 beta, log10 t; column xi1^2 d/dtau for the vanishing time).
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm
from b1_mpcheck import mp_equilibrated_svals

BB = (-2.0, 2.0)
TB = (-3.0, 3.0)


def s_all(beta, tR, tl, metric="eq", pad=1.0):
    w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, pad=pad, nmax=60000)
    M = zm.real_rel_jac_np(w, 1.0, beta, tR, tl, 1.0, 0.0, "code", False, "log", wts)
    if metric == "eq":
        M = M / np.linalg.norm(M, axis=0)
    S = np.linalg.svd(M, compute_uv=False)
    return S, ok, len(w)


def point(x, stratum, sign):
    if stratum == "interior":
        lb, ltl, d = x
        tl = 10.0 ** ltl
        tR = tl * np.exp(d * sign)
        return 10.0 ** lb, tR, tl
    lb, lt = x
    return (10.0 ** lb, 10.0 ** lt, 0.0) if stratum == "tl=0" else (10.0 ** lb, 0.0, 10.0 ** lt)


def penalty(x, stratum, sign, c):
    beta, tR, tl = point(x, stratum, sign)
    pen = 0.0
    for v, (a, b) in ((np.log10(beta), BB),):
        pen += max(a - v, 0) ** 2 + max(v - b, 0) ** 2
    for t in (tR, tl):
        if t > 0:
            v = np.log10(t)
            pen += max(TB[0] - v, 0) ** 2 + max(v - TB[1], 0) ** 2
    if stratum == "interior":
        pen += max(c - x[2], 0) ** 2
    return pen


def objective(x, stratum, sign, c, metric):
    pen = penalty(x, stratum, sign, c)
    if pen > 0.0:                      # hard constraints: s4 -> 0 on the resonance would beat any soft penalty
        return 10.0 + 1e3 * pen
    beta, tR, tl = point(x, stratum, sign)
    S, ok, n = s_all(beta, tR, tl, metric)
    val = np.log10(max(S[-1], 1e-300))
    if not ok:
        val += 5.0
    return val + 1e3 * pen


def classify(x, stratum, sign, c):
    beta, tR, tl = point(x, stratum, sign)
    at_box = abs(np.log10(beta) - BB[0]) < 0.05 or abs(np.log10(beta) - BB[1]) < 0.05
    for t in (tR, tl):
        if t > 0:
            at_box = at_box or abs(np.log10(t) - TB[0]) < 0.05 or abs(np.log10(t) - TB[1]) < 0.05
    if stratum == "interior" and x[2] < c + 0.02:
        return "resonance"
    return "box" if at_box else "interior"


def run(stratum, c, metric, nstart, seed):
    rng = np.random.default_rng(seed)
    out = []
    for k in range(nstart):
        sign = 1.0 if k % 2 == 0 else -1.0
        if stratum == "interior":
            while True:
                x0 = np.array([rng.uniform(*BB), rng.uniform(*TB), rng.uniform(c, 6.0)])
                if penalty(x0, stratum, sign, c) == 0:
                    break
        else:
            x0 = np.array([rng.uniform(*BB), rng.uniform(*TB)])
        r = minimize(objective, x0, args=(stratum, sign, c, metric), method="Nelder-Mead",
                     options={"xatol": 1e-4, "fatol": 1e-5, "maxiter": 1500})
        beta, tR, tl = point(r.x, stratum, sign)
        out.append({"beta": beta, "tR": tR, "tl": tl, "log10_s4": float(r.fun),
                    "class": classify(r.x, stratum, sign, c)})
    return out


if __name__ == "__main__":
    t0 = time.time()
    nstart = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    res = {}
    for metric in ("eq", "raw"):
        for c in (0.1, 1.0):
            res[f"interior_c{c}_{metric}"] = run("interior", c, metric, nstart, seed=7)
        for stratum in ("tl=0", "tR=0"):
            res[f"{stratum}_{metric}"] = run(stratum, 0.0, metric, nstart // 2, seed=11)
    summary = {}
    for key, r in res.items():
        cls = [q["class"] for q in r]
        inter = sorted([q for q in r if q["class"] == "interior"], key=lambda q: q["log10_s4"])
        summary[key] = {"counts": {k: cls.count(k) for k in set(cls)},
                        "global_min": min(r, key=lambda q: q["log10_s4"]),
                        "best_interior": inter[:3]}
        print(key, summary[key]["counts"], "global min", summary[key]["global_min"])
        for q in inter[:3]:
            print("    interior:", q)
    checks = []
    for key, r in res.items():
        inter = sorted([q for q in r if q["class"] == "interior"], key=lambda q: q["log10_s4"])[:2]
        for q in inter:
            w, wts, ok = zm.resolved_grid(1.0, q["beta"], q["tR"], q["tl"], 1.0, pad=1.0, nmax=60000)
            if len(w) > 4000:
                continue
            ts = "log" if (q["tR"] > 0 and q["tl"] > 0) else "log"
            sv, svd, n = mp_equilibrated_svals(q["beta"], q["tR"], q["tl"], dps=50, wlist=w,
                                               equilibrate=key.endswith("eq"), tscale=ts)
            # mp_equilibrated_svals with wlist uses unit weights: recompute double with unit weights
            checks.append({"key": key, "point": q, "s_mp50_unitweights": sv, "s_double_unitweights": svd})
            print("mp check", key, "beta=%.4g tR=%.4g tl=%.4g" % (q["beta"], q["tR"], q["tl"]),
                  "s4 mp50 = %.6e  double = %.6e" % (sv[-1], svd[-1]))
    out = {"summary": summary, "runs": res, "mp_checks": checks, "seconds": time.time() - t0}
    with open(os.path.join(HERE, "data", "b1_adversarial.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("seconds", out["seconds"])
