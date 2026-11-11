"""B1: exact rank on the strata, random points, 50-digit arithmetic, raw singular values (no equilibration).

Strata (xi1 = bs = 1): Cattaneo tl = 0 < tR; GK without relaxation tR = 0 < tl; matched effusivity beta = 1
off resonance; high-frequency impedance matching beta = sqrt(tR/tl) (rho -> 0 as omega -> inf);
Fourier resonance tR = tl (beta != 1); invisible film {tR = tl, beta = 1}.
The vanishing time uses the column xi1^2 d/dtau (one-sided derivative; the map is analytic across tau = 0
on the imaginary axis).  Grid: characteristic frequencies +- 1 decade, at most ~120 points.
Exact rank r  <=>  s_{r+1}/s_1 at the level of 10^-50 while s_r/s_1 is not.
"""
import sys, os, json
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


def raw_svals(beta, tR, tl, dps=50, nmax=120):
    w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, pad=1.0, nmax=60000)
    if len(w) > nmax:
        idx = np.unique(np.linspace(0, len(w) - 1, nmax).astype(int))
        w = w[idx]
    mp.mp.dps = dps
    M = zm.real_rel_jac_mp([mp.mpf(x) for x in w], mp.mpf(1), mp.mpf(beta), mp.mpf(tR), mp.mpf(tl), 1, 0,
                           "code", False, "log", None)
    sv = zm.svals_mp(M)
    return [float(x / sv[0]) for x in sv]


if __name__ == "__main__":
    rng = np.random.default_rng(29)
    out = {}
    def draw():
        return 10 ** rng.uniform(-1.5, 1.5), 10 ** rng.uniform(-2, 2), 10 ** rng.uniform(-2, 2)
    strata = {
        "cattaneo tl=0": lambda b, t1, t2: (b, t1, 0.0),
        "gk tR=0": lambda b, t1, t2: (b, 0.0, t2),
        "beta=1 off resonance": lambda b, t1, t2: (1.0, t1, t2),
        "beta=sqrt(tR/tl)": lambda b, t1, t2: (float(np.sqrt(t1 / t2)), t1, t2),
        "resonance tR=tl, beta!=1": lambda b, t1, t2: (b, t1, t1),
        "invisible film tR=tl, beta=1": lambda b, t1, t2: (1.0, t1, t1),
    }
    for name, f in strata.items():
        rows = []
        for k in range(30):
            b, t1, t2 = draw()
            beta, tR, tl = f(b, t1, t2)
            sv = raw_svals(beta, tR, tl)
            rows.append([beta, tR, tl] + sv)
        rows = np.array(rows)
        rec = {"n": len(rows), "min_s2": float(rows[:, 4].min()), "min_s3": float(rows[:, 5].min()),
               "min_s4": float(rows[:, 6].min()), "max_s3": float(rows[:, 5].max()), "max_s4": float(rows[:, 6].max()),
               "rows": rows.tolist()}
        out[name] = rec
        print("%-30s s3/s1 in [%.2e, %.2e], s4/s1 in [%.2e, %.2e]" % (name, rec["min_s3"], rec["max_s3"],
                                                                     rec["min_s4"], rec["max_s4"]))
    with open(os.path.join(HERE, "data", "b1_strata_mp.json"), "w") as fh:
        json.dump(out, fh, indent=1)
