"""B1: scaling of the smallest singular values near the degenerate strata.

(a) Fourier resonance: tR = tbar e^{d/2}, tl = tbar e^{-d/2}, d -> 0.  Expect s4 ~ |d|^k; fit k.
(b) Invisible film: beta = 1 + eps and resonance (d = 0) or both small: s3, s4 versus (eps, d).
(c) Exact strata at 50 digits: rank 3 on resonance (beta != 1), rank 2 on {resonance, beta = 1}.
Metrics: raw singular values of the log-parameter relative Jacobian with weights d ln(omega)
(inverse standard deviations at unit relative noise per e-fold of frequency), and the column-
equilibrated ones (collinearity).  Window: every characteristic frequency +- pad decades.
"""
import sys, os, json
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm
from b1_mpcheck import mp_equilibrated_svals


def svals(beta, tR, tl, pad=1.0, ppd=24, variant="code"):
    w, wts = zm.wide_grid(1.0, beta, tR, tl, 1.0, pad=pad, ppd=ppd)
    M = zm.real_rel_jac_np(w, 1.0, beta, tR, tl, 1.0, 0.0, variant, False, "log", wts)
    Sraw = np.linalg.svd(M, compute_uv=False)
    Mn = M / np.linalg.norm(M, axis=0)
    Sn = np.linalg.svd(Mn, compute_uv=False)
    return Sraw, Sn


out = {"a_resonance": [], "b_invisible": [], "c_exact": []}
ds = 10.0 ** np.arange(-1, -9.5, -1)
for beta in (0.1, 0.5, 2.0, 10.0):
    for tbar in (1e-2, 1.0, 1e2):
        rec = {"beta": beta, "tbar": tbar, "d": ds.tolist(), "s4raw": [], "s4n": [], "s3raw": []}
        for d in ds:
            Sraw, Sn = svals(beta, tbar * np.exp(d / 2), tbar * np.exp(-d / 2))
            rec["s4raw"].append(float(Sraw[-1])); rec["s3raw"].append(float(Sraw[-2]))
            rec["s4n"].append(float(Sn[-1]))
        k = np.polyfit(np.log(ds[2:]), np.log(rec["s4raw"][2:]), 1)[0]
        rec["slope_s4raw"] = float(k)
        rec["coef_s4raw_over_d"] = float(rec["s4raw"][-1] / ds[-1])
        out["a_resonance"].append(rec)
        print(f"(a) beta={beta:5g} tbar={tbar:6g}: s4raw/|d| at d=1e-9: {rec['coef_s4raw_over_d']:.4e}, "
              f"slope {k:.4f}, s3raw={rec['s3raw'][-1]:.3e}")

# (b) invisible film: beta = 1 + eps, d small
for tbar in (1e-2, 1.0, 1e2):
    grid = []
    for eps in (1e-1, 1e-2, 1e-3, 1e-4):
        for d in (1e-1, 1e-2, 1e-3, 1e-4):
            Sraw, Sn = svals(1.0 + eps, tbar * np.exp(d / 2), tbar * np.exp(-d / 2))
            grid.append([eps, d, float(Sraw[-2]), float(Sraw[-1])])
    g = np.array(grid)
    # fit log s = a log eps + c log d + const for s3 and s4
    X = np.column_stack([np.log(g[:, 0]), np.log(g[:, 1]), np.ones(len(g))])
    c3 = np.linalg.lstsq(X, np.log(g[:, 2]), rcond=None)[0]
    c4 = np.linalg.lstsq(X, np.log(g[:, 3]), rcond=None)[0]
    out["b_invisible"].append({"tbar": tbar, "grid[eps,d,s3raw,s4raw]": g.tolist(),
                               "fit_s3_exponents[eps,d]": c3[:2].tolist(),
                               "fit_s4_exponents[eps,d]": c4[:2].tolist()})
    print(f"(b) tbar={tbar:g}: s3 ~ eps^{c3[0]:.3f} d^{c3[1]:.3f};  s4 ~ eps^{c4[0]:.3f} d^{c4[1]:.3f}")
    for row in g:
        print("     eps=%.0e d=%.0e  s3raw=%.3e s4raw=%.3e" % tuple(row))

# (c) exact strata at 50 digits.  Raw singular values (no equilibration: an identically zero column,
# as d/dln xi1 on the invisible-film set, would be blown up to O(1) noise by equilibration).
for (beta, tR, tl) in ((2.0, 0.7, 0.7), (0.3, 5.0, 5.0), (1.0, 0.7, 0.7), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                       (1.0, 0.7, 0.3), (1.0, 0.7, 0.0), (1.0, 0.0, 0.7)):
    sv, svd, n = mp_equilibrated_svals(beta, tR, tl, pad=1.0, ppd=8, dps=50, tscale="lin",
                                       equilibrate=False)
    out["c_exact"].append({"beta": beta, "tR": tR, "tl": tl, "svals_raw_mp50": sv})
    print(f"(c) beta={beta} tR={tR} tl={tl}: raw s (50 digits, lin tau columns) =",
          ["%.3e" % x for x in sv])

# (d) beta = 1 exactly, d -> 0 (approach to the invisible-film set along matched effusivity)
out["d_beta1"] = []
for tbar in (1e-2, 1.0, 1e2):
    rec = {"tbar": tbar, "d": [], "s": []}
    for d in 10.0 ** np.arange(-1, -7.5, -1):
        Sraw, Sn = svals(1.0, tbar * np.exp(d / 2), tbar * np.exp(-d / 2))
        rec["d"].append(float(d)); rec["s"].append(Sraw.tolist())
    s = np.array(rec["s"]); dd = np.array(rec["d"])
    rec["slope_s3"] = float(np.polyfit(np.log(dd[2:]), np.log(s[2:, 2]), 1)[0])
    rec["slope_s4"] = float(np.polyfit(np.log(dd[2:]), np.log(s[2:, 3]), 1)[0])
    out["d_beta1"].append(rec)
    print(f"(d) beta=1, tbar={tbar:g}: s3 ~ d^{rec['slope_s3']:.3f}, s4 ~ d^{rec['slope_s4']:.3f}; "
          f"at d=1e-7: s3={s[-1,2]:.3e}, s4={s[-1,3]:.3e}")
# (e) on resonance, beta = 1 + eps -> 1: s3 ~ eps^k
out["e_res_eps"] = []
for tbar in (1e-2, 1.0, 1e2):
    rec = {"tbar": tbar, "eps": [], "s3": []}
    for eps in 10.0 ** np.arange(-1, -7.5, -1):
        Sraw, Sn = svals(1.0 + eps, tbar, tbar)
        rec["eps"].append(float(eps)); rec["s3"].append(float(Sraw[2]))
    rec["slope"] = float(np.polyfit(np.log(rec["eps"][2:]), np.log(rec["s3"][2:]), 1)[0])
    out["e_res_eps"].append(rec)
    print(f"(e) resonance tbar={tbar:g}: s3 ~ eps^{rec['slope']:.3f}; s3/eps at 1e-7 = {rec['s3'][-1]/1e-7:.4e}")

with open(os.path.join(HERE, "data", "b1_resonance_scaling.json"), "w") as fh:
    json.dump(out, fh, indent=1)
