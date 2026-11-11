"""B2: refinement of the best distinct-parameter pairs found by b2_injectivity.py.

For each selected (target, band, stratum) minimum theta':
  1. re-minimise from theta' with tighter tolerances (does the misfit decrease further or stall?);
  2. recompute the misfit in 40-digit arithmetic (round-off excluded);
  3. evaluate the same theta' on the full band and on bands widened by 1 and 2 decades on each side
     (an exact counterexample would keep misfit ~ 0 on every band).
All with the physical (transmission-line) impedance, variant "code".
"""
import sys, os, json
import numpy as np
import mpmath as mp
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm
from b2_injectivity import TARGETS, bands, residual_and_jac, unpack, distance

d = json.load(open(os.path.join(HERE, "data", "b2_injectivity.json")))


def rms(tp, w, wts, Zt):
    Z = zm.Z_np(w, *tp)
    return float(np.sqrt(np.sum(wts * np.abs(Z / Zt - 1) ** 2) / np.sum(wts)))


def rms_mp(tp, theta, w, wts, dps=40):
    mp.mp.dps = dps
    num = mp.mpf(0); den = mp.mpf(0)
    for wi, wt in zip(w, wts):
        z1 = zm.Z_mp(mp.mpf(wi), *[mp.mpf(v) for v in tp])
        z0 = zm.Z_mp(mp.mpf(wi), *[mp.mpf(v) for v in theta])
        num += mp.mpf(wt) * abs(z1 / z0 - 1) ** 2
        den += mp.mpf(wt)
    return float(mp.sqrt(num / den))


def stratum_of(tp):
    if tp[2] > 0 and tp[3] > 0:
        return "gk"
    if tp[2] > 0:
        return "cat"
    if tp[3] > 0:
        return "gk0"
    return "fourier"


picks = [("T1_gk", "fdtr_thin", "cat"), ("T1_gk", "fdtr_mid", "cat"), ("T2_gk_tl>tR", "fdtr_mid", "gk0"),
         ("T4_beta1", "fdtr_mid", "gk"), ("T4_beta1", "full", "gk"), ("T5_near_res", "fdtr_mid", "cat"),
         ("T5_near_res", "full", "gk"), ("T1_gk", "full", "gk"), ("T8_gk0", "fdtr_thin", "gk")]
out = []
for tname, bname, st in picks:
    theta = TARGETS[tname]
    B = bands(theta)
    w, wts = B[bname]
    Zt = zm.Z_np(w, *theta)
    sw = np.sqrt(wts / np.sum(wts))
    far = d[tname][bname][st]["best_far"]
    if not far:
        continue
    tp = far[0]["theta_p"]
    stp = stratum_of(tp)
    x0 = np.log([v for v in tp if v > 0])
    f = lambda x: residual_and_jac(x, stp, w, sw, Zt)[0]
    g = lambda x: residual_and_jac(x, stp, w, sw, Zt)[1]
    r = least_squares(f, x0, jac=g, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    tp2 = unpack(r.x, stp)
    m_band = rms(tp2, w, wts, Zt)
    m_band_mp = rms_mp(tp2, theta, w[::max(1, len(w) // 60)], wts[::max(1, len(w) // 60)])
    wf, wtf = B["full"]
    m_full = rms(tp2, wf, wtf, zm.Z_np(wf, *theta))
    widened = []
    for k in (1, 2):
        lo, hi = w[0] / 10 ** k, w[-1] * 10 ** k
        ww, wwt, ok = zm.resolved_grid(*theta, 1.0, wmin=lo, wmax=hi, nmax=60000)
        widened.append(rms(tp2, ww, wwt, zm.Z_np(ww, *theta)))
    rec = {"target": tname, "theta": list(theta), "band": bname, "band_limits": [float(w[0]), float(w[-1])],
           "theta_p": [float(v) for v in tp2], "dist": distance(tp2, theta),
           "misfit_band_initial": far[0]["misfit"], "misfit_band_refined": m_band,
           "misfit_band_40digits": m_band_mp, "misfit_full_band": m_full,
           "misfit_band_widened_1dec": widened[0], "misfit_band_widened_2dec": widened[1]}
    out.append(rec)
    print("%-12s %-10s %-4s theta'=%s dist=%.2f | band %.3e -> %.3e (40 dig %.3e) | widened 1 dec %.2e, 2 dec %.2e | full %.2e" % (
        tname, bname, stp, np.array2string(np.array(tp2), precision=4), rec["dist"], far[0]["misfit"], m_band,
        m_band_mp, widened[0], widened[1], m_full))
with open(os.path.join(HERE, "data", "b2_refine_pairs.json"), "w") as fh:
    json.dump(out, fh, indent=1)
