"""Practical near-counterexamples in finite bands: explicit pairs theta != theta' with small misfit.

(a) Profile of the minimal misfit along a prescribed separation: for theta' constrained to
    ln xi1' = ln xi1 + D (D fixed), minimise the RMS relative misfit over the other three parameters.
    Done for the matched-effusivity film (beta = 1), a generic film (beta = 2), in the FDTR-like mid band
    [1e-4, 2]/xi1^2 and in the full band; the curve misfit(D) is the practical fibre.
(b) Band-limited near-resonance law: s4 versus d = ln(tR/tl) inside the mid band.
(c) Low-frequency gauge tau -> tau + c: misfit of theta' = (xi1, b, tR + c, tl + c) versus c and omega_b tau.
"""
import sys, os, json
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


def rms_misfit(tp, w, wts, Zt):
    Z = zm.Z_np(w, *tp)
    return float(np.sqrt(np.sum(wts * np.abs(Z / Zt - 1) ** 2) / np.sum(wts)))


def profile_xi(theta, w, wts, Ds, rng, nstart=12):
    Zt = zm.Z_np(w, *theta)
    sw = np.sqrt(wts / np.sum(wts))
    out = []
    for D in Ds:
        xi_p = theta[0] * np.exp(D)
        def res(x):
            b, tR, tl = np.exp(x)
            r = (zm.Z_np(w, xi_p, b, tR, tl) / Zt - 1) * sw
            return np.concatenate([r.real, r.imag])
        best = None
        for k in range(nstart):
            x0 = np.log([theta[1], max(theta[2], 1e-3), max(theta[3], 1e-3)]) + rng.normal(0, 1.5, 3)
            try:
                r = least_squares(res, x0, method="trf", xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=2000,
                                  bounds=(np.log([1e-3, 1e-6, 1e-6]), np.log([1e3, 1e5, 1e5])))
            except Exception:
                continue
            m = float(np.sqrt(2 * r.cost))
            if best is None or m < best[0]:
                best = (m, np.exp(r.x).tolist())
        out.append({"D": float(D), "xi1_p": float(xi_p), "misfit": best[0], "b_tR_tl": best[1]})
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(31)
    out = {}
    Ds = np.array([-1.5, -1.0, -0.5, -0.2, 0.2, 0.5, 1.0])
    for name, theta in (("beta1", (1.0, 1.0, 0.5, 0.1)), ("beta2", (1.0, 2.0, 0.5, 0.1)),
                        ("beta1_cattaneo_small_t", (1.0, 1.0, 0.05, 0.0))):
        for bname, (lo, hi) in (("mid", (1e-4, 2.0)), ("full", (None, None))):
            if lo is None:
                w, wts, ok = zm.resolved_grid(*theta, 1.0, pad=1.0)
            else:
                w, wts, ok = zm.resolved_grid(*theta, 1.0, wmin=lo, wmax=hi)
            prof = profile_xi(theta, w, wts, Ds, rng)
            out[f"{name}_{bname}"] = prof
            print(f"(a) {name:24s} {bname:4s}: " + "  ".join("D=%+.1f:%.1e" % (q["D"], q["misfit"]) for q in prof))
    # (b) band-limited near resonance
    rows = []
    for d in 10.0 ** np.arange(-1, -6.5, -1):
        theta = (1.0, 2.0, 0.05 * np.exp(d / 2), 0.05 * np.exp(-d / 2))
        w, wts, ok = zm.resolved_grid(*theta, 1.0, wmin=1e-4, wmax=2.0)
        M = zm.real_rel_jac_np(w, *theta, 1.0, 0.0, "code", False, "log", wts)
        S = np.linalg.svd(M, compute_uv=False) / np.sqrt(np.log(2e4))
        rows.append([d, *S.tolist()])
    rows = np.array(rows)
    k = np.polyfit(np.log(rows[1:, 0]), np.log(rows[1:, 4]), 1)[0]
    out["band_near_resonance"] = {"rows[d,s1..s4]": rows.tolist(), "slope_s4": float(k)}
    print("(b) mid band, beta=2, tbar=0.05: s4 ~ d^%.3f; s4/d = %.3e; s3 = %.3e" % (k, rows[-1, 4] / rows[-1, 0], rows[-1, 3]))
    # (c) low-frequency gauge tau -> tau + c
    rows = []
    for tau_scale in (1.0, 0.1, 0.01):
        theta = (1.0, 2.0, 0.5 * tau_scale, 0.1 * tau_scale)
        w, wts, ok = zm.resolved_grid(*theta, 1.0, wmin=1e-4, wmax=2.0)
        Zt = zm.Z_np(w, *theta)
        for c in (0.1, 0.5, 1.0):
            cc = c * theta[3]
            tp = (1.0, 2.0, theta[2] + cc, theta[3] + cc)
            rows.append([tau_scale, c, rms_misfit(tp, w, wts, Zt)])
            print("(c) omega_b tR = %.2f, shift c = %.1f tl: misfit = %.3e" % (2.0 * theta[2], c, rows[-1][-1]))
    out["shift_gauge"] = rows
    with open(os.path.join(HERE, "data", "practical_pairs.json"), "w") as fh:
        json.dump(out, fh, indent=1)
