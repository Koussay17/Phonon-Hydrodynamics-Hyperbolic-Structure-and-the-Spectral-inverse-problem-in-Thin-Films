"""B-2 / B-4, step 6: numerical test of the global (fibre) statements at function level.

For targets on the resonance line and at generic points, multistart least squares (log parameters, 60 frequencies
omega in [1e-2, 1e2], relative complex misfit) searches for a second parameter point with the same response:
  bare film, T observable: local (control: the resonance line is a continuum of exact twins), entropic, kinetic
  eps = 0.084 and eps = -0.5; transducer stack (6 unknowns): local (control) and entropic.
Reported per case: the best misfit among local minima farther than 1e-3 (log distance) from the target, and the
largest log distance among minima with misfit < 1e-10 (exact twins). Theorems C1/T2 predict exact twins only for
the local controls. Double precision. Output: b2_06_global_search.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import numpy as np
from scipy.optimize import least_squares

ws = np.logspace(-2, 2, 60)
XA, BA = 0.2, 3.0


def bare(th, w, eps, R=0.0):
    xi, b, tR, tl = th[:4]
    p = 1j * w
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(p)
    L = (R + 1 / S) / (1 + eps * tl * p)
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    return (1 + r * E) / ((1 - r * E) * B)


def stack(th, w, eps):
    xi, b, tR, tl, R1, R2 = th
    p = 1j * w
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(p)
    N = 1 + eps * tl * p
    L = (R2 + 1 / S) / N
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    W = R1 + N * (1 + r * E) / ((1 - r * E) * B)
    At, Bt = XA * np.sqrt(p), BA * np.sqrt(p)
    return (np.cosh(At) * W + np.sinh(At) / Bt) / (Bt * np.sinh(At) * W + np.cosh(At))


def run(model, th_star, eps, nstart, rng):
    th_star = np.asarray(th_star, float)
    f = (lambda th: bare(th, ws, eps)) if model == "bare" else (lambda th: stack(th, ws, eps))
    target = f(th_star)
    l0 = np.log(th_star)

    def res(lq):
        d = (f(np.exp(lq)) - target) / np.abs(target)
        return np.concatenate([d.real, d.imag])
    best_far, far_exact = np.inf, 0.0
    minima = []
    for _ in range(nstart):
        x0 = l0 + 1.2 * rng.standard_normal(len(l0))
        try:
            sol = least_squares(res, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=1500)
        except Exception:
            continue
        if not np.all(np.isfinite(sol.x)) or np.max(np.abs(sol.x)) > 30:
            continue
        mis = float(np.sqrt(np.mean(sol.fun ** 2)))
        dist = float(np.max(np.abs(sol.x - l0)))
        minima.append((mis, dist, [float(v) for v in np.exp(sol.x)]))
        if dist > 1e-3:
            best_far = min(best_far, mis)
        if mis < 1e-10:
            far_exact = max(far_exact, dist)
    minima.sort()
    return dict(target=list(th_star), eps=eps, best_misfit_far_from_target=best_far,
                largest_logdist_of_exact_minima=far_exact,
                best_far_minima=[m for m in minima if m[1] > 1e-3][:3])


rng = np.random.default_rng(11)
cases = {
    "bare/local Delta (control)": ("bare", (0.8, 2.0, 0.4, 0.4), 0.0),
    "bare/entropic Delta": ("bare", (0.8, 2.0, 0.4, 0.4), 1.0),
    "bare/entropic Delta b=bs": ("bare", (0.8, 1.0, 0.4, 0.4), 1.0),
    "bare/entropic generic": ("bare", (0.8, 2.0, 0.5, 0.1), 1.0),
    "bare/kinetic eps=0.084 Delta": ("bare", (0.8, 2.0, 0.4, 0.4), 0.084),
    "bare/kinetic eps=-0.5 Delta": ("bare", (0.8, 2.0, 0.4, 0.4), -0.5),
    "stack/local Delta (control)": ("stack", (0.8, 2.0, 0.4, 0.4, 0.3, 0.2), 0.0),
    "stack/entropic Delta": ("stack", (0.8, 2.0, 0.4, 0.4, 0.3, 0.2), 1.0),
    "stack/entropic Delta R2=0": ("stack", (0.8, 2.0, 0.4, 0.4, 0.3, 1e-12), 1.0),
    "stack/entropic generic": ("stack", (0.8, 2.0, 0.5, 0.1, 0.3, 0.2), 1.0),
}
OUT = {}
for name, (model, th, eps) in cases.items():
    OUT[name] = run(model, th, eps, 100, rng)
    print(f"{name:32s} best far misfit={OUT[name]['best_misfit_far_from_target']:.2e} "
          f"exact-twin logdist={OUT[name]['largest_logdist_of_exact_minima']:.2e}", flush=True)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_06_global_search.json"), "w") as fh:
        json.dump(OUT, fh, indent=1, default=float)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b2_06_global_search.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, default=float)
