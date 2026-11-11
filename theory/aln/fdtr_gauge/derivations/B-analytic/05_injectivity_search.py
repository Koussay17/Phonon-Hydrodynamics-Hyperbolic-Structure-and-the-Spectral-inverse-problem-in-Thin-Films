"""Workstream B: adversarial search for two distinct parameter points with the same response.

(1) Function-level data (60 frequencies, 1e-2 .. 1e2 in units xi_1^2 = 1, b_s = 1): for target
    points on several strata, a multistart least-squares fit (log-parameters, 400 starts) looks for
    local minima with tiny misfit far from the target. Off the resonance none should exist.
(2) Two-frequency data (4 real equations, 4 unknowns): multistart root finding collects the
    distinct exact roots; several distinct roots would show that N = 2 gives local but not global
    identifiability (numerical observation only).
(3) The same search with the interface resistance R (5 parameters).
Run:  python 05_injectivity_search.py > 05_injectivity_search.out
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(20261111)
BS = 1.0


def Z(w, xi, b, tR, tl, R=0.0):
    p = 1j*np.asarray(w)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*np.sqrt(p*g)
    B = b*np.sqrt(p/g)
    S = BS*np.sqrt(p)
    L = R + 1/S
    r = (B*L - 1)/(B*L + 1)
    E = np.exp(-2*A)
    return (1/B)*(1 + r*E)/(1 - r*E)


def resid(x, w, data, nR):
    th = np.exp(x)
    z = Z(w, *th[:4], R=(th[4] if nR else 0.0))
    d = np.log(z) - np.log(data)          # relative amplitude and phase misfit
    return np.concatenate([d.real, d.imag])


def search(theta, w, nR=False, nstart=400, tag=""):
    data = Z(w, *theta[:4], R=(theta[4] if nR else 0.0))
    x0 = np.log(np.array(theta))
    found = []
    for k in range(nstart):
        xs = x0 + rng.normal(0, 1.5, size=len(x0))
        sol = least_squares(resid, xs, args=(w, data, nR), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        cost = np.sqrt(np.mean(sol.fun**2))
        found.append((cost, np.exp(sol.x)))
    found.sort(key=lambda t: t[0])
    th = np.array(theta)
    far = [(c, x) for c, x in found if np.max(np.abs(np.log(x) - np.log(th))) > 1e-3]
    near = [(c, x) for c, x in found if np.max(np.abs(np.log(x) - np.log(th))) <= 1e-3]
    print(f"[{tag}] target {np.round(th, 4)}: {len(near)} starts converged to the target (best rms {near[0][0] if near else float('nan'):.2e});"
          f" best rms among minima farther than 1e-3 in log-distance: {far[0][0] if far else float('nan'):.3e}")
    for c, x in far[:3]:
        extra = "" if not nR else f", R={x[4]:.4g}"
        print(f"      rms {c:.3e} at xi={x[0]:.6g}, b={x[1]:.6g}, tau_R={x[2]:.6g}, tau_l={x[3]:.6g}{extra};"
              f"  |tau_R - tau_l|/(tau_R + tau_l) = {abs(x[2]-x[3])/(x[2]+x[3]):.3e}")
    return found


wf = np.logspace(-2, 2, 60)
print("=== (1) function-level data, 4 parameters ===")
for tag, th in [("interior", (1.0, 2.0, 0.7, 0.3)), ("b = b_s", (1.0, 1.0, 0.7, 0.3)),
                ("near-Cattaneo tau_l = 1e-3", (1.0, 2.0, 0.7, 1e-3)), ("near tau_R = 0", (1.0, 2.0, 1e-3, 0.5)),
                ("tau_l > tau_R", (0.8, 0.5, 0.2, 0.9)), ("resonance (expected: continuum)", (1.0, 2.0, 0.5, 0.5))]:
    search(th, wf, tag=tag)

print("\n=== (2) two-frequency data: distinct exact roots ===")
for tag, th, w2 in [("interior", (1.0, 2.0, 0.7, 0.3), np.array([0.5, 4.0])),
                    ("b = b_s", (1.0, 1.0, 0.7, 0.3), np.array([0.5, 4.0])),
                    ("AlN-like, 10 and 100 MHz", (1.0, 2.70, 0.532, 0.266), np.array([0.118, 1.18]))]:
    found = search(th, w2, nstart=600, tag="N=2 " + tag)
    roots = []
    for c, x in found:
        if c < 1e-11 and not any(np.max(np.abs(np.log(x) - np.log(y))) < 1e-6 for y in roots):
            roots.append(x)
    print(f"      distinct roots with rms < 1e-11: {len(roots)}")
    for x in roots[:8]:
        print("        ", np.array2string(x, precision=8))

print("\n=== (3) function-level data, 5 parameters (interface resistance) ===")
for tag, th in [("interior R > 0", (1.0, 2.0, 0.7, 0.3, 0.3)), ("tau_R small, R > 0", (1.0, 2.0, 1e-3, 0.5, 0.3)),
                ("b = b_s, R > 0", (1.0, 1.0, 0.7, 0.3, 0.3)), ("small R", (1.0, 2.0, 0.7, 0.3, 1e-3))]:
    search(th, wf, nR=True, tag=tag)
