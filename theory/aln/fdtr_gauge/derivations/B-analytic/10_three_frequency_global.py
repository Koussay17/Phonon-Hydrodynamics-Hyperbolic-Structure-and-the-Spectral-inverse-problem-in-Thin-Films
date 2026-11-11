"""Workstream B: global uniqueness from finitely many frequencies (numerical exploration only).

For N = 2 the certified twin of 09 shows that local identifiability does not imply uniqueness of the
two-frequency inverse problem. Here a multistart search (wide log-normal starts, 800 per case) collects
the distinct exact roots (rms < 1e-11) of the N = 2 and N = 3 real equations for several targets.
A result "one root" is evidence, not proof, of uniqueness within the explored region.
Run:  python 10_three_frequency_global.py > 10_three_frequency_global.out
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(31415)


def Z(w, xi, b, tR, tl, bs=1.0):
    p = 1j*np.asarray(w)
    g = (1 + tR*p)/(1 + tl*p)
    A = xi*np.sqrt(p*g); B = b*np.sqrt(p/g); S = bs*np.sqrt(p)
    r = (B - S)/(B + S); E = np.exp(-2*A)
    return (1/B)*(1 + r*E)/(1 - r*E)


def roots(theta, w, nstart=300, spread=2.0):
    data = Z(w, *theta)

    def res(x):
        d = np.log(Z(w, *np.exp(x))) - np.log(data)
        return np.concatenate([d.real, d.imag])
    x0 = np.log(theta)
    found = []
    for _ in range(nstart):
        sol = least_squares(res, x0 + rng.normal(0, spread, 4), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        if np.sqrt(np.mean(sol.fun**2)) < 1e-11:
            x = np.exp(sol.x)
            if not any(np.max(np.abs(np.log(x) - np.log(y))) < 1e-6 for y in found):
                found.append(x)
    return found


cases = [("interior", (1.0, 2.0, 0.7, 0.3)), ("b = b_s", (1.0, 1.0, 0.7, 0.3)),
         ("tau_l > tau_R", (0.8, 0.5, 0.2, 0.9)), ("AlN-like", (1.0, 2.70, 0.532, 0.266))]
# N = 2 with (0.5, 4) and 800 starts was run first; its partial record is 10a_two_frequency_roots_partial.out.
for wlabel, w in [("N=3 (0.5, 4, 1.7)", np.array([0.5, 4.0, 1.7])), ("N=3 (0.2, 1.1, 6)", np.array([0.2, 1.1, 6.0])),
                  ("N=2 (0.2, 1.1)", np.array([0.2, 1.1]))]:
    print(f"=== {wlabel} ===")
    for label, th in cases:
        fr = roots(np.array(th), w)
        print(f"  {label}: {len(fr)} distinct exact root(s)")
        for x in fr:
            print("      ", np.array2string(x, precision=7))
