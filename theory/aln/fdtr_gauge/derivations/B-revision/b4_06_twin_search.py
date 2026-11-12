"""B-4, step 6: global uniqueness with finitely many frequencies (known transducer, unknown R1, R2; 6 unknowns).

Multistart root collection for Z_tot(omega_j; theta) = Z_tot(omega_j; theta*), j = 1..N (2N real equations),
unknowns theta = (xi, b, tau_R, tau_l, R1, R2) > 0 in log coordinates, N = 3 (square system) and N = 4
(overdetermined). Levenberg-Marquardt (scipy least_squares) from wide log-normal starts; a root is accepted
at relative residual < 1e-11 and distinct roots are those farther than 1e-5 in log distance from the target.
Double precision is enough to locate roots; candidate twins are refined and certified separately
(b4_07_certify_twin.py). Local and entropic conventions. Nondimensional: b_s = 1, transducer xi_A = 0.2, b_A = 3.
Output: b4_06_twin_search.json
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import json
import numpy as np
from scipy.optimize import least_squares

XA, BA = 0.2, 3.0


def ztot(th, w, conv):
    xi, b, tR, tl, R1, R2 = th
    p = 1j * w
    g = (1 + tR * p) / (1 + tl * p)
    A = xi * np.sqrt(p * g)
    B = b * np.sqrt(p / g)
    S = np.sqrt(p)
    N = 1 + (tl * p if conv == "entropic" else 0)
    L = (R2 + 1 / S) / N
    E = np.exp(-2 * A)
    r = (B * L - 1) / (B * L + 1)
    Zf = N * (1 + r * E) / ((1 - r * E) * B)        # N * Mobius(M; L2/N)
    W = R1 + Zf
    At, Bt = XA * np.sqrt(p), BA * np.sqrt(p)
    return (np.cosh(At) * W + np.sinh(At) / Bt) / (Bt * np.sinh(At) * W + np.cosh(At))


def residual(lq, target, ws, conv):
    th = np.exp(lq)
    z = ztot(th, ws, conv)
    d = (z - target) / np.abs(target)
    return np.concatenate([d.real, d.imag])


def search(th_star, ws, conv, nstart, rng, spread=1.5):
    ws = np.asarray(ws, float)
    target = ztot(np.asarray(th_star, float), ws, conv)
    l0 = np.log(th_star)
    roots = []
    for _ in range(nstart):
        x0 = l0 + spread * rng.standard_normal(6)
        try:
            sol = least_squares(residual, x0, args=(target, ws, conv), method="lm", xtol=1e-15, ftol=1e-15,
                                gtol=1e-15, max_nfev=1500)
        except Exception:
            continue
        if not np.all(np.isfinite(sol.x)) or np.max(np.abs(sol.x)) > 40:
            continue
        rr = np.max(np.abs(sol.fun))
        if rr < 1e-11:
            if all(np.max(np.abs(sol.x - q)) > 1e-5 for q in roots):
                roots.append(sol.x.copy())
    out = []
    for q in roots:
        out.append(dict(theta=[float(v) for v in np.exp(q)], logdist_to_target=float(np.max(np.abs(q - l0))),
                        is_target=bool(np.max(np.abs(q - l0)) < 1e-5)))
    return out


rng = np.random.default_rng(20261112)
targets = {
    "generic": (0.8, 2.0, 0.5, 0.1, 0.3, 0.2),
    "b=bs": (0.8, 1.0, 0.5, 0.1, 0.3, 0.2),
    "AlN-like": (1.0, 2.07, 0.03, 0.012, 0.05, 0.1),
}
triples = {"N3a": (0.5, 4.0, 1.7), "N3b": (0.2, 1.1, 6.0)}
quads = {"N4a": (0.5, 4.0, 1.7, 0.2)}
OUT = {}
for conv in ("local", "entropic"):
    for tn, th in targets.items():
        for wn, ws in list(triples.items()) + list(quads.items()):
            key = f"{conv}:{tn}:{wn}"
            roots = search(th, ws, conv, 150 if wn.startswith("N3") else 80, rng)
            OUT[key] = dict(frequencies=list(ws), target=list(th), roots=roots,
                            n_distinct_roots=len(roots), n_twins=sum(1 for r in roots if not r["is_target"]))
            print(key, "roots:", len(roots), "twins:", OUT[key]["n_twins"],
                  [[round(v, 5) for v in r["theta"]] for r in roots if not r["is_target"]][:3], flush=True)
            with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_06_twin_search.json"), "w") as fh:
                json.dump(OUT, fh, indent=1)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_06_twin_search.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
