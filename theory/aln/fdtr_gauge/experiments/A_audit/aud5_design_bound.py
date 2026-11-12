"""Audit, part 5: how much could the frequency design itself lower w*?

For a fixed number n = 80 of measurements (same per-point noise), allocated with weights xi_i >= 0, sum xi_i = 1,
over the candidate frequencies of the band (the 80 log-spaced points, repetitions allowed), the marginal variance of
ln tau_R is c^T M(xi)^-1 c / n, M(xi) = sum_i xi_i A_i^T A_i (A_i: whitened 2 x p block of frequency i). c-optimality
is convex in xi; it is minimised here by exponentiated-gradient (mirror) descent on the simplex, with the
Kiefer-Wolfowitz-type equivalence check max_i d_i / phi <= 1 + tol at the optimum (d_i = |A_i M^-1 c|^2).
For N0 the optimum is closed form: all n points at the frequency of largest per-point information.
Jacobians: the independent model (out/aud_julia_core.json).
Usage: python aud5_design_bound.py -> out/aud5_design_bound.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
JL = json.loads((OUT / "aud_julia_core.json").read_text())
F = np.array(JL["freq80"])
PN = ["ln_tauR", "ln_lamz", "ln_lamr", "ln_C", "ln_R1", "ln_R2", "ln_lamAu", "ln_CAu", "ln_dAu", "ln_dAlN",
      "ln_lamS", "ln_CS", "ln_wscale"]
SETS = {"N0": ["ln_tauR"], "N1": ["ln_tauR", "g", "phi0", "td_ns"]}
SETS["N2"] = SETS["N1"] + ["ln_R1", "ln_R2"]
SETS["N3"] = SETS["N2"] + ["ln_lamz", "ln_C"]


def cvec(a):
    a = np.asarray(a, float)
    return a[..., 0] + 1j * a[..., 1]


def blocks(J, names, sa, sp_deg):
    cols = []
    for n in names:
        if n == "g":
            cols.append(np.ones(F.size, complex))
        elif n == "phi0":
            cols.append(np.full(F.size, 1j))
        elif n == "td_ns":
            cols.append(-2j * np.pi * F * 1e-9)
        else:
            cols.append(J[:, PN.index(n)])
    D = np.stack(cols, axis=1)
    return np.stack([D.real / sa, D.imag / np.radians(sp_deg)], axis=1)      # (nf, 2, p)


def c_opt(A, n=80, iters=400000, tol=1e-6):
    """Exponentiated-gradient descent in log weights; returns sigma, equivalence gap, iterations, support.
    Convexity gives phi* >= phi (1 - gap), gap = max_i d_i / phi - 1 (certificate of near-optimality)."""
    nf, _, p = A.shape
    G = np.einsum("fki,fkj->fij", A, A)
    c = np.zeros(p); c[0] = 1.0
    s = np.sqrt(np.einsum("fii->i", G).clip(1e-300)) / np.sqrt(nf)      # column scaling (conditioning only)
    A_s = A / s
    G = np.einsum("fki,fkj->fij", A_s, A_s)
    c_s = c / s
    lw = np.full(nf, -np.log(nf))
    for it in range(iters):
        xi = np.exp(lw - lw.max()); xi /= xi.sum()
        M = np.einsum("f,fij->ij", xi, G)
        u = np.linalg.solve(M, c_s)
        phi = float(c_s @ u)
        d = (np.einsum("fki,i->fk", A_s, u) ** 2).sum(axis=1)          # |A_i M^-1 c|^2 = -d phi / d xi_i
        r = d / phi
        gap = float(r.max() - 1.0)
        if gap < tol:
            break
        lw = lw + (0.3 / r.max()) * (r - 1.0)
        lw = np.maximum(lw, lw.max() - 700.0)
    sigma = float(np.sqrt(phi / n))
    sigma_lower = float(np.sqrt(max(phi * (1.0 - gap), 0.0) / n))
    support = [(float(F[i]), float(xi[i])) for i in np.argsort(-xi)[:8] if xi[i] > 1e-4]
    return sigma, gap, it, support, sigma_lower


res = {}
for Lk in ("100", "300", "500"):
    J = np.stack([cvec(col) for col in JL["cases"][Lk]["J_all"]], axis=1)
    for N in SETS:
        for sp in (0.01, 0.1):
            A = blocks(J, SETS[N], 0.01, sp)
            # log-spaced design of the branch (all weights 1/80)
            Fm = np.einsum("fki,fkj->ij", A, A)
            s_ls = float(np.sqrt(np.linalg.inv(Fm)[0, 0]))
            if N == "N0":
                info = (A[:, :, 0] ** 2).sum(axis=1)
                s_opt = float(1 / np.sqrt(80 * info.max()))
                gap, it, sup, s_low = 0.0, 0, [(float(F[info.argmax()]), 1.0)], s_opt
            else:
                s_opt, gap, it, sup, s_low = c_opt(A)
            res[f"{Lk}|{N}|{sp}"] = {"sigma_logspaced80": s_ls, "sigma_c_optimal80": s_opt, "sigma_lower_bound": s_low,
                                     "gain": s_ls / s_opt, "equivalence_gap": gap, "iterations": it, "support_top": sup}
            print(Lk, N, sp, "log-spaced %.4g  c-optimal %.4g (>= %.4g)  gain %.2f  gap %.1e  it %d  support %s" % (
                s_ls, s_opt, s_low, s_ls / s_opt, gap, it, [(round(f / 1e6, 3), round(w, 3)) for f, w in sup[:5]]),
                flush=True)
(OUT / "aud5_design_bound.json").write_text(json.dumps(res, indent=1))
