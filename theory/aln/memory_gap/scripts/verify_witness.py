"""Independent verification of a witness g: four evaluations of K, N, tau and constraint residuals.

(a) dense Cholesky of C + gamma e_hat e_hat^T (float64);
(b) dense eigendecomposition of C restricted to H (float64), spectral sums, lambda_min on H;
(c) matrix-free conjugate gradients with the event action y -> sum g a (a^T y), projected on H;
(d) ball arithmetic (python-flint arb, rigorous enclosures) if n <= nmax_arb.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import scipy.linalg as sla
from scipy.sparse.linalg import LinearOperator, cg

sys.dont_write_bytecode = True


def chol(geom, g, col):
    R = geom.response(g, col)
    return {"K": R["K"], "N": R["N"], "tau": R["tau"]}


def eig(geom, g, col):
    C = geom.C(g)
    w, U = np.linalg.eigh(C)
    c = U.T @ geom.b[:, col]
    # drop the energy direction (largest overlap with e_hat)
    ov = np.abs(U.T @ geom.ehat)
    k0 = int(np.argmax(ov))
    keep = np.ones(len(w), bool)
    keep[k0] = False
    K = float(np.sum(c[keep] ** 2 / w[keep]))
    N = float(np.sum(c[keep] ** 2 / w[keep] ** 2))
    return {"K": K, "N": N, "tau": N / K, "lambda_min_H": float(w[keep].min()), "lambda_max": float(w.max()),
            "cond_H": float(w.max() / w[keep].min()), "energy_mode_eig": float(w[k0]), "energy_overlap": float(ov[k0])}


def cgsolve(geom, g, col, rtol=1e-14):
    A = geom.A.tocsr()
    AT = A.T.tocsr()
    eh = geom.ehat

    def mv(y):
        y = y - eh * (eh @ y)
        z = A @ (g * (AT @ y))
        return z - eh * (eh @ z)

    n = geom.n
    op = LinearOperator((n, n), matvec=mv, dtype=float)
    b = geom.b[:, col] - eh * (eh @ geom.b[:, col])
    x, info1 = cg(op, b, rtol=rtol, atol=0.0, maxiter=50 * n)
    u, info2 = cg(op, x, rtol=rtol, atol=0.0, maxiter=50 * n)
    K = float(b @ x)
    N = float(x @ x)
    res = float(np.linalg.norm(mv(x) - b) / np.linalg.norm(b))
    return {"K": K, "N": N, "tau": N / K, "cg_info": [int(info1), int(info2)], "rel_residual": res}


def arbcert(geom, g, col, prec=128):
    from certify import evaluate, summary
    out = summary(evaluate(geom, g, col=col, prec=prec))
    return {k: v for k, v in out.items()}


def verify(geom, g, col=0, r_ref=None, K0=None, nmax_arb=400, prec=128):
    t0 = time.time()
    rec = {"n": geom.n, "m": geom.m, "col": col}
    rec["min_g"] = float(np.min(g))
    rec["chol"] = chol(geom, g, col)
    rec["eig"] = eig(geom, g, col)
    rec["cg"] = cgsolve(geom, g, col)
    if r_ref is not None:
        r = geom.W @ g
        rec["diag_resid_rel_max"] = float(np.max(np.abs(r - r_ref) / r_ref))
    if K0 is not None:
        rec["K_resid_rel"] = rec["chol"]["K"] / K0 - 1
    t = rec["chol"]["tau"]
    rec["tau_spread_rel"] = max(abs(rec["eig"]["tau"] / t - 1), abs(rec["cg"]["tau"] / t - 1))
    if geom.n <= nmax_arb:
        try:
            rec["arb"] = arbcert(geom, g, col, prec)
            rec["arb_vs_chol_rel"] = abs(rec["arb"]["tau"]["mid"] / t - 1)
        except Exception as exc:  # noqa: BLE001
            rec["arb"] = {"error": repr(exc)}
    rec["seconds"] = time.time() - t0
    return rec


if __name__ == "__main__":
    print("library module")
