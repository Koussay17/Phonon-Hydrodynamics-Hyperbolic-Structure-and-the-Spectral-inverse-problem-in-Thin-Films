"""Rigorous evaluation of K, N, tau and the diagonal residual for a witness g (ball arithmetic).

The witness g (float64, hence an exact dyadic rational) and the model data (energies, kT,
velocities, stoichiometric coefficients) are converted exactly to arb balls; the entropy scales
d = sqrt(n(1+n)), n = 1/expm1(eps/kT), are evaluated with rigorous enclosures; C(g) is assembled
in arb_mat and the linear system (C + gamma e_hat e_hat^T) x = b is solved with arb_mat.solve,
which returns rigorous enclosures (or raises if the matrix cannot be certified nonsingular).
Since b is orthogonal to e (enclosed), x = C^+ b exactly for admissible g; e^T x is enclosed and
reported. Output: midpoints and radii of K, N, tau and of max relative |W g - r_ref|.
"""
from __future__ import annotations

import numpy as np
from flint import arb, arb_mat, ctx


def arb_vec(v):
    return [arb(float(x)) for x in v]


def model_balls(geom, prec=200):
    ctx.prec = prec
    kT = arb(float(geom.kT))
    eps = [arb(float(e)) for e in geom.eps]
    dsc = []
    for e in eps:
        nocc = 1 / (arb.expm1(e / kT))
        dsc.append((nocc * (1 + nocc)).sqrt())
    return eps, dsc


def assemble(geom, g, eps, dsc, col=0):
    n, m = geom.n, geom.m
    rows = [[arb(0) for _ in range(n)] for _ in range(n)]
    diag = [arb(0) for _ in range(n)]
    for a in range(m):
        ga = arb(float(g[a]))
        if float(g[a]) == 0.0:
            continue
        ids, vals = [], []
        for s in range(3):
            mu = int(geom.events[a, s])
            if mu < 0:
                continue
            ids.append(mu)
            vals.append(arb(float(geom.coef[a, s])) / dsc[mu])
        for i, vi in zip(ids, vals):
            for j, vj in zip(ids, vals):
                rows[i][j] += ga * vi * vj
    C = arb_mat(rows)
    e = [dsc[i] * eps[i] for i in range(n)]
    b = [dsc[i] * eps[i] * arb(float(geom.vel[i, col])) for i in range(n)]
    return C, e, b


def evaluate(geom, g, col=0, prec=200, r_ref=None):
    eps, dsc = model_balls(geom, prec)
    C, e, b = assemble(geom, g, eps, dsc, col)
    n = geom.n
    ee = sum((x * x for x in e), arb(0))
    enorm = ee.sqrt()
    eh = [x / enorm for x in e]
    gam = sum((C[i, i] for i in range(n)), arb(0)) / n
    Creg = arb_mat([[C[i, j] + gam * eh[i] * eh[j] for j in range(n)] for i in range(n)])
    x = Creg.solve(arb_mat([[bi] for bi in b]))
    xs = [x[i, 0] for i in range(n)]
    K = sum((b[i] * xs[i] for i in range(n)), arb(0))
    N = sum((xs[i] * xs[i] for i in range(n)), arb(0))
    tau = N / K
    be = sum((b[i] * e[i] for i in range(n)), arb(0))
    ex = sum((e[i] * xs[i] for i in range(n)), arb(0))
    out = {"K": K, "N": N, "tau": tau, "b_dot_e": be, "e_dot_x": ex}
    if r_ref is not None:
        diag = [C[i, i] for i in range(n)]
        res = max(abs(float((diag[i] - arb(float(r_ref[i]))).mid())) / float(r_ref[i]) for i in range(n))
        rad = max(float(diag[i].rad()) / float(r_ref[i]) for i in range(n))
        out["diag_resid_rel"] = res
        out["diag_rad_rel"] = rad
    return out


def summary(out):
    s = {}
    for k, v in out.items():
        if hasattr(v, "mid"):
            s[k] = {"mid": float(v.mid()), "rad": float(v.rad())}
        else:
            s[k] = v
    return s
