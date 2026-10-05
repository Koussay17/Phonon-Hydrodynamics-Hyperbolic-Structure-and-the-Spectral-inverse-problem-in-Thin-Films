"""Exact Ben-Tal--Teboulle constant M = max_I |y_I| by batched enumeration of all (d-1)-subsets.

y_I solves a_alpha^T y = 0 (alpha in I), b^T y = 1 on H = e^perp (reduced coordinates).
Singular subsets are detected by the condition number of B_I (a clear gap separates exactly
singular subsets, cond ~ 1e15-1e17, from nonsingular ones); the gap is reported.
Second method: the maximising flat is recomputed from its rational kernel vector phi
(integer stoichiometry), |y_I| = |P_H D phi| / |b^T D phi|, in exact rational arithmetic for phi
and 50-digit arithmetic for the transcendental Bose scales.
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time

import numpy as np

from memgap import reduced_basis


def enumerate_M(geom, col=0, batch=200_000, cond_cut=1e11):
    Q = reduced_basis(geom)
    Ared = Q.T @ geom.A.toarray()
    bred = Q.T @ geom.b[:, col]
    d, m = Ared.shape
    k = d - 1
    total = math.comb(m, k)
    it = itertools.combinations(range(m), k)
    best, bestI = 0.0, None
    nbasis = 0
    cond_hist = np.zeros(40, dtype=np.int64)      # log10 cond histogram
    t0 = time.time()
    rhs = np.zeros(d); rhs[-1] = 1.0
    while True:
        chunk = list(itertools.islice(it, batch))
        if not chunk:
            break
        I = np.array(chunk, dtype=np.int64)               # (B, k)
        Bm = np.empty((len(I), d, d))
        Bm[:, :k, :] = np.transpose(Ared[:, I], (1, 2, 0))  # rows a_alpha^T
        Bm[:, k, :] = bred
        s = np.linalg.svd(Bm, compute_uv=False)
        cond = s[:, 0] / np.maximum(s[:, -1], 1e-300)
        lc = np.clip(np.floor(np.log10(cond)).astype(int), 0, 39)
        np.add.at(cond_hist, lc, 1)
        ok = cond < cond_cut
        if ok.any():
            y = np.linalg.solve(Bm[ok], np.broadcast_to(rhs, (ok.sum(), d))[..., None])[..., 0]
            nrm = np.linalg.norm(y, axis=1)
            nbasis += int(ok.sum())
            j = int(np.argmax(nrm))
            if nrm[j] > best:
                best = float(nrm[j])
                bestI = I[ok][j].tolist()
    return {"M": best, "I": bestI, "n_bases": nbasis, "n_subsets": total, "d": d, "m": m,
            "cond_hist_log10": cond_hist.tolist(), "seconds": time.time() - t0}


def verify_flat_exact(geom, I, col=0, dps=50):
    """Recompute |y_I| from the rational kernel of the integer stoichiometric matrix s_I."""
    import sympy
    import mpmath
    n = geom.n
    S = np.zeros((n, len(I)), dtype=object)
    for jj, a in enumerate(I):
        for slot in range(3):
            mu = geom.events[a, slot]
            if mu >= 0:
                c = geom.coef[a, slot]
                assert float(c).is_integer(), "exact verification needs integer stoichiometry"
                S[mu, jj] += int(round(c))
    Sm = sympy.Matrix(S.tolist())
    ker = (Sm.T).nullspace()
    assert len(ker) == 2, f"kernel dimension {len(ker)}"
    # energies are integers c_b |kappa| in units of 2 pi / N (1D exact model)
    N = geom.labels["N"]
    unit = 2 * math.pi / N
    eint = [sympy.Integer(int(round(e / unit))) for e in geom.eps]
    epsv = sympy.Matrix(eint)
    # choose phi in ker not proportional to eps
    K = sympy.Matrix.hstack(*ker)
    # component of kernel orthogonal (Euclidean) to eps is fine: any phi with phi not || eps
    cand = [v for v in ker]
    phi = None
    for v in cand:
        if sympy.Matrix.hstack(v, epsv).rank() == 2:
            phi = v
            break
    assert phi is not None
    mpmath.mp.dps = dps
    kT = mpmath.mpf(geom.kT)
    epsf = [mpmath.mpf(int(e)) * mpmath.mpf(unit) for e in eint]
    # use the same float unit as the geometry for consistency
    epsf = [mpmath.mpf(float(e)) for e in geom.eps]
    dsc = [mpmath.sqrt(1 / mpmath.expm1(e / kT) * (1 + 1 / mpmath.expm1(e / kT))) for e in epsf]
    Dphi = [dsc[i] * mpmath.mpf(sympy.Rational(phi[i]).p) / mpmath.mpf(sympy.Rational(phi[i]).q) for i in range(n)]
    e = [dsc[i] * epsf[i] for i in range(n)]
    vel = [mpmath.mpf(float(v)) for v in geom.vel[:, col]]
    bvec = [dsc[i] * epsf[i] * vel[i] for i in range(n)]
    ee = mpmath.fsum(x * x for x in e)
    pe = mpmath.fsum(Dphi[i] * e[i] for i in range(n))
    PH = [Dphi[i] - pe / ee * e[i] for i in range(n)]
    num = mpmath.sqrt(mpmath.fsum(x * x for x in PH))
    den = abs(mpmath.fsum(bvec[i] * Dphi[i] for i in range(n)))
    return {"phi": [str(sympy.Rational(x)) for x in phi], "norm_y": float(num / den), "norm_y_str": mpmath.nstr(num / den, 30)}


if __name__ == "__main__":
    from debye_events import build
    d = int(sys.argv[1]); N = int(sys.argv[2])
    geom = build(d, N)
    res = enumerate_M(geom)
    g0 = geom.gphys
    K0 = geom.response(g0)["K"]
    res["K0"] = K0
    res["BT_bound_tau"] = res["M"] ** 2 * K0
    print(json.dumps({k: v for k, v in res.items() if k != "cond_hist_log10"}, indent=1))
    print("cond histogram (log10 bins):", res["cond_hist_log10"])
    if geom.labels["exact"] and res["I"] is not None:
        ver = verify_flat_exact(geom, res["I"])
        res["exact_check"] = ver
        print("exact rational-kernel check: |y_I| =", ver["norm_y_str"], " float:", res["M"])
    with open(f"../results/bt_exact_d{d}_N{N}.json", "w") as f:
        json.dump(res, f, indent=1)
