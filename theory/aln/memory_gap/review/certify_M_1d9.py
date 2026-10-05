"""Independent certification of the Ben-Tal--Teboulle constant M for the exact 1D N = 9 Debye model
(07-experiments.md E1: M = 16284.0067 by float enumeration with a condition-number cut at 1e11).

Gap in the original argument: a subset I is declared singular when cond([a_I^T; b^T]) >= 1e11. An exactly
singular subset and a subset with a genuine but enormous |y_I| (>~ 1e15) are indistinguishable in float64, so
the float enumeration by itself only certifies the maximum over well-conditioned subsets.

Method here (P3 route; independent of the original reduced-coordinate solve):
  For each 14-subset I of the 26 events, S_I is a 16 x 14 INTEGER matrix (entries 0, +-1, +2).
  (i) Integer rank.  If rank S_I = 14, det(S_I^T S_I) is a positive integer, so prod_i sigma_i^2 >= 1 and
      sigma_14 >= L_I := 1 / prod_{i<=13} sigma_i.  A computed sigma_14 < 1e-9 with sigma_14 / L_I << 1
      certifies rank < 14 (then v_I = 0: the subset contributes nothing).  sigma_14 > 1e-9 certifies rank 14.
  (ii) Rank 14: ker S_I^T = span{eps, phi_I}, phi_I orthogonal to eps; v_I^T b != 0 <=> beta_I = b^T D phi_I != 0
      and |y_I| = |P_H D phi_I| / |beta_I|  (P3).
  (iii) Time reversal T (k -> -k within a branch) fixes eps and D and b is exactly odd.  If T phi_I = +phi_I then
      beta_I = 0 exactly.  Rigour of the float test: phi_I is parallel to the integer Cramer vector c_I of the
      15 x 16 integer matrix [S_I^T; eps^T]; by Hadamard |c_I| <= 4 * 5^7 * |eps| =: Hmax.  If T c != c then
      |T c - c| >= 1, so for the unit vector |T phi - phi| >= 1/Hmax (~1.8e-7).  The float test uses 1e-9.
  (iv) For every rank-14 subset that is not T-even, |y_I| is evaluated in float (well conditioned); the largest
      values are recomputed exactly (sympy integer kernel + 60-digit Bose scales).
Run: python -B certify_M_1d9.py   (writes certify_M_1d9.json)
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

import numpy as np  # noqa: E402
import mpmath  # noqa: E402
import sympy  # noqa: E402

from debye_events import build  # noqa: E402


def integer_stoich(geom):
    n, m = geom.n, geom.m
    S = np.zeros((n, m), dtype=np.int64)
    for a in range(m):
        for slot in range(3):
            mu = geom.events[a, slot]
            if mu >= 0:
                c = geom.coef[a, slot]
                assert float(c).is_integer()
                S[mu, a] += int(round(c))
    return S


def exact_ynorm(S, I, eint, T, geom, N):
    n = S.shape[0]
    SI = sympy.Matrix(S[:, list(I)].tolist())
    rows = SI.T.col_join(sympy.Matrix([[int(v) for v in eint]]))
    ker = rows.nullspace()
    assert len(ker) == 1
    phi = ker[0]
    den = sympy.ilcm(*[sympy.fraction(x)[1] for x in phi])
    phi = (phi * den).applyfunc(sympy.Integer)
    gg = sympy.igcd(*[int(x) for x in phi])
    phi = phi / gg
    is_even = all(phi[int(T[i])] == phi[i] for i in range(n))
    mpmath.mp.dps = 60
    kT = mpmath.mpf(geom.kT)
    unit = 2 * mpmath.pi / N
    eps_mp = [mpmath.mpf(int(v)) * unit for v in eint]
    dsc = [mpmath.sqrt(mpmath.exp(x / kT)) / (mpmath.exp(x / kT) - 1) for x in eps_mp]
    vel = geom.vel[:, 0]
    b = [dsc[i] * eps_mp[i] * int(vel[i]) for i in range(n)]
    e = [dsc[i] * eps_mp[i] for i in range(n)]
    Dphi = [dsc[i] * int(phi[i]) for i in range(n)]
    beta = mpmath.fsum(b[i] * Dphi[i] for i in range(n))
    pe = mpmath.fsum(Dphi[i] * e[i] for i in range(n)) / mpmath.fsum(x * x for x in e)
    PH = [Dphi[i] - pe * e[i] for i in range(n)]
    nPH = mpmath.sqrt(mpmath.fsum(x * x for x in PH))
    return {"phi": [str(x) for x in phi], "even": bool(is_even),
            "ynorm": mpmath.nstr(nPH / abs(beta), 30) if beta != 0 else None}


def main():
    t0 = time.time()
    geom = build(1, 9)
    n, m = geom.n, geom.m
    S = integer_stoich(geom)
    N = geom.labels["N"]
    unit = 2 * math.pi / N
    eint = np.rint(geom.eps / unit).astype(np.int64)
    assert np.allclose(eint * unit, geom.eps)
    assert np.all(S.T @ eint == 0), "integer energy conservation"
    mk, mb = geom.labels["modes_k"], geom.labels["modes_b"]
    key = {(int(b), tuple(k)): i for i, (b, k) in enumerate(zip(mb, map(tuple, mk)))}
    T = np.array([key[(int(b), tuple(-k))] for b, k in zip(mb, mk)], dtype=np.int64)
    b = geom.b[:, 0]
    Dd = geom.dscale
    e = geom.e
    assert np.array_equal(b[T], -b) and np.array_equal(Dd[T], Dd) and np.array_equal(eint[T], eint)
    Hmax = 4 * 5 ** 7 * float(np.linalg.norm(eint))
    even_tol = 1e-9
    assert even_tol < 0.01 / Hmax

    epsn = eint / np.linalg.norm(eint)
    k = n - 2
    total = math.comb(m, k)
    it = itertools.combinations(range(m), k)
    batch = 100_000
    st = {"n_subsets": total, "rank_deficient": 0, "rank14": 0, "rank14_T_even": 0, "rank14_not_even": 0,
          "max_sigma14_deficient": 0.0, "min_sigma14_rank14": np.inf,
          "max_sigma14_over_L_deficient": 0.0, "Hadamard_Hmax": Hmax, "even_threshold": even_tol,
          "max_dist_even_flagged": 0.0, "min_dist_even_not_flagged": np.inf,
          "max_rel_beta_even": 0.0, "min_rel_beta_not_even": np.inf}
    hist_logy = np.zeros(40, dtype=np.int64)
    top = []                                   # (ynorm, subset)
    while True:
        chunk = list(itertools.islice(it, batch))
        if not chunk:
            break
        I = np.array(chunk, dtype=np.int64)
        SI = np.transpose(S[:, I], (1, 0, 2)).astype(float)
        U, s, _ = np.linalg.svd(SI, full_matrices=True)
        s14 = s[:, -1]
        with np.errstate(divide="ignore"):
            L = 1.0 / np.prod(s[:, :-1], axis=1)
        defi = s14 < 1e-9
        st["rank_deficient"] += int(defi.sum())
        if defi.any():
            st["max_sigma14_deficient"] = max(st["max_sigma14_deficient"], float(s14[defi].max()))
            ratio = np.where(np.isfinite(L[defi]), s14[defi] / L[defi], 0.0)
            st["max_sigma14_over_L_deficient"] = max(st["max_sigma14_over_L_deficient"], float(ratio.max()))
        ok = ~defi
        if not ok.any():
            continue
        st["rank14"] += int(ok.sum())
        st["min_sigma14_rank14"] = min(st["min_sigma14_rank14"], float(s14[ok].min()))
        N2 = U[ok][:, :, k:]
        c = np.einsum("bij,i->bj", N2, epsn)
        phi = N2[:, :, 0] * c[:, 1:2] - N2[:, :, 1] * c[:, 0:1]
        phi /= np.linalg.norm(phi, axis=1)[:, None]
        dev = np.linalg.norm(phi[:, T] - phi, axis=1)
        even = dev < even_tol
        if even.any():
            st["max_dist_even_flagged"] = max(st["max_dist_even_flagged"], float(dev[even].max()))
        if (~even).any():
            st["min_dist_even_not_flagged"] = min(st["min_dist_even_not_flagged"], float(dev[~even].min()))
        st["rank14_T_even"] += int(even.sum())
        st["rank14_not_even"] += int((~even).sum())
        Dphi = phi * Dd[None, :]
        beta = Dphi @ b
        relb = np.abs(beta) / (np.linalg.norm(b) * np.linalg.norm(Dphi, axis=1))
        if even.any():
            st["max_rel_beta_even"] = max(st["max_rel_beta_even"], float(relb[even].max()))
        ne = ~even
        if ne.any():
            st["min_rel_beta_not_even"] = min(st["min_rel_beta_not_even"], float(relb[ne].min()))
            PH = Dphi[ne] - np.outer(Dphi[ne] @ e / (e @ e), e)
            yn = np.linalg.norm(PH, axis=1) / np.abs(beta[ne])
            np.add.at(hist_logy, np.clip(np.floor(np.log10(yn)).astype(int) + 5, 0, 39), 1)
            Ine = I[ok][ne]
            kk = min(20, len(yn))
            idx = np.argpartition(-yn, kk - 1)[:kk]
            for j in idx:
                top.append((float(yn[j]), Ine[j].tolist()))
            top.sort(key=lambda t: -t[0])
            top = top[:40]
    st["seconds_enumeration"] = time.time() - t0
    st["hist_log10_ynorm_offset5"] = hist_logy.tolist()
    # exact recomputation of the distinct top values
    seen, exact = set(), []
    for yv, Isub in top:
        kf = round(yv, 6)
        if kf in seen:
            continue
        seen.add(kf)
        ex = exact_ynorm(S, Isub, eint, T, geom, N)
        exact.append({"float": yv, "subset": Isub, **ex})
        if len(exact) >= 6:
            break
    out = {"stats": st, "top_exact": exact,
           "M_certified": exact[0]["ynorm"], "M_campaign_float": 16284.006724417613,
           "M_campaign_rational_kernel": "16284.0067235629307297011295459",
           "n_bases_campaign": 3566333, "seconds_total": time.time() - t0}
    (HERE / "certify_M_1d9.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
