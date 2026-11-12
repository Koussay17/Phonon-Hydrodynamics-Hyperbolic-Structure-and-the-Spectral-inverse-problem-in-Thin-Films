"""cx_02: 50-digit function-level ranks on random points of every stratum of Tables 0.2 / 0.3 (and more).

Jacobian of log Z w.r.t. the parameters at 40 frequencies, log-spaced in [1e-3, 1e3]; columns normalised.
Numerical rank = number of singular values > 1e-30 (the deficient ones sit at ~1e-48 at 50 digits).
Each stratum: several random points drawn independently of the revision branch's test points.
Predicted kernels (from the claims) are applied to the unnormalised Jacobian.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
OM = [mp.mpf(10) ** (mp.mpf(-3) + 6 * mp.mpf(i) / 39) for i in range(40)]
rng = np.random.default_rng(20261112)
TH = mp.mpf("1e-30")


def U(a, b):
    return float(rng.uniform(a, b))


def bnot1():
    # b away from b_s = 1
    return U(0.3, 0.8) if rng.random() < 0.5 else U(1.25, 3.0)


def ks(b):
    b = mp.mpf(b)
    return b / (b * b - 1)          # b_s = 1: k_s = b b_s^2/(b^2 - b_s^2)


def run(spec, th, pred_rank, pred_kernels, tag):
    J = cx.realjac_mp(spec, th, OM, rel=True)
    # exactly-zero columns (|col| < 1e-35 max|col|) are not normalised: they are kernel directions
    norms = [mp.sqrt(mp.fsum(J[i, k] ** 2 for i in range(J.rows))) for k in range(J.cols)]
    big = max(norms)
    zc = [k for k in range(J.cols) if norms[k] < mp.mpf("1e-35") * big]
    Jn = J.copy()
    for k in range(J.cols):
        for i in range(J.rows):
            Jn[i, k] = J[i, k] / norms[k] if k not in zc else mp.mpf(0)
    S = mp.svd_r(Jn, compute_uv=False)
    sv = sorted([S[i] for i in range(len(S))], reverse=True)
    rank = sum(1 for x in sv if x > TH)
    kres = [float(cx.kernel_residual_mp(J, c)) for c in pred_kernels]
    ok = (rank == pred_rank) and all(r < 1e-35 for r in kres)
    rec = dict(tag=tag, spec=spec.label(), theta=[float(x) for x in th], rank=rank, pred=pred_rank,
               sv_min_retained=float(sv[rank - 1]) if rank > 0 else None,
               sv_max_deficient=float(sv[rank]) if rank < len(sv) else None,
               kernel_residuals=kres, ok=ok)
    print(("OK  " if ok else "FAIL"), tag, spec.label(), "rank", rank, "pred", pred_rank,
          "retained>=%.2e" % rec["sv_min_retained"],
          "deficient<=%.2e" % (rec["sv_max_deficient"] or 0), "kres", ["%.1e" % r for r in kres], flush=True)
    return rec


records = []
NP = 3   # points per stratum

# ---------------------------------------------------------------- stack, local
sp = cx.Spec("stack", "tied", 0.0, 0.0)
for i in range(NP):
    tR, tl = U(0.05, 1), U(0.05, 1)
    records.append(run(sp, [U(.3, 2), bnot1(), tR, tl, U(.05, 1), U(.05, 1)], 6, [], "L off-Delta"))
    records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1), U(.05, 1), 0], 6, [], "L off-Delta R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, U(.05, 1), U(.05, 1), U(.05, 1), 0], 6, [], "L off-Delta b=bs R2=0"))
    records.append(run(sp, [U(.3, 2), bnot1(), 0.0, U(.05, 1), 0.0, U(.05, 1)], 6, [], "L tR=0 R1=0"))
    records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), 0.0, U(.05, 1), 0.0], 6, [], "L tl=0 R2=0"))
    t = U(0.05, 1)
    records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), U(.05, 1)], 5, [[0, 0, 1, 1, 0, 0]], "L Delta R2>0"))
    b = bnot1()
    records.append(run(sp, [U(.3, 2), b, t, t, U(.05, 1), 0], 4, [[0, 0, 1, 1, 0, 0], [ks(b), 0, 0, 0, 0, 1]],
                       "L Delta R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, t, t, U(.05, 1), 0], 4, [[0, 0, 1, 1, 0, 0], [1, 0, 0, 0, 0, 0]],
                       "L D2 R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, t, t, U(.05, 1), U(.05, 1)], 5, [[0, 0, 1, 1, 0, 0]], "L D2 R2>0"))
    b = bnot1()
    records.append(run(sp, [U(.3, 2), b, 0, 0, U(.05, 1), 0], 4, [[0, 0, 1, 1, 0, 0], [ks(b), 0, 0, 0, 0, 1]],
                       "L Pi R2=0"))

# ---------------------------------------------------------------- stack, tied eps != 0 (same at both faces)
for eps in [1.0, 0.5, 0.084, -0.5, 2.5, -1.0, 1 / 3, 2.0, -3.0, 1e-3]:
    sp = cx.Spec("stack", "tied", eps, eps)
    for i in range(2):
        t = U(0.05, 1)
        records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1), U(.05, 1), U(.05, 1)], 6, [], "T interior"))
        records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), U(.05, 1)], 6, [], "T Delta R2>0"))
        records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), 0], 6, [], "T Delta R2=0"))
        records.append(run(sp, [U(.3, 2), 1.0, t, t, U(.05, 1), 0], 6, [], "T D2 R2=0"))
        records.append(run(sp, [U(.3, 2), bnot1(), 0, 0, U(.05, 1), U(.05, 1)], 6, [], "T Pi R2>0"))
        records.append(run(sp, [U(.3, 2), 1.0, 0, 0, U(.05, 1), U(.05, 1)], 6, [], "T Pi b=bs R2>0"))
        b = bnot1()
        records.append(run(sp, [U(.3, 2), b, 0, 0, U(.05, 1), 0], 5, [[ks(b), 0, 0, 0, 0, 1]], "T Pi R2=0"))
        records.append(run(sp, [U(.3, 2), 1.0, 0, 0, U(.05, 1), 0], 4, [[1, 0, 0, 0, 0, 0], [0, 0, 1 - 2 * mp.mpf(eps), 1, 0, 0]],
                           "T Pi b=bs R2=0"))
        records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), 0, U(.05, 1), 0], 6, [], "T tl=0 R2=0"))
        records.append(run(sp, [U(.3, 2), 1.0, U(.05, 1), 0, U(.05, 1), 0], 6, [], "T tl=0 b=bs R2=0"))
        records.append(run(sp, [U(.3, 2), bnot1(), 0, U(.05, 1), U(.05, 1), 0], 6, [], "T tR=0 R2=0"))
        tl = U(0.05, 1)
        if eps > 0:
            records.append(run(sp, [U(.3, 2), bnot1(), eps * tl, tl, U(.05, 1), 0], 6, [], "T tR=eps*tl R2=0"))

# ---------------------------------------------------------------- stack, free c_n
sp = cx.Spec("stack", "free")
for i in range(NP):
    t = U(0.05, 1)
    records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1), U(.05, 1), U(.05, 1), U(.05, 1)], 7, [], "F off-Delta"))
    records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1), U(-.5, 1), U(.05, 1), 0], 7, [],
                       "F off-Delta R2=0 (tn any sign)"))
    records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), U(.05, 1), U(.05, 1)], 6, [[0, 0, 1, 1, 0, 0, 0]], "F Delta"))
    records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), U(.05, 1), 0], 6, [[0, 0, 1, 1, 0, 0, 0]], "F Delta R2=0"))
    records.append(run(sp, [U(.3, 2), bnot1(), t, t, 0, U(.05, 1), U(.05, 1)], 6, [[0, 0, 1, 1, 0, 0, 0]], "F Delta tn=0"))
    b = bnot1()
    records.append(run(sp, [U(.3, 2), b, t, t, 0, U(.05, 1), 0], 5, [[0, 0, 1, 1, 0, 0, 0], [ks(b), 0, 0, 0, 0, 0, 1]],
                       "F Delta tn=0 R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, t, t, 0, U(.05, 1), 0], 5, [[0, 0, 1, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0]],
                       "F D2 tn=0 R2=0"))
    records.append(run(sp, [U(.3, 2), b, 0, 0, 0, U(.05, 1), 0], 5, [[0, 0, 1, 1, 0, 0, 0], [ks(b), 0, 0, 0, 0, 0, 1]],
                       "F Pi tn=0 R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, 0, 0, 0, U(.05, 1), 0], 4,
                       [[1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0], [0, 0, -2, 0, 1, 0, 0]], "F Pi tn=0 R2=0 b=bs"))
    tR, tl = U(0.05, 1), U(0.05, 1)
    records.append(run(sp, [U(.3, 2), bnot1(), tR, tl, tR, U(.05, 1), 0], 7, [], "F tn=tR R2=0"))
    records.append(run(sp, [U(.3, 2), bnot1(), tR, tl, tl, U(.05, 1), 0], 7, [], "F tn=tl R2=0"))
    records.append(run(sp, [U(.3, 2), bnot1(), 0, tl, tl, U(.05, 1), 0], 7, [], "F tn=tl tR=0 R2=0"))
    records.append(run(sp, [U(.3, 2), 1.0, 0, tl, 0, U(.05, 1), 0], 7, [], "F tR=tn=0 b=bs R2=0"))

# ---------------------------------------------------------------- stack, mixed face coefficients
for (e0, e) in [(0.3, 0.7), (-0.4, 1.3), (2.0, -1.0), (0.5, 0.0), (0.0, 0.5), (0.5, 1.0)]:
    sp = cx.Spec("stack", "tied", e0, e)
    for i in range(2):
        t = U(0.05, 1)
        records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1), U(.05, 1), U(.05, 1)], 6, [], "M interior"))
        records.append(run(sp, [U(.3, 2), bnot1(), t, t, U(.05, 1), U(.05, 1)], 6, [], "M Delta R2>0"))
        pr = 5 if e == 0 else 6
        b = bnot1()
        records.append(run(sp, [U(.3, 2), b, t, t, U(.05, 1), 0], pr, [[ks(b), 0, 0, 0, 0, 1]] if e == 0 else [],
                           "M Delta R2=0"))
        records.append(run(sp, [U(.3, 2), b, 0, 0, U(.05, 1), 0], 5, [[ks(b), 0, 0, 0, 0, 1]], "M Pi R2=0"))
        records.append(run(sp, [U(.3, 2), 1.0, 0, 0, U(.05, 1), 0], 5, [[1, 0, 0, 0, 0, 0]], "M Pi b=bs R2=0"))
        records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), 0, U(.05, 1), 0], 6, [], "M tl=0 R2=0"))

# exotic strata E_0 (entropic top, local bottom) and E_eps (entropic top, eps<0 bottom), stack
for i in range(2):
    xi = mp.mpf(U(0.4, 1.5))
    b = mp.mpf(U(0.2, 0.9))
    t = 4 * xi * xi * (1 - b * b)
    R = 2 * xi * (1 - b * b) / b
    sp = cx.Spec("stack", "tied", 1.0, 0.0)
    records.append(run(sp, [xi, b, t, t, U(.05, 1), R], 5, [[1, 0, -t / xi, t / xi, 0, (b * b - 1) / b]], "E0 stack"))
    e = -mp.mpf(U(0.2, 2.0))
    b = mp.sqrt(1 - e)
    t = 4 * abs(e) * xi * xi / (1 + abs(e))
    R = 2 * abs(e) * xi / b
    sp = cx.Spec("stack", "tied", 1.0, e)  # e is an mpf
    records.append(run(sp, [xi, b, t, t, U(.05, 1), R], 5, [[1, 0, -t / xi, t / xi, 0, -e / b]], "E_eps stack"))

# ---------------------------------------------------------------- bare film, read-out T(0)
for eps in [0.0, 1.0, 0.084, -0.5, 2.5, 0.5]:
    for Rk in [None, 0.0, 0.37]:
        sp = cx.Spec("bare", "tied", 0.0, eps, R_known=Rk)
        nR = 1 if Rk is None else 0
        for i in range(2):
            t = U(0.05, 1)
            Rv = [U(.05, 1)] if Rk is None else []
            R0 = [0.0] if Rk is None else []
            full = 4 + nR
            # interior
            records.append(run(sp, [U(.3, 2), bnot1(), U(.05, 1), U(.05, 1)] + Rv, full, [], "bare interior"))
            # on Delta with R>0 (unknown) / known
            if eps == 0:
                records.append(run(sp, [U(.3, 2), bnot1(), t, t] + Rv, full - 1, [[0, 0, 1, 1] + [0] * nR],
                                   "bare L Delta R>0/known"))
            else:
                records.append(run(sp, [U(.3, 2), bnot1(), t, t] + Rv, full, [], "bare T Delta R>0/known"))
            # unknown R at R = 0
            if Rk is None:
                b = bnot1()
                if eps == 0:
                    records.append(run(sp, [U(.3, 2), b, t, t, 0.0], 3, [[0, 0, 1, 1, 0], [ks(b), 0, 0, 0, 1]],
                                       "bare L Delta R=0 unknown"))
                    records.append(run(sp, [U(.3, 2), 1.0, t, t, 0.0], 3, [[0, 0, 1, 1, 0], [1, 0, 0, 0, 0]],
                                       "bare L D2 R=0 unknown"))
                else:
                    records.append(run(sp, [U(.3, 2), b, t, t, 0.0], 5, [], "bare T Delta R=0 unknown"))
                    records.append(run(sp, [U(.3, 2), b, 0, 0, 0.0], 4, [[ks(b), 0, 0, 0, 1]], "bare T Pi R=0 unknown"))
                    records.append(run(sp, [U(.3, 2), 1.0, 0, 0, 0.0], 4, [[1, 0, 0, 0, 0]], "bare T Pi b=bs R=0 unknown"))
            if Rk == 0.0:
                if eps == 0:
                    records.append(run(sp, [U(.3, 2), 1.0, t, t], 2, [[0, 0, 1, 1], [1, 0, 0, 0]], "bare L D2 R=0 known"))
                else:
                    records.append(run(sp, [U(.3, 2), 1.0, t, t], 4, [], "bare T D2 R=0 known"))
                    records.append(run(sp, [U(.3, 2), bnot1(), 0, 0], 4, [], "bare T Pi R=0 known"))
                    records.append(run(sp, [U(.3, 2), 1.0, 0, 0], 3, [[1, 0, 0, 0]], "bare T Pi b=bs R=0 known"))

nfail = sum(1 for r in records if not r["ok"])
print("points", len(records), "failures", nfail)
with open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_02_rank_tables.json", "w") as f:
    json.dump(dict(npoints=len(records), nfail=nfail, records=records), f, indent=1)
