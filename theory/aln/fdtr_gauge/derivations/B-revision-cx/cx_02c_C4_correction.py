"""cx_02c: check of the audit correction D-1 (Proposition C4, read-out theta_e = (1 + tau_l p) T(0), bottom face eps tau_l)
on Pi = {tau_R = tau_l = 0} with b = b_s and R = 0: claimed kernel span{e_xi, (0,0,-1,1[,0])} for eps = 1 (rank 2 of 4 with R
known, 3 of 5 with R unknown) and e_xi alone for eps != 1 (rank 3 / 4).  50 digits, 40 frequencies, random xi; also the
neighbouring strata D2 (t > 0) for eps = 1, where full rank is claimed."""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
OM = [mp.mpf(10) ** (mp.mpf(-3) + 6 * mp.mpf(i) / 39) for i in range(40)]
rng = np.random.default_rng(1212)
recs = []


def rank_and_kernels(spec, th, kers):
    J = cx.realjac_mp(spec, th, OM)
    norms = [mp.sqrt(mp.fsum(J[i, k] ** 2 for i in range(J.rows))) for k in range(J.cols)]
    big = max(norms)
    Jn = J.copy()
    for k in range(J.cols):
        for i in range(J.rows):
            Jn[i, k] = J[i, k] / norms[k] if norms[k] > mp.mpf("1e-35") * big else mp.mpf(0)
    S = mp.svd_r(Jn, compute_uv=False)
    sv = sorted([S[i] for i in range(len(S))], reverse=True)
    rank = sum(1 for x in sv if x > mp.mpf("1e-30"))
    kr = [float(cx.kernel_residual_mp(J, c)) for c in kers]
    return rank, [float(x) for x in sv], kr


for eps in ["1", "0.5", "0", "-1", "2"]:
    e = mp.mpf(eps)
    for Rk in ("known", "unknown"):
        for trial in range(2):
            xi = mp.mpf(float(rng.uniform(0.4, 1.6)))
            if Rk == "known":
                sp = cx.Spec("bare", "tied", 1, e, R_known=mp.mpf(0))
                th = [xi, 1, 0, 0]
                kers = [[1, 0, 0, 0], [0, 0, -1, 1]]
                pred = 2 if e == 1 else 3
            else:
                sp = cx.Spec("bare", "tied", 1, e, R_known=None)
                th = [xi, 1, 0, 0, 0]
                kers = [[1, 0, 0, 0, 0], [0, 0, -1, 1, 0]]
                pred = 3 if e == 1 else 4
            rank, sv, kr = rank_and_kernels(sp, th, kers)
            ok = rank == pred and kr[0] < 1e-35 and ((kr[1] < 1e-35) if e == 1 else (kr[1] > 1e-20))
            recs.append(dict(eps=eps, R=Rk, xi=float(xi), rank=rank, pred=pred, kernel_residuals=kr, ok=ok,
                             sv_min_retained=sv[rank - 1]))
            print(("OK  " if ok else "FAIL"), f"eps={eps} R {Rk}: rank {rank} (claimed after D-1: {pred}); residual of e_xi {kr[0]:.1e}, "
                  f"of (0,0,-1,1) {kr[1]:.1e}; smallest retained s {sv[rank - 1]:.2e}", flush=True)
# D2 with t > 0, eps = 1, R = 0 known: claimed full rank 4
for trial in range(3):
    xi = mp.mpf(float(rng.uniform(0.4, 1.6))); t = mp.mpf(float(rng.uniform(0.05, 1)))
    sp = cx.Spec("bare", "tied", 1, 1, R_known=mp.mpf(0))
    rank, sv, kr = rank_and_kernels(sp, [xi, 1, t, t], [])
    recs.append(dict(eps="1", R="known", stratum="D2 t>0", rank=rank, pred=4, ok=rank == 4, sv_min_retained=sv[rank - 1]))
    print(("OK  " if rank == 4 else "FAIL"), f"eps=1 D2 t={float(t):.3f} R=0 known: rank {rank} (claimed 4), smallest s {sv[rank - 1]:.2e}",
          flush=True)
json.dump(recs, open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_02c_C4_correction.json", "w"),
          indent=1)
