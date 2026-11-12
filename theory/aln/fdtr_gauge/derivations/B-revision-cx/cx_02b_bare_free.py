"""cx_02b: Table 0.2, free c_n rows (bare film, read-out T(0)), 50 digits, random points; plus bare film with the
theta_e read-out (eps0 = 1) on the claimed strata of Proposition C4 (a)-(d)."""
import sys as _s
_s.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle as _t
_t.wait_slot(2)          # light wait before any numerical import
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
OM = [mp.mpf(10) ** (mp.mpf(-3) + 6 * mp.mpf(i) / 39) for i in range(40)]
rng = np.random.default_rng(99)
U = lambda a, b: mp.mpf(float(rng.uniform(a, b)))
bn1 = lambda: U(0.3, 0.8) if rng.random() < 0.5 else U(1.25, 3.0)
ks = lambda b: b / (b * b - 1)
recs = []


def run(spec, th, pred, kers, tag):
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
    ok = rank == pred and all(r < 1e-35 for r in kr)
    recs.append(dict(tag=tag, spec=spec.label(), theta=[float(x) for x in th], rank=rank, pred=pred, kres=kr, ok=ok,
                     smin_retained=float(sv[rank - 1])))
    print("OK  " if ok else "FAIL", tag, spec.label(), "rank", rank, "pred", pred, "retained>=%.2e" % float(sv[rank - 1]),
          ["%.0e" % r for r in kr], flush=True)


for i in range(3):
    t = U(0.05, 1)
    # free c_n, R known (0 and >0) and unknown
    for Rk in (mp.mpf(0), mp.mpf("0.37"), None):
        sp = cx.Spec("bare", "free", 0.0, 0.0, R_known=Rk)
        nR = 1 if Rk is None else 0
        Rv = [U(.05, 1)] if Rk is None else []
        run(sp, [U(.3, 2), bn1(), U(.05, 1), U(.05, 1), U(-.5, 1)] + Rv, 5 + nR, [], "free off-Delta")
        run(sp, [U(.3, 2), bn1(), t, t, U(.05, 1)] + Rv, 4 + nR, [[0, 0, 1, 1, 0] + [0] * nR], "free Delta tn!=0")
        b = bn1()
        if Rk is None:
            run(sp, [U(.3, 2), b, t, t, 0, U(.05, 1)], 5, [[0, 0, 1, 1, 0, 0]], "free Delta tn=0 R>0 unknown")
            run(sp, [U(.3, 2), b, t, t, 0, 0], 4, [[0, 0, 1, 1, 0, 0], [ks(b), 0, 0, 0, 0, 1]], "free Delta tn=0 R=0 unknown")
            run(sp, [U(.3, 2), 1, t, t, 0, 0], 4, [[0, 0, 1, 1, 0, 0], [1, 0, 0, 0, 0, 0]], "free D2 tn=0 R=0 unknown")
        elif Rk == 0:
            run(sp, [U(.3, 2), b, t, t, 0], 4, [[0, 0, 1, 1, 0]], "free Delta tn=0 R=0 known")
            run(sp, [U(.3, 2), 1, t, t, 0], 3, [[0, 0, 1, 1, 0], [1, 0, 0, 0, 0]], "free D2 tn=0 R=0 known")
        else:
            run(sp, [U(.3, 2), 1, t, t, 0], 4, [[0, 0, 1, 1, 0]], "free D2 tn=0 R>0 known")
    # theta_e read-out (eps0 = 1), Proposition C4
    for e in (mp.mpf(0), mp.mpf("0.5"), mp.mpf(1), mp.mpf(-1), mp.mpf("2.5")):
        sp = cx.Spec("bare", "tied", 1, e, R_known=None)
        run(sp, [U(.3, 2), bn1(), U(.05, 1), U(.05, 1), U(.05, 1)], 5, [], "theta_e off-Delta")
        run(sp, [U(.3, 2), bn1(), t, t, U(.05, 1)], 5, [], "theta_e Delta R>0 (off exotic)")
        b = bn1()
        if e == 0:
            run(sp, [U(.3, 2), b, t, t, 0], 4, [[ks(b), 0, 0, 0, 1]], "theta_e eps=0 Delta R=0")
            xi = U(.4, 1.5); bb = U(.2, .9); tt = 4 * xi * xi * (1 - bb * bb); RR = 2 * xi * (1 - bb * bb) / bb
            run(sp, [xi, bb, tt, tt, RR], 4, [[1, 0, -tt / xi, tt / xi, (bb * bb - 1) / bb]], "theta_e E_0")
        else:
            run(sp, [U(.3, 2), b, t, t, 0], 5, [], "theta_e Delta R=0 (eps!=0)")
            if e < 0:
                xi = U(.4, 1.5); bb = mp.sqrt(1 - e); tt = 4 * abs(e) * xi * xi / (1 + abs(e)); RR = 2 * abs(e) * xi / bb
                run(sp, [xi, bb, tt, tt, RR], 4, [[1, 0, -tt / xi, tt / xi, -e / bb]], "theta_e E_eps")
        run(sp, [U(.3, 2), b, 0, 0, 0], 4, [[ks(b), 0, 0, 0, 1]], "theta_e Pi R=0")
        # R known = 0 on D2: claimed degenerate (kernel e_xi) iff eps = 0 or t = 0
        spk = cx.Spec("bare", "tied", 1, e, R_known=mp.mpf(0))
        run(spk, [U(.3, 2), 1, t, t], 3 if e == 0 else 4, [[1, 0, 0, 0]] if e == 0 else [], "theta_e D2 R=0 known")
        # audit correction D-1: at eps = 1 the kernel is span{e_xi, (0,0,-1,1)} (rank 2)
        run(spk, [U(.3, 2), 1, 0, 0], 2 if e == 1 else 3, [[1, 0, 0, 0], [0, 0, -1, 1]] if e == 1 else [[1, 0, 0, 0]],
            "theta_e Pi b=bs R=0 known")

nf = sum(1 for r in recs if not r["ok"])
print("points", len(recs), "failures", nf)
json.dump(dict(npoints=len(recs), nfail=nf, records=recs),
          open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_02b_bare_free.json", "w"), indent=1)
