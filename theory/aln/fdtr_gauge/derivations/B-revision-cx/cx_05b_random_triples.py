"""cx_05b: the branch's own targets with random frequency triples (entropic convention): multistart for
three-frequency twins.  Shows whether the absence of entropic twins at (1/2, 4, 1.7) and (0.2, 1.1, 6) is a
property of the target or of the two triples.
Usage: python cx_05b_random_triples.py <target: generic|bbs|alnlike> <ntriples> <starts>
"""
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
import cxtwin as ct

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
tname = sys.argv[1]
NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NS = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
conv = sys.argv[4] if len(sys.argv) > 4 else "entropic"
eps = {"entropic": 1.0, "local": 0.0, "kin": 0.084}[conv]
tgt = {"generic": [0.8, 2.0, 0.5, 0.1, 0.3, 0.2], "bbs": [0.8, 1.0, 0.5, 0.1, 0.3, 0.2],
       "alnlike": [1.0, 2.07, 0.03, 0.012, 0.05, 0.1]}[tname]
tgt = np.array(tgt)
sp = cx.Spec("stack", "tied", eps, eps, xiA=0.2, bA=3.0)
rng = np.random.default_rng(55)
lt = np.log(tgt)
out = dict(target=tgt.tolist(), conv=conv, triples=[])
for t in range(NTR):
    om = np.sort(10 ** rng.uniform(-1.3, 1.3, 3))
    om = np.round(om, 4)
    Zt = cx.Z_np(sp, tgt, om)
    Y0 = lt[None, :] + rng.uniform(-3, 3, (NS, 6))
    Y, res = ct.batch_lm(sp, Y0, om, Zt, iters=250)
    roots = []
    for y, r in zip(Y, res):
        if r < 1e-11 and np.max(np.abs(y - lt)) > 1e-5 and not any(np.max(np.abs(y - q)) < 1e-6 for q in roots):
            roots.append(y)
    twins = []
    mp.mp.dps = 60
    for y in roots:
        thm, nF = ct.refine_mp(sp, np.exp(y), [repr(float(w)) for w in om], [mp.mpf(repr(float(x))) for x in tgt], dps=60)
        if nF < mp.mpf("1e-45") and min(thm) > 0:
            Ztm = [cx.Z_mp(sp, [mp.mpf(repr(float(x))) for x in tgt], mp.mpf(repr(float(w)))) for w in om]
            sv = ct.jac_cond_mp(sp, thm, [repr(float(w)) for w in om], Ztm)
            twins.append(dict(theta=[mp.nstr(x, 30) for x in thm], residual_60=float(nF), smin_over_smax=float(sv[-1] / sv[0])))
    out["triples"].append(dict(omega=om.tolist(), n_twins=len(twins), twins=twins))
    print(f"triple {om.tolist()}: twins {len(twins)}", [tw["theta"][:6] for tw in twins][:2], flush=True)
    json.dump(out, open(DIR + f"cx_05b_random_triples_{tname}_{conv}.json", "w"), indent=1)
print("triples with twins:", sum(1 for t in out["triples"] if t["n_twins"] > 0), "of", NTR)
