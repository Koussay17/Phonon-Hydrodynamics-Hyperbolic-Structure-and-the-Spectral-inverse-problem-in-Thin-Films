"""cx_14: bare film, read-out T(0), entropic bottom face (eps = 1), R known: two frequencies on the resonance Delta.
50-digit 4 x 4 real Jacobian at random frequency pairs; normalised s_min, and the same for eps = 0.084, 1e-2 (lifting
of the resonance proportional to eps?) and for the local convention (rank 3 expected)."""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
rng = np.random.default_rng(12)
out = []
for eps in [1.0, 0.084, 0.01, 0.0]:
    for Rk in [0.0, 0.3]:
        sp = cx.Spec("bare", "tied", 0.0, eps, R_known=Rk)
        vals = []
        for i in range(6):
            t = float(rng.uniform(0.05, 1)); th = [float(rng.uniform(0.4, 1.6)), float(rng.uniform(0.4, 2.5)), t, t]
            for k in range(5):
                om = [mp.mpf(float(x)) for x in 10 ** rng.uniform(-1.5, 1.5, 2)]
                J = cx.realjac_mp(sp, th, om)
                sv, _ = cx.svd_mp(J)
                vals.append(float(sv[-1] / sv[0]))
        rec = dict(eps=eps, R_known=Rk, smin_min=min(vals), smin_median=float(np.median(vals)), smin_max=max(vals))
        out.append(rec)
        print(f"eps={eps} R={Rk}: two-frequency normalised s_min on Delta: min {min(vals):.2e} median {np.median(vals):.2e} max {max(vals):.2e}",
              flush=True)
json.dump(out, open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_14_two_freq_entropic.json", "w"), indent=1)
