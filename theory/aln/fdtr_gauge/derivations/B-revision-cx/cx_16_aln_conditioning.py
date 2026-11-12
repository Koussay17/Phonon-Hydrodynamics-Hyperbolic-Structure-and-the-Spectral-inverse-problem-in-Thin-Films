"""cx_16: conditioning of the Au/AlN/sapphire note-16 stack (1D) in the FDTR band, 50 digits.

(1) three-frequency Jacobians (columns d log Z / d log th, th = (xi, b, tR, tl, R1, R2)) at 300 random band triples
    (10 kHz - 200 MHz): best and median s_min/s_max -> how far "three generic frequencies suffice" is from usable;
(2) 80-frequency Jacobian: singular values and the weakest direction, both conventions, tau_l = 0, 10, 25 ps.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

mp.mp.dps = 50
rng = np.random.default_rng(16)
out = {}
for tl_ps in (0.0, 10.0, 25.0):
    th, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=tl_ps * 1e-12)
    th = list(th)
    if tl_ps == 0.0:
        th[3] = 0.0
    for conv, eps in (("local", 0.0), ("entropic", 1.0)):
        sp = cx.Spec("stack", "tied", eps, eps, xiA=xiA, bA=bA)
        # log-parameter scaling (tl = 0: raw column scaled by the 10 ps value)
        scale = [x if x > 0 else 10e-12 / t0 for x in th]

        def J_of(fr):
            om = [mp.mpf(repr(float(2 * np.pi * f * t0))) for f in fr]
            J = cx.realjac_mp(sp, th, om, rel=True)
            for i in range(J.rows):
                for k in range(J.cols):
                    J[i, k] *= scale[k]
            return J
        vals = []
        for t in range(300):
            fr = np.sort(10 ** rng.uniform(4, np.log10(2e8), 3))
            S = mp.svd_r(J_of(fr), compute_uv=False)
            sv = sorted([S[i] for i in range(len(S))], reverse=True)
            vals.append(float(sv[-1] / sv[0]))
        fr80 = np.logspace(4, np.log10(2e8), 80)
        J80 = J_of(fr80)
        U, S, V = mp.svd_r(J80)
        sv = [S[i] for i in range(len(S))]
        order = sorted(range(len(sv)), key=lambda i: -sv[i])
        weakest = [float(V[order[-1], k]) for k in range(6)]
        key = f"{conv} tau_l={tl_ps:g}ps"
        out[key] = dict(three_freq_best=max(vals), three_freq_median=float(np.median(vals)),
                        sv80=[float(sv[i] / sv[order[0]]) for i in order], weakest_direction_logparams=weakest)
        print(f"{key}: 3-frequency s_min/s_max best {max(vals):.2e} median {np.median(vals):.2e}; 80-frequency s/s1 =",
              ["%.1e" % (sv[i] / sv[order[0]]) for i in order], "weakest dir (xi,b,tR,tl,R1,R2) =", np.round(weakest, 3).tolist(),
              flush=True)
json.dump(out, open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_16_aln_conditioning.json", "w"), indent=1)
