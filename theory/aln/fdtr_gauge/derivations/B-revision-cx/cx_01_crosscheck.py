"""cx_01: consistency of the independent model.

(a) exact forward-mode derivatives vs 80-digit central differences;
(b) numpy vs mpmath values;
(c) independent model vs the twin `src/twin_fdtr.py` (k = 0) in SI units: Au/AlN/sapphire stack and bare
    AlN/sapphire, conventions local / entropic / kinetic (eps = 0.3, -0.5, 2), observable T and theta_e.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
sys.path.insert(0, r"C:/Users/Koussay/these/src")
import cxmodel as cx

out = {}

# (a) derivatives
mp.mp.dps = 80
rng = np.random.default_rng(1)
worst = 0
for spec in [cx.Spec("stack", "tied", 0, 0), cx.Spec("stack", "tied", 1, 1), cx.Spec("stack", "tied", 0.3, -0.7),
             cx.Spec("stack", "free"), cx.Spec("bare", "tied", 0, 1), cx.Spec("bare", "tied", 1, 0.5),
             cx.Spec("bare", "free", 0, 0)]:
    for trial in range(3):
        th = [mp.mpf(float(x)) for x in rng.uniform(0.1, 2.0, spec.n)]
        for w in [mp.mpf("0.01"), mp.mpf("1.3"), mp.mpf("40")]:
            z, d = cx.ZJ_mp(spec, th, w)
            h = mp.mpf("1e-30")
            for k in range(spec.n):
                tp = list(th); tm = list(th)
                tp[k] += h; tm[k] -= h
                fd = (cx.Z_mp(spec, tp, w) - cx.Z_mp(spec, tm, w)) / (2 * h)
                err = abs(fd - d[k]) / (abs(d[k]) + abs(z) * mp.mpf("1e-10"))
                worst = max(worst, err)
out["a_dual_vs_fd_worst_rel"] = float(worst)
print("(a) dual vs FD worst rel err", float(worst))

# (b) numpy vs mpmath
mp.mp.dps = 30
worst = 0
for spec in [cx.Spec("stack", "tied", 1, 1), cx.Spec("bare", "tied", 0, 0.084), cx.Spec("stack", "free")]:
    th = rng.uniform(0.1, 2.0, spec.n)
    om = np.logspace(-3, 3, 13)
    zn, Jn = cx.ZJ_np(spec, th, om)
    for i, w in enumerate(om):
        zm, dm = cx.ZJ_mp(spec, th, w)
        worst = max(worst, abs(complex(zm) - zn[i]) / abs(zn[i]))
        for k in range(spec.n):
            worst = max(worst, abs(complex(dm[k]) - Jn[k, i]) / (abs(Jn[k, i]) + 1e-300))
out["b_numpy_vs_mp_worst_rel"] = float(worst)
print("(b) numpy vs mp worst rel", worst)

# (c) twin comparison
import twin_fdtr as tw

lam, C = 188.6, 2.41e6
au = tw.Material(150.0, 2.49e6)
sap = tw.Material(35.0, 3.03e6)
a = lam / C
bs = np.sqrt(35.0 * 3.03e6)
freqs = np.logspace(4, 11, 29)
om = 2 * np.pi * freqs
res = {}
for alpha in (1 / 3, 2.0):
    for (tR, tl) in [(25e-12, 10e-12), (25e-12, 25e-12), (0.0, 20e-12)]:
        ell = np.sqrt(tl * a / (1 + alpha))
        aln = tw.Material(lam, C, tau_r=tR, ell=ell, alpha=alpha)
        for conv, eps in [("local", 0.0), ("entropic", 1.0), ("kin+0.3", 0.3), ("kin-0.5", -0.5), ("kin2", 2.0)]:
            if conv == "local":
                bnd = tw.GKBoundary(temperature="local")
            elif conv == "entropic":
                bnd = tw.GKBoundary(temperature="entropic")
            else:
                bnd = tw.GKBoundary(temperature="kinetic", c_n=eps * tl / C, c_t=0.0)
            for R1, R2 in [(0.0, 0.0), (1e-8, 2e-8)]:
                st = tw.Stack((tw.Layer(au, 80e-9, R1), tw.Layer(aln, 500e-9, R2, bnd)), sap)
                zt = tw.kernel(freqs, 0.0, st)
                spec = cx.Spec("stack", "tied", eps, eps, xiA=80e-9 / np.sqrt(150 / 2.49e6),
                               bA=np.sqrt(150 * 2.49e6), bs=bs)
                th = (500e-9 / np.sqrt(a), np.sqrt(lam * C), tR, tl, R1, R2)
                zm = cx.Z_np(spec, th, om)
                e1 = float(np.max(np.abs(zm - zt) / np.abs(zt)))
                # bare film, observable T and theta_e
                sb = tw.Stack((tw.Layer(aln, 500e-9, R2, bnd),), sap)
                zbt = tw.kernel(freqs, 0.0, sb)
                zbe = tw.kernel(freqs, 0.0, sb, observable="theta_e")
                specT = cx.Spec("bare", "tied", 0.0, eps, R_known=R2, bs=bs)
                specE = cx.Spec("bare", "tied", 1.0, eps, R_known=R2, bs=bs)
                thb = th[:4]
                e2 = float(np.max(np.abs(cx.Z_np(specT, thb, om) - zbt) / np.abs(zbt)))
                e3 = float(np.max(np.abs(cx.Z_np(specE, thb, om) - zbe) / np.abs(zbe)))
                key = f"alpha={alpha:.3g} tR={tR:g} tl={tl:g} {conv} R=({R1:g},{R2:g})"
                res[key] = (e1, e2, e3)
worst = np.max(np.array(list(res.values())), axis=0)
out["c_twin_vs_model_worst_rel(stack,bare_T,bare_theta_e)"] = [float(x) for x in worst]
out["c_detail"] = res
print("(c) twin vs model worst (stack, bare T, bare theta_e):", worst)

with open(r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/cx_01_crosscheck.json", "w") as f:
    json.dump(out, f, indent=1)
