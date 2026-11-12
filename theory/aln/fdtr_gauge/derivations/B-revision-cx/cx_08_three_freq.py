"""cx_08: three-frequency local rank (Theorem T4(iv), Conjecture T4').

(1) Independent computation of the low-frequency Laurent coefficients c_k of the film block
    W - R1 - 1/(b_s s) = sum_k c_k s^k  by trapezoidal Cauchy integrals in s (mpmath, exact parameter
    derivatives by forward mode), and of G'_12345 = det[grad c_1..c_5] (columns xi, b, tR, tl, R2);
    comparison with the closed forms (P for local, Q for entropic) of the reduction record.
(2) Points ON the hypersurface {G' = 0} off Delta (local: P = 0; entropic: Q = 0): smallest normalised
    singular value of the 3-frequency Jacobian J_3 at 40 random triples, 50 digits.  A point where every
    triple gives s_min = 0 would refute Conjecture T4'.
(3) The leading-order formula det J_3 ~ sigma^15 t(lam,mu) G': ratio check at small sigma.
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
mp.mp.dps = 50
out = {}


def W_dual(spec, th, s):
    """Film block W - R1 (bare film block with N0 = top factor), at complex s, dual in (xi,b,tR,tl,R2)."""
    xi, b, tR, tl, R2 = th
    sp = cx.Spec("bare", "tied", spec.eps0, spec.eps, R_known=None, bs=1)
    d = cx.seed([mp.mpf(x) for x in (xi, b, tR, tl, R2)], mp.mpc(1), mp.mpc(0))
    return cx.response(sp, d, s, cx.MP_OPS)


def laurent_grads(spec, th, K=6, rho=mp.mpf("0.25"), M=96):
    """c_1..c_K and their gradients, by the trapezoidal rule on |s| = rho."""
    cs = [mp.mpc(0)] * (K + 1)
    gs = [[mp.mpc(0)] * 5 for _ in range(K + 1)]
    for j in range(M):
        s = rho * mp.expjpi(mp.mpf(2 * j) / M)
        w = W_dual(spec, th, s)
        v = w.v - 1 / s
        for k in range(K + 1):
            sk = s ** (-k)
            cs[k] += v * sk / M
            for i in range(5):
                gs[k][i] += w.d[i] * sk / M
    return cs, gs


def Gprime(spec, th):
    xi, b = mp.mpf(th[0]), mp.mpf(th[1])
    rho = mp.mpf("0.3") * (mp.atanh(min(1 / b, mp.mpf("0.999"))) if b > 1 else mp.pi / 2) / xi
    rho = min(rho, mp.mpf("0.15"))
    cs, gs = laurent_grads(spec, th, rho=rho, M=160)
    cs2, gs2 = laurent_grads(spec, th, rho=rho * mp.mpf("0.6"), M=160)
    dif = max(abs(gs[k][i] - gs2[k][i]) / (abs(gs[k][i]) + mp.mpf("1e-30")) for k in range(1, 6) for i in range(5))
    if dif > mp.mpf("1e-25"):
        print("   WARNING Cauchy radius check: rel diff %.1e" % float(dif))
    Mx = mp.matrix([[gs[k][i].real for i in range(5)] for k in range(1, 6)])
    imag = max(abs(gs[k][i].imag) for k in range(1, 6) for i in range(5))
    return mp.det(Mx), cs, imag


def P_local(b, R, tl, tR, xi):
    return (90*b*R*tl**2 - 90*b*R*tl*tR + 15*b**4*R**4*xi - 30*b**2*R**2*tl*xi - 45*b**2*tl**2*xi + 30*b**2*R**2*tR*xi
            + 45*tl*tR*xi + 45*b**2*tl*tR*xi - 45*tR**2*xi + 24*b**3*R**3*xi**2 - 30*b**5*R**3*xi**2 - 6*b*R*tl*xi**2
            + 60*b**3*R*tl*xi**2 + 6*b*R*tR*xi**2 - 60*b**3*R*tR*xi**2 + 11*b**2*R**2*xi**3 - 31*b**4*R**2*xi**3
            - 3*tl*xi**3 + 3*b**2*tl*xi**3 + 3*tR*xi**3 - 3*b**2*tR*xi**3 + 2*b*R*xi**4 - 4*b**3*R*xi**4 + 2*b**5*R*xi**4)


def Q_ent(b, R, tl, tR, xi):
    return (-1350*b*R*tl**2*tR + 1350*b*R*tl*tR**2 - 45*b**4*R**4*tl*xi + 630*b**4*R**2*tl**2*xi + 675*b**4*tl**3*xi
            + 225*b**4*R**4*tR*xi - 450*b**2*R**2*tl*tR*xi - 450*b**4*R**2*tl*tR*xi - 675*tl**2*tR*xi
            + 2025*b**2*tl**2*tR*xi - 675*b**4*tl**2*tR*xi + 450*b**2*R**2*tR**2*xi + 1350*tl*tR**2*xi
            - 2025*b**2*tl*tR**2*xi - 675*tR**3*xi + 90*b**5*R**5*xi**2 - 90*b**3*R**3*tl*xi**2 + 540*b**5*R**3*tl*xi**2
            + 630*b**3*R*tl**2*xi**2 + 450*b**5*R*tl**2*xi**2 + 360*b**3*R**3*tR*xi**2 - 450*b**5*R**3*tR*xi**2
            - 90*b*R*tl*tR*xi**2 + 540*b**3*R*tl*tR*xi**2 - 450*b**5*R*tl*tR*xi**2 - 900*b**3*R*tR**2*xi**2
            + 150*b**4*R**4*xi**3 - 15*b**2*R**2*tl*xi**3 + 675*b**4*R**2*tl*xi**3 + 45*b**2*tl**2*xi**3
            + 45*b**4*tl**2*xi**3 + 135*b**2*R**2*tR*xi**3 - 435*b**4*R**2*tR*xi**3 - 45*tl*tR*xi**3
            + 180*b**2*tl*tR*xi**3 + 45*b**4*tl*tR*xi**3 + 45*tR**2*xi**3 + 45*b**2*tR**2*xi**3 + 90*b**3*R**3*xi**4
            + 30*b**5*R**3*xi**4 + 150*b**3*R*tl*xi**4 - 30*b**5*R*tl*xi**4 + 30*b*R*tR*xi**4 + 30*b**3*R*tR*xi**4
            + 60*b**5*R*tR*xi**4 + 24*b**2*R**2*xi**5 + 24*b**4*R**2*xi**5 + 2*b*R*xi**6 - 2*b**5*R*xi**6)


rng = np.random.default_rng(5)
spL = cx.Spec("stack", "tied", 0, 0, xiA=0.2, bA=3.0)
spE = cx.Spec("stack", "tied", 1, 1, xiA=0.2, bA=3.0)

# (1) closed forms vs numerical G'
cmp = []
for i in range(6):
    xi, b, tR, tl, R = [mp.mpf(float(x)) for x in (rng.uniform(.4, 1.5), rng.uniform(.4, 2.5), rng.uniform(.05, .8),
                                                   rng.uniform(.05, .8), rng.uniform(0, .8))]
    GL, cs, im = Gprime(spL, (xi, b, tR, tl, R))
    GE, _, imE = Gprime(spE, (xi, b, tR, tl, R))
    fL = 2 * (tl - tR) * xi ** 7 * P_local(b, R, tl, tR, xi) / (135 * b ** 2)
    fE = -2 * xi ** 7 * Q_ent(b, R, tl, tR, xi) / (2025 * b ** 2)
    cmp.append(dict(point=[float(x) for x in (xi, b, tR, tl, R)], G_local_num=float(GL), G_local_formula=float(fL),
                    G_ent_num=float(GE), G_ent_formula=float(fE), imag_parts=float(max(im, imE))))
    print("G' local num %.12e formula %.12e | entropic num %.12e formula %.12e | imag %.1e" %
          (float(GL), float(fL), float(GE), float(fE), float(max(im, imE))), flush=True)
out["closed_form_check"] = cmp


# (2) points on {G' = 0}
def triples_smin(spec, th, ntr=40):
    vals = []
    for t in range(ntr):
        om = [mp.mpf(float(x)) for x in 10 ** rng.uniform(-1.5, 1.5, 3)]
        J = cx.realjac_mp(spec, th, om, rel=True)
        sv, _ = cx.svd_mp(J)
        vals.append(float(sv[-1] / sv[0]))
    return vals


hyp = []
for conv, sp, fun in (("local", spL, P_local), ("entropic", spE, Q_ent)):
    for i in range(8):
        # pick (xi, b, tl, R2, R1) and solve for tR on the hypersurface (polynomial in tR)
        for attempt in range(50):
            xi, b, tl, R = [mp.mpf(float(x)) for x in (rng.uniform(.4, 1.5), rng.uniform(.4, 2.5), rng.uniform(.02, .8),
                                                        rng.choice([0.0, float(rng.uniform(0.01, 0.6))]))]
            f = lambda x: fun(b, R, tl, x, xi)
            # scan tR in (0, 3) for sign changes
            xs = [mp.mpf(k) / 200 * 3 for k in range(1, 201)]
            ys = [f(x) for x in xs]
            roots = []
            for k in range(len(xs) - 1):
                if ys[k] * ys[k + 1] < 0:
                    roots.append(mp.findroot(f, (xs[k], xs[k + 1]), solver="anderson"))
            roots = [r for r in roots if abs(r - tl) > 0.02]
            if roots:
                break
        tR = roots[0]
        R1 = mp.mpf(float(rng.uniform(0.05, 0.6)))
        th = [xi, b, tR, tl, R1, R]
        G, _, _ = Gprime(sp, (xi, b, tR, tl, R))
        # function-level rank (40 frequencies)
        OM = [mp.mpf(10) ** (mp.mpf(-3) + 6 * mp.mpf(k) / 39) for k in range(40)]
        sv40, _ = cx.svd_mp(cx.realjac_mp(sp, th, OM))
        sm = triples_smin(sp, th)
        rec = dict(conv=conv, theta=[float(x) for x in th], G_prime_numeric=float(G), fun_level_smin=float(sv40[-1] / sv40[0]),
                   triples_smin_min=min(sm), triples_smin_median=float(np.median(sm)), triples_smin_max=max(sm))
        hyp.append(rec)
        print(conv, "on {G'=0}:", [mp.nstr(x, 6) for x in th], "G'num %.1e" % float(G),
              "fun-level smin %.2e" % float(sv40[-1] / sv40[0]), "3-freq smin min/med/max %.2e %.2e %.2e" %
              (min(sm), np.median(sm), max(sm)), flush=True)
out["on_hypersurface"] = hyp

# (3) leading-order check of det J_3 at small sigma (generic point, local)
th = [mp.mpf(x) for x in ("0.8", "2", "0.5", "0.1", "0.3", "0.2")]
lam_, mu_ = mp.mpf(2), mp.mpf(3)
G, _, _ = Gprime(spL, th[:4] + [th[5]])
tfac = (lam_ - 1) ** 2 * lam_ * (1 + lam_ ** 2) * (lam_ - mu_) ** 2 * (mu_ - 1) ** 2 * mu_ * (1 + mu_ ** 2) * (lam_ ** 2 + mu_ ** 2) / (2 * mp.sqrt(2))
lead = []
for sig in [mp.mpf("1e-2"), mp.mpf("3e-3"), mp.mpf("1e-3")]:
    om = [sig ** 2, (lam_ * sig) ** 2, (mu_ * sig) ** 2]
    # Jacobian of W (not Z): use the film block plus R1 column
    rows_re, rows_im = [], []
    for w in om:
        s = mp.sqrt(mp.mpc(0, w))
        wd = W_dual(spL, th[:4] + [th[5]], s)
        d = wd.d
        full = [d[0], d[1], d[2], d[3], mp.mpc(1), d[4]]
        rows_re.append([x.real for x in full]); rows_im.append([x.imag for x in full])
    D = mp.det(mp.matrix(rows_re + rows_im))
    ratio = D / (sig ** 15 * tfac * G)
    lead.append(dict(sigma=float(sig), det=float(D), ratio_to_leading=float(ratio)))
    print("sigma %.0e det %.4e ratio det/(sigma^15 t G') = %.6f" % (float(sig), float(D), float(ratio)), flush=True)
out["leading_order"] = lead
with open(DIR + "cx_08_three_freq.json", "w") as f:
    json.dump(out, f, indent=1)
