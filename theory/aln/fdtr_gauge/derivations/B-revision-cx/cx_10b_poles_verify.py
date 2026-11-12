"""cx_10b: (1) verification of right-half-plane poles for eps < 0 directly on Z (40 digits): |Z(p)| near the
located real pole p0 grows like 1/|p - p0|, and Z changes sign across p0 on the real p axis;
(2) thick-film asymptotics: pole of the bare film at p ~ (1 + b/b_s)/(|eps| tau_l) as xi -> infinity
(tau_R = tau_l, R = 0); (3) passivity on the imaginary axis for eps > 0: min over omega of Re Z/|Z|
(bare film read-out T(0) and stack) for random parameters; (4) |1/Z(i omega)| minima (poles on the axis).
"""
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
mp.mp.dps = 40
out = {}


def Z_at_p(spec, th, p):
    s = mp.sqrt(mp.mpc(p))
    return cx.response(spec, [mp.mpf(x) for x in th], s, cx.MP_OPS)


# (1) direct verification on Z: poles located as zeros of the cleared denominator D (bare: B L tanh A + N;
# stack: as in cx_10) on the real p axis, then |Z| must scale like 1/|p - p0|.
def D_mp(kind, prm, s, xiA, bA):
    xi, b, tR, tl, tn, to, R1, R2 = prm
    p = s * s
    g = (1 + tR * p) / (1 + tl * p); sg = mp.sqrt(g)
    A = xi * s * sg; B = b * s / sg; N = 1 + tn * p; N0 = 1 + to * p; L = R2 + 1 / s; T = mp.tanh(A)
    if kind == "bare":
        return B * L * T + N
    TA = mp.tanh(xiA * s)
    return bA * s * TA * (R1 * (B * L * T + N) + N0 * (L + N * T / B)) + (B * L * T + N)

ver = []
for kind, eps, th in [("bare", -0.5, (0.6477, 2.098, 0.288, 0.288, 0.0)), ("bare", -1.0, (0.8298, 0.3996, 1.2365, 1.2365, 0.0)),
                      ("stack", -0.5, (1.1908, 0.5106, 0.1825, 0.1825, 0.0153, 0.0)), ("stack", -0.1, None)]:
    if th is None:
        thv, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=10e-12)
        th = thv
        label = "AlN note-16 stack eps=-0.1"
    else:
        xiA, bA = 0.2, 3.0
        label = f"{kind} eps={eps}"
    sp = cx.Spec(kind, "tied", 0.0 if kind == "bare" else eps, eps, xiA=xiA, bA=bA, R_known=None)
    th = [mp.mpf(float(x)) for x in th]
    if kind == "bare":
        prm = (th[0], th[1], th[2], th[3], eps * th[3], 0, 0, th[4])
    else:
        prm = (th[0], th[1], th[2], th[3], eps * th[3], eps * th[3], th[4], th[5])
    prm = [mp.mpf(float(x)) for x in prm]
    ss = [mp.mpf(10) ** (mp.mpf(k) / 200) for k in range(-200, 801)]
    vals = [D_mp(kind, prm, s, mp.mpf(xiA), mp.mpf(bA)) for s in ss]
    poles = []
    for k in range(len(ss) - 1):
        if vals[k] * vals[k + 1] < 0:
            s0 = mp.findroot(lambda s: D_mp(kind, prm, s, mp.mpf(xiA), mp.mpf(bA)), (ss[k], ss[k + 1]), solver="anderson")
            p0 = s0 * s0
            za = abs(Z_at_p(sp, th, p0 * (1 + mp.mpf("1e-8"))))
            zb = abs(Z_at_p(sp, th, p0 * (1 + mp.mpf("1e-12"))))
            zc = abs(Z_at_p(sp, th, p0 * (1 - mp.mpf("1e-12"))))
            poles.append(dict(p0=float(p0), absZ_at_1e8=float(za), absZ_at_1e12=float(zb), absZ_at_m1e12=float(zc),
                              scaling_ratio=float(zb / za)))
    ver.append(dict(case=label, theta=[float(x) for x in th], real_poles=poles))
    print(label, "real RHP poles:", [("p0=%.6g" % q["p0"], "|Z|(+1e-8)=%.3e" % q["absZ_at_1e8"],
                                      "|Z|(+1e-12)=%.3e" % q["absZ_at_1e12"], "ratio %.4g" % q["scaling_ratio"]) for q in poles], flush=True)
out["direct_verification"] = ver

# (2) thick-film asymptotics
asy = []
for eps in [-0.2, -0.5, -1.0]:
    for b in [0.5, 2.0]:
        t = mp.mpf("0.3")
        pred = (1 + b) / (abs(eps) * t)
        sp = cx.Spec("bare", "tied", 0.0, eps, R_known=0.0)
        for xi in [2, 5, 10]:
            th = (xi, b, t, t)
            f = lambda p: (1 / Z_at_p(sp, th, p)).real
            p0 = mp.findroot(f, pred * mp.mpf("1.0001"), solver="secant")
            asy.append(dict(eps=eps, b=b, tau=float(t), xi=xi, pole=float(p0), predicted=float(pred)))
            print(f"thick film eps={eps} b={b} xi={xi}: pole p0 = {float(p0):.6g}, prediction (1+b/bs)/(|eps| tau) = {float(pred):.6g}",
                  flush=True)
out["thick_film_asymptotics"] = asy

# (3) passivity on the axis and (4) near-poles on the axis, eps >= 0
rng = np.random.default_rng(3)
om = np.logspace(-4, 4, 4001)
pas = {}
for kind in ("bare", "stack"):
    for eps in [0.0, 0.084, 0.3, 0.5, 0.9, 1.0, 1.2, 2.5, 5.0, 20.0]:
        worst, worst_th, mininv = 1.0, None, np.inf
        for t in range(200):
            xi = np.exp(rng.uniform(np.log(.2), np.log(5))); b = np.exp(rng.uniform(np.log(.2), np.log(5)))
            tR = np.exp(rng.uniform(np.log(1e-3), np.log(3))); tl = np.exp(rng.uniform(np.log(1e-3), np.log(3)))
            R1 = np.exp(rng.uniform(np.log(1e-2), np.log(2))); R2 = 0.0 if t % 3 == 0 else np.exp(rng.uniform(np.log(1e-2), np.log(2)))
            if kind == "bare":
                sp = cx.Spec("bare", "tied", 0.0, eps, R_known=None); th = (xi, b, tR, tl, R2)
            else:
                sp = cx.Spec("stack", "tied", eps, eps, xiA=0.2, bA=3.0); th = (xi, b, tR, tl, R1, R2)
            with np.errstate(all="ignore"):
                z = cx.Z_np(sp, th, om)
            m = np.min(z.real / np.abs(z))
            if m < worst:
                worst, worst_th = m, th
            mininv = min(mininv, np.min(np.abs(1 / z) * np.sqrt(om)))
        pas[f"{kind} eps={eps:g}"] = dict(min_ReZ_over_absZ=float(worst), at=[float(x) for x in worst_th],
                                          min_abs_invZ_times_sqrt_omega=float(mininv))
        print(f"{kind} eps={eps:g}: min Re Z/|Z| over omega and 200 parameter sets = {worst:.4f}; min |sqrt(w)/Z| = {mininv:.3e}",
              flush=True)
out["passivity_axis"] = pas
with open(DIR + "cx_10b_poles_verify.json", "w") as f:
    json.dump(out, f, indent=1)
