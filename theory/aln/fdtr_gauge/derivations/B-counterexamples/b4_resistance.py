"""B4: interface resistance R between film and substrate; rank of theta5 = (xi1, b, tR, tl, R) -> Z.

Zs = R + 1/(bs sqrt p);  Z = (Zs + tanh A/B)/(1 + B Zs tanh A).  Units xi1 = bs = 1 at the point
(R in units of xi1/bs, i.e. of 1/(bs sqrt p) at p = 1/xi1^2).  Columns: d/dln xi1, d/dln b,
tau d/dtau (or xi1^2 d/dtau at tau = 0), R d/dR (or (xi1/b) d/dR at R = 0).
Tests
 1. random full-band scan of s5 (equilibrated and raw), 50-digit re-evaluation of the smallest;
 2. exact strata at 50 digits: resonance, matched effusivity, R = 0, combinations;
 3. near-resonance scaling of s5;
 4. limits R -> infinity (decoupled substrate), R -> 0, thick and thin film;
 5. N = 3 frequencies: ball-arithmetic certificate of a nonzero 5x5 minor at random points;
 6. adversarial minimisation of s5 away from resonance (hard constraint |ln tR/tl| >= 1).
"""
import sys, os, json, time
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from flint import acb, arb, arb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


LAST_OK = [True]


def svals5(beta, tR, tl, R, pad=1.0, metric="eq", wlist=None):
    ok = True
    if wlist is None:
        w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, R=R, pad=pad, nmax=60000)
        # extend the window to the resistance frequency (R bs)^-2 (where R ~ 1/(bs sqrt p))
        if R > 0:
            fR = 1.0 / R ** 2
            lo, hi = min(w[0], fR / 10 ** pad), max(w[-1], fR * 10 ** pad)
            w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, R=R, wmin=lo, wmax=hi, nmax=60000)
    else:
        w = np.asarray(wlist); wts = np.ones_like(w)
    LAST_OK[0] = ok                      # False: thermal-wave resonances not resolved (grid capped)
    M = zm.real_rel_jac_np(w, 1.0, beta, tR, tl, 1.0, R, "code", True, "log", wts)
    if metric == "eq":
        M = M / np.linalg.norm(M, axis=0)
    return np.linalg.svd(M, compute_uv=False), w, wts


def svals5_mp(beta, tR, tl, R, w, wts, dps=50, equilibrate=True):
    mp.mp.dps = dps
    M = zm.real_rel_jac_mp([mp.mpf(x) for x in w], mp.mpf(1), mp.mpf(beta), mp.mpf(tR), mp.mpf(tl), 1,
                           mp.mpf(R), "code", True, "log", [mp.mpf(x) for x in wts])
    if equilibrate:
        for j in range(M.cols):
            nrm = mp.sqrt(mp.fsum(M[i, j] ** 2 for i in range(M.rows)))
            if nrm > 0:
                for i in range(M.rows):
                    M[i, j] /= nrm
    return [float(x) for x in zm.svals_mp(M)]


def jac5_acb(w, xi1, b, tR, tl, R, bs=1):
    p = acb(0, w)
    s = p.sqrt()
    r = ((1 + tR * p) / (1 + tl * p)).sqrt()
    A = xi1 * s * r
    B = b * s / r
    E = (-2 * A).exp()
    ome = -((-2 * A).expm1())
    y = B * (R + 1 / (bs * s))
    rho = (y - 1) / (y + 1)
    D1 = 1 - rho * E
    Z = (1 + rho * E) / (B * D1)
    u = -4 * A * rho * E / (B * D1 ** 2)
    v = -ome * ((1 + E) * (1 + y ** 2) + 2 * y * ome) / (B * D1 ** 2 * (y + 1) ** 2)
    gR = p / (2 * (1 + tR * p))
    gl = p / (2 * (1 + tl * p))
    dR = 4 * E / (D1 ** 2 * (y + 1) ** 2)
    return Z, [u, v, gR * (u - v), -gl * (u - v), dR]


def certify3(th, ws, prec=200):
    ctx.prec = prec
    xi1, b, tR, tl, R = (arb(str(x)) for x in th)
    rows = []
    for w in ws:
        Z, J = jac5_acb(arb(str(w)), xi1, b, tR, tl, R)
        c = [arb(1), arb(1), tR, tl, R]
        jr = [J[k] * c[k] / Z for k in range(5)]
        rows.append([x.real for x in jr])
        rows.append([x.imag for x in jr])
    best = None
    # try the 6 minors obtained by deleting one row
    for drop in range(6):
        sub = [rows[i] for i in range(6) if i != drop]
        det = arb_mat(sub).det()
        if not det.contains(0):
            return det, True, drop
        best = det
    return best, False, None


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(17)
    out = {}
    # 1. random scan
    rows = []
    for k in range(3000):
        beta = 10 ** rng.uniform(-2, 2); tR = 10 ** rng.uniform(-3, 3); tl = 10 ** rng.uniform(-3, 3)
        R = 10 ** rng.uniform(-3, 3)
        if k % 6 == 5:
            tl = 0.0
        S, w, wts = svals5(beta, tR, tl, R)
        okk = LAST_OK[0]
        Sr, _, _ = svals5(beta, tR, tl, R, metric="raw")
        rows.append([beta, tR, tl, R, S[-1], Sr[-1], len(w), float(okk)])
    rows = np.array(rows)
    nunres = int((rows[:, 7] == 0).sum())
    rows = rows[rows[:, 7] == 1]
    order = np.argsort(rows[:, 4])
    out["scan"] = {"n": len(rows), "n_unresolved_dropped": nunres, "median_s5_eq": float(np.median(rows[:, 4])),
                   "min_s5_eq": float(rows[order[0], 4]), "worst10": rows[order[:10]].tolist()}
    print("1. scan: %d resolved points (%d unresolved dropped): median s5(eq) = %.3e, min = %.3e at %s" % (
        len(rows), nunres, out["scan"]["median_s5_eq"], rows[order[0], 4], rows[order[0], :4]))
    checks = []
    for i in order[:6]:
        beta, tR, tl, R = rows[i, :4]
        S, w, wts = svals5(beta, tR, tl, R)
        if len(w) > 3000:
            idx = np.unique(np.linspace(0, len(w) - 1, 3000).astype(int))
            w, wts = w[idx], wts[idx] * len(w) / len(idx)
            S, _, _ = svals5(beta, tR, tl, R, wlist=w)
            S = None
        smp = svals5_mp(beta, tR, tl, R, w, wts)
        Md = zm.real_rel_jac_np(np.asarray(w), 1.0, beta, tR, tl, 1.0, R, "code", True, "log", np.asarray(wts))
        Md = Md / np.linalg.norm(Md, axis=0)
        sd = np.linalg.svd(Md, compute_uv=False)
        checks.append({"point": rows[i, :4].tolist(), "n_omega": len(w), "s5_double_same_grid": float(sd[-1]),
                       "s5_mp50": smp[-1]})
        print("   worst point beta=%.3g tR=%.3g tl=%.3g R=%.3g (%d nodes): s5 double %.4e, mp50 %.4e" %
              (beta, tR, tl, R, len(w), sd[-1], smp[-1]))
    out["scan"]["mp_checks"] = checks
    # 2. exact strata (raw, 50 digits)
    strata = []
    for (beta, tR, tl, R, label) in ((2.0, 0.5, 0.1, 0.3, "generic"),
                                     (2.0, 0.5, 0.5, 0.3, "resonance"),
                                     (1.0, 0.5, 0.5, 0.3, "resonance, beta=1, R>0"),
                                     (1.0, 0.5, 0.5, 0.0, "resonance, beta=1, R=0"),
                                     (2.0, 0.5, 0.1, 0.0, "R=0"),
                                     (1.0, 0.5, 0.1, 0.0, "beta=1, R=0"),
                                     (2.0, 0.5, 0.0, 0.3, "Cattaneo"),
                                     (2.0, 0.0, 0.3, 0.3, "tR=0"),
                                     (2.0, 0.0, 0.0, 0.3, "Fourier film + R")):
        w, wts, ok = zm.resolved_grid(1.0, beta, tR, tl, 1.0, R=R, pad=1.0)
        w = w[::max(1, len(w) // 120)]; wts = np.ones_like(w)
        sv = svals5_mp(beta, tR, tl, R, w, wts, dps=50, equilibrate=False)
        mp.mp.dps = 50
        strata.append({"label": label, "point": [beta, tR, tl, R], "s_raw_mp50": sv})
        print("2. %-26s s/s1 (50 digits) =" % label, ["%.2e" % (x / sv[0]) for x in sv])
    out["strata"] = strata
    # 3. near-resonance scaling
    res = []
    for d in 10.0 ** np.arange(-1, -8.5, -1):
        S, w, wts = svals5(2.0, 0.5 * np.exp(d / 2), 0.5 * np.exp(-d / 2), 0.3, metric="raw")
        res.append([d, float(S[-1]), float(S[-2])])
    res = np.array(res)
    slope = float(np.polyfit(np.log(res[2:, 0]), np.log(res[2:, 1]), 1)[0])
    out["near_resonance"] = {"rows[d,s5,s4]": res.tolist(), "slope_s5": slope}
    print("3. near resonance (beta=2, tbar=0.5, R=0.3): s5 ~ d^%.4f, s5/d = %.3e" % (slope, res[-1, 1] / res[-1, 0]))
    # 4. limits in R
    lim = []
    for R in (1e-4, 1e-2, 1.0, 1e2, 1e4, 1e6):
        S, w, wts = svals5(2.0, 0.5, 0.1, R, metric="raw")
        Se, _, _ = svals5(2.0, 0.5, 0.1, R, metric="eq")
        lim.append([R, float(S[-1]), float(Se[-1])])
        print("4. R=%.0e: s5 raw %.3e, eq %.3e" % (R, S[-1], Se[-1]))
    out["R_limits"] = lim
    # 5. N = 3 certificates
    cert = []
    for k in range(100):
        th = (1.0, float("%.5g" % 10 ** rng.uniform(-1.5, 1.5)), float("%.5g" % 10 ** rng.uniform(-2, 2)),
              float("%.5g" % 10 ** rng.uniform(-2, 2)), float("%.5g" % 10 ** rng.uniform(-2, 2)))
        ws = [float("%.5g" % x) for x in 10 ** rng.uniform(-2, 2, 3)]
        det, okc, drop = certify3(th, ws)
        cert.append({"theta5": th, "ws": ws, "certified": okc, "det": str(det), "dropped_row": drop})
    out["N3_certificates"] = {"n": len(cert), "certified": sum(c["certified"] for c in cert),
                              "failures": [c for c in cert if not c["certified"]][:5]}
    print("5. N = 3: rank 5 certified (ball arithmetic) at %d/%d random points" %
          (out["N3_certificates"]["certified"], len(cert)))
    # 6. adversarial (eq metric), hard constraint |d| >= 1, box beta in [1e-2,1e2], times, R in [1e-3, 1e3]
    def obj(x, sign):
        lb, ltl, d, lR = x
        if not (-2 <= lb <= 2 and -3 <= ltl <= 3 and 1.0 <= d <= 9 and -3 <= lR <= 3):
            return 10.0
        tl = 10 ** ltl; tR = tl * np.exp(sign * d)
        if not (1e-3 <= tR <= 1e3):
            return 10.0
        S, w, wts = svals5(10 ** lb, tR, tl, 10 ** lR)
        return float(np.log10(max(S[-1], 1e-300)))
    adv = []
    for k in range(24):
        sign = 1.0 if k % 2 else -1.0
        while True:
            x0 = np.array([rng.uniform(-2, 2), rng.uniform(-3, 3), rng.uniform(1, 6), rng.uniform(-3, 3)])
            if obj(x0, sign) < 5:
                break
        r = minimize(obj, x0, args=(sign,), method="Nelder-Mead",
                     options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 800})
        lb, ltl, d, lR = r.x
        tl = 10 ** ltl; tR = tl * np.exp(sign * d)
        at_box = (abs(abs(lb) - 2) < 0.05 or abs(abs(ltl) - 3) < 0.05 or abs(abs(np.log10(tR)) - 3) < 0.05
                  or abs(abs(lR) - 3) < 0.05)
        cls = "resonance" if d < 1.02 else ("box" if at_box else "interior")
        adv.append({"beta": 10 ** lb, "tR": tR, "tl": tl, "R": 10 ** lR, "log10_s5": float(r.fun), "class": cls})
    adv.sort(key=lambda z: z["log10_s5"])
    out["adversarial"] = adv
    for z in adv[:5]:
        print("6.", z)
    print("6. classes:", {c: [z["class"] for z in adv].count(c) for c in ("interior", "box", "resonance")})
    inter = [z for z in adv if z["class"] == "interior"]
    if inter:
        print("6. best interior:", inter[0])
    out["seconds"] = time.time() - t0
    with open(os.path.join(HERE, "data", "b4_resistance.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("seconds", out["seconds"])
