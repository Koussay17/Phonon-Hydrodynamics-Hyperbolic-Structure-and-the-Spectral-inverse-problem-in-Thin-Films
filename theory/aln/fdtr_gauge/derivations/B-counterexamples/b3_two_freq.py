"""B3: rank of the 4x4 real Jacobian built from N = 2 frequencies.

M2(theta; w1, w2) = [Re J(w1); Im J(w1); Re J(w2); Im J(w2)],  J = d ln Z / d(ln xi1, ln b, ln tR, ln tl).
Units xi1 = bs = 1 at the point (frequencies in units of 1/xi1^2).

Structural facts used (derived in the report):
  rank M2 < 4 for ALL pairs  <=>  (case 1) full-band rank < 4, or
                                  (case 2) Re Phi(i w) is constant in w, where
                                  Phi(p) = (1 + tR p)(1 + tl p) q(p)/p,  q = (u + v)/(2(u - v)).
Tests
 1. random (theta, w1, w2): Hadamard ratio h = |det M2| / prod ||columns||, and s4/s1, at 50 digits;
 2. ball-arithmetic certificates (python-flint, 200 bits) that det M2 != 0 at random theta and pairs;
 3. case-2 defect: spread of Re Phi(i w) over w (zero would break B3 at that theta);
 4. explicit singular pair: a root of det M2 along w2 for fixed w1 (non-generic pairs exist);
 5. adversarial: minimise over theta (|ln tR/tl| >= c) the best-pair h over a frequency grid.
"""
import sys, os, json, time
import numpy as np
import mpmath as mp
from scipy.optimize import minimize
from flint import acb, arb, arb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm


def M2_mp(th, w1, w2):
    xi1, b, tR, tl = th
    rows = []
    for w in (w1, w2):
        Z, J = zm.jac_mp(mp.mpf(w), mp.mpf(xi1), mp.mpf(b), mp.mpf(tR), mp.mpf(tl), 1, 0, "code")
        c = [1, 1, mp.mpf(tR), mp.mpf(tl)]
        jr = [J[k] * c[k] / Z for k in range(4)]
        rows.append([x.real for x in jr])
        rows.append([x.imag for x in jr])
    return mp.matrix(rows)


def hadamard(M):
    det = mp.det(M)
    cn = 1
    for j in range(M.cols):
        cn *= mp.sqrt(mp.fsum(M[i, j] ** 2 for i in range(M.rows)))
    return abs(det) / cn, det


def M2_np(th, ws):
    xi1, b, tR, tl = th
    Z, J = zm.jac_np(np.asarray(ws, float), xi1, b, tR, tl, 1.0, 0.0)
    return J[:, :4] * np.array([1, 1, tR, tl]) / Z[:, None]     # (n, 4) complex rows


def best_pair_h(th, ws):
    """max over pairs (i < j) of the Hadamard ratio of [Re Ji; Im Ji; Re Jj; Im Jj]."""
    J = M2_np(th, ws)
    n = len(ws)
    best = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            M = np.array([J[i].real, J[i].imag, J[j].real, J[j].imag])
            h = abs(np.linalg.det(M)) / np.prod(np.linalg.norm(M, axis=0))
            best = max(best, h)
    return best


# ---- ball arithmetic -------------------------------------------------------------------------
def jac_acb(w, xi1, b, tR, tl, bs=1):
    p = acb(0, w)
    s = p.sqrt()
    r = ((1 + tR * p) / (1 + tl * p)).sqrt()
    A = xi1 * s * r
    B = b * s / r
    E = (-2 * A).exp()
    ome = -((-2 * A).expm1())
    y = B / (bs * s)
    rho = (y - 1) / (y + 1)
    D1 = 1 - rho * E
    Z = (1 + rho * E) / (B * D1)
    u = -4 * A * rho * E / (B * D1 ** 2)
    v = -ome * ((1 + E) * (1 + y ** 2) + 2 * y * ome) / (B * D1 ** 2 * (y + 1) ** 2)
    gR = p / (2 * (1 + tR * p))
    gl = p / (2 * (1 + tl * p))
    return Z, [u, v, gR * (u - v), -gl * (u - v)]


def certify(th, w1, w2, prec=200):
    ctx.prec = prec
    xi1, b, tR, tl = (arb(str(x)) for x in th)   # decimal strings: exact rationals near the floats
    rows = []
    for w in (w1, w2):
        Z, J = jac_acb(arb(str(w)), xi1, b, tR, tl)
        c = [arb(1), arb(1), tR, tl]
        jr = [J[k] * c[k] / Z for k in range(4)]
        rows.append([x.real for x in jr])
        rows.append([x.imag for x in jr])
    det = arb_mat(rows).det()
    return det, (not det.contains(0))


def rephi_spread(th, ws):
    xi1, b, tR, tl = th
    Z, J = zm.jac_np(np.asarray(ws, float), xi1, b, tR, tl)
    u, v = J[:, 0], J[:, 1]
    q = (u + v) / (2 * (u - v))
    p = 1j * np.asarray(ws)
    Phi = (1 + tR * p) * (1 + tl * p) * q / p
    re = Phi.real
    return float(np.std(re) / (np.std(Phi.imag) + abs(np.mean(re)) + 1e-300)), re


def rand_theta(rng, box=(2.0, 3.0)):
    beta = 10 ** rng.uniform(-box[0], box[0])
    tR = 10 ** rng.uniform(-box[1], box[1])
    tl = 10 ** rng.uniform(-box[1], box[1])
    return (1.0, beta, tR, tl)


if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(3)
    out = {}
    mp.mp.dps = 50
    # 1. random sampling
    H, S = [], []
    worst = None
    for k in range(1500):
        th = rand_theta(rng)
        w1, w2 = 10 ** rng.uniform(-3, 3, 2)
        M = M2_mp(th, w1, w2)
        h, det = hadamard(M)
        sv = mp.svd_r(M, compute_uv=False)
        sv = sorted([abs(x) for x in sv], reverse=True)
        H.append(float(h)); S.append(float(sv[-1] / sv[0]))
        if worst is None or h < worst[0]:
            worst = (float(h), th, float(w1), float(w2), float(sv[-1] / sv[0]))
    H = np.array(H); S = np.array(S)
    out["random"] = {"n": len(H), "h_min": float(H.min()), "h_median": float(np.median(H)),
                     "h_quantiles_1_10_50": np.quantile(H, [0.01, 0.1, 0.5]).tolist(),
                     "s4_over_s1_min": float(S.min()), "s4_over_s1_median": float(np.median(S)),
                     "worst": worst}
    print("1.", json.dumps(out["random"]))
    # 2. certificates
    cert = []
    for k in range(200):
        th = rand_theta(rng, box=(1.5, 2.0))
        w1, w2 = 10 ** rng.uniform(-2, 2, 2)
        th_r = tuple(float("%.6g" % x) for x in th)
        w1, w2 = float("%.6g" % w1), float("%.6g" % w2)
        det, ok = certify(th_r, w1, w2)
        cert.append({"theta": th_r, "w1": w1, "w2": w2, "det": str(det), "certified_nonzero": ok})
    nok = sum(c["certified_nonzero"] for c in cert)
    out["certificates"] = {"n": len(cert), "certified": nok, "examples": cert[:5],
                           "failures": [c for c in cert if not c["certified_nonzero"]]}
    print(f"2. ball-arithmetic det != 0 certified for {nok}/{len(cert)} random (theta, w1, w2)")
    # 3. case-2 defect
    ws = np.logspace(-3, 3, 400)
    spreads = []
    for k in range(300):
        th = rand_theta(rng)
        sp, re = rephi_spread(th, ws)
        spreads.append((sp, th))
    spreads.sort(key=lambda z: z[0])
    out["case2_spread"] = {"min": spreads[0][0], "argmin": spreads[0][1],
                           "median": float(np.median([z[0] for z in spreads]))}
    print("3. case-2 spread of Re Phi: min %.3e at %s, median %.3e" %
          (spreads[0][0], spreads[0][1], out["case2_spread"]["median"]))
    # 4. explicit singular pair: root of det along w2 at fixed w1
    th = (1.0, 2.0, 0.5, 0.1)
    w1 = 0.3
    f = lambda lw2: mp.det(M2_mp(th, w1, mp.exp(lw2)))
    grid = np.linspace(np.log(1e-3), np.log(1e3), 241)
    vals = [f(mp.mpf(g)) for g in grid]
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] * vals[i + 1] < 0:
            rt = mp.findroot(f, (mp.mpf(grid[i]), mp.mpf(grid[i + 1])), solver="anderson")
            roots.append(float(mp.exp(rt)))
    out["singular_pairs"] = {"theta": th, "w1": w1, "w2_roots": roots}
    print("4. theta", th, "w1", w1, "singular w2 (det M2 = 0):", roots)
    if roots:
        M = M2_mp(th, w1, roots[0])
        sv = sorted([abs(x) for x in mp.svd_r(M, compute_uv=False)], reverse=True)
        print("   s/s1 at the first root (50 digits):", [float(x / sv[0]) for x in sv])
        out["singular_pairs"]["svals_at_root"] = [float(x / sv[0]) for x in sv]
    # 5. adversarial over theta: best-pair Hadamard ratio on a 25-point grid
    wg = np.logspace(-3, 3, 25)
    def obj(x, c, sign):
        lb, ltl, d = x
        if not (-2 <= lb <= 2 and -3 <= ltl <= 3 and c <= d <= 8):
            return 10.0
        tl = 10 ** ltl; tR = tl * np.exp(sign * d)
        if not (1e-3 <= tR <= 1e3):
            return 10.0
        return np.log10(max(best_pair_h((1.0, 10 ** lb, tR, tl), wg), 1e-300))
    adv = []
    for k in range(24):
        c = 0.1 if k < 12 else 1.0
        sign = 1.0 if k % 2 else -1.0
        while True:
            x0 = np.array([rng.uniform(-2, 2), rng.uniform(-3, 3), rng.uniform(c, 5)])
            if obj(x0, c, sign) < 5:
                break
        r = minimize(obj, x0, args=(c, sign), method="Nelder-Mead",
                     options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 600})
        lb, ltl, d = r.x
        adv.append({"c": c, "beta": 10 ** lb, "tl": 10 ** ltl, "tR": 10 ** ltl * np.exp(sign * d),
                    "log10_best_pair_h": float(r.fun)})
    adv.sort(key=lambda z: z["log10_best_pair_h"])
    out["adversarial_best_pair"] = adv
    for z in adv[:6]:
        print("5.", z)
    out["seconds"] = time.time() - t0
    with open(os.path.join(HERE, "data", "b3_two_freq.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("seconds", out["seconds"])
