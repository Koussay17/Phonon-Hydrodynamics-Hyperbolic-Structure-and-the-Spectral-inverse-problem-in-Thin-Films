"""cx_09b: Krawczyk certification of twin PAIRS with N frequencies (square system).

Unknowns: the free coordinates of th (those not fixed in cx_06) and all of th'; equations
F = [Re; Im](Z(omega_j; th)/Z(omega_j; th') - 1), j = 1..N; the fixed coordinates of th are set to their decimal values
(rounded to 12 digits, exact as decimals).  Newton at 110 digits, then K(X) inside int(X) with python-flint balls.
Usage: python cx_09b_certify_pairs.py <cx_06 json> <npairs> <out json>
"""
import json
import sys
sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle
throttle.wait_slot(2)          # light wait before any numerical import
import flint
from flint import arb, acb, arb_mat
flint.ctx.prec = 400
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

src = json.load(open(sys.argv[1]))
NP = int(sys.argv[2])
eps = src["eps"]
om_s = [repr(float(w)) for w in src["omega"]]
OPS = cx.flint_ops()


def F_mp(spm, th_full):
    rows_re, rows_im, Fr, Fi = [], [], [], []
    for w in om_s:
        w = mp.mpf(w)
        z1, d1 = cx.ZJ_mp(spm, th_full[:6], w)
        z2, d2 = cx.ZJ_mp(spm, th_full[6:], w)
        f = z1 / z2 - 1
        g = [x / z2 for x in d1] + [-z1 / z2 ** 2 * x for x in d2]
        Fr.append(f.real); Fi.append(f.imag)
        rows_re.append([x.real for x in g]); rows_im.append([x.imag for x in g])
    return Fr + Fi, rows_re + rows_im


def certify(pair):
    fixed = pair["fixed"]
    free = [i for i in range(12) if i not in fixed]
    mp.mp.dps = 110
    spm = cx.Spec("stack", "tied", mp.mpf(eps), mp.mpf(eps), xiA=mp.mpf("0.2"), bA=mp.mpf(3))
    th = [mp.mpf(x) for x in pair["theta"]] + [mp.mpf(x) for x in pair["theta_twin"]]
    for i in fixed:                       # exact decimal values of the fixed target coordinates
        th[i] = mp.mpf(mp.nstr(th[i], 12))
    for it in range(60):
        Fv, J = F_mp(spm, th)
        Jm = mp.matrix([[J[r][i] for i in free] for r in range(len(J))])
        nF = mp.sqrt(mp.fsum(x * x for x in Fv))
        if nF < mp.mpf(10) ** (-100):
            break
        dx = mp.lu_solve(Jm, mp.matrix([-x for x in Fv]))
        for k, i in enumerate(free):
            th[i] += dx[k]
    mid = [mp.nstr(x, 100) for x in th]
    # Krawczyk in flint
    flint.ctx.prec = 400
    spf = cx.Spec("stack", "tied", arb(repr(float(eps))) if eps not in (0.0, 1.0) else arb(int(eps)),
                  arb(repr(float(eps))) if eps not in (0.0, 1.0) else arb(int(eps)), xiA=arb("0.2"), bA=arb(3))
    rad = "1e-40"

    def FJ_flint(vals, with_rad):
        X = [acb(arb(v, rad)) if (with_rad and i in free) else acb(arb(v)) for i, v in enumerate(vals)]
        Fr, Fi, Jr, Ji = [], [], [], []
        for w in om_s:
            s = acb(0, arb(w)).sqrt()
            z1 = cx.response(spf, cx.seed(X[:6], acb(1), acb(0)), s, OPS)
            z2 = cx.response(spf, cx.seed(X[6:], acb(1), acb(0)), s, OPS)
            f = z1.v / z2.v - 1
            g = [d / z2.v for d in z1.d] + [-(z1.v / (z2.v * z2.v)) * d for d in z2.d]
            Fr.append(f.real); Fi.append(f.imag)
            Jr.append([g[i].real for i in free]); Ji.append([g[i].imag for i in free])
        return Fr + Fi, Jr + Ji
    Fm, _ = FJ_flint(mid, False)
    _, JX = FJ_flint(mid, True)
    n = len(free)
    JXm = arb_mat(JX)
    Jmid = arb_mat([[arb(JXm[i, k].mid()) for k in range(n)] for i in range(n)])
    Y = Jmid.inv()
    Y = arb_mat([[arb(Y[i, k].mid()) for k in range(n)] for i in range(n)])
    I = arb_mat([[arb(1) if i == k else arb(0) for k in range(n)] for i in range(n)])
    m = arb_mat([[arb(mid[i])] for i in free])
    Kv = m - Y * arb_mat([[x] for x in Fm]) + (I - Y * JXm) * arb_mat([[arb(0, rad)] for _ in range(n)])
    r = arb(rad)
    devs = [abs(Kv[k, 0] - m[k, 0]).upper() for k in range(n)]
    ok = all(d < r for d in devs)
    sep = max(abs(mp.log(th[i] / th[6 + i])) for i in range(6))
    w = mp.mpf("1.3")
    mis = abs(cx.Z_mp(spm, th[:6], w) / cx.Z_mp(spm, th[6:], w) - 1)
    return dict(certified=bool(ok), K_dev_over_radius=float(max(devs) / r.mid()), newton_residual=float(nF),
                theta=[mp.nstr(x, 30) for x in th[:6]], theta_twin=[mp.nstr(x, 30) for x in th[6:]],
                fixed_target_coordinates=fixed, max_log_sep=float(sep), mismatch_at_1p3=float(mis))


res = []
pairs = [p for p in src["pairs"] if p["smin_over_smax"] > 1e-20 and p["admissible"]][:NP]
for p in pairs:
    r = certify(p)
    res.append(r)
    print("certified" if r["certified"] else "NOT certified", "K-dev/r %.1e" % r["K_dev_over_radius"],
          "sep %.3f mis(1.3) %.2e" % (r["max_log_sep"], r["mismatch_at_1p3"]), [x[:10] for x in r["theta"]], "<->",
          [x[:10] for x in r["theta_twin"]], flush=True)
json.dump(dict(source=sys.argv[1], omega=src["omega"], eps=eps, results=res), open(sys.argv[3], "w"), indent=1)
