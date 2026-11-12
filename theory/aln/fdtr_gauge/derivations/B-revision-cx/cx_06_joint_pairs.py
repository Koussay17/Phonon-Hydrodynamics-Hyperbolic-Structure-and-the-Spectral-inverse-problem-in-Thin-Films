"""cx_06: joint search for twin PAIRS (th, th') with equal data at N frequencies, N = 3..7.

Pairs with equal N-frequency data form, by dimension count, a set of dimension 12 - 2N (if nonempty).
We fix 12 - 2N coordinates of th (random values) and solve the square system of 2N real equations in
the remaining 2N unknowns (2N - 6 coordinates of th, and th') by batched Levenberg-Marquardt in log
variables from starts with th' far from th.  N = 6: all 12 coordinates free.  N = 7: overdetermined
least squares (14 equations, 12 unknowns): any exact pair would contradict the 2r+1 heuristic.
Roots are refined at 60 digits and the square Jacobian is checked for regularity.

Usage: python cx_06_joint_pairs.py <conv: entropic|local|kin> <N> [nbatches]
"""
import json
import sys
import time
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
conv = sys.argv[1]
NF = int(sys.argv[2])
NB = int(sys.argv[3]) if len(sys.argv) > 3 else 10
BOX = sys.argv[4] if len(sys.argv) > 4 else "wide"
LOB, HIB = (-10.0, 5.0) if BOX == "wide" else (np.log(0.05), np.log(10.0))
EPS = {"entropic": 1.0, "local": 0.0, "kin": 0.084}[conv]
sp = cx.Spec("stack", "tied", EPS, EPS, xiA=0.2, bA=3.0)
OMS = {3: [0.5, 4.0, 1.7], 4: [0.5, 4.0, 1.7, 0.2], 5: [0.5, 4.0, 1.7, 0.2, 9.0],
       6: [0.5, 4.0, 1.7, 0.2, 9.0, 0.05], 7: [0.5, 4.0, 1.7, 0.2, 9.0, 0.05, 25.0]}
om = np.array(OMS[NF])
nfix = max(0, 12 - 2 * NF)
rng = np.random.default_rng(100 + NF + int(10 * EPS))


def ev(Y):
    """Y (m, 12) log variables -> residual (m, 2NF), Jacobian (m, 12, 2NF)."""
    T1, T2 = np.exp(Y[:, :6]), np.exp(Y[:, 6:])
    Z1, J1 = cx.ZJ_np_batch(sp, T1, om)
    Z2, J2 = cx.ZJ_np_batch(sp, T2, om)
    F = Z1 / Z2 - 1
    G1 = J1 / Z2[:, None, :] * T1[:, :, None]
    G2 = -(Z1 / Z2 ** 2)[:, None, :] * J2 * T2[:, :, None]
    G = np.concatenate([G1, G2], axis=1)
    return np.concatenate([F.real, F.imag], axis=1), np.concatenate([G.real, G.imag], axis=2)


def lm(Y, free, iters=300):
    m = Y.shape[0]
    lam = np.full(m, 1e-2)
    with np.errstate(all="ignore"):
        F, G = ev(Y)
        G = G[:, free, :]
        cost = np.sum(F * F, axis=1)
        nfree = int(np.sum(free))
        for it in range(iters):
            JTJ = np.einsum("mik,mjk->mij", G, G)
            JTF = np.einsum("mik,mk->mi", G, F)
            dia = np.einsum("mii->mi", JTJ)
            A = JTJ + lam[:, None, None] * (np.eye(nfree)[None] * (dia[:, :, None] + 1e-12))
            try:
                d = -np.linalg.solve(A, JTF[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                break
            d = np.clip(d, -1.5, 1.5)
            Yn = Y.copy()
            Yn[:, free] = np.clip(Y[:, free] + d, LOB, HIB)
            Fn, Gn = ev(Yn)
            cn = np.sum(Fn * Fn, axis=1)
            cn[~np.isfinite(cn)] = np.inf
            acc = cn < cost
            Y[acc] = Yn[acc]
            F[acc] = Fn[acc]
            G[acc] = Gn[acc][:, free, :]
            cost[acc] = cn[acc]
            lam[acc] = np.maximum(lam[acc] / 3, 1e-12)
            lam[~acc] = np.minimum(lam[~acc] * 4, 1e8)
    return Y, np.sqrt(cost)


def refine_pair_mp(y, free, dps=60, maxit=40):
    mp.mp.dps = dps
    th = [mp.mpf(float(v)) for v in np.exp(y)]
    idx = [i for i in range(12) if free[i]]
    for it in range(maxit):
        Fr, Jr = [], []
        rows_re, rows_im, fre, fim = [], [], [], []
        for w in om:
            w = mp.mpf(repr(float(w)))
            z1, d1 = cx.ZJ_mp(sp, th[:6], w)
            z2, d2 = cx.ZJ_mp(sp, th[6:], w)
            f = z1 / z2 - 1
            g = [x / z2 for x in d1] + [-z1 / z2 ** 2 * x for x in d2]
            fre.append(f.real); fim.append(f.imag)
            rows_re.append([g[i].real for i in idx]); rows_im.append([g[i].imag for i in idx])
        Fv = fre + fim
        Jm = mp.matrix(rows_re + rows_im)
        nF = mp.sqrt(mp.fsum(x * x for x in Fv))
        if nF < mp.mpf(10) ** (-dps + 8):
            break
        try:
            if Jm.rows == Jm.cols:
                dx = mp.lu_solve(Jm, mp.matrix([-x for x in Fv]))
            else:
                dx = mp.qr_solve(Jm, mp.matrix([-x for x in Fv]))[0]
        except ZeroDivisionError:
            # singular square Jacobian: least-squares step via SVD pseudo-inverse
            U_, S_, V_ = mp.svd_r(Jm)
            smax = max(S_[i] for i in range(len(S_)))
            y_ = [mp.fsum(U_[r, i] * (-Fv[r]) for r in range(Jm.rows)) for i in range(len(S_))]
            z_ = [y_[i] / S_[i] if S_[i] > smax * mp.mpf(10) ** (-dps + 10) else mp.mpf(0) for i in range(len(S_))]
            dx = [mp.fsum(V_[i, k] * z_[i] for i in range(len(S_))) for k in range(Jm.cols)]
        for k, i in enumerate(idx):
            th[i] += dx[k]
    S = mp.svd_r(Jm, compute_uv=False)
    sv = sorted([S[i] for i in range(len(S))], reverse=True)
    return th, nF, sv[-1] / sv[0]


t0 = time.time()
pairs = []
nb = 1000
for batch in range(NB):
    # fixed coordinates of th (target) and free mask
    Y = np.empty((nb, 12))
    base = np.log(np.array([0.8, 2.0, 0.5, 0.1, 0.3, 0.2]))
    Y[:, :6] = base[None, :] + rng.uniform(-1.5, 1.5, (nb, 6))
    Y[:, 6:] = Y[:, :6] + rng.choice([-1, 1], (nb, 6)) * rng.uniform(0.2, 1.5, (nb, 6))
    free = np.ones(12, bool)
    if nfix > 0:
        fixed_idx = rng.choice(6, nfix, replace=False)
        free[fixed_idx] = False
    Y, res = lm(Y, free)
    if NF < 7:
        ok = np.isfinite(res) & (res < 1e-11)
    else:
        ok = np.isfinite(res) & (res < 1e-9)
    _, Gok = ev(Y[ok])
    svs = np.linalg.svd(Gok[:, free, :], compute_uv=False) if np.any(ok) else np.zeros((0, 1))
    regular = (svs[:, -1] / svs[:, 0] > 1e-9) if np.any(ok) else np.zeros(0, bool)
    for y, reg in zip(Y[ok], regular):
        sep = np.max(np.abs(y[:6] - y[6:]))
        if sep < 1e-4 or not reg:
            continue
        if BOX != "wide" and (np.any(y[free] <= LOB + 1e-6) or np.any(y[free] >= HIB - 1e-6)):
            continue
        if any(np.max(np.abs(y - p["y"])) < 1e-6 for p in pairs):
            continue
        pairs.append(dict(y=y, free=free.copy(), res=float(np.linalg.norm(ev(y[None, :])[0]))))
    print(f"batch {batch}: converged {int(np.sum(ok))}, distinct separated pairs so far {len(pairs)}, "
          f"min residual {np.nanmin(res):.2e}, {time.time() - t0:.0f}s", flush=True)
    if NF == 7:
        print("   N=7 best residuals among separated:", np.sort(res[np.max(np.abs(Y[:, :6] - Y[:, 6:]), axis=1) > 1e-3])[:5])

# drop pairs on the resonance continuum (theta on Delta and the same (xi, b, R1, R2)): the known function-level fibre
def _on_continuum(y):
    return abs(y[2] - y[3]) < 1e-6 and abs(y[8] - y[9]) < 1e-6 and max(abs(y[i] - y[6 + i]) for i in (0, 1, 4, 5)) < 1e-6
n_cont = sum(1 for p in pairs if _on_continuum(p["y"]))
pairs = [p for p in pairs if not _on_continuum(p["y"])]
print(f"pairs on the resonance continuum removed: {n_cont}; remaining {len(pairs)}", flush=True)
out = []
for p in pairs[:40]:
    th, nF, cnd = refine_pair_mp(p["y"], p["free"])
    sep = max(abs(mp.log(th[i] / th[6 + i])) for i in range(6))
    # mismatch at a frequency not in the data
    w = mp.mpf("1.3")
    mis = abs(cx.Z_mp(sp, th[:6], w) / cx.Z_mp(sp, th[6:], w) - 1)
    rec = dict(theta=[mp.nstr(x, 30) for x in th[:6]], theta_twin=[mp.nstr(x, 30) for x in th[6:]],
               fixed=[int(i) for i in np.where(~p["free"])[0]], residual_60=float(nF), smin_over_smax=float(cnd),
               max_log_sep=float(sep), admissible=bool(min(th) > 0), mismatch_at_1p3=float(mis))
    out.append(rec)
    print("pair res60 %.1e cond %.1e sep %.3f mis(1.3) %.2e" % (float(nF), float(cnd), float(sep), float(mis)),
          [mp.nstr(x, 8) for x in th[:6]], "<->", [mp.nstr(x, 8) for x in th[6:]], flush=True)
with open(DIR + f"cx_06_joint_{conv}_N{NF}" + ("" if BOX == "wide" else "_moderate") + ".json", "w") as f:
    json.dump(dict(conv=conv, eps=EPS, omega=list(om), nfixed=nfix, n_pairs_double=len(pairs), n_continuum=n_cont, pairs=out), f, indent=1)
