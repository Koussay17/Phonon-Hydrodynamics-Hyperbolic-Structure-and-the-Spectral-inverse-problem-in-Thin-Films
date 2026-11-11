"""Task 5. Interface resistance R between film and substrate: theta = (xi1, b, tau_R, tau_l, R),
Z_sub = R + 1/(bs sqrt p).

(R1) 50-digit SVD of the real (2N x 5) and complex (N x 5) Jacobians, N = 1, 2, 3, 4, 5, 8, at
     random rational points; derivatives (a) SymPy symbolic, (b) mpmath numerical differentiation
     of the quadrupole product film x [[1, R], [0, 1]] closed by 1/(bs sqrt p).
(R2) Arb certificates: real N = 2 -> some 4 x 4 minor != 0, real N = 3 -> some 5 x 5 minor != 0,
     complex N = 5 -> det != 0.
(R3) strata: tau_R = tau_l with R > 0; tau_R = tau_l with R = 0 (b != bs and b = bs);
     R = 0 off the diagonal; b = bs off the diagonal; certified lower bounds and the predicted
     null vectors (B-symbolic.md, section 6), checked at 50 digits.
(R4) grid search in the dimensionless variables (beta, T_R, T_l, Rhat = R bs/xi1), 65 frequencies.
Output: results/05_resistance.txt
"""
import itertools
import json
import os
import random
import sys
from fractions import Fraction as Fr
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "2"
import numpy as np
import mpmath as mp
from flint import acb, arb_mat, acb_mat, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


mp.mp.dps = 50
ctx.prec = 320
rng = random.Random(555)


def rnd_log(lo, hi, den=1000):
    uu = rng.uniform(lo, hi)
    return Fr(max(1, round(10 ** uu * den)), den)


jsym = M.sympy_jacobian_mp(with_R=True)

# ------------------------------------------------------------------ (R1)
pts = []
for k in range(16):
    xi1 = rnd_log(-1, 0.7)
    T = xi1 ** 2
    b = rnd_log(-1, 1)
    pt = dict(xi1=xi1, b=b, bs=rnd_log(-1, 1), tau_R=rnd_log(-2, 1.3) * T, tau_l=rnd_log(-2, 1.3) * T,
              R=rnd_log(-2, 1) * xi1 / b)
    if k >= 12:
        pt["tau_l"] = Fr(0)          # Cattaneo film with interface resistance
    pt["omegas"] = sorted(rnd_log(-2, 2) / T for _ in range(8))
    pts.append(pt)
SUB = {1: [3], 2: [1, 5], 3: [0, 3, 6], 4: [0, 2, 5, 7], 5: [0, 2, 3, 5, 7], 8: list(range(8))}
log("(R1) theta = (xi1, b, tau_R, tau_l, R); singular values s_i/s_1 of the equilibrated Jacobians, 50 digits")
summary = {}
worst_ab = mp.mpf(0)
for ip, pt in enumerate(pts):
    th = (pt["xi1"], pt["b"], pt["tau_R"], pt["tau_l"], pt["R"])
    rows = []
    for om in pt["omegas"]:
        pv = mp.mpc(0, M._mpf(om))
        ra = jsym(pv, th[0], th[1], th[2], th[3], pt["bs"], th[4])
        rb = M.jac_fd_mp(pv, th, pt["bs"], model="quad", dps_inner=100)
        sc = max(abs(x) for x in ra)
        worst_ab = max(worst_ab, max(abs(x - y) for x, y in zip(ra, rb)) / sc)
        rows.append(ra)
    out = []
    for N, sel in SUB.items():
        sr = M.singular_values(M.equilibrate(M.real_jacobian([rows[i] for i in sel])))
        sc_ = M.singular_values(M.equilibrate(M.complex_jacobian([rows[i] for i in sel])), complex_=True)
        rr = sum(1 for x in sr if x / sr[0] > mp.mpf("1e-30"))
        rc = sum(1 for x in sc_ if x / sc_[0] > mp.mpf("1e-30"))
        summary.setdefault(N, []).append((rr, rc, sr[-1] / sr[0], sc_[-1] / sc_[0]))
        out.append(f"N={N}: real {rr} (min {mp.nstr(sr[-1]/sr[0], 2)}), complex {rc} (min {mp.nstr(sc_[-1]/sc_[0], 2)})")
    log(f"  [{ip:02d}] tau_l={'0' if pt['tau_l'] == 0 else 'pos'}: " + "; ".join(out))
log(f"  max rel. difference symbolic vs numerical derivatives: {mp.nstr(worst_ab, 3)}")
for N in SUB:
    v = summary[N]
    log(f"  N={N}: real ranks {sorted(set(x[0] for x in v))}, complex ranks {sorted(set(x[1] for x in v))},"
        f" min s_min/s1 real {mp.nstr(min(x[2] for x in v), 3)}, complex {mp.nstr(min(x[3] for x in v), 3)}")

# ------------------------------------------------------------------ (R2)
log("")
log("(R2) Arb certificates (forward-mode AD in ball arithmetic)")


def best_minor(rows_, k, complex_=False):
    best = None
    for ri in itertools.combinations(range(len(rows_)), k):
        for ci in itertools.combinations(range(len(rows_[0])), k):
            sub = [[rows_[i][j] for j in ci] for i in ri]
            d = (acb_mat(sub) if complex_ else arb_mat(sub)).det()
            mag = abs(d).mid() if complex_ else abs(d.mid())
            if best is None or mag > best[0]:
                best = (mag, d, ci)
    return best[1], best[2]


def realrows(cr):
    out = []
    for r in cr:
        out.append([c.real for c in r])
        out.append([c.imag for c in r])
    return out


allok = True
for ip, pt in enumerate(pts):
    th = (pt["xi1"], pt["b"], pt["tau_R"], pt["tau_l"], pt["R"])
    cr = {i: M.jac_ball(pt["omegas"][i], th, pt["bs"]) for i in sorted(set(SUB[2] + SUB[3] + SUB[5]))}
    d2, _ = best_minor(realrows([cr[i] for i in SUB[2]]), 4)
    d3, _ = best_minor(realrows([cr[i] for i in SUB[3]]), 5)
    d5 = acb_mat([list(cr[i]) for i in SUB[5]]).det()
    ok = (not d2.contains(0)) and (not d3.contains(0)) and (not d5.contains(0))
    allok &= ok
    log(f"  [{ip:02d}] real N=2 4x4 minor {d2.mid().str(5, radius=False)} +/- {d2.rad().str(2, radius=False)};"
        f" real N=3 5x5 minor {d3.mid().str(5, radius=False)} +/- {d3.rad().str(2, radius=False)};"
        f" complex N=5 |det| {abs(d5).mid().str(5, radius=False)}  {'ok' if ok else 'NOT CERTIFIED'}")
log(f"  all certificates valid: {allok}")

# ------------------------------------------------------------------ (R3)
log("")
log("(R3) strata with R (50-digit SVD of the real Jacobian at 8 frequencies, natural scaling; Arb lower bounds)")


def strat(name, th, bs, r_exp, nulls):
    oms = [Fr(1, 10), Fr(1, 3), Fr(1), Fr(2), Fr(5), Fr(12), Fr(30), Fr(80)]
    oms = [o / (th[0] * th[0]) for o in oms]
    rows = [jsym(mp.mpc(0, M._mpf(o)), *th[:4], bs, th[4]) for o in oms]
    Jr = M.real_jacobian(rows)
    sv = M.singular_values(Jr)
    res = []
    for nv in nulls:
        vv = [M._mpf(x) for x in nv]
        r = max(abs(sum(rows[i][j] * vv[j] for j in range(5))) for i in range(len(rows)))
        sc = max(abs(x) for rr in rows for x in rr)
        res.append(mp.nstr(r / sc, 3))
    cr = [M.jac_ball(o, th, bs) for o in oms[:4]]
    d, ci = best_minor(realrows(cr), r_exp)
    ok = not d.contains(0)
    log(f"  {name}: s/s1 = [{', '.join(mp.nstr(x / sv[0], 3) for x in sv)}]; certified {r_exp}x{r_exp} minor"
        f" {d.mid().str(4, radius=False)} +/- {d.rad().str(2, radius=False)} {'ok' if ok else 'NOT CERTIFIED'};"
        f" null-vector residuals {res}")
    return ok


ok = True
for k in range(2):
    x1, bb, bsv, tau = rnd_log(-0.5, 0.5), rnd_log(-0.5, 0.5), rnd_log(-0.5, 0.5), rnd_log(-1, 0.5)
    Rv = rnd_log(-1, 0) * x1 / bb
    ok &= strat("tau_R = tau_l, R > 0 (rank 4)", (x1, bb, tau, tau, Rv), bsv, 4, [(0, 0, 1, 1, 0)])
    c5 = -(1 - (bb / bsv) ** 2) / bb
    ok &= strat("tau_R = tau_l, R = 0, b != bs (rank 3)", (x1, bb, tau, tau, Fr(0)), bsv, 3,
                [(0, 0, 1, 1, 0), (1, 0, 0, 0, c5)])
    ok &= strat("tau_R = tau_l, R = 0, b = bs (rank 3)", (x1, bb, tau, tau, Fr(0)), bb, 3,
                [(0, 0, 1, 1, 0), (1, 0, 0, 0, 0)])
    ok &= strat("tau_R = tau_l = 0 (Fourier), R = 0 (rank 3)", (x1, bb, Fr(0), Fr(0), Fr(0)), bsv, 3,
                [(0, 0, 1, 1, 0), (1, 0, 0, 0, c5)])
    ok &= strat("tau_R = tau_l, R > 0, b = bs (rank 4)", (x1, bb, tau, tau, Rv), bb, 4, [(0, 0, 1, 1, 0)])
    ok &= strat("R = 0, tau_R != tau_l (rank 5)", (x1, bb, tau, tau / 3, Fr(0)), bsv, 5, [])
    ok &= strat("b = bs, R > 0, tau_R != tau_l (rank 5)", (x1, bb, tau, tau / 3, Rv), bb, 5, [])
log(f"  all stratum certificates valid: {ok}")

# ------------------------------------------------------------------ (R4)
log("")
log("(R4) grid search: beta in 10^[-2,2] (9), T_R, T_l in {0} U 10^[-3,3] (13), Rhat = R bs/xi1 in {0} U 10^[-3,2] (7);")
log("     65 frequencies P in [1e-4, 1e4]; real Jacobian 130 x 5, equilibrated columns (zero columns cannot occur here")
log("     except dZ/dxi1 on {T_R = T_l, R = 0, beta = 1}, which is excluded from the statistics below and covered by the certified strata of R3)")
def columns(P, beta, TR, TL, Rh):
    p = 1j * P
    rho = np.sqrt((1 + TR * p) / (1 + TL * p))
    sq = np.sqrt(p)
    A = sq * rho
    B = beta * sq / rho
    w = B * (Rh + 1 / sq)
    t = np.tanh(A)
    den = B * (1 + w * t) ** 2
    U = A * (1 - t ** 2) * (1 - w ** 2) / den
    V = -t * (1 + 2 * w * t + w ** 2) / den
    W = U - V
    Y = (1 - t ** 2) / (1 + w * t) ** 2
    return np.stack([U, V / beta, p / (2 * (1 + TR * p)) * W, -p / (2 * (1 + TL * p)) * W, Y], axis=-1)


betas = 10.0 ** np.linspace(-2, 2, 9)
Ts = np.concatenate([[0.0], 10.0 ** np.linspace(-3, 3, 13)])
Rhs = np.concatenate([[0.0], 10.0 ** np.linspace(-3, 2, 6)])
Ps = 10.0 ** np.linspace(-4, 4, 65)
BB, TR, TL, RH = np.meshgrid(betas, Ts, Ts, Rhs, indexing="ij")
with np.errstate(all="ignore"):
    J = columns(Ps, BB[..., None], TR[..., None], TL[..., None], RH[..., None])
    Jr = np.concatenate([J.real, J.imag], axis=-2)
    nrm = np.linalg.norm(Jr, axis=-2, keepdims=True)
    nrm[nrm == 0] = 1
    S = np.linalg.svd(Jr / nrm, compute_uv=False)
    S = S / S[..., :1]
diag = TR == TL
for name, m in [("off the diagonal, R > 0", (~diag) & (RH > 0)), ("off the diagonal, R = 0", (~diag) & (RH == 0))]:
    i = np.unravel_index(np.argmin(np.where(m, S[..., 4], np.inf)), BB.shape)
    log(f"  {name}: min s5/s1 = {S[..., 4][i]:.3e} at beta={BB[i]:.3g}, T_R={TR[i]:.3g}, T_l={TL[i]:.3g}, Rhat={RH[i]:.3g}")
m = diag & (RH > 0)
log(f"  on the diagonal, R > 0: max s5/s1 = {np.max(np.where(m, S[..., 4], 0)):.2e};"
    f" min s4/s1 = {np.min(np.where(m, S[..., 3], np.inf)):.3e}")
m = diag & (RH == 0) & ~np.isclose(BB, 1)
log(f"  on the diagonal, R = 0, beta != 1: max s4/s1 = {np.max(np.where(m, S[..., 3], 0)):.2e};"
    f" min s3/s1 = {np.min(np.where(m, S[..., 2], np.inf)):.3e}")
(HERE / "results" / "05_resistance.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
