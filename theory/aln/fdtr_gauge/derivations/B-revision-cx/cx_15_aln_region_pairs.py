"""cx_15: twin pairs for targets in the AlN note-16 neighbourhood.

Target th = (1, b_AlN, tR, tl, R1, R2) with xi = 1 and b = b_AlN/b_s fixed (materials known), (tR, tl) in
[0.25, 4] x the note-16 values, (R1, R2) in [0.3, 3] x the note-16 values; twin th' free (positive); Au transducer and
N frequencies in the FDTR band (10 kHz - 200 MHz).  Underdetermined system (2N equations, 4 + 6 unknowns) solved by
batched LM in log variables from starts with th' far from th; a converged regular pair with th inside the box is a
target of the AlN neighbourhood with an N-frequency twin.  Refined at 60 digits.
Usage: python cx_15_aln_region_pairs.py <conv> <N> <nbatches>
"""
import sys as _s
_s.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle as _t
_t.wait_slot(2)          # light wait before any numerical import
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
NB = int(sys.argv[3]) if len(sys.argv) > 3 else 6
eps = 1.0 if conv == "entropic" else 0.0
th0, (xiA, bA), t0 = cx.aln_stack_nd(tau_l=10e-12)
sp = cx.Spec("stack", "tied", eps, eps, xiA=xiA, bA=bA)
fMHz = {3: [1.0, 20.0, 150.0], 4: [0.3, 3.0, 30.0, 150.0], 5: [0.1, 1.0, 10.0, 60.0, 180.0]}[NF]
om = 2 * np.pi * np.array(fMHz) * 1e6 * t0
lt0 = np.log(np.array(th0))
LOT = lt0.copy(); HIT = lt0.copy()
LOT[2:4] += np.log(0.25); HIT[2:4] += np.log(4.0)
LOT[4:6] += np.log(0.3); HIT[4:6] += np.log(3.0)
LOW, HIW = np.log(1e-4), np.log(1e3)
free = np.array([False, False, True, True, True, True] + [True] * 6)
rng = np.random.default_rng(2024 + NF + int(eps))


def ev(Y):
    T1, T2 = np.exp(Y[:, :6]), np.exp(Y[:, 6:])
    Z1, J1 = cx.ZJ_np_batch(sp, T1, om)
    Z2, J2 = cx.ZJ_np_batch(sp, T2, om)
    F = Z1 / Z2 - 1
    G1 = J1 / Z2[:, None, :] * T1[:, :, None]
    G2 = -(Z1 / Z2 ** 2)[:, None, :] * J2 * T2[:, :, None]
    G = np.concatenate([G1, G2], axis=1)
    return np.concatenate([F.real, F.imag], axis=1), np.concatenate([G.real, G.imag], axis=2)


def clip(Y):
    Y = Y.copy()
    Y[:, :6] = np.clip(Y[:, :6], LOT, HIT)
    Y[:, 6:] = np.clip(Y[:, 6:], LOW, HIW)
    Y[:, 0] = lt0[0]; Y[:, 1] = lt0[1]
    return Y


t1 = time.time()
found = []
for batch in range(NB):
    m = 1000
    Y = np.empty((m, 12))
    Y[:, :6] = LOT + (HIT - LOT) * rng.uniform(0, 1, (m, 6))
    Y[:, 6:] = Y[:, :6] + rng.choice([-1, 1], (m, 6)) * rng.uniform(0.3, 2.5, (m, 6))
    Y = clip(Y)
    lam = np.full(m, 1e-2)
    with np.errstate(all="ignore"):
        F, G = ev(Y); G = G[:, free, :]
        cost = np.sum(F * F, axis=1)
        for it in range(300):
            JTJ = np.einsum("mik,mjk->mij", G, G); JTF = np.einsum("mik,mk->mi", G, F)
            dia = np.einsum("mii->mi", JTJ)
            d = -np.linalg.solve(JTJ + lam[:, None, None] * (np.eye(10)[None] * (dia[:, :, None] + 1e-12)), JTF[:, :, None])[:, :, 0]
            Yn = Y.copy(); Yn[:, free] = Y[:, free] + np.clip(d, -1.0, 1.0); Yn = clip(Yn)
            Fn, Gn = ev(Yn); cn = np.sum(Fn * Fn, axis=1); cn[~np.isfinite(cn)] = np.inf
            acc = cn < cost
            Y[acc] = Yn[acc]; F[acc] = Fn[acc]; G[acc] = Gn[acc][:, free, :]; cost[acc] = cn[acc]
            lam[acc] = np.maximum(lam[acc] / 3, 1e-12); lam[~acc] = np.minimum(lam[~acc] * 4, 1e8)
    res = np.sqrt(cost)
    ok = np.isfinite(res) & (res < 1e-12)
    for y in Y[ok]:
        sep = np.max(np.abs(y[:6] - y[6:]))
        inside = np.all(y[:6] > LOT + 1e-6) and np.all(y[:6] < HIT - 1e-6) or True
        on_bound_twin = np.any(y[6:] <= LOW + 1e-6) or np.any(y[6:] >= HIW - 1e-6)
        if sep < 1e-3 or on_bound_twin:
            continue
        if any(np.max(np.abs(y - q)) < 1e-6 for q in found):
            continue
        found.append(y)
    print(f"batch {batch}: converged {int(ok.sum())}, distinct separated pairs {len(found)}, min res {np.nanmin(res):.1e} "
          f"[{time.time() - t1:.0f}s]", flush=True)

# refine (square system: fix 2 more target coordinates, the times of the target, at their found values)
mp.mp.dps = 60
out = []
for y in found[:25]:
    th = [mp.mpf(float(v)) for v in np.exp(y)]
    idx = [4, 5] + list(range(6, 12)) if NF == 4 else ([2, 3, 4, 5] + list(range(6, 12)) if NF == 5 else list(range(6, 12)))
    for it in range(40):
        rows_re, rows_im, fre, fim = [], [], [], []
        for w in om:
            w = mp.mpf(repr(float(w)))
            z1, d1 = cx.ZJ_mp(sp, th[:6], w); z2, d2 = cx.ZJ_mp(sp, th[6:], w)
            f = z1 / z2 - 1; g = [x / z2 for x in d1] + [-z1 / z2 ** 2 * x for x in d2]
            fre.append(f.real); fim.append(f.imag)
            rows_re.append([g[i].real for i in idx]); rows_im.append([g[i].imag for i in idx])
        Fv = fre + fim
        nF = mp.sqrt(mp.fsum(x * x for x in Fv))
        Jm = mp.matrix(rows_re + rows_im)
        if nF < mp.mpf("1e-52"):
            break
        try:
            dx = mp.lu_solve(Jm, mp.matrix([-x for x in Fv])) if Jm.rows == Jm.cols else mp.qr_solve(Jm, mp.matrix([-x for x in Fv]))[0]
        except ZeroDivisionError:
            break
        for k, i in enumerate(idx):
            th[i] += dx[k]
    S = mp.svd_r(Jm, compute_uv=False); sv = sorted([S[i] for i in range(len(S))], reverse=True)
    adm = min(th) > 0
    ratio = [float(th[i] / mp.mpf(float(th0[i]))) for i in range(6)]
    w = mp.mpf(repr(float(om[0] * 3.1)))
    mis = abs(cx.Z_mp(sp, th[:6], w) / cx.Z_mp(sp, th[6:], w) - 1)
    rec = dict(theta=[mp.nstr(x, 25) for x in th[:6]], theta_twin=[mp.nstr(x, 25) for x in th[6:]], residual_60=float(nF),
               smin_over_smax=float(sv[-1] / sv[0]), admissible=bool(adm), target_over_note16=ratio,
               tau_ps_target=[float(th[2] * t0 * 1e12), float(th[3] * t0 * 1e12)],
               tau_ps_twin=[float(th[8] * t0 * 1e12), float(th[9] * t0 * 1e12)], mismatch_off_data=float(mis))
    out.append(rec)
    print("pair res %.1e cond %.1e adm %s target/note16 %s twin tau (ps) %s  mis %.1e" %
          (rec["residual_60"], rec["smin_over_smax"], adm, np.round(ratio, 3).tolist(), np.round(rec["tau_ps_twin"], 2).tolist(),
           rec["mismatch_off_data"]), flush=True)
json.dump(dict(conv=conv, N=NF, f_MHz=fMHz, omega=list(om), note16=list(th0), transducer=(xiA, bA), pairs=out),
          open(DIR + f"cx_15_aln_region_{conv}_N{NF}.json", "w"), indent=1)
