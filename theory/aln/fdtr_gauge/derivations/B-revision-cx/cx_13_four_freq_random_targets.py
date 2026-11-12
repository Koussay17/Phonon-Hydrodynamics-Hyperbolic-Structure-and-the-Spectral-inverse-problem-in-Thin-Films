"""cx_13: Conjecture T4'' (four frequencies, generic targets).

(1) Random targets: for each of NT random targets (log-uniform box), batched LM multistart for exact roots of the
    8 real equations Z(omega_j; th') = Z(omega_j; th*), j = 1..4, th' != th*.  An exact root (relative residual
    < 1e-11, refined to < 1e-45 at 60 digits) for a random target would contradict T4''.
(2) Exceptional targets: for twin pairs found by cx_06 at N = 4, the rank of the 8 x 12 Jacobian of the pair
    equations (full rank 8 => the pair set is locally a 4-manifold) and the rank of the projection of its tangent
    space onto the target coordinates (<= 4 => the targets with twins form a set of dimension <= 4 locally).
Usage: python cx_13_four_freq_random_targets.py <conv> <NT> <starts per target>
"""
import sys as _s
_s.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import throttle as _t
_t.wait_slot(2)          # light wait before any numerical import
import json
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx")
import cxmodel as cx
import cxtwin as ct

DIR = r"D:/ResearchLab/orchestration/campaigns/20261111-fdtr-gauge-memory/derivations/B-revision-cx/"
conv = sys.argv[1]
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 30
NS = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
eps = 1.0 if conv == "entropic" else 0.0
sp = cx.Spec("stack", "tied", eps, eps, xiA=0.2, bA=3.0)
om = np.array([0.5, 4.0, 1.7, 0.2])
rng = np.random.default_rng(404)
base = np.log(np.array([0.8, 2.0, 0.5, 0.1, 0.3, 0.2]))
out = dict(conv=conv, omega=list(om), random_targets=[], exceptional=[])
best_overall = np.inf
for t in range(NT):
    tgt = np.exp(base + rng.uniform(-1.5, 1.5, 6))
    Zt = cx.Z_np(sp, tgt, om)
    lt = np.log(tgt)
    Y0 = lt[None, :] + rng.uniform(-3, 3, (NS, 6))
    Y, res = ct.batch_lm(sp, Y0, om, Zt, iters=250)
    far = np.max(np.abs(Y - lt[None, :]), axis=1) > 1e-4
    bf = float(np.min(res[far])) if np.any(far) else np.inf
    best_overall = min(best_overall, bf)
    exact = [y for y, r, f in zip(Y, res, far) if f and r < 1e-11]
    out["random_targets"].append(dict(target=tgt.tolist(), best_residual_distinct=bf, n_exact=len(exact),
                                      exact=[np.exp(y).tolist() for y in exact[:5]]))
    print(f"target {t}: best relative residual over distinct roots {bf:.2e}, exact roots {len(exact)}", flush=True)
out["best_residual_over_all_random_targets"] = best_overall

# (2) local structure of the exceptional set
try:
    pairs = json.load(open(DIR + f"cx_06_joint_{conv}_N4.json"))["pairs"]
except FileNotFoundError:
    pairs = []
mp.mp.dps = 50
for p in pairs[:10]:
    th = [mp.mpf(x) for x in p["theta"]] + [mp.mpf(x) for x in p["theta_twin"]]
    if p["smin_over_smax"] < 1e-20:
        continue
    rows_re, rows_im = [], []
    for w in om:
        w = mp.mpf(repr(float(w)))
        z1, d1 = cx.ZJ_mp(sp, th[:6], w)
        z2, d2 = cx.ZJ_mp(sp, th[6:], w)
        g = [x / z2 for x in d1] + [-z1 / z2 ** 2 * x for x in d2]
        rows_re.append([x.real for x in g]); rows_im.append([x.imag for x in g])
    Jm = mp.matrix(rows_re + rows_im)                     # 8 x 12
    S_ = mp.svd_r(Jm, compute_uv=False)
    sv = [S_[i] for i in range(len(S_))]
    rank = sum(1 for x in sv if x > mp.mpf("1e-30") * max(sv))
    # tangent space = null space of Jm: eigenvectors of Jm^T Jm with (numerically) zero eigenvalues
    E_, Q_ = mp.eigsy(Jm.T * Jm)
    order = sorted(range(12), key=lambda i: E_[i])
    null_idx = order[:12 - rank]
    N = mp.matrix([[Q_[j, i] for j in range(12)] for i in null_idx])               # (12-rank) x 12
    P = mp.matrix([[N[i, j] for j in range(6)] for i in range(N.rows)])
    Sp = mp.svd_r(P, compute_uv=False)
    spv = sorted([Sp[i] for i in range(len(Sp))], reverse=True)
    prank = sum(1 for x in spv if x > mp.mpf("1e-30"))
    out["exceptional"].append(dict(theta=p["theta"][:3], jac_rank=rank, tangent_dim=12 - rank, projected_rank=prank,
                                   sv_min=float(min(sv) / max(sv))))
    print(f"exceptional pair: Jacobian rank {rank} of 8, tangent dim {12 - rank}, projection rank on targets {prank}", flush=True)
json.dump(out, open(DIR + f"cx_13_four_freq_{conv}.json", "w"), indent=1)
