"""Task 2, numerical search for degeneracies of the function-level Jacobian.

Scaling: Z(p; xi1, b, tau_R, tau_l; bs) = (xi1/bs) Zhat(p xi1^2; b/bs, tau_R/xi1^2, tau_l/xi1^2),
so the rank of the function-level Jacobian depends only on (beta, T_R, T_l) = (b/bs, tau_R/xi1^2,
tau_l/xi1^2) (B-symbolic.md, section 3.1).  We set xi1 = bs = 1 and scan a log grid
  beta in 10^[-2, 2] (17 values, contains beta = 1),  T_R, T_l in {0} U 10^[-3, 3] (25 values)
with 65 frequencies P = omega xi1^2 log-spaced in [1e-4, 1e4] (real Jacobian 130 x 4, columns
equilibrated).  Columns from the verified closed forms [U/xi1, V/b, gR (U-V), -gl (U-V)].
Reported: s4/s1 (rank-4 indicator) off and near the diagonal, its scaling with |T_R - T_l|,
s3/s1 on the diagonal, the boundary strata, and the N = 2 determinant (Hadamard ratio, best of
8 random frequency pairs) on the whole grid.  The lowest off-diagonal points are recomputed at
50 digits with mpmath.
Output: results/02_sigma_numeric.txt, results/02_grid.npz
"""
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "2"
import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdtr1d as M  # noqa: E402

lines = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    lines.append(s)


def columns(P, beta, TR, TL, xi1=1.0, bs=1.0, R=0.0):
    """Jacobian columns (complex) at p = iP; broadcasting over all inputs."""
    p = 1j * P
    rho = np.sqrt((1 + TR * p) / (1 + TL * p))
    sq = np.sqrt(p)
    A = xi1 * sq * rho
    B = beta * sq / rho
    Zs = R + 1 / (bs * sq)
    w = B * Zs
    t = np.tanh(A)
    den = B * (1 + w * t) ** 2
    U = A * (1 - t ** 2) * (1 - w ** 2) / den
    V = -t * (1 + 2 * w * t + w ** 2) / den
    W = U - V
    gR = p / (2 * (1 + TR * p))
    gL = p / (2 * (1 + TL * p))
    Y = (1 - t ** 2) / (1 + w * t) ** 2
    return np.stack([U / xi1, V / beta, gR * W, -gL * W, Y], axis=-1)


def svals(J):
    """J complex (..., nP, k) -> normalised singular values of the equilibrated real (..., 2nP, k)."""
    Jr = np.concatenate([J.real, J.imag], axis=-2)
    nrm = np.linalg.norm(Jr, axis=-2, keepdims=True)
    nrm[nrm == 0] = 1.0
    s = np.linalg.svd(Jr / nrm, compute_uv=False)
    return s / s[..., :1]


betas = 10.0 ** np.linspace(-2, 2, 17)
Ts = np.concatenate([[0.0], 10.0 ** np.linspace(-3, 3, 25)])
Ps = 10.0 ** np.linspace(-4, 4, 65)
BB, TR, TL = np.meshgrid(betas, Ts, Ts, indexing="ij")
J = columns(Ps[None, None, None, :], BB[..., None], TR[..., None], TL[..., None])[..., :4]
S = svals(J)                                   # (17, 26, 26, 4)
s4, s3 = S[..., 3], S[..., 2]
np.savez_compressed(HERE / "results" / "02_grid.npz", betas=betas, Ts=Ts, Ps=Ps, S=S)

diag = TR == TL
both_pos = (TR > 0) & (TL > 0)
with np.errstate(divide="ignore"):
    dlog = np.where(both_pos, np.abs(np.log10(np.where(both_pos, TR / np.where(TL > 0, TL, 1), 1))), np.inf)
log("Function-level Jacobian on the grid (xi1 = bs = 1), real 130 x 4, equilibrated columns.")
log(f"grid points: {S[..., 0].size};  diagonal points: {diag.sum()}")
log("")
log("1. Off the diagonal")
for d0 in (0.25, 0.5, 1.0, 2.0):
    m = (~diag) & (dlog >= d0 - 1e-12)
    i = np.unravel_index(np.argmin(np.where(m, s4, np.inf)), s4.shape)
    log(f"   |log10(T_R/T_l)| >= {d0}: min s4/s1 = {s4[i]:.3e} at beta={BB[i]:.3g}, T_R={TR[i]:.3g}, T_l={TL[i]:.3g}")
for name, m in [("T_l = 0 < T_R (Cattaneo)", (TL == 0) & (TR > 0)), ("T_R = 0 < T_l", (TR == 0) & (TL > 0))]:
    i = np.unravel_index(np.argmin(np.where(m, s4, np.inf)), s4.shape)
    log(f"   {name}: min s4/s1 = {s4[i]:.3e} at beta={BB[i]:.3g}, T_R={TR[i]:.3g}, T_l={TL[i]:.3g}")
m = (~diag) & (np.isclose(BB, 1.0))
i = np.unravel_index(np.argmin(np.where(m & (dlog >= 0.25), s4, np.inf)), s4.shape)
log(f"   beta = 1 (b = bs), |log10(T_R/T_l)| >= 0.25: min s4/s1 = {s4[i]:.3e} at T_R={TR[i]:.3g}, T_l={TL[i]:.3g}")
log("")
log("2. On the diagonal T_R = T_l (rank <= 3 by the identity col3 + col4 = 0)")
log(f"   max s4/s1 on the diagonal: {np.max(np.where(diag, s4, 0)):.2e}")
mb = diag & ~np.isclose(BB, 1.0)
i = np.unravel_index(np.argmin(np.where(mb, s3, np.inf)), s3.shape)
log(f"   beta != 1: min s3/s1 = {s3[i]:.3e} at beta={BB[i]:.3g}, T={TR[i]:.3g}")
mb1 = diag & np.isclose(BB, 1.0)
log(f"   beta = 1, equilibrated columns: max s3/s1 = {np.max(np.where(mb1, s3, 0)):.2e}  <- artefact: the xi1 column is")
log("      identically zero there (U = 0), and equilibration rescales its round-off (~1e-16) to unit norm.")
# natural scaling (no equilibration): [xi1 dZ/dxi1, b dZ/db, dZ/dT_R, dZ/dT_l]
Jn = J[np.isclose(BB, 1.0) & diag]                       # (n, nP, 4) complex, beta = 1 so b dZ/db = dZ/db
Jnr = np.concatenate([Jn.real, Jn.imag], axis=-2)
sn = np.linalg.svd(Jnr, compute_uv=False)
sn = sn / sn[:, :1]
log(f"   beta = 1, natural scaling (no equilibration): max s3/s1 = {sn[:, 2].max():.2e}, min s2/s1 = {sn[:, 1].min():.3e}  (rank 2)")
log("")
log("3. Scaling of s4/s1 with the distance to the diagonal, T_l = T (1 + eps), T_R = T")
for beta in (0.1, 1.0, 10.0):
    for T in (0.01, 1.0, 100.0):
        eps = 10.0 ** np.arange(-1, -7, -1.0)
        Jd = columns(Ps[None, :], beta, T, T * (1 + eps[:, None]))[..., :4]
        sd = svals(Jd)[:, 3]
        slope = np.polyfit(np.log10(eps), np.log10(sd), 1)[0]
        log(f"   beta={beta:5g} T={T:7g}: s4/s1 = {', '.join(f'{x:.2e}' for x in sd)}  (eps = 1e-1..1e-6), slope {slope:.3f}")
log("")
log("4. N = 2 frequencies: Hadamard ratio |det J_real| / prod ||col||, best of 8 random pairs P1, P2 in 10^[-2, 2]")
rng = np.random.default_rng(11)
pairs = 10.0 ** rng.uniform(-2, 2, size=(8, 2))
best = np.zeros(BB.shape)
for P1, P2 in pairs:
    Jp = columns(np.array([P1, P2])[None, None, None, :], BB[..., None], TR[..., None], TL[..., None])[..., :4]
    Jr = np.concatenate([Jp.real, Jp.imag], axis=-2)              # (..., 4, 4)
    nrm = np.prod(np.linalg.norm(Jr, axis=-2), axis=-1)
    had = np.abs(np.linalg.det(Jr)) / np.where(nrm > 0, nrm, 1)
    best = np.maximum(best, had)
for d0 in (0.25, 1.0):
    m = (~diag) & (dlog >= d0 - 1e-12)
    i = np.unravel_index(np.argmin(np.where(m, best, np.inf)), best.shape)
    log(f"   |log10(T_R/T_l)| >= {d0}: min over grid of best ratio = {best[i]:.3e} at beta={BB[i]:.3g}, T_R={TR[i]:.3g}, T_l={TL[i]:.3g}")
for name, m in [("T_l = 0 < T_R", (TL == 0) & (TR > 0)), ("T_R = 0 < T_l", (TR == 0) & (TL > 0))]:
    i = np.unravel_index(np.argmin(np.where(m, best, np.inf)), best.shape)
    log(f"   {name}: min best ratio = {best[i]:.3e} at beta={BB[i]:.3g}, T_R={TR[i]:.3g}, T_l={TL[i]:.3g}")
log(f"   on the diagonal: max best ratio = {np.max(np.where(diag, best, 0)):.2e}")
np.savez_compressed(HERE / "results" / "02_grid_N2.npz", best=best)

# ------------------------------------------------------------------ 50-digit recomputation
log("")
log("5. 50-digit recomputation (mpmath, SymPy derivatives) at the 6 lowest off-diagonal grid points (|log10 ratio| >= 0.5)")
mp.mp.dps = 50
jsym = M.sympy_jacobian_mp(with_R=False)
m = (~diag) & (dlog >= 0.5)
order = np.argsort(np.where(m, s4, np.inf), axis=None)[:6]
for flat in order:
    i = np.unravel_index(flat, s4.shape)
    rows = [jsym(mp.mpc(0, mp.mpf(float(P))), 1, float(BB[i]), float(TR[i]), float(TL[i]), 1, 0) for P in Ps]
    sv = M.singular_values(M.equilibrate(M.real_jacobian(rows)))
    log(f"   beta={BB[i]:.3g} T_R={TR[i]:.3g} T_l={TL[i]:.3g}: double s4/s1 = {s4[i]:.4e}; 50-digit s4/s1 = {mp.nstr(sv[3]/sv[0], 5)}")

(HERE / "results" / "02_sigma_numeric.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
