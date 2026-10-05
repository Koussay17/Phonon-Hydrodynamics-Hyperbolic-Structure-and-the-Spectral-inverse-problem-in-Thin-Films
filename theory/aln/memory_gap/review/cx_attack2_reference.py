"""Attack 2: is the event operator C(g_ref) consistent with phono3py's physical operator Omega'?

Omega' = D + C1 - (C0 + C2) J (phono3py channels, degeneracy-averaged, symmetrised), rebuilt here from
results/aln/physop_m553_s0.1.npz exactly as described in physop_check_local.py, restricted to the 897
non-zero modes. Compared with C_ev = C(g_ref) (orbit-averaged event rates):
  * entrywise: diagonal, off-diagonal, Frobenius; low spectrum on H
  * K, kappa, tau_mem under three treatments of Omega's (non-conserved) energy direction:
      (a) projected: P Omega' P on H = e^perp (e from the event model, D eps);
      (b) raw inverse of Omega' on the full 897-space (b orthogonal to e);
      (c) spectral: drop the eigenvector with the largest overlap with e.
  * the slowest-mode anatomy of both operators.
Output: review/cx_attack2_reference.json
"""
import json
import sys

import numpy as np
import scipy.linalg as sla

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, PHYSOP, FOURPI, basis_H  # noqa: E402

M = AlN()
gref = M.g_ref()
Cev = M.C(gref)
Z = np.load(PHYSOP)
C0, C1, C2 = Z["C0"], Z["C1"], Z["C2"]
gam = Z["gamma"].ravel()
f = Z["freqs"]
Jmap = Z["Jmap"]
nq, nb = f.shape
n = nq * nb
perm = (Jmap[:, None] * nb + np.arange(nb)[None, :]).ravel()
Om = C1 - (C0 + C2)[:, perm]
# degeneracy averaging as phono3py (blocks of degenerate bands at each q)
blocks = []
for q in range(nq):
    fq = f[q]
    j = 0
    while j < nb:
        k = j + 1
        while k < nb and abs(fq[k] - fq[j]) < 1e-4:
            k += 1
        if k - j > 1:
            blocks.append(q * nb + np.arange(j, k))
        j = k
for idx in blocks:
    Om[idx, :] = Om[idx, :].mean(axis=0)
for idx in blocks:
    Om[:, idx] = Om[:, idx].mean(axis=1)[:, None]
Om[np.arange(n), np.arange(n)] += gam
Om = 0.5 * (Om + Om.T)
v = M.valid
Om = Om[np.ix_(v, v)]
out = {"n_degenerate_blocks": len(blocks)}

# entrywise comparison (in the same entropy coordinates; phono3py's symmetrised variables)
d_ev, d_om = np.diag(Cev), np.diag(Om)
off = ~np.eye(M.n, dtype=bool)
out["diag_rel_median"] = float(np.median(np.abs(d_ev / d_om - 1)))
out["diag_rel_max"] = float(np.max(np.abs(d_ev / d_om - 1)))
out["offdiag_rel_frobenius"] = float(np.linalg.norm((Cev - Om)[off]) / np.linalg.norm(Om[off]))
out["total_rel_frobenius"] = float(np.linalg.norm(Cev - Om) / np.linalg.norm(Om))
# energy residuals
e = M.e
out["energy_residual_Omega"] = float(np.linalg.norm(Om @ e) / np.linalg.norm(d_om * e))
out["energy_residual_Cev"] = float(np.linalg.norm(Cev @ e) / np.linalg.norm(d_ev * e))

kf = M.kappa_factor()
Q = basis_H(M.ehat)


def anatomy(A, label):
    rec = {}
    AH = Q.T @ A @ Q
    AH = 0.5 * (AH + AH.T)
    w, U = np.linalg.eigh(AH)
    rec["lam_H_lowest6"] = [float(x) for x in w[:6]]
    for c in (0, 2):
        bh = U.T @ (Q.T @ M.b[:, c])
        K = float(np.sum(bh ** 2 / w)); N = float(np.sum(bh ** 2 / w ** 2))
        rec[f"projected_{c}"] = {"kappa_WmK": kf * K, "tau_ps": N / K / FOURPI}
        # raw full-space solve
        x = np.linalg.solve(A, M.b[:, c])
        K2 = float(M.b[:, c] @ x); N2 = float(x @ x)
        rec[f"raw_{c}"] = {"kappa_WmK": kf * K2, "tau_ps": N2 / K2 / FOURPI}
        # spectral: drop eigvec with largest overlap with e
        wf, Uf = np.linalg.eigh(A)
        k0 = int(np.argmax(np.abs(Uf.T @ M.ehat)))
        keep = np.arange(len(wf)) != k0
        bf = Uf.T @ M.b[:, c]
        K3 = float(np.sum(bf[keep] ** 2 / wf[keep])); N3 = float(np.sum(bf[keep] ** 2 / wf[keep] ** 2))
        rec[f"spectral_{c}"] = {"kappa_WmK": kf * K3, "tau_ps": N3 / K3 / FOURPI,
                                "dropped_eig": float(wf[k0]), "dropped_overlap_with_e": float(abs(Uf[:, k0] @ M.ehat)),
                                "b_overlap_dropped": float(bf[k0] ** 2 / np.sum(bf ** 2))}
        # share of K and N in slow modes (projected)
        relax = 1 / w / FOURPI
        Kk = bh ** 2 / w; Nk = bh ** 2 / w ** 2
        rec[f"projected_{c}"]["K_share_relax_gt_50ps"] = float(Kk[relax > 50].sum() / Kk.sum())
    return rec


out["Omega_prime"] = anatomy(Om, "Om")
out["event_operator"] = anatomy(Cev, "Cev")
# operator distance in the metric that matters: relative difference of C^+ b
for c in (0, 2):
    x_ev = M.moments(gref, (c,))[c]["x"]
    AH = Q.T @ Om @ Q
    x_om = Q @ np.linalg.solve(0.5 * (AH + AH.T), Q.T @ M.b[:, c])
    out[f"rel_diff_x_{c}"] = float(np.linalg.norm(x_ev - x_om) / np.linalg.norm(x_om))
print(json.dumps(out, indent=1))
(CAMP / "review" / "cx_attack2_reference.json").write_text(json.dumps(out, indent=1))
