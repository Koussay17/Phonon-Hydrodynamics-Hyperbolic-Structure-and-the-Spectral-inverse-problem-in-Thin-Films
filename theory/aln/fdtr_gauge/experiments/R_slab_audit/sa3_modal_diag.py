"""SA3: diagnostics of the audited modal solver (rs_core.slab_paired) and of its premises.

(1) Parity map (own code, sa_core operator): grid map 15 = (qx, qy, -qz); frequency symmetry; v_z antisymmetry;
    invariance of the assembled operator under the mode permutation; on 9x9x6 the zone-boundary plane q_z = 1/2 is
    mapped onto itself: count of modes there with v_z != 0 (the paired reduction requires v_z = 0 on self-mapped
    modes).
(2) Conditioning (audited code, REF 7x7x5): column-normalised eigenvector matrix Phi of M = v^-1 A_o v^-1 A_e,
    the boundary matrices P+ and P- of _bc_solve at 300/500 nm (2-norm condition numbers), spread of lambda,
    diffusive eigenvalue before/after the polish against the closed-form expectation mu_d ~ p C/K (1 + delta p).
(3) Diffusive-mode sensitivity: two-port with the unpolished geev eigenpair vs polished (relative change), i.e. how
    much the polish matters in the FDTR band.
Output: out/sa3_modal_diag.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sa_core as sc  # noqa: E402
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

sys.path.insert(0, str(HERE.parent / "R_slab"))
import rs_core as rc  # noqa: E402

T0 = time.perf_counter()
RES = HERE / "out" / "sa3_modal_diag.json"


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


res = {}
# ---------------------------------------------------------------- (1) parity on 7x7x5 (own operator) and 9x9x6 plane
G = sc.load_modes("m775")
C = sc.operator(G, "film2", "REF")
addr, msh, R = G["addr"], G["meshdims"], G["rot_maps"][15]
tgt = np.stack([addr[:, 0] % msh[0], addr[:, 1] % msh[1], (-addr[:, 2]) % msh[2]], 1)
ok_map = bool(np.array_equal(addr[R] % msh, tgt))
nb = G["nb"]
full = np.arange(G["valid"].size)
permf = R[full // nb] * nb + full % nb
perm = G["new"][permf[G["valid"]]]
vz, w = G["vz"], G["w"]
mov = np.abs(vz) >= 1e-12
res["parity_m775"] = {"map15_is_qz_reflection": ok_map,
                      "freq_sym_max_rel": float(np.max(np.abs(w[perm] - w) / w)),
                      "vz_antisym_max_rel": float(np.max(np.abs(vz[perm] + vz)) / np.max(np.abs(vz))),
                      "C_invariance_max_rel": float(np.max(np.abs(C[np.ix_(perm, perm)] - C)) / np.max(np.abs(C))),
                      "n_selfmapped": int(np.sum(perm == np.arange(perm.size))),
                      "n_selfmapped_moving": int(np.sum((perm == np.arange(perm.size)) & mov)),
                      "n_local": int(np.sum(~mov))}
log("parity m775", res["parity_m775"])
del C
Z = np.load(sc.EV["m996"])
fr = Z["freqs"]
valid = (fr > float(Z["cutoff"])).ravel()
vz9 = Z["gv"].reshape(-1, 3)[valid][:, 2] * 100.0
qz = np.repeat(Z["addr"][:, 2] % Z["mesh"][2], fr.shape[1])[valid]
plane = qz == Z["mesh"][2] // 2
res["m996_qz_half_plane"] = {"n_modes_plane": int(plane.sum()),
                             "n_moving_on_plane": int(np.sum(plane & (np.abs(vz9) >= 1e-12))),
                             "max_abs_vz_on_plane": float(np.max(np.abs(vz9[plane]))),
                             "median_abs_vz_on_plane": float(np.median(np.abs(vz9[plane]))),
                             "n_qz0_local": int(np.sum((qz == 0) & (np.abs(vz9) < 1e-12))),
                             "n_modes": int(valid.sum())}
# degenerate pairs on the plane: v_z of the two partners
fq = fr.reshape(-1)[valid]
qidx = np.repeat(np.arange(fr.shape[0]), fr.shape[1])[valid]
sums, mags, nonpair = [], [], 0
for q in np.unique(qidx[plane]):
    ii = np.where((qidx == q))[0]
    ff = fq[ii]
    used = np.zeros(ii.size, bool)
    for a in range(ii.size):
        if used[a]:
            continue
        partners = [b for b in range(ii.size) if b != a and not used[b] and abs(ff[b] - ff[a]) < 1e-6]
        if len(partners) == 1:
            b = partners[0]
            used[a] = used[b] = True
            sums.append(abs(vz9[ii[a]] + vz9[ii[b]]))
            mags.append(max(abs(vz9[ii[a]]), abs(vz9[ii[b]])))
        else:
            nonpair += 1
res["m996_qz_half_plane"].update({"n_degenerate_pairs": len(sums), "n_unpaired_modes": nonpair,
                                  "max_abs_vsum_pair": float(max(sums)) if sums else None,
                                  "median_abs_v_pair": float(np.median(mags)) if mags else None,
                                  "max_abs_v_pair": float(max(mags)) if mags else None})
log("m996 plane", res["m996_qz_half_plane"])
RES.write_text(json.dumps(res, indent=1, default=float))
# ---------------------------------------------------------------- (2), (3) conditioning of the audited modal solver
Gr = rc.geometry("m775")
Cr = rc.operator(Gr, "film2", "REF")
Bk = rc.blocks(Gr, Cr)
kb = rc.kinetic_times(Gr, Cr)
del Cr
rows = []
for f in (1e2, 1e4, 1.5e6, 2e8):
    p = 2j * np.pi * f
    outs, info, modes = rc.slab_paired(Bk, p, (300e-9, 500e-9), return_modes=True)
    lam, Phi, Psi = modes["lam"], modes["Phi"], modes["Psi"]
    Pc = Phi / np.linalg.norm(Phi, axis=0)
    sv = sla.svdvals(Pc)
    row = {"f": f, "cond_Phi_colnorm": float(sv[0] / sv[-1]), "lam_min": float(np.abs(lam).min()),
           "lam_max": float(np.abs(lam).max()), "mu_d_geev": str(info["mu_d_geev"]),
           "mu_d_polished": str(info["mu_d_polished"]), "polish_shift_rel": info["polish_shift_rel"],
           "eig_backward_max": info["eig_backward_max"]}
    # closed-form first-order expectation of the diffusive eigenvalue: p C/K (1 + delta p)
    mu_exp = p * Gr["C_tot"] / kb["K"] * (1 + kb["delta"] * p)
    row["mu_d_rel_dev_from_first_order"] = float(abs(info["mu_d_polished"] / mu_exp - 1))
    for d in (300e-9, 500e-9):
        E = np.exp(-lam * d)
        Pp = (Phi + Psi) + (Phi - Psi) * E[None, :]
        Pm = (Phi + Psi) - (Phi - Psi) * E[None, :]
        for tag, M_ in (("Pplus", Pp), ("Pminus", Pm)):
            Mc = M_ / np.linalg.norm(M_, axis=0)
            s_ = sla.svdvals(Mc)
            row[f"cond_{tag}_{d * 1e9:.0f}"] = float(s_[0] / s_[-1])
    # (3) unpolished vs polished two-port
    outs_u, info_u = rc.slab_paired(Bk, p, (300e-9, 500e-9), polish=False)
    for (Yp, Thp, _), (Yu, Thu, _), d in zip(outs, outs_u, (300, 500)):
        row[f"polish_effect_Y_{d}"] = float(np.max(np.abs(Yp - Yu)) / np.max(np.abs(Yp)))
        row[f"polish_effect_ImY11_rel_{d}"] = float(abs((Yp[0, 0] - Yu[0, 0]).imag) / abs(Yp[0, 0].imag))
    rows.append(row)
    log(row)
    res["conditioning_REF_m775"] = rows
    RES.write_text(json.dumps(res, indent=1, default=float))
