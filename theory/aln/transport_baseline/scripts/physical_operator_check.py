"""Even-sector (energy) test of phono3py's reducible collision matrix for AlN.

Hypothesis tested (thesis repository, operator-parity check of phono3py 4.5.0):
phono3py's reducible matrix is Chaput's conductivity-equivalent operator
    Omega  = D + C0 + C1 + C2,
while the physical linearized (symmetrized-variable) operator is
    Omega' = D + C1 - (C0 + C2) J,
with D = diag(Gamma) (phono3py linewidth), Ck the off-diagonal contribution of the
triplet channel k (g0: d(w0-w1-w2), g1: d(w0+w1-w2), g2: d(w0-w1+w2), triplet
q0+q1+q2=G, column index = q1) and J: (q, j) -> (-q, j).
Consequences: Omega - Omega' = (C0+C2)(I+J); the two agree on odd vectors; the
energy vector e_(q j) = f_qj / sinh(h f_qj / 2 k_B T) (even) is annihilated by
Omega' up to discretization error but not by Omega.

Method: for every grid point of the full mesh, phono3py's CollisionMatrix
(reducible mode) is run; its interaction strengths pp and the three channel weights
are read (phono3py stores g[0]=g0, g[1]=g1-g2, g[2]=g0+g1+g2, so g1 and g2 are
recovered exactly), channel rows are rebuilt and checked against phono3py's own row
(C0+C1+C2), then dense matrices are formed and the vectors tested.  Also checks
symmetry of the raw matrices and parity consistency on random odd/even vectors.

usage: python physical_operator_check.py --mesh 9 9 5 [--sigma 0.1] [--lang Rust]
       [--no-r0-average] [--T 300] [--json out.json] [--save-npz out.npz]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=int, nargs=3, required=True)
    ap.add_argument("--sigma", type=float, default=None, help="Gaussian width (THz); default tetrahedron")
    ap.add_argument("--lang", default="Rust")
    ap.add_argument("--no-r0-average", action="store_true")
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--json", default=None)
    ap.add_argument("--save-npz", default=None)
    a = ap.parse_args()
    t0 = time.perf_counter()

    from phono3py.cui.load import load
    from phono3py.phonon3.collision_matrix import CollisionMatrix
    from phonopy.phonon.grid import get_grid_point_from_address
    from phonopy.physical_units import get_physical_units
    u = get_physical_units()

    common.verify_inputs()
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=[3, 3, 2],
               phonon_supercell_matrix=[5, 5, 3], fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"),
               is_nac=True, symmetrize_fc=False, make_r0_average=not a.no_r0_average,
               log_level=0, lang=a.lang)
    ph3.mesh_numbers = a.mesh
    ph3.sigmas = [a.sigma]
    ph3.init_phph_interaction()
    ph3.run_phonon_solver()
    itr = ph3.phph_interaction
    bz = itr.bz_grid
    D = bz.D_diag
    N = int(np.prod(D))
    freqs_bz, _, _ = itr.get_phonons()
    nb = freqs_bz.shape[1]
    n = N * nb
    T = a.T
    cutoff = itr.cutoff_frequency

    # GR-grid helpers
    grg2bzg = bz.grg2bzg
    addr = bz.addresses[grg2bzg]               # GR addresses (BZ representative)
    Jmap = np.array([get_grid_point_from_address(-ad, D) for ad in addr], dtype=int)
    assert np.array_equal(Jmap[Jmap], np.arange(N))
    f = freqs_bz[grg2bzg]                      # (N, nb) frequencies on GR grid
    fJ_err = float(np.abs(f[Jmap] - f).max())   # time-reversal check w(-q)=w(q)

    col = CollisionMatrix(itr, rot_grid_points=None, lang=a.lang)
    C = np.zeros((3, n, n))
    gam = np.zeros((N, nb))       # degeneracy-averaged linewidth (what phono3py puts on the diagonal)
    gam_raw = np.zeros((N, nb))   # per-band linewidth before degeneracy averaging
    row_check = 0.0
    neg_weight_min = 0.0
    for g0 in range(N):
        gp0 = int(grg2bzg[g0])
        col.set_grid_point(gp0)
        col.set_sigma(a.sigma, sigma_cutoff=None)
        col.run_integration_weights()
        col.run_interaction(is_full_pp=False)
        col.temperature = T
        col.run()
        gam[g0] = col.imag_self_energy
        gam_raw[g0] = np.array(col._imag_self_energy)
        ref_row = np.array(col.get_collision_matrix())          # (nb, N, nb)
        pp = col._pp_strength                                    # (ntp, nb0, nb1, nb2)
        g = col._g                                               # (3, ntp, nb0, nb1, nb2)
        ch0 = g[0]
        ch1 = 0.5 * (g[2] - g[0] + g[1])
        ch2 = 0.5 * (g[2] - g[0] - g[1])
        neg_weight_min = min(neg_weight_min, float(ch1.min()), float(ch2.min()), float(ch0.min()))
        gp2tp, tp2s, swapped = col._get_gp2tp_map()
        conv = col._unit_conversion
        inv_sinh_all = np.zeros_like(freqs_bz)
        with np.errstate(over="ignore", divide="ignore"):
            s = np.where(freqs_bz > cutoff, np.sinh(freqs_bz * u.THzToEv / (2 * u.KB * T)), -1.0)
            inv_sinh_all = np.where(s > 0, 1.0 / np.where(s > 0, s, 1.0), 0.0)
        for g1 in range(N):
            ti = gp2tp[g1]
            isn = inv_sinh_all[tp2s[g1]]                         # third phonon bands
            if swapped[g1]:
                # stored triplet (q0, q2, q1): our q1 band k sits in the last slot
                c0 = np.einsum("jlk,jlk,l->jk", pp[ti], ch0[ti], isn)
                c1 = np.einsum("jlk,jlk,l->jk", pp[ti], ch2[ti], isn)   # roles of g1,g2 exchange
                c2 = np.einsum("jlk,jlk,l->jk", pp[ti], ch1[ti], isn)
            else:
                c0 = np.einsum("jkl,jkl,l->jk", pp[ti], ch0[ti], isn)
                c1 = np.einsum("jkl,jkl,l->jk", pp[ti], ch1[ti], isn)
                c2 = np.einsum("jkl,jkl,l->jk", pp[ti], ch2[ti], isn)
            rows = slice(g0 * nb, (g0 + 1) * nb)
            cols = slice(g1 * nb, (g1 + 1) * nb)
            C[0, rows, cols] = conv * c0
            C[1, rows, cols] = conv * c1
            C[2, rows, cols] = conv * c2
        mine = (C[0, g0 * nb:(g0 + 1) * nb] + C[1, g0 * nb:(g0 + 1) * nb]
                + C[2, g0 * nb:(g0 + 1) * nb]).reshape(nb, N, nb)
        row_check = max(row_check, float(np.abs(mine - ref_row).max() / max(np.abs(ref_row).max(), 1e-300)))

    # ---------------------------------------------------------------- operators
    # memory-lean: matrices are built one at a time; P (degeneracy averaging) applied blockwise
    from phonopy.phonon.degeneracy import degenerate_sets
    perm = (Jmap[:, None] * nb + np.arange(nb)[None, :]).ravel()   # (J x)[(q,j)] = x[(-q,j)]
    blocks = []
    for q in range(N):
        for dset in degenerate_sets(f[q]):
            if len(dset) > 1:
                blocks.append(q * nb + np.array(dset))

    def P_rows_cols(M):
        for idx in blocks:
            M[idx, :] = M[idx, :].mean(axis=0)
        for idx in blocks:
            M[:, idx] = M[:, idx].mean(axis=1)[:, None]
        return M

    def Pvec(v):
        v = v.copy()
        for idx in blocks:
            v[idx] = v[idx].mean()
        return v

    valid = (f > cutoff).ravel()
    with np.errstate(over="ignore"):
        sh = np.where(f > cutoff, np.sinh(f * u.THzToEv / (2 * u.KB * T)), 1.0)
    e = np.where(f > cutoff, f / sh, 0.0).ravel()                   # energy vector (symmetrized vars)
    De = gam_raw.ravel() * e

    def rel(v):
        return float(np.linalg.norm(v[valid]) / np.linalg.norm(De[valid]))

    def asym(M):
        return float(np.abs(M - M.T).max() / np.abs(M).max())

    def jcomm(M):
        return float(np.abs(M[np.ix_(perm, perm)] - M).max() / np.abs(M).max())

    rng = np.random.default_rng(20261003)
    x = rng.standard_normal(n)
    x_odd = Pvec(0.5 * (x - x[perm]))
    x_even = Pvec(0.5 * (x + x[perm]))
    res = {
        "mesh": a.mesh, "sigma": a.sigma, "lang": a.lang, "make_r0_average": not a.no_r0_average, "T": T,
        "n": n, "time_reversal_freq_maxdiff_THz": fJ_err,
        "channel_rebuild_vs_phono3py_row_maxrel": row_check,
        "min_channel_weight": neg_weight_min,
        "gamma_raw_vs_avg_maxrel": float(np.abs(gam_raw - gam).max() / gam.max()),
        "timereversal_gamma_raw_maxrel": float(np.abs(gam_raw[Jmap] - gam_raw).max() / gam_raw.max()),
        "timereversal_gamma_avg_maxrel": float(np.abs(gam[Jmap] - gam).max() / gam.max()),
        "asym_C0": asym(C[0]), "asym_C1": asym(C[1]), "asym_C2": asym(C[2]),
    }
    idx_d = np.arange(n)
    # (a) phono3py-type operator, raw rows + raw linewidth
    M = C[0] + C[1] + C[2]
    M[idx_d, idx_d] += gam_raw.ravel()
    Me_O = M @ e
    res["energy_residual_Omega_raw"] = rel(Me_O)
    res["energy_rayleigh_Omega_raw"] = float(e @ Me_O / (e @ De))
    res["asym_Omega_raw"] = asym(M)
    Mx_odd_O, Mx_even_O = M @ x_odd, M @ x_even
    # (b) physical operator, raw
    M = C[1] - (C[0] + C[2])[:, perm]
    M[idx_d, idx_d] += gam_raw.ravel()
    Me_P = M @ e
    res["energy_residual_OmegaPrime_raw"] = rel(Me_P)
    res["energy_rayleigh_OmegaPrime_raw"] = float(e @ Me_P / (e @ De))
    res["asym_OmegaPrime_raw"] = asym(M)
    res["identity_(Omega-OmegaPrime)e_vs_2(C0+C2)e"] = float(
        np.linalg.norm((Me_O - Me_P) - 2 * ((C[0] + C[2]) @ e)) / np.linalg.norm(De))
    res["odd_sector_|(Omega-OmegaPrime)x|/|Omega x|"] = float(
        np.linalg.norm(Mx_odd_O - M @ x_odd) / np.linalg.norm(Mx_odd_O))
    res["even_sector_|(Omega-OmegaPrime)x|/|Omega x|"] = float(
        np.linalg.norm(Mx_even_O - M @ x_even) / np.linalg.norm(Mx_even_O))
    # (c) physical operator with phono3py-style degeneracy averaging and averaged linewidth
    M[idx_d, idx_d] -= gam_raw.ravel()
    P_rows_cols(M)
    M[idx_d, idx_d] += gam.ravel()
    res["energy_residual_OmegaPrime_degavg"] = rel(M @ e)
    res["asym_OmegaPrime_degavg"] = asym(M)
    res["JMJ_rel_OmegaPrime_degavg"] = jcomm(M)
    y = M @ x_odd
    res["parity_leak_OmegaPrime_degavg_odd_to_even"] = float(np.linalg.norm(0.5 * (y + y[perm])) / np.linalg.norm(y))
    M = 0.5 * (M + M.T)                                            # phono3py-style (A + A^T)/2
    res["energy_residual_symOmegaPrime_degavg"] = rel(M @ e)
    res["energy_rayleigh_symOmegaPrime_degavg"] = float(e @ (M @ e) / (e @ De))
    w = np.linalg.eigvalsh(M)
    res["eig_symOmegaPrime_degavg_smallest6"] = [float(v) for v in np.sort(w)[:6]]
    res["eig_symOmegaPrime_degavg_num_below_-1e-10"] = int((w < -1e-10).sum())
    res["eig_symOmegaPrime_degavg_max"] = float(w.max())
    # (d) phono3py-type operator degeneracy averaged + symmetrized (what phono3py assembles)
    M = C[0] + C[1] + C[2]
    P_rows_cols(M)
    M[idx_d, idx_d] += gam.ravel()
    res["asym_Omega_degavg"] = asym(M)
    res["JMJ_rel_Omega_degavg"] = jcomm(M)
    M = 0.5 * (M + M.T)
    res["energy_residual_symOmega_degavg"] = rel(M @ e)
    res["energy_rayleigh_symOmega_degavg"] = float(e @ (M @ e) / (e @ De))
    w = np.linalg.eigvalsh(M)
    res["eig_symOmega_degavg_smallest6"] = [float(v) for v in np.sort(w)[:6]]
    res["eig_symOmega_degavg_max"] = float(w.max())
    del M
    res["gamma_max"] = float(gam.max())
    res["seconds"] = time.perf_counter() - t0
    print(json.dumps(res, indent=1))
    if a.json:
        common.write_json(a.json, res)
    if a.save_npz:
        np.savez_compressed(a.save_npz, C0=C[0], C1=C[1], C2=C[2], gamma=gam, freqs=f, Jmap=Jmap, e=e)


if __name__ == "__main__":
    main()
