"""Out-of-core LBTE for AlN: phono3py collision rows -> disk -> own assembly -> PCG.

Purpose
-------
(1) Second, independent solution route for the LBTE thermal conductivity: per-grid-
    point collision rows and linewidths come from phono3py's CollisionMatrix kernel,
    but the global assembly (stabilizer diagonal, star weights, degeneracy averaging,
    symmetrization) is re-implemented here and the linear system is solved by
    projected Jacobi-preconditioned conjugate gradients instead of phono3py's dense
    eigendecomposition pseudo-inverse.  Validated bit-for-bit (assembled matrix) and
    to ~1e-12 (kappa) against phono3py's dense LBTE at 9x9x5.
(2) Memory: the assembled matrix lives in a .npy memmap on disk, so meshes whose
    dense matrix does not fit in RAM (31x31x17: n = 31104, 7.2 GiB) can be solved.

Conventions reproduced (phono3py 4.5.0 irreducible "kappa-star" representation)
-----------------------------------------------------------------------------
Index (i, b, a): irreducible grid point i, band b, Cartesian a; n = n_ir*nb*3.
  A_raw[(i b a),(j c d)] = sum_{R: q' = R q_j} Omega_{(q_i b),(q' c)} R_ad   (3-phonon, THz)
                         + delta_ij delta_bc sum_{R in stab(q_i)} Gamma_ib R_ad
  A_w = w_i w_j A_raw,   w_i = sqrt(|star(q_i)| / |G|),  |G| = 24 (6mm x time reversal)
  A   = sym(P A_w P),    P = average over degenerate bands, sym(M) = (M + M^T)/2
  X_(i b a) = w_i v_iba f_ib THzToEv / (4 k_B T^2 sinh(f_ib THzToEv / 2 k_B T))
  Y   = A^+ X  -> projected PCG on range(Pi), Pi = P (bands) x S_i (stabilizer average)
  kappa = conv k_B T^2 / N sum_{i,b} sum_R [(R x_ib)(R y_ib)^T + transpose],
  conv = phono3py get_unit_to_WmK()/V.  Same with Y = X/Gamma gives kappa_RTA.
Optional isotope scattering: Gamma_iso (from a phono3py RTA --isotope kappa file) is
added to the diagonal blocks as w_i^2 Gamma_iso sum_{R in stab} R (inside Pi).

Stages (resumable; outputs in --out)
------------------------------------
  rows     : rows for all irreducible points at all --temps; writes A_T<k>.npy (A_w),
             gamma.npy, meta.npz, progress.json.
  assemble : optional copy of A_T<k> before averaging (--keep-raw k), then in-place
             degeneracy averaging + symmetrization.
  solve    : PCG per temperature (optionally with --gamma-iso-from); solve_T<k>[_iso].json
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def build_phono3py(args):
    from phono3py.cui.load import load
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"),
               supercell_matrix=common.SUPERCELL_FC3,
               phonon_supercell_matrix=common.SUPERCELL_FC2,
               fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"),
               born_filename=str(common.DATASET / "BORN"), is_nac=True,
               symmetrize_fc=False, make_r0_average=not args.no_r0_average,
               log_level=0, lang=args.lang)
    ph3.mesh_numbers = args.mesh
    ph3.sigmas = [args.sigma]
    ph3.init_phph_interaction()
    ph3.run_phonon_solver()
    return ph3


def stage_rows(args):
    from phono3py.conductivity.velocity_solvers import GroupVelocitySolver
    from phono3py.phonon3.collision_matrix import CollisionMatrix
    from phonopy.phonon.grid import get_grid_points_by_rotations, get_ir_grid_points

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rec_path = out / "run.json"
    if not rec_path.exists():
        common.write_json(rec_path, {
            "argv": sys.argv, "settings": vars(args), "inputs": common.verify_inputs(),
            "software": common.software_versions(), "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "script_sha256": common.sha256(__file__), "dataset_dir": str(common.DATASET)})
    else:
        common.verify_inputs()
    ph3 = build_phono3py(args)
    itr = ph3.phph_interaction
    bz = itr.bz_grid
    freqs_all, _, _ = itr.get_phonons()
    ir_grg, ir_weights, _ = get_ir_grid_points(bz)
    gps = np.array(bz.grg2bzg[ir_grg], dtype="int64")
    rots = bz.rotations
    rcart = np.array(bz.rotations_cartesian)
    rot_gps = np.array([get_grid_points_by_rotations(gp, bz, reciprocal_rotations=rots) for gp in gps],
                       dtype="int64")
    n_ir, n_rot = rot_gps.shape
    nb = freqs_all.shape[1]
    n = n_ir * nb * 3
    w = np.array([math.sqrt(len(np.unique(r))) for r in rot_gps]) / math.sqrt(n_rot)
    temps = np.array(args.temps, dtype=float)

    prog_path = out / "progress.json"
    if prog_path.exists():
        prog = json.loads(prog_path.read_text())
        assert prog["mesh"] == list(args.mesh) and prog["temps"] == list(map(float, temps)), "settings changed"
        done = set(prog["done"])
        gamma = np.load(out / "gamma_partial.npy")
    else:
        prog = {"mesh": list(args.mesh), "temps": list(map(float, temps)), "done": [], "n": n, "n_ir": n_ir}
        done = set()
        gamma = np.zeros((len(temps), n_ir, nb))
    mms = []
    for k in range(len(temps)):
        p = out / f"A_T{k}.npy"
        mms.append(np.lib.format.open_memmap(p, mode=("r+" if p.exists() else "w+"),
                                             dtype=np.float64, shape=(n, n)))
    gv_solver = GroupVelocitySolver(itr, is_kappa_star=True, gv_delta_q=args.gv_delta_q)
    gv = np.array([gv_solver.compute(int(gp)).group_velocities for gp in gps])
    np.savez(out / "meta.npz", gps=gps, ir_weights=ir_weights, w=w, rot_gps=rot_gps, rcart=rcart,
             freqs=freqs_all[gps], gv=gv, temps=temps, mesh=np.array(args.mesh),
             qpoints=bz.addresses[gps] @ np.linalg.inv(np.diag(bz.D_diag)).T,
             volume=ph3.primitive.volume, cutoff_frequency=itr.cutoff_frequency,
             sigma=np.nan if args.sigma is None else args.sigma)
    col = CollisionMatrix(itr, rot_grid_points=rot_gps, is_kappa_star=True, lang=args.lang)
    t0 = time.perf_counter()
    for i, gp in enumerate(gps):
        if i in done:
            continue
        col.set_grid_point(int(gp))
        col.set_sigma(args.sigma, sigma_cutoff=None)
        col.run_integration_weights()
        col.run_interaction(is_full_pp=False)
        for k, T in enumerate(temps):
            col.temperature = float(T)
            col.run()
            gam = np.array(col.imag_self_energy, dtype=float)
            row = np.array(col.get_collision_matrix(), dtype=float)  # (nb, 3, n_ir, nb, 3)
            for r, rgp in zip(rcart, rot_gps[i]):
                if rgp == gp:
                    for ll in range(nb):
                        row[ll, :, i, ll, :] += gam[ll] * r
            row *= w[i] * w[None, None, :, None, None]
            mms[k][i * nb * 3:(i + 1) * nb * 3, :] = row.reshape(nb * 3, n)
            gamma[k, i] = gam
        done.add(i)
        if (len(done) % 50 == 0) or len(done) == n_ir:
            for mm in mms:
                mm.flush()
            np.save(out / "gamma_partial.npy", gamma)
            prog["done"] = sorted(done)
            prog_path.write_text(json.dumps(prog))
            print(f"[rows] {len(done)}/{n_ir} ir points, {time.perf_counter() - t0:.0f}s", flush=True)
    for mm in mms:
        mm.flush()
    np.save(out / "gamma.npy", gamma)
    prog["done"] = sorted(done)
    prog["rows_complete"] = True
    prog["rows_seconds"] = prog.get("rows_seconds", 0.0) + time.perf_counter() - t0
    prog_path.write_text(json.dumps(prog))


def degeneracy_projectors(freqs):
    from phonopy.phonon.degeneracy import degenerate_sets
    n_ir, nb = freqs.shape
    P = np.zeros((n_ir, nb, nb))
    for i in range(n_ir):
        for dset in degenerate_sets(freqs[i]):
            idx = np.array(dset)
            P[i][np.ix_(idx, idx)] = 1.0 / len(idx)
    return P


def stage_assemble(args):
    out = Path(args.out)
    meta = np.load(out / "meta.npz")
    prog = json.loads((out / "progress.json").read_text())
    assert prog.get("rows_complete"), "rows stage incomplete"
    freqs = meta["freqs"]
    n_ir, nb = freqs.shape
    n = n_ir * nb * 3
    P = degeneracy_projectors(freqs)
    bs = args.block
    for k in range(len(meta["temps"])):
        flag = out / f"assembled_T{k}.flag"
        if flag.exists():
            continue
        if args.keep_raw is not None and k in args.keep_raw:
            raw = out / f"A_T{k}_prephono3pyavg_presym.npy"
            if not raw.exists():
                shutil.copyfile(out / f"A_T{k}.npy", raw)
        mm = np.lib.format.open_memmap(out / f"A_T{k}.npy", mode="r+")
        t0 = time.perf_counter()
        asym_num = 0.0
        asym_den = 0.0
        for i0 in range(0, n_ir, bs):
            i1 = min(n_ir, i0 + bs)
            blk = np.array(mm[i0 * nb * 3:i1 * nb * 3, :]).reshape(i1 - i0, nb, 3, n_ir, nb, 3)
            blk = np.einsum("iab,ibxjcy->iaxjcy", P[i0:i1], blk, optimize=True)
            blk = np.einsum("iaxjcy,jdc->iaxjdy", blk, P, optimize=True)
            mm[i0 * nb * 3:i1 * nb * 3, :] = blk.reshape((i1 - i0) * nb * 3, n)
        mm.flush()
        t1 = time.perf_counter()
        edges = list(range(0, n_ir, bs)) + [n_ir]
        maxabs = 0.0
        for a in range(len(edges) - 1):
            ra = slice(edges[a] * nb * 3, edges[a + 1] * nb * 3)
            for b in range(a, len(edges) - 1):
                rb = slice(edges[b] * nb * 3, edges[b + 1] * nb * 3)
                Aab = np.array(mm[ra, rb])
                Aba = Aab if a == b else np.array(mm[rb, ra])
                asym_num = max(asym_num, float(np.abs(Aab - Aba.T).max()))
                maxabs = max(maxabs, float(np.abs(Aab).max()), float(np.abs(Aba).max()))
                S = 0.5 * (Aab + Aba.T)
                mm[ra, rb] = S
                if a != b:
                    mm[rb, ra] = S.T
        mm.flush()
        del mm
        info = {"deg_avg_s": t1 - t0, "sym_s": time.perf_counter() - t1,
                "pre_symmetrization_max_asym": asym_num, "max_abs": maxabs,
                "pre_symmetrization_max_asym_rel": asym_num / maxabs}
        flag.write_text(json.dumps(info))
        print(f"[assemble] T index {k}: {info}", flush=True)


def x_vector(meta, T):
    from phonopy.physical_units import get_physical_units
    u = get_physical_units()
    f = meta["freqs"]
    cutoff = float(meta["cutoff_frequency"])
    with np.errstate(over="ignore"):
        s = np.where(f > cutoff, np.sinh(f * u.THzToEv / (2 * u.KB * T)), -1.0)
    inv_s = np.where(s > 0, 1.0 / np.where(s > 0, s, 1.0), 0.0)
    ff = f * u.THzToEv * inv_s / (4 * u.KB * T**2)
    return meta["gv"] * meta["w"][:, None, None] * ff[:, :, None]


def kappa_from_xy(meta, X, Y, T):
    from phono3py.conductivity.utils import get_unit_to_WmK
    from phonopy.physical_units import get_physical_units
    u = get_physical_units()
    conv = get_unit_to_WmK() / float(meta["volume"])
    R = meta["rcart"]
    RX = np.einsum("rab,ijb->rija", R, X)
    RY = np.einsum("rab,ijb->rija", R, Y)
    K = np.einsum("rija,rijb->ab", RX, RY)
    K = K + K.T
    return K * conv * u.KB * T**2 / int(np.prod(meta["mesh"]))


def stabilizer_projectors(meta):
    gps, rot_gps, R = meta["gps"], meta["rot_gps"], meta["rcart"]
    S = np.zeros((len(gps), 3, 3))
    stab_sum = np.zeros((len(gps), 3, 3))
    for i, gp in enumerate(gps):
        m = rot_gps[i] == gp
        S[i] = R[m].mean(axis=0)
        stab_sum[i] = R[m].sum(axis=0)
    return S, stab_sum


def stage_solve(args):
    out = Path(args.out)
    meta = dict(np.load(out / "meta.npz"))
    gamma = np.load(out / "gamma.npy")
    freqs = meta["freqs"]
    n_ir, nb = freqs.shape
    n = n_ir * nb * 3
    P = degeneracy_projectors(freqs)
    S, stab_sum = stabilizer_projectors(meta)
    w = meta["w"]
    g_iso = None
    tag = ""
    if args.gamma_iso_from:
        import h5py
        with h5py.File(args.gamma_iso_from, "r") as f:
            g_iso = np.array(f["gamma_isotope"][()]).reshape(n_ir, nb)
            fr = f["frequency"][()]
        assert np.abs(fr - freqs).max() < 1e-8, "isotope file frequencies differ"
        tag = "_iso"

    def proj(v):
        V = v.reshape(n_ir, nb, 3)
        V = np.einsum("iab,ibx->iax", P, V)
        V = np.einsum("ixy,iay->iax", S, V)
        return V.ravel()

    temps_sel = range(len(meta["temps"])) if args.temp_index is None else args.temp_index
    for k in temps_sel:
        T = float(meta["temps"][k])
        assert (out / f"assembled_T{k}.flag").exists(), "assemble first"
        mm = np.lib.format.open_memmap(out / f"A_T{k}.npy", mode="r")
        t0 = time.perf_counter()
        diag = np.empty(n)
        for r0 in range(0, n, args.rows_per_chunk):
            r1 = min(n, r0 + args.rows_per_chunk)
            blk = np.asarray(mm[r0:r1, r0:r1])
            diag[r0:r1] = np.diagonal(blk)
        iso_blocks = None
        if g_iso is not None:
            iso_blocks = (w**2)[:, None, None, None] * g_iso[:, :, None, None] * stab_sum[:, None, :, :]
            diag = diag + np.einsum("ibxx->ibx", iso_blocks).ravel()

        def mv(x):
            y = np.empty(n)
            for r0 in range(0, n, args.rows_per_chunk):
                r1 = min(n, r0 + args.rows_per_chunk)
                y[r0:r1] = np.asarray(mm[r0:r1, :]) @ x
            if iso_blocks is not None:
                y += np.einsum("ibxy,iby->ibx", iso_blocks, x.reshape(n_ir, nb, 3)).ravel()
            return y

        active = diag > args.active_rel * np.abs(diag).max()
        M_inv = np.where(active, 1.0 / np.where(active, diag, 1.0), 0.0)
        X = x_vector(meta, T)
        b_raw = X.ravel()
        b = proj(b_raw)
        x = np.zeros(n)
        r = b.copy()
        z = proj(M_inv * r)
        p = z.copy()
        rz = r @ z
        bnorm = np.linalg.norm(b)
        hist = []
        alphas, betas = [], []
        min_curv = np.inf
        for it in range(1, args.maxiter + 1):
            Ap = proj(mv(p))
            pAp = p @ Ap
            min_curv = min(min_curv, pAp / (p @ p))
            alpha = rz / pAp
            alphas.append(alpha)
            x += alpha * p
            r -= alpha * Ap
            hist.append(float(np.linalg.norm(r) / bnorm))
            if it % 25 == 0:
                print(f"   T={T:g} pcg it {it}: rel residual {hist[-1]:.3e}", flush=True)
            if hist[-1] < args.tol:
                break
            z = proj(M_inv * r)
            rz_new = r @ z
            betas.append(rz_new / rz)
            p = z + (rz_new / rz) * p
            rz = rz_new
        # Lanczos tridiagonal from PCG coefficients -> Ritz values of the preconditioned operator
        kL = len(alphas)
        Tm = np.zeros((kL, kL))
        for j in range(kL):
            Tm[j, j] = 1.0 / alphas[j] + (betas[j - 1] / alphas[j - 1] if j > 0 else 0.0)
            if j + 1 < kL:
                Tm[j, j + 1] = Tm[j + 1, j] = np.sqrt(betas[j]) / alphas[j]
        ritz = np.linalg.eigvalsh(Tm)
        r_true = b - mv(x)
        Y = x.reshape(n_ir, nb, 3)
        kap = kappa_from_xy(meta, X, Y, T)
        g = gamma[k] + (g_iso if g_iso is not None else 0.0)
        valid = (freqs > float(meta["cutoff_frequency"]))[:, :, None] & (g[:, :, None] > 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            Yr = np.where(valid, X / g[:, :, None], 0.0)
        kap_rta = kappa_from_xy(meta, X, Yr, T)
        res = {"T": T, "isotope": g_iso is not None, "n": n, "active": int(active.sum()),
               "pcg_iterations": len(hist), "pcg_final_rel_residual": hist[-1],
               "true_rel_residual": float(np.linalg.norm(r_true) / bnorm),
               "true_rel_residual_projected": float(np.linalg.norm(proj(r_true)) / bnorm),
               "rhs_fraction_removed_by_projection": float(np.linalg.norm(b_raw - b) / np.linalg.norm(b_raw)),
               "kappa_xx": float(kap[0, 0]), "kappa_yy": float(kap[1, 1]), "kappa_zz": float(kap[2, 2]),
               "kappa_offdiag_max": float(np.abs(kap - np.diag(np.diag(kap))).max()),
               "kappa_RTA_xx": float(kap_rta[0, 0]), "kappa_RTA_zz": float(kap_rta[2, 2]),
               "min_curvature_pAp_over_pp_THz": float(min_curv),
               "ritz_preconditioned_min": float(ritz.min()), "ritz_preconditioned_max": float(ritz.max()),
               "seconds": time.perf_counter() - t0, "residual_history": hist}
        np.save(out / f"Y_T{k}{tag}.npy", Y)
        (out / f"solve_T{k}{tag}.json").write_text(json.dumps(res, indent=2))
        print(json.dumps({kk: v for kk, v in res.items() if kk != "residual_history"}), flush=True)
        del mm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["rows", "assemble", "solve", "all"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--mesh", type=int, nargs=3)
    ap.add_argument("--temps", type=float, nargs="+", default=[300.0])
    ap.add_argument("--sigma", type=float, default=None, help="Gaussian width THz (default tetrahedron)")
    ap.add_argument("--lang", choices=["C", "Rust"], default="Rust")
    ap.add_argument("--no-r0-average", action="store_true")
    ap.add_argument("--gv-delta-q", type=float, default=None)
    ap.add_argument("--block", type=int, default=32)
    ap.add_argument("--keep-raw", type=int, nargs="*", default=None,
                    help="temperature indices whose pre-averaging/pre-symmetrization matrix is kept")
    ap.add_argument("--rows-per-chunk", type=int, default=2048)
    ap.add_argument("--active-rel", type=float, default=1e-12)
    ap.add_argument("--tol", type=float, default=1e-10)
    ap.add_argument("--maxiter", type=int, default=3000)
    ap.add_argument("--temp-index", type=int, nargs="*", default=None)
    ap.add_argument("--gamma-iso-from", default=None, help="phono3py kappa hdf5 with gamma_isotope")
    args = ap.parse_args()
    t0 = time.perf_counter()
    if args.stage in ("rows", "all"):
        stage_rows(args)
    if args.stage in ("assemble", "all"):
        stage_assemble(args)
    if args.stage in ("solve", "all"):
        stage_solve(args)
    print(f"[done] {args.stage} in {time.perf_counter() - t0:.0f}s; memory {common.memory_info()}", flush=True)


if __name__ == "__main__":
    main()
