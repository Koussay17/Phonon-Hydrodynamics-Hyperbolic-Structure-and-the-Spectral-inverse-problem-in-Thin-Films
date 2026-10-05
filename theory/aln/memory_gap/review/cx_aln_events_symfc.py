# REVIEW COPY of scripts/aln_events.py with symmetrize_fc=True (selection-rule test); output to review/cx_runs only.
"""Export the AlN three-phonon event geometry (physical-operator convention) from phono3py 4.5.0.

Run with the phono3py environment (D:\\ResearchLab\\envs\\aln-phono3py-4.5.0).

For every full-grid parent mode p = (q0, j0) and every pair (q1, j1), (q2, l) of phono3py's
triplet q0 + q1 + q2 = G, the decay channel g0 = delta_sigma(w0 - w1 - w2) defines the event
    p -> a + b,  a = (-q1, j1),  b = (-q2, l).
In phono3py's physical operator Omega' = D + C1 - (C0 + C2) J (thesis campaign
20261002-aln-transport-baseline) this event contributes -conv*pp*g0/sinh(x_b/2) to
Omega'[p, a]; matching g * z_p z_a with z = D^{-1} s, D = diag(1/(2 sinh(x/2))) gives the rate
    g = conv * pp * g0 / (4 sinh(x_p/2) sinh(x_a/2) sinh(x_b/2))      (a != b),
    g = conv * pp * g0 / (8 sinh(x_p/2) sinh(x_a/2)^2), s = -e_p + 2 e_a   (a == b),
which is symmetric in the three participants. Each unordered event is emitted once
(canonical a <= b). Events are kept if |w0 - w1 - w2| <= nsig*sigma and pp*g0 > 0.
Conserving stoichiometry (declared): s = -l_p e_p + l_d (e_a + e_b), l_p = sqrt((w_a+w_b)/w_p),
l_d = 1/l_p, so s^T w = 0 exactly; the rank-one event operator is PSD and conserves energy.
Units: g in THz (phono3py collision-matrix units); physical rates are 4*pi*g in 1/ps.

usage: python aln_events.py --mesh 5 5 3 [--sigma 0.1] [--nsig 4] [--T 300] --out file.npz
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

BASE = Path(r"D:/ResearchLab/orchestration/campaigns/20261004-memory-gap/scripts")
sys.path.insert(0, str(BASE))
sys.dont_write_bytecode = True
import aln_common as common  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=int, nargs=3, required=True)
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--nsig", type=float, default=4.0)
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--lang", default="Rust")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.perf_counter()

    from phono3py.cui.load import load
    from phono3py.phonon3.collision_matrix import CollisionMatrix
    from phono3py.conductivity.velocity_solvers import GroupVelocitySolver
    from phono3py.conductivity.utils import get_unit_to_WmK
    from phonopy.phonon.grid import get_grid_point_from_address, get_grid_points_by_rotations
    from phonopy.physical_units import get_physical_units
    u = get_physical_units()

    inputs = common.verify_inputs()
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=common.SUPERCELL_FC3,
               phonon_supercell_matrix=common.SUPERCELL_FC2, fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"),
               is_nac=True, symmetrize_fc=True, make_r0_average=True, log_level=0, lang=a.lang)
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
    T = a.T
    cutoff = itr.cutoff_frequency
    grg2bzg = bz.grg2bzg
    bzg2grg = bz.bzg2grg
    addr = bz.addresses[grg2bzg]
    Jmap = np.array([get_grid_point_from_address(-ad, D) for ad in addr], dtype=int)
    f = freqs_bz[grg2bzg]                                   # (N, nb)
    xk = f * u.THzToEv / (u.KB * T)
    with np.errstate(over="ignore"):
        sh = np.where(f > cutoff, np.sinh(xk / 2), np.inf)   # sinh(x/2); inf excludes zero modes

    # symmetry maps on the GR grid (rotations incl. time reversal as provided by phono3py)
    rots = bz.rotations
    rot_maps = []
    imgs = np.array([get_grid_points_by_rotations(int(grg2bzg[g]), bz, reciprocal_rotations=rots)
                     for g in range(N)], dtype=np.int64)    # (N, nrot) BZ indices
    imgs_gr = bzg2grg[imgs]
    for k in range(imgs_gr.shape[1]):
        rot_maps.append(imgs_gr[:, k])
    rot_maps = np.array(rot_maps)
    has_tr = any(np.array_equal(m, Jmap) for m in rot_maps)

    gv_solver = GroupVelocitySolver(itr, is_kappa_star=False)
    gv = np.array([gv_solver.compute(int(grg2bzg[g])).group_velocities for g in range(N)])  # (N, nb, 3)

    col = CollisionMatrix(itr, rot_grid_points=None, lang=a.lang)
    ev_p, ev_a, ev_b, ev_g, ev_dl, ev_G = [], [], [], [], [], []
    gam_raw = np.zeros((N, nb))
    tol = a.nsig * a.sigma
    max_rep_freq_mismatch = 0.0
    for g0 in range(N):
        gp0 = int(grg2bzg[g0])
        col.set_grid_point(gp0)
        col.set_sigma(a.sigma, sigma_cutoff=None)
        col.run_integration_weights()
        col.run_interaction(is_full_pp=False)
        col.temperature = T
        col.run()
        gam_raw[g0] = np.array(col._imag_self_energy)
        pp = col._pp_strength
        ch0 = col._g[0]
        gp2tp, tp2s, swapped = col._get_gp2tp_map()
        conv = col._unit_conversion
        for g1 in range(N):
            ti = gp2tp[g1]
            # phono3py's tp2s[g1] is the third phonon of the (little-group) representative
            # triplet; only its frequency is rotation invariant. The event needs the actual
            # third phonon q2 = -q0 - q1 (mod G):
            g2 = int(get_grid_point_from_address(-(addr[g0] + addr[g1]), D))
            g2rep = int(bzg2grg[int(tp2s[g1])])
            fmis = float(np.abs(f[g2] - f[g2rep]).max())
            max_rep_freq_mismatch = max(max_rep_freq_mismatch, fmis)
            if swapped[g1]:
                P = pp[ti].transpose(0, 2, 1)
                W0 = ch0[ti].transpose(0, 2, 1)
            else:
                P = pp[ti]
                W0 = ch0[ti]
            # detuning in THz: w0 - w1 - w2
            dl = f[g0][:, None, None] - f[g1][None, :, None] - f[g2][None, None, :]
            E = conv * P * W0
            ok = (np.abs(dl) <= tol) & (E > 0)
            ok &= (f[g0][:, None, None] > cutoff) & (f[g1][None, :, None] > cutoff) & (f[g2][None, None, :] > cutoff)
            if not ok.any():
                continue
            j0, j1, l = np.nonzero(ok)
            pa = Jmap[g1] * nb + j1
            pb = Jmap[g2] * nb + l
            keep = pa <= pb
            j0, j1, l, pa, pb = j0[keep], j1[keep], l[keep], pa[keep], pb[keep]
            Ev = E[j0, j1, l]
            sp, sa, sb = sh[g0, j0], sh[Jmap[g1], j1], sh[Jmap[g2], l]
            rep = pa == pb
            gval = np.where(rep, Ev / (8 * sp * sa * sa), Ev / (4 * sp * sa * sb))
            # BZ-representative test for normal/umklapp (reporting only)
            a0 = bz.addresses[gp0]; a1 = bz.addresses[int(grg2bzg[g1])]; a2 = bz.addresses[int(grg2bzg[g2])]
            Gvec = a0 + a1 + a2
            ev_p.append(g0 * nb + j0); ev_a.append(pa); ev_b.append(pb); ev_g.append(gval)
            ev_dl.append(dl[j0, j1, l]); ev_G.append(np.repeat(int(np.any(Gvec != 0)), len(j0)))
        if g0 % 10 == 0:
            print(f"parent q {g0}/{N}  events so far {sum(len(x) for x in ev_p)}  {time.perf_counter()-t0:.0f}s", flush=True)
    ev_p = np.concatenate(ev_p); ev_a = np.concatenate(ev_a); ev_b = np.concatenate(ev_b)
    ev_g = np.concatenate(ev_g); ev_dl = np.concatenate(ev_dl); ev_G = np.concatenate(ev_G)

    np.savez_compressed(a.out, ev_p=ev_p, ev_a=ev_a, ev_b=ev_b, ev_g=ev_g, ev_dl=ev_dl, ev_umk=ev_G,
                        freqs=f, gv=gv, Jmap=Jmap, rot_maps=rot_maps, gam_raw=gam_raw, addr=addr,
                        mesh=np.array(a.mesh), sigma=a.sigma, nsig=a.nsig, T=T, cutoff=cutoff,
                        volume=ph3.primitive.volume, unit_to_WmK=get_unit_to_WmK(),
                        THzToEv=u.THzToEv, KB=u.KB, lattice=ph3.primitive.cell)
    meta = {"mesh": a.mesh, "sigma": a.sigma, "nsig": a.nsig, "T": T, "n_modes": N * nb,
            "n_events": int(len(ev_p)), "n_umklapp": int(ev_G.sum()), "n_repeated": int((ev_a == ev_b).sum()),
            "time_reversal_in_rotations": bool(has_tr), "n_rot": int(len(rot_maps)),
            "max_freq_mismatch_actual_vs_representative_q2_THz": max_rep_freq_mismatch,
            "inputs": inputs, "software": common.software_versions(), "seconds": time.perf_counter() - t0}
    Path(a.out).with_suffix(".json").write_text(json.dumps(meta, indent=1, default=str))
    print(json.dumps({k: v for k, v in meta.items() if k not in ("inputs", "software")}, indent=1))


if __name__ == "__main__":
    main()
