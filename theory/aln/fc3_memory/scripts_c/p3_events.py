"""phono3py-route evaluation of |Phi_alpha|^2 on a FIXED three-phonon event list (C: mirror copy; identical
logic to the D: original). Run with a phono3py 4.5.0 environment.

For every event of the memory-gap list (p -> a + b, a = (-q1, j1), b = (-q2, l), canonical a <= b) the parent's
CollisionMatrix machinery is re-run (as aln_events.py) and the script records
    tp (BZ indices of phono3py's representative triplet), bands (in that triplet's order),
    pp (interaction strength |Phi|^2 * Interaction.unit_conversion), g0w (Gaussian weight, 1/THz),
so that g = conv * pp * g0w / (4 s_p s_a s_b) (distinct daughters) or / (8 s_p s_a^2) (repeated daughter).
--static also writes phonons on the BZ grid, smallest vectors, all_shortest, supercell.
usage: python -B p3_events.py --mesh 5 5 3 --events EVENTS.npz --outdir DIR [--r0avg 1] [--static]
                              [--fc3 REF | --fc3 name=path.npy ...] [--tag TAG]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("RAYON_NUM_THREADS", "4")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import aln_common as common  # noqa: E402


def load_ph3(r0avg: bool, mesh, sigma, lang="Rust"):
    from phono3py.cui.load import load
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=common.SUPERCELL_FC3,
               phonon_supercell_matrix=common.SUPERCELL_FC2, fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"),
               is_nac=True, symmetrize_fc=False, make_r0_average=bool(r0avg), log_level=0, lang=lang)
    ph3.mesh_numbers = mesh
    ph3.sigmas = [sigma]
    return ph3


def event_rows(ph3, EV, sigma, T, lang="Rust", want_static=False):
    from phono3py.phonon3.collision_matrix import CollisionMatrix
    ph3.init_phph_interaction()
    ph3.run_phonon_solver()
    itr = ph3.phph_interaction
    bz = itr.bz_grid
    D = bz.D_diag
    N = int(np.prod(D))
    freqs_bz, _, _ = itr.get_phonons()
    nb = freqs_bz.shape[1]
    grg2bzg = bz.grg2bzg
    Jmap = EV["Jmap"]
    ev_p, ev_a, ev_b = EV["ev_p"], EV["ev_a"], EV["ev_b"]
    m = len(ev_p)
    g0_of, j0_of = ev_p // nb, ev_p % nb
    g1_of, j1_of = Jmap[ev_a // nb], ev_a % nb
    l_of = ev_b % nb
    tp = -np.ones((m, 3), dtype=np.int64)
    bands = -np.ones((m, 3), dtype=np.int64)
    pp_out = np.full(m, np.nan)
    g0w = np.full(m, np.nan)
    order = np.argsort(g0_of, kind="stable")
    bounds = np.searchsorted(g0_of[order], np.arange(N + 1))
    col = CollisionMatrix(itr, rot_grid_points=None, lang=lang)
    conv = None
    for g0 in range(N):
        idx = order[bounds[g0]:bounds[g0 + 1]]
        col.set_grid_point(int(grg2bzg[g0]))
        col.set_sigma(sigma, sigma_cutoff=None)
        col.run_integration_weights()
        col.run_interaction(is_full_pp=False)
        col.temperature = T
        conv = col._unit_conversion
        if len(idx) == 0:
            continue
        pp = col._pp_strength
        ch0 = col._g[0]
        gp2tp, tp2s, swapped = col._get_gp2tp_map()
        trip, _, _, _ = itr.get_triplets_at_q()
        g1 = g1_of[idx]
        ti = gp2tp[g1]
        assert (ti >= 0).all()
        sw = swapped[g1].astype(bool)
        j0 = j0_of[idx]
        b1 = np.where(sw, l_of[idx], j1_of[idx])
        b2 = np.where(sw, j1_of[idx], l_of[idx])
        tp[idx] = trip[ti]
        bands[idx, 0], bands[idx, 1], bands[idx, 2] = j0, b1, b2
        pp_out[idx] = pp[ti, j0, b1, b2]
        g0w[idx] = ch0[ti, j0, b1, b2]
    out = {"tp": tp, "bands": bands, "pp": pp_out, "g0w": g0w, "conv": float(conv),
           "unit_conversion": float(itr.unit_conversion_factor)}
    if want_static:
        prim = ph3.primitive
        svecs, multi = prim.get_smallest_vectors()
        fz, ez, _ = itr.get_phonons()
        out["static"] = dict(freqs_bz=np.array(fz), eigv_bz=np.array(ez), bz_addresses=np.array(bz.addresses),
                             D_diag=np.array(D), grg2bzg=np.array(grg2bzg), bzg2grg=np.array(bz.bzg2grg),
                             gp_Gamma=int(bz.gp_Gamma), svecs=np.array(svecs), multi=np.array(multi),
                             p2s=np.array(prim.p2s_map), s2p=np.array(prim.s2p_map), masses=np.array(prim.masses),
                             prim_scaled=np.array(prim.scaled_positions), prim_cell=np.array(prim.cell),
                             all_shortest=np.array(itr.all_shortest), cutoff_frequency=float(itr.cutoff_frequency),
                             sc_cell=np.array(ph3.supercell.cell), sc_scaled=np.array(ph3.supercell.scaled_positions),
                             sc_numbers=np.array(ph3.supercell.numbers))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=int, nargs=3, required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--r0avg", type=int, default=1)
    ap.add_argument("--static", action="store_true")
    ap.add_argument("--fc3", action="append", default=None)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = common.verify_inputs()
    EV = np.load(a.events)
    ph3 = load_ph3(a.r0avg, a.mesh, a.sigma)
    fc3_ref = np.array(ph3.fc3)
    log = {"inputs": inputs, "software": common.software_versions(), "args": vars(a), "runs": []}
    tag = f"{a.tag}_" if a.tag else ""
    for k, spec in enumerate(a.fc3 or ["REF"]):
        t0 = time.perf_counter()
        if spec == "REF":
            name, fc3 = "REF", fc3_ref
        else:
            name, path = spec.split("=", 1)
            fc3 = np.load(path)
        ph3.fc3 = np.ascontiguousarray(fc3, dtype="double")
        res = event_rows(ph3, EV, a.sigma, a.T, want_static=(a.static and k == 0))
        static = res.pop("static", None)
        fn = outdir / f"pp_{tag}r0{a.r0avg}_{name}.npz"
        np.savez_compressed(fn, **res)
        if static is not None:
            np.savez_compressed(outdir / f"static_{tag}r0{a.r0avg}.npz", **static)
        rec = {"name": name, "file": str(fn), "seconds": time.perf_counter() - t0}
        log["runs"].append(rec)
        print(json.dumps(rec), flush=True)
    common.write_json(outdir / f"p3_events_{tag}r0{a.r0avg}.json", log)


if __name__ == "__main__":
    main()
