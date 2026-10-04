"""Reproducible phono3py 4.5.0 driver for wurtzite AlN lattice thermal conductivity.

Inputs: hash-pinned Phonon Olympics files (POSCAR, BORN, fc2.hdf5 [5x5x3 supercell],
fc3.hdf5 [3x3x2 supercell]); see common.py.  Each invocation creates a NEW directory
<campaign>/runs/<run_id>/ (refuses to overwrite), writes phono3py outputs there and a
run.json record (settings, versions, input hashes, grid sizes, timings, peak memory,
kappa tensors, output-file hashes).

Version-compatibility switches (phono3py 2.1.0 was used by the Phonon Olympics team):
  --no-r0-average     make_r0_average=False ("rough backward compatibility with v2.x")
  --gv-delta-q 1e-5   finite-difference group velocities (v4.0.x and earlier with NAC)

Examples
--------
  python run_kappa.py --run-id rta-m9x9x5 --mesh 9 9 5 --method rta --temps 300
  python run_kappa.py --run-id lbte-m9x9x5 --mesh 9 9 5 --method lbte --temps 300 \
      --write-collision --dump-assembled
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-id", required=True)
    p.add_argument("--mesh", type=int, nargs=3, required=True)
    p.add_argument("--method", choices=["rta", "lbte"], required=True)
    p.add_argument("--sigmas", nargs="+", default=["tetra"],
                   help="'tetra' (linear tetrahedron) or Gaussian widths in THz. "
                        "Do not mix 'tetra' with Gaussian widths in one run.")
    p.add_argument("--sigma-cutoff", type=float, default=None)
    p.add_argument("--temps", type=float, nargs="+", default=[300.0])
    p.add_argument("--isotope", action="store_true", help="natural-abundance isotope scattering")
    p.add_argument("--gv-delta-q", type=float, default=None,
                   help="finite-difference group velocity step (1/Angstrom); default analytic")
    p.add_argument("--lang", choices=["C", "Rust"], default="Rust")
    p.add_argument("--is-N-U", action="store_true", help="RTA only: split gamma into N and U")
    p.add_argument("--write-gamma", action="store_true", help="RTA only: per-grid-point files")
    p.add_argument("--write-collision", action="store_true",
                   help="LBTE only: per-grid-point raw collision rows (before assembly)")
    p.add_argument("--dump-assembled", action="store_true",
                   help="LBTE only: save the assembled (weighted, degeneracy-averaged, "
                        "symmetrized) collision matrix just before diagonalization")
    p.add_argument("--reducible", action="store_true", help="LBTE only: full-grid collision matrix")
    p.add_argument("--pinv-solver", type=int, default=0)
    p.add_argument("--pinv-cutoff", type=float, default=None)
    p.add_argument("--is-full-pp", action="store_true")
    p.add_argument("--symmetrize-fc", action="store_true",
                   help="apply phono3py symmetrization to the fc2/fc3 read from file")
    p.add_argument("--nac-gcutoff-factor", type=float, default=None,
                   help="multiply the Gonze-Lee reciprocal cutoff G_cutoff (Lambda fixed)")
    p.add_argument("--no-nac", action="store_true")
    p.add_argument("--no-r0-average", action="store_true",
                   help="make_r0_average=False (phono3py v2.x-like fc3 real-to-reciprocal phase)")
    p.add_argument("--log-level", type=int, default=1)
    p.add_argument("--note", default="")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    timer = common.Timer()
    run_dir = common.CAMPAIGN / "runs" / args.run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    os.chdir(run_dir)

    import numpy as np
    import phono3py
    from phono3py.cui.load import load

    assert phono3py.__version__ == "4.5.0", phono3py.__version__
    inputs = common.verify_inputs()
    record = {
        "run_id": args.run_id,
        "argv": sys.argv,
        "note": args.note,
        "campaign_dir": str(common.CAMPAIGN),
        "dataset_dir": str(common.DATASET),
        "script": {"path": str(Path(__file__).resolve()), "sha256": common.sha256(__file__),
                   "common_sha256": common.sha256(Path(__file__).with_name("common.py"))},
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "software": common.software_versions(),
        "inputs": inputs,
        "settings": vars(args).copy(),
        "status": "running",
    }
    common.write_json(run_dir / "run.json", record)

    sigmas = [None if s.lower() in ("tetra", "none", "thm") else float(s) for s in args.sigmas]
    if None in sigmas and len(sigmas) > 1:
        raise SystemExit("Mixing tetrahedron with Gaussian widths in one run reuses pp "
                         "computed with tetrahedron zero-skipping; run them separately.")

    def _load(log_level):
        return load(unitcell_filename=str(common.DATASET / "POSCAR"),
                    supercell_matrix=common.SUPERCELL_FC3,
                    phonon_supercell_matrix=common.SUPERCELL_FC2,
                    fc3_filename=str(common.DATASET / "fc3.hdf5"),
                    fc2_filename=str(common.DATASET / "fc2.hdf5"),
                    born_filename=None if args.no_nac else str(common.DATASET / "BORN"),
                    is_nac=not args.no_nac,
                    symmetrize_fc=False,
                    make_r0_average=not args.no_r0_average,
                    log_level=log_level,
                    lang=args.lang)

    ph3 = _load(args.log_level)
    if args.nac_gcutoff_factor is not None:
        # Read default Gonze-Lee parameters, then rebuild with an enlarged reciprocal
        # cutoff at fixed Lambda (same control as the interaction pilot).
        ph3.mesh_numbers = [1, 1, 1]
        ph3.init_phph_interaction()
        dm = ph3.dynamical_matrix
        dm.run(np.array([1.0, 0.0, 0.0]) / 3)
        _, _, g_cut, g_list, lam = dm.Gonze_nac_dataset
        record["nac_default"] = {"G_cutoff": float(g_cut), "Lambda": float(lam), "G_count": int(len(g_list))}
        ph3 = _load(args.log_level)
        ph3.nac_params = dict(ph3.nac_params, G_cutoff=float(g_cut) * args.nac_gcutoff_factor,
                              Lambda=float(lam))
    if args.symmetrize_fc:
        ph3.symmetrize_fc3()
        ph3.symmetrize_fc2()
        record["settings"]["fc_symmetrized_after_read"] = True
    nac = ph3.nac_params
    record["nac_params"] = None if nac is None else {
        "born": np.asarray(nac["born"]).tolist(),
        "dielectric": np.asarray(nac["dielectric"]).tolist(),
        "factor": nac.get("factor"), "method": nac.get("method"),
        "G_cutoff": nac.get("G_cutoff"), "Lambda": nac.get("Lambda")}
    record["masses"] = np.asarray(ph3.primitive.masses).tolist()
    record["primitive_volume_A3"] = float(ph3.primitive.volume)
    timer.mark("load")

    ph3.mesh_numbers = args.mesh
    ph3.sigmas = sigmas
    if args.sigma_cutoff is not None:
        ph3.sigma_cutoff = args.sigma_cutoff
    ph3.init_phph_interaction()
    if not args.no_nac:
        dm = ph3.dynamical_matrix
        dm.run(np.array([1.0, 0.0, 0.0]) / 3)
        try:
            _, _, g_cut, g_list, lam = dm.Gonze_nac_dataset
            record["nac_effective"] = {"G_cutoff": float(g_cut), "Lambda": float(lam),
                                       "G_count": int(len(g_list)), "class": type(dm).__name__}
        except Exception as exc:  # noqa: BLE001
            record["nac_effective"] = {"error": repr(exc), "class": type(dm).__name__}
    bz = ph3.grid
    from phonopy.phonon.grid import get_ir_grid_points
    ir_grg, ir_w, _ = get_ir_grid_points(bz)
    record["grid"] = {"mesh": list(map(int, bz.D_diag)), "num_mesh": int(np.prod(bz.D_diag)),
                      "num_ir_grid_points": int(len(ir_grg)), "num_bz_grid_points": int(len(bz.addresses)),
                      "num_rotations": int(len(bz.rotations))}
    print("Grid:", record["grid"], flush=True)
    timer.mark("init_interaction")
    ph3.run_phonon_solver()
    freqs, _, _ = ph3.get_phonon_data()
    record["frequency_THz"] = {"min": float(freqs.min()), "max": float(freqs.max())}
    timer.mark("phonons")
    record["memory_after_setup"] = common.memory_info()
    common.write_json(run_dir / "run.json", record)

    dumped = []
    if args.method == "lbte" and args.dump_assembled:
        import h5py
        import phono3py.conductivity.collision_matrix_kernel as cmk
        _orig = cmk.diagonalize_collision_matrix
        _orig_weights = cmk.IrreducibleCollisionMatrixKernel._get_collision_weights
        captured = {}

        def _weights_hook(self):
            w = _orig_weights(self)
            captured["weights"] = np.array(w)
            captured["rotations_cartesian"] = np.array(self._rotations_cartesian)
            captured["rot_grid_points"] = np.array(self._rot_grid_points)
            return w

        def _diag_hook(collision_matrices, i_sigma=None, i_temp=None, pinv_solver=0, log_level=0):
            fn = run_dir / f"assembled-colmat-s{i_sigma}-t{i_temp}.hdf5"
            t0 = time.perf_counter()
            with h5py.File(fn, "w") as f:
                f.create_dataset("collision_matrix", data=collision_matrices[i_sigma, i_temp])
                for k, v in captured.items():
                    f.create_dataset(k, data=v)
                f.attrs["i_sigma"] = i_sigma
                f.attrs["i_temp"] = i_temp
                f.attrs["description"] = (
                    "phono3py 4.5.0 assembled collision matrix exactly as passed to the "
                    "eigensolver: raw rows + gamma*R on stabilizer diagonal, multiplied by "
                    "w_i*w_j (w=sqrt(star size/|G|)), degeneracy-averaged rows and columns, "
                    "then (A+A^T)/2. Units THz (ordinary frequency, half-linewidth convention).")
            dumped.append({"file": fn.name, "seconds": time.perf_counter() - t0})
            print(f"[dump] wrote {fn.name} in {time.perf_counter() - t0:.1f}s", flush=True)
            return _orig(collision_matrices, i_sigma=i_sigma, i_temp=i_temp,
                         pinv_solver=pinv_solver, log_level=log_level)

        cmk.diagonalize_collision_matrix = _diag_hook
        cmk.IrreducibleCollisionMatrixKernel._get_collision_weights = _weights_hook

    kw = dict(temperatures=args.temps, is_isotope=args.isotope, write_kappa=True,
              gv_delta_q=args.gv_delta_q, is_full_pp=args.is_full_pp, log_level=args.log_level)
    if args.method == "rta":
        ph3.run_thermal_conductivity(is_LBTE=False, is_N_U=args.is_N_U,
                                     write_gamma=args.write_gamma, **kw)
    else:
        ph3.run_thermal_conductivity(is_LBTE=True, write_collision=args.write_collision,
                                     is_reducible_collision_matrix=args.reducible,
                                     pinv_solver=args.pinv_solver, pinv_cutoff=args.pinv_cutoff,
                                     **kw)
    timer.mark("kappa")
    tc = ph3.thermal_conductivity
    kappa = np.asarray(tc.kappa)
    res = {"temperatures": np.asarray(tc.temperatures).tolist(),
           "sigmas": [s if s is not None else "tetra" for s in sigmas],
           "kappa_WmK[sigma][T][xx,yy,zz,yz,xz,xy]": kappa.tolist()}
    if args.method == "lbte":
        res["kappa_RTA_WmK[sigma][T][xx,yy,zz,yz,xz,xy]"] = np.asarray(tc.kappa_RTA).tolist()
        ev = tc.collision_eigenvalues
        if ev is not None:
            ev = np.asarray(ev)
            res["collision_eigenvalues_summary"] = [[{
                "min": float(ev[i, j].min()), "max": float(ev[i, j].max()),
                "num_abs_below_1e-8": int((np.abs(ev[i, j]) < 1e-8).sum()),
                "num_negative_below_-1e-8": int((ev[i, j] < -1e-8).sum()),
                "size": int(ev.shape[-1])} for j in range(ev.shape[1])] for i in range(ev.shape[0])]
    record["result"] = res
    record["timing_s"] = timer.marks
    record["memory_end"] = common.memory_info()
    record["dumped_assembled"] = dumped
    outs = {}
    for pth in sorted(run_dir.glob("*.hdf5")):
        outs[pth.name] = {"bytes": pth.stat().st_size, "sha256": common.sha256(pth)}
    record["outputs"] = outs
    record["status"] = "finished"
    record["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    common.write_json(run_dir / "run.json", record)
    print("RESULT", res, flush=True)
    print("TIMING", timer.marks, flush=True)
    print("MEMORY", record["memory_end"], flush=True)


if __name__ == "__main__":
    main()
