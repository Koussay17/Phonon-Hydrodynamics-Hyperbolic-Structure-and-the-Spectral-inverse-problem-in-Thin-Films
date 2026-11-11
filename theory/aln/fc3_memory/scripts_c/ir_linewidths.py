"""E25 (T1): three-phonon linewidths of long-wavelength acoustic modes at fixed q, converged in the integration mesh.

For each target q = address / mesh (addresses given on a reference division so that the same physical q is on every
integration mesh), phono3py 4.5.0 RTA linewidths gamma(q, band) at 300 K (half width, THz) with the linear
tetrahedron method (default) or a Gaussian sigma, NAC as in the transport baseline. Records frequencies, group
velocities and gamma of all 12 bands, so the low-frequency law gamma ~ f^n of each acoustic branch and its
convergence in the integration mesh can be read off.
usage: python -B ir_linewidths.py --mesh N1 N2 N3 --q a/b,c/d,e/f [--q ...] --out out.json [--sigma s] [--threads n]
       [--symmetrize-fc]
       q components as fractions, e.g. --q 1/48,0,0 --q 0,0,1/24 ; each must lie on the mesh.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--mesh", type=int, nargs=3, required=True)
ap.add_argument("--q", action="append", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--sigma", type=float, default=None)
ap.add_argument("--threads", type=int, default=12)
ap.add_argument("--symmetrize-fc", action="store_true",
                help="impose translational (sum-rule) and permutation symmetry on fc2/fc3 after reading (phono3py)")
args = ap.parse_args()
os.environ["OMP_NUM_THREADS"] = str(args.threads)
os.environ["RAYON_NUM_THREADS"] = str(args.threads)
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import phono3py  # noqa: E402
from phono3py.cui.load import load  # noqa: E402

import aln_common as common  # noqa: E402

t0 = time.perf_counter()
assert phono3py.__version__ == "4.5.0"
inputs = common.verify_inputs()
mesh = np.array(args.mesh)
targets, addrs = [], []
for qs in args.q:
    fr = [Fraction(x) for x in qs.split(",")]
    a = [f * int(m) for f, m in zip(fr, mesh)]
    if any(x.denominator != 1 for x in a):
        raise SystemExit(f"q {qs} is not on mesh {mesh.tolist()}")
    targets.append(qs)
    addrs.append([int(x) for x in a])
ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=common.SUPERCELL_FC3,
           phonon_supercell_matrix=common.SUPERCELL_FC2, fc3_filename=str(common.DATASET / "fc3.hdf5"),
           fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"), is_nac=True,
           symmetrize_fc=False, make_r0_average=True, log_level=0, lang="Rust")
if args.symmetrize_fc:
    ph3.symmetrize_fc3()
    ph3.symmetrize_fc2()
ph3.mesh_numbers = mesh.tolist()
ph3.sigmas = [args.sigma]
ph3.init_phph_interaction()
bz = ph3.grid
gps = np.array(bz.get_indices_from_addresses(np.array(addrs, dtype=np.int64)), dtype=np.int64).ravel()
ph3.run_phonon_solver()
t1 = time.perf_counter()
ph3.run_thermal_conductivity(temperatures=[300.0], grid_points=gps, write_kappa=False, log_level=0)
tc = ph3.thermal_conductivity
gam = np.asarray(tc.gamma)[0, 0]
fq = np.asarray(tc.frequencies)
gv = np.asarray(tc.group_velocities)
gp_out = np.asarray(tc.grid_points)
recs = []
for i, (qs, a) in enumerate(zip(targets, addrs)):
    j = int(np.nonzero(gp_out == gps[i])[0][0])
    recs.append({"q": qs, "address": a, "grid_point": int(gps[i]),
                 "bz_address": np.asarray(bz.addresses[gps[i]]).tolist(),
                 "frequency_THz": fq[j].tolist(), "gamma_THz": gam[j].tolist(),
                 "group_velocity_THzA": gv[j].tolist()})
out = {"mesh": mesh.tolist(), "symmetrize_fc": args.symmetrize_fc, "sigma": args.sigma, "T": 300.0, "inputs": inputs,
       "software": {"phono3py": phono3py.__version__}, "threads": args.threads,
       "seconds_setup": t1 - t0, "seconds_kappa": time.perf_counter() - t1, "points": recs}
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
Path(args.out).write_text(json.dumps(out, indent=1))
print(json.dumps({"mesh": mesh.tolist(), "n": len(recs), "seconds_kappa": out["seconds_kappa"]}))
for r in recs:
    print(r["q"], "f", np.round(r["frequency_THz"][:3], 4).tolist(), "gamma", ["%.3e" % g for g in r["gamma_THz"][:3]])
