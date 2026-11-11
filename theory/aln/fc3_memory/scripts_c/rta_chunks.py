"""E27: resumable full-mesh RTA (phono3py 4.5.0, 300 K, linear tetrahedron, NAC as the transport baseline) in chunks of
irreducible grid points. Each chunk is written as chunk_<i>.npz (grid points, weights, q-points, frequencies, gamma,
group velocities, mode_kappa); existing chunks are skipped, so the run survives interruptions. 'assemble' combines
the chunks: kappa, tau_mem(RTA) = sum mode_kappa tau / sum mode_kappa (tau = 1/(4 pi gamma)), the part above
frequency cutoffs, and the film-regularised tau_mem(L) with 1/tau_L = 1/tau + |v_z|/L (as ir_film.py).
Checked on 9x9x5: chunked kappa 231.1818 / 212.4176 W/mK equals the single-run value.
usage: python -B rta_chunks.py run --mesh N1 N2 N3 --outdir DIR [--chunk 64] [--threads 12]
       python -B rta_chunks.py assemble --outdir DIR
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("mode", choices=["run", "assemble"])
ap.add_argument("--mesh", type=int, nargs=3)
ap.add_argument("--outdir", required=True)
ap.add_argument("--chunk", type=int, default=64)
ap.add_argument("--threads", type=int, default=12)
a = ap.parse_args()
os.environ["OMP_NUM_THREADS"] = str(a.threads)
os.environ["RAYON_NUM_THREADS"] = str(a.threads)
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

out = Path(a.outdir)
out.mkdir(parents=True, exist_ok=True)

if a.mode == "run":
    import phono3py
    from phono3py.cui.load import load
    from phonopy.phonon.grid import get_ir_grid_points

    import aln_common as common
    assert phono3py.__version__ == "4.5.0"
    inputs = common.verify_inputs()
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=common.SUPERCELL_FC3,
               phonon_supercell_matrix=common.SUPERCELL_FC2, fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"), is_nac=True,
               symmetrize_fc=False, make_r0_average=True, log_level=0, lang="Rust")
    ph3.mesh_numbers = a.mesh
    ph3.sigmas = [None]
    ph3.init_phph_interaction()
    ph3.run_phonon_solver()
    bz = ph3.grid
    ir, w, _ = get_ir_grid_points(bz)
    gps = np.asarray(bz.grg2bzg[ir])
    nch = (len(gps) + a.chunk - 1) // a.chunk
    meta = {"mesh": a.mesh, "n_ir": int(len(gps)), "sum_weights": int(w.sum()), "chunk": a.chunk, "n_chunks": nch,
            "inputs": inputs, "phono3py": phono3py.__version__}
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    for i in range(nch):
        f = out / f"chunk_{i:04d}.npz"
        if f.exists():
            continue
        t0 = time.perf_counter()
        sel = gps[i * a.chunk:(i + 1) * a.chunk]
        ph3.run_thermal_conductivity(temperatures=[300.0], grid_points=sel, write_kappa=False, log_level=0)
        tc = ph3.thermal_conductivity
        tmp = out / f"chunk_{i:04d}.tmp.npz"
        np.savez(tmp, gp=np.asarray(tc.grid_points), w=np.asarray(tc.grid_weights), q=np.asarray(tc.qpoints),
                 f=np.asarray(tc.frequencies), gamma=np.asarray(tc.gamma)[0, 0], gv=np.asarray(tc.group_velocities),
                 mk=np.asarray(tc.mode_kappa)[0, 0])
        os.replace(tmp, f)
        print(json.dumps({"chunk": i, "of": nch, "seconds": time.perf_counter() - t0}), flush=True)
    (out / "done").write_text("done")
else:
    meta = json.loads((out / "meta.json").read_text())
    fs = sorted(out.glob("chunk_*.npz"))
    fs = [f for f in fs if ".tmp" not in f.name]
    assert len(fs) == meta["n_chunks"], (len(fs), meta["n_chunks"])
    D = {k: np.concatenate([np.load(f)[k] for f in fs]) for k in ("w", "f", "gamma", "gv", "mk")}
    assert int(D["w"].sum()) == meta["sum_weights"]
    ok = (D["gamma"] > 0) & (D["f"] > 1e-3)
    tau = np.where(ok, 1 / (4 * np.pi * np.where(ok, D["gamma"], 1)), 0.0)
    vz = np.abs(D["gv"][..., 2]) * 100.0
    res = {"mesh": meta["mesh"], "f_min_THz": float(D["f"][D["f"] > 1e-3].min())}
    nq = int(np.prod(meta["mesh"]))
    for c, nm in ((0, "x"), (2, "z")):
        k = np.where(ok, D["mk"][..., c], 0.0)
        K, N = k.sum(), (k * tau).sum()
        res[f"kappa_sum_{nm}"] = float(K)
        res[f"kappa_{nm}"] = float(K / meta["sum_weights"])          # W/mK (sum of mode_kappa / N_q)
        res[f"tau_mem_{nm}"] = float(N / K)
        res[f"tau_above_{nm}"] = {str(fc): float((k * tau)[D["f"] > fc].sum() / K) for fc in (1.5, 2.0, 3.0, 4.0)}
        wgt = np.where(ok, k / np.where(ok, tau, 1), 0.0)
        res[f"film_vz_{nm}"] = {}
        for L in (30, 100, 300, 1000, 3000, 10000):
            tl = np.where(ok, 1 / (1 / np.where(ok, tau, 1) + vz / (L * 1e-9) * 1e-12), 0.0)
            res[f"film_vz_{nm}"][str(L)] = float((wgt * tl ** 2).sum() / (wgt * tl).sum())
    (out / "assembled.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
