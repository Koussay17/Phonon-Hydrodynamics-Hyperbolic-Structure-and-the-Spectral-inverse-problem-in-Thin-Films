"""Localize the C-vs-Rust backend difference in phono3py 4.5.0 (AlN, 9x9x5 mesh).

Observed: RTA kappa differs by 1.3e-4 (xx) / 3.5e-4 (zz) between lang='C' and
lang='Rust' with identical settings (FD group velocities), individual linewidths
by up to 1.6 %.  This script separates (1) phonon solver, (2) ph-ph interaction
strength, (3) tetrahedron integration weights, (4) linewidth accumulation, by
feeding identical inputs to both backends.

usage: python diag_c_vs_rust.py [--mesh 9 9 5] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from phono3py.cui.load import load  # noqa: E402
from phono3py.phonon3.triplets import get_triplets_integration_weights  # noqa: E402
from phonopy.phonon.grid import get_ir_grid_points  # noqa: E402


def setup(lang, mesh, r0avg=True):
    ph = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=[3, 3, 2],
              phonon_supercell_matrix=[5, 5, 3], fc3_filename=str(common.DATASET / "fc3.hdf5"),
              fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"),
              is_nac=True, symmetrize_fc=False, make_r0_average=r0avg, log_level=0, lang=lang)
    ph.mesh_numbers = mesh
    ph.init_phph_interaction()
    ph.run_phonon_solver()
    return ph


def degenerate_projector_diff(fa, ea, fb, eb, tol=1e-4):
    """max || P_a - P_b || over degenerate blocks at each grid point (eigvecs columns)."""
    worst = 0.0
    for gp in range(fa.shape[0]):
        f = fa[gp]
        i = 0
        while i < len(f):
            j = i + 1
            while j < len(f) and abs(f[j] - f[i]) < tol:
                j += 1
            A = ea[gp][:, i:j]
            B = eb[gp][:, i:j]
            worst = max(worst, float(np.abs(A @ A.conj().T - B @ B.conj().T).max()))
            i = j
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=int, nargs=3, default=[9, 9, 5])
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    out = {"mesh": a.mesh}
    phC = setup("C", a.mesh)
    phR = setup("Rust", a.mesh)
    fC, eC, addrC = phC.get_phonon_data()
    fR, eR, addrR = phR.get_phonon_data()
    assert np.array_equal(addrC, addrR)
    out["phonon_freq_maxabs_diff_THz"] = float(np.abs(fC - fR).max())
    out["phonon_degenerate_projector_maxdiff"] = degenerate_projector_diff(fC, eC, fR, eR)
    # Inject C phonons into the Rust object so that both use identical phonons.
    phR.phph_interaction.set_phonon_data(fC, eC, addrC)
    itC, itR = phC.phph_interaction, phR.phph_interaction
    bz = itC.bz_grid
    ir_grg, _, _ = get_ir_grid_points(bz)
    rows = []
    for i_ir in range(len(ir_grg)):
        gp = int(bz.grg2bzg[ir_grg[i_ir]])
        itC.set_grid_point(gp)
        itR.set_grid_point(gp)
        tC, wC = itC.get_triplets_at_q()[:2]
        tR, wR = itR.get_triplets_at_q()[:2]
        same_triplets = bool(np.array_equal(tC, tR) and np.array_equal(wC, wR))
        itC.run()
        itR.run()
        ppC, ppR = itC.interaction_strength, itR.interaction_strength
        fpts = fC[gp]
        gC, _ = get_triplets_integration_weights(itC, fpts, None, lang="C")
        gR, _ = get_triplets_integration_weights(itC, fpts, None, lang="Rust")
        gP, _ = get_triplets_integration_weights(itC, fpts, None, lang="Python")
        rows.append({
            "ir_index": i_ir, "bz_gp": gp, "same_triplets": same_triplets,
            "pp_maxabs_diff_over_max": float(np.abs(ppC - ppR).max() / np.abs(ppC).max()),
            "g_C_vs_Rust_maxabs": float(np.abs(gC - gR).max()),
            "g_C_vs_Py_maxabs": float(np.abs(gC - gP).max()),
            "g_Rust_vs_Py_maxabs": float(np.abs(gR - gP).max()),
            "g_max": float(np.abs(gC).max()),
        })
    out["per_ir"] = rows
    out["summary"] = {k: max(r[k] for r in rows) for k in rows[0] if k.endswith(("_over_max", "_maxabs"))}
    out["all_same_triplets"] = all(r["same_triplets"] for r in rows)
    print(json.dumps(out["summary"], indent=2))
    print("phonon freq diff", out["phonon_freq_maxabs_diff_THz"], "projector diff", out["phonon_degenerate_projector_maxdiff"])
    worst = sorted(rows, key=lambda r: -r["g_C_vs_Rust_maxabs"])[:5]
    print(json.dumps(worst, indent=1))
    if a.json:
        common.write_json(a.json, out)


if __name__ == "__main__":
    main()
