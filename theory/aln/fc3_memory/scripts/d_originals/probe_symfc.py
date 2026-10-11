"""Probe: fc3 layout of the hash-pinned AlN inputs, BZ-grid multiplicities, symfc fc3 basis sizes.

Run with the phono3py environment. Output: results/probe_symfc.json
usage: python -B probe_symfc.py [cutoffs ...]
"""
from __future__ import annotations

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

OUT = HERE.parent / "results"


def main():
    cutoffs = [None] + [float(x) for x in sys.argv[1:]]
    from phono3py.cui.load import load
    from symfc.basis_sets import FCBasisSetO3
    from symfc.utils.utils import SymfcAtoms

    inputs = common.verify_inputs()
    ph3 = load(unitcell_filename=str(common.DATASET / "POSCAR"), supercell_matrix=common.SUPERCELL_FC3,
               phonon_supercell_matrix=common.SUPERCELL_FC2, fc3_filename=str(common.DATASET / "fc3.hdf5"),
               fc2_filename=str(common.DATASET / "fc2.hdf5"), born_filename=str(common.DATASET / "BORN"),
               is_nac=True, symmetrize_fc=False, make_r0_average=True, log_level=0, lang="Rust")
    fc3 = np.asarray(ph3.fc3)
    sc = ph3.supercell
    prim = ph3.primitive
    out = {"inputs": inputs, "fc3_shape": list(fc3.shape), "fc3_norm": float(np.linalg.norm(fc3)),
           "fc3_absmax": float(np.abs(fc3).max()), "n_satom": len(sc), "p2s": prim.p2s_map.tolist(),
           "s2p_first20": prim.s2p_map[:20].tolist(), "masses": prim.masses.tolist(),
           "prim_scaled_positions": prim.scaled_positions.tolist(), "prim_cell": prim.cell.tolist()}
    svecs, multi = prim.get_smallest_vectors()
    out["svecs_shape"] = list(svecs.shape)
    out["multi_shape"] = list(multi.shape)
    out["multi_gt1_count"] = int((multi[:, :, 0] > 1).sum())
    # fraction of exactly-zero compact elements
    out["fc3_zero_fraction"] = float((fc3 == 0).mean())
    # permutation-symmetry and ASR drift of the compact fc3 (phono3py helper)
    try:
        from phono3py.phonon3.fc3 import get_drift_fc3
        out["drift_fc3"] = [float(x) for x in get_drift_fc3(fc3, primitive=prim)]
    except Exception as e:  # pragma: no cover
        out["drift_fc3_error"] = repr(e)

    ph3.mesh_numbers = [5, 5, 3]
    ph3.sigmas = [0.1]
    ph3.init_phph_interaction()
    itr = ph3.phph_interaction
    bz = itr.bz_grid
    out["bz_n_addresses"] = int(len(bz.addresses))
    out["gr_n"] = int(np.prod(bz.D_diag))
    out["all_shortest_count"] = int(itr.all_shortest.sum())
    out["unit_conversion_interaction"] = float(itr.unit_conversion_factor)

    satoms = SymfcAtoms(cell=sc.cell, scaled_positions=sc.scaled_positions, numbers=sc.numbers)
    res = []
    for c in cutoffs:
        t0 = time.time()
        bs = FCBasisSetO3(satoms, cutoff=c, use_mkl=False, log_level=0).run()
        t1 = time.time()
        Cn = bs._n_a_compression_matrix
        B = bs.basis_set
        rec = {"cutoff": c, "seconds": t1 - t0, "n_a_compr_shape": list(Cn.shape), "n_a_compr_nnz": int(Cn.nnz),
               "basis_shape": list(B.shape), "p2s_symfc": np.asarray(bs.p2s_map).tolist()}
        # orthonormality of the compression matrix columns
        G = (Cn.T @ Cn).toarray() if Cn.shape[1] < 20000 else None
        if G is not None:
            rec["CnTCn_dev_from_I"] = float(np.abs(G - np.eye(G.shape[0])).max())
        rec["BTB_dev_from_I"] = float(np.abs(B.T @ B - np.eye(B.shape[1])).max())
        # projection of the reference compact fc3 onto the symmetric (+ASR) span
        ccm = bs.compact_compression_matrix          # (n_a*N*N*27, n_compr)
        f = fc3[prim.p2s_map].reshape(-1) if fc3.shape[0] == fc3.shape[1] else fc3.reshape(-1)
        # least squares in compressed coordinates (columns of ccm are orthogonal with equal norms)
        colnorm2 = np.asarray(ccm.multiply(ccm).sum(axis=0)).ravel()
        rec["ccm_colnorm2_minmax"] = [float(colnorm2.min()), float(colnorm2.max())]
        eta = (ccm.T @ f) / colnorm2                 # projection onto symmetric (no ASR) span
        f_sym = ccm @ eta
        rec["rel_resid_sym_noASR"] = float(np.linalg.norm(f - f_sym) / np.linalg.norm(f))
        # with ASR: eta restricted to span(B) (B orthonormal in compressed coordinates); metric diag(colnorm2)
        # all colnorm2 equal => orthogonal projection is B B^T eta
        eta_asr = B @ (B.T @ eta)
        f_asr = ccm @ eta_asr
        rec["rel_resid_sym_ASR"] = float(np.linalg.norm(f - f_asr) / np.linalg.norm(f))
        res.append(rec)
        print(json.dumps(rec), flush=True)
    out["symfc"] = res
    common.write_json(OUT / "probe_symfc.json", out)


if __name__ == "__main__":
    main()
