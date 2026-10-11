"""Symmetric fc3 space (symfc 1.7.3, 3x3x2 supercell) and its cluster-orbit structure. phono3py environment.

compact fc3 = ccm @ eta, eta = B @ c (ccm: compact compression matrix, columns orthogonal with |col|^2 = 1/18;
B orthonormal; sum rule included in B). Saves basis_ccm.npz, basis_B.npy, basis_c0.npy (coordinates of the
reference fc3), fc3_ref_compact.npy, cluster_orbits.npz (column -> cluster-orbit block, block reference norm and
largest pair distance; each column lies in one distance shell).
usage: python -B basis_p3.py <outdir>
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import aln_common as common  # noqa: E402
from p3_events import load_ph3  # noqa: E402


def symfc_basis(ph3, cutoff=None):
    from symfc.basis_sets import FCBasisSetO3
    from symfc.utils.utils import SymfcAtoms
    sc = ph3.supercell
    bs = FCBasisSetO3(SymfcAtoms(cell=sc.cell, scaled_positions=sc.scaled_positions, numbers=sc.numbers),
                      cutoff=cutoff, use_mkl=False, log_level=0).run()
    assert list(np.asarray(bs.p2s_map)) == list(ph3.primitive.p2s_map)
    return sp.csr_matrix(bs.compact_compression_matrix), np.asarray(bs.basis_set)


def pair_distances(cell, pos):
    d = np.full((len(pos), len(pos)), np.inf)
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                diff = pos[None, :, :] - pos[:, None, :] + np.array([i, j, k])[None, None, :]
                d = np.minimum(d, np.linalg.norm(diff @ cell, axis=2))
    return d


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    inputs = common.verify_inputs()
    ph3 = load_ph3(True, [5, 5, 3], 0.1)
    p2s = ph3.primitive.p2s_map
    f0 = np.array(ph3.fc3)[p2s].copy()
    ccm, B = symfc_basis(ph3, None)
    eta0 = 18.0 * (ccm.T @ f0.ravel())
    c0 = B.T @ eta0
    resid = float(np.linalg.norm(ccm @ (B @ c0) - f0.ravel()) / np.linalg.norm(f0))
    sp.save_npz(out / "basis_ccm.npz", ccm)
    np.save(out / "basis_B.npy", B)
    np.save(out / "basis_c0.npy", c0)
    np.save(out / "fc3_ref_compact.npy", f0)
    # cluster orbits
    sc = ph3.supercell
    d = pair_distances(np.array(sc.cell), np.array(sc.scaled_positions))
    npa, nsat = f0.shape[0], f0.shape[1]
    dmax = np.zeros((npa, nsat, nsat))
    for a, pa in enumerate(p2s):
        dmax[a] = np.maximum(np.maximum(d[pa][:, None], d[pa][None, :]), d)
    dm_el = np.repeat(dmax.ravel(), 27)
    trip_el = np.repeat(np.arange(npa * nsat * nsat), 27)
    csc = ccm.tocsc()
    ncol = csc.shape[1]
    parent = np.arange(ncol)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner, spread, dcol = {}, np.zeros(ncol), np.zeros(ncol)
    for k in range(ncol):
        rows = csc.indices[csc.indptr[k]:csc.indptr[k + 1]]
        dv = dm_el[rows]
        spread[k], dcol[k] = dv.max() - dv.min(), dv.mean()
        for t in np.unique(trip_el[rows]):
            if t in owner:
                a_, b_ = find(k), find(owner[t])
                if a_ != b_:
                    parent[a_] = b_
            else:
                owner[t] = k
    roots = np.array([find(k) for k in range(ncol)])
    _, blk = np.unique(roots, return_inverse=True)
    nblk = int(blk.max()) + 1
    eta_ref = B @ c0
    refnorm = np.sqrt(np.bincount(blk, weights=eta_ref ** 2, minlength=nblk))
    bdmax = np.array([dcol[blk == b].mean() for b in range(nblk)])
    np.savez(out / "cluster_orbits.npz", col_block=blk, block_dmax=bdmax, block_refnorm=refnorm, col_dmax=dcol)
    meta = {"inputs": inputs, "software": common.software_versions(), "ccm_shape": list(ccm.shape),
            "basis_shape": list(B.shape), "proj_resid_rel": resid, "n_cluster_orbits": nblk,
            "max_shell_spread_within_column": float(spread.max()),
            "refnorm_minmax": [float(refnorm.min()), float(refnorm.max())]}
    common.write_json(out / "basis_meta.json", meta)
    print(json.dumps({k: v for k, v in meta.items() if k not in ("inputs", "software")}, indent=1))


if __name__ == "__main__":
    main()
