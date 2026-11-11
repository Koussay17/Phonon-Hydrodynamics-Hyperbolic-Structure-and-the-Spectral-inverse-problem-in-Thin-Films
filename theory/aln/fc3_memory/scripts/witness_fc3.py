"""Build explicit compact fc3 tensors for optimisation witnesses and describe where the change sits.

fc3 = fc3_0 + ccm @ B @ (|c0| T y)   (T = identity for the full space, embed_T_r*.npy for cutoff subspaces).
Anatomy: every compact element (a, J, K) gets the largest of its three minimum-image pair distances; the
squared norms of fc3_0 and of the change are accumulated per distance shell, and the relative change of the
nearest-neighbour cluster block is reported.
usage: python -B witness_fc3.py <resultdir> <outdir> <witness_y.npy> [subspace_r or none] [name]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp


def pair_distances(S):
    cell = S["sc_cell"]
    pos = S["sc_scaled"]
    n = len(pos)
    d = np.full((n, n), np.inf)
    rng = range(-2, 3)
    for i in rng:
        for j in rng:
            for k in rng:
                diff = pos[None, :, :] - pos[:, None, :] + np.array([i, j, k])[None, None, :]
                dist = np.linalg.norm(diff @ cell, axis=2)
                d = np.minimum(d, dist)
    return d


def anatomy(f0, df, S):
    d = pair_distances(S)
    p2s = S["p2s"]
    npa, nsat = f0.shape[0], f0.shape[1]
    dmax = np.zeros((npa, nsat, nsat))
    for a, pa in enumerate(p2s):
        dmax[a] = np.maximum(np.maximum(d[pa][:, None], d[pa][None, :]), d)
    e0 = (f0.reshape(npa, nsat, nsat, -1) ** 2).sum(-1)
    e1 = (df.reshape(npa, nsat, nsat, -1) ** 2).sum(-1)
    edges = [0, 0.01, 2.0, 3.2, 3.5, 4.0, 4.5, 5.0, 6.0, 20.0]
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (dmax >= lo) & (dmax < hi)
        rows.append({"shell_A": [lo, hi], "n_triplets": int(m.sum()), "ref_norm2_frac": float(e0[m].sum() / e0.sum()),
                     "delta_norm2_frac": float(e1[m].sum() / e1.sum()),
                     "rel_change_in_shell": float(np.sqrt(e1[m].sum() / max(e0[m].sum(), 1e-300)))})
    return rows


def main():
    rdir, outdir, ypath = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sub = sys.argv[4] if len(sys.argv) > 4 else "none"
    name = sys.argv[5] if len(sys.argv) > 5 else ypath.stem.replace("_y", "")
    outdir.mkdir(parents=True, exist_ok=True)
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    c0 = np.load(rdir / "basis_c0.npy")
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    y = np.load(ypath)
    if sub != "none":
        y = np.load(rdir / f"embed_T_r{float(sub):.2f}.npy") @ y
    c = np.linalg.norm(c0) * y
    df = (ccm @ (B @ c)).reshape(f0.shape)
    fc3 = f0 + df
    np.save(outdir / f"fc3_{name}.npy", fc3)
    S = np.load(rdir / "static_r01.npz")
    out = {"name": name, "eps_rel_frob": float(np.linalg.norm(df) / np.linalg.norm(f0)),
           "max_abs_change_over_max_abs_ref": float(np.abs(df).max() / np.abs(f0).max()),
           "anatomy": anatomy(f0, df, S)}
    (outdir / f"anatomy_{name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
