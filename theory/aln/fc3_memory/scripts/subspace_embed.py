"""Embedding of the cutoff-restricted symmetric fc3 spaces into the full symmetric + ASR coordinates.

For cutoff r, symfc gives ccm_r, B_r (symmetric + ASR fc3 supported on triplets with all pair distances <= r).
These tensors are symmetric + ASR in the full space, so c = T_r c_r with T_r = 18 B^T ccm^T ccm_r B_r
(7752 x n_r). Checks: containment residual and isometry (T_r^T T_r = I). Saves T_r for restricted problems.
Run with the phono3py environment.
usage: python -B subspace_embed.py <resultdir> r1 r2 ...
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import aln_common as common  # noqa: E402
from p3_events import load_ph3  # noqa: E402
from variants_p3 import symfc_basis  # noqa: E402


def main():
    rdir = Path(sys.argv[1])
    radii = [float(x) for x in sys.argv[2:]]
    common.verify_inputs()
    ph3 = load_ph3(True, [5, 5, 3], 0.1)
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    out = {}
    for r in radii:
        ccm_r, B_r = symfc_basis(ph3, r)
        Xr = ccm_r @ B_r                                  # compact tensors (n_compact x n_r)
        T = 18.0 * (B.T @ (ccm.T @ Xr))                   # (7752 x n_r)
        back = ccm @ (B @ T)
        rec = {"n_r": int(B_r.shape[1]), "containment_resid_rel": float(np.linalg.norm(back - Xr) / np.linalg.norm(Xr)),
               "isometry_dev": float(np.abs(T.T @ T - np.eye(T.shape[1])).max())}
        np.save(rdir / f"embed_T_r{r:.2f}.npy", T)
        out[f"{r:.2f}"] = rec
        print(r, json.dumps(rec), flush=True)
    (rdir / "subspace_embed.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
