"""Review R9: rebuild an fc3 tensor from optimiser coordinates of the review runs (scaled u -> c -> eta -> fc3) and
check the per-orbit bound DIRECTLY on the rebuilt tensor (projection 18 ccm^T df, independent of W / M^-1/2), the
symmetry residual and the acoustic sum rule. Writes review/work/fc3_<name>.npy (campaign files untouched).
usage: python -B r9_build.py <rdir> <sdir> <tag> <u.npy> <name> <outdir>
"""
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

rdir, sdir, tag, upath, name, outdir = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], Path(sys.argv[6])
pre = f"scaled_{tag}_" if tag else "scaled_"
Mmh = np.load(sdir / f"{pre}Mmh.npy", mmap_mode="r")
u = np.load(upath)
dc = Mmh @ u                           # c - c0 (scaled coordinates: Phi = Phi0 + G M^-1/2 u)
B = np.load(rdir / "basis_B.npy", mmap_mode="r")
ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
c0 = np.load(rdir / "basis_c0.npy")
f0 = np.load(rdir / "fc3_ref_compact.npy")
CO = np.load(rdir / "cluster_orbits.npz")
deta = B @ dc
df = (ccm @ deta).reshape(f0.shape)
fc3 = f0 + df
np.save(outdir / f"fc3_{name}.npy", fc3)
# independent re-projection of the tensor change
deta2 = 18.0 * (ccm.T @ df.ravel())
rel = np.sqrt(np.bincount(CO["col_block"], weights=deta2 ** 2)) / CO["block_refnorm"]
out = {"name": name, "frob_rel": float(np.linalg.norm(df) / np.linalg.norm(f0)),
       "orbit_rel_max_from_tensor": float(rel.max()),
       "reproj_resid": float(np.linalg.norm(ccm @ deta2 - df.ravel()) / np.linalg.norm(df)),
       "asr_rel": float(np.abs(fc3.sum(axis=2)).max() / np.abs(fc3).max())}
if tag:
    epsv = np.load(sdir / f"{pre}epsvec.npy") if (sdir / f"{pre}epsvec.npy").exists() else None
    if epsv is not None:
        out["orbit_rel_over_eps_max_from_tensor"] = float(np.max(rel / epsv))
(outdir / f"build_{name}.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
