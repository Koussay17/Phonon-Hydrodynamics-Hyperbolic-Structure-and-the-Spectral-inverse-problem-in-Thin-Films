"""Well-conditioned coordinates for per-cluster-orbit relative bounds.

D = per-column reference norm of the column's cluster orbit (cluster_orbits.npz), M = B^T D^-2 B,
u = M^(1/2) c  =>  w = D^-1 B c = W u with W = D^-1 B M^-1/2 isometric (W^T W = I): the per-orbit relative change
of fc3 is ||(W u)_o|| and the constraint 'every cluster orbit changes by at most eps relative to its reference'
reads ||(W u)_o|| <= eps. Amplitudes: Phi = Phi0 + G M^-1/2 u. Frobenius size: |M^-1/2 u| / |c0|.
Saves scaled_Gr.npy, scaled_Gi.npy, scaled_W.npy, scaled_Mmh.npy and a check json.
Weighted variant (calibrated class): --eps <npz with 'rel'> k floor tag replaces each orbit's reference norm by
refnorm_o * eps_o, eps_o = max(k rel_o, floor); then every per-orbit constraint reads ||(W u)_o|| <= 1. Files are
saved as scaled_<tag>_*.npy.
usage: python -B orbit_scaling.py <resultdir> [--eps <alamode_orbits.npz> <k> <floor> <tag>]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np


def main():
    rdir = Path(sys.argv[1])
    pre, epsb = "scaled_", None
    if len(sys.argv) > 2 and sys.argv[2] == "--eps":
        epsb = np.maximum(float(sys.argv[4]) * np.load(sys.argv[3])["rel"], float(sys.argv[5]))
        pre = f"scaled_{sys.argv[6]}_"
    t0 = time.perf_counter()
    B = np.load(rdir / "basis_B.npy")
    CO = np.load(rdir / "cluster_orbits.npz")
    refb = CO["block_refnorm"] if epsb is None else CO["block_refnorm"] * epsb
    Dcol = refb[CO["col_block"]]
    BD = B / Dcol[:, None]
    M = BD.T @ BD
    lam, Q = np.linalg.eigh(M)
    Mmh = (Q * lam ** -0.5) @ Q.T
    W = BD @ Mmh
    Z = np.load(rdir / "linearization.npz")
    Gu = Z["G"] @ Mmh
    np.save(rdir / f"{pre}Gr.npy", np.ascontiguousarray(Gu.real))
    np.save(rdir / f"{pre}Gi.npy", np.ascontiguousarray(Gu.imag))
    np.save(rdir / f"{pre}W.npy", W)
    np.save(rdir / f"{pre}Mmh.npy", Mmh)
    rng = np.random.default_rng(0)
    v = rng.standard_normal(W.shape[1])
    out = {"M_eig_minmax": [float(lam.min()), float(lam.max())], "cond_M": float(lam.max() / lam.min()),
           "isometry_check": float(abs(np.linalg.norm(W @ v) / np.linalg.norm(v) - 1)),
           "WtW_diag_dev_sample": float(np.abs((W[:, :50].T @ W[:, :50]) - np.eye(50)).max()),
           "seconds": time.perf_counter() - t0}
    (rdir / f"{pre.rstrip('_')}_orbit_scaling.json" if epsb is not None else rdir / "orbit_scaling.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
