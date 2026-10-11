"""Review R4 (mesh variant, R4_MESH env, e.g. "7 7 5"): independent route for the witnesses through phono3py's own LBTE machinery (reducible collision
matrix on the full 5x5x3 mesh, ALL triplets with Gaussian weights, no 4-sigma event cut, phono3py's own assembly
and symmetrisation (Omega + Omega^T)/2). Shares only fc2, fc3 and the mesh with the campaign's event operator.

For each fc3: phono3py LBTE kappa (its own solver), and from the captured symmetrised collision matrix Omega and
phono3py's X vectors: K = X^T Omega^+ X, N = |Omega^+ X|^2, tau = N/K (pseudo-inverse dropping |eig| < 1e-8 as
phono3py does). Ratios witness / reference are unit free; the reference value fixes the unit (compared with the
campaign's 34.25 / 27.45 ps). Also the phono3py RTA linewidths gamma (self-energy convention) per mode,
degeneracy-averaged, relative to the reference.
phono3py env. usage: python -B r4_p3_lbte.py <outjson> [--sigma-cutoff X] name=path.npy ...  (name REF = dataset fc3)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("RAYON_NUM_THREADS", "2")
sys.dont_write_bytecode = True
CAMP = Path(r"C:\Users\Koussay\ResearchLab\orchestration\campaigns\20261005-fc3-memory\scripts")
sys.path.insert(0, str(CAMP))
import numpy as np  # noqa: E402

import aln_common as common  # noqa: E402
from p3_events import load_ph3  # noqa: E402
from phono3py.conductivity import collision_matrix_kernel as cmk  # noqa: E402

STASH = {}
_orig_prep = cmk.ReducibleCollisionMatrixKernel._prepare_collision_matrix
_orig_getX = cmk.CollisionMatrixKernel._get_X


def _prep(self):
    w = _orig_prep(self)
    STASH["Omega"] = np.array(self._collision_matrix, copy=True)
    STASH["freqs"] = np.array(self._frequencies, copy=True)
    return w


def _getX(self, i_temp, weights):
    X = _orig_getX(self, i_temp, weights)
    STASH["X"] = np.array(X, copy=True)
    return X


cmk.ReducibleCollisionMatrixKernel._prepare_collision_matrix = _prep
cmk.CollisionMatrixKernel._get_X = _getX


def degsets(freqs, tol=1e-6):
    grp = -np.ones(freqs.size, int)
    g = 0
    for q in range(freqs.shape[0]):
        f = freqs[q]
        order = np.argsort(f)
        st = 0
        for k in range(1, len(f) + 1):
            if k == len(f) or f[order[k]] - f[order[k - 1]] > tol:
                grp[q * len(f) + order[st:k]] = g
                g += 1
                st = k
    return grp


MESH = [int(v) for v in os.environ.get("R4_MESH", "5 5 3").split()]


def main():
    outp = Path(sys.argv[1])
    args = sys.argv[2:]
    sigma_cutoff = None
    if args and args[0] == "--sigma-cutoff":
        sigma_cutoff = float(args[1])
        args = args[2:]
    common.verify_inputs()
    ph3 = load_ph3(True, MESH, 0.1)
    if sigma_cutoff is not None:
        ph3.sigma_cutoff = sigma_cutoff
    fc3_ref = np.array(ph3.fc3)
    res = {"sigma_cutoff": sigma_cutoff, "runs": {}}
    for spec in args:
        name, path = spec.split("=", 1)
        fc3 = fc3_ref if path == "REF" else np.load(path)
        t0 = time.time()
        ph3.fc3 = np.ascontiguousarray(fc3, dtype="double")
        ph3.init_phph_interaction()
        STASH.clear()
        ph3.run_thermal_conductivity(is_LBTE=True, temperatures=[300.0], is_reducible_collision_matrix=True,
                                     log_level=0)
        tc = ph3.thermal_conductivity
        kappa = np.array(tc.kappa)[0, 0]           # (6,) xx yy zz yz xz xy
        Om = STASH["Omega"][0, 0]
        n = int(np.sqrt(Om.size))
        Om = Om.reshape(n, n)
        X = STASH["X"].reshape(n, 3)
        w, V = np.linalg.eigh(Om)
        keep = np.abs(w) > 1e-8
        rec = {"kappa_p3_xx_zz": [float(kappa[0]), float(kappa[2])], "n_pinv_dropped": int((~keep).sum()),
               "lam_min_kept": float(np.abs(w[keep]).min()), "asym": float(np.abs(Om - Om.T).max() / np.abs(Om).max())}
        for c, nm in ((0, "x"), (2, "z")):
            bc = V[:, keep].T @ X[:, c]
            K = float(np.sum(bc ** 2 / w[keep]))
            N = float(np.sum(bc ** 2 / w[keep] ** 2))
            rec[f"K_{nm}"], rec[f"N_{nm}"], rec[f"tauN/K_{nm}"] = K, N, N / K
            rec[f"tauCS_{nm}"] = K / float(X[:, c] @ X[:, c])
        gam = np.array(tc.gamma)[0, 0]             # (num_ir_gp, num_band), IR grid points of phono3py
        f_ir = STASH["freqs"][np.array(tc.grid_points)]
        rec["gamma_shape"] = list(gam.shape)
        np.save(outp.with_name(f"{outp.stem}_gamma_{name}.npy"), gam)
        np.save(outp.with_name(f"{outp.stem}_freqs_ir.npy"), f_ir)
        rec["seconds"] = time.time() - t0
        res["runs"][name] = rec
        print(name, json.dumps(rec), flush=True)
    ref = res["runs"].get("REF")
    if ref is not None:
        f_ir = np.load(outp.with_name(f"{outp.stem}_freqs_ir.npy"))
        g0 = np.load(outp.with_name(f"{outp.stem}_gamma_REF.npy"))
        grp = degsets(f_ir)
        da = lambda a: (np.bincount(grp, weights=a.ravel()) / np.bincount(grp))[grp]
        valid = (f_ir.ravel() > 1e-3) & (da(g0) > 1e-12)
        for name, rec in res["runs"].items():
            rec["ratio_kappa_xx_zz"] = [rec["kappa_p3_xx_zz"][0] / ref["kappa_p3_xx_zz"][0],
                                        rec["kappa_p3_xx_zz"][1] / ref["kappa_p3_xx_zz"][1]]
            rec["ratio_tau_x_z"] = [rec["tauN/K_x"] / ref["tauN/K_x"], rec["tauN/K_z"] / ref["tauN/K_z"]]
            rec["ratio_K_x_z"] = [rec["K_x"] / ref["K_x"], rec["K_z"] / ref["K_z"]]
            g = np.load(outp.with_name(f"{outp.stem}_gamma_{name}.npy"))
            rr = (da(g) / np.where(da(g0) > 0, da(g0), 1) - 1)[valid]
            rec["gamma_rel_change_degavg"] = {"min": float(rr.min()), "max": float(rr.max()),
                                              "rms": float(np.sqrt(np.mean(rr ** 2))),
                                              "frac_within_10pct": float(np.mean(np.abs(rr) <= 0.1))}
    outp.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
