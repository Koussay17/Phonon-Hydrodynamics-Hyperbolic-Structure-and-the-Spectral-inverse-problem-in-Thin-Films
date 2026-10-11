"""First-order (local) analysis of the fc3 class at the reference fc3.

y = c/|c0| (relative coordinates, Frobenius metric). Data map D(y) = (log r_mu (mode-orbit reps), log K_xx, log K_zz).
J_D = dD/dy at y = 0 (rows: 117 lifetimes + 2 kappa components at 5x5x3); singular spectrum = how many fc3
directions the static data see and how well. The first-order slope of log tau along data-preserving directions,
    s = | P_null(J_D) grad_y log tau |,
gives tau_max(rho) ~ tau_0 exp(s rho) for small rho with the data held fixed to first order (no slack).
Also for the short-range subspace (cutoff 3.5 A) and per mesh.
usage: python -B local_analysis.py <events.npz> <resultdir> [subspace ...]
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "1")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

from optimize_fc3 import FC3Model  # noqa: E402


def analyse(model):
    ev = model.evaluate(np.zeros(model.n_var), cols=(0, 2), grad=True)
    Phi = ev["Phi"]
    out = {"n_var": int(model.n_var)}
    rows = []
    for mu in range(len(model.r0)):
        w = np.zeros(len(model.r0)); w[mu] = 1.0 / model.r0[mu]
        rows.append(model.chain(model.dr_dorb(w), Phi))
    for c in (0, 2):
        rows.append(model.chain(ev["dK"][c], Phi) / ev["K"][c])
    J = np.array(rows)
    U, s, Vt = np.linalg.svd(J, full_matrices=False)
    rank = int((s > 1e-10 * s[0]).sum())
    out["J_data_shape"] = list(J.shape)
    out["J_data_rank"] = rank
    out["J_data_sv_quantiles(min,10%,50%,90%,max)"] = np.quantile(s, [0, 0.1, 0.5, 0.9, 1]).tolist()
    V = Vt[:rank].T
    for c in (0, 2):
        Dorb = ev["dN"][c] / ev["N"][c] - ev["dK"][c] / ev["K"][c]
        gt = model.chain(Dorb, Phi)
        gp = gt - V @ (V.T @ gt)
        out[f"grad_logtau_{'xyz'[c]}_norm"] = float(np.linalg.norm(gt))
        out[f"slope_s_{'xyz'[c]}_data_fixed"] = float(np.linalg.norm(gp))
        out[f"fraction_of_grad_in_data_nullspace_{'xyz'[c]}"] = float(np.linalg.norm(gp) / np.linalg.norm(gt))
    return out


def main():
    ev_path, rdir = sys.argv[1], Path(sys.argv[2])
    subs = [None] + [float(x) for x in sys.argv[3:]]
    res = {}
    for sub in subs:
        model = FC3Model(ev_path, rdir, subspace=sub)
        key = "full" if sub is None else f"cut{sub:.2f}"
        res[key] = analyse(model)
        print(key, json.dumps(res[key]), flush=True)
    (rdir / "local_analysis.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
