"""Task 5: independent re-evaluation of fc3 witnesses.

Inputs: explicit witness fc3 (compact .npy, from witness_fc3.py) and its phono3py-route pp (p3_events.py --fc3).
Routes compared for kappa, tau_mem, lifetimes:
  (A) linear amplitude model + block averaging (optimisation model, from the witness y),
  (B) phono3py Rust interaction strengths on all events, orbit-averaged rates (memory-gap convention),
  (C) independent NumPy engine on all events with the explicit fc3, block-averaged per orbit;
and three linear-algebra routes for tau on (B): Cholesky, dense eigendecomposition on H, conjugate gradients.
usage: python -B verify_witness.py <events.npz> <resultdir> <witness_dir> <name> <pp.npz> [subspace or none]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI, summary  # noqa: E402
from amp_engine import AmpEngine  # noqa: E402


def tau_routes(M, g, col):
    C = M.ops.C(g)
    b = M.geom.b[:, col]
    e = M.geom.ehat
    # dense eigendecomposition on H
    w, V = np.linalg.eigh(C)
    ke = np.argmax(np.abs(V.T @ e))                      # eigenvector closest to the energy invariant
    keep = np.ones(len(w), bool); keep[ke] = False
    bc = V[:, keep].T @ b
    K_eig = float(np.sum(bc ** 2 / w[keep]))
    N_eig = float(np.sum(bc ** 2 / w[keep] ** 2))
    # conjugate gradients on H (projected operator)
    P = lambda v: v - e * (e @ v)
    A = spla.LinearOperator(C.shape, matvec=lambda v: P(C @ P(v)) + (np.trace(C) / len(C)) * e * (e @ v))
    x, info = spla.cg(A, b, rtol=1e-13, maxiter=20000)
    K_cg = float(b @ x); N_cg = float(x @ x)
    return {"tau_eig_ps": N_eig / K_eig / FOURPI, "tau_cg_ps": N_cg / K_cg / FOURPI, "cg_info": int(info),
            "K_eig": K_eig, "K_cg": K_cg, "lam_min_H": float(np.min(w[keep])), "cond_H": float(w[keep].max() / w[keep].min())}


def main():
    ev_path, rdir, wdir, name, pp_path = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4], sys.argv[5]
    sub = sys.argv[6] if len(sys.argv) > 6 else "none"
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    Z = np.load(rdir / "linearization.npz")
    out = {"name": name}
    # (B) phono3py route
    pp = np.load(pp_path)["pp"]
    gB = M.rates(pp)
    mB = M.metrics(gB)
    out["B_phono3py_orbitavg"] = summary(mB)
    for col in (0, 2):
        out[f"B_tau_routes_{'xyz'[col]}"] = tau_routes(M, gB, col)
    # (C) engine, explicit fc3, block-averaged
    fc3 = np.load(wdir / f"fc3_{name}.npy")
    eng = AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=float(Z["uc"]))
    phi = eng.phi_events(fc3, Z["rows_tp"], Z["rows_b"])
    row_orb, blk, W_orb, allowed = Z["row_orb"], Z["blk"], Z["W_orb"], Z["allowed"]
    g_orb = np.where(allowed, 0.0, Z["g_ref_orb"])
    g_orb += np.bincount(row_orb, weights=W_orb[row_orb] * float(Z["uc"]) * np.abs(phi) ** 2 / blk[row_orb],
                         minlength=M.n_orb)
    mC = M.metrics(g_orb[M.orb])
    out["C_engine_blockavg"] = summary(mC)
    # engine vs Rust on all events (raw pp)
    Pref = np.load(rdir / "pp_r01_REF.npz")
    phi_all = eng.phi_events(fc3, Pref["tp"], Pref["bands"])
    pe = np.abs(phi_all) ** 2 * float(Z["uc"])
    big = pp > 1e-9 * pp.max()
    out["engine_vs_rust_pp_relmax"] = float(np.max(np.abs(pe[big] - pp[big]) / pp[big]))
    # (A) linear model from y
    yfile = wdir / f"{name}_y.npy"
    if yfile.exists():
        from optimize_fc3 import FC3Model
        model = FC3Model(ev_path, rdir, subspace=None if sub == "none" else float(sub))
        ev = model.evaluate(np.load(yfile))
        out["A_model"] = {"tau_x": ev["tau"][0] / FOURPI, "tau_z": ev["tau"][2] / FOURPI,
                          "kappa_x": M.kf * ev["K"][0], "kappa_z": M.kf * ev["K"][2],
                          "max_lifetime_dev_vs_base": float(np.max(np.abs(ev["r"] / model.r0 - 1)))}
        # lifetimes of the phono3py route relative to the phono3py-route reference; per-mode values of
        # degenerate partners are gauge dependent in the orbit-average convention, so also report the
        # degeneracy-averaged deviation (modes at the same q with |df| < 1e-6 THz)
        r_ref = M.metrics(M.rates(Pref["pp"]))["r"]
        out["B_max_lifetime_dev_vs_orbitavg_ref"] = float(np.max(np.abs(mB["r"] / r_ref - 1)))
        Lb = M.geom.labels
        qidx = (np.arange(Lb["Nq"] * Lb["nb"]) // Lb["nb"])[Lb["valid"]]
        f = M.geom.eps
        grp = -np.ones(len(f), int)
        gcount = 0
        for q in np.unique(qidx):
            ii = np.nonzero(qidx == q)[0]
            ii = ii[np.argsort(f[ii])]
            st = 0
            for k in range(1, len(ii) + 1):
                if k == len(ii) or f[ii[k]] - f[ii[k - 1]] > 1e-6:
                    grp[ii[st:k]] = gcount
                    gcount += 1
                    st = k
        da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
        out["B_max_lifetime_dev_degavg"] = float(np.max(np.abs(da(mB["r"]) / da(r_ref) - 1)))
    (wdir / f"verify_{name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
