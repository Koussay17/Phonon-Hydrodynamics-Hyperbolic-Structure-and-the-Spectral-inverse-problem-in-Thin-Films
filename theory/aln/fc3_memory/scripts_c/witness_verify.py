"""Build and independently verify optimisation witnesses (task 5), C: mirror.

build: fc3 = fc3_0 + ccm B (|c0| y) -> results/witness/fc3_<name>.npy, plus per-shell / per-cluster relative change.
check: route B = phono3py Rust interaction strengths for all events (pp file from p3_events.py --fc3), orbit-averaged
       rates; tau by Cholesky, dense eigendecomposition on H and conjugate gradients; lifetimes compared after
       averaging over degenerate partners; route A = linear model (optimisation model) at y.
usage: python -B witness_verify.py build <events> <rdir> <y.npy> <name>
       python -B witness_verify.py check <events> <rdir> <name> <pp.npz> <y.npy>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI, summary  # noqa: E402


def build(rdir, ypath, name):
    wdir = rdir / "witness"
    wdir.mkdir(exist_ok=True)
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    c0 = np.load(rdir / "basis_c0.npy")
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    CO = np.load(rdir / "cluster_orbits.npz")
    y = np.load(ypath)
    if name.startswith("orb") or name.startswith("ms_orb"):   # scaled coordinates u -> y = M^-1/2 u / |c0|
        y = (np.load(rdir / "scaled_Mmh.npy") @ y) / np.linalg.norm(c0)
        np.save(wdir / f"{name}_yunscaled.npy", y)
    elif name.startswith("cal"):            # eps-weighted coordinates of the calibrated class (tag cal<k>)
        y = (np.load(rdir / f"scaled_{name.split('_')[0]}_Mmh.npy") @ y) / np.linalg.norm(c0)
        np.save(wdir / f"{name}_yunscaled.npy", y)
    deta = np.linalg.norm(c0) * (B @ y)
    df = (ccm @ deta).reshape(f0.shape)
    np.save(wdir / f"fc3_{name}.npy", f0 + df)
    eta0 = B @ c0
    edges = [0, 0.01, 2.0, 3.2, 4.0, 5.0, 6.0, 20.0]
    shells = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (CO["col_dmax"] >= lo) & (CO["col_dmax"] < hi)
        shells.append({"shell_A": [lo, hi], "ref_norm2_frac": float((eta0[m] ** 2).sum() / (eta0 ** 2).sum()),
                       "rel_change": float(np.sqrt((deta[m] ** 2).sum() / (eta0[m] ** 2).sum()))})
    rel_orb = np.sqrt(np.bincount(CO["col_block"], weights=deta ** 2) / CO["block_refnorm"] ** 2)
    out = {"name": name, "eps_frob": float(np.linalg.norm(df) / np.linalg.norm(f0)),
           "max_orbit_rel_change": float(rel_orb.max()), "shells": shells}
    (wdir / f"anatomy_{name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


def tau_routes(M, g, col):
    C = M.ops.C(g)
    b, e = M.geom.b[:, col], M.geom.ehat
    w, V = np.linalg.eigh(C)
    keep = np.ones(len(w), bool)
    keep[np.argmax(np.abs(V.T @ e))] = False
    bc = V[:, keep].T @ b
    Ke, Ne = float(np.sum(bc ** 2 / w[keep])), float(np.sum(bc ** 2 / w[keep] ** 2))
    P = lambda v: v - e * (e @ v)
    A = spla.LinearOperator(C.shape, matvec=lambda v: P(C @ P(v)) + (np.trace(C) / len(C)) * e * (e @ v))
    x, info = spla.cg(A, b, rtol=1e-13, maxiter=50000)
    return {"tau_eig_ps": Ne / Ke / FOURPI, "tau_cg_ps": float(x @ x) / float(b @ x) / FOURPI, "cg_info": int(info),
            "cond_H": float(w[keep].max() / w[keep].min())}


def degavg_groups(M):
    L = M.geom.labels
    qidx = (np.arange(L["Nq"] * L["nb"]) // L["nb"])[L["valid"]]
    f = M.geom.eps
    grp = -np.ones(len(f), int)
    g = 0
    for q in np.unique(qidx):
        ii = np.nonzero(qidx == q)[0]
        ii = ii[np.argsort(f[ii])]
        st = 0
        for k in range(1, len(ii) + 1):
            if k == len(ii) or f[ii[k]] - f[ii[k - 1]] > 1e-6:
                grp[ii[st:k]] = g
                g += 1
                st = k
    return grp


def check(ev_path, rdir, name, pp_path, ypath):
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    pp_ref = np.load(rdir / "pp_r01_REF.npz")["pp"]
    gB = M.rates(np.load(pp_path)["pp"])
    mB = M.metrics(gB)
    out = {"name": name, "B_phono3py_orbitavg": summary(mB)}
    for c in (0, 2):
        out[f"B_tau_routes_{'xyz'[c]}"] = tau_routes(M, gB, c)
    grp = degavg_groups(M)
    da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
    r_ref = M.metrics(M.rates(pp_ref), cols=(), want_lmin=False)["r"]
    out["B_max_lifetime_dev_raw"] = float(np.max(np.abs(mB["r"] / r_ref - 1)))
    out["B_max_lifetime_dev_degavg"] = float(np.max(np.abs(da(mB["r"]) / da(r_ref) - 1)))
    m0 = M.metrics(M.rates(pp_ref))
    out["B_kappa_ratio_x"] = mB["kappa_x"] / m0["kappa_x"]
    out["B_kappa_ratio_z"] = mB["kappa_z"] / m0["kappa_z"]
    from optimize_fc3 import FC3Model
    model = FC3Model(ev_path, rdir)
    ev = model.evaluate(np.load(ypath))
    out["A_model"] = {"tau_x": ev["tau"][0] / FOURPI, "tau_z": ev["tau"][2] / FOURPI,
                      "kappa_x": M.kf * ev["K"][0], "kappa_z": M.kf * ev["K"][2]}
    (rdir / "witness" / f"verify_{name}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build(Path(sys.argv[3]), sys.argv[4], sys.argv[5])
    else:
        check(sys.argv[2], Path(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])
