"""Numerical audit, part 8: first-order optimality (KKT residual) of the E20 joint-mesh witnesses.

Every inner L-BFGS-B solve of E20 stopped at maxiter. This script measures how far the reported best feasible
iterate u* is from a first-order KKT point of   max log tau_obj(u)  s.t.  h_i(u) <= 0, with the constraint set of
joint_mesh.py (mode-orbit lifetimes +-10 % and kappa_xx, kappa_zz +-1 % on 5x5x3 and 7x7x5, per-cluster-orbit
bounds ||(W u)_o|| <= eps), using the campaign's own linearised models (joint_mesh.MeshModel; scaled matrices
memory-mapped). Active set: constraints within a band of their bound (|s| >= dr - a_r, |dK| >= dK - a_K,
||(W u)_o|| / eps >= 1 - a_o). KKT residual: rho = min_{lambda >= 0} ||grad f - J_act^T lambda|| / ||grad f||
(non-negative least squares). rho ~ 0: first-order stationary for that active set; rho ~ 1: the gradient has a large
component that no active constraint opposes, i.e. a first-order ascent direction exists.
Also: the first-order gain available along the projected direction d = (grad f - J^T lambda), normalised to the
step that keeps every active constraint to first order: reported as the directional derivative per unit |d|.
usage: python -B na_kkt.py <run name: cal1nz_x_max | cal1nz_z_max | orb01_x_max | orb01_z_max> [out.json]
       python -B na_kkt.py custom <u.npy> <tag> <eps> <col> <objective mesh> <constraint meshes, comma sep.> <out.json>
"""
import json
import os
import sys
import time
from pathlib import Path

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_k, "4")
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy.optimize import nnls  # noqa: E402

HERE = Path(__file__).resolve().parent
SCRIPTS = Path(r"C:/Users/Koussay/ResearchLab/orchestration/campaigns/20261005-fc3-memory/scripts")
sys.path.insert(0, str(SCRIPTS))
sys.dont_write_bytecode = True
_np_load = np.load


def _load(path, *a, **k):          # memory-map the large .npy matrices instead of reading them into memory
    if str(path).endswith(".npy") and "mmap_mode" not in k:
        k["mmap_mode"] = "r"
    return _np_load(path, *a, **k)


np.load = _load
import joint_mesh as J  # noqa: E402


def chain_batch(md, Dmat, Phi, chunk=2000):
    """Columns of d/du of sum_orb Dmat[orb, k] g_orb(u):  2 Re(...) as in FC3Model.chain, for k columns at once."""
    V = (md.w_row[:, None] * Dmat[md.row_orb]) * np.conj(Phi)[:, None]
    out = np.zeros((Dmat.shape[1], md.Gr.shape[1]))
    for s in range(0, md.Gr.shape[0], chunk):
        e = min(s + chunk, md.Gr.shape[0])
        out += V[s:e].real.T @ np.asarray(md.Gr[s:e]) - V[s:e].imag.T @ np.asarray(md.Gi[s:e])
    return 2.0 * out


def main():
    run = sys.argv[1]
    t0 = time.perf_counter()
    if run == "custom":
        upath, tag, eps, col, obj, cm, out_path = sys.argv[2:9]
        eps, col, out_path = float(eps), int(col), Path(out_path)
        H = {"constraint_meshes": cm.split(",")}
        u = np.asarray(_np_load(upath))
        run = Path(upath).stem
    else:
        out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"kkt_{run}.json"
        H = json.load(open(J.RES / "joint" / f"{run}.json"))
        tag, eps, col, obj = H["tag"], float(H["eps"]), int(H["col"]), H["objective_mesh"]
        u = np.asarray(_np_load(J.RES / "joint" / f"{run}_u.npy"))
    pre = "scaled_" if tag == "orb" else f"scaled_{tag}_"
    models = {m: J.MeshModel(*J.MESH[m], tag, J.RES / "m553") for m in sorted(set(H["constraint_meshes"]) | {obj})}
    Wmat = np.asarray(_np_load(J.RES / "m553" / f"{pre}W.npy"))
    col_block = _np_load(J.RES / "m553" / "cluster_orbits.npz")["col_block"]
    nb = int(col_block.max()) + 1
    dr, dK = 0.1, 0.01
    R = {"run": run, "tag": tag, "eps": eps, "col": col, "objective_mesh": obj}
    grads, kinds = [], []
    gf = None
    for m, md in models.items():
        ev = md.evaluate(u, grad=True)
        s = ev["r"] / md.r0 - 1.0
        kk = np.array([ev["K"][0] / md.K0[0] - 1.0, ev["K"][2] / md.K0[2] - 1.0])
        R[f"{m}_tau_ratio"] = [ev["tau"][0] / md.base["tau"][0], ev["tau"][2] / md.base["tau"][2]]
        R[f"{m}_max_abs_s"] = float(np.max(np.abs(s)))
        R[f"{m}_kappa_dev"] = kk.tolist()
        if m == obj:
            Df = ev["dN"][col] / ev["N"][col] - ev["dK"][col] / ev["K"][col]
            gf = md.chain(Df, ev["Phi"])
        # lifetime constraints: d s_mu / d g_orb = (Wr S^T)[mu, orb] / r0_mu
        orb = md.M.orb
        S = sp.csr_matrix((np.ones(len(orb)), (orb, np.arange(len(orb)))), shape=(md.n_orb, len(orb)))
        Wro = (md.Wr @ S.T).tocsc()                                   # (mode-orbit rows x event orbits)
        act = np.nonzero(np.abs(s) >= dr - 2e-3)[0]
        R[f"{m}_n_active_lifetimes"] = int(len(act))
        if len(act):
            Wm = np.zeros((len(s), len(act)))
            Wm[act, np.arange(len(act))] = np.sign(s[act]) / md.r0[act]
            Dmat = np.asarray(Wro.T @ Wm)
            grads.append(chain_batch(md, Dmat, ev["Phi"]))
            kinds += [f"{m}_life"] * len(act)
        for i, c in enumerate((0, 2)):
            if abs(kk[i]) >= dK - 3e-4:
                Dk = (np.sign(kk[i]) * ev["dK"][c] / md.K0[c])[:, None]
                grads.append(chain_batch(md, Dk, ev["Phi"]))
                kinds.append(f"{m}_kappa{'xyz'[c]}")
        print(m, "evaluated", round(time.perf_counter() - t0), flush=True)
    Wu = Wmat @ u
    rel2 = np.bincount(col_block, weights=Wu * Wu, minlength=nb)
    ratio = np.sqrt(rel2) / eps
    acto = np.nonzero(ratio >= 1 - 2e-3)[0]
    R["n_active_orbit_bounds"] = int(len(acto))
    if len(acto):
        Mo = np.zeros((len(Wu), len(acto)))
        for k, o in enumerate(acto):
            sel = col_block == o
            Mo[sel, k] = 2 * Wu[sel] / eps ** 2
        grads.append((Wmat.T @ Mo).T)
        kinds += ["orbit"] * len(acto)
    Jt = np.vstack(grads).T                                            # (n_var x n_active)
    Jn = Jt / np.linalg.norm(Jt, axis=0)
    lam, res = nnls(Jn, gf, maxiter=50 * Jn.shape[1])
    lsq = np.linalg.lstsq(Jn, gf, rcond=None)[0]
    R["n_var"] = int(len(u))
    R["n_active_total"] = int(Jt.shape[1])
    R["kkt_residual_nnls_rel"] = float(res / np.linalg.norm(gf))
    R["kkt_residual_unsigned_lsq_rel"] = float(np.linalg.norm(gf - Jn @ lsq) / np.linalg.norm(gf))
    kinds = np.array(kinds)
    R["multiplier_mass_by_kind"] = {k: float(lam[kinds == k].sum() / max(lam.sum(), 1e-300)) for k in np.unique(kinds)}
    d = gf - Jn @ lam
    R["projected_ascent_rate"] = float(gf @ d / np.linalg.norm(d)) if np.linalg.norm(d) > 0 else 0.0
    R["grad_f_norm"] = float(np.linalg.norm(gf))
    # sensitivity of rho to the width of the active band (wider band = more constraints counted as active)
    R["t_total"] = time.perf_counter() - t0
    out_path.write_text(json.dumps(R, indent=1))
    print(json.dumps(R, indent=1), flush=True)


if __name__ == "__main__":
    main()
