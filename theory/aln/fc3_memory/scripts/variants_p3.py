"""Task 2 (b, c) on the phono3py route: fc3 variants evaluated with phono3py's Rust kernel on the fixed event list.

Symmetric fc3 space: symfc 1.7.3 FCBasisSetO3 on the 3x3x2 supercell (space group, index permutation,
acoustic sum rule). compact fc3 = ccm @ eta, eta = B @ c; ccm columns orthogonal with |col|^2 = 1/18, B
orthonormal, so the compact Frobenius norm is |c| / sqrt(18) for eta in span(B). The reference fc3 lies in
this space to 4.5e-15 (probe_symfc.json), with coordinates c0.

Variants (relative size eps = ||fc3 - fc3_0||_F / ||fc3_0||_F, compact form):
  cut_zero_r : phono3py cutoff_fc3_by_zero(r) (any pair distance > r -> 0), not re-symmetrised (ASR broken)
  cut_sym_r  : cut_zero_r projected orthogonally onto the symfc cutoff-r symmetric + ASR subspace
  iso_e_s    : c = c0 + eps |c0| u, u uniform on the unit sphere of the 7752-dim symmetric+ASR space
  mul_e_s    : delta = P_sym(xi * fc3_0), xi iid N(0,1) per compact element, rescaled to eps
  scale_e    : (1 + eps) fc3_0 (exact check: tau -> tau / (1 + eps)^2, kappa -> kappa / (1 + eps)^2)
Output: <outdir>/variants_<family>.npz (pp per variant, names, eps) and variants_meta.json;
        basis files basis_ccm.npz, basis_B.npy, basis_c0.npy (for linearize.py).
usage: python -B variants_p3.py --mesh 5 5 3 --events EV.npz --outdir DIR [--cutoffs ...] [--eps ...] [--seeds N]
                                [--save-fc3 NAME ...] [--no-basis-save]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("RAYON_NUM_THREADS", "4")
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import aln_common as common  # noqa: E402
from p3_events import load_ph3, event_rows  # noqa: E402


def symfc_basis(ph3, cutoff=None):
    from symfc.basis_sets import FCBasisSetO3
    from symfc.utils.utils import SymfcAtoms
    sc = ph3.supercell
    sa = SymfcAtoms(cell=sc.cell, scaled_positions=sc.scaled_positions, numbers=sc.numbers)
    bs = FCBasisSetO3(sa, cutoff=cutoff, use_mkl=False, log_level=0).run()
    ccm = sp.csr_matrix(bs.compact_compression_matrix)
    B = np.asarray(bs.basis_set)
    assert list(np.asarray(bs.p2s_map)) == list(ph3.primitive.p2s_map)
    return ccm, B


def project(ccm, B, f):
    """Orthogonal projection (compact Frobenius metric) of compact fc3 f onto span(ccm @ B)."""
    eta = 18.0 * (ccm.T @ f.ravel())
    c = B.T @ eta
    return (ccm @ (B @ c)).reshape(f.shape), c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", type=int, nargs=3, required=True)
    ap.add_argument("--events", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sigma", type=float, default=0.1)
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--cutoffs", type=float, nargs="*", default=[2.0, 3.2, 3.5, 4.0, 4.5, 5.0])
    ap.add_argument("--eps", type=float, nargs="*", default=[0.01, 0.05, 0.10, 0.20])
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--save-fc3", nargs="*", default=[])
    ap.add_argument("--no-basis-save", action="store_true")
    ap.add_argument("--families", nargs="*", default=["cut", "iso", "mul", "scale"])
    a = ap.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    t00 = time.perf_counter()
    inputs = common.verify_inputs()
    EV = np.load(a.events)
    ph3 = load_ph3(True, a.mesh, a.sigma)
    p2s = ph3.primitive.p2s_map
    f0 = np.array(ph3.fc3)[p2s].copy()                       # compact reference fc3
    n0 = np.linalg.norm(f0)
    ccm, B = symfc_basis(ph3, None)
    f0p, c0 = project(ccm, B, f0)
    meta = {"inputs": inputs, "software": common.software_versions(), "args": vars(a),
            "fc3_compact_norm": float(n0), "proj_resid_rel": float(np.linalg.norm(f0p - f0) / n0),
            "basis_shape": list(B.shape), "ccm_shape": list(ccm.shape), "c0_norm": float(np.linalg.norm(c0)),
            "c0_norm_over_sqrt18_fnorm": float(np.linalg.norm(c0) / np.sqrt(18) / n0), "variants": {}}
    if not a.no_basis_save:
        sp.save_npz(out / "basis_ccm.npz", ccm)
        np.save(out / "basis_B.npy", B)
        np.save(out / "basis_c0.npy", c0)
        np.save(out / "fc3_ref_compact.npy", f0)

    def run(name, fc3c, family, store):
        t0 = time.perf_counter()
        ph3.fc3 = np.ascontiguousarray(fc3c)
        res = event_rows(ph3, EV, a.sigma, a.T)
        store["names"].append(name)
        store["pp"].append(res["pp"])
        eps = float(np.linalg.norm(fc3c - f0) / n0)
        store["eps"].append(eps)
        rec = {"family": family, "eps_rel_frob": eps, "seconds": time.perf_counter() - t0}
        meta["variants"][name] = rec
        if name in a.save_fc3:
            np.save(out / f"fc3_{name}.npy", fc3c)
        print(name, json.dumps(rec), flush=True)

    fams = {}
    if "cut" in a.families:
        from phono3py.phonon3.fc3 import cutoff_fc3_by_zero
        st = fams.setdefault("cut", {"names": [], "pp": [], "eps": []})
        for r in a.cutoffs:
            fz = f0.copy()
            cutoff_fc3_by_zero(fz, ph3.supercell, r, p2s_map=p2s)
            run(f"cut_zero_{r:.2f}", fz, "cut_zero", st)
            ccm_r, B_r = symfc_basis(ph3, r)
            fs, _ = project(ccm_r, B_r, fz)
            meta["variants"].setdefault("_cut_basis_sizes", {})[f"{r:.2f}"] = list(B_r.shape)
            run(f"cut_sym_{r:.2f}", fs, "cut_sym", st)
    if "iso" in a.families or "mul" in a.families:
        for fam in ("iso", "mul"):
            if fam not in a.families:
                continue
            st = fams.setdefault(fam, {"names": [], "pp": [], "eps": []})
            for eps in a.eps:
                for s in range(a.seeds):
                    rng = np.random.default_rng(1000 * s + int(round(eps * 1000)) + (0 if fam == "iso" else 500000))
                    if fam == "iso":
                        u = rng.standard_normal(B.shape[1])
                        dc = u / np.linalg.norm(u) * eps * np.linalg.norm(c0)
                        d = (ccm @ (B @ dc)).reshape(f0.shape)
                    else:
                        xi = rng.standard_normal(f0.shape)
                        d, _ = project(ccm, B, xi * f0)
                        d *= eps * n0 / np.linalg.norm(d)
                    run(f"{fam}_{eps:.2f}_s{s}", f0 + d, fam, st)
    if "scale" in a.families:
        st = fams.setdefault("scale", {"names": [], "pp": [], "eps": []})
        run("scale_0.05", 1.05 * f0, "scale", st)
    for fam, st in fams.items():
        np.savez_compressed(out / f"variants_{fam}.npz", names=np.array(st["names"]), pp=np.array(st["pp"]),
                            eps=np.array(st["eps"]))
    meta["seconds_total"] = time.perf_counter() - t00
    common.write_json(out / f"variants_meta_{'_'.join(a.families)}.json", meta)


if __name__ == "__main__":
    main()
