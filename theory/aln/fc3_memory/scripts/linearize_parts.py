"""Exact linearisation of the event amplitudes, split into row chunks for large meshes (7x7x5); same rows, same
engine and same formulas as linearize.py.

setup : rows (allowed orbit representatives x degenerate blocks), Phi0, metadata -> lin_meta.npz, linearization.json
part  : G rows [s, e) of chunk i of n -> lin_G_<i>.npy (complex, rows x 7752)
merge : for each scaling tag (orbit_scaling.py files in <scaledir>), Gs = G M^-1/2 (scaled coordinates u, see
        optimize_fc3.FC3ModelScaled) -> scaled_<tag>_Gr.npy / _Gi.npy in <rdir>; checks Phi0 = G c0.
usage: python -B linearize_parts.py setup <events.npz> <rdir>
       python -B linearize_parts.py part <events.npz> <rdir> <i> <n>
       python -B linearize_parts.py merge <rdir> <n> <scaledir> <tag> [<tag> ...]     (tag "orb" = untagged files)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(k, "2")
import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from amp_engine import AmpEngine  # noqa: E402
from fc3class import FC3Class, summary  # noqa: E402

TOL = 1e-6


def engine(rdir):
    P = np.load(rdir / "pp_r01_REF.npz")
    uc = float(P["unit_conversion"])
    return AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=uc), uc, P


def setup(ev_path, rdir):
    t0 = time.perf_counter()
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    eng, uc, P = engine(rdir)
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    c0 = np.load(rdir / "basis_c0.npy")
    diag = {"W_pp_vs_evg_relmax": float(np.max(np.abs(M.W * P["pp"] - M.geom.gphys) / M.geom.gphys))}
    g_ref_orb = np.bincount(M.orb, weights=M.geom.gphys, minlength=M.n_orb) / M.osize
    diag["reference_orbit_avg"] = summary(M.metrics(g_ref_orb[M.orb], want_lmin=False))
    allowed = g_ref_orb > 1e-9 * g_ref_orb.max()
    reps = np.array([np.nonzero(M.orb == o)[0][0] for o in range(M.n_orb)])
    fb = eng.f
    rows_tp, rows_b, row_orb = [], [], []
    blk = np.zeros(M.n_orb, dtype=np.int64)
    for o in np.nonzero(allowed)[0]:
        a = reps[o]
        t = M.tp[a]
        Ds = [np.nonzero(np.abs(fb[t[k]] - fb[t[k]][M.bands[a, k]]) < TOL)[0] for k in range(3)]
        for j0 in Ds[0]:
            for j1 in Ds[1]:
                for j2 in Ds[2]:
                    rows_tp.append(t); rows_b.append((j0, j1, j2)); row_orb.append(o); blk[o] += 1
    rows_tp, rows_b, row_orb = np.array(rows_tp), np.array(rows_b), np.array(row_orb)
    W_orb = M.W[reps]
    Phi0 = eng.phi_events(f0, rows_tp, rows_b)
    pp_blk = np.bincount(row_orb, weights=np.abs(Phi0) ** 2 * uc, minlength=M.n_orb)
    g_blk_orb = np.where(allowed, W_orb * pp_blk / np.maximum(blk, 1), g_ref_orb)
    diag.update({"n_orbits": int(M.n_orb), "n_allowed": int(allowed.sum()), "n_rows": int(len(row_orb)),
                 "reference_block_avg": summary(M.metrics(g_blk_orb[M.orb], want_lmin=False)),
                 "seconds": time.perf_counter() - t0})
    np.savez(rdir / "lin_meta.npz", Phi0=Phi0, row_orb=row_orb, blk=blk, W_orb=W_orb, allowed=allowed, reps=reps,
             rows_tp=rows_tp, rows_b=rows_b, g_ref_orb=g_ref_orb, g_blk_orb=g_blk_orb, uc=uc, c0=c0)
    (rdir / "linearization.json").write_text(json.dumps(diag, indent=1))
    print(json.dumps(diag, indent=1))


def part(rdir, i, n):
    t0 = time.perf_counter()
    eng, _, _ = engine(rdir)
    Z = np.load(rdir / "lin_meta.npz")
    rows_tp, rows_b = Z["rows_tp"], Z["rows_b"]
    nrows = len(rows_tp)
    s0, e0 = (nrows * i) // n, (nrows * (i + 1)) // n
    ccmT = sp.load_npz(rdir / "basis_ccm.npz").tocsr().T.tocsr()
    B = np.load(rdir / "basis_B.npy")
    G = np.zeros((e0 - s0, B.shape[1]), dtype=complex)
    bs = 24
    for s in range(s0, e0, bs):
        e = min(s + bs, e0)
        L = eng.event_functionals(rows_tp[s:e], rows_b[s:e])
        G[s - s0:e - s0] = ((ccmT @ L.real.T).T + 1j * (ccmT @ L.imag.T).T) @ B
        if ((s - s0) // bs) % 50 == 0:
            print(f"part {i}: rows {s - s0}/{e0 - s0} {time.perf_counter() - t0:.0f}s", flush=True)
    np.save(rdir / f"lin_G_{i}.npy", G)
    print(json.dumps({"part": i, "rows": [int(s0), int(e0)], "seconds": time.perf_counter() - t0}))


def merge(rdir, n, scaledir, tags):
    t0 = time.perf_counter()
    Z = np.load(rdir / "lin_meta.npz")
    G = np.concatenate([np.load(rdir / f"lin_G_{i}.npy") for i in range(n)])
    rel = float(np.max(np.abs(G @ Z["c0"] - Z["Phi0"]) / np.abs(Z["Phi0"])))
    info = {"n_rows": int(G.shape[0]), "Phi0_from_G_relmax": rel, "tags": {}}
    for tag in tags:
        pre = "scaled_" if tag == "orb" else f"scaled_{tag}_"
        Mmh = np.load(scaledir / f"{pre}Mmh.npy")
        np.save(rdir / f"{pre}Gr.npy", np.ascontiguousarray(G.real @ Mmh))
        np.save(rdir / f"{pre}Gi.npy", np.ascontiguousarray(G.imag @ Mmh))
        info["tags"][tag] = {"Mmh": str(scaledir / f"{pre}Mmh.npy")}
        print(tag, time.perf_counter() - t0, flush=True)
    info["seconds"] = time.perf_counter() - t0
    (rdir / "linearization_merge.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "setup":
        setup(sys.argv[2], Path(sys.argv[3]))
    elif mode == "part":
        part(Path(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    else:
        merge(Path(sys.argv[2]), int(sys.argv[3]), Path(sys.argv[4]), sys.argv[5:])
