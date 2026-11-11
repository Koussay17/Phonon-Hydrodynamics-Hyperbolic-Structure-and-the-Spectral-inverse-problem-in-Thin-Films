"""Chunked, parallel version of linearize.py for larger meshes (7x7x5).

Same definitions as linearize.py (rows = allowed-orbit representatives and their degenerate blocks; exact linear
amplitude model Phi = Phi0 + G c in the 7752-dim symmetric + ASR fc3 space). G is written as real/imaginary
float64 parts per worker chunk: lin_Gr_<k>.npy, lin_Gi_<k>.npy; metadata in linearization.npz (no 'G' key,
'parts' gives the chunk count). The symfc basis files are mesh independent and read from <basisdir>.
usage: python -B linearize_big.py <events.npz> <resultdir> <basisdir> <nworkers> [deg_tol_THz]
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[k] = "1"
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

_W = {}


def _winit(rdir, bdir, uc):
    from amp_engine import AmpEngine
    _W["eng"] = AmpEngine(Path(rdir) / "static_r01.npz", r0avg=True, unit_conversion=uc)
    _W["ccmT"] = sp.load_npz(Path(bdir) / "basis_ccm.npz").tocsr().T.tocsr()
    _W["B"] = np.load(Path(bdir) / "basis_B.npy")
    _W["rdir"] = Path(rdir)


def _wchunk(args):
    k, tp, bands = args
    eng, ccmT, B = _W["eng"], _W["ccmT"], _W["B"]
    n = len(tp)
    Gr = np.zeros((n, B.shape[1]))
    Gi = np.zeros((n, B.shape[1]))
    bs = 24
    t0 = time.time()
    for s in range(0, n, bs):
        L = eng.event_functionals(tp[s:s + bs], bands[s:s + bs])
        Ge_r = (ccmT @ L.real.T).T
        Ge_i = (ccmT @ L.imag.T).T
        Gr[s:s + bs] = Ge_r @ B
        Gi[s:s + bs] = Ge_i @ B
    np.save(_W["rdir"] / f"lin_Gr_{k}.npy", Gr)
    np.save(_W["rdir"] / f"lin_Gi_{k}.npy", Gi)
    return k, n, time.time() - t0


def main():
    ev_path, rdir, bdir, nw = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), int(sys.argv[4])
    tol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-6
    t0 = time.perf_counter()
    from amp_engine import AmpEngine
    from fc3class import FC3Class, summary
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    P = np.load(rdir / "pp_r01_REF.npz")
    uc = float(P["unit_conversion"])
    eng = AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=uc)
    f0 = np.load(bdir / "fc3_ref_compact.npy")
    fb = eng.f
    diag = {}
    gaps = np.diff(np.sort(fb, axis=1), axis=1).ravel()
    diag["freq_gaps_counts"] = {"<1e-9": int((gaps < 1e-9).sum()), "1e-9..1e-6": int(((gaps >= 1e-9) & (gaps < 1e-6)).sum()),
                                "1e-6..1e-3": int(((gaps >= 1e-6) & (gaps < 1e-3)).sum()), "tol": tol}
    g_ref_orb = np.bincount(M.orb, weights=M.geom.gphys, minlength=M.n_orb) / M.osize
    allowed = g_ref_orb > 1e-9 * g_ref_orb.max()
    first = np.full(M.n_orb, -1)
    for i in range(len(M.orb) - 1, -1, -1):
        first[M.orb[i]] = i
    reps = first
    tp, bands = M.tp, M.bands
    rows_tp, rows_b, row_orb = [], [], []
    blk = np.zeros(M.n_orb, dtype=np.int64)
    for o in np.nonzero(allowed)[0]:
        a = reps[o]
        t = tp[a]
        Ds = [np.nonzero(np.abs(fb[t[k]] - fb[t[k]][bands[a, k]]) < tol)[0] for k in range(3)]
        cnt = 0
        for j0 in Ds[0]:
            for j1 in Ds[1]:
                for j2 in Ds[2]:
                    rows_tp.append(t); rows_b.append((j0, j1, j2)); row_orb.append(o); cnt += 1
        blk[o] = cnt
    rows_tp = np.array(rows_tp); rows_b = np.array(rows_b); row_orb = np.array(row_orb)
    nrows = len(row_orb)
    diag.update({"n_orbits": int(M.n_orb), "n_allowed_orbits": int(allowed.sum()), "n_rows": int(nrows),
                 "n_mode_orbits": int(M.tie["n_mode_orb"]),
                 "block_size_hist": {int(k): int(v) for k, v in zip(*np.unique(blk[allowed], return_counts=True))}})
    W_orb = M.W[reps]
    Phi0 = eng.phi_events(f0, rows_tp, rows_b)
    pp_blk = np.bincount(row_orb, weights=np.abs(Phi0) ** 2 * uc, minlength=M.n_orb)
    g_blk_orb = np.where(allowed, W_orb * pp_blk / np.maximum(blk, 1), g_ref_orb)
    rel = np.abs(g_blk_orb[allowed] - g_ref_orb[allowed]) / g_ref_orb[allowed]
    diag["blockavg_vs_orbitavg_rate_rel"] = {"max": float(rel.max()), "median": float(np.median(rel)),
                                             "n_gt_1e-6": int((rel > 1e-6).sum())}
    diag["reference_orbit_avg"] = summary(M.metrics(g_ref_orb[M.orb]))
    diag["reference_block_avg"] = summary(M.metrics(g_blk_orb[M.orb]))
    print(json.dumps(diag, indent=1), flush=True)
    c0 = np.load(bdir / "basis_c0.npy")
    # parallel G
    bounds = np.linspace(0, nrows, nw * 4 + 1).astype(int)
    jobs = [(k, rows_tp[bounds[k]:bounds[k + 1]], rows_b[bounds[k]:bounds[k + 1]]) for k in range(len(bounds) - 1)]
    with Pool(nw, initializer=_winit, initargs=(str(rdir), str(bdir), uc)) as pool:
        for k, n, dt in pool.imap_unordered(_wchunk, jobs):
            print(f"chunk {k}: {n} rows in {dt:.0f}s  (elapsed {time.perf_counter() - t0:.0f}s)", flush=True)
    # consistency: G c0 == Phi0 on the first chunk
    Gr0 = np.load(rdir / "lin_Gr_0.npy"); Gi0 = np.load(rdir / "lin_Gi_0.npy")
    n0 = Gr0.shape[0]
    diag["Phi0_from_G_relmax_chunk0"] = float(np.max(np.abs(Gr0 @ c0 + 1j * (Gi0 @ c0) - Phi0[:n0]) / np.abs(Phi0[:n0])))
    np.savez(rdir / "linearization.npz", Phi0=Phi0, row_orb=row_orb, blk=blk, W_orb=W_orb, allowed=allowed, reps=reps,
             rows_tp=rows_tp, rows_b=rows_b, g_ref_orb=g_ref_orb, g_blk_orb=g_blk_orb, uc=uc, c0=c0,
             parts=len(jobs))
    diag["seconds_total"] = time.perf_counter() - t0
    (rdir / "linearization.json").write_text(json.dumps(diag, indent=1))
    print(json.dumps({k: v for k, v in diag.items() if not k.startswith("reference")}, indent=1))


if __name__ == "__main__":
    main()
