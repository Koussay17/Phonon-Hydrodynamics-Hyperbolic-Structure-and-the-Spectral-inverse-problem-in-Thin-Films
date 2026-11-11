"""Exact linearisation of the event amplitudes in the symmetric fc3 space (C: mirror copy).

fc3 = fc3_0 + ccm B c  =>  Phi_m(c) = Phi0_m + G_m c exactly, G = L_rows ccm B.
Rows: representatives of ALLOWED event orbits (reference orbit rate > 1e-9 of max) and their degenerate blocks
D(j0) x D(j1) x D(j2) (tolerance 1e-6 THz). fc3-class orbit rate = W_O uc mean_block |Phi|^2 (gauge invariant).
Also re-checks the reference (orbit-averaged, memory-gap convention) on the regenerated data.
usage: python -B linearize.py <events.npz> <resultdir>
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from amp_engine import AmpEngine  # noqa: E402
from fc3class import FC3Class, summary  # noqa: E402


def main():
    ev_path, rdir = sys.argv[1], Path(sys.argv[2])
    tol = 1e-6
    t0 = time.perf_counter()
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    P = np.load(rdir / "pp_r01_REF.npz")
    uc = float(P["unit_conversion"])
    eng = AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=uc)
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    ccmT = sp.load_npz(rdir / "basis_ccm.npz").tocsr().T.tocsr()
    B = np.load(rdir / "basis_B.npy")
    c0 = np.load(rdir / "basis_c0.npy")
    diag = {"W_pp_vs_evg_relmax": float(np.max(np.abs(M.W * P["pp"] - M.geom.gphys) / M.geom.gphys))}
    # engine versus Rust on all events (reference fc3)
    phi_all = eng.phi_events(f0, P["tp"], P["bands"])
    big = P["pp"] > 1e-9 * P["pp"].max()
    diag["engine_vs_rust_relmax"] = float(np.max(np.abs(np.abs(phi_all[big]) ** 2 * uc - P["pp"][big]) / P["pp"][big]))
    g_ref_orb = np.bincount(M.orb, weights=M.geom.gphys, minlength=M.n_orb) / M.osize
    diag["reference_orbit_avg"] = summary(M.metrics(g_ref_orb[M.orb]))
    allowed = g_ref_orb > 1e-9 * g_ref_orb.max()
    reps = np.array([np.nonzero(M.orb == o)[0][0] for o in range(M.n_orb)])
    fb = eng.f
    rows_tp, rows_b, row_orb = [], [], []
    blk = np.zeros(M.n_orb, dtype=np.int64)
    for o in np.nonzero(allowed)[0]:
        a = reps[o]
        t = M.tp[a]
        Ds = [np.nonzero(np.abs(fb[t[k]] - fb[t[k]][M.bands[a, k]]) < tol)[0] for k in range(3)]
        for j0 in Ds[0]:
            for j1 in Ds[1]:
                for j2 in Ds[2]:
                    rows_tp.append(t); rows_b.append((j0, j1, j2)); row_orb.append(o); blk[o] += 1
    rows_tp, rows_b, row_orb = np.array(rows_tp), np.array(rows_b), np.array(row_orb)
    nrows = len(row_orb)
    W_orb = M.W[reps]
    Phi0 = eng.phi_events(f0, rows_tp, rows_b)
    pp_blk = np.bincount(row_orb, weights=np.abs(Phi0) ** 2 * uc, minlength=M.n_orb)
    g_blk_orb = np.where(allowed, W_orb * pp_blk / np.maximum(blk, 1), g_ref_orb)
    diag.update({"n_orbits": int(M.n_orb), "n_allowed": int(allowed.sum()), "n_rows": int(nrows),
                 "reference_block_avg": summary(M.metrics(g_blk_orb[M.orb]))})
    print(json.dumps(diag, indent=1), flush=True)
    G = np.zeros((nrows, B.shape[1]), dtype=complex)
    bs = 24
    for s in range(0, nrows, bs):
        L = eng.event_functionals(rows_tp[s:s + bs], rows_b[s:s + bs])
        G[s:s + bs] = ((ccmT @ L.real.T).T + 1j * (ccmT @ L.imag.T).T) @ B
        if (s // bs) % 25 == 0:
            print(f"rows {s}/{nrows} {time.perf_counter() - t0:.0f}s", flush=True)
    diag["Phi0_from_G_relmax"] = float(np.max(np.abs(G @ c0 - Phi0) / np.abs(Phi0)))
    np.savez(rdir / "linearization.npz", G=G, Phi0=Phi0, row_orb=row_orb, blk=blk, W_orb=W_orb, allowed=allowed,
             reps=reps, rows_tp=rows_tp, rows_b=rows_b, g_ref_orb=g_ref_orb, g_blk_orb=g_blk_orb, uc=uc, c0=c0)
    diag["seconds_total"] = time.perf_counter() - t0
    (rdir / "linearization.json").write_text(json.dumps(diag, indent=1))
    print(json.dumps({k: diag[k] for k in ("Phi0_from_G_relmax", "seconds_total")}))


if __name__ == "__main__":
    main()
