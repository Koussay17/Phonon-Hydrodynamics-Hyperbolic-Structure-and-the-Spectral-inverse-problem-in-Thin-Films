"""Exact linearisation of the event amplitudes in the symmetric fc3 space (fc3-class rates as functions of c).

Phi is linear in fc3, so for symmetric perturbations fc3 = fc3_0 + ccm B c (symfc basis, see variants_p3.py)
    Phi_m(c) = Phi0_m + G_m c   (exact, no truncation),  G = Lrows @ ccm @ B.
Rows m: for every ALLOWED event orbit O (reference orbit rate > 1e-9 of the maximum) the representative event
and all band triples of its degenerate block D(j0) x D(j1) x D(j2) at the same representative triplet.
The fc3-class orbit rate is the block mean (gauge invariant, constant on orbits for symmetric fc3):
    g_O(c) = W_O * unit_conversion * mean_{m in block(O)} |Phi_m(c)|^2 .
Selection-rule orbits (reference rate <= 1e-9 of max; exactly zero for symmetric fc3) keep their reference
(numerically zero) rates.
Outputs (npz): G (rows x n_c, complex128), Phi0 (rows), row_orb, blk, W_orb, allowed orbit ids, rep events,
and diagnostics (json): block-averaged versus orbit-averaged reference, rank of the rate Jacobian.
usage: python -B linearize.py <events.npz> <resultdir> [deg_tol_THz]
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


def degenerate_set(freqs_q, j, tol):
    return np.nonzero(np.abs(freqs_q - freqs_q[j]) < tol)[0]


def main():
    ev_path, rdir = sys.argv[1], Path(sys.argv[2])
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    t0 = time.perf_counter()
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    P = np.load(rdir / "pp_r01_REF.npz")
    uc = float(P["unit_conversion"])
    eng = AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=uc)
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    diag = {}
    # frequency gap statistics (degeneracy tolerance)
    fb = eng.f
    gaps = np.diff(np.sort(fb, axis=1), axis=1).ravel()
    diag["freq_gaps_counts"] = {"<1e-9": int((gaps < 1e-9).sum()), "1e-9..1e-6": int(((gaps >= 1e-9) & (gaps < 1e-6)).sum()),
                                "1e-6..1e-3": int(((gaps >= 1e-6) & (gaps < 1e-3)).sum()),
                                "1e-3..1e-2": int(((gaps >= 1e-3) & (gaps < 1e-2)).sum()), "tol": tol}
    g_ref_orb = np.bincount(M.orb, weights=M.geom.gphys, minlength=M.n_orb) / M.osize
    allowed = g_ref_orb > 1e-9 * g_ref_orb.max()
    reps = np.array([np.nonzero(M.orb == o)[0][0] for o in range(M.n_orb)])
    tp, bands = M.tp, M.bands
    rows_tp, rows_b, row_orb = [], [], []
    blk = np.zeros(M.n_orb, dtype=np.int64)
    for o in np.nonzero(allowed)[0]:
        a = reps[o]
        t = tp[a]
        Ds = [degenerate_set(fb[t[k]], bands[a, k], tol) for k in range(3)]
        cnt = 0
        for j0 in Ds[0]:
            for j1 in Ds[1]:
                for j2 in Ds[2]:
                    rows_tp.append(t)
                    rows_b.append((j0, j1, j2))
                    row_orb.append(o)
                    cnt += 1
        blk[o] = cnt
    rows_tp = np.array(rows_tp)
    rows_b = np.array(rows_b)
    row_orb = np.array(row_orb)
    nrows = len(row_orb)
    diag.update({"n_orbits": int(M.n_orb), "n_allowed_orbits": int(allowed.sum()), "n_rows": int(nrows),
                 "block_size_hist": {int(k): int(v) for k, v in zip(*np.unique(blk[allowed], return_counts=True))}})
    W_orb = M.W[reps]
    # check W constant on orbits
    Wo_mean = np.bincount(M.orb, weights=M.W, minlength=M.n_orb) / M.osize
    diag["W_orbit_rel_spread_max"] = float(np.max(np.abs(M.W - Wo_mean[M.orb]) / Wo_mean[M.orb]))
    # reference amplitudes and block-averaged rates
    Phi0 = eng.phi_events(f0, rows_tp, rows_b)
    pp_rows = np.abs(Phi0) ** 2 * uc
    pp_blk = np.bincount(row_orb, weights=pp_rows, minlength=M.n_orb)
    g_blk_orb = np.where(allowed, W_orb * pp_blk / np.maximum(blk, 1), g_ref_orb)
    rel = np.abs(g_blk_orb[allowed] - g_ref_orb[allowed]) / g_ref_orb[allowed]
    diag["blockavg_vs_orbitavg_rate_rel"] = {"max": float(rel.max()), "median": float(np.median(rel)),
                                             "n_gt_1e-6": int((rel > 1e-6).sum())}
    m_orb = M.metrics(g_ref_orb[M.orb])
    m_blk = M.metrics(g_blk_orb[M.orb])
    diag["reference_orbit_avg"] = summary(m_orb)
    diag["reference_block_avg"] = summary(m_blk)
    print(json.dumps(diag, indent=1), flush=True)
    # linear functionals -> G
    t1 = time.perf_counter()
    nc = B.shape[1]
    Geta = np.zeros((nrows, ccm.shape[1]), dtype=complex)
    ccmT = ccm.T.tocsr()
    bs = 24
    for s in range(0, nrows, bs):
        L = eng.event_functionals(rows_tp[s:s + bs], rows_b[s:s + bs])        # (bs, n_compact)
        Geta[s:s + bs] = (ccmT @ L.real.T).T + 1j * (ccmT @ L.imag.T).T
        if s == 0:
            # functional check: L @ f0 == Phi0
            diag["functional_check_relmax"] = float(np.max(np.abs(L @ f0.ravel() - Phi0[:bs]) / np.abs(Phi0[:bs])))
        if (s // bs) % 20 == 0:
            print(f"rows {s}/{nrows}  {time.perf_counter() - t1:.0f}s", flush=True)
    G = Geta @ B                                                                # (nrows, nc)
    diag["seconds_G"] = time.perf_counter() - t1
    # check: Phi0 == Geta @ eta0 (fc3_0 is in the symmetric space)
    c0 = np.load(rdir / "basis_c0.npy")
    diag["Phi0_from_G_relmax"] = float(np.max(np.abs(G @ c0 - Phi0) / np.abs(Phi0)))
    # rate Jacobian at c = 0 (allowed orbits) and its rank
    Jrow = 2 * (np.conj(Phi0)[:, None] * G).real * uc                         # d|Phi_m|^2/dc (nrows, nc)
    J = np.zeros((M.n_orb, nc))
    np.add.at(J, row_orb, Jrow)
    J = (W_orb / np.maximum(blk, 1))[:, None] * J
    Ja = J[allowed] / g_ref_orb[allowed][:, None]                              # relative rate Jacobian
    sv = np.linalg.svd(Ja, compute_uv=False)
    diag["rel_rate_jacobian"] = {"shape": list(Ja.shape), "sv_max": float(sv[0]), "sv_min": float(sv[-1]),
                                 "rank_1e-10": int((sv > 1e-10 * sv[0]).sum()), "rank_1e-6": int((sv > 1e-6 * sv[0]).sum()),
                                 "sv_quantiles": np.quantile(sv, [0, 0.01, 0.1, 0.5, 0.9, 1]).tolist()}
    GR = np.vstack([G.real, G.imag])
    svG = np.linalg.svd(GR, compute_uv=False)
    diag["amplitude_map_real_rank"] = {"shape": list(GR.shape), "rank_1e-10": int((svG > 1e-10 * svG[0]).sum()),
                                       "sv_max": float(svG[0]), "sv_min": float(svG[-1])}
    # per-row operator norm of c -> Phi_m (2 x nc real map) and fragility |Phi0| / (||G_m|| |c0|)
    a11 = (G.real ** 2).sum(1); a22 = (G.imag ** 2).sum(1); a12 = (G.real * G.imag).sum(1)
    opn = np.sqrt(0.5 * (a11 + a22) + np.sqrt(0.25 * (a11 - a22) ** 2 + a12 ** 2))
    frag = np.abs(Phi0) / (opn * np.linalg.norm(c0))
    diag["fragility_rho_star_quantiles(0,0.001,0.01,0.1,0.5,0.9,1)"] = np.quantile(frag, [0, 0.001, 0.01, 0.1, 0.5, 0.9, 1]).tolist()
    np.savez(rdir / "linearization.npz", G=G, Phi0=Phi0, row_orb=row_orb, blk=blk, W_orb=W_orb, allowed=allowed,
             reps=reps, rows_tp=rows_tp, rows_b=rows_b, g_ref_orb=g_ref_orb, g_blk_orb=g_blk_orb, uc=uc, c0=c0,
             opn=opn, Jrel_sv=sv)
    diag["seconds_total"] = time.perf_counter() - t0
    (rdir / "linearization.json").write_text(json.dumps(diag, indent=1))
    print(json.dumps({k: v for k, v in diag.items() if k not in ("reference_orbit_avg", "reference_block_avg")}, indent=1))


if __name__ == "__main__":
    main()
