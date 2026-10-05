"""Attack 2b: mesh dependence. 7x7x5 event geometry (exported with the campaign's aln_events.py into
review/cx_runs/events_m775_s0.1.npz): reference tau, tau_CS, tau_RTA, kappa, lambda_min(C_ref|H)
(=> rigorous factor-F ceiling F/lambda_min), share of numerically-forbidden events, and comparison with 5x5x3.
Output: review/cx_attack2_mesh775.json
"""
import json, sys, time
import numpy as np
import scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI, basis_H
t0 = time.time()
out = {}
for label, path in (("m553", None), ("m775", CAMP / "review" / "cx_runs" / "events_m775_s0.1.npz")):
    M = AlN() if path is None else AlN(path)
    o = M.orbits()
    gref = M.g_ref()
    forb = gref < 1e-12 * gref.max()
    r = M.diag(gref)
    C = M.C(gref)
    Q = basis_H(M.ehat)
    CH = Q.T @ C @ Q
    w, U = np.linalg.eigh(0.5 * (CH + CH.T))
    rec = {"n": M.n, "m": M.m, "n_ev_orb": int(o["n_ev_orb"]), "n_mode_orb": int(o["n_mode_orb"]),
           "frac_forbidden": float(forb.mean()), "lam_min_H": float(w[0]),
           "F1_ceiling_ps": float(1 / w[0] / FOURPI), "min_freq_THz": float(M.eps.min())}
    kf = M.kappa_factor()
    for c in (0, 2):
        bh = U.T @ (Q.T @ M.b[:, c])
        K = float(np.sum(bh ** 2 / w)); N = float(np.sum(bh ** 2 / w ** 2))
        b = M.b[:, c]
        KR = float(np.sum(b * b / r)); NR = float(np.sum(b * b / r ** 2))
        rec[c] = {"kappa_WmK": kf * K, "tau_ref_ps": N / K / FOURPI, "tau_CS_ps": K / float(b @ b) / FOURPI,
                  "tau_RTA_ps": NR / KR / FOURPI, "K_over_KRTA": K / KR,
                  "F_needed_for_ratio_100_over_tauCS": 100 * (K / float(b @ b)) * w[0]}
    # orbit-average sanity: raw vs averaged kappa
    out[label] = rec
    print(label, json.dumps(rec), time.time() - t0, flush=True)
(CAMP / "review" / "cx_attack2_mesh775.json").write_text(json.dumps(out, indent=1))
