"""Attack 3 (part 1): independent re-computation of the AlN reference and of the key witnesses.

Checks: event count, kernel, kappa and tau_mem of the reference; decode witnesses; orbit constancy
under independently computed orbits; lifetime and K residuals; tau by two routes (Cholesky,
eigendecomposition on an explicit basis of H); spectral anatomy of each witness (lambda_min on H,
weight of the slowest modes in K and N); fraction of suppressed / enhanced events.
Output: review/cx_verify_reference.json
"""
import json
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, their_ev_orb, load_witness, FOURPI, CAMP  # noqa: E402

t0 = time.time()
M = AlN()
o = M.orbits()
print("n", M.n, "m", M.m, "event orbits", o["n_ev_orb"], "mode orbits", o["n_mode_orb"], "umklapp", M.umk.mean())
gref = M.g_ref()
r = M.diag(gref)
out = {"n": M.n, "m": M.m, "n_ev_orb": int(o["n_ev_orb"]), "n_mode_orb": int(o["n_mode_orb"])}

# energy conservation of the event vectors
Ae = M.AT @ M.e
out["max_|a.e|/(|a||e|)"] = float(np.max(np.abs(Ae) / (np.sqrt(np.asarray(M.W.sum(axis=0)).ravel()) * np.linalg.norm(M.e))))
out["b.e/(|b||e|)"] = [float(M.b[:, c] @ M.e / np.linalg.norm(M.b[:, c]) / np.linalg.norm(M.e)) for c in range(3)]

kf = M.kappa_factor()
mo = M.moments(gref, cols=(0, 1, 2))
me = M.moments_eig(gref, cols=(0, 1, 2))
ref = {}
for c in (0, 1, 2):
    b = M.b[:, c]
    KR = float(np.sum(b * b / r))
    NR = float(np.sum(b * b / r ** 2))
    ref[c] = dict(K=mo[c]["K"], kappa_WmK=kf * mo[c]["K"], tau_ps_chol=mo[c]["tau_ps"], tau_ps_eig=me[c]["tau_ps"],
                  tau_CS_ps=mo[c]["K"] / float(b @ b) / FOURPI, K_RTA=KR, K_over_KRTA=mo[c]["K"] / KR,
                  tau_RTA_ps=NR / KR / FOURPI, b2=float(b @ b))
out["reference"] = ref
out["reference_lam_min_H"] = me["lam_min_H"]
out["reference_cond_H"] = me["cond_H"]
w, U, Q = me["w"], me["U"], me["Q"]
print("reference", json.dumps({c: {k: v for k, v in ref[c].items()} for c in ref}, indent=1))
print("lam_min_H(ref) =", w[:6])

# isotypic restriction: smallest eigenvalue of C_ref among eigenvectors with non-negligible b overlap
for c in (0, 2):
    bh = U.T @ (Q.T @ M.b[:, c])
    wt = bh ** 2 / np.sum(bh ** 2)
    sel = wt > 1e-20
    out[f"reference_lam_min_coupled_{c}"] = float(w[sel].min())
    out[f"reference_first_coupled_weights_{c}"] = [[float(w[i]), float(wt[i])] for i in np.nonzero(sel)[0][:5]]

# raw (not orbit averaged) rates
mo_raw = M.moments(M.g_raw, cols=(0, 2))
out["raw_rates"] = {c: dict(kappa=kf * mo_raw[c]["K"], tau_ps=mo_raw[c]["tau_ps"]) for c in (0, 2)}

# ------------------------------------------------------------------ witnesses
ev_orb_theirs, geom = their_ev_orb()
assert np.array_equal(geom.events[:, 0], M.P) and np.array_equal(geom.events[:, 1], M.A)
names = ["full_c0_max_08_vertex", "full_c0_max_09_vertex", "full_c0_max_04_logn3.0", "full_c2_max_01_logn0.5",
         "full_c0_min_02_logn1.0", "full_c2_min_06_vertex"]
K0 = {c: ref[c]["K"] for c in (0, 2)}
wit = {}
for nm in names:
    g = load_witness(nm, ev_orb_theirs)
    rec = {}
    # orbit constancy under MY orbits
    ga = M.orbit_average(g)
    rec["orbit_const_relmax"] = float(np.max(np.abs(ga - g) / np.maximum(np.abs(g), 1e-300)))
    rec["min_g"] = float(g.min())
    rd = M.diag(g)
    rec["diag_res_relmax"] = float(np.max(np.abs(rd - r) / r))
    mw = M.moments(g, cols=(0, 2))
    ew = M.moments_eig(g, cols=(0, 2))
    for c in (0, 2):
        rec[f"K{c}_rel"] = mw[c]["K"] / K0[c] - 1
        rec[f"tau{c}_ps_chol"] = mw[c]["tau_ps"]
        rec[f"tau{c}_ps_eig"] = ew[c]["tau_ps"]
    rec["lam_min_H"] = ew["lam_min_H"]
    rec["cond_H"] = ew["cond_H"]
    # spectral anatomy for the objective column
    c = 0 if "_c0_" in nm else 2
    w2, U2, Q2 = ew["w"], ew["U"], ew["Q"]
    bh = U2.T @ (Q2.T @ M.b[:, c])
    Kk = bh ** 2 / w2
    Nk = bh ** 2 / w2 ** 2
    order = np.argsort(w2)
    rec["slowest_modes"] = [dict(lam=float(w2[i]), relax_ps=float(1 / w2[i] / FOURPI), K_share=float(Kk[i] / Kk.sum()),
                                 N_share=float(Nk[i] / Nk.sum()), overlap2=float(bh[i] ** 2 / np.sum(bh ** 2)))
                            for i in order[:5]]
    # share of K carried by modes with relaxation time > 1 ns, > 100 ps
    relax_ps = 1 / w2 / FOURPI
    for thr in (100.0, 1000.0, 1e4, 1e5):
        rec[f"K_share_relax_gt_{int(thr)}ps"] = float(Kk[relax_ps > thr].sum() / Kk.sum())
        rec[f"N_share_relax_gt_{int(thr)}ps"] = float(Nk[relax_ps > thr].sum() / Nk.sum())
    # rate ratios
    ratio = g / gref
    rec["frac_events_below_1e-6_ref"] = float(np.mean(ratio < 1e-6))
    for F in (2, 10, 100, 1000):
        rec[f"frac_below_1/{F}"] = float(np.mean(ratio < 1.0 / F))
        rec[f"frac_above_{F}"] = float(np.mean(ratio > F))
    rec["ratio_max"] = float(ratio.max())
    rec["ratio_min"] = float(ratio.min())
    # rate-weighted share: fraction of total reference rate carried by suppressed events
    rec["ref_rate_share_suppressed_1e-6"] = float(gref[ratio < 1e-6].sum() / gref.sum())
    rec["umklapp_share_suppressed"] = float(M.umk[ratio < 1e-6].mean()) if np.any(ratio < 1e-6) else None
    wit[nm] = rec
    print(nm, json.dumps({k: v for k, v in rec.items() if k != "slowest_modes"}, indent=None))
    print("   slowest:", rec["slowest_modes"][:3])
out["witnesses"] = wit
out["seconds"] = time.time() - t0
(CAMP / "review" / "cx_verify_reference.json").write_text(json.dumps(out, indent=1, default=float))
print("done", time.time() - t0)
