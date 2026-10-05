"""Verify review witnesses (cx_runs/<id>_u.npy) independently of the optimiser:
rebuild g from u (forbidden events exactly zero if the job says so), check lifetimes and K residuals
with the Cholesky and eigen routes, lambda_min on H, slow-mode K share, suppression statistics,
K(z) dispersion, accumulation, and the kappa / lifetimes implied at 100-800 K (T-independent |Phi|^2).
usage: python cx_verify_newwitness.py id1 id2 ...
Output: review/cx_verify_newwitness.json (appended / updated per id)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI  # noqa: E402

RUNS = CAMP / "review" / "cx_runs"
OUTF = CAMP / "review" / "cx_verify_newwitness.json"
M = AlN()
o = M.orbits()
gref = M.g_ref()
r = M.diag(gref)
forb = gref < 1e-12 * gref.max()
mo_ref = M.moments(gref, (0, 2))
K0 = {c: mo_ref[c]["K"] for c in (0, 2)}
w0 = gref / M.Bf
MT = {T: AlN(T=T) for T in (100.0, 200.0, 500.0, 800.0)}
refT = {T: MM.moments(w0 * MM.Bf, (0, 2)) for T, MM in MT.items()}
res = json.loads(OUTF.read_text()) if OUTF.exists() else {}
for jid in sys.argv[1:]:
    job = json.loads((RUNS / f"{jid}.json").read_text())["job"]
    u = np.load(RUNS / f"{jid}_u.npy")
    base = gref * (~forb) if job.get("forbid_zero") else gref
    g = base * np.exp(u[o["ev_orb"]])
    rec = {"forbid_zero": bool(job.get("forbid_zero", False)), "start": job.get("start")}
    rec["forbidden_rate_max"] = float(np.max(g[forb])) if job.get("forbid_zero") else None
    rec["diag_res_relmax"] = float(np.max(np.abs(M.diag(g) / r - 1)))
    mc = M.moments(g, (0, 2))
    me = M.moments_eig(g, (0, 2))
    for c in (0, 2):
        rec[f"K{c}_rel"] = mc[c]["K"] / K0[c] - 1
        rec[f"tau{c}_ps_chol"] = mc[c]["tau_ps"]
        rec[f"tau{c}_ps_eig"] = me[c]["tau_ps"]
    rec["lam_min_H"] = me["lam_min_H"]
    rec["cond_H"] = me["cond_H"]
    c = job["col"]
    w, U, Q = me["w"], me["U"], me["Q"]
    bh = U.T @ (Q.T @ M.b[:, c])
    Kk = bh ** 2 / w
    Nk = bh ** 2 / w ** 2
    relax = 1 / w / FOURPI
    rec["K_share_relax_gt_1000ps"] = float(Kk[relax > 1000].sum() / Kk.sum())
    rec["N_share_relax_gt_1000ps"] = float(Nk[relax > 1000].sum() / Nk.sum())
    rec["slowest_relax_ps"] = float(relax.max())
    ratio = g / gref
    al = ~forb
    rec["frac_allowed_events_below_1e-6_ref"] = float(np.mean(ratio[al] < 1e-6))
    rec["frac_allowed_events_below_1e-2_ref"] = float(np.mean(ratio[al] < 1e-2))
    rec["min_ratio_allowed"] = float(ratio[al].min())
    rec["F_needed_lower_bound_(tau*lam_ref)"] = float(mc[c]["tau_ps"] * FOURPI * 9.764279680497229e-4)
    # dispersion and accumulation
    C = M.C(g)
    b = M.b[:, c]
    for lab, zp in {"1/10ns": 1e-4, "1/ns": 1e-3, "1/100ps": 1e-2}.items():
        z = zp / FOURPI
        xz = np.linalg.solve(C + z * np.eye(M.n), b)
        xr = np.linalg.solve(M.C(gref) + z * np.eye(M.n), b)
        rec[f"K({lab})/K_ref({lab})"] = float(b @ xz) / float(b @ xr)
    x = mc[c]["x"]
    xr = mo_ref[c]["x"]
    for fc in (5.0, 10.0):
        mk = M.eps < fc
        rec[f"acc<{fc:g}THz (ref)"] = [float(b[mk] @ x[mk] / K0[c]), float(b[mk] @ xr[mk] / K0[c])]
    # other temperatures
    ww = g / M.Bf
    for T, MM in MT.items():
        gT = ww * MM.Bf
        m2 = MM.moments(gT, (0, 2))
        rT = MM.diag(gT) / MM.diag(w0 * MM.Bf)
        rec[f"T{int(T)}"] = {"kappa_x_ratio": m2[0]["K"] / refT[T][0]["K"], "kappa_z_ratio": m2[2]["K"] / refT[T][2]["K"],
                             "lifetime_ratio_range": [float(rT.min()), float(rT.max())]}
    res[jid] = rec
    print(jid, json.dumps(rec), flush=True)
OUTF.write_text(json.dumps(res, indent=1))
