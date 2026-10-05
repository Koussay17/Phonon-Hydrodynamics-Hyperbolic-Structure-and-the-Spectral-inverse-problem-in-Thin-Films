"""Verify multi-temperature / extra-data witnesses independently: rebuild rates at every temperature
from u (T-independent |Phi|^2 delta), check every constrained datum, report tau at 300 K and the
slow-mode anatomy. usage: python cx_verify_multiT.py id1 id2 ...   Output: review/cx_verify_multiT.json
"""
import json, sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI
RUNS = CAMP / "review" / "cx_runs"
OUTF = CAMP / "review" / "cx_verify_multiT.json"
M = AlN(); o = M.orbits(); gref = M.g_ref(); forb = gref < 1e-12 * gref.max(); w0 = gref / M.Bf
res = json.loads(OUTF.read_text()) if OUTF.exists() else {}
cache = {300.0: M}
for jid in sys.argv[1:]:
    job = json.loads((RUNS / f"{jid}.json").read_text())["job"]
    u = np.load(RUNS / f"{jid}_u.npy")
    temps = job.get("temps", [300.0])
    lifeT = job.get("lifetimes_T", temps)
    rec = {"temps": temps, "lifetimes_T": lifeT, "extras": job.get("extras")}
    for T in temps:
        if T not in cache:
            cache[T] = AlN(T=T)
        MT = cache[T]
        gT_ref = w0 * MT.Bf
        base = gT_ref * (~forb) if job.get("forbid_zero") else gT_ref
        g = base * np.exp(u[o["ev_orb"]])
        m_ref = MT.moments(gT_ref, (0, 2)); m = MT.moments(g, (0, 2))
        rr = MT.diag(g) / MT.diag(gT_ref)
        rec[f"T{int(T)}"] = {"K_x_rel": m[0]["K"] / m_ref[0]["K"] - 1, "K_z_rel": m[2]["K"] / m_ref[2]["K"] - 1,
                             "lifetime_ratio_range": [float(rr.min()), float(rr.max())],
                             "tau_x_ps": m[0]["tau_ps"], "tau_z_ps": m[2]["tau_ps"],
                             "tau_ref_x_ps": m_ref[0]["tau_ps"], "tau_ref_z_ps": m_ref[2]["tau_ps"]}
        if T == 300.0:
            me = MT.moments_eig(g, (0, 2))
            c = job["col"]
            w, U, Q = me["w"], me["U"], me["Q"]
            bh = U.T @ (Q.T @ MT.b[:, c])
            Kk = bh ** 2 / w; Nk = bh ** 2 / w ** 2; relax = 1 / w / FOURPI
            rec["lam_min_H"] = me["lam_min_H"]; rec["cond_H"] = me["cond_H"]
            rec["tau_eig_ps"] = me[c]["tau_ps"]
            rec["K_share_relax_gt_1000ps"] = float(Kk[relax > 1000].sum() / Kk.sum())
            rec["slowest_relax_ps"] = float(relax.max())
            rec["K_share_slowest"] = float(Kk[np.argmax(relax)] / Kk.sum())
    res[jid] = rec
    print(jid, json.dumps(rec), flush=True)
OUTF.write_text(json.dumps(res, indent=1))
