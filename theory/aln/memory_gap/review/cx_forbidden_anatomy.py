"""Attack 1b (anatomy): reliance of the witnesses on events whose reference vertex is numerically zero.

'Forbidden' = reference rate < 1e-12 * max (the reference histogram has an empty gap between 1e-15 and
1e-9 relative: these are selection-rule zeros / fc3 symmetry noise, symmetrize_fc=False).
For each witness: share of every mode's lifetime carried by forbidden events, the effect on tau of
zeroing them (data then violated), and a classification of forbidden events by the q-points involved.
Output: review/cx_forbidden_anatomy.json
"""
import json, sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, their_ev_orb, load_witness, CAMP

M = AlN(); o = M.orbits(); gref = M.g_ref(); r = M.diag(gref)
forb = gref < 1e-12 * gref.max()
out = {"n_forbidden": int(forb.sum()), "frac_forbidden": float(forb.mean()),
       "n_forbidden_orbits": int(len(np.unique(o["ev_orb"][forb]))),
       "orbit_consistency": bool(np.all([forb[o["ev_orb"] == k].all() or (~forb[o["ev_orb"] == k]).all()
                                        for k in np.unique(o["ev_orb"][forb])]))}
# where are the forbidden events? q-point of parent / daughters (high-symmetry?)
Nq, nb = M.Nq, M.nb
qidx = np.nonzero(M.valid)[0] // nb
addr = M.Z["addr"]
def is_hs(q):  # high-symmetry: fixed by more than identity+TR of the 24 operations
    rot = M.Z["rot_maps"]
    return int(np.sum(rot[:, q] == q))
stab = np.array([is_hs(q) for q in range(Nq)])
sP, sA, sB = stab[qidx[M.P]], stab[qidx[M.A]], stab[qidx[M.B]]
out["stabilizer_size_parent_forbidden_hist"] = {int(k): int(v) for k, v in zip(*np.unique(sP[forb], return_counts=True))}
out["stabilizer_size_parent_allowed_hist"] = {int(k): int(v) for k, v in zip(*np.unique(sP[~forb], return_counts=True))}
out["max_stabilizer_any_participant_forbidden_hist"] = {int(k): int(v) for k, v in zip(*np.unique(np.maximum(np.maximum(sP, sA), sB)[forb], return_counts=True))}
out["max_stabilizer_any_participant_allowed_hist"] = {int(k): int(v) for k, v in zip(*np.unique(np.maximum(np.maximum(sP, sA), sB)[~forb], return_counts=True))}
# raw (not orbit-averaged) phono3py rates for forbidden events
out["raw_rate_forbidden_max_rel"] = float(M.g_raw[forb].max() / M.g_raw.max())
out["raw_rate_allowed_min_rel"] = float(M.g_raw[~forb].min() / M.g_raw.max())

ev_orb_theirs, _ = their_ev_orb()
K0 = {c: M.moments(gref, (0, 2))[c]["K"] for c in (0, 2)}
for nm in ["full_c0_max_08_vertex", "full_c0_max_09_vertex", "full_c0_max_04_logn3.0", "full_c2_max_01_logn0.5",
           "full_c0_min_02_logn1.0", "full_c2_min_06_vertex"]:
    g = load_witness(nm, ev_orb_theirs)
    d_all = M.diag(g)
    d_forb = M.W @ (g * forb)
    share = d_forb / d_all
    g2 = g * (~forb)
    rec = {"lifetime_share_forbidden_max": float(share.max()), "lifetime_share_forbidden_median": float(np.median(share)),
           "n_modes_share_gt_1pct": int(np.sum(share > 0.01)), "n_modes_share_gt_10pct": int(np.sum(share > 0.1)),
           "total_rate_share_forbidden": float(g[forb].sum() / g.sum()),
           "n_forbidden_active_(g>1e-6*max)": int(np.sum(g[forb] > 1e-6 * g.max()))}
    try:
        mo2 = M.moments(g2, (0, 2))
        rec["zeroed: tau_x_ps"] = mo2[0]["tau_ps"]; rec["zeroed: tau_z_ps"] = mo2[2]["tau_ps"]
        rec["zeroed: K_x rel"] = mo2[0]["K"] / K0[0] - 1; rec["zeroed: K_z rel"] = mo2[2]["K"] / K0[2] - 1
        rec["zeroed: lifetime residual max"] = float(np.max(np.abs(M.diag(g2) / r - 1)))
        lam = M.moments_eig(g2, (0, 2))["lam_min_H"]
        rec["zeroed: lam_min_H"] = lam
    except Exception as exc:
        rec["zeroed"] = repr(exc)
    out[nm] = rec
    print(nm, json.dumps(rec))
(CAMP / "review" / "cx_forbidden_anatomy.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if not k.startswith("full_")}, indent=1))
