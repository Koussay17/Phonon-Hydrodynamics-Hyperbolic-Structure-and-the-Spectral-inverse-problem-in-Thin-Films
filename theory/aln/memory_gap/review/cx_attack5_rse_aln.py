"""Attack 5b: does phono3py's RTA linewidth Gamma differ from r = diag C on modes with repeated-daughter
events (self-coupling counted once vs twice)? And does K >= K_RTA/3 hold with the phono3py linewidth?
Output: review/cx_attack5_rse_aln.json
"""
import json, sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP
M = AlN()
graw = M.g_raw
r = M.diag(graw)
gam = M.Z["gam_raw"].ravel()[M.valid]
# repeated-daughter contribution to diag at the daughter
Dinv2 = 1.0 / M.D ** 2
rep = M.rep
contrib_rep = np.bincount(M.A[rep], weights=graw[rep] * (M.sA[rep] ** 2) * Dinv2[M.A[rep]], minlength=M.n)
share = contrib_rep / r
rel = (r - gam) / gam
sel = share > 1e-3
out = {"n_modes_with_repeated_share_gt_1e-3": int(sel.sum()),
       "median_rel_dev_all": float(np.median(np.abs(rel))),
       "rel_dev_vs_half_repeated_share_selected": [[float(share[i]), float(rel[i])] for i in np.nonzero(sel)[0][:20]],
       "corr_rel_dev_vs_share_selected": float(np.corrcoef(share[sel], rel[sel])[0, 1]) if sel.sum() > 2 else None}
# fit rel ~ alpha * share on selected modes
if sel.sum() > 2:
    alpha = float(np.sum(share[sel] * rel[sel]) / np.sum(share[sel] ** 2))
    out["fit_rel_dev_equals_alpha_times_share_alpha"] = alpha
# L3 with phono3py Gamma
gref = M.g_ref()
mo = M.moments(gref, (0, 2))
for c in (0, 2):
    b = M.b[:, c]
    out[f"K_over_KRTA_gamma_{c}"] = float(mo[c]["K"] / np.sum(b * b / gam))
    out[f"K_over_KRTA_diag_{c}"] = float(mo[c]["K"] / np.sum(b * b / M.diag(gref)))
print(json.dumps(out, indent=1))
(CAMP / "review" / "cx_attack5_rse_aln.json").write_text(json.dumps(out, indent=1))
