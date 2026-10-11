"""Referee check: for the lowest-|q| modes, compare the change of x = C^+ b with the RTA expectation b/r.
If x changes by a large factor while r = diag C stays within the lifetime tolerance, the tau change is carried
by off-diagonal (in-scattering / drag) couplings rather than by the diagonal rates."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from referee_check import Mesh, ROOT
cases = [("m775", ROOT/"m775x"/"events_m775_s0.1.npz", ROOT/"m775x"/"pp_r01_REF.npz",
          [("cal1nz_J_c2_max(E20)", "z", ROOT/"m775x"/"pp_r01_cal1nz_J_c2_max.npz"),
           ("orb0.1_c2_max(E15)", "z", ROOT/"m775x"/"pp_r01_orb0.1_c2_max.npz")]),
         ("m553", ROOT/"m553"/"events_m553_s0.1.npz", ROOT/"m553"/"pp_r01_REF.npz",
          [("cal1_c0_max(E18)", "x", ROOT/"m553"/"witness"/"pp_cal_r01_cal1_c0_max.npz"),
           ("cal1_c2_max(E18)", "z", ROOT/"m553"/"witness"/"pp_cal_r01_cal1_c2_max.npz")])]
out = {}
for mesh, ev, ppref, wits in cases:
    M = Mesh(ev); f = M.eps
    g0, g0w = M.rates(ppref); m0 = M.metrics(g0, want_x=True)
    for nm, ax, pth in wits:
        c = 0 if ax == "x" else 2
        g, _ = M.rates(pth, g0w_ref=g0w); m = M.metrics(g, want_x=True)
        x0, x1, b = m0[f"x_{ax}"], m[f"x_{ax}"], M.b[:, c]
        N0 = (x0**2).sum()
        # modes ranked by their reference contribution to N
        top = np.argsort(-(x0**2))[:8]
        rows = []
        for i in top:
            rows.append({"f": round(float(f[i]), 3), "N_share_ref": round(float(x0[i]**2/N0), 4),
                         "x_ratio": round(float(x1[i]/x0[i]), 3),
                         "r_ratio": round(float(m["r"][i]/m0["r"][i]), 3),
                         "x_ref_over_RTA": round(float(x0[i]/(b[i]/m0["r"][i])), 3),
                         "x_wit_over_RTA": round(float(x1[i]/(b[i]/m["r"][i])), 3)})
        out[f"{mesh}:{nm}"] = rows
        print(mesh, nm, ax)
        for r in rows: print("   ", r)
Path("referee_drag.json").write_text(json.dumps(out, indent=1))
