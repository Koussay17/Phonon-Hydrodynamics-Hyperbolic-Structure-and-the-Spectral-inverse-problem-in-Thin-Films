"""Referee check (continued): diagonal versus off-diagonal mechanism, per witness and mesh, summarised over the
modes that carry the top 50 % of the reference N: N-weighted mean of log(x_ratio), split into the diagonal part
-log(r_ratio) (RTA expectation) and the remainder (off-diagonal / drag)."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from referee_check import Mesh, ROOT
cases = [
 ("m553", ROOT/"m553"/"events_m553_s0.1.npz", ROOT/"m553"/"pp_r01_REF.npz", [
   ("cal1_c0_max(E18)", "x", ROOT/"m553"/"witness"/"pp_cal_r01_cal1_c0_max.npz"),
   ("cal1_c2_max(E18)", "z", ROOT/"m553"/"witness"/"pp_cal_r01_cal1_c2_max.npz"),
   ("orb0.1_c0_max(E15)", "x", ROOT/"m553"/"witness"/"pp_orb_r01_orb0.1_c0_max.npz"),
   ("orb0.1_c2_max(E15)", "z", ROOT/"m553"/"witness"/"pp_orb_r01_orb0.1_c2_max.npz"),
   ("alamode", "z", ROOT/"m553"/"alamode"/"pp_alm_r01_alamode.npz")]),
 ("m775", ROOT/"m775x"/"events_m775_s0.1.npz", ROOT/"m775x"/"pp_r01_REF.npz", [
   ("cal1nz_J_c0_max(E20)", "x", ROOT/"m775x"/"pp_r01_cal1nz_J_c0_max.npz"),
   ("cal1nz_J_c2_max(E20)", "z", ROOT/"m775x"/"pp_r01_cal1nz_J_c2_max.npz"),
   ("orb0.1_J_c0_max(E20)", "x", ROOT/"m775x"/"pp_r01_orb0.1_J_c0_max.npz"),
   ("orb0.1_J_c2_max(E20)", "z", ROOT/"m775x"/"pp_r01_orb0.1_J_c2_max.npz"),
   ("orb0.1_c2_max(E15) transferred", "z", ROOT/"m775x"/"pp_r01_orb0.1_c2_max.npz"),
   ("cal1_c2_max(E18) transferred", "z", ROOT/"m775x"/"pp_r01_cal1_c2_max.npz")]),
 ("m996", ROOT/"m996"/"events_m996_s0.1.npz", ROOT/"m996"/"pp_r01_REF.npz", [
   ("orb0.1_J_c0_max(E20) validation", "x", ROOT/"m996"/"joint"/"pp_r01_orb0.1_J_c0_max.npz"),
   ("orb0.1_J_c2_max(E20) validation", "z", ROOT/"m996"/"joint"/"pp_r01_orb0.1_J_c2_max.npz"),
   ("cal1nz_J_c0_max(E20) validation", "x", ROOT/"m996"/"joint"/"pp_r01_cal1nz_J_c0_max.npz"),
   ("alamode", "z", ROOT/"m996"/"pp_r01_alamode.npz")])]
out = {}
for mesh, ev, ppref, wits in cases:
    M = Mesh(ev)
    g0, g0w = M.rates(ppref); m0 = M.metrics(g0, want_x=True)
    for nm, ax, pth in wits:
        c = 0 if ax == "x" else 2
        g, _ = M.rates(pth, g0w_ref=g0w); m = M.metrics(g, want_x=True)
        x0, x1 = m0[f"x_{ax}"], m[f"x_{ax}"]
        w = x0**2 / (x0**2).sum()
        order = np.argsort(-w); cum = np.cumsum(w[order]); sel = order[: int(np.searchsorted(cum, 0.5)) + 1]
        ok = sel[(x0[sel] != 0) & (x1[sel] / x0[sel] > 0)]
        ww = w[ok] / w[ok].sum()
        lx = np.log(x1[ok] / x0[ok]); ld = -np.log(m["r"][ok] / m0["r"][ok])
        rec = {"tau_ratio": round(m[f"tau_{ax}"] / m0[f"tau_{ax}"], 3), "n_modes_top50pctN": int(len(sel)),
               "f_range": [round(float(M.eps[sel].min()), 2), round(float(M.eps[sel].max()), 2)],
               "mean_log_x_ratio": round(float(ww @ lx), 3), "diagonal_part": round(float(ww @ ld), 3),
               "offdiag_part": round(float(ww @ (lx - ld)), 3),
               "min_r_ratio_in_set": round(float((m["r"][sel] / m0["r"][sel]).min()), 3)}
        out[f"{mesh}:{nm}"] = rec
        print(mesh, nm, ax, rec, flush=True)
Path("referee_drag2.json").write_text(json.dumps(out, indent=1))
