"""Referee check: spectral location of the tau change of each witness, per mesh.

For tensor W on mesh M: N_W - N_ref and K_W - K_ref decomposed by mode frequency (x = C^+ b, full operator),
expressed as fractions of N_ref and K_ref, plus the change carried by the modes with f < 1.15 f_min
(lowest |q| acoustic star of that mesh).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from referee_check import Mesh, ROOT  # noqa: E402

BINS = [0, 3, 5, 8, 12, 40]
CASES = {
    "m553": (ROOT / "m553" / "events_m553_s0.1.npz", ROOT / "m553" / "pp_r01_REF.npz",
             {"cal1_c0_max(E18)": ("x", ROOT / "m553" / "witness" / "pp_cal_r01_cal1_c0_max.npz"),
              "cal1_c2_max(E18)": ("z", ROOT / "m553" / "witness" / "pp_cal_r01_cal1_c2_max.npz"),
              "orb0.1_c2_max(E15)": ("z", ROOT / "m553" / "witness" / "pp_orb_r01_orb0.1_c2_max.npz"),
              "alamode": ("z", ROOT / "m553" / "alamode" / "pp_alm_r01_alamode.npz")}),
    "m775": (ROOT / "m775x" / "events_m775_s0.1.npz", ROOT / "m775x" / "pp_r01_REF.npz",
             {"cal1nz_J_c0_max(E20)": ("x", ROOT / "m775x" / "pp_r01_cal1nz_J_c0_max.npz"),
              "cal1nz_J_c2_max(E20)": ("z", ROOT / "m775x" / "pp_r01_cal1nz_J_c2_max.npz"),
              "orb0.1_J_c2_max(E20)": ("z", ROOT / "m775x" / "pp_r01_orb0.1_J_c2_max.npz"),
              "orb0.1_c2_max(E15)": ("z", ROOT / "m775x" / "pp_r01_orb0.1_c2_max.npz"),
              "cal1_c0_max(E18)": ("x", ROOT / "m775x" / "pp_r01_cal1_c0_max.npz")}),
    "m996": (ROOT / "m996" / "events_m996_s0.1.npz", ROOT / "m996" / "pp_r01_REF.npz",
             {"cal1nz_J_c0_max(E20)": ("x", ROOT / "m996" / "joint" / "pp_r01_cal1nz_J_c0_max.npz"),
              "cal1nz_J_c2_max(E20)": ("z", ROOT / "m996" / "joint" / "pp_r01_cal1nz_J_c2_max.npz"),
              "orb0.1_J_c0_max(E20)": ("x", ROOT / "m996" / "joint" / "pp_r01_orb0.1_J_c0_max.npz"),
              "orb0.1_J_c2_max(E20)": ("z", ROOT / "m996" / "joint" / "pp_r01_orb0.1_J_c2_max.npz"),
              "alamode": ("z", ROOT / "m996" / "pp_r01_alamode.npz")}),
}

out = {}
for mesh, (ev, ppref, wits) in CASES.items():
    M = Mesh(ev)
    f = M.eps
    g0, g0w = M.rates(ppref)
    m0 = M.metrics(g0, want_x=True)
    star = f < 1.15 * f.min()
    out[mesh] = {"f_min": float(f.min())}
    for nm, (ax, pth) in wits.items():
        g, _ = M.rates(pth, g0w_ref=g0w)
        m = M.metrics(g, want_x=True)
        c = 0 if ax == "x" else 2
        b = M.b[:, c]
        x0, x1 = m0[f"x_{ax}"], m[f"x_{ax}"]
        N0, N1 = x0 * x0, x1 * x1
        K0, K1 = b * x0, b * x1
        rows = []
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            s = (f >= lo) & (f < hi)
            rows.append({"f": [lo, hi], "n": int(s.sum()), "dN_over_Nref": float((N1[s] - N0[s]).sum() / N0.sum()),
                         "dK_over_Kref": float((K1[s] - K0[s]).sum() / K0.sum())})
        rec = {"axis": ax, "tau_ratio": m[f"tau_{ax}"] / m0[f"tau_{ax}"],
               "N_ratio": float(N1.sum() / N0.sum()), "K_ratio": float(K1.sum() / K0.sum()),
               "dN_lowest_star_over_Nref": float((N1[star] - N0[star]).sum() / N0.sum()),
               "bins": rows}
        out[mesh][nm] = rec
        print(mesh, nm, ax, "tau %.3f N %.3f K %.3f | dN(lowest star)/Nref %.3f | "
              % (rec["tau_ratio"], rec["N_ratio"], rec["K_ratio"], rec["dN_lowest_star_over_Nref"])
              + " ".join("[%g-%g] dN %.3f dK %.3f" % (r["f"][0], r["f"][1], r["dN_over_Nref"], r["dK_over_Kref"])
                         for r in rows), flush=True)
Path(Path(__file__).resolve().parent / "referee_anatomy.json").write_text(json.dumps(out, indent=1))
