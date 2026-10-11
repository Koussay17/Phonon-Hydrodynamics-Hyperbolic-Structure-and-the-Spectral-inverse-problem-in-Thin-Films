"""Numerical audit, part 6: gauge-invariant (degenerate-block) averaging of the event rates from the pp files.

The campaign's reported numbers use the orbit average of per-event rates, which is not gauge invariant for events
touching degenerate bands (09 item 2). The linearised model uses the block average: for the triplet of degenerate
groups (D0 at q0, D1 at q1, D2 at q2), every event of the block gets W_event x mean over the ORDERED band triples of
|Phi|^2. This script builds the block average directly from the pp file of any tensor on any mesh:
  key  = (group of p, unordered pair of daughter groups);
  sum  = sum over the events of the block of pp, counting an event with distinct daughters inside one group twice
         (ordered pairs) when both daughters belong to the same group;
  size = |D0| |D1| |D2| (ordered band triples); missing partner events (pp = 0, not exported) count as zeros;
  g    = W_event x sum / size.
Validation: reproduces the campaign's block-average references at 5x5x3 and 7x7x5 (linearization.json).
Then: block versus orbit convention for the reference and the witnesses at every mesh.
usage: python -B na_blockavg.py m553|m775|m996 [out.json]
"""
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import na_common as nc  # noqa: E402
from na_transfer import CASES, degeneracy_groups  # noqa: E402

CAMPAIGN_BLOCK_REF = {"m553": {"tau": [34.2530, 27.4410], "kappa": [214.2999, 201.9039]},
                      "m775": {"tau": [52.60872171535434, 50.62596472213647], "kappa": [243.02706782834366, 272.30194463964125]}}


def block_rates(op, grp, sizes, pp, g0w, conv):
    gp, ga, gb = grp[op.P], grp[op.A], grp[op.B]
    lo, hi = np.minimum(ga, gb), np.maximum(ga, gb)
    G = int(grp.max()) + 1
    key = (gp * G + lo) * G + hi
    uniq, inv = np.unique(key, return_inverse=True)
    same = (ga == gb) & (op.A != op.B)
    wsum = np.bincount(inv, weights=np.where(same, 2.0, 1.0) * pp, minlength=len(uniq))
    size = (sizes[gp] * sizes[ga] * sizes[gb]).astype(float)
    W = op.rates_from_pp(np.ones_like(pp), g0w, conv)
    return W * wsum[inv] / size, len(uniq)


def response(op, g):
    C = op.C(g)
    rr = nc.chol_response(op, C, np.trace(C) / op.n)
    return [rr["x"]["tau_ps"], rr["z"]["tau_ps"]], [rr["x"]["kappa"], rr["z"]["kappa"]]


def main():
    mk = sys.argv[1]
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"blockavg_{mk}.json"
    t0 = time.perf_counter()
    op = nc.Operator(nc.MESHES[mk][0])
    grp = degeneracy_groups(op)
    sizes = np.bincount(grp)
    cs = CASES[mk]
    names = ["REF"] + cs["names"]
    R = {"mesh": op.mesh, "n_degenerate_modes": int(np.sum(sizes[grp] > 1)), "tensors": {}}
    base = {}
    for nm in names:
        if nm == "REF":
            path = nc.MESHES[mk][1] + "/pp_r01_REF.npz"
        elif "/" in nm:
            path = cs["dir"] + "/" + nm.split("/")[0] + "/pp_r01_" + nm.split("/")[1] + ".npz"
        else:
            path = cs["dir"] + f"/pp_r01_{nm}.npz"
        P = np.load(path)
        conv = float(P["conv"])
        g_orb = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], conv))
        g_blk, nblk = block_rates(op, grp, sizes, P["pp"], P["g0w"], conv)
        sym = float(np.max(np.abs(op.group_average(g_blk) - g_blk)) / np.max(g_blk))
        t_o, k_o = response(op, g_orb)
        t_b, k_b = response(op, g_blk)
        rec = {"orbit": {"tau": t_o, "kappa": k_o}, "block": {"tau": t_b, "kappa": k_b}, "n_blocks": nblk,
               "block_rates_symmetry_defect": sym,
               "block_minus_orbit_rel": {"tau": [t_b[i] / t_o[i] - 1 for i in range(2)], "kappa": [k_b[i] / k_o[i] - 1 for i in range(2)]}}
        if nm == "REF":
            base = rec
            if mk in CAMPAIGN_BLOCK_REF:
                cb = CAMPAIGN_BLOCK_REF[mk]
                rec["block_vs_campaign_block_rel"] = {"tau": [t_b[i] / cb["tau"][i] - 1 for i in range(2)],
                                                      "kappa": [k_b[i] / cb["kappa"][i] - 1 for i in range(2)]}
        else:
            rec["ratio_orbit"] = [t_o[i] / base["orbit"]["tau"][i] for i in range(2)]
            rec["ratio_block"] = [t_b[i] / base["block"]["tau"][i] for i in range(2)]
            rec["ratio_block_minus_orbit"] = [rec["ratio_block"][i] - rec["ratio_orbit"][i] for i in range(2)]
            rec["kappa_ratio_orbit"] = [k_o[i] / base["orbit"]["kappa"][i] for i in range(2)]
            rec["kappa_ratio_block"] = [k_b[i] / base["block"]["kappa"][i] for i in range(2)]
        R["tensors"][nm.split("/")[-1]] = rec
        print(nm, json.dumps({k: rec[k] for k in rec if k in ("block_minus_orbit_rel", "block_vs_campaign_block_rel",
                                                              "ratio_block_minus_orbit", "block_rates_symmetry_defect")}),
              round(time.perf_counter() - t0), flush=True)
    out_path.write_text(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
