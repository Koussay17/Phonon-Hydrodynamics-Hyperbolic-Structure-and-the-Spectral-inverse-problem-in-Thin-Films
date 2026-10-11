"""Numerical audit, part 10: 5x5x3 per-cluster witnesses never evaluated on 9x9x6 (E15 orb 0.1 max x/z, review
continuations c1_x/z_cont02). phono3py 4.5.0 Rust amplitudes on the 9x9x6 event list (p3_events.py, outputs in the
temporary folder), operator of na_common.py, orbit and gauge-invariant block conventions.
usage: python -B na_untested.py"""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import na_common as nc
from na_transfer import degeneracy_groups
from na_blockavg import block_rates, response
T = Path(r"C:/Users/Koussay/AppData/Local/Temp/na_fc3/m996w")
op = nc.Operator(nc.MESHES["m996"][0])
grp = degeneracy_groups(op)
sizes = np.bincount(grp)
da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
P0c = np.load(nc.MESHES["m996"][1] + "/pp_r01_REF.npz")
P0 = np.load(T / "pp_r01_REF.npz")
R = {"REF_pp_recomputed_vs_campaign_maxrel": float(np.max(np.abs(P0["pp"] - P0c["pp"]) / np.maximum(P0c["pp"], 1e-300)))}
conv = float(P0["conv"])
g0o = op.group_average(op.rates_from_pp(P0["pp"], P0["g0w"], conv))
g0b = block_rates(op, grp, sizes, P0["pp"], P0["g0w"], conv)[0]
t0o, k0o = response(op, g0o)
t0b, k0b = response(op, g0b)
r0 = op.diagC(g0o)
for nm in ("orb0.1_c0_max", "orb0.1_c2_max", "c1_x_cont02", "c1_z_cont02"):
    P = np.load(T / f"pp_r01_{nm}.npz")
    go = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], conv))
    gb = block_rates(op, grp, sizes, P["pp"], P["g0w"], conv)[0]
    to, ko = response(op, go)
    tb, kb = response(op, gb)
    lr = da(op.diagC(go)) / da(r0) - 1
    R[nm] = {"tau_ratio_orbit": [to[i] / t0o[i] for i in range(2)], "tau_ratio_block": [tb[i] / t0b[i] for i in range(2)],
             "kappa_ratio_orbit": [ko[i] / k0o[i] for i in range(2)], "kappa_ratio_block": [kb[i] / k0b[i] for i in range(2)],
             "lifetime_degavg_minmax": [float(lr.min()), float(lr.max())]}
    print(nm, json.dumps({k: [round(v, 4) for v in R[nm][k]] for k in R[nm]}), flush=True)
print("REF recompute", R["REF_pp_recomputed_vs_campaign_maxrel"])
(HERE / "na_out" / "untested_m996.json").write_text(json.dumps(R, indent=1))
