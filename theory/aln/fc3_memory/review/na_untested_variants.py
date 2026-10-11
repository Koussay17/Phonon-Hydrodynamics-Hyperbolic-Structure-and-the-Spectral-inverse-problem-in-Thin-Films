"""Numerical audit, part 10b: robustness of the 9x9x6 tau ratios of the untested 5x5x3 per-cluster witnesses (event cut
3 / 3.5 sigma, smearing 0.08 / 0.09 THz) and spectral anatomy (lambda_min, share of N in the lowest eigenpairs).
usage: python -B na_untested_variants.py"""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import na_common as nc
from na_transfer import analyse, variants
T = Path(r"C:/Users/Koussay/AppData/Local/Temp/na_fc3/m996w")
op = nc.Operator(nc.MESHES["m996"][0])
P0 = np.load(T / "pp_r01_REF.npz")
conv = float(P0["conv"])
g0 = op.group_average(op.rates_from_pp(P0["pp"], P0["g0w"], conv))
A0, V0 = analyse(op, g0), variants(op, g0)
R = {}
for nm in ("orb0.1_c2_max", "c1_z_cont02", "orb0.1_c0_max", "c1_x_cont02"):
    P = np.load(T / f"pp_r01_{nm}.npz")
    g = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], conv))
    A, V = analyse(op, g), variants(op, g)
    R[nm] = {"tau_ratio": [A["tau"][i] / A0["tau"][i] for i in range(2)], "lam_min_H_ratio": A["lam_min_H"] / A0["lam_min_H"],
             "cond_H": A["cond_H"], "N_frac_lowest5": A["N_frac_lowest5"], "N_frac_lowest20": A["N_frac_lowest20"],
             "N_top_pair": A["N_top_pair"], "variants_tau_ratio": {k: [V[k][0] / V0[k][0], V[k][1] / V0[k][1]] for k in V},
             "variants_kappa_ratio": {k: [V[k][2] / V0[k][2], V[k][3] / V0[k][3]] for k in V}}
    print(nm, json.dumps({k: R[nm][k] for k in ("tau_ratio", "lam_min_H_ratio", "N_frac_lowest5", "variants_tau_ratio")}), flush=True)
R["REF"] = {"lam_min_H": A0["lam_min_H"], "N_frac_lowest5": A0["N_frac_lowest5"], "N_frac_lowest20": A0["N_frac_lowest20"]}
(HERE / "na_out" / "untested_m996_variants.json").write_text(json.dumps(R, indent=1))
