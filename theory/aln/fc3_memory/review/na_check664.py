"""Cross-check of the lean solver (na_bigmesh.py) on 6x6x4 with the dense eigendecomposition route of na_common."""
import sys
sys.dont_write_bytecode = True
import numpy as np
import na_common as nc
op = nc.Operator(r"C:/Users/Koussay/AppData/Local/Temp/na_fc3/events_m664_s0.1.npz")
g = op.group_average(op.ev_g)
E = nc.eig_response(op, op.C(g))
print("tau", E["cols"]["x"]["tau_ps"], E["cols"]["z"]["tau_ps"], "kappa", E["cols"]["x"]["kappa"], E["cols"]["z"]["kappa"])
print("lam_min_H", E["lam_min_H"], "cond_H", E["cond_H"], "null overlap", E["null_overlap"], "lam_null", E["lam_null"])
print("N frac lowest", {a: E["cols"][a]["N_frac_lowest"] for a in "xz"})
