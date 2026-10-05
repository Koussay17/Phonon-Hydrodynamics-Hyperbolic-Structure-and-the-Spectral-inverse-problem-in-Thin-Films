"""Attack 1c (direct): compare event rates exported with symmetrize_fc=False (campaign) and True.
Selection-rule zeros polluted by fc3 asymmetry noise should drop further with symmetrised fc3,
while allowed rates change little. Output: review/cx_selection_symfc_compare.json
"""
import json
import numpy as np
base = r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap"
A = np.load(base + r"\results\aln\events_m553_s0.1.npz")
B = np.load(base + r"\review\cx_runs\events_m553_symfc.npz")
def keys(Z):
    return Z["ev_p"] * 10**8 + np.minimum(Z["ev_a"], Z["ev_b"]) * 10**4 + np.maximum(Z["ev_a"], Z["ev_b"])
ka, kb = keys(A), keys(B)
ia, ib = np.argsort(ka), np.argsort(kb)
assert np.array_equal(ka[ia], kb[ib])
ga, gb = A["ev_g"][ia], B["ev_g"][ib]
gmax = ga.max()
forb = ga < 1e-12 * gmax
out = {"n_events": int(len(ga)), "n_forbidden_unsym": int(forb.sum()),
       "forbidden_unsym_log10rel_quantiles": np.quantile(np.log10(ga[forb] / gmax), [0, .25, .5, .75, 1]).tolist(),
       "forbidden_sym_log10rel_quantiles": np.quantile(np.log10(np.maximum(gb[forb], 1e-300) / gmax), [0, .25, .5, .75, 1]).tolist(),
       "n_forbidden_sym_above_1e-12": int(np.sum(gb[forb] > 1e-12 * gmax)),
       "n_allowed_unsym_becoming_forbidden_sym": int(np.sum((~forb) & (gb < 1e-12 * gmax))),
       "allowed_rel_change_median": float(np.median(np.abs(gb[~forb] / ga[~forb] - 1))),
       "allowed_rel_change_max": float(np.max(np.abs(gb[~forb] / ga[~forb] - 1)))}
# the 1e-20 cluster vs the 1e-30 cluster
mid = forb & (ga > 1e-27 * gmax)
low = forb & (ga <= 1e-27 * gmax)
out["cluster_1e-25_to_1e-15: n"] = int(mid.sum())
out["cluster_1e-25_to_1e-15: sym log10rel median"] = float(np.median(np.log10(np.maximum(gb[mid], 1e-300) / gmax)))
out["cluster_<1e-27: n"] = int(low.sum())
out["cluster_<1e-27: sym log10rel median"] = float(np.median(np.log10(np.maximum(gb[low], 1e-300) / gmax)))
print(json.dumps(out, indent=1))
json.dump(out, open(base + r"\review\cx_selection_symfc_compare.json", "w"), indent=1)
