"""Transport and memory moments of the hybrid model fc2 (phono3py Olympics) + fc3 (ALAMODE Olympics, folded and
projected, alamode_compare.py) on the fixed 5x5x3 event list; comparison with the reference (orbit-average
convention, lifetimes compared after degeneracy averaging).
usage: python -B alamode_transport.py <events.npz> <resultdir(m553)> <pp_alamode.npz> <out.json>
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from fc3class import FC3Class, summary  # noqa: E402
from witness_verify import degavg_groups, tau_routes  # noqa: E402

ev, rdir, ppA, outp = sys.argv[1], Path(sys.argv[2]), sys.argv[3], sys.argv[4]
M = FC3Class(ev, rdir / "pp_r01_REF.npz")
g0 = M.rates(np.load(rdir / "pp_r01_REF.npz")["pp"])
gA = M.rates(np.load(ppA)["pp"])
m0, mA = M.metrics(g0), M.metrics(gA)
grp = degavg_groups(M)
da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
lr = da(mA["r"]) / da(m0["r"]) - 1
al = g0 > 1e-9 * g0.max()
ratio = gA[al] / g0[al]
out = {"reference": summary(m0), "alamode_fc3": summary(mA),
       "routes_x": tau_routes(M, gA, 0), "routes_z": tau_routes(M, gA, 2),
       "kappa_ratio": [mA["kappa_x"] / m0["kappa_x"], mA["kappa_z"] / m0["kappa_z"]],
       "tau_ratio": [mA["tau_x"] / m0["tau_x"], mA["tau_z"] / m0["tau_z"]],
       "tauCS_ratio": [mA["tauCS_x"] / m0["tauCS_x"], mA["tauCS_z"] / m0["tauCS_z"]],
       "lifetime_rel_change_degavg": {"min": float(lr.min()), "max": float(lr.max()),
                                      "rms": float(np.sqrt(np.mean(lr ** 2))),
                                      "frac_within_10pct": float(np.mean(np.abs(lr) <= 0.1))},
       "event_rate_ratio_quantiles(0.1%,1%,50%,99%,99.9%)": np.quantile(ratio, [1e-3, 1e-2, .5, .99, .999]).tolist(),
       "selection_rule_events_alamode_max_rel": float(gA[~al].max() / gA.max())}
Path(outp).write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
