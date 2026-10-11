"""E19: mesh transfer of 5x5x3 witness tensors, evaluated with the campaign's own operator (fc3class.FC3Class) on the
7x7x5 event set. Amplitudes from p3_events.py --fc3 on the same event list (results/m775x).
usage: python -B mesh_transfer_check.py <events_m775.npz> <dir with pp_r01_*.npz> <out.json> name ...
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class  # noqa: E402
from witness_verify import degavg_groups  # noqa: E402

ev, d, outp, names = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4:]
M = FC3Class(ev, d / "pp_r01_REF.npz")
grp = degavg_groups(M)
da = lambda r: (np.bincount(grp, weights=r) / np.bincount(grp))[grp]
m0 = M.metrics(M.rates(np.load(d / "pp_r01_REF.npz")["pp"]))
FOURPI = 4 * np.pi
res = {"REF": {"kappa": [m0["kappa_x"], m0["kappa_z"]], "tau_ps": [m0["tau_x"], m0["tau_z"]]}}
for nm in names:
    m = M.metrics(M.rates(np.load(d / f"pp_r01_{nm}.npz")["pp"]))
    lr = da(m["r"]) / da(m0["r"]) - 1
    res[nm] = {"kappa_ratio": [m["kappa_x"] / m0["kappa_x"], m["kappa_z"] / m0["kappa_z"]],
               "tau_ratio": [m["tau_x"] / m0["tau_x"], m["tau_z"] / m0["tau_z"]],
               "lifetime_degavg_minmax": [float(lr.min()), float(lr.max())]}
    print(nm, json.dumps(res[nm]), flush=True)
print("REF", res["REF"])
outp.write_text(json.dumps(res, indent=1))
