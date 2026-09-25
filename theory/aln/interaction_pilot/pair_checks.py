"""Reciprocity and broadening controls, no convergence claim."""
from pathlib import Path
import h5py,numpy as np,json
ROOT=Path(__file__).resolve().parent
RUN=ROOT/"runs/pair-m333"
rows={}
for gp in [1,2]:
 for sigma in [.05,.1,.2]:
  with h5py.File(RUN/f"kappa-m333-g{gp}-s{sigma}.hdf5") as f:
   rows[gp,sigma]=f["gamma"][0]
result={"pair_grid_points":[1,2],"sigma_THz":[.05,.1,.2],
 "reciprocity_max_scaled_error_by_sigma":{str(s):float(np.max(abs(rows[1,s]-rows[2,s]))/np.max(abs(rows[1,s]))) for s in [.05,.1,.2]},
 "gp1_gamma_THz":{str(s):rows[1,s].tolist() for s in [.05,.1,.2]},
 "sigma_005_over_01_range":[float(np.min(rows[1,.05]/rows[1,.1])),float(np.max(rows[1,.05]/rows[1,.1]))],
 "sigma_02_over_01_range":[float(np.min(rows[1,.2]/rows[1,.1])),float(np.max(rows[1,.2]/rows[1,.1]))],
 "interpretation":"width sensitivity on one coarse grid, not a grid convergence study or a physical uncertainty interval"}
result["checks"]={"reciprocity":all(x<1e-10 for x in result["reciprocity_max_scaled_error_by_sigma"].values()),
 "positive_rates":all(bool(np.all(x>0)) for x in rows.values())}
(ROOT/"pair-checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result,indent=2))
assert all(result["checks"].values())
