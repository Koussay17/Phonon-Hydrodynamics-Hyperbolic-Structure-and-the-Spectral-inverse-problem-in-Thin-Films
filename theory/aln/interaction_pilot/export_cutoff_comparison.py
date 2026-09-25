"""Re-export selected interactions under a fixed-Lambda reciprocal-cutoff sweep."""
from pathlib import Path
import os,json,shutil,hashlib,subprocess,sys
import numpy as np,h5py
ROOT=Path(__file__).resolve().parent
base=json.loads((ROOT/"cutoff-control.json").read_text())[0]
source=ROOT/"runs/gp1-m333"
template=(ROOT/"pilot_run.py").read_text()
rows=[]
for factor in (1.,1.5,2.):
 work=ROOT/"runs"/f"export-cutoff-{factor}"
 work.mkdir(exist_ok=False)
 run=work/"runs/gp1-m333";run.mkdir(parents=True)
 for name in ("POSCAR","BORN","fc2.hdf5","fc3.hdf5"):
  shutil.copy2(source/name,run/name)
 inject=f'ph.nac_params=dict(ph.nac_params, G_cutoff={base["G_cutoff"]*factor!r}, Lambda={base["Lambda"]!r})\n'
 script=template.replace("ph.mesh_numbers=[3,3,3]",inject+"ph.mesh_numbers=[3,3,3]")
 (work/"pilot_run.py").write_text(script)
 env=dict(os.environ,OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="2",RAYON_NUM_THREADS="2")
 result=subprocess.run([sys.executable,"-B",str(work/"pilot_run.py")],capture_output=True,env=env,timeout=180)
 (work/"export.log").write_bytes(result.stdout+b"\n"+result.stderr)
 if result.returncode: raise RuntimeError(str(work))
 if factor>1:
  for check in ("inspect_pilot.py","reconstruct_gamma.py"):
   shutil.copy2(ROOT/check,work/check)
   checked=subprocess.run([sys.executable,"-B",str(work/check)],capture_output=True,env=env,timeout=180)
   (work/(check+".log")).write_bytes(checked.stdout+b"\n"+checked.stderr)
   if checked.returncode: raise RuntimeError(f"{check} failed in {work}")
 with h5py.File(run/"phonon-m333.hdf5") as f:
  freq=f["frequency"][:];address=f["grid_address"][:]
 pairs=[(i,j) for i in range(len(address)) for j in range(i) if np.all((address[i]-address[j])%3==0)]
 defect=max(float(np.max(abs(freq[i]-freq[j]))) for i,j in pairs)
 with h5py.File(run/"pp-m333-g1-s0.1.hdf5") as f: pp=f["pp"][:];triplets=f["triplet"][:]
 with h5py.File(run/"kappa-m333-g1-s0.1.hdf5") as f: gamma=f["gamma"][:]
 rows.append(dict(factor=factor,G_cutoff=base["G_cutoff"]*factor,Lambda=base["Lambda"],
  duplicate_pair_count=len(pairs),duplicate_frequency_max_THz=defect,
  original_tolerance_pass=defect<=1e-10,run=str(run),gamma=gamma.tolist()))
 if factor==1.:
  pp0=pp.copy();g0=gamma.copy();t0=triplets.copy()
 else:
  assert np.array_equal(t0,triplets)
  rows[-1].update(pp_max_scaled_change_from_default=float(np.max(abs(pp-pp0))/np.max(pp0)),
   gamma_max_scaled_change_from_default=float(np.max(abs(gamma-g0))/np.max(g0)),
   gamma_max_branch_relative_change_from_default=float(np.max(abs((gamma-g0)/g0))))
 if factor==1.5: pp15=pp.copy();g15=gamma.copy()
 if factor==2.:
  rows[-1].update(pp_max_scaled_change_from_1_5=float(np.max(abs(pp-pp15))/np.max(pp15)),
   gamma_max_scaled_change_from_1_5=float(np.max(abs(gamma-g15))/np.max(g15)),
   gamma_max_branch_relative_change_from_1_5=float(np.max(abs((gamma-g15)/g15))))
 print(json.dumps(rows[-1]),flush=True)
(ROOT/"export-cutoff-comparison.json").write_text(json.dumps(rows,indent=2)+"\n")
