"""Diagnostic controls for BZ-equivalent phonon frequencies."""
from pathlib import Path
import os,json,time
import numpy as np
import phono3py
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT/"runs/gp1-m333")
rows=[]
for nac in [False,True]:
 for backend in ["Rust","C"]:
  t=time.time()
  try:
   ph=phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
       phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",fc3_filename="fc3.hdf5",
       born_filename="BORN" if nac else None,is_nac=nac,symmetrize_fc=False,lang=backend,log_level=0)
   ph.mesh_numbers=[3,3,3];ph.init_phph_interaction();ph.run_phonon_solver()
   f,v,a=ph.get_phonon_data()
   seen={};diff=[]
   for i,key in enumerate(map(tuple,a%3)):
    if key in seen:
     j=seen[key]
     diff.append(float(np.max(np.abs(f[i]-f[j]))))
    else:seen[key]=i
   np.savez_compressed(ROOT/f"harmonic-{backend}-nac{int(nac)}.npz",frequency=f,address=a)
   rows.append({"backend":backend,"NAC":nac,"actual_NAC_loaded":ph.nac_params is not None,"max_duplicate_difference_THz":max(diff),
       "min_frequency_THz":float(f.min()),"seconds":time.time()-t})
  except Exception as e:
   rows.append({"backend":backend,"NAC":nac,"error":type(e).__name__+": "+str(e)})
(ROOT/"harmonic-controls.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
print(json.dumps(rows,indent=2))
