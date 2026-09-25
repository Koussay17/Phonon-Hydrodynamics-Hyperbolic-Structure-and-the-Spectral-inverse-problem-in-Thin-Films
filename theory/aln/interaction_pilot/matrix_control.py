"""Four-matrix periodic-gauge diagnostic requested by numerical reviewer."""
from pathlib import Path
import json,os,numpy as np
import phono3py
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT/"runs/gp1-m333")
q0=np.array([-1.,-1.,1.])/3;q1=np.array([-1.,2.,1.])/3
rows=[]
for nac in [False,True]:
 ph=phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
     phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",fc3_filename="fc3.hdf5",
     born_filename="BORN" if nac else None,is_nac=nac,symmetrize_fc=False,lang="Rust")
 ph.mesh_numbers=[3,3,3];ph.init_phph_interaction()
 dm=ph.dynamical_matrix
 matrices=[];residual=[];frequencies=[]
 for q in [q0,q1]:
  dm.run(q)
  d=dm.dynamical_matrix.copy();matrices.append(d)
  ev,v=np.linalg.eigh(d)
  residual.append(float(np.linalg.norm(d@v-v*ev)/np.linalg.norm(d)))
  frequencies.append(np.sign(ev)*np.sqrt(abs(ev))*ph.unit_conversion_factor)
 pos=ph.primitive.scaled_positions
 phase=np.repeat(np.exp(2j*np.pi*(pos@(q1-q0))),3)
 d0,d1=matrices
 aligned=phase[:,None]*d1*phase.conj()[None,:]
 opposite=phase.conj()[:,None]*d1*phase[None,:]
 rows.append({"NAC":nac,"actual_NAC_loaded":ph.nac_params is not None,
   "matrix_class":type(dm).__name__,"nac_method":getattr(dm,"nac_method",None),
   "eigenpair_relative_residuals":residual,
   "frequency_max_difference_THz":float(np.max(abs(frequencies[0]-frequencies[1]))),
   "gauge_aligned_matrix_relative_difference":float(np.linalg.norm(aligned-d0)/np.linalg.norm(d0)),
   "opposite_gauge_relative_difference":float(np.linalg.norm(opposite-d0)/np.linalg.norm(d0))})
 np.savez_compressed(ROOT/f"matrix-control-nac{int(nac)}.npz",d0=d0,d1=d1,phase=phase,q0=q0,q1=q1)
(ROOT/"matrix-control.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
print(json.dumps(rows,indent=2))
