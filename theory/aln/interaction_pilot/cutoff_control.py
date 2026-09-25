"""Hold Gaussian splitting parameter fixed while enlarging reciprocal cutoff."""
from pathlib import Path
import os,json,numpy as np,phono3py
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT/"runs/gp1-m333")
def load():
 return phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
  phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",fc3_filename="fc3.hdf5",
  born_filename="BORN",is_nac=True,symmetrize_fc=False,lang="Rust")
ph=load();ph.mesh_numbers=[3,3,3];ph.init_phph_interaction()
dm=ph.dynamical_matrix
dm.run(np.array([-1.,-1.,1.])/3)
_,_,cutoff,_,splitting=dm.Gonze_nac_dataset
rows=[]
for factor in [1.,1.25,1.5,2.]:
 ph=load()
 params=dict(ph.nac_params)
 params.update(G_cutoff=cutoff*factor,Lambda=splitting)
 ph.nac_params=params
 ph.mesh_numbers=[3,3,3];ph.init_phph_interaction()
 dm=ph.dynamical_matrix
 ds=[]
 for q in [np.array([-1.,-1.,1.])/3,np.array([-1.,2.,1.])/3]:
  dm.run(q);ds.append(dm.dynamical_matrix.copy())
 phase=np.repeat(np.exp(2j*np.pi*ph.primitive.scaled_positions[:,1]),3)
 ev=[np.linalg.eigvalsh(d) for d in ds]
 fs=[np.sign(e)*np.sqrt(abs(e))*ph.unit_conversion_factor for e in ev]
 _,_,actual_cutoff,gs,actual_lambda=dm.Gonze_nac_dataset
 rows.append(dict(cutoff_factor=factor,G_cutoff=actual_cutoff,Lambda=actual_lambda,
  G_count=len(gs),matrix_relative_defect=float(np.linalg.norm(phase[:,None]*ds[1]*phase.conj()[None,:]-ds[0])/np.linalg.norm(ds[0])),
  frequency_max_difference_THz=float(np.max(abs(fs[0]-fs[1])))))
 np.savez_compressed(ROOT/f"cutoff-control-{factor}.npz",d0=ds[0],d1=ds[1],phase=phase)
(ROOT/"cutoff-control.json").write_text(json.dumps(rows,indent=2)+"\n")
print(json.dumps(rows,indent=2))
