"""One-triplet complex contraction from pinned FC3 and saved phonons."""
from pathlib import Path
import json,os,hashlib
import numpy as np,h5py,phono3py
from phono3py.phonon3.real_to_reciprocal import RealToReciprocal
ROOT=Path(__file__).resolve().parents[1]
BASE=Path(r"D:\ResearchLab\orchestration\campaigns\20260924-105051-aln-interaction-pilot")
runs={v:BASE/"runs"/f"export-cutoff-{v}"/"runs/gp1-m333" for v in (1.5,2.0)}
os.chdir(runs[1.5])
ph=phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
 phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",fc3_filename="fc3.hdf5",
 born_filename="BORN",is_nac=True,symmetrize_fc=False,lang="Rust")
ph.mesh_numbers=[3,3,3];ph.init_phph_interaction()
interaction=ph.phph_interaction
arrays={}
for factor,run in runs.items():
 with h5py.File(run/"phonon-m333.hdf5") as f:
  freq=f["frequency"][:];eig=f["eigenvector"][:];addr=f["grid_address"][:]
 with h5py.File(run/"pp-m333-g1-s0.1.hdf5") as f:
  trip=f["triplet"][4];pp=f["pp"][4]
 arrays[factor]=(freq[trip],eig[trip],addr[trip],pp)
assert np.array_equal(arrays[1.5][2],arrays[2.][2])
r2r=RealToReciprocal(ph.fc3,ph.primitive,np.array([3,3,3]),
 make_r0_average=interaction._make_r0_average,
 all_shortest=interaction._all_shortest if interaction._make_r0_average else None)
r2r.run(arrays[1.5][2]);fc=r2r.get_fc3_reciprocal()
mass=np.sqrt(ph.primitive.masses)
amplitudes={};bare={};rows=[]
for factor,(freq,eig,addr,pp) in arrays.items():
 e=[v.reshape(4,3,12)/mass[:,None,None] for v in eig]
 raw=np.einsum("ijklmn,ila,jmb,knc->abc",fc,*e,optimize=True)
 amp=raw/np.sqrt(np.einsum("a,b,c->abc",*freq))*np.sqrt(interaction._unit_conversion)
 prediction=abs(amp)**2
 amplitudes[factor]=amp;bare[factor]=raw
 rows.append(dict(factor=factor,pp_max_scaled_error=float(np.max(abs(prediction-pp))/np.max(pp)),
  pp_L1_relative_error=float(np.sum(abs(prediction-pp))/np.sum(pp))))
ea=arrays[1.5][1];eb=arrays[2.][1]
u=[a.conj().T@b for a,b in zip(ea,eb)]
predicted=np.einsum("ijk,ia,jb,kc->abc",bare[1.5],*u,optimize=True)
f=arrays[2.][0]
predicted_amp=predicted/np.sqrt(np.einsum("a,b,c->abc",*f))*np.sqrt(interaction._unit_conversion)
block_u=[]
for v,frequencies in zip(u,arrays[1.5][0]):
 grouped=np.zeros_like(v);used=set()
 for j in range(12):
  if j in used: continue
  group=np.flatnonzero(abs(frequencies-frequencies[j])<1e-10)
  used.update(map(int,group))
  left,sv,right=np.linalg.svd(v[np.ix_(group,group)])
  grouped[np.ix_(group,group)]=left@right
 block_u.append(grouped)
block_prediction=np.einsum("ijk,ia,jb,kc->abc",amplitudes[1.5],*block_u,optimize=True)
result=dict(scope="one selected AlN triplet, shared real-to-reciprocal source, independent complex einsum contraction; not collision action",
 triplet_index=4,triplet_grid_rows=trip.tolist(),make_r0_average=interaction._make_r0_average,
 symmetrize_fc3q=interaction._symmetrize_fc3q,unit_conversion=interaction._unit_conversion,
 comparisons=rows,block_only_transport_relative_frobenius=float(np.linalg.norm(block_prediction-amplitudes[2.])/np.linalg.norm(amplitudes[2.])),
 off_block_overlap_norms=[float(np.linalg.norm(a-b)) for a,b in zip(u,block_u)],
 complex_transport_relative_frobenius=float(np.linalg.norm(predicted_amp-amplitudes[2.])/np.linalg.norm(amplitudes[2.])),
 unitarity_residuals=[float(np.linalg.norm(v.conj().T@v-np.eye(12))) for v in u])
np.savez_compressed(ROOT/"experiments/complex-triplet.npz",amplitude_a=amplitudes[1.5],amplitude_b=amplitudes[2.],predicted_b=predicted_amp,
 u0=u[0],u1=u[1],u2=u[2],freq_a=arrays[1.5][0],freq_b=arrays[2.][0],block_prediction=block_prediction,block_u0=block_u[0],block_u1=block_u[1],block_u2=block_u[2])
result["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(ROOT/"experiments/complex-triplet.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
