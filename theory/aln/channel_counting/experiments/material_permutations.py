"""Selected material tensor permutation audit, no kinetic/channel assignment."""
from pathlib import Path
import os,json,itertools,hashlib,argparse,importlib.metadata
import numpy as np,h5py,phono3py
from phono3py.phonon3.real_to_reciprocal import RealToReciprocal
from frequency_blocks import frequency_blocks
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--pilot-root",type=Path,required=True)
p.add_argument("--output-dir",type=Path,required=True)
a=p.parse_args();root=a.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
run=a.pilot_root.resolve()/"runs/export-cutoff-2.0/runs/gp1-m333"
meta=json.loads((run.parents[1]/"pilot-run.json").read_text())
for name,item in {**meta["inputs"],**meta["outputs"]}.items():
 assert hashlib.sha256((run/name).read_bytes()).hexdigest()==item["sha256"],name
assert importlib.metadata.version("phono3py")=="4.5.0"
os.chdir(run)
ph=phono3py.load(unitcell_filename="POSCAR",supercell_matrix=[3,3,2],
 phonon_supercell_matrix=[5,5,3],fc2_filename="fc2.hdf5",fc3_filename="fc3.hdf5",
 born_filename="BORN",is_nac=True,symmetrize_fc=False,lang="Rust")
ph.mesh_numbers=[3,3,3];ph.init_phph_interaction();interaction=ph.phph_interaction
with h5py.File("phonon-m333.hdf5") as f:
 frequencies=f["frequency"][:];vectors=f["eigenvector"][:];addresses=f["grid_address"][:]
with h5py.File("pp-m333-g1-s0.1.hdf5") as f:trip=f["triplet"][4];reference_pp=f["pp"][4]
r2r=RealToReciprocal(ph.fc3,ph.primitive,np.array([3,3,3]),
 make_r0_average=interaction._make_r0_average,
 all_shortest=interaction._all_shortest if interaction._make_r0_average else None)
mass=np.sqrt(ph.primitive.masses)
def amplitude(points):
 f=frequencies[points]
 if np.any(f<=1e-4):raise ValueError("selected audit excludes cutoff/zero modes")
 e=[v.reshape(4,3,12)/mass[:,None,None] for v in vectors[points]]
 r2r.run(addresses[points])
 raw=np.einsum("ijklmn,ila,jmb,knc->abc",r2r.get_fc3_reciprocal(),*e,optimize=True)
 return raw/np.sqrt(np.einsum("a,b,c->abc",*f))*np.sqrt(interaction._unit_conversion)
rows=[];tensors=[]
for perm in itertools.permutations(range(3)):
 v=amplitude(trip[list(perm)])
 aligned=v.transpose(np.argsort(perm));tensors.append(aligned)
 if not rows:reference=aligned.copy()
 rows.append({"permutation":perm,"amplitude_relative_frobenius":float(np.linalg.norm(aligned-reference)/np.linalg.norm(reference)),
              "pp_global_scaled_difference":float(np.max(abs(abs(aligned)**2-reference_pp))/np.max(reference_pp))})
 print(rows[-1],flush=True)
reverse=[]
for q in addresses[trip]:
 candidates=np.flatnonzero(np.all(addresses==-q,axis=1))
 reverse.append(int(candidates[0]) if len(candidates) else None)
record={"scope":"one all-incoming triplet, exact saved address reversals only; not Hamiltonian/channel normalization",
 "triplet":trip.tolist(),"addresses":addresses[trip].tolist(),"permutations":rows,
 "sum_vs_six_reference_relative":float(np.linalg.norm(sum(tensors)-6*reference)/np.linalg.norm(6*reference)),
 "exact_reverse_rows":reverse,"unit_conversion":interaction._unit_conversion,
 "make_r0_average":interaction._make_r0_average,"symmetrize_fc3q":interaction._symmetrize_fc3q}
payload={f"aligned_{i}":v for i,v in enumerate(tensors)}
if all(v is not None for v in reverse):
 rev=np.array(reverse);negative=amplitude(rev)
 sewing=[vectors[q].T@vectors[r] for q,r in zip(trip,rev)]
 # Undo frequency normalization before a general basis transformation.
 raw_conjugate=reference.conj()*np.sqrt(np.einsum("a,b,c->abc",*frequencies[trip]))
 predicted=np.einsum("ijk,ia,jb,kc->abc",raw_conjugate,*sewing,optimize=True)
 predicted/=np.sqrt(np.einsum("a,b,c->abc",*frequencies[rev]))
 block_sewing=[];block_groups=[]
 for u,fa,fb in zip(sewing,frequencies[trip],frequencies[rev]):
  groups=frequency_blocks(fa);other=frequency_blocks(fb)
  if [g.tolist() for g in groups]!=[g.tolist() for g in other]:
   raise ValueError("reversed points have different numerical block partitions")
  b=np.zeros_like(u)
  for g in groups:
   left,sv,right=np.linalg.svd(u[np.ix_(g,g)])
   if sv.min()<1-1e-8: raise ValueError("time-reversed eigenspaces are not aligned")
   b[np.ix_(g,g)]=left@right
  block_sewing.append(b);block_groups.append([g.tolist() for g in groups])
 block_predicted=np.einsum("ijk,ia,jb,kc->abc",raw_conjugate,*block_sewing,optimize=True)
 block_predicted/=np.sqrt(np.einsum("a,b,c->abc",*frequencies[rev]))
 record["time_reversal"]={"amplitude_relative_frobenius":float(np.linalg.norm(predicted-negative)/np.linalg.norm(negative)),
 "block_only_relative_frobenius":float(np.linalg.norm(block_predicted-negative)/np.linalg.norm(negative)),
 "block_groups":block_groups,"frequency_max_difference_THz":float(np.max(abs(frequencies[trip]-frequencies[rev]))),
 "unitarity_residuals":[float(np.linalg.norm(u.conj().T@u-np.eye(12))) for u in sewing]}
 payload.update(reversed_amplitude=negative,reversed_prediction=predicted,
                sewing0=sewing[0],sewing1=sewing[1],sewing2=sewing[2],block_prediction=block_predicted,
 block_sewing0=block_sewing[0],block_sewing1=block_sewing[1],block_sewing2=block_sewing[2],
 original_frequencies=frequencies[trip],reversed_frequencies=frequencies[rev])
f0,f1,f2=frequencies[trip]
mask=reference_pp>np.max(reference_pp)*1e-12
detuning=f0[:,None,None]-f1[None,:,None]-f2[None,None,:]
record["energy_admission"]={"orientation":"first leg parent, reversed other two legs; no quadrature conversion",
 "pp_relative_floor":1e-12,"minimum_absolute_decay_detuning_THz":float(np.min(abs(detuning[mask]))),
 "number_above_floor":int(mask.sum()),"number_with_absolute_detuning_below_1e_10_THz":int(np.count_nonzero(mask&(abs(detuning)<1e-10)))}
record["helper_sha256"]=hashlib.sha256((Path(__file__).parent/"frequency_blocks.py").read_bytes()).hexdigest()
record["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(root/"permutation-audit.json").write_text(json.dumps(record,indent=2)+"\n")
np.savez_compressed(root/"permutation-audit.npz",**payload)
print(json.dumps(record,indent=2))
