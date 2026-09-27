"""Selected material tensor permutation audit, no kinetic/channel assignment."""
from pathlib import Path
import os,json,itertools,hashlib,argparse,importlib.metadata
import numpy as np,h5py,phono3py
from phono3py.phonon3.real_to_reciprocal import RealToReciprocal
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
 predicted=np.einsum("ijk,ia,jb,kc->abc",reference.conj(),*sewing,optimize=True)
 # Destination frequency normalization corrects only stored roundoff differences.
 predicted*=np.sqrt(np.einsum("a,b,c->abc",*frequencies[trip])/np.einsum("a,b,c->abc",*frequencies[rev]))
 record["time_reversal"]={"amplitude_relative_frobenius":float(np.linalg.norm(predicted-negative)/np.linalg.norm(negative)),
 "frequency_max_difference_THz":float(np.max(abs(frequencies[trip]-frequencies[rev]))),
 "unitarity_residuals":[float(np.linalg.norm(u.conj().T@u-np.eye(12))) for u in sewing]}
 payload.update(reversed_amplitude=negative,reversed_prediction=predicted,
                sewing0=sewing[0],sewing1=sewing[1],sewing2=sewing[2])
record["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(root/"permutation-audit.json").write_text(json.dumps(record,indent=2)+"\n")
np.savez_compressed(root/"permutation-audit.npz",**payload)
print(json.dumps(record,indent=2))
