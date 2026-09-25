"""Independent array/accounting checks; not material convergence."""
from pathlib import Path
import json,numpy as np,h5py
ROOT=Path(__file__).resolve().parent
RUN=ROOT/"runs/gp1-m333"
def read(name):
    with h5py.File(RUN/name) as f:return {k:f[k][()] for k in f}
phon=read("phonon-m333.hdf5"); pp=read("pp-m333-g1-s0.1.hdf5")
detail=read("gamma_detail-m333-g1-s0.1.hdf5"); row=read("kappa-m333-g1-s0.1.hdf5")
mesh=phon["mesh"]; addr=phon["grid_address"];freq=phon["frequency"]
keys=[tuple(a%mesh) for a in addr]
lookup={}
for i,key in enumerate(keys):lookup.setdefault(key,i)
representatives=np.array(list(lookup.values()))
rev=np.array([lookup[tuple((-addr[i])%mesh)] for i in representatives])
quotient={key:i for i,key in enumerate(lookup)}
involution=np.array([quotient[tuple((-np.array(key))%mesh)] for key in lookup])
reconstructed=(detail["gamma_detail"].sum(axis=(3,4))*detail["weight"][None,:,None]).sum(axis=1)
scale=np.max(np.abs(row["gamma"]))
triplet=pp["triplet"]; full=pp["triplet_all"]
mismatch=addr[triplet].sum(axis=1)
delta=freq[triplet[:,0],:,None,None]-freq[triplet[:,1],None,:,None]-freq[triplet[:,2],None,None,:]
active=pp["pp"]>0
eigs=phon["eigenvector"]
identity=np.eye(freq.shape[1])
result={
 "scope":"array/serialization/accounting pilot, no material operator or continuum convergence",
 "grid":{"bz_entries":len(addr),"unique_regular_points":len(lookup),"expected_regular_points":int(np.prod(mesh)),
         "involution":bool(np.array_equal(involution[involution],np.arange(len(involution)))),
         "reversal_fixed_points":int(np.sum(involution==np.arange(len(involution)))),
         "reversal_frequency_max_abs_THz":float(np.max(np.abs(freq[representatives]-freq[rev]))),
         "duplicate_frequency_max_abs_THz":float(max(np.max(np.abs(freq[i]-freq[lookup[key]])) for i,key in enumerate(keys)))},
 "triplets":{"reduced":len(triplet),"full":len(full),"weight_sum":int(pp["weight"].sum()),
             "all_fixed_target":bool(np.all(full[:,0]==1)),
             "sum_mod_mesh_max":int(np.max(np.abs(mismatch%mesh))),
             "full_sum_mod_mesh_max":int(np.max(np.abs(addr[full].sum(axis=1)%mesh))),
             "unique_full_second_q":len({keys[i] for i in full[:,1]})},
 "arrays":{"pp_shape":list(pp["pp"].shape),"pp_min":float(pp["pp"].min()),
           "pp_max":float(pp["pp"].max()),"pp_nonzero_count":int(active.sum()),
           "all_pp_finite":bool(np.isfinite(pp["pp"]).all()),
           "phonon_eigenvector_orthogonality_max":float(np.max(np.abs(eigs.conj().transpose(0,2,1)@eigs-identity)))},
 "rates":{"weighted_detail_max_abs_error_THz":float(np.max(np.abs(reconstructed-row["gamma"]))),
          "weighted_detail_max_scaled_error":float(np.max(np.abs(reconstructed-row["gamma"]))/scale),
          "N_plus_U_max_abs_error_THz":float(np.max(np.abs(row["gamma_N"]+row["gamma_U"]-row["gamma"]))),
          "gamma_min_THz":float(row["gamma"].min()),"gamma_max_THz":float(row["gamma"].max())},
 "energy_mismatch":{"orientation":"f0-f1-f2 only; not all channels and not an event-weighted physical statistic",
                    "among_nonzero_pp_min_abs_THz":float(np.abs(delta[active]).min()),
                    "among_nonzero_pp_max_abs_THz":float(np.abs(delta[active]).max()),
                    "fraction_nonzero_pp_with_abs_delta_above_1e_8_THz":float(np.mean(np.abs(delta[active])>1e-8))}}
checks={
 "grid_cardinality":len(lookup)==int(np.prod(mesh)),
 "reversal_involution":result["grid"]["involution"],
 "reversal_fixed_count":result["grid"]["reversal_fixed_points"]==1,
 "reciprocal_frequencies":result["grid"]["reversal_frequency_max_abs_THz"]<1e-10,
 "BZ_duplicate_frequencies":result["grid"]["duplicate_frequency_max_abs_THz"]<1e-10,
 "weights_cover_mesh":int(pp["weight"].sum())==int(np.prod(mesh)),
 "momentum_mod_grid":result["triplets"]["sum_mod_mesh_max"]==0 and result["triplets"]["full_sum_mod_mesh_max"]==0,
 "full_triplet_coverage":result["triplets"]["unique_full_second_q"]==int(np.prod(mesh)),
 "positive_finite_pp":bool(np.isfinite(pp["pp"]).all() and pp["pp"].min()>=0),
 "phonon_orthonormality":result["arrays"]["phonon_eigenvector_orthogonality_max"]<1e-10,
 "detail_reconstructs_gamma":result["rates"]["weighted_detail_max_scaled_error"]<1e-12,
 "NU_reconstructs_gamma":result["rates"]["N_plus_U_max_abs_error_THz"]<1e-12*scale}
checks={k:bool(v) for k,v in checks.items()}
result["checks"]=checks
(ROOT/"pilot-inspection.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result,indent=2))
assert all(checks.values()),checks
