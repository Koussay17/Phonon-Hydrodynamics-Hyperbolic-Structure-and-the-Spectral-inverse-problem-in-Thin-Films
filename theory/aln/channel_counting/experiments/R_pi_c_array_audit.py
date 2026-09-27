"""Independent PI/C saved-array audit; no D branch, solver, or material run."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import hashlib,itertools,json,math
import numpy as np
import h5py
ROOT=Path(__file__).resolve().parents[1]
EXP=ROOT/'experiments'
load=lambda p:json.loads(Path(p).read_text())
prov=load(EXP/'material-provenance.json');material=load(EXP/'permutation-audit.json')
channel=load(EXP/'channel-sewing.json');cf=load(EXP/'C_fock_check.json');c2=load(EXP/'C2_closure_check.json')
array_path=Path(prov['array_path'])
z=dict(np.load(array_path,allow_pickle=False));v=z['aligned_0'];vr=z['reversed_amplitude']
s=[z[f'block_sewing{i}'] for i in range(3)]

def transform(t,rotations):
    value=t.copy()
    for axis,u in enumerate(rotations):value=np.moveaxis(np.tensordot(value,u,axes=(axis,0)),-1,axis)
    return value

def metric(x,y,power_floor=1e-12):
    delta=abs(x-y);mask=abs(y)**2>power_floor*np.max(abs(y)**2)
    return {'relative_frobenius':float(np.linalg.norm(delta)/np.linalg.norm(y)),
        'absolute_max':float(delta.max()),'entries_above_relative_power_floor':int(mask.sum()),
        'relative_power_floor':power_floor,'maximum_entry_relative_above_floor':float(np.max(delta[mask]/abs(y[mask])))}

def channels(t,tr,ss):
    return 6*transform(t,[np.eye(t.shape[0]),ss[1].conj(),ss[2].conj()]),6*transform(tr,[ss[0].conj().T,np.eye(tr.shape[1]),np.eye(tr.shape[2])])

forward,inverse=channels(v,vr,s)
result={'scope':'Independent numerical attack of PI/C only; no D self-approval, no material run',
 'hash_checks':{'material_array':hashlib.sha256(array_path.read_bytes()).hexdigest()==prov['array_sha256'],
  'channel_array_reference':channel['material']['array_sha256']==prov['array_sha256']},
 'permutations':[metric(z[f'aligned_{i}'],v) for i in range(6)],
 'channel':metric(inverse,forward.conj()),'full_reverse_saved_prediction':metric(z['reversed_prediction'],vr),
 'block_reverse_saved_prediction':metric(z['block_prediction'],vr),
 'block_sewing_singular_values':[(np.linalg.svd(q,compute_uv=False).min().item(),np.linalg.svd(q,compute_uv=False).max().item()) for q in s],
 'block_unitarity_frobenius':[float(np.linalg.norm(q.conj().T@q-np.eye(len(q)))) for q in s]}
for name,rec in [('material_permutations.py',material),('channel_sewing_check.py',channel),('C_fock_check.py',cf),('C2_closure_check.py',c2)]:
 result['hash_checks'][name]=hashlib.sha256((EXP/name).read_bytes()).hexdigest()==rec['script_sha256']
run=Path(prov['previous_pilot'])/'runs/export-cutoff-2.0/runs/gp1-m333'
with h5py.File(run/'phonon-m333.hdf5') as q:addr=q['grid_address'][:];frequencies=q['frequency'][:];eig=q['eigenvector'][:]
with h5py.File(run/'pp-m333-g1-s0.1.hdf5') as q:pp=q['pp'][4];trip=q['triplet'][4]
reverse=np.array(material['exact_reverse_rows'])
result['address_and_input_checks']={'exact_negation':bool(np.array_equal(addr[reverse],-addr[trip])),
 'triplet_matches_record':trip.tolist()==material['triplet'],
 'saved_frequencies_match':bool(np.array_equal(frequencies[trip],z['original_frequencies']) and np.array_equal(frequencies[reverse],z['reversed_frequencies'])),
 'raw_sewing_recomputation_residuals':[float(np.linalg.norm(eig[i].T@eig[j]-z[f'sewing{k}'])) for k,(i,j) in enumerate(zip(trip,reverse))]}
f0,f1,f2=z['original_frequencies'];det=f0[:,None,None]-f1[None,:,None]-f2[None,None,:]
mask=pp>1e-12*pp.max()
result['energy_admission']={'above_floor_count':int(mask.sum()),'minimum_absolute_detuning_THz':float(abs(det[mask]).min()),
 'count_below_1e_10_THz':int(np.sum(mask&(abs(det)<1e-10)))}
# Fresh independent gauges, restricted to the recorded numerical blocks.
rng=np.random.default_rng(9272026)
def gauge(groups,n):
    u=np.zeros((n,n),complex)
    for g in groups:
        t=rng.normal(size=(len(g),len(g)))+1j*rng.normal(size=(len(g),len(g)))
        q,_=np.linalg.qr(t);u[np.ix_(g,g)]=q
    return u
result['independent_gauges']=[]
for trial in range(3):
    aa=[gauge(g,12) for g in material['time_reversal']['block_groups']]
    bb=[gauge(g,12) for g in material['time_reversal']['block_groups']]
    vv=transform(v,aa);vv_r=transform(vr,bb)
    ss=[a.T@q@b for a,q,b in zip(aa,s,bb)]
    f,g=channels(vv,vv_r,ss)
    cov=transform(forward,[aa[0],bb[1].conj(),bb[2].conj()])
    wrong=[a.conj().T@q@b for a,q,b in zip(aa,s,bb)]
    wf,wg=channels(vv,vv_r,wrong)
    result['independent_gauges'].append({'trial':trial,'channel_conjugation_relative':metric(g,f.conj())['relative_frobenius'],
        'forward_covariance_relative':metric(f,cov)['relative_frobenius'],
        'wrong_adjoint_in_sewing_relative':metric(wg,wf.conj())['relative_frobenius']})
# Untruncated ladder-path enumeration, independent of C's sparse Kronecker matrices.
def apply_field(states,index,reverse):
    out={}
    for state,coef in states.items():
        if state[index]:
            final=list(state);final[index]-=1;key=tuple(final)
            out[key]=out.get(key,0j)+coef*math.sqrt(state[index])
        final=list(state);final[reverse[index]]+=1;key=tuple(final)
        out[key]=out.get(key,0j)+coef*math.sqrt(state[reverse[index]]+1)
    return out

def h_element(initial,final,indices,reverse,duplicate=False):
    total=0j
    for incoming,coef in ((indices,1+2j),(tuple(reverse[i] for i in indices),1-2j)):
        permutations=list(itertools.permutations(incoming))
        if not duplicate:permutations=set(permutations)
        for order in permutations:
            states={tuple(initial):1+0j}
            for index in reversed(order):states=apply_field(states,index,reverse)
            total+=coef/6*states.get(tuple(final),0j)
    return total
result['C_fock_untruncated_samples']={}
for name in ('distinct','repeated'):
    data=cf[name];names=data['names'];reverse=[names.index(n[4:] if n.startswith('bar_') else 'bar_'+n) for n in names]
    errors=[];conjerrors=[]
    for sample in data['samples']:
        value=h_element(sample['initial'],sample['final'],tuple(data['incoming_tuple']),reverse)
        errors.append(abs(value-complex(*sample['matrix_amplitude'])))
        back=h_element(sample['final'],sample['initial'],tuple(data['incoming_tuple']),reverse)
        conjerrors.append(abs(back-value.conjugate()))
    result['C_fock_untruncated_samples'][name]={'maximum_difference_from_saved_sparse_sample':max(errors),'maximum_reverse_conjugation_error':max(conjerrors),'sample_count':len(errors)}
    if name=='repeated':
        sample=data['samples'][0]
        correct=h_element(sample['initial'],sample['final'],tuple(data['incoming_tuple']),reverse)
        wrong=h_element(sample['initial'],sample['final'],tuple(data['incoming_tuple']),reverse,True)
        result['C_fock_untruncated_samples'][name]['duplicate_permutation_squared_ratio']=abs(wrong/correct)**2
# Closed geometric tails from memorylessness, not finite state enumeration.
def rawmom(r,k):
    n=r/(1-r)
    return [F(1),n,2*n*n+n,6*n**3+6*n*n+n][k]

def partialmom(r,cut,k):
    tail=r**(cut+1)*sum(F(math.comb(k,j))*(cut+1)**(k-j)*rawmom(r,j) for j in range(k+1))
    return rawmom(r,k)-tail
result['C_thermal_closed_tail_check']=[]
for case in cf['thermal']['geometric_partial_sums']:
    r=F(case['geometric_ratio_exact']);n=r/(1-r)
    errors=[]
    for row in case['convergence']:
        cut=row['cutoff'];m=[partialmom(r,cut,k) for k in range(3)]
        a=m[2]-m[1];c=m[2]+3*m[1]+2*m[0]
        erra=1-a/(2*n*n);errc=1-c/(2*(n+1)**2)
        errors.append(max(abs(float(erra)-row['annihilation_relative_error']),abs(float(errc)-row['creation_relative_error'])))
    result['C_thermal_closed_tail_check'].append({'ratio':str(r),'maximum_difference_from_reported_error':max(errors),
        'final_annihilation_relative_tail':float(erra),'final_creation_relative_tail':float(errc),'final_omitted_mass':float(r**(cut+1))})
polys={
 'flux':{(1,1):4,(1,0):2,(0,2):-1,(0,1):1},
 'parent_factorial2_derivative':{(2,1):-8,(2,0):-4,(1,2):4,(1,1):4,(1,0):4},
 'daughter_factorial2_derivative':{(1,2):24,(1,1):8,(1,0):4,(0,3):-4,(0,2):10,(0,1):-6},
 'parent_daughter_derivative':{(2,1):8,(2,0):4,(1,2):-10,(1,1):-4,(1,0):-4,(0,3):1,(0,2):-3,(0,1):2}}
polys['parent_derivative']={k:-v for k,v in polys['flux'].items()};polys['daughter_derivative']={k:2*v for k,v in polys['flux'].items()}
def expectation(poly,mp,md):return sum(F(coef)*mp[i]*md[j] for (i,j),coef in poly.items())
rp,rd=F(1,4),F(1,3);targets={k:expectation(p,[rawmom(rp,j) for j in range(4)],[rawmom(rd,j) for j in range(4)]) for k,p in polys.items()}
result['C2_independent_polynomial_tail_check']={'exact_targets_match':all(str(targets[k])==val for k,val in c2['rational_partial_sums']['targets_exact'].items()),'cutoffs':[]}
for row in c2['rational_partial_sums']['convergence']:
    cp,cd=row['max_parent'],row['max_daughter'];mp=[partialmom(rp,cp,j) for j in range(4)];md=[partialmom(rd,cd,j) for j in range(4)]
    vals={k:expectation(p,mp,md) for k,p in polys.items()};errors={k:abs(vals[k]-targets[k]) for k in vals}
    result['C2_independent_polynomial_tail_check']['cutoffs'].append({'cutoffs':[cp,cd],
        'all_partial_sums_match_exactly':all(vals[k]==F(v) for k,v in row['partial_generator_expectations'].items()),
        'max_error_matches_exactly':max(errors.values())==F(row['maximum_absolute_error_exact']),
        'maximum_absolute_error':float(max(errors.values())),
        'maximum_error_over_omitted_probability':float(max(errors.values())/F(row['omitted_probability_exact']))})
J=targets['flux'];b=F(1,2);daughter_tangent=4*b*targets['daughter_derivative']
result['C2_tangency']={'flux':str(J),'daughter_moment_derivative':str(targets['daughter_factorial2_derivative']),
 'geometric_tangent':str(daughter_tangent),'defect':str(targets['daughter_factorial2_derivative']-daughter_tangent),
 'claimed_2_times_2b_plus1_times_J':str(2*(2*b+1)*J)}
result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(EXP/'R_pi_c_array_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
