"""Independent saved-array numerical attack of PI/C artifacts; no material run."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,itertools,sys
import numpy as np
import h5py
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent
PILOT=Path(r'D:\ResearchLab\orchestration\campaigns\20260924-105051-aln-interaction-pilot')
z=dict(np.load(OUT/'complex-triplet.npz'))
a,b=z['amplitude_a'],z['amplitude_b']
fa,fb=z['freq_a'],z['freq_b']
u=[z[f'u{i}'] for i in range(3)]
meta=json.loads((OUT/'complex-triplet.json').read_text())
provenance=json.loads((OUT/'complex-triplet-provenance.json').read_text())

def transport(tensor, rotations):
    value=tensor.copy()
    for axis,rotation in enumerate(rotations):
        value=np.moveaxis(np.tensordot(value,rotation,axes=(axis,0)),-1,axis)
    return value

def summary_difference(x,y):
    d=np.abs(x-y)
    scale=float(np.max(np.abs(y)))
    record={'absolute_max':float(d.max()),'global_max_scaled':float(d.max()/scale),
            'relative_frobenius':float(np.linalg.norm(d)/np.linalg.norm(y))}
    for floor in (1e-6,1e-12):
        mask=np.abs(y)>floor*scale
        record[f'max_entry_relative_above_{floor:g}_global_floor']=float(np.max(d[mask]/np.abs(y[mask])))
    return record

def pi_groups_and_polar(frequencies, overlap, tolerance):
    grouped=np.zeros_like(overlap);used=set();groups=[];svs=[]
    for j in range(len(frequencies)):
        if j in used:continue
        group=np.flatnonzero(abs(frequencies-frequencies[j])<tolerance)
        used.update(map(int,group));groups.append(group.tolist())
        left,s,right=np.linalg.svd(overlap[np.ix_(group,group)])
        svs.extend(s.tolist());grouped[np.ix_(group,group)]=left@right
    return grouped,groups,svs


weight_a=np.sqrt(fa[0,:,None,None]*fa[1,None,:,None]*fa[2,None,None,:])
weight_b=np.sqrt(fb[0,:,None,None]*fb[1,None,:,None]*fb[2,None,None,:])
full=transport(a*weight_a,u)/weight_b
results={'scope':'PI/C saved arrays only; no material recalculation or D self-review',
 'npz_shapes':{k:list(v.shape) for k,v in z.items()},
 'all_arrays_finite':all(bool(np.all(np.isfinite(v))) for v in z.values()),
 'minimum_frequency':float(min(fa.min(),fb.min())),
 'full_transport':summary_difference(full,b),
 'saved_full_prediction_recomputed_difference':summary_difference(full,z['predicted_b']),
 'hash_checks':{},'native_export_checks':[],'clustering':[]}
for name,digest in provenance['outputs'].items():
    results['hash_checks'][name]=hashlib.sha256((OUT/name).read_bytes()).hexdigest()==digest
for factor,amp,freq in ((1.5,a,fa),(2.0,b,fb)):
    run=PILOT/'runs'/f'export-cutoff-{factor}'/'runs/gp1-m333'
    with h5py.File(run/'pp-m333-g1-s0.1.hdf5') as h:
        pp=h['pp'][4];trip=h['triplet'][4]
    with h5py.File(run/'phonon-m333.hdf5') as h:
        native_f=h['frequency'][trip];native_e=h['eigenvector'][trip]
    results['native_export_checks'].append({'factor':factor,'triplet':trip.tolist(),
        'squared_amplitude':summary_difference(abs(amp)**2,pp),
        'frequency_matches_npz':bool(np.array_equal(freq,native_f)),
        'pp_file_hash_matches_provenance':hashlib.sha256((run/'pp-m333-g1-s0.1.hdf5').read_bytes()).hexdigest()==provenance['inputs'][str(factor)]['pp-m333-g1-s0.1.hdf5'],
        'phonon_file_hash_matches_provenance':hashlib.sha256((run/'phonon-m333.hdf5').read_bytes()).hexdigest()==provenance['inputs'][str(factor)]['phonon-m333.hdf5']})
    if factor==1.5:ea=native_e
    else:eb=native_e
results['native_overlap_recomputed_residuals']=[float(np.linalg.norm(v-aa.conj().T@bb)) for v,aa,bb in zip(u,ea,eb)]
for tol in (1e-12,1e-10,1e-8):
    rotations=[];legs=[]
    for leg,(freq,ov) in enumerate(zip(fa,u)):
        rot,groups,svs=pi_groups_and_polar(freq,ov,tol);rotations.append(rot)
        coverage=np.zeros(len(freq),int);mask=np.zeros_like(ov,dtype=bool)
        for group in groups:coverage[group]+=1;mask[np.ix_(group,group)]=True
        outside=np.abs(freq[:,None]-freq[None,:])[~mask]
        spans_b=[float(np.ptp(fb[leg,group])) for group in groups]
        legs.append({'groups':groups,'coverage':coverage.tolist(),'min_overlap_singular_value':min(svs),
            'max_overlap_singular_value':max(svs),'max_frequency_span_a':max(float(np.ptp(freq[g])) for g in groups),
            'max_frequency_span_b':max(spans_b),'minimum_interblock_gap':float(outside.min()),
            'pure_offblock_frobenius':float(np.linalg.norm(ov[~mask])),
            'full_overlap_minus_polar_frobenius':float(np.linalg.norm(ov-rot)),
            'polar_unitarity_residual':float(np.linalg.norm(rot.conj().T@rot-np.eye(len(freq))))})
    results['clustering'].append({'tolerance':tol,'legs':legs,'transport':summary_difference(transport(a,rotations),b)})
# Exact source rule attacked with a near-degenerate chain, not asserted in AlN.
tol=1e-10;f=np.array([1,1+.75*tol,1+1.5*tol]+list(range(2,11)),float)
v=np.eye(12,dtype=complex)
v[:3,:3]=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
rot,groups,sv=pi_groups_and_polar(f,v,tol)
coverage=np.zeros(12,int)
for g in groups:coverage[g]+=1
results['synthetic_clustering_counterexample']={'frequencies':f.tolist(),'groups':groups,'coverage':coverage.tolist(),
 'input_unitarity_residual':float(np.linalg.norm(v.conj().T@v-np.eye(12))),
 'output_unitarity_residual':float(np.linalg.norm(rot.conj().T@rot-np.eye(12))),
 'output_singular_values':np.linalg.svd(rot,compute_uv=False).tolist(),
 'scope':'Nontransitive overlapping-group failure of source algorithm; distinct from selected AlN triplet.'}
# C: independent matrix-function evaluation via the exact two nonzero eigenvalues.
c=json.loads((OUT/'C_lumpability_results.json').read_text())
W=np.diag([6/25,6/25,2,3/4]);B=np.array([[1,1,0,0],[0,0,1,0],[0,0,0,1.]])
H=W@B.T@np.linalg.inv(B@W@B.T);Q=np.eye(4)-H@B
incidence=np.array([[-1,0],[0,-1],[1,1],[1,1.]])
results['C_independent_checks']={}
for name,rates in (('equal',np.array([2.,2.])),('unequal',np.array([3.,1.]))):
    M=(3/5)*incidence@np.diag(rates)@incidence.T
    L=M@np.linalg.inv(W);G=B@L@H
    disc=np.sqrt((72/5)**2-4*(47/4)*np.prod(rates))
    lambdas=np.array([(72/5-disc)/2,(72/5+disc)/2])
    checks=[]
    for old in [r for r in c['dynamics'] if r['model']==name]:
        t=old['t'];y=np.expm1(-lambdas*t)/lambdas
        c2=(y[1]-y[0])/(lambdas[1]-lambdas[0]);c1=y[0]-c2*lambdas[0]
        exact_projection=B@(np.eye(4)+c1*L+c2*(L@L))@H
        closed=np.eye(3)+np.expm1(-(47/5)*t)/(47/5)*G
        error=np.linalg.norm(exact_projection-closed,2)/np.linalg.norm(exact_projection,2)
        checks.append({'t':t,'independent_relative_propagator_error':float(error),
            'difference_from_reported_error':float(abs(error-old['propagator_relative_error'])),
            'maximum_reported_state_difference':float(np.max(np.abs(exact_projection@np.array([.02,0,0])-old['total_actual'])))})
    results['C_independent_checks'][name]={'nonzero_eigenvalues':lambdas.tolist(),
        'closure_defect_frobenius':float(np.linalg.norm(B@L@Q)),
        'hidden_derivative':(-B@L@np.array([1.,-1,0,0])).tolist(),
        'second_derivative_defect_frobenius':float(np.linalg.norm(B@L@L@H-G@G)),
        'propagation':checks}
results['second_memory_files_available']=[p.name for p in OUT.glob('*memory*')]
results['diagnostic_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(OUT/'R-numerical-audit.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8')
print(json.dumps(results,indent=2))

