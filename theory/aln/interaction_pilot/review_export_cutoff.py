"""Read saved exports only; compare squared-interaction blocks and eigenbases."""
from pathlib import Path
import json,itertools,numpy as np,h5py
ROOT=Path(__file__).resolve().parent
records=json.loads((ROOT/'export-cutoff-comparison.json').read_text())
data={}
for record in records:
 f=record['factor'];run=Path(record['run'])
 def read(name):
  with h5py.File(run/name) as z:return {k:z[k][()] for k in z}
 data[f]={'phonon':read('phonon-m333.hdf5'),'pp':read('pp-m333-g1-s0.1.hdf5'),'gamma':read('kappa-m333-g1-s0.1.hdf5')['gamma'],
          'metadata':json.loads((run.parents[1]/'pilot-run.json').read_text())}
a,b=data[1.5],data[2.0];fa=a['phonon']['frequency'];fb=b['phonon']['frequency'];va=a['phonon']['eigenvector'];vb=b['phonon']['eigenvector']
pa=a['pp']['pp'];pb=b['pp']['pp'];trip=a['pp']['triplet'];delta=np.abs(pb-pa);scale=float(pa.max());maximum=np.unravel_index(np.argmax(delta),delta.shape)
result={'scope':'Saved-array test, not a new material calculation','input_hashes_identical':a['metadata']['inputs']==b['metadata']['inputs'],
 'addresses_identical':bool(np.array_equal(a['phonon']['grid_address'],b['phonon']['grid_address'])),
 'triplets_identical':bool(np.array_equal(trip,b['pp']['triplet'])),'weights_identical':bool(np.array_equal(a['pp']['weight'],b['pp']['weight'])),
 'raw_pp_max_scaled_change':float(delta.max()/scale),'raw_pp_L1_relative_change':float(delta.sum()/pa.sum()),
 'gamma_max_branch_relative_change':float(np.max(np.abs((b['gamma']-a['gamma'])/a['gamma']))),
 'maximum_pp_difference':{'index':[int(x) for x in maximum],'q_rows':trip[maximum[0]].tolist(),'pa':float(pa[maximum]),'pb':float(pb[maximum]),
     'frequencies_a':[float(fa[trip[maximum[0],leg],maximum[leg+1]]) for leg in range(3)]},'threshold_studies':[]}

def clusters(q,tol):
 means=(fa[q]+fb[q])/2
 inactive=np.where((fa[q]<1e-4)&(fb[q]<1e-4))[0].tolist()
 groups=[inactive] if inactive else []
 active=[i for i in range(len(means)) if i not in inactive]
 for i in active:
  if groups and groups[-1][-1] not in inactive and abs(means[i]-means[groups[-1][-1]])<tol:groups[-1].append(i)
  else:groups.append([i])
 return groups
for tol in (1e-12,1e-10,1e-8):
 group=[clusters(q,tol) for q in range(len(fa))];maps=[];overlaps=[]
 for q,gs in enumerate(group):
  m={i:g for g in gs for i in g};maps.append(m)
  for g in gs:
   if len(g)>1 and np.min(fa[q,g])>=1e-4:
    aa=va[q][:,g];bb=vb[q][:,g];o=aa.conj().T@bb;svals=np.linalg.svd(o,compute_uv=False)
    off=np.abs(o)**2;np.fill_diagonal(off,0)
    overlaps.append({'q_row':q,'address':a['phonon']['grid_address'][q].tolist(),'bands':g,
      'frequency_span_a_THz':float(np.ptp(fa[q,g])),'frequency_span_b_THz':float(np.ptp(fb[q,g])),
      'max_offdiagonal_overlap_squared':float(off.max()),'minimum_overlap_singular_value':float(svals.min()),
      'projector_spectral_difference':float(np.linalg.norm(aa@aa.conj().T-bb@bb.conj().T,2))})
 blocka=[];blockb=[];worst=None
 for t,qs in enumerate(trip):
  for gg in itertools.product(*(group[q] for q in qs)):
   idx=np.ix_(*gg);sa=float(pa[t][idx].sum());sb=float(pb[t][idx].sum());blocka.append(sa);blockb.append(sb)
   if worst is None or abs(sb-sa)>worst['absolute_change']:worst={'triplet_index':t,'bands':gg,'sum_a':sa,'sum_b':sb,'absolute_change':abs(sb-sa)}
 blocka=np.array(blocka);blockb=np.array(blockb);bd=np.abs(blocka-blockb);selected=blocka>1e-12*blocka.max()
 significant=np.argwhere(delta>1e-10*scale);not_degenerate=[]
 for ix in significant:
  t,i,j,k=map(int,ix)
  if all(len(maps[int(q)][bb])==1 for q,bb in zip(trip[t],(i,j,k))):not_degenerate.append([t,i,j,k])
 result['threshold_studies'].append({'frequency_cluster_tolerance_THz':tol,'block_count':len(blocka),
  'block_max_change_scaled_by_max_block':float(bd.max()/blocka.max()),'block_L1_relative_change':float(bd.sum()/blocka.sum()),
  'max_block_relative_change_above_declared_floor':float(np.max(bd[selected]/blocka[selected])),
  'relative_block_floor':1e-12,'maximum_block_change':worst,
  'raw_entries_above_1e_10_global_scale':len(significant),'significant_entries_with_only_nondegenerate_legs':not_degenerate,
  'rotating_active_clusters':[v for v in overlaps if v['max_offdiagonal_overlap_squared']>1e-5],
  'max_active_cluster_projector_difference':max(v['projector_spectral_difference'] for v in overlaps),
  'cluster_of_each_leg_at_maximum':[maps[int(q)][bb] for q,bb in zip(trip[maximum[0]],maximum[1:])]})
print(json.dumps(result,indent=2))
