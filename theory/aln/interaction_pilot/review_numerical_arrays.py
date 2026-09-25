"""Read-only inspection of saved harmonic controls; no phonon solve."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parent
out={'scope':'Saved-array comparisons only; no force-constant or phonon recalculation','controls':[],'backend_comparisons':[]}
arrays={}
for nac in (0,1):
 for backend in ('Rust','C'):
  with np.load(ROOT/f'harmonic-{backend}-nac{nac}.npz') as z:
   f=z['frequency'];a=z['address']
  arrays[backend,nac]=(f,a)
  seen={};pairs=[]
  for i,addr in enumerate(a):
   key=tuple(addr%3)
   if key in seen:
    j=seen[key];difference=np.abs(f[i]-f[j]);b=int(np.argmax(difference))
    pairs.append({'indices':[j,i],'addresses':[a[j].tolist(),a[i].tolist()],
       'band_zero_based':b,'frequencies_THz':[float(f[j,b]),float(f[i,b])],'difference_THz':float(difference[b])})
   else:seen[key]=i
  worst=max(pairs,key=lambda p:p['difference_THz']);minimum=np.unravel_index(np.argmin(f),f.shape)
  out['controls'].append({'backend':backend,'NAC':bool(nac),'worst_duplicate':worst,
       'negative_count':int(np.sum(f<0)),'negative_non_Gamma_count':int(np.sum((f<0)&np.any(a%3!=0,axis=1)[:,None])),
       'minimum':{'index':[int(v) for v in minimum],'address':a[minimum[0]].tolist(),'frequency_THz':float(f[minimum])}})
 fr,ar=arrays['Rust',nac];fc,ac=arrays['C',nac]
 same=np.array_equal(ar,ac);d=np.abs(fr-fc);notgamma=np.any(ar%3!=0,axis=1)
 out['backend_comparisons'].append({'NAC':bool(nac),'address_arrays_identical':bool(same),
     'max_abs_frequency_difference_THz':float(d.max()),'max_abs_non_Gamma_difference_THz':float(d[notgamma].max()),
     'max_abs_signed_square_difference_THz2':float(np.max(np.abs(fr*np.abs(fr)-fc*np.abs(fc))))})
(ROOT/'review-numerical-data.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
