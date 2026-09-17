"""Bulk diagonal-RTA harmonic conductivity; not a simulated FDTR signal."""
import argparse,json,sys
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
a=parser.parse_args()
sys.path.insert(0,str(a.repo))
from src.spectral import mode_heat_capacity,rta_moments
with np.load(a.repo/'theory/aln/rao_300K_modes.npz') as z:
    c=mode_heat_capacity(z['omega_rad_s'],300)*z['degeneracy'][:,None]/z['degeneracy'].sum()/z['cell_volume_m3']
    v=z['velocity_m_s'];r=z['total_rate_s']
    tau=np.divide(1.,r,out=np.zeros_like(r),where=r>0)
frequencies=[1e3,1e5,1e6,1e7,5e7,1e8,2e8,1e9]
result={'scope':'uniform bulk diagonal RTA conductivity at 300 K; no spatial streaming, no FDTR boundary response, no N/U inference','harmonic_convention':'exp(+i omega t)', 'directions':{}}
for name,v2 in [('basal',(v[:,:,0]**2+v[:,:,1]**2)/2),('c_axis',v[:,:,2]**2)]:
    m=rta_moments(c,v2,r); contributions=m.pop('contributions');rows=[]
    for f in frequencies:
        s=2j*np.pi*f
        exact=np.sum(contributions/(1+s*tau))
        static=m['kappa']/(1+s*m['tau_static_s'])
        memory=m['kappa']/(1+s*m['tau_memory_s'])
        rows.append({'frequency_Hz':f,'kappa_real_W_mK':float(exact.real),'kappa_imag_W_mK':float(exact.imag),'phase_degree':float(np.angle(exact,deg=True)),
                     'static_pole_relative_complex_error':float(abs(static-exact)/abs(exact)),
                     'memory_pole_relative_complex_error':float(abs(memory-exact)/abs(exact))})
    result['directions'][name]={'moments':m,'response':rows}
(a.repo/'theory/aln/bulk_response.json').write_text(json.dumps(result,indent=2)+'\n')
for name,d in result['directions'].items():print(name,d['response'][-2])
