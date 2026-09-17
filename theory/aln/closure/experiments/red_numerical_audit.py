"""Independent red-team checks; writes red_numerical_* only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,platform,sys
import mpmath as mp
import numpy as np
OUT=Path(__file__).resolve().parent
DATA=Path(__file__).resolve().parents[2]/'rao_300K_modes.npz'
S=[[1 if (i&j).bit_count()%2==0 else -1 for j in range(8)] for i in range(8)]
SF=np.array(S,dtype=float)

def exact_matrix(t):
    d=[0,t*t,2,2-t*t,1,1,1,1]
    return [[(sum(d[k]*S[i][k]*S[j][k] for k in range(8))+t*(S[i][1]*S[j][2]+S[i][2]*S[j][1]-S[i][4]*S[j][7]-S[i][7]*S[j][4]))/8 for j in range(8)] for i in range(8)]

def rational_solve(a,b):
    a=[list(row)+[b[i]] for i,row in enumerate(a)];n=len(b)
    for k in range(n):
        p=max(range(k,n),key=lambda i:abs(a[i][k]));a[k],a[p]=a[p],a[k]
        for i in range(k+1,n):
            v=a[i][k]/a[k][k]
            for j in range(k,n+1):a[i][j]-=v*a[k][j]
    x=[F(0)]*n
    for i in reversed(range(n)):x[i]=(a[i][n]-sum(a[i][j]*x[j] for j in range(i+1,n)))/a[i][i]
    return x

def algebra_and_conditioning():
    result={'exact':[],'double':[],'precision':[],'perturbation':[]}
    e=np.ones(8)/np.sqrt(8);b=SF[:,2]/np.sqrt(8);h=SF/np.sqrt(8)
    for ts in ('0.5','0.001','1e-8','1e-20'):
        t=F(ts);c=exact_matrix(t)
        x=rational_solve([[v+F(1,8) for v in row] for row in c],[F(row[2]) for row in S])
        dc=sum(S[i][2]*x[i] for i in range(8))/8;m2=sum(v*v for v in x)/8
        result['exact'].append({'t':ts,'dc':str(dc),'memory':str(m2/dc),
          'zero_residual':all(sum(c[i][j]*x[j] for j in range(8))==S[i][2] for i in range(8)),
          'fixed_diagonal':all(c[i][i]==1 for i in range(8)),'energy_conservation':all(sum(row)==0 for row in c)})
    for ts in ('0.001','1e-5','1e-6','1e-7','1e-8','1e-9','0.9999999999999999'):
        t=float(ts);a=np.diag([0,t*t,2,2-t*t,1,1,1,1]);a[1,2]=a[2,1]=t;a[4,7]=a[7,4]=-t
        target=SF@a@SF.T/8;nearest=np.array([[float(v) for v in row] for row in exact_matrix(F(ts))])
        with mp.workdps(80):
            tt=mp.mpf(ts);lp=(2+tt*tt+mp.sqrt(4+tt**4))/2;lm=tt*tt/lp
            condition=float(lp/min(lm,1-tt))
        for name,c in [('target_transform',target),('exact_entries_rounded_once',nearest)]:
            row={'t':ts,'assembly':name,'exact_condition':condition,'assembly_difference_norm':float(np.linalg.norm(c-nearest,2))}
            lifted=c+np.outer(e,e)
            try:
                x=np.linalg.solve(lifted,b);dc=float(b@x);tm=float(x@x/dc)
                row.update(dc=dc,memory=tm,dc_error=abs(dc-1),memory_relative_error=abs(tm/(1+1/t**2)-1),
                  unlifted_residual=float(np.linalg.norm(c@x-b)),
                  normwise_backward_error=float(np.linalg.norm(lifted@x-b)/(np.linalg.norm(lifted,2)*np.linalg.norm(x)+np.linalg.norm(b))))
            except np.linalg.LinAlgError as exc:row['solve_error']=str(exc)
            values,vectors=np.linalg.eigh(h[:,1:].T@c@h[:,1:]);weights=(vectors.T@(h[:,1:].T@b))**2
            dc=float(np.sum(weights/values));tm=float(np.sum(weights/values**2)/dc)
            row['spectral']={'least_eigenvalue':float(values.min()),'dc':dc,'memory':tm,'memory_relative_error':abs(tm/(1+1/t**2)-1)}
            p=.3;kp=(p+t*t)/(p*p+(2+t*t)*p+t*t)
            row['p_point3_resolvent_error']=abs(float(b@np.linalg.solve(c+p*np.eye(8),b))-kp)
            result['double'].append(row)
    for ts in ('1e-20','1e-40'):
        cc=exact_matrix(F(ts))
        for dps in (50,80,110):
            with mp.workdps(dps):
                t=mp.mpf(ts);c=mp.matrix([[mp.mpf(v.numerator)/v.denominator for v in row] for row in cc]);bb=mp.matrix([row[2] for row in S])
                row={'t':ts,'decimal_precision':dps}
                try:
                    x=mp.lu_solve(c+mp.ones(8,8)/8,bb);dc=(bb.T*x)[0]/8;tm=(x.T*x)[0]/(8*dc)
                    row.update(dc=mp.nstr(dc,25),memory=mp.nstr(tm,25),dc_error=mp.nstr(abs(dc-1),12),
                      memory_relative_error=mp.nstr(abs(tm/(1+1/t**2)-1),12),residual=mp.nstr(mp.norm(c*x-bb)/mp.norm(bb),12))
                except Exception as exc:row['error']=str(exc)
                result['precision'].append(row)
    for beta in (F(1),F(-49,100)):
        t=F('1e-8');delta=beta*t*t;dc=(t*t+delta)/(t*t+2*delta)
        tm=(t*t+(t*t+delta)**2)/((t*t+2*delta)*(t*t+delta))
        result['perturbation'].append({'t':str(t),'perturbation_spectral_norm':str(abs(delta)),'dc':str(dc),
          'memory_ratio':str(tm/(1+1/t**2)),'psd':t*t+2*delta>0,
          'description':'Add delta to Walsh A11 and subtract delta from A33: diagonal, energy, parity and b preserved; DC changes.'})
    return result

def bulk_data():
    mp.mp.dps=75
    kb=mp.mpf('1.380649e-23');hbar=mp.mpf('6.62607015e-34')/(2*mp.pi)
    with np.load(DATA,allow_pickle=False) as z:
        w=z['omega_rad_s'];v=z['velocity_m_s'];r=z['total_rate_s'];deg=z['degeneracy'];vol=mp.mpf(float(z['cell_volume_m3']))
    data={'basal':[],'c_axis':[]};ng=int(sum(deg))
    for i in range(w.shape[0]):
        for j in range(w.shape[1]):
            ww=mp.mpf(float(w[i,j]));rate=mp.mpf(float(r[i,j]));vv=[mp.mpf(float(u)) for u in v[i,j]]
            xx=hbar*ww/(kb*300);heat=kb if xx==0 else kb*(xx/(2*mp.sinh(xx/2)))**2
            cap=heat*int(deg[i])/ng/vol
            if rate<=0:
                assert all(u==0 for u in vv)
                continue
            tau=1/rate;freq=float(w[i,j]/(2*np.pi*1e12));vz=float(abs(v[i,j,2]))
            for key,v2 in [('basal',(vv[0]**2+vv[1]**2)/2),('c_axis',vv[2]**2)]:data[key].append((cap*v2,tau,freq,vz))
    published=json.loads((OUT/'pi_bulk_response.json').read_text())
    result={'shape':list(w.shape),'degeneracy_sum':ng,'mp_dps':75,'constants':'Exact SI h and kB; hbar=h/(2pi)','directions':{}}
    for key,rows in data.items():
        m=[mp.fsum(a*t**n for a,t,_,_ in rows) for n in range(4)];tm=m[2]/m[1];ts=m[1]/m[0]
        pub=published['directions'][key]
        out={'moments':[mp.nstr(u,25) for u in m],'memory_ps':float(tm*1e12),'static_ps':float(ts*1e12),
          'published_dc_relative_error':float(abs(mp.mpf(pub['moments']['kappa'])/m[1]-1)),
          'published_memory_relative_error':float(abs(mp.mpf(pub['moments']['tau_memory_s'])/tm-1)),'response':[],'infrared':{}}
        for cutoff in (1,2):
            low=[row for row in rows if 0<row[2]<=cutoff]
            for n in (1,2):out['infrared'][f'M{n}_below_{cutoff}_THz_fraction']=float(mp.fsum(a*t**n for a,t,_,_ in low)/m[n])
        if key=='basal':
            zero=[row for row in rows if row[3]==0]
            out['grazing_zero_dc_fraction']=float(mp.fsum(a*t for a,t,_,_ in zero)/m[1])
            out['grazing_zero_M2_fraction']=float(mp.fsum(a*t*t for a,t,_,_ in zero)/m[2])
            out['grazing_zero_floor']=float(mp.fsum(a*t for a,t,_,_ in zero))
        for pubrow in pub['response']:
            f=pubrow['frequency_Hz'];s=2j*mp.pi*f
            exact=mp.fsum(a/(1/t+s) for a,t,_,_ in rows);pole=m[1]/(1+s*tm);static=m[1]/(1+s*ts)
            err=abs(pole-exact)/abs(exact)
            val=mp.mpc(pubrow['kappa_real_W_mK'],pubrow['kappa_imag_W_mK'])
            t64=np.array([float(t) for _,t,_,_ in rows]);wt=np.array([float(a*t) for a,t,_,_ in rows]);kap=math.fsum(wt.tolist())
            mean=math.fsum((wt*t64).tolist())/kap;ss=2j*np.pi*f
            gap=ss**2/(1+ss*mean)**2*np.sum(wt*(t64-mean)**2/(1+ss*t64));stable=abs(gap)/float(abs(exact))
            out['response'].append({'f_Hz':f,'real':float(exact.real),'imag':float(exact.imag),'response_relative_error':float(abs(val-exact)/abs(exact)),
              'memory_pole_relative_error':float(err),'reported_error_relative_discrepancy':float(abs(mp.mpf(pubrow['memory_pole_relative_complex_error'])/err-1)),
              'stable_identity_relative_discrepancy':float(abs(mp.mpf(stable)/err-1)),
              'memory_pole_dynamic_correction_error':float(abs(pole-exact)/abs(exact-m[1])),
              'static_pole_dynamic_correction_error':float(abs(static-exact)/abs(exact-m[1]))})
        result['directions'][key]=out
    return result

if __name__=='__main__':
    result={'environment':{'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'mpmath':mp.__version__},
      'input_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'D':algebra_and_conditioning()}
    out=OUT/'red_numerical_results.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print('Saved D checks',flush=True)
    result['bulk']=bulk_data();out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(str(out),flush=True)
