"""Independent Bose-event Jacobian/limit checks; no prototype imports.

Gamma=1 defines an illustrative event scale, never an AlN material rate.
Only files named D_* in this campaign's experiments folder are written.
"""
from pathlib import Path
import json,math,platform,sys,hashlib
import mpmath as mp
import numpy as np
OUT=Path(__file__).resolve().parent

def incidence(n,event):
    s=[0]*n
    s[event[0]]-=1;s[event[1]]+=1;s[event[2]]+=1
    return mp.matrix(s)

def flux(n,event,form='products'):
    p,a,b=event
    if form=='products':return n[p]*(1+n[a])*(1+n[b])-(1+n[p])*n[a]*n[b]
    return n[p]+n[p]*n[a]+n[p]*n[b]-n[a]*n[b]

def quantities(text,event):
    x=mp.matrix([mp.mpf(t) for t in text]);n=mp.matrix([1/mp.expm1(z) for z in x]);r=mp.matrix([mp.sqrt(z*(1+z)) for z in n])
    s=incidence(len(x),event);p,a,b=event
    g=mp.matrix(len(x),1);g[p]+=1+n[a]+n[b];g[a]+=n[p]-n[b];g[b]+=n[p]-n[a]
    v=mp.matrix([s[i]/r[i] for i in range(len(x))]);w=mp.matrix([g[i]*r[i] for i in range(len(x))])
    c=-v*w.T
    q=n[p]*(1+n[a])*(1+n[b]);outer=q*v*v.T
    energy=mp.matrix([x[i]*r[i] for i in range(len(x))])
    return x,n,r,s,g,c,outer,energy,q

def jacobian_fd(n,r,s,event,eta,form='products'):
    size=len(n);c=mp.matrix(size)
    for j in range(size):
        step=eta*n[j]/r[j];npop=n.copy();nneg=n.copy();npop[j]+=r[j]*step;nneg[j]-=r[j]*step
        df=(flux(npop,event,form)-flux(nneg,event,form))/(2*step)
        for i in range(size):c[i,j]=-s[i]*df/r[i]
    return c

def rel(a,b):return mp.norm(a-b)/mp.norm(b)
def dec(x):return mp.nstr(x,14)

def precision_tests():
    out=[]
    cases=[(['1.8','.7','1.1'],(0,1,2)),(['3e-30','1e-30','2e-30'],(0,1,2)),(['1000','400','600'],(0,1,2)),(['1.4','.7'],(0,1,1))]
    for xs,event in cases:
        for dps in (50,80,120):
            with mp.workdps(dps):
                x,n,r,s,g,c,outer,e,q=quantities(xs,event)
                fd=jacobian_fd(n,r,s,event,mp.mpf('1e-6'))
                row={'x':xs,'event':event,'dps':dps,'fd_norm_relative_error':dec(rel(fd,c)),
                     'fd_max_entry_relative_error':dec(max(abs((fd[i,j]-c[i,j])/c[i,j]) for i in range(len(n)) for j in range(len(n)) if c[i,j])),
                     'analytic_vs_resonant_outer_error':dec(rel(c,outer)),
                     'normalized_energy_residual':dec(mp.norm(c.T*e)/(mp.norm(c)*mp.norm(e))),
                     'symmetry_defect':dec(mp.norm(c-c.T)/mp.norm(c)),
                     'nonzero_eigenvalue_trace':dec(mp.trace(c) if hasattr(mp,'trace') else sum(c[i,i] for i in range(len(n)))),
                     'minimum_physical_entropy_step':dec(min(n[i]/r[i] for i in range(len(n))))}
                out.append(row)
    return out

def directional_convergence():
    out=[]
    with mp.workdps(90):
        for xs,event,coeff in [(['1.8','.7','1.1'],(0,1,2),[1,mp.mpf('.3'),mp.mpf('-.4')]),(['1.4','.7'],(0,1,1),[1,mp.mpf('.3')])]:
            x,n,r,s,g,c,outer,e,q=quantities(xs,event)
            u=mp.matrix([n[i]/r[i]*coeff[i] for i in range(len(n))]);expected=c*u;previous=None
            for exponent in range(1,13):
                eta=mp.power(10,-exponent);plus=n+eta*mp.matrix([r[i]*u[i] for i in range(len(n))]);minus=2*n-plus
                fp=(flux(plus,event)-flux(n,event))/eta
                fc=(flux(plus,event)-flux(minus,event))/(2*eta)
                actionp=mp.matrix([-s[i]*fp/r[i] for i in range(len(n))]);actionc=mp.matrix([-s[i]*fc/r[i] for i in range(len(n))])
                er=rel(actionp,expected)
                row={'x':xs,'eta':dec(eta),'forward_error':dec(er),'central_error':dec(rel(actionc,expected)),
                     'min_relative_population':dec(min(plus[i]/n[i] for i in range(len(n))))}
                if previous is not None:row['observed_forward_order']=dec(mp.log(previous/er)/mp.log(10))
                out.append(row);previous=er
    return out

def float_small_limits():
    out=[]
    for ep in ('1','1e-2','1e-6','1e-10','1e-14','1e-18','1e-30'):
        with mp.workdps(100):
            eps=mp.mpf(ep);xs=[str(3*eps),str(eps),str(2*eps)];x,n,r,s,g,c,outer,e,q=quantities(xs,(0,1,2))
            ref=np.array(c.tolist(),float)
        xx=np.array([3.,1.,2.])*float(ep);nn=1/np.expm1(xx);rr=np.sqrt(nn*(1+nn));ss=np.array([-1.,1.,1.])
        for form in ('products','quadratic'):
            def ff(pop):
                p,a,b=pop
                return p*(1+a)*(1+b)-(1+p)*a*b if form=='products' else p+p*a+p*b-a*b
            base=ff(nn);rhs=-ss*base/rr
            for eta in (1e-2,1e-4,1e-6,1e-8):
                cc=np.zeros((3,3))
                for j in range(3):
                    step=eta*nn[j]/rr[j];plus=nn.copy();minus=nn.copy();plus[j]+=step*rr[j];minus[j]-=step*rr[j]
                    cc[:,j]=-ss/rr*(ff(plus)-ff(minus))/(2*step)
                out.append({'epsilon':ep,'form':form,'eta':eta,'jacobian_relative_error':float(np.linalg.norm(cc-ref)/np.linalg.norm(ref)),
                     'equilibrium_entropy_drift_over_operator_norm':float(np.linalg.norm(rhs)/np.linalg.norm(ref)),
                     'all_fd_entries_zero':bool(np.all(cc==0))})
    return out

def limits():
    out={'small':[],'large':[],'zero_energy':[]}
    with mp.workdps(120):
        for ep in ('1e-2','1e-6','1e-12','1e-30'):
            eps=mp.mpf(ep)
            for rep in (False,True):
                xs=[str(2*eps),str(eps)] if rep else [str(3*eps),str(eps),str(2*eps)]
                event=(0,1,1) if rep else (0,1,2)
                *_,c,outer,e,q=quantities(xs,event)
                lam=sum(c[i,i] for i in range(c.rows));target=mp.mpf(4) if rep else mp.mpf(7)/3
                out['small'].append({'epsilon':ep,'repeated':rep,'epsilon_times_rate':dec(eps*lam),'asymptotic_coefficient':dec(target),'relative_asymptotic_error':dec(abs(eps*lam/target-1))})
        for xs in (['20','8','12'],['200','80','120'],['1000','400','600']):
            x,n,r,s,g,c,outer,e,q=quantities(xs,(0,1,2));lam=sum(c[i,i] for i in range(c.rows))
            xx=np.array(xs,float)
            with np.errstate(all='ignore'):
                nn=1/np.expm1(xx);rr=np.sqrt(nn*(1+nn));qq=nn[0]*(1+nn[1])*(1+nn[2]);vv=np.array([-1,1,1])/rr;cc=qq*np.outer(vv,vv)
            out['large'].append({'x':xs,'high_precision_rate':dec(lam),'parent_diagonal':dec(c[0,0]),'min_population':dec(min(n)),
                 'max_common_entropy_step':dec(min(n[i]/r[i] for i in range(len(n)))),
                 'binary64_population':nn.tolist(),'binary64_outer_all_finite':bool(np.isfinite(cc).all()),
                 'fixed_entropy_step_1e-6_preserves_population':all(mp.mpf('1e-6')<n[i]/r[i] for i in range(len(n)))})
        for ep in ('1e-2','1e-6','1e-12','1e-30'):
            eps=mp.mpf(ep);x,n,r,s,g,c,outer,e,q=quantities([str(1+eps),str(eps),'1'],(0,1,2));lam=sum(c[i,i] for i in range(c.rows))
            out['zero_energy'].append({'daughter_x':ep,'x_times_rate':dec(eps*lam),'rate':dec(lam)})
        try:mp.mpf(1)/mp.expm1(0)
        except Exception as exc:out['exact_zero_exception']=type(exc).__name__
    return out

def detuning():
    out=[]
    with mp.workdps(100):
        for ds in ('0','.1','.001','.000001','-.1'):
            delta=mp.mpf(ds);x,n,r,s,g,c,outer,e,q=quantities([str(mp.mpf('1.8')+delta),'.7','1.1'],(0,1,2))
            reverse=(1+n[0])*n[1]*n[2];f=flux(n,(0,1,2));sym=(c+c.T)/2;values=mp.eigsy(sym,eigvals_only=True)
            energydrift=(x.T*s)[0]*f
            out.append({'detuning':ds,'equilibrium_flux_over_reverse':dec(f/reverse),
                'expected_expm1_negative_detuning':dec(mp.expm1(-delta)),
                'relative_symmetry_defect':dec(mp.norm(c-c.T)/mp.norm(c)),
                'minimum_symmetric_eigenvalue':dec(values[0]),'minimum_over_maximum_symmetric_eigenvalue':dec(values[0]/values[values.rows-1]),
                'left_energy_residual':dec(mp.norm(c.T*e)/(mp.norm(c)*mp.norm(e))),
                'equilibrium_energy_drift_over_reverse':dec(energydrift/reverse),
                'resonant_outer_relative_error':dec(rel(outer,c))})
    return out

def repeated_and_reciprocal():
    with mp.workdps(90):
        x,n,r,s,g,c,outer,e,q=quantities(['1.4','.7'],(0,1,1))
        wrong=mp.matrix([-1,1]);wrongouter=q*mp.matrix([wrong[i]/r[i] for i in range(2)])*mp.matrix([wrong[i]/r[i] for i in range(2)]).T
        repeated={'correct_incidence':[int(v) for v in s],'correct_jacobian':[[dec(z) for z in row] for row in c.tolist()],
            'wrong_unit_incidence_energy_residual':dec(mp.norm(wrongouter*e)/(mp.norm(wrongouter)*mp.norm(e))),
            'correct_daughter_diagonal_over_Q_by_d':dec(c[1,1]/(q/r[1]**2)),
            'wrong_unit_incidence_relative_matrix_error':dec(rel(wrongouter,c))}
        x,n,r,s,g,c,outer,e,q=quantities(['1.8','.7','1.1'],(0,1,2))
        pair=mp.matrix(6);parity=mp.matrix(6)
        for i in range(3):
            parity[i,i+3]=parity[i+3,i]=1
            for j in range(3):pair[i,j]=pair[i+3,j+3]=c[i,j]
        ee=mp.matrix(list(e)+list(e));test=mp.matrix([mp.mpf(t) for t in ('1','.3','-.8','.1','-.2','.7')])
        pairout={'relative_parity_commutator':dec(mp.norm(pair*parity-parity*pair)/mp.norm(pair)),
           'relative_energy_residual':dec(mp.norm(pair*ee)/(mp.norm(pair)*mp.norm(ee))),
           'chosen_vector_dissipation':dec((test.T*pair*test)[0]),'rank_exact_from_disjoint_events':2,'kernel_dimension':4,
           'daughter_permutation_flux_difference':dec(flux(n,(0,1,2))-flux(n,(0,2,1))),
           'counting_both_ordered_daughter_permutations_factor':2}
        return {'repeated':repeated,'reciprocal_pair':pairout}

if __name__=='__main__':
    target=OUT/'D_bose_jacobian_results.json'
    results={'environment':{'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'mpmath':mp.__version__},
      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gamma':1,'retained_failure':'Initial run used negative indexing of an mpmath matrix, returning zero in an eigenvalue ratio. Corrected to the explicit last row; no physical inference used that failed ratio.',
      'precision_tests':precision_tests(),'directional_convergence':directional_convergence()}
    target.write_text(json.dumps(results,indent=2,allow_nan=False)+'\n',encoding='utf-8');print('Saved differentiation and precision checks',flush=True)
    results.update(binary64_small= float_small_limits(),asymptotic_limits=limits(),detuning=detuning(),counting_and_reciprocity=repeated_and_reciprocal())
    target.write_text(json.dumps(results,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(str(target),flush=True)

