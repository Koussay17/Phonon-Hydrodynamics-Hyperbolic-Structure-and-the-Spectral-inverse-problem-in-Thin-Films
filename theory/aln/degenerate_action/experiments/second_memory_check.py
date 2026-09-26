"""Small-matrix retained resolvent and exact memory check.
Writes second_memory_check.json only. No first-pass or repository edits.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
sys.dont_write_bytecode=True
import numpy as np
import scipy
import scipy.linalg as la
from scipy.integrate import quad_vec
import mpmath as mp

OUT=Path(__file__).resolve().parent
REPO=Path('C:/Users/Koussay/these/src/collision_events.py')
EXPECTED_HASHES={'77196934d0e3ee209e7121f31655af48e074864d7fff5668cfb6f07b1ea99e6c',
                 '1848dc1b248174bdabe14fc5cef10d3073be73a346f5c068753f5623eaec884f'}  # CRLF and canonical git LF
W=np.diag([6/25,6/25,2,3/4])
SW=np.diag(np.sqrt(np.diag(W)))
B=np.array([[1.,1,0,0],[0,0,1,0],[0,0,0,1]])
H=np.array([[.5,0,0],[.5,0,0],[0,1,0],[0,0,1]])
E=np.array([[1/np.sqrt(2),0,0],[1/np.sqrt(2),0,0],[0,1,0],[0,0,1]])
q=np.array([1.,-1,0,0])/np.sqrt(2)
nu=np.array([[-1.,0,1,1],[0,-1,1,1]])
c=np.array([-5/np.sqrt(12),1/np.sqrt(2),2/np.sqrt(3)])
r=c/la.norm(c)
P0=np.eye(3)-np.outer(r,r)
Sigma=B@W@B.T
Sigma_sqrt=np.diag(np.sqrt(np.diag(Sigma)))
rates={'equal':[2.,2.],'unequal':[3.,1.]}

def norm(a):
    return float(la.norm(a,2))

def pair(z):
    return [float(np.real(z)),float(np.imag(z))]

def make(gamma):
    M=sum((.6*gamma[i]*np.outer(nu[i],nu[i]) for i in range(2)),np.zeros((4,4)))
    L=M@np.linalg.inv(W)
    C=np.linalg.solve(SW,M)@np.linalg.inv(SW)
    A=E.T@C@E
    K=q@C@E
    d=float(q@C@q)
    return M,L,C,A,K,d

def resolvents():
    rows=[]
    spectra=[]
    samples=[complex(s,w) for s in [.001,.1,1.] for w in [0.,.4,3.]]+[complex(1e-6)]
    for name,gamma in rates.items():
        M,L,C,A,K,d=make(gamma)
        nullvectors=SW@np.array([[1.,1],[1,1],[1,0],[0,1]])
        spectra.append(dict(model=name,full_eigenvalues=la.eigvalsh(C).tolist(),
            retained_eigenvalues=la.eigvalsh(A).tolist(),D=d,K=K.tolist(),
            K_squared_norm=float(K@K),zero_mode_residual=norm(C@nullvectors)/norm(nullvectors),
            predicted_two_zero_modes=True,
            K_outer_relation_error=norm(np.outer(K,K)-(0 if name=='equal' else 1.25)*A)))
        for z in samples:
            full=z*np.eye(4)+C
            inverse=la.solve(full,np.eye(4))
            retained=E.T@inverse@E
            schur=z*np.eye(3)+A-np.outer(K,K)/(z+d)
            exact=la.solve(schur,np.eye(3))
            local=la.solve(z*np.eye(3)+A,np.eye(3))
            hidden_direct=E.T@inverse@q
            hidden_schur=-la.solve(schur,K)/(z+d)
            pop_retained=B@la.solve(z*np.eye(4)+L,H)
            transformed=Sigma_sqrt@retained@np.linalg.inv(Sigma_sqrt)
            active=r@retained@r
            active_local=1/(z+47/5)
            rows.append(dict(model=name,z=pair(z),condition_full=float(np.linalg.cond(full)),
                schur_relative_error=norm(retained-exact)/norm(retained),
                normalized_full_residual=norm(full@inverse-np.eye(4))/(norm(full)*norm(inverse)+1),
                hidden_source_absolute_error=norm(hidden_direct-hidden_schur),
                hidden_source_norm=norm(hidden_direct),
                population_coordinate_relative_error=norm(pop_retained-transformed)/norm(pop_retained),
                conserved_resolvent_residue_error=norm(z*retained@P0-P0),
                projected_local_full_relative_error=norm(retained-local)/norm(retained),
                active_response=pair(active),
                active_local_relative_error=float(abs(active-active_local)/abs(active))))
    return dict(spectra=spectra,rows=rows)

def time_memory():
    rows=[]
    for name,gamma in rates.items():
        _,_,C,A,K,d=make(gamma)
        initial=E@r+q  # Both retained and hidden initial components.
        for t in [.05,.2,1.]:
            yt=la.expm(-C*t)@initial
            a=E.T@yt
            derivative=-E.T@C@yt
            def integrand(s):
                past=E.T@(la.expm(-C*s)@initial)
                return K*np.exp(-d*(t-s))*(K@past)
            integral,err=quad_vec(integrand,0,t,epsabs=1e-12,epsrel=1e-12)
            source=K*np.exp(-d*t)
            residual=derivative+A@a-integral+source
            scale=norm(derivative)+norm(A@a)+norm(integral)+norm(source)+1
            rows.append(dict(model=name,t=t,quadrature_error_estimate=float(err),
                normalized_memory_equation_residual=norm(residual)/scale,
                hidden_initial_term_norm=norm(source),
                residual_if_hidden_term_omitted=norm(derivative+A@a-integral)))
    return rows

def high_precision():
    mp.mp.dps=80
    sw=mp.diag([mp.sqrt(mp.mpf(6)/25),mp.sqrt(mp.mpf(6)/25),mp.sqrt(2),mp.sqrt(mp.mpf(3)/4)])
    em=mp.matrix([[1/mp.sqrt(2),0,0],[1/mp.sqrt(2),0,0],[0,1,0],[0,0,1]])
    qm=mp.matrix([1,-1,0,0])/mp.sqrt(2)
    cm=mp.matrix([-5/mp.sqrt(12),1/mp.sqrt(2),2/mp.sqrt(3)])
    rm=cm/mp.norm(cm)
    pzero=mp.eye(3)-rm*rm.T
    vv=[mp.matrix([-1,0,1,1]),mp.matrix([0,-1,1,1])]
    rows=[]
    for name,gamma in rates.items():
        C=mp.zeros(4)
        for i in range(2):
            dd=sw**-1*vv[i]
            C+=mp.mpf(3)/5*int(gamma[i])*dd*dd.T
        A=em.T*C*em
        K=qm.T*C*em
        D=(qm.T*C*qm)[0]
        for zs in ['1e-3','1e-6','1e-10']:
            z=mp.mpf(zs)
            full=mp.matrix(4,3)
            for col in range(3):
                solved=mp.lu_solve(z*mp.eye(4)+C,em[:,col])
                for i in range(4):
                    full[i,col]=solved[i]
            retained=em.T*full
            exact=(z*mp.eye(3)+A-K.T*K/(z+D))**-1
            active=(rm.T*retained*rm)[0]
            analytic=1/(z+mp.mpf(47)/5) if name=='equal' else (z+5)/(z*z+mp.mpf(72)/5*z+mp.mpf(141)/4)
            local=1/(z+mp.mpf(47)/5)
            rows.append(dict(model=name,z=zs,
                dense_schur_relative_frobenius_error=float(mp.norm(retained-exact)/mp.norm(retained)),
                active_scalar_relative_error=float(abs(active-analytic)/abs(analytic)),
                conserved_residue_error=float(mp.norm(z*retained*pzero-pzero)),
                active_local_relative_error=float(abs(active-local)/abs(active)),
                full_local_spectral_relative_error=float(abs(analytic-local)/max(1/z,abs(analytic)))))
    return dict(decimal_digits=80,rows=rows)

def repo_comparison():
    sha=hashlib.sha256(REPO.read_bytes()).hexdigest()
    if sha not in EXPECTED_HASHES:
        return dict(status='skipped: source changed since inspected version',sha256=sha)
    spec=importlib.util.spec_from_file_location('_memory_repo_events',REPO)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    x=np.log(np.array([6.,6.,2.,3.]))
    events=np.array([[0,2,3],[1,2,3]])
    rows=[]
    mp.mp.dps=80
    xm=[mp.mpf(float(v)) for v in x]
    pops=[1/mp.expm1(v) for v in xm]
    sm=mp.diag([mp.sqrt(p*(1+p)) for p in pops])
    for name,gamma in rates.items():
        op=module.DecayEventOperator(x,events,np.array(gamma))
        got=np.column_stack([op.action(np.eye(4)[:,i]) for i in range(4)])
        _,_,exact,_,_,_=make(gamma)
        jac=mp.zeros(4)
        rhs=mp.zeros(4,1)
        for i,g in enumerate(gamma):
            vec=mp.matrix([-1,0,1,1] if i==0 else [0,-1,1,1])
            grad=mp.zeros(4,1)
            grad[i]=1+pops[2]+pops[3]
            grad[2]=pops[i]-pops[3]
            grad[3]=pops[i]-pops[2]
            jac+=int(g)*vec*grad.T
            flux=pops[i]*(1+pops[2])*(1+pops[3])-(1+pops[i])*pops[2]*pops[3]
            rhs+=int(g)*vec*flux
        actual=-sm**-1*jac*sm
        rows.append(dict(model=name,operator_relative_error=norm(got-exact)/norm(exact),
            stored_detuning=op.detailed_balance_log_residual.tolist(),
            exact_binary_detuning=mp.nstr(xm[0]-xm[2]-xm[3],30),
            true_binary_energy_jacobian_asymmetry_relative_frobenius=float(mp.norm(actual-actual.T)/mp.norm(actual)),
            true_binary_energy_equilibrium_rhs_norm=float(mp.norm(rhs))))
    return dict(status='compared',source=str(REPO),sha256_before=sha,
        sha256_after=hashlib.sha256(REPO.read_bytes()).hexdigest(),energy_values=x.tolist(),
        default_resonance_tolerance=1e-12,rows=rows)

def main():
    result=dict(metadata=dict(python=platform.python_version(),python_executable=sys.executable,
        repo_root=str(REPO.parent.parent),output_dir=str(OUT),numpy=np.__version__,
        scipy=scipy.__version__,mpmath=mp.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='Four-mode scalar reversible model; no material or quantum-generator validation'),
        resolvents=resolvents(),time_memory=time_memory(),high_precision=high_precision(),
        repo=repo_comparison())
    rr=result['resolvents']['rows']
    hp=result['high_precision']['rows']
    checks=dict(
        double_schur_identity=max(r['schur_relative_error'] for r in rr)<1e-7,
        population_coordinate_change=max(r['population_coordinate_relative_error'] for r in rr)<1e-7,
        high_precision_identity=max(r['dense_schur_relative_frobenius_error'] for r in hp)<1e-60,
        high_precision_active_formula=max(r['active_scalar_relative_error'] for r in hp)<1e-60,
        memory_equation=max(r['normalized_memory_equation_residual'] for r in result['time_memory'])<1e-12,
        equal_rate_positive_control=max(r['active_local_relative_error'] for r in hp if r['model']=='equal')<1e-60,
        unequal_relaxing_response_error=hp[-1]['active_local_relative_error']>.2499,
        zero_modes_mask_error=hp[-1]['full_local_spectral_relative_error']<1e-10,
        hidden_initial_term_required=max(r['residual_if_hidden_term_omitted'] for r in result['time_memory'] if r['model']=='unequal')>1)
    if result['repo']['status']=='compared':
        checks['repository_action_matches']=max(r['operator_relative_error'] for r in result['repo']['rows'])<1e-13
        checks['repository_source_unchanged']=result['repo']['sha256_before']==result['repo']['sha256_after']
    result['checks']=checks
    (OUT/'second_memory_check.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    summary=dict(result)
    summary['resolvents']=dict(spectra=result['resolvents']['spectra'],
        maxima={key:max(r[key] for r in rr) for key in ['schur_relative_error','normalized_full_residual',
        'hidden_source_absolute_error','population_coordinate_relative_error','condition_full']})
    print(json.dumps(summary,indent=2))
    assert all(checks.values()),checks

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=REPO.parent.parent,
        help='Repository root containing src/collision_events.py; the inspected source hash is pinned.')
    parser.add_argument('--output-dir',type=Path,default=OUT,
        help='Directory for second_memory_check.json; defaults to the script directory.')
    args=parser.parse_args()
    REPO=args.repo_root.expanduser().resolve()/'src'/'collision_events.py'
    OUT=args.output_dir.expanduser().resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    main()

