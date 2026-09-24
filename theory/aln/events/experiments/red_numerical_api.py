"""Small independent tests of the frozen event reference; no repository edits."""
from pathlib import Path
import importlib.util,json,hashlib,platform,warnings
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
FILE=ROOT/'review-inputs/collision_events.py'
spec=importlib.util.spec_from_file_location('frozen_collision_events',FILE)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
Operator=module.DecayEventOperator

def dense(op):return np.column_stack([op.action(v) for v in np.eye(len(op.energy_over_kbt))])
def outcome(fun):
    try:
        result=fun()
        return {'accepted':True,'result':result}
    except Exception as exc:return {'accepted':False,'exception':type(exc).__name__,'message':str(exc)}

def event_outcome(x,event,tol=1e-12):
    op=Operator(x,[event],[1.],resonance_rtol=tol)
    return {'relative_residual':op.relative_energy_residual.tolist(),'reported_detuning':op.detailed_balance_log_residual.tolist(),
            'factor':op.factor.toarray().tolist(),'collision':dense(op).tolist()}

out={'environment':{'python':platform.python_version(),'numpy':np.__version__,'mpmath':mp.__version__},
     'frozen_source_sha256':hashlib.sha256(FILE.read_bytes()).hexdigest(),'resonance_checks':[]}
for x,events,tol in [([1e16,1.,1e16],(0,1,2),1e-12),([1e16,1.,1e16],(0,2,1),1e-12),
                      ([3.,1.,2.],(0,1,2),0.),([3.,1.,2.],(0,2,1),0.),
                      ([1e5,5e-12,1e5],(0,1,2),1e-12)]:
    with mp.workdps(100):
        xx=[mp.mpf(z) for z in x];detuning=xx[events[0]]-xx[events[1]]-xx[events[2]]
    out['resonance_checks'].append({'x':x,'event':events,'tolerance':tol,'exact_input_detuning':mp.nstr(detuning,30),
                                  **outcome(lambda:event_outcome(x,events,tol))})

# High-precision derivative from the population polynomial, using the actual
# stored binary64 energy inputs. No equilibrium outer-product identity used.
with mp.workdps(100):
    x=[mp.mpf(1e16),mp.mpf(1),mp.mpf(1e16)];n=[1/mp.expm1(v) for v in x];d=[mp.sqrt(v*(1+v)) for v in n]
    s=[-1,1,1];g=[1+n[1]+n[2],n[0]-n[2],n[0]-n[1]]
    truth=mp.matrix([[-s[i]*g[j]*d[j]/d[i] for j in range(3)] for i in range(3)])
    actual=mp.matrix(dense(Operator([1e16,1.,1e16],[[0,1,2]],[1.])).tolist())
    out['accepted_detuned_derivative']={'true_collision':[[mp.nstr(z,18) for z in row] for row in truth.tolist()],
          'matrix_relative_error':mp.nstr(mp.norm(actual-truth)/mp.norm(truth),18),
          'forward_reverse_flux_ratio':mp.nstr(mp.e,18),'equilibrium_flux_over_reverse':mp.nstr(mp.e-1,18)}

op=Operator([3.,1.,2.],[[0,1,2]],[1.]);linear=op.as_linear_operator()
out['sparse_api']={}
for name,fun in [('nan_vector',lambda:linear@np.array([np.nan,0.,0.])),
                 ('nan_two_columns',lambda:linear@np.array([[np.nan,1.],[0.,0.],[0.,0.]]))]:
    def check(fun=fun):
        z=np.asarray(fun());return {'shape':list(z.shape),'all_finite':bool(np.isfinite(z).all())}
    out['sparse_api'][name]=outcome(check)

# Test tolerance effect without attributing numerical equality to bit equality.
c=dense(op)
try:
    np.testing.assert_allclose(1.00000005*c,c,atol=0)
    out['allclose_atol_zero_allows_relative_change_5e_8']=True
except AssertionError:out['allclose_atol_zero_allows_relative_change_5e_8']=False

# Retain the stated moderate-input limitation: no attempt to improve the RHS.
small=Operator([3e-18,1e-18,2e-18],[[0,1,2]],[1.]);pop=1/np.expm1(small.energy_over_kbt);scale=small.entropy_scale
h=1e-4;direction=np.array([1.,.3,-.4])*pop/scale
rhsdiff=(small.population_rhs(pop+h*scale*direction)-small.population_rhs(pop-h*scale*direction))/(2*h*scale)
out['rhs_small_energy']={'finite_difference_relative_error':float(np.linalg.norm(rhsdiff+small.action(direction))/np.linalg.norm(small.action(direction))),
    'computed_derivative_all_zero':bool(np.all(rhsdiff==0)),'base_state_rhs_norm':float(np.linalg.norm(small.population_rhs(pop))),
    'scope':'Documented extreme-occupation limitation, not claimed as an in-contract moderate-input defect.'}

path=ROOT/'experiments/red_numerical_api_results.json'
path.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2,allow_nan=False))
