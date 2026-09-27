"""Independent bounded algebra checks; no material data or candidate code imported."""
import json
from pathlib import Path
import sympy as s
from sympy.functions.combinatorial.numbers import stirling
X,Y,p,a,c,k,A=s.symbols('X Y p a c k A', positive=True)
wf=c*X*(Y+1)*(Y+2)
wr=c*(X+1)*Y*(Y-1)
def gen(f):
    return s.expand(wf*(f.subs({X:X-1,Y:Y+2},simultaneous=True)-f)+wr*(f.subs({X:X+1,Y:Y-2},simultaneous=True)-f))
def raw(n,z):
    return sum(stirling(n,j,kind=2)*s.factorial(j)*z**j for j in range(n+1))
def av(f):
    return s.expand(sum(coef*raw(i,p)*raw(j,a) for (i,j),coef in s.Poly(s.expand(f),X,Y).terms()))
checks={}
def zero(name,expr):
    r=s.simplify(expr)
    assert r==0,(name,r)
    checks[name]=str(r)
F=p*(1+2*a)-a*a
J=2*c*F
zero('mean_flux',av(wf-wr)-J)
zero('parent_mean',av(gen(X))+J)
zero('daughter_mean',av(gen(Y))-2*J)
zero('daughter_factorial_derivative',av(gen(Y*(Y-1)))-2*(6*a+1)*J)
zero('factorial_tangency_defect',av(gen(Y*(Y-1)))-4*a*av(gen(Y))-2*(2*a+1)*J)
zero('covariance_defect',av(gen(X*Y))-a*av(gen(X))-p*av(gen(Y))-2*(p-a)*J)
zero('first_order_factorial_defect',(2*(2*a+1)*J).subs({a:1,p:s.Rational(1,3)+s.Symbol('eps')})-36*c*s.Symbol('eps'))
P=A*A/(1+2*A)
wp=P*(1+P); wa=A*(1+A); Q=P*(1+A)**2
r=s.Matrix([1,-2]); W=s.diag(wp,wa)
L=k*s.Matrix([F]).jacobian([p,a])
L=r*L.subs({p:P,a:A})
for i in range(2):
    for j in range(2):
        zero('jacobian_factor_%d%d'%(i,j),(L-k*Q*r*r.T*W.inv())[i,j])
energy=s.Matrix([2,1])
for i,v in enumerate(energy.T*L): zero('energy_left_kernel_%d'%i,v)
for i,v in enumerate(L*W*energy): zero('energy_right_kernel_%d'%i,v)
zero('nonzero_eigenvalue',s.trace(L)-k*(1+8*A+8*A*A)/(1+2*A))
zero('sample_eigenvalue',s.trace(L).subs({A:1,k:2})-s.Rational(34,3))
S=[s.Matrix([[1,s.I],[s.I,1]])/s.sqrt(2),s.Matrix([[1,1],[-1,1]])/s.sqrt(2),s.diag(1,s.I)]
for q in range(3): assert s.simplify(S[q]*S[q].conjugate().T)==s.eye(2)
V={(i,j,l):s.Integer(1+4*i+2*j+l)+s.I*(2+3*i-j+2*l) for i in range(2) for j in range(2) for l in range(2)}
Vm={(aa,b,cc):s.simplify(sum(s.conjugate(V[i,j,l])*S[0][i,aa]*S[1][j,b]*S[2][l,cc] for i in range(2) for j in range(2) for l in range(2))) for aa in range(2) for b in range(2) for cc in range(2)}
wrong=[]
for i in range(2):
    for b in range(2):
        for cc in range(2):
            gf=6*sum(V[i,j,l]*s.conjugate(S[1][j,b])*s.conjugate(S[2][l,cc]) for j in range(2) for l in range(2))
            gi=6*sum(Vm[aa,b,cc]*s.conjugate(S[0][i,aa]) for aa in range(2))
            zero('sewing_%d%d%d'%(i,b,cc),gi-s.conjugate(gf))
            bad=6*sum(Vm[aa,b,cc]*S[0][i,aa] for aa in range(2))
            wrong.append(s.simplify(bad-s.conjugate(gf)))
assert any(v!=0 for v in wrong)
T=s.Matrix([[1,s.I],[s.I,2]])
pair_strength=sum(s.simplify(s.Abs(6*T[i,j]/s.sqrt(1+int(i==j)))**2) for i in range(2) for j in range(i,2))
frob=sum(s.Abs(T[i,j])**2 for i in range(2) for j in range(2))
zero('symmetric_pair_strength',pair_strength-18*frob)
hb,pp,Df=s.symbols('hbar_ps pp Df',positive=True)
Cgamma=18*s.pi/(2*s.pi*hb)**2
zero('distinct_unit_conversion',8*s.pi*Cgamma*pp*Df-36*pp*Df/hb**2)
zero('repeated_unit_conversion',4*s.pi*Cgamma*pp*Df-18*pp*Df/hb**2)
# Exact non-geometric state with the thermal means used in B2.
zero('same_means_counterexample',av(wf-wr).subs({p:s.Rational(1,3),a:1}))
fixed_daughter=s.simplify((wf-wr).subs({X:s.Rational(1,3),Y:1}))
zero('same_means_fixed_daughter_flux',fixed_daughter-2*c)
result={'status':'PASS','scope':'Exact symbolic finite algebra only; no material run, source audit, novelty assessment, or global quantum semigroup construction.','zero_residual_checks':checks,'negative_control_missing_parent_conjugation':[str(v) for v in wrong],'symmetric_pair_strength':str(pair_strength),'same_means_fixed_daughter_flux':str(fixed_daughter),'sympy_version':s.__version__}
out=Path(__file__).with_suffix('.json')
out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':result['status'],'exact_zero_checks':len(checks),'wrong_conjugation_detected':True,'output':str(out)},indent=2))
