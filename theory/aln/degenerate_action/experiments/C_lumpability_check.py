"""Exact algebra and finite linear dynamics for block lumpability.
No other branch, repository implementation, or external data are read.
"""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import scipy
import scipy.linalg as la
import sympy as s

OUT=Path(__file__).resolve().parent
W=s.diag(s.Rational(6,25),s.Rational(6,25),2,s.Rational(3,4))
B=s.Matrix([[1,1,0,0],[0,0,1,0],[0,0,0,1]])
nu=[s.Matrix([-1,0,1,1]),s.Matrix([0,-1,1,1])]
Sigma=B*W*B.T
H=W*B.T*Sigma.inv()
P=H*B
Q=s.eye(4)-P
h=s.Matrix([1,-1,0,0])
u=s.Matrix([-1,1,1])
F=s.Rational(3,5)
rates={"equal":[2,2],"unequal":[3,1]}
models={}
for name,gamma in rates.items():
    M=F*sum((gamma[i]*nu[i]*nu[i].T for i in range(2)),s.zeros(4))
    L=M*W.inv()
    G=B*L*H
    models[name]={"M":M,"L":L,"G":G}
G=models["equal"]["G"]
Lm=models["unequal"]["L"]
mzero=lambda z: z==s.zeros(*z.shape)
n0=s.Matrix([s.Rational(1,5),s.Rational(1,5),1,s.Rational(1,2)])
eta=s.Rational(1,100)
population=n0+eta*h
rhs=s.zeros(4,1)
for i,gamma in enumerate(rates["unequal"]):
    flux=population[i]*(1+population[2])*(1+population[3])-(1+population[i])*population[2]*population[3]
    rhs+=gamma*nu[i]*flux
rotation=s.Matrix([[s.cos(s.pi/12),s.sin(s.pi/12)],[-s.sin(s.pi/12),s.cos(s.pi/12)]])
rotation_res=(rotation*s.Matrix([s.sqrt(2),s.sqrt(2)])-s.Matrix([s.sqrt(3),1])).applyfunc(s.simplify)
energy=s.Matrix([s.log(2)+s.log(3),s.log(2)+s.log(3),s.log(2),s.log(3)])
checks={
    "same_compressed_G":mzero(models["unequal"]["G"]-G),
    "equal_exactly_lumpable":mzero(B*models["equal"]["L"]-G*B),
    "unequal_not_lumpable":not mzero(B*Lm-G*B),
    "hidden_same_block_totals":mzero(B*h),
    "hidden_aggregate_derivative":mzero(-B*Lm*h-5*u),
    "equal_occupation_leakage":mzero(Q*Lm*H-h*s.Matrix([[s.Rational(5,4),-s.Rational(3,10),-s.Rational(4,5)]])),
    "second_derivative_defect_identity":mzero(B*Lm**2*H-G**2-B*Lm*Q*Lm*H),
    "nonlinear_positive_state_same_derivative":mzero(B*rhs-eta*5*u),
    "positive_perturbed_populations":all(x>0 for x in population),
    "rotation_preserves_strengths":mzero(rotation_res),
    "both_conserve_energy":all(mzero(energy.T*model["M"]) for model in models.values())
}
numeric=[]
bn=np.array(B,float)
hn=np.array(H,float)
gn=np.array(G,float)
N0=np.array([.02,0.,0.])
for name,model in models.items():
    ln=np.array(model["L"],float)
    for t in [.001,.01,.1,.5]:
        projected=bn@la.expm(-ln*t)@hn
        closed=la.expm(-gn*t)
        actualN=projected@N0
        closedN=closed@N0
        numeric.append(dict(model=name,t=t,
            propagator_relative_error=float(la.norm(projected-closed,2)/la.norm(projected,2)),
            total_actual=actualN.tolist(),total_compressed=closedN.tolist(),
            selected_initial_state_relative_error=float(la.norm(actualN-closedN)/la.norm(actualN))))
checks["numeric_equal_closure"]=max(r["propagator_relative_error"] for r in numeric if r["model"]=="equal")<1e-13
checks["numeric_unequal_failure"]=max(r["propagator_relative_error"] for r in numeric if r["model"]=="unequal")>1e-3

def entries(a):
    return [[str(a[i,j]) for j in range(a.cols)] for i in range(a.rows)]

data=dict(metadata=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
    sympy=s.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope="Declared finite reversible scalar event model; no quantum generator or material input"),
    exact=dict(W=entries(W),B=entries(B),Sigma=entries(Sigma),H=entries(H),G=entries(G),
        mean_generator=entries(s.diag(s.Rational(1,2),1,1)*G*s.diag(2,1,1)),
        unequal_closure_defect=entries(B*Lm-G*B),
        hidden_aggregate_derivative=entries(-B*Lm*h),
        equal_occupation_leakage=entries(Q*Lm*H),
        second_derivative_defect=entries(B*Lm**2*H-G**2),
        nonlinear_block_derivative=entries(B*rhs)),
    dynamics=numeric,checks=checks)
(OUT/"C_lumpability_results.json").write_text(json.dumps(data,indent=2),encoding="utf-8")
print(json.dumps(data,indent=2))
assert all(checks.values()),checks

