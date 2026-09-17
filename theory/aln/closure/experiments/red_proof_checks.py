"""Independent exact proof probes; no campaign implementation is imported."""
from pathlib import Path
import json, sys
import sympy as s

out = {"python": sys.version, "sympy": s.__version__}
t, p, x = s.symbols("t p x", positive=True)
H2 = s.Matrix([[1, 1], [1, -1]])
S = s.kronecker_product(H2, H2, H2)
A = s.diag(0, t*t, 2, 2-t*t, 1, 1, 1, 1)
A[1,2] = A[2,1] = t
A[4,7] = A[7,4] = -t
C = s.simplify(S*A*S.T/8)
e, b = s.ones(8,1), S[:,2]
P = s.zeros(8)
for i in range(8):
    P[i,i ^ 3] = 1

def zero(M):
    return all(s.cancel(v) == 0 for v in M)

assert zero(S.T*S-8*s.eye(8))
assert zero(C*e)
assert all(s.cancel(C[i,i]-1) == 0 for i in range(8))
assert zero(P*C-C*P) and zero(P*b+b)
sol0 = -S[:,1]/t + b
assert zero(C*sol0-b) and s.cancel((e.T*sol0)[0]) == 0
K0 = s.simplify((b.T*sol0)[0]/8)
M2 = s.simplify((sol0.T*sol0)[0]/8)
assert K0 == 1 and s.cancel(M2-1-t**-2) == 0
d = p*p + (2+t*t)*p + t*t
solp = (-t*S[:,1] + (p+t*t)*b)/d
assert zero((C+p*s.eye(8))*solp-b)
J = S[:,2].row_join(S[:,5]).row_join(S[:,6])
Jsol = sol0.row_join(S[:,5]).row_join(S[:,6])
assert zero(C*Jsol-J)
assert zero(J.T*Jsol/8-s.eye(3))
assert zero(Jsol.T*Jsol/8-s.diag(1+t**-2,1,1))
cp = C.charpoly()
actual_poly = cp.as_expr().subs(cp.gen,x)
expected_poly = x*(x-(2-t*t))*(x-1)**2*((x-1)**2-t*t)*(x*x-(2+t*t)*x+t*t)
assert s.factor(actual_poly-expected_poly) == 0
out["D"] = {"K0": str(K0), "M2": str(M2),
    "response": str(s.simplify((b.T*solp)[0]/8)),
    "physical_characteristic_polynomial": str(s.factor(actual_poly)),
    "physical_diagonal_energy_parity_tensor_residuals": "identically zero"}

A1 = s.Matrix([[0,1,-1,0],[1,0,0,-1],[-1,0,0,1],[0,-1,1,0]])
A2 = s.Matrix([[0,1,0,-1],[1,0,-1,0],[0,-1,0,1],[-1,0,1,0]])
joint = A1.col_join(A2)
assert joint.rank() == 3 and joint.nullspace() == [s.ones(4,1)]
f = s.Matrix([3,-1,-1,-1])
first = -(f.T*A1*f)[0]
second = 2*s.Rational(3,4)*(A1*f).dot(A1*f)
assert first == 0 and second == 48
out["B1"] = {"four_coordinate_common_kernel": "span(1,1,1,1)",
    "first_derivative_zero_test": str(first), "second_derivative_same_test": str(second)}

Dm = s.diag(1,2,3,4)
H = s.Matrix([[1,0],[1,1],[1,-1],[1,2]])
Q = s.eye(4)-H*(H.T*H).inv()*H.T
CD = Dm-Dm*H*(H.T*Dm*H).inv()*H.T*Dm
G = Q*Dm.inv()*Q
assert zero(CD*G-Q) and zero(G*CD-Q)
assert zero(CD*G*CD-CD) and zero(G*CD*G-G)
assert zero(G-G.T) and zero(CD-CD.T)
out["B_surrogate"] = {"rational_Moore_Penrose_conditions": "all exact",
    "diagonal": [str(CD[i,i]) for i in range(4)]}

z = s.Rational(1,7)+s.I*s.Rational(2,9)
U, B, W = s.Rational(1,3), s.Rational(2,5), -s.Rational(1,2)
full = s.Matrix([[z+2+s.I*U, 1+s.I*B],[1+s.I*B,z+1+s.I*W]])
Schur = z+2+s.I*U-(1+s.I*B)**2/(z+1+s.I*W)
static = z+1+s.I*(U-2*B)+B**2
assert s.cancel(full.inv()[0,0]-1/Schur) == 0
assert s.cancel(static-(z+s.Rational(29,25)-s.I*s.Rational(7,15))) == 0
out["A"] = {"exact_inverse_Schur_identity": "zero residual", "static_sign_check": str(static)}

q = s.symbols("q", positive=True)
r = q*q*(1+q*s.sin(1/q))
dr = s.diff(r,q)
assert s.simplify(dr-(2*q+3*q*q*s.sin(1/q)-q*s.cos(1/q))) == 0
ratio = (q*q/dr)/(s.sqrt(r)/2)
even = s.simplify(ratio.subs(q,1/(10*s.pi)))
odd = s.simplify(ratio.subs(q,1/(11*s.pi)))
assert even == 2 and odd == s.Rational(2,3)
out["C_density_counterexample"] = {"rate": str(r), "derivative": str(dr),
    "even_sequence_density_ratio": str(even), "odd_sequence_density_ratio": str(odd),
    "satisfies_value_asymptotic": True, "pointwise_density_asymptotic": False,
    "moment_thresholds_and_cumulative_tail_conclusions": "not refuted"}

w, tau = [s.Integer(1),s.Integer(2)], [s.Integer(1),s.Integer(3)]
s0 = s.Integer(2)
k0 = sum(wi*ti for wi,ti in zip(w,tau))
k = sum(wi*ti/(1+s0*ti) for wi,ti in zip(w,tau))
ts = k0/sum(w)
tm = sum(wi*ti*ti for wi,ti in zip(w,tau))/k0
lower, upper = k0/(1+s0*tm), k0/(1+s0*ts)
assert lower < k < upper
out["PI"] = {"lower": str(lower),"actual": str(k),"upper": str(upper)}

path = Path(__file__).with_suffix(".json")
path.write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
