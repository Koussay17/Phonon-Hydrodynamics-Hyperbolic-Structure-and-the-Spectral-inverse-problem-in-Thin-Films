"""Independent bounded proof probes of the frozen event implementation."""
from pathlib import Path
from fractions import Fraction
import importlib.util, json, sys, hashlib
import numpy as np
import mpmath as mp
import sympy as sp

root = Path(__file__).resolve().parents[1]
source = root / "review-inputs" / "collision_events.py"
spec = importlib.util.spec_from_file_location("frozen_event_proof_target", source)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
Op = mod.DecayEventOperator

def dense(op):
    return np.column_stack([op.action(v) for v in np.eye(len(op.energy_over_kbt))])

energies = np.array([1e16, 1., 1e16])
op = Op(energies, [[0,1,2]], [1.])
c = dense(op)
exact_delta = Fraction.from_float(energies[0]) - Fraction.from_float(energies[1]) - Fraction.from_float(energies[2])
assert exact_delta == -1 and op.detailed_balance_log_residual[0] == 0
try:
    Op(energies, [[0,2,1]], [1.])
    swap_status = "accepted"
except ValueError as exc:
    swap_status = str(exc)
assert swap_status != "accepted"
with mp.workdps(100):
    xm = list(map(mp.mpf,energies))
    n = [1/mp.expm1(a) for a in xm]
    d = [mp.sqrt(a*(1+a)) for a in n]
    grad = [1+n[1]+n[2], n[0]-n[2], n[0]-n[1]]
    stoich = [-1,1,1]
    negative_jacobian = mp.matrix([[-stoich[i]*grad[j]*d[j]/d[i] for j in range(3)] for i in range(3)])
    expected = np.array(negative_jacobian.tolist(),dtype=float)
    ratio = float(mp.exp(1))

# Exact resonant algebra, independently differentiating the quadratic flux.
n0,n1,n2,G,u,v = sp.symbols("n0 n1 n2 G u v", positive=True)
flux = G*(n0*(1+n1+n2)-n1*n2)
N = [u*v/(1-u*v),u/(1-u),v/(1-v)]
sub = dict(zip([n0,n1,n2],N))
activity = G*N[0]*(1+N[1])*(1+N[2])
residuals = [sp.factor(sp.diff(flux,ni).subs(sub) + activity*si/(Ni*(1+Ni)))
             for ni,si,Ni in zip([n0,n1,n2],[-1,1,1],N)]
assert residuals == [0,0,0]

# Self-reciprocal distinct daughter event and exact odd invisibility.
self_good = Op([1.,1.,2.],[[2,0,1]],[1.])
self_double = Op([1.,1.,2.],[[2,0,1],[2,1,0]],[1.,1.])
cg,cd = dense(self_good),dense(self_double)
odd = np.array([1.,-1.,0.])/np.sqrt(2)
np.testing.assert_allclose(cd,2*cg,atol=1e-15)
assert np.linalg.norm((cd-cg)@odd)<1e-15
z,k = .25,.3
V = np.diag([1.,-1.,0.])
e = self_good.energy_over_kbt*self_good.entropy_scale
j = V@e
hg = j@np.linalg.solve(z*np.eye(3)+cg+1j*k*V,j)
hd = j@np.linalg.solve(z*np.eye(3)+cd+1j*k*V,j)
assert abs(hg-hd)>1e-3
out = {
    "python":sys.version,"numpy":np.__version__,"mpmath":mp.__version__,"sympy":sp.__version__,
    "frozen_source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
    "resonant_gradient_symbolic_residuals":list(map(str,residuals)),
    "admission_counterexample": {
        "energies":energies.tolist(),"event":[0,1,2],"gamma":1.,
        "exact_represented_input_detuning":str(exact_delta),
        "reported_detuning":op.detailed_balance_log_residual.tolist(),
        "reported_relative_residual":op.relative_energy_residual.tolist(),
        "true_Bose_forward_reverse_ratio":ratio,
        "swapped_daughters_status":swap_status,
        "returned_C":c.tolist(),"true_negative_entropy_jacobian_100_digits":expected.tolist(),
        "maximum_entry_error":float(np.max(np.abs(c-expected))),
        "relative_Frobenius_error":float(np.linalg.norm(c-expected)/np.linalg.norm(expected))},
    "self_reciprocal_counting": {
        "doubling_matrix_residual":float(np.linalg.norm(cd-2*cg)),
        "odd_change_norm":float(np.linalg.norm((cd-cg)@odd)),
        "finite_k_current_good":[float(hg.real),float(hg.imag)],
        "finite_k_current_double":[float(hd.real),float(hd.imag)]}
}
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2),encoding="utf-8")
print(json.dumps(out,indent=2))
