"""Independent attack of frozen event code; no repository imports or edits."""
import hashlib
import importlib.util
import json
import platform
from pathlib import Path
import numpy as np
import scipy
import mpmath as mp

BASE = Path(__file__).resolve().parents[1]
FROZEN = BASE / "review-inputs" / "collision_events.py"
spec = importlib.util.spec_from_file_location("red_counter_frozen_events", FROZEN)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
Op = module.DecayEventOperator

def matrix(op):
    return op.factor.T.toarray() @ op.factor.toarray()

def nonlinear_entropy_generator(x, gamma=1):
    # Differentiate the quadratic population flux independently of the factor.
    with mp.workdps(100):
        xx = [mp.mpf(float(a)) for a in x]
        n = [1 / mp.expm1(a) for a in xx]
        d = [mp.sqrt(a*(1+a)) for a in n]
        grad = [gamma*(1+n[1]+n[2]), gamma*(n[0]-n[2]),
                gamma*(n[0]-n[1])]
        signs = [-1, 1, 1]
        c = mp.matrix([[-signs[i]*grad[j]*d[j]/d[i]
                        for j in range(3)] for i in range(3)])
        return np.array(c.tolist(), dtype=float)

out = {
    "python": platform.python_version(), "numpy": np.__version__,
    "scipy": scipy.__version__, "mpmath": mp.__version__,
    "frozen_sha256": hashlib.sha256(FROZEN.read_bytes()).hexdigest(),
    "scope": "Frozen implementation only; synthetic inputs; no AlN data.",
}
x = np.array([float(2**53+2), 1., float(2**53)])
op = Op(x, [[0, 1, 2]], [1.])
cm = matrix(op)
ct = nonlinear_entropy_generator(x)
with mp.workdps(100):
    exact_delta = mp.mpf(float(x[0]))-mp.mpf(float(x[1]))-mp.mpf(float(x[2]))
try:
    Op(x, [[0, 2, 1]], [1.])
    swapped = {"accepted": True}
except ValueError as exc:
    swapped = {"accepted": False, "error": str(exc)}
out["lost_detuning"] = {
    "x": x.tolist(), "accepted_event": [0,1,2],
    "reported_relative_residual": op.relative_energy_residual.tolist(),
    "reported_thermal_detuning": op.detailed_balance_log_residual.tolist(),
    "exact_stored_float_detuning": str(exact_delta),
    "thermal_forward_reverse_ratio": float(mp.exp(-exact_delta)),
    "swapped_daughters": swapped,
    "factor_generator": cm.tolist(),
    "nonlinear_entropy_generator_100digits": ct.tolist(),
    "relative_frobenius_error": float(np.linalg.norm(cm-ct)/np.linalg.norm(ct)),
    "reference_nonsymmetry_norm": float(np.linalg.norm(ct-ct.T)),
    "parent_large_daughter_entry_model": float(cm[0,2]),
    "parent_large_daughter_entry_reference": float(ct[0,2]),
}
assert exact_delta == 1 and not swapped["accepted"]
assert op.detailed_balance_log_residual[0] == 0
assert np.linalg.norm(cm-ct)/np.linalg.norm(ct) > .4

# Negative control: exactly resonant large energies retain the contract.
xr = np.array([float(2**53+2), 2., float(2**53)])
opr = Op(xr, [[0,1,2]], [1.])
cr = nonlinear_entropy_generator(xr)
out["exact_resonance_control"] = {
    "relative_frobenius_error": float(np.linalg.norm(matrix(opr)-cr)/np.linalg.norm(cr)),
    "thermal_detuning": opr.detailed_balance_log_residual.tolist(),
}
assert out["exact_resonance_control"]["relative_frobenius_error"] < 1e-13

# Counting attack: structural and homogeneous odd checks miss a duplicate.
good = Op([1.,1.,2.], [[2,0,1]], [1.])
double = Op([1.,1.,2.], [[2,0,1],[2,1,0]], [1.,1.])
cg, cd = matrix(good), matrix(double)
odd = np.array([1.,-1.,0.])/np.sqrt(2)
energy = good.energy_over_kbt*good.entropy_scale
velocity = np.diag([1.,-1.,0.])
current = velocity@energy
z, k = .25, .3
hg = current@np.linalg.solve(z*np.eye(3)+cg+1j*k*velocity,current)
hd = current@np.linalg.solve(z*np.eye(3)+cd+1j*k*velocity,current)
out["duplicate_self_reciprocal"] = {
    "odd_action_change_norm": float(np.linalg.norm((cd-cg)@odd)),
    "full_change_norm": float(np.linalg.norm(cd-cg)),
    "energy_residual_good": float(np.linalg.norm(cg@energy)),
    "energy_residual_double": float(np.linalg.norm(cd@energy)),
    "finite_k_current_good": [float(hg.real),float(hg.imag)],
    "finite_k_current_double": [float(hd.real),float(hd.imag)],
    "finite_k_relative_change": float(abs(hd-hg)/abs(hg)),
}
assert np.linalg.norm((cd-cg)@odd) < 1e-14
assert abs(hd-hg)/abs(hg) > .007

# Sparse complex action at moderate scales: attempted contract attack fails.
y = np.array([.2+1j*.4, -.1+1j*.2, .3-1j*.5])
out["sparse_complex_control_error"] = float(np.linalg.norm(
    good.action(y)-cg@y))
assert out["sparse_complex_control_error"] < 1e-14
target = Path(__file__).with_suffix(".json")
target.write_text(json.dumps(out, indent=2)+"\n", encoding="utf-8")
print(json.dumps(out, indent=2))

