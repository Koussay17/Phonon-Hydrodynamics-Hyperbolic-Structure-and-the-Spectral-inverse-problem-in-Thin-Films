"""Write the inputs of the Julia audit model: frequency grid (bitwise the numpy grid of the branch), stack values,
and the AlN base values as re-assembled in aud1_base_values.py (not the branch file).
Usage: python aud0_inputs.py -> out/aud_inputs.json
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
b = json.loads((OUT / "aud1_base_values.json").read_text())
ref, conv, est = b["event_REF"], b["rta"]["m646435"], b["L500"]
base = {}
for L in ("30", "100", "300"):
    base[L] = {"tauR": ref[L]["tau_z"] * 1e-12, "lamz": ref[L]["kappa_z"], "lamr": conv["kappa_x"][L]}
base["500"] = {"tauR": est["tau_z"]["ratio_method_500"] * 1e-12, "lamz": est["kappa_z"]["ratio_method_500"],
               "lamr": conv["kappa_x"]["500"]}
base["500_rta"] = {"tauR": conv["tau_z"]["500"] * 1e-12, "lamz": conv["kappa_z"]["500"], "lamr": conv["kappa_x"]["500"]}
f80 = np.geomspace(1e4, 2e8, 80)
inp = {
    "freq80": [float(x) for x in f80],
    "stack": {"lamAu": 150.0, "CAu": 2.49e6, "dAu": 80e-9, "R1": 1e-8, "R2": 2e-8, "C": 2.41e6,
              "lamS": 35.0, "CS": 3.03e6, "wp": 10e-6, "wr": 8e-6},
    "base": base,
    "w_admissible": {"30": 0.028, "100": 0.037, "300": 0.052, "500": [0.055, 0.08]},
}
(OUT / "aud_inputs.json").write_text(json.dumps(inp, indent=1))
print(json.dumps(base, indent=1))


# Julia-readable copy of the same inputs (exact: repr of every float)
def jl(v):
    if isinstance(v, dict):
        return "Dict{String,Any}(" + ", ".join(f'"{k}" => {jl(x)}' for k, x in v.items()) + ")"
    if isinstance(v, (list, tuple)):
        return "Any[" + ", ".join(jl(x) for x in v) + "]"
    if isinstance(v, str):
        return '"' + v + '"'
    return repr(float(v))


(OUT / "aud_inputs.jl").write_text("const INP = " + jl(inp) + "\n")
