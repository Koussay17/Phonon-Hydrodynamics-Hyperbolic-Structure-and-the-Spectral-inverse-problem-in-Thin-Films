"""Review R4b: phono3py LBTE internal-option variants for REF and two witnesses (is_kappa_star, reducible vs IR
solver, RTA), to locate the origin of the witness kappa-ratio discrepancy with the event operator.
phono3py env. usage: python -B r4b_p3_variants.py <outjson> name=path ...
"""
import json, os, sys, time
from pathlib import Path
os.environ.setdefault("OMP_NUM_THREADS", "2"); os.environ.setdefault("RAYON_NUM_THREADS", "2")
sys.dont_write_bytecode = True
sys.path.insert(0, r"C:\Users\Koussay\ResearchLab\orchestration\campaigns\20261005-fc3-memory\scripts")
import numpy as np
import aln_common as common
from p3_events import load_ph3

outp = Path(sys.argv[1]); specs = sys.argv[2:]
common.verify_inputs()
ph3 = load_ph3(True, [5, 5, 3], 0.1)
fc3_ref = np.array(ph3.fc3)
variants = {"LBTE_red_star": dict(is_LBTE=True, is_reducible_collision_matrix=True, is_kappa_star=True),
            "LBTE_red_nostar": dict(is_LBTE=True, is_reducible_collision_matrix=True, is_kappa_star=False),
            "LBTE_ir": dict(is_LBTE=True, is_reducible_collision_matrix=False),
            "RTA": dict(is_LBTE=False)}
res = {}
for spec in specs:
    name, path = spec.split("=", 1)
    fc3 = fc3_ref if path == "REF" else np.load(path)
    ph3.fc3 = np.ascontiguousarray(fc3, dtype="double")
    ph3.init_phph_interaction()
    for vn, kw in variants.items():
        t0 = time.time()
        ph3.run_thermal_conductivity(temperatures=[300.0], log_level=0, **kw)
        k = np.array(ph3.thermal_conductivity.kappa)[0, 0]
        res.setdefault(vn, {})[name] = [float(k[0]), float(k[2])]
        print(vn, name, k[0], k[2], "%.1fs" % (time.time() - t0), flush=True)
for vn, d in res.items():
    if "REF" in d:
        for name, v in d.items():
            d[name] = {"kappa": v, "ratio": [v[0] / d["REF"][0] if isinstance(d["REF"], list) else None, None]}
outp.write_text(json.dumps(res, indent=1))
