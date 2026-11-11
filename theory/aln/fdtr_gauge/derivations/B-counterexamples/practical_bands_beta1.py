"""Thin-film exponents at matched effusivity (b = b_s): numerical check of the symbolic hierarchy
(s01_thin_hierarchy.py predicts cumulative-rank orders 1, 3, 4, 5, i.e. s_k ~ (omega_b xi1^2)^{1/2, 3/2, 2, 5/2})."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from practical_bands import band_svd, fit
out = {}
for name, th in (("beta1_gk", (1.0, 1.0, 0.5, 0.1)), ("beta1_cattaneo", (1.0, 1.0, 0.5, 0.0))):
    rows = []
    for wb in 10.0 ** np.arange(-5, -2.49, 0.5):
        S, Vt = band_svd(th, wb)
        rows.append([wb, *S.tolist()])
    rows = np.array(rows)
    out[name] = {"rows": rows.tolist(), "exponents": [fit(rows[:, 0], rows[:, 1 + k])[0] for k in range(4)]}
    print(name, "exponents s1..s4 vs omega_b xi1^2:", ["%.3f" % e for e in out[name]["exponents"]])
json.dump(out, open(os.path.join(HERE, "data", "practical_bands_beta1.json"), "w"), indent=1)
