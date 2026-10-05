"""Arithmetic check of the vacuity claims about the cone bound (13-final-report.md section 1, bullet 3;
07-experiments.md E2): "M^2 K0 >= 1e8 tau_ref in every geometry computed" and P6 "1D N = 9: 7.3e10".
Only the campaign's own result files are used (results/bt_local_search.json, results/bt_exact_d1_N9.json,
results/scan_summary.json).
Run: python -B check_bt_numbers.py  (writes check_bt_numbers.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"


def main():
    rows = json.loads((RES / "bt_local_search.json").read_text())
    out = {"per_geometry": []}
    for r in rows:
        ratio = r["M_lower"] ** 2 * r["K0"] / r["tau_ref"]
        out["per_geometry"].append({
            "geometry": f"d={r['d']} N={r['N']}", "M_lower": r["M_lower"], "K0": r["K0"], "tau_ref": r["tau_ref"],
            "sqrt_tau_ref_over_K0": (r["tau_ref"] / r["K0"]) ** 0.5,
            "M2K0_over_tau_ref_certified_lower": ratio, "claim_ge_1e8_supported": ratio >= 1e8})
    out["claim_holds_for_all"] = all(x["claim_ge_1e8_supported"] for x in out["per_geometry"])
    out["min_ratio"] = min(x["M2K0_over_tau_ref_certified_lower"] for x in out["per_geometry"])
    ex = json.loads((RES / "bt_exact_d1_N9.json").read_text())
    # P6 value M^2 |b|^2 for 1D N = 9 needs |b|^2: recover it from tau_CS = K0/|b|^2 is not stored; recompute.
    sys.path.insert(0, str(HERE.parent / "scripts"))
    from debye_events import build  # noqa: E402
    g = build(1, 9)
    b2 = float(g.b[:, 0] @ g.b[:, 0])
    out["P6_1D_N9"] = {"M": ex["M"], "b2": b2, "M2_b2": ex["M"] ** 2 * b2, "claimed": 7.3e10,
                       "tau_CS": ex["K0"] / b2}
    # the 1D N = 15 and 21 pivoting lower bounds are far below the exact 1D N = 9 value: the pivoting search
    # is not monotone in N and these M_lower are weak
    (HERE / "check_bt_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
