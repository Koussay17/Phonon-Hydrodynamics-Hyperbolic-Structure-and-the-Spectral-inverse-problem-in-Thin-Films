"""Extract the Phonon Olympics AlN benchmark numbers from the author spreadsheets.

Sources (commit 0640f07735059be9717a7565c2a0f22dc0da7a17, hashes in sources/olympics_phono3py/SHA256SUMS):
  phono3py team: AlN_kappa_withT_qgrid_313117.xlsx (kappa(T), 31x31x17, tetrahedron,
                 no isotope), AlN_Kappa_qgrid_convergence_check.xlsx (300 K mesh series),
                 AlN_summary_final.xlsx (settings), Kappa_timing_info.xlsx.
  ShengBTE team: Thermal conductivity vs Temperature.xlsx (columns RTA_a,b,c / iter_a,b,c).
Writes results/olympics_benchmark.json.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
from xlsx_dump import read_sheets  # noqa: E402

SRC = common.CAMPAIGN / "sources" / "olympics_phono3py"


def main():
    out = {"source_commit": "0640f07735059be9717a7565c2a0f22dc0da7a17",
           "units": "W/(m K)", "files": {}}
    for fn in ["AlN_kappa_withT_qgrid_313117.xlsx", "AlN_Kappa_qgrid_convergence_check.xlsx",
               "ShengBTE__Thermal_conductivity_vs_Temperature.xlsx", "Kappa_timing_info.xlsx"]:
        p = SRC / fn
        if not p.exists():
            p = SRC / ("Convergence_tests__" + fn)
        if not p.exists():
            p = SRC / ("Timing_info__" + fn)
        out["files"][fn] = {"path": str(p), "sha256": common.sha256(p)}

    # phono3py kappa(T) at 31x31x17: columns T | full xx yy zz | (blank) | RTA xx yy zz
    sh = read_sheets(Path(out["files"]["AlN_kappa_withT_qgrid_313117.xlsx"]["path"]))["Sheet1"]
    rows = []
    for r, c in sorted(sh.items()):
        T = c.get(1)
        if isinstance(T, float):
            rows.append({"T": T, "lbte_xx": c.get(2), "lbte_zz": c.get(4), "rta_xx": c.get(6), "rta_zz": c.get(8)})
    out["phono3py_kappa_T_31x31x17"] = rows

    # phono3py 300 K mesh series
    sh = read_sheets(Path(out["files"]["AlN_Kappa_qgrid_convergence_check.xlsx"]["path"]))["300K"]
    mesh_rows = []
    for r, c in sorted(sh.items()):
        if isinstance(c.get(1), float):
            mesh_rows.append({"mesh": [int(c[1]), int(c[2]), int(c[3])], "lbte_xx": c.get(4), "lbte_zz": c.get(5),
                              "rta_xx": c.get(7), "rta_zz": c.get(8)})
    out["phono3py_mesh_series_300K"] = mesh_rows

    # timing (hours, 128 cores, phono3py 2.1.0)
    sh = read_sheets(Path(out["files"]["Kappa_timing_info.xlsx"]["path"]))["Sheet1"]
    out["phono3py_lbte_timing_hours_128cores"] = [
        {"mesh": [int(c[1]), int(c[2]), int(c[3])], "hours": c[4]}
        for r, c in sorted(sh.items()) if isinstance(c.get(1), float)]

    # ShengBTE kappa(T): find header row with T
    sheets = read_sheets(Path(out["files"]["ShengBTE__Thermal_conductivity_vs_Temperature.xlsx"]["path"]))
    sheng = []
    for name, sh in sheets.items():
        for r, c in sorted(sh.items()):
            vals = [c.get(i) for i in range(max(c) + 1)]
            nums = [v for v in vals if isinstance(v, float)]
            if len(nums) >= 7 and 10 <= nums[0] <= 1000:
                sheng.append({"T": nums[0], "rta_a": nums[1], "rta_b": nums[2], "rta_c": nums[3],
                              "iter_a": nums[4], "iter_b": nums[5], "iter_c": nums[6]})
    out["shengbte_kappa_T"] = sheng
    out["paper_table_300K"] = {
        "reference": "McGaughey et al., J. Appl. Phys. 138, 135108 (2025), AlN table, 300 K (in-plane/cross-plane)",
        "phono3py": {"rta": [253, 232], "lbte": [285, 271]},
        "ShengBTE": {"rta": [271, 251], "lbte": [298, 291]},
        "ALAMODE": {"rta": [282, 263]}}
    dest = common.CAMPAIGN / "results" / "olympics_benchmark.json"
    common.write_json(dest, out)
    print(json.dumps({k: (v if not isinstance(v, list) else v[:3]) for k, v in out.items()}, indent=1)[:4000])


if __name__ == "__main__":
    main()
