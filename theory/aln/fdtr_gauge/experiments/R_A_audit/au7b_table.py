"""AU7b: final margin table for the report (0.01 and 0.1 deg; note-16 and favourable; 300 and 500 nm).
Readings: (A) campaign: first-order regularised GK film, campaign half-range; (B) + geometric mesh continuation of the
range; (C) + k^2-truncation transfer (au5b, exact regularised RTA kernel; N1/N2 only); (D) kinetic slab: sigma at the
slab-fitted base point, visible half-range = campaign half-range x witness-specific range transfer (au8, 9x9x6), sets
N2lam (lam_z free, C known) and N3; (E) = (D) + mesh continuation. Output: out/au7b_table.json"""
import json
from pathlib import Path
OUT = Path(__file__).resolve().parent / "out"
s7 = json.loads((OUT / "au7_summary.json").read_text())["margins"]
a8 = json.loads((OUT / "au8_witness_slab.json").read_text())["m996"]
a4 = json.loads((OUT / "au4_delta_range.json").read_text())["table"]
rows = []
for L in (300, 500):
    t = a4[f"{float(L):.6g}"]["all"]
    fmesh = t["x12_geometric"]["halfrange_delta_over_tau0"] / t["robust"]["halfrange_delta_over_tau0"]
    for R in ("note16", "favourable"):
        m = s7[f"{L}|{R}"]
        Tw = {"N2lam": a8[f"{L}|note16"]["summary"]["range_transfer_N2lam"], "N3": a8[f"{L}|note16"]["summary"]["range_transfer_N3"]}
        for sp in ("0.01", "0.1"):
            r = {"L": L, "stack": R, "sigma_phi": sp}
            for N in ("N1", "N2", "N2lam", "N3"):
                x = m[f"{N}|{sp}"]
                r[f"A_{N}"] = x["m_campaign"]
                r[f"B_{N}"] = x["m_campaign"] / fmesh
                if x["T_ktrunc"] is not None and N in ("N1", "N2"):
                    r[f"C_{N}"] = x["m_campaign"] / fmesh / x["T_ktrunc"]
                if N in ("N2lam", "N3"):
                    D = x["sigma_delta_at_slab_base_ps"] / (x["halfrange_campaign_ps"] * Tw[N])
                    r[f"D_{N}"] = D
                    r[f"E_{N}"] = D / fmesh
            rows.append(r)
            print(L, R, sp, {k: round(v, 2) for k, v in r.items() if k[:2] in ("A_", "B_", "C_", "D_", "E_")})
(OUT / "au7b_table.json").write_text(json.dumps({"rows": rows}, indent=1))
