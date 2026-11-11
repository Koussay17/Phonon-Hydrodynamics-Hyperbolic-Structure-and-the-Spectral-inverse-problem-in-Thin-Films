"""Rewrite out/a4_wstar_table.csv from out/a4_fisher.json (no recomputation).

w*_1sigma = sigma(ln tau_R) and w*_2sigma = 2 sigma(ln tau_R) as relative changes of tau_R (H is linear in tau_R,
A5/A7). Usage: python a4_csv.py
"""
import json

from common import OUT, BANDS, SIG_PHI_DEG

d = json.loads((OUT / "a4_fisher.json").read_text())
rows = ["L_nm,variant,set,band,sigma_phi_deg,sigma_ln_tauR,wstar_1sigma_rel,wstar_2sigma_rel,cond,flag,indep_rel_diff"]
for key, c in d["cases"].items():
    name = key.split("_", 1)[1]
    for N, byband in c["sets"].items():
        for band in BANDS:
            for s in SIG_PHI_DEG:
                r = byband[band][str(s)]
                sg = r["sigma"]
                dd = c.get("indep_rel_diff", {}).get(f"{N}|{band}|{s}")
                rows.append(f"{c['L_nm']},{name},{N},{band},{s},{sg if sg is not None else 'NA'},"
                            f"{sg if sg is not None else 'NA'},{2 * sg if sg is not None else 'NA'},"
                            f"{r['cond']:.4g},{r['flag'].replace(',', ';')},{dd if dd is not None else ''}")
(OUT / "a4_wstar_table.csv").write_text("\n".join(rows) + "\n")
print(len(rows) - 1, "rows")
