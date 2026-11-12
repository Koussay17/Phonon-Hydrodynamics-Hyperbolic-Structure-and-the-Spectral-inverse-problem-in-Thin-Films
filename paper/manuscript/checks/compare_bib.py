"""Compare volume, first page and year of every refs.bib entry with its registry record."""
import json, re
from pathlib import Path
here = Path(__file__).resolve().parent
bib = (here.parent / "refs.bib").read_text(encoding="utf-8")
recs = json.loads((here / "doi_records.json").read_text(encoding="utf-8"))
issues = 0
for m in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", bib, flags=re.S):
    key, body = m.group(1).strip(), m.group(2)
    f = dict((k.lower(), v.strip()) for k, v in re.findall(r"(\w+)\s*=\s*\{(.*?)\}\s*,?\s*$", body, flags=re.M))
    r = recs[key]
    if "message" not in r.get("data", {}):
        print(f"{key:26s} (DataCite record; year {r['data']['data']['attributes'].get('publicationYear')}, bib {f.get('year')})")
        continue
    msg = r["data"]["message"]
    yrs = set()
    for k in ("published-print", "published-online", "issued"):
        if k in msg and msg[k].get("date-parts") and msg[k]["date-parts"][0][0]:
            yrs.add(str(msg[k]["date-parts"][0][0]))
    vol = msg.get("volume", "")
    page = msg.get("page", msg.get("article-number", "")).split("-")[0]
    bpage = f.get("pages", "").split("-")[0]
    prob = []
    if f.get("year") not in yrs:
        prob.append(f"year bib {f.get('year')} vs registry {sorted(yrs)}")
    if vol and f.get("volume") and f.get("volume") != vol:
        prob.append(f"volume bib {f.get('volume')} vs {vol}")
    if page and bpage and bpage != page:
        prob.append(f"page bib {bpage} vs {page}")
    if prob:
        issues += 1
        print(f"{key:26s} " + "; ".join(prob))
print(f"entries with discrepancies: {issues}")
