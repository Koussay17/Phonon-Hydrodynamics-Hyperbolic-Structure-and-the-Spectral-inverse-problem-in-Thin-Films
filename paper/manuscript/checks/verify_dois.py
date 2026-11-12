"""Retrieve the registry record of every DOI in ../refs.bib.

Crossref is queried for journal and book DOIs; DataCite for arXiv DOIs (prefix 10.48550).
The raw records are written to doi_records.json and a compact table, one line per bibliography key,
to doi_records.txt, so that each entry can be compared field by field with its registry record.
The script exits with a nonzero status if any DOI does not resolve.

Usage (from this directory):  python -B verify_dois.py
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BIB = Path(__file__).resolve().parent.parent / "refs.bib"


def bib_dois():
    text = BIB.read_text(encoding="utf-8")
    pairs = []
    for m in re.finditer(r"@\w+\{([^,]+),(.*?)\n\}", text, flags=re.S):
        key, body = m.group(1).strip(), m.group(2)
        d = re.search(r"doi\s*=\s*\{([^}]*)\}", body)
        pairs.append((key, d.group(1).strip() if d else None))
    return pairs


def fetch(doi):
    if doi.lower().startswith("10.48550/"):
        url = "https://api.datacite.org/dois/" + urllib.parse.quote(doi, safe="")
    else:
        url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/:")
    req = urllib.request.Request(url, headers={"User-Agent": "bibliography-check/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=40) as r:
            return {"status": r.status, "data": json.loads(r.read().decode("utf-8"))}
    except Exception as exc:  # record the failure, do not hide it
        return {"status": "error", "error": repr(exc)}


def summarise(key, doi, rec):
    if doi is None:
        return f"{key}\tNO DOI FIELD"
    if rec.get("status") != 200:
        return f"{key}\t{doi}\tNOT FOUND\t{rec.get('error', rec.get('status'))}"
    d = rec["data"]
    if "message" in d:  # Crossref
        m = d["message"]
        title = " ".join(m.get("title", [""]))
        cont = " ".join(m.get("container-title", [""]))
        auth = "; ".join(
            f"{a.get('family', a.get('name', ''))}, {a.get('given', '')}".strip(", ")
            for a in m.get("author", m.get("editor", [])))
        year = None
        for k in ("published-print", "published-online", "issued", "published"):
            if k in m and m[k].get("date-parts") and m[k]["date-parts"][0][0]:
                year = m[k]["date-parts"][0][0]
                break
        return "\t".join(str(x) for x in (
            key, doi, auth, title, cont, m.get("volume", ""), m.get("issue", ""),
            m.get("page", m.get("article-number", "")), year, m.get("type", ""), m.get("publisher", "")))
    a = d["data"]["attributes"]  # DataCite
    auth = "; ".join(c.get("name", "") for c in a.get("creators", []))
    title = " ".join(t.get("title", "") for t in a.get("titles", []))
    return "\t".join(str(x) for x in (key, doi, auth, title, a.get("publisher", ""), "", "", "",
                                      a.get("publicationYear", ""), "preprint", "DataCite"))


def main():
    records, lines = {}, []
    for key, doi in bib_dois():
        rec = fetch(doi) if doi else {"status": "no doi"}
        records[key] = {"doi": doi, **rec}
        lines.append(summarise(key, doi, rec))
        time.sleep(0.4)
    here = Path(__file__).resolve().parent
    sys.path.insert(0, str(here))
    from trim_records import trim  # keep only the bibliographic fields of each record
    (here / "doi_records.json").write_text(json.dumps(trim(records), indent=1, ensure_ascii=False), encoding="utf-8")
    header = "key\tdoi\tauthors\ttitle\tcontainer\tvolume\tissue\tpages\tyear\ttype\tpublisher\n"
    (here / "doi_records.txt").write_text(header + "\n".join(lines) + "\n", encoding="utf-8")
    bad = [l for l in lines if "NOT FOUND" in l or "NO DOI" in l]
    print(f"{len(lines) - len(bad)} of {len(lines)} bibliography DOIs resolved")
    for l in bad:
        print(l)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
