"""Reduce saved registry records to their bibliographic fields.

The raw Crossref/DataCite answers contain reference lists, licences and funder data that are not needed to
check a bibliography entry. This keeps, for each record: DOI, type, title, container title, authors, volume,
issue, pages or article number, dates, publisher and ISBN. Applied in place to the JSON files given on the
command line; verify_dois.py calls trim() before writing.

Usage (from this directory):  python -B trim_records.py doi_records.json doi_records_initial.json title_search_records.json
"""
import json
import sys

KEEP = ("DOI", "type", "title", "container-title", "author", "editor", "volume", "issue", "page",
        "article-number", "issued", "published-print", "published-online", "publisher", "ISBN")


def trim_message(m):
    out = {k: m[k] for k in KEEP if k in m}
    for who in ("author", "editor"):
        if who in out:
            out[who] = [{k: a[k] for k in ("family", "given", "name") if k in a} for a in out[who]]
    return out


def trim(obj):
    """Recursively replace Crossref 'message' payloads and DataCite attributes by their trimmed form."""
    if isinstance(obj, dict):
        if "message" in obj and isinstance(obj["message"], dict):
            return {**{k: v for k, v in obj.items() if k != "message"}, "message": trim_message(obj["message"])}
        if "attributes" in obj and isinstance(obj["attributes"], dict):
            a = obj["attributes"]
            keep = {k: a[k] for k in ("doi", "titles", "creators", "publisher", "publicationYear", "types") if k in a}
            return {"attributes": keep}
        if "DOI" in obj and "title" in obj:  # item of a Crossref search result
            return trim_message(obj)
        return {k: trim(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [trim(v) for v in obj]
    return obj


if __name__ == "__main__":
    for path in sys.argv[1:]:
        data = json.load(open(path, encoding="utf-8"))
        json.dump(trim(data), open(path, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
        print("trimmed", path)
