"""Bibliographic search on Crossref for references whose DOI is not known in advance."""
import json, sys, time, urllib.parse, urllib.request
QUERIES = {
    "krapez2019": "Linear, trigonometric and hyperbolic profiles of thermal effusivity in the Liouville space and related quadrupoles",
    "krapez2023": "Krapez Comment on Simultaneous density and thermal conductivity depth profile reconstructions",
    "maillet2000": "Thermal quadrupoles solving the heat equation through integral transforms Maillet",
    "kovacs2024": "Kovacs Heat equations beyond Fourier from heat waves to thermal metamaterials",
    "cattaneo1958": "Cattaneo Sur une forme de l'equation de la chaleur eliminant le paradoxe d'une propagation instantanee",
    "vernotte1958": "Vernotte Les paradoxes de la theorie continue de l'equation de la chaleur",
    "golden1983": "Golden Papanicolaou Bounds for effective parameters of heterogeneous media by analytic continuation",
}
out = {}
for key, q in QUERIES.items():
    url = "https://api.crossref.org/works?rows=4&query.bibliographic=" + urllib.parse.quote(q)
    req = urllib.request.Request(url, headers={"User-Agent": "bibliography-check/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=40) as r:
            items = json.loads(r.read().decode())["message"]["items"]
    except Exception as exc:
        print(key, "ERROR", exc); continue
    out[key] = items
    print("==", key)
    for it in items:
        au = "; ".join(f"{a.get('family','')}, {a.get('given','')}" for a in it.get("author", [])[:3])
        yr = (it.get("issued", {}).get("date-parts") or [[None]])[0][0]
        print("  ", it.get("DOI"), "|", au, "|", " ".join(it.get("title", [""]))[:110], "|",
              " ".join(it.get("container-title", [""]))[:50], it.get("volume", ""), it.get("page", it.get("article-number", "")), yr)
    time.sleep(0.5)
json.dump(out, open("title_search_records.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
