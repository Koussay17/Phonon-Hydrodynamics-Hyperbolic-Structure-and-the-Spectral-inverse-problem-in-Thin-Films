"""Draft BibTeX entries from the initial registry records (doi_records_initial.json, first DOI list) and search_titles.py.

The output (refs_draft.bib) is reviewed by hand before being copied to ../refs.bib: title capitalization
is protected, LaTeX accents are restored and book metadata missing from the registry is completed only
from the registry records themselves.
"""
import json

recs = json.load(open("doi_records_initial.json", encoding="utf-8"))
titles = json.load(open("title_search_records.json", encoding="utf-8"))
extra = {}
for key, items in titles.items():
    for it in items:
        extra[it["DOI"].lower()] = {"status": 200, "data": {"message": it}}

KEYS = {
    "Cahill2004": "10.1063/1.1819431", "Schmidt2009": "10.1063/1.3212673",
    "KohCahill2007": "10.1103/PhysRevB.76.075207", "Regner2013": "10.1038/ncomms2630",
    "WilsonCahill2014": "10.1038/ncomms6075", "Minnich2012": "10.1103/PhysRevLett.109.205901",
    "Forghani2016": "10.1103/PhysRevB.94.155439", "Beardo2020": "10.1103/PhysRevB.101.075303",
    "HennessyMyers2021": "10.1007/978-3-030-64272-3_2", "CamachoDeLaRosa2025": "10.1063/5.0257299",
    "KrapezRigollet2017": "10.1063/1.4997733", "Krapez2019": "10.1016/j.ijthermalsci.2018.10.018",
    "Krapez2023": "10.1063/5.0151465", "Cheng2020": "10.1103/PhysRevMaterials.4.044602",
    "McGaughey2025": "10.1063/5.0289819", "MaLiLuo2014": "10.1103/PhysRevB.90.035203",
    "JosephPreziosi1989": "10.1103/RevModPhys.61.41", "GuyerKrumhansl1966a": "10.1103/PhysRev.148.766",
    "GuyerKrumhansl1966b": "10.1103/PhysRev.148.778", "Callaway1959": "10.1103/PhysRev.113.1046",
    "HardyAlbers1974": "10.1103/PhysRevB.10.3546", "GuoWang2022": "10.1016/j.ijthermalsci.2021.107178",
    "Sendra2021": "10.1103/PhysRevB.103.L140301", "Sendra2022": "10.1103/PhysRevB.106.155301",
    "GuoWang2015": "10.1016/j.physrep.2015.07.003", "Struchtrup2025": "10.1007/978-3-031-93918-1_27",
    "LebonDauby1990": "10.1103/PhysRevA.42.4710", "GurtinPipkin1968": "10.1007/BF00281373",
    "Both2016": "10.1515/jnet-2015-0035", "Kovacs2018": "10.1016/j.ijheatmasstransfer.2018.06.082",
    "Kovacs2024": "10.1016/j.physrep.2023.11.001", "CepellottiMarzari2016": "10.1103/PhysRevX.6.041013",
    "Simoncelli2020": "10.1103/PhysRevX.10.011019", "Chaput2013": "10.1103/PhysRevLett.110.265506",
    "Fugallo2013": "10.1103/PhysRevB.88.045430", "Ziman2001": "10.1093/acprof:oso/9780198507796.001.0001",
    "Zwanzig1961": "10.1103/PhysRev.124.983", "Mori1965": "10.1143/PTP.33.423",
    "MalviyaRavichandran2023": "10.1103/PhysRevB.108.155201", "ZengDong2026": "10.1103/hx28-4ksc",
    "WeiKuo1969": "10.1021/i160029a019", "LiRabitz1989": "10.1016/0009-2509(89)85014-6",
    "Aoki1968": "10.1109/TAC.1968.1098900", "Buchholz1994": "10.2307/3215235",
    "BenTalTeboulle1990": "10.1016/0024-3795(90)90395-S", "Berg1986": "10.1016/0024-3795(86)90307-1",
    "Stewart1989": "10.1016/0024-3795(89)90594-6", "OLeary1990": "10.1016/0024-3795(90)90056-I",
    "Todd1990": "10.1287/opre.38.6.1006", "HankeNeumann1993": "10.1016/0024-3795(93)90223-B",
    "VavasisYe1996": "10.1007/BF02592148", "Tuncel1999": "10.1007/s101070050087",
    "Saunderson2012": "10.1137/120872516", "Akhiezer2020": "10.1137/1.9781611976397",
    "Schmudgen2017": "10.1007/978-3-319-64546-9", "Vermeersch2015a": "10.1103/PhysRevB.91.085202",
    "Vermeersch2015b": "10.1103/PhysRevB.91.085203", "Herring1954": "10.1103/PhysRev.95.954",
    "Chen2025": "10.48550/arXiv.2502.19533", "Raue2009": "10.1093/bioinformatics/btp358",
    "Golden1983": "10.1007/bf01216179",
}


def rec_for(doi):
    r = recs.get(doi) or extra.get(doi.lower())
    if r is None or r.get("status") != 200:
        return None
    return r["data"]


out = []
for key, doi in KEYS.items():
    d = rec_for(doi)
    if d is None:
        out.append(f"% {key}: NO REGISTRY RECORD for {doi}\n")
        continue
    if "message" in d:
        m = d["message"]
        authors = " and ".join(
            (f"{a.get('family', '')}, {a.get('given', '')}" if a.get("given") else a.get("family", a.get("name", "")))
            for a in m.get("author", []))
        year = None
        for k in ("published-print", "issued", "published-online"):
            if k in m and m[k].get("date-parts") and m[k]["date-parts"][0][0]:
                year = m[k]["date-parts"][0][0]
                break
        typ = m.get("type")
        title = " ".join(m.get("title", [""]))
        cont = m.get("container-title", [""])
        fields = {"author": authors, "title": title, "doi": doi, "year": year}
        if typ == "journal-article":
            et = "article"
            fields["journal"] = cont[0] if cont else ""
            fields["volume"] = m.get("volume", "")
            fields["pages"] = m.get("page", m.get("article-number", ""))
        elif typ == "book-chapter":
            et = "incollection"
            fields["booktitle"] = cont[-1] if cont else ""
            fields["series"] = cont[0] if len(cont) > 1 else ""
            fields["pages"] = m.get("page", "")
            fields["publisher"] = m.get("publisher", "")
        else:
            et = "book"
            fields["publisher"] = m.get("publisher", "")
            if cont:
                fields["series"] = cont[0]
    else:
        a = d["data"]["attributes"]
        et = "misc"
        fields = {"author": " and ".join(c.get("name", "") for c in a.get("creators", [])),
                  "title": " ".join(t.get("title", "") for t in a.get("titles", [])),
                  "doi": doi, "year": a.get("publicationYear", ""), "howpublished": "arXiv preprint"}
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items() if v not in (None, ""))
    out.append(f"@{et}{{{key},\n{body}\n}}\n")
open("refs_draft.bib", "w", encoding="utf-8").write("\n".join(out))
print(f"{len(KEYS)} entries drafted")
