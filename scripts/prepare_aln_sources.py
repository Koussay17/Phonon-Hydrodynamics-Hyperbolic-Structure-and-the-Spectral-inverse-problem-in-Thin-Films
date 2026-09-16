"""Fetch and convert published AlN inputs; external source files stay in .build.

No DFT calculation is run. Rao et al. data are CC BY 4.0; the compact NPZ is
an explicitly attributed conversion, not newly computed ab initio data.
"""
from pathlib import Path
import hashlib
import io
import json
import urllib.request
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / ".build/research"
OUT = ROOT / "theory/aln"
RAO_URL = "https://data.mendeley.com/public-files/datasets/w9hg2mnnwy/files/c91f4a36-8294-4e63-99ee-a2d1528d21ba/file_downloaded"
RAO_SHA = "ee51e6c7e7c51105504bed98298f96826fcd56ac5542c686f6955d33496d0709"
COMMIT = "0640f07735059be9717a7565c2a0f22dc0da7a17"
OLY_PATH = "Aluminum Nitride/ShengBTE/AlN_final result/5.thermal conductivity/Thermal conductivity vs Temperature.xlsx"
OLY_URL = "https://raw.githubusercontent.com/McGaughey-Lab/Phonon-Olympics/" + COMMIT + "/" + urllib.parse.quote(OLY_PATH)


def fetch(name, url, expected=None):
    path = CACHE/name
    if not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        path.write_bytes(urllib.request.urlopen(request, timeout=60).read())
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if expected and digest != expected:
        raise ValueError("source checksum mismatch: " + name)
    return raw, digest


def read_temperature_table(raw):
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        strings = ["".join(node.itertext()) for node in
                   ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns)]
        root = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in root.findall(".//m:row", ns):
            cells = {}
            for cell in row.findall("m:c", ns):
                value = cell.find("m:v", ns)
                if value is None:
                    continue
                value = strings[int(value.text)] if cell.get("t") == "s" else value.text
                col = "".join(c for c in cell.get("r") if c.isalpha())
                cells[col] = value
            if not cells:
                continue
            if not rows:
                assert cells["B"] == "RTA-a (W/mK)" and cells["G"] == "Iter-c(W/mK)"
                rows.append(None)
                continue
            if "A" in cells:
                rows.append([None if cells[c] == "NA" else float(cells[c])
                             for c in "ABCDEFG"])
    return rows[1:]


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    raw, digest = fetch("rao_2025.zip", RAO_URL, RAO_SHA)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        def get(name):
            return np.loadtxt(io.BytesIO(z.read("ShengBTE file/AlN/"+name)))
        q = get("BTE.qpoints")
        omega = get("BTE.omega")*1e12
        nq, nb = omega.shape
        assert nb == 12 and nq == 793
        velocity = get("BTE.v").reshape(nb, nq, 3).transpose(1, 0, 2)*1e3
        rate_files = {"total_rate_s": "T300K/BTE.w",
                      "anharmonic_rate_s": "T300K/BTE.w_anharmonic",
                      "isotope_rate_s": "BTE.w_isotopic"}
        rates = {}
        for key, name in rate_files.items():
            values = get(name)
            np.testing.assert_allclose(values[:, 0].reshape(nb, nq).T*1e12,
                                       omega, rtol=1e-9, atol=1.)
            rates[key] = values[:, 1].reshape(nb, nq).T*1e12
        reciprocal = get("BTE.ReciprocalLatticeVectors")
        volume = (2*np.pi)**3 / abs(np.linalg.det(reciprocal))*1e-27
        degeneracy = q[:, 2].astype(int)
        # Verify irreducible degeneracies against the full mesh.
        full = get("BTE.qpoints_full")
        counts = np.bincount(full[:, 1].astype(int)-1, minlength=nq)
        np.testing.assert_array_equal(counts, degeneracy)
        np.testing.assert_allclose(rates["total_rate_s"],
            rates["anharmonic_rate_s"]+rates["isotope_rate_s"], rtol=2e-8, atol=1e3)
        refs = {"C_J_m3K": float(get("T300K/BTE.cv")),
                "kappa_RTA_W_mK": get("BTE.KappaTensorVsT_RTA")[1:10].reshape(3,3).tolist(),
                "kappa_iterative_W_mK": get("BTE.KappaTensorVsT_CONV")[1:10].reshape(3,3).tolist()}
    np.savez_compressed(OUT/"rao_300K_modes.npz", omega_rad_s=omega,
                       velocity_m_s=velocity, degeneracy=degeneracy,
                       q_fractional=q[:,3:6], cell_volume_m3=volume, **rates)
    table_raw, table_sha = fetch("olympics_kappa.xlsx", OLY_URL, "b7c07d2170e9047250f9fcd1da4a188f4086a452012746364e042dff1309cf67")
    table = {"source_url": OLY_URL, "source_commit": COMMIT, "sha256": table_sha,
             "columns": ["T_K","RTA_a","RTA_b","RTA_c","iter_a","iter_b","iter_c"],
             "units": "W/(m K) except T_K", "rows": read_temperature_table(table_raw),
             "scope": "Published bulk calculation, not our rerun. CONTROL disables isotope scattering. 20 K unavailable; low-T convergence not certified here."}
    (OUT/"olympics_temperature.json").write_text(json.dumps(table, indent=2)+"\n", encoding="utf-8")
    manifest = {"rao": {"doi": "10.17632/w9hg2mnnwy.1", "url": RAO_URL,
                       "archive_sha256": digest, "license": "CC BY 4.0",
                       "temperature_K": 300, "irreducible_q": nq, "branches": nb,
                       "full_q": int(degeneracy.sum()), "references": refs,
                       "normal_umklapp_split_available": False,
                       "conversion": "rad/ps -> rad/s, km/s -> m/s, ps^-1 -> s^-1. q index first then branch. Original arrays retained, no fitted rates."},
                "olympics": {"doi": "10.1063/5.0289819", "url": OLY_URL,
                             "sha256": table_sha, "commit": COMMIT},
                "missing": ["mode-resolved N/U separation versus temperature",
                            "full collision operator for dynamic GK closure",
                            "mesh convergence for the Rao dataset and derived higher moments",
                            "sample defect populations and boundary properties"]}
    (OUT/"sources.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(manifest,indent=2))


if __name__ == "__main__":
    main()
