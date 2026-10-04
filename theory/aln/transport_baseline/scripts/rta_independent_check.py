"""Independent re-summation of the RTA lattice thermal conductivity.

Reads only primitive per-mode arrays saved by phono3py (frequency, group velocity,
gamma = imaginary self-energy / half linewidth, star weights) plus the POSCAR
lattice, and recomputes

    kappa_ab = 1/(N V) * sum_{q in BZ, j} C_qj v_a v_b tau_qj,
    tau_qj   = 1 / (2 * 2*pi * Gamma_qj)      (Gamma in THz, ordinary frequency),
    C_qj     = k_B x^2 e^x/(e^x-1)^2,          x = h f_qj/(k_B T),

WITHOUT phono3py's rotation machinery: for the wurtzite point group 6mm, the star sum
of v (x) v over the arms of an irreducible point equals w * diag((vx^2+vy^2)/2,
(vx^2+vy^2)/2, vz^2) (3-fold axis => in-plane isotropy of rank-2 tensors).
Heat capacities are recomputed from the frequencies (not read from the file) with two
constant sets:
  * "phonopy" constants (as used by phono3py 4.5.0 internally: CODATA-2006-era
     k_B, h and e = 1.60217733e-19 C) -> must reproduce phono3py to round-off;
  * exact SI-2019 constants -> differs by the documented constant-set offset.
Modes with f < 1e-4 THz (phono3py cutoff_frequency) are excluded, as in phono3py.

usage: python rta_independent_check.py <kappa-*.hdf5> [<POSCAR>] [--json out.json]
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np

CONST = {
    "phonopy": {"kB_eV": 8.617338256808316e-05, "h_eVs": 4.13566733e-15, "e_C": 1.60217733e-19},
    "SI2019": {"kB_eV": 1.380649e-23 / 1.602176634e-19, "h_eVs": 6.62607015e-34 / 1.602176634e-19,
               "e_C": 1.602176634e-19},
}
CUTOFF_THZ = 1e-4


def poscar_volume(path: Path) -> float:
    lines = Path(path).read_text().splitlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    return abs(float(np.linalg.det(lat)))


def mode_cv(freq_THz: np.ndarray, T: float, c: dict) -> np.ndarray:
    """Mode heat capacity in eV/K; zero below the frequency cutoff."""
    out = np.zeros_like(freq_THz)
    m = freq_THz >= CUTOFF_THZ
    x = c["h_eVs"] * freq_THz[m] * 1e12 / (c["kB_eV"] * T)
    ex = np.exp(x)
    out[m] = c["kB_eV"] * x**2 * ex / (ex - 1.0) ** 2
    return out


def independent_kappa(freq, gv, gamma_T, weights, T, N, V_A3, c):
    cv = mode_cv(freq, T, c)                      # (nir, nb) eV/K
    valid = freq >= CUTOFF_THZ
    with np.errstate(divide="ignore", invalid="ignore"):
        tau_ps = np.where(valid & (gamma_T > 0), 1.0 / (4.0 * math.pi * gamma_T), 0.0)  # ps
    vx2y2 = 0.5 * (gv[..., 0] ** 2 + gv[..., 1] ** 2)   # (THz A)^2
    vz2 = gv[..., 2] ** 2
    w = weights[:, None]
    # SI: cv*e [J/K] * v^2*(1e2)^2 [m^2/s^2] * tau*1e-12 [s] / (N V*1e-30 [m^3])
    fac = c["e_C"] * 1e4 * 1e-12 / (N * V_A3 * 1e-30)
    kxx = fac * np.sum(w * cv * vx2y2 * tau_ps)
    kzz = fac * np.sum(w * cv * vz2 * tau_ps)
    mode_xx = fac * w * cv * vx2y2 * tau_ps
    mode_zz = fac * w * cv * vz2 * tau_ps
    return kxx, kzz, cv, mode_xx, mode_zz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kappa_file")
    ap.add_argument("poscar", nargs="?", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    kf = Path(a.kappa_file)
    poscar = Path(a.poscar) if a.poscar else None
    if poscar is None:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import common
        poscar = common.DATASET / "POSCAR"
    V = poscar_volume(poscar)
    with h5py.File(kf, "r") as f:
        freq = f["frequency"][()]
        gv = f["group_velocity"][()]
        gamma = f["gamma"][()]
        weights = f["weight"][()].astype(float)
        temps = f["temperature"][()]
        mesh = f["mesh"][()]
        cv_file = f["heat_capacity"][()]
        g_iso = f["gamma_isotope"][()] if "gamma_isotope" in f else None
        is_lbte = "kappa_RTA" in f
        k_ref = f["kappa_RTA"][()] if is_lbte else f["kappa"][()]
        mode_ref = f["mode_kappa_RTA"][()] if is_lbte else (f["mode_kappa"][()] if "mode_kappa" in f else None)
    N = int(np.prod(mesh))
    assert int(weights.sum()) == N, (weights.sum(), N)
    out = {"kappa_file": str(kf), "poscar": str(poscar), "volume_A3": V, "mesh": mesh.tolist(),
           "num_ir": int(len(weights)), "source": "kappa_RTA (LBTE file)" if is_lbte else "kappa (RTA file)",
           "isotope_rate_included": g_iso is not None,
           "rows": []}
    for it, T in enumerate(temps):
        if T <= 0:
            continue
        row = {"T": float(T), "phono3py_xx": float(k_ref[it, 0]), "phono3py_yy": float(k_ref[it, 1]),
               "phono3py_zz": float(k_ref[it, 2]),
               "phono3py_offdiag_max": float(np.abs(k_ref[it, 3:]).max())}
        for name, c in CONST.items():
            g_tot = gamma[it] + (g_iso.reshape(gamma[it].shape) if g_iso is not None else 0.0)
            kxx, kzz, cv, mxx, mzz = independent_kappa(freq, gv, g_tot, weights, float(T), N, V, c)
            row[f"indep_{name}_xx"] = kxx
            row[f"indep_{name}_zz"] = kzz
            row[f"rel_{name}_xx"] = (kxx - k_ref[it, 0]) / k_ref[it, 0]
            row[f"rel_{name}_zz"] = (kzz - k_ref[it, 2]) / k_ref[it, 2]
            m = cv_file[it] > 0
            row[f"cv_maxrel_{name}"] = float(np.max(np.abs(cv[m] - cv_file[it][m]) / cv_file[it][m]))
            if name == "phonopy" and mode_ref is not None:
                # per-ir-point star sums: compare (xx+yy)/2 and zz mode by mode
                ref_xy = 0.5 * (mode_ref[it, :, :, 0] + mode_ref[it, :, :, 1]) / N
                ref_zz = mode_ref[it, :, :, 2] / N
                scale = max(np.abs(ref_xy).max(), 1e-300)
                row["mode_maxabs_diff_over_max_xx"] = float(np.abs(mxx - ref_xy).max() / scale)
                row["mode_maxabs_diff_over_max_zz"] = float(np.abs(mzz - ref_zz).max() / max(np.abs(ref_zz).max(), 1e-300))
        out["rows"].append(row)
    txt = json.dumps(out, indent=2)
    print(txt)
    if a.json:
        Path(a.json).write_text(txt + "\n")


if __name__ == "__main__":
    main()
