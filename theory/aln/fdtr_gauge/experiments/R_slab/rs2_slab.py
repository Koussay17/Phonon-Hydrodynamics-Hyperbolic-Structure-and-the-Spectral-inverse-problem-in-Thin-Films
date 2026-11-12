"""RS2: two-ports of the full-operator slab with diffuse black walls, and of its RTA reduction, for REF and the
witnesses (no regulator in the slab). One eigendecomposition per frequency serves all thicknesses.

usage: python -B rs2_slab.py <mesh> <group> <ops> <freqset> <name> [name ...]
  mesh   : m775 | m553
  group  : film2 (REF + 12 E33 optima, 7x7x5) | filmw (REF + 10 film witnesses, 7x7x5 and 5x5x3)
  ops    : comma list of full, rta
  freqset: sub (every third FDTR frequency + the last one + 100 Hz, 1 kHz, 2 kHz; for AAA interpolation)
           | all (the 80 FDTR frequencies + 100 Hz, 1 kHz, 2 kHz)
  (append "-nowalls" to the freqset to compute only black walls for the full operator)
Walls: full operator - black (audit convention), specular r = 0.5, 0.9, diffuse transmission t = 0.5, 0.1 (rs_core.
wall_operator, same at both faces); rta - black only.
Per (tensor, op): out/rs2_<mesh>_<group>_<op>_<freqset>_<name>.npz with f, Y[wall, d, f, 2, 2], Th[...] and the
per-frequency diagnostics; static film times (regularised full operator QDQ and E-phys completions; RTA V1) in the
companion JSON.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs_core as rc  # noqa: E402

FREQ = np.geomspace(1e4, 2e8, 80)                   # au5_fisher1d.FREQ
DS = (300e-9, 500e-9)
EXTRA = np.array([1e2, 1e3, 2e3])
WALLS_FULL = [("black",), ("specular", 0.5), ("specular", 0.9), ("diffuse", 0.5), ("diffuse", 0.1)]
WALLS_RTA = [("black",)]


def freqs(tag):
    """'all': 80 FDTR frequencies; 'sub' (= 'sub3') / 'subK': every K-th FDTR frequency and the last one; always with
    the three low frequencies 100 Hz, 1 kHz, 2 kHz (static and first-order two-port coefficients)."""
    if tag == "all":
        return np.concatenate([EXTRA, FREQ])
    k = 3 if tag == "sub" else int(tag[3:])
    sub = np.unique(np.r_[np.arange(0, 80, k), 79])
    return np.concatenate([EXTRA, FREQ[sub]])


def rta_film(cap, vz, tau, L):
    """RTA (V1, one local temperature incl. the regulator) film times on unbinned modes (au8_witness_slab.rta_film)."""
    te = 1.0 / (1.0 / tau + np.abs(vz) / L)
    S = cap.sum()
    lam = (cap * vz * vz * te).sum()
    s = lam / S
    tauR = (cap * vz * vz * te ** 2).sum() / lam
    taul = (cap * te * (vz * vz * te - s) ** 2).sum() / (S * s * s)
    return {"lam": float(lam), "tauR": float(tauR), "taul": float(taul), "delta": float(tauR - taul)}


def main():
    mesh, group, ops, ftag = sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4]
    walls_on = not ftag.endswith("-nowalls")
    ftag = ftag.replace("-nowalls", "")
    names = sys.argv[5:]
    t0 = time.perf_counter()
    G = rc.geometry(mesh)
    fr = freqs(ftag)
    print(mesh, group, ops, ftag, "modes", G["n"], "pairs", G["ip"].size, "local", G["i0"].size, "frequencies", fr.size,
          flush=True)
    for nm in names:
        C = rc.operator(G, group, nm)
        meta_p = OUT / f"rs2_{mesh}_{group}_static_{nm}.json"
        meta = {"name": nm, "mesh": mesh, "group": group, "C_tot": G["C_tot"]}
        kb = rc.kinetic_times(G, C)
        meta["bulk_full"] = kb
        r = np.diag(C)
        for L in DS:
            for var in ("QDQ", "phys"):
                meta[f"reg_full_{var}_{L * 1e9:.0f}"] = rc.kinetic_times(G, rc.regulated(G, C, L, var))
            meta[f"reg_rta_V1_{L * 1e9:.0f}"] = rta_film(G["cap"], G["vz"], 1.0 / r, L)
        meta["bulk_rta"] = rc.kinetic_times(G, rc.rta_operator(G, C))
        meta_p.write_text(json.dumps(meta, indent=1))
        for op in ops:
            Ch = C if op == "full" else rc.rta_operator(G, C)
            Bk = rc.blocks(G, Ch)
            walls = WALLS_FULL if (op == "full" and walls_on) else WALLS_RTA
            Y = np.zeros((len(walls), len(DS), fr.size, 2, 2), complex)
            Th = np.zeros_like(Y)
            dg = {k: np.zeros((len(walls), len(DS), fr.size)) for k in ("energy_balance_rel", "bc_residual_rel", "pivot_ratio_Pminus")}
            ei = {k: np.zeros(fr.size) for k in ("eig_backward_max", "eig_backward_diffusive", "polish_shift_rel",
                                                 "schur_energy_correction_rel", "lam_min_abs", "lam_max_abs")}
            for k, f in enumerate(fr):
                outs, info = rc.slab_paired(Bk, 2j * np.pi * f, DS, walls=walls)
                for w_, ow in enumerate(outs):
                    for i, (Y_, Th_, d_) in enumerate(ow):
                        Y[w_, i, k], Th[w_, i, k] = Y_, Th_
                        for q in dg:
                            dg[q][w_, i, k] = d_[q]
                for q in ei:
                    ei[q][k] = info[q]
            np.savez(OUT / f"rs2_{mesh}_{group}_{op}_{ftag}_{nm}.npz", f=fr, d=np.array(DS), Y=Y, Th=Th, sym=json.dumps(Bk["sym"]),
                     walls=json.dumps(walls),
                     **{f"dg_{q}": v for q, v in dg.items()}, **{f"ei_{q}": v for q, v in ei.items()})
            print(nm, op, "done", f"max energy balance {dg['energy_balance_rel'].max():.1e}",
                  f"eig backward max {ei['eig_backward_max'].max():.1e}", f"polish shift max {ei['polish_shift_rel'].max():.1e}",
                  f"G0(300,500) {Y[0, 0, 0, 0, 0].real:.4e} {Y[0, 1, 0, 0, 0].real:.4e}", f"[{time.perf_counter() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
