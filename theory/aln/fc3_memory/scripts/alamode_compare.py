"""Calibration of the per-cluster size eps with an independent workflow: the Phonon Olympics ALAMODE AlN cubic force
constants (cubic_15nn_ndata200_ols.xml, ALM 1.4.2, 4x4x3 supercell, own DFT setup) mapped onto the phono3py 3x3x2
supercell of this campaign and projected onto the same symmetric + ASR basis. phono3py environment.

Mapping (every step checked in the output json):
  * ALAMODE entries: pair1 = (primitive atom, xyz), pair2/pair3 = (supercell atom, xyz, image cell); stored with
    pair2 <= pair3 only -> the (pair2, pair3) swap is added. Image cell c (1-based) = shift number c-1 of the ALM
    ordering: 0 -> (0,0,0), then (ix, iy, iz) in -1..1 (ix slowest), (0,0,0) skipped. Checked: the stored image
    is a shortest image of the pair (fraction reported).
  * Polarity: the ALAMODE cell has N at +u c above Al, the phono3py cell at -u c; positions are mapped by z -> -z
    and tensors by R = diag(1,1,-1) (factor (-1)^(number of z indices)); equivalent to inversion up to the 6_3
    screw axis, a symmetry of P6_3mc.
  * Each atom gets an integer primitive label (lattice vector n, basis atom k); the ALAMODE triplet
    (k1, 0; k2, n2; k3, n3) is added to the phono3py compact fc3 at (p2s[k1], atom(n_p2s(k1) + n2, k2),
    atom(n_p2s(k1) + n3, k3)) with lattice vectors modulo (3, 3, 2): the 3x3x2-periodic folding of the
    longer-range ALAMODE force constants (what a 3x3x2 displacement calculation would see).
  * Units: Ry/bohr^3 -> eV/A^3 (13.605693123 / 0.529177210903^3).
Outputs: per-cluster-orbit relative difference ||(eta_A - eta_0)_o|| / ||eta_0,o|| (the eps metric of task 3),
per-shell summaries, Frobenius size, sum-rule and symmetry residuals, and fc3_alamode_folded.npy (ASR-projected,
for p3_events.py --fc3).
usage: python -B alamode_compare.py <alamode_xml> <resultdir(m553)> <outdir> [map]
map: mirror (default; z -> -z, factor (-1)^(#z)), inversion (r -> -r + c/2, factor -1), mirror_nosign (control,
deliberately wrong: positions mirrored, components not transformed).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

RY_BOHR3_TO_EV_A3 = 13.605693123 / 0.529177210903 ** 3
SHIFTS = np.array([[0, 0, 0]] + [[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)
                                 if (i, j, k) != (0, 0, 0)], float)
SC_A = np.array([4, 4, 3])           # ALAMODE supercell (checked against the lattice vectors below)
SC_P = np.array([3, 3, 2])           # phono3py fc3 supercell


def parse_alamode(path):
    txt = Path(path).read_text()
    lat = np.array([[float(x) for x in re.search(rf"<a{i}>(.*?)</a{i}>", txt).group(1).split()] for i in (1, 2, 3)])
    pos = np.array([[float(x) for x in m.group(2).split()] for m in
                    re.finditer(r'<pos index="(\d+)" element="\w+">(.*?)</pos>', txt)])
    elem = [m.group(1) for m in re.finditer(r'<pos index="\d+" element="(\w+)">', txt)]
    tran1 = [int(m.group(2)) - 1 for m in re.finditer(r'<map tran="1" atom="(\d+)">(\d+)</map>', txt)]
    blk = txt[txt.index("<ANHARM3>"):txt.index("</ANHARM3>")]
    pat = re.compile(r'<FC3 pair1="(\d+) (\d+)" pair2="(\d+) (\d+) (\d+)" pair3="(\d+) (\d+) (\d+)">(\S+)</FC3>')
    rows = np.array([[int(g) for g in m.groups()[:8]] + [0] for m in pat.finditer(blk)], dtype=np.int64)
    vals = np.array([float(m.group(9)) for m in pat.finditer(blk)])
    blk2 = txt[txt.index("<HARMONIC>"):txt.index("</HARMONIC>")]
    pat2 = re.compile(r'<FC2 pair1="(\d+) (\d+)" pair2="(\d+) (\d+) (\d+)">(\S+)</FC2>')
    fc2 = [(int(a), int(b), int(c), int(d), int(e), float(v)) for a, b, c, d, e, v in pat2.findall(blk2)]
    return lat, pos, elem, tran1, rows[:, :8], vals, fc2


def labels(frac_prim, basis_frac, tol=2e-3):
    """Integer lattice vector n and basis index k with frac = n + basis[k] (tolerance in fractional units)."""
    d = frac_prim[:, None, :] - basis_frac[None, :, :]
    r = np.round(d)
    err = np.abs(d - r).max(axis=2)
    k = err.argmin(axis=1)
    assert np.all(err[np.arange(len(k)), k] < tol), float(err.min(axis=1).max())
    return r[np.arange(len(k)), k].astype(int), k


def main():
    xml, rdir, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    mode = sys.argv[4] if len(sys.argv) > 4 else "mirror"
    Rd, tr = {"mirror": (np.array([1, 1, -1]), np.zeros(3)), "mirror_nosign": (np.array([1, 1, -1]), np.zeros(3)),
              "inversion": (np.array([-1, -1, -1]), np.array([0, 0, 0.5]))}[mode]
    Rsign = np.ones(3) if mode == "mirror_nosign" else Rd
    outdir = outdir if mode == "mirror" else outdir / mode
    outdir.mkdir(parents=True, exist_ok=True)
    out = {"xml": str(xml), "map": mode}
    lat, pos, elem, tran1, rows, vals, fc2 = parse_alamode(xml)
    out["n_atoms_alamode"], out["n_fc3_entries"] = len(pos), len(vals)
    prim_lat_A = lat / SC_A[:, None]
    out["alamode_prim_a_c_A"] = [float(np.linalg.norm(prim_lat_A[0]) * 0.529177210903),
                                 float(np.linalg.norm(prim_lat_A[2]) * 0.529177210903)]

    # --- phono3py side (hash-checked inputs) ---
    from p3_events import load_ph3
    ph3 = load_ph3(True, [5, 5, 3], 0.1)
    sc, prim = ph3.supercell, ph3.primitive
    p2s = np.asarray(prim.p2s_map)
    P_basis = np.asarray(prim.scaled_positions)
    P_frac = np.asarray(sc.scaled_positions) * SC_P            # supercell frac -> primitive frac
    nP, kP = labels(P_frac, P_basis)
    lut = {(tuple(np.mod(n, SC_P)), k): i for i, (n, k) in enumerate(zip(nP, kP))}
    out["phono3py_prim_a_c_A"] = [float(np.linalg.norm(prim.cell[0])), float(np.linalg.norm(prim.cell[2]))]

    # --- ALAMODE labels, polarity map z -> -z ---
    A_frac = pos * SC_A
    A_basis_raw = A_frac[tran1]
    A_basis_m = A_basis_raw * Rd + tr
    # basis matching ALAMODE primitive atom -> phono3py basis atom by mapped position (mod lattice)
    nb, kb = labels(A_basis_m, P_basis)
    out["alamode_to_phono3py_basis"] = kb.tolist()
    assert sorted(kb.tolist()) == [0, 1, 2, 3]
    elemP = [s for s in prim.symbols]
    assert all(elem[tran1[a]] == elemP[kb[a]] for a in range(4)), "species mismatch"
    A_frac_m = A_frac * Rd + tr

    def alabel(atom, cell):
        f = A_frac_m[atom] + (SHIFTS[cell] * SC_A) * Rd
        n, k = labels(f[None, :], P_basis)
        return n[0], k[0]

    # shortest-image check for the stored images
    cart = lambda f: f @ np.asarray(prim.cell)
    n_chk, n_short = 0, 0
    rng = np.random.default_rng(0)
    for e in rng.choice(len(vals), 4000, replace=False):
        i1 = tran1[rows[e, 0] - 1]
        for (at, cl) in ((rows[e, 2] - 1, rows[e, 4] - 1), (rows[e, 5] - 1, rows[e, 7] - 1)):
            d_st = np.linalg.norm(cart(A_frac[at] + SHIFTS[cl] * SC_A - A_frac[i1]))
            d_all = [np.linalg.norm(cart(A_frac[at] + s * SC_A - A_frac[i1])) for s in SHIFTS]
            n_chk += 1
            n_short += d_st <= min(d_all) + 1e-6
    out["stored_image_is_shortest_fraction"] = n_short / n_chk

    # --- fold into phono3py compact fc3 ---
    natP = len(sc)
    F = np.zeros((4, natP, natP, 3, 3, 3))
    keys = set()
    for e in range(len(vals)):
        keys.add(tuple(rows[e]))
    swapped_listed = sum(1 for r in keys if (r[0], r[1], r[5], r[6], r[7], r[2], r[3], r[4]) in keys
                         and (r[2], r[3], r[4]) != (r[5], r[6], r[7]))
    out["entries_with_swapped_partner_listed"] = int(swapped_listed)
    cache = {}
    for e in range(len(vals)):
        if vals[e] == 0.0:
            continue
        a1, x1, j, xj, cj, k, xk, ck = rows[e]
        a = a1 - 1
        pa = int(kb[a])
        for key in ((j - 1, cj - 1), (k - 1, ck - 1)):
            if key not in cache:
                cache[key] = alabel(*key)
        nJ, kJ = cache[(j - 1, cj - 1)]
        nK, kK = cache[(k - 1, ck - 1)]
        n0 = nb[a]                                    # lattice vector of ALAMODE primitive atom a (mapped)
        i0 = p2s[pa]
        nI = nP[i0]
        J = lut[(tuple(np.mod(nI + nJ - n0, SC_P)), kJ)]
        K = lut[(tuple(np.mod(nI + nK - n0, SC_P)), kK)]
        sgn = Rsign[x1 - 1] * Rsign[xj - 1] * Rsign[xk - 1]
        v = sgn * vals[e] * RY_BOHR3_TO_EV_A3
        pi = list(p2s).index(i0)
        F[pi, J, K, x1 - 1, xj - 1, xk - 1] += v
        if (j, xj, cj) != (k, xk, ck) and swapped_listed == 0:
            F[pi, K, J, x1 - 1, xk - 1, xj - 1] += v
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    np.save(outdir / "fc3_alamode_raw.npy", F)
    out["asr_rel_alamode_folded"] = float(np.abs(F.sum(axis=2)).max() / np.abs(F).max())
    out["asr_rel_reference"] = float(np.abs(f0.sum(axis=2)).max() / np.abs(f0).max())
    out["perm_jk_rel_alamode"] = float(np.abs(F - F.transpose(0, 2, 1, 3, 5, 4)).max() / np.abs(F).max())

    # --- projection onto the symmetric (+ASR) basis ---
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    c0 = np.load(rdir / "basis_c0.npy")
    CO = np.load(rdir / "cluster_orbits.npz")
    blk, refn, bdmax = CO["col_block"], CO["block_refnorm"], CO["block_dmax"]
    fA = F.ravel()
    etaA = 18.0 * (ccm.T @ fA)
    out["sym_resid_rel_noASR"] = float(np.linalg.norm(ccm @ etaA - fA) / np.linalg.norm(fA))
    cA = B.T @ etaA
    fA_asr = ccm @ (B @ cA)
    out["sym_ASR_resid_rel"] = float(np.linalg.norm(fA_asr - fA) / np.linalg.norm(fA))
    np.save(outdir / "fc3_alamode_folded.npy", fA_asr.reshape(f0.shape))
    eta0, etaA2 = B @ c0, B @ cA
    nb_ = len(refn)
    dn = np.sqrt(np.bincount(blk, weights=(etaA2 - eta0) ** 2, minlength=nb_))
    an = np.sqrt(np.bincount(blk, weights=etaA2 ** 2, minlength=nb_))
    rel = dn / refn
    out["frobenius_rel"] = float(np.linalg.norm(cA - c0) / np.linalg.norm(c0))
    out["norm_ratio_alamode_over_ref"] = float(np.linalg.norm(cA) / np.linalg.norm(c0))
    w2 = refn ** 2 / np.sum(refn ** 2)
    shells = [(0, 0.01), (0.01, 2.0), (2.0, 3.2), (3.2, 4.0), (4.0, 5.0), (5.0, 6.0), (6.0, 20.0)]
    out["shells"] = []
    for lo, hi in shells:
        m = (bdmax >= lo) & (bdmax < hi)
        if not m.any():
            continue
        out["shells"].append({"shell_A": [lo, hi], "n_orbits": int(m.sum()), "ref_norm2_frac": float(w2[m].sum()),
                              "rel_change_shell": float(np.sqrt(np.sum(dn[m] ** 2) / np.sum(refn[m] ** 2))),
                              "orbit_rel_median": float(np.median(rel[m])), "orbit_rel_max": float(rel[m].max()),
                              "norm_ratio_shell": float(np.sqrt(np.sum(an[m] ** 2) / np.sum(refn[m] ** 2)))})
    order = np.argsort(-refn)
    out["top_orbits"] = [{"orbit": int(o), "dmax_A": float(bdmax[o]), "ref_norm": float(refn[o]),
                          "ref_norm2_frac": float(w2[o]), "rel_change": float(rel[o]),
                          "norm_ratio": float(an[o] / refn[o])} for o in order[:12]]
    for q in (0.5, 0.9, 0.99):
        out[f"orbit_rel_quantile_{q}"] = float(np.quantile(rel, q))
    out["orbit_rel_weighted_rms"] = float(np.sqrt(np.sum(w2 * rel ** 2)))
    np.savez(outdir / "alamode_orbits.npz", rel=rel, refnorm=refn, alanorm=an, dmax=bdmax, cA=cA)
    (outdir / "alamode_compare.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k not in ("shells", "top_orbits")}, indent=1))
    for s in out["shells"]:
        print(s)
    for t in out["top_orbits"]:
        print(t)


if __name__ == "__main__":
    main()
