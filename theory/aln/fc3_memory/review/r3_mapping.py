"""Review R3: independent re-derivation of the ALAMODE -> phono3py fc3 folding (no code shared with
alamode_compare.py): ElementTree parsing, Cartesian displacement vectors, mirror applied to vectors, nearest-atom
matching in the phono3py 3x3x2 supercell, per-axis sign product. Also: polarity of both structures, image
multiplicity statistics of the ALAMODE FC3 entries, physical cluster diameters of the ALAMODE clusters
(cutoff), and consistent (single-embedding) diameters of the 456 phono3py cluster orbits versus the
independent-shortest-image d_max used for the shell labels.
Research env. usage: python -B r3_mapping.py <xml> <rdir m553> <outdir>
"""
import itertools
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import scipy.sparse as sp

BOHR = 0.529177210903
RY = 13.605693123
xmlp, rdir, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
outdir.mkdir(parents=True, exist_ok=True)
out = {}

# ---------------- ALAMODE structure (ElementTree)
root = ET.parse(xmlp).getroot()
latA = np.array([[float(x) for x in root.find(f"Structure/LatticeVector/a{i}").text.split()] for i in (1, 2, 3)]) * BOHR
posA, elA = [], []
for p in root.find("Structure/Position"):
    posA.append([float(x) for x in p.text.split()])
    elA.append(p.attrib["element"])
posA = np.array(posA)
cartA = posA @ latA
nA = len(posA)
tran1 = {}
for mp in root.find("Symmetry/Translations"):
    if mp.attrib["tran"] == "1":
        tran1[int(mp.attrib["atom"])] = int(mp.text) - 1
# image ordering of ALM: index 0 = (0,0,0), then ia, ja, ka in -1..1 (ia slowest), (0,0,0) skipped
imgs = [(0, 0, 0)] + [t for t in itertools.product((-1, 0, 1), repeat=3) if t != (0, 0, 0)]
imgs = np.array(imgs, float)

# polarity: for each Al, the N with the same in-plane position at |dz| < 0.45 c (axial bond)
def axial(cart, el, lat, c):
    res = []
    for i in range(len(cart)):
        if el[i] != "Al":
            continue
        best = None
        for j in range(len(cart)):
            if el[j] != "N":
                continue
            for s in itertools.product((-1, 0, 1), repeat=3):
                d = cart[j] + np.array(s) @ lat - cart[i]
                if np.hypot(d[0], d[1]) < 1e-3 and abs(d[2]) < 0.45 * c:
                    best = d[2]
        res.append(best)
    return np.array(res, float)

cA = np.linalg.norm(latA[2]) / 3
axA = axial(cartA, elA, latA, cA)
S = np.load(rdir / "static_r01.npz")
latP, scaP, p2s = S["sc_cell"], S["sc_scaled"], S["p2s"]
primP = S["prim_cell"]
cartP = scaP @ latP
elP = ["Al" if m > 20 else "N" for m in np.repeat(S["masses"], 1)[S["s2p"]]] if False else None
numsP = S["sc_numbers"]
elP = ["Al" if n == 13 else "N" for n in numsP]
cP = np.linalg.norm(primP[2])
axP = axial(cartP, elP, latP, cP)
out["axial_bond_dz_A"] = {"alamode_min_max": [float(np.nanmin(axA)), float(np.nanmax(axA))],
                          "phono3py_min_max": [float(np.nanmin(axP)), float(np.nanmax(axP))]}

# ---------------- FC3 entries
primA = latA / np.array([4, 4, 3])[:, None]
inv_primA = np.linalg.inv(primA)
R = np.diag([1.0, 1.0, -1.0])
rows = []
vals = []
for fc in root.find("ForceConstants/ANHARM3"):
    a1 = fc.attrib["pair1"].split()
    a2 = fc.attrib["pair2"].split()
    a3 = fc.attrib["pair3"].split()
    rows.append([int(a1[0]), int(a1[1]), int(a2[0]), int(a2[1]), int(a2[2]), int(a3[0]), int(a3[1]), int(a3[2])])
    vals.append(float(fc.text))
rows = np.array(rows)
vals = np.array(vals) * RY / BOHR ** 3
out["n_entries"] = len(vals)
i_sc = np.array([tran1[r] for r in rows[:, 0]])
rJ = cartA[rows[:, 2] - 1] + imgs[rows[:, 4] - 1] @ latA
rK = cartA[rows[:, 5] - 1] + imgs[rows[:, 7] - 1] @ latA
ri = cartA[i_sc]
dJ, dK = rJ - ri, rK - ri
dJK = rK - rJ
diam = np.maximum(np.maximum(np.linalg.norm(dJ, axis=1), np.linalg.norm(dK, axis=1)), np.linalg.norm(dJK, axis=1))
nz = vals != 0
out["alamode_cluster_diameter_max_A"] = float(diam[nz].max())
out["alamode_pair_distance_max_A"] = float(max(np.linalg.norm(dJ[nz], axis=1).max(), np.linalg.norm(dK[nz], axis=1).max()))
out["alamode_diameter_quantiles_nonzero"] = np.quantile(diam[nz], [0.5, 0.9, 0.99, 1]).tolist()
# image multiplicity: same (pair1, atomJ, xyzJ, atomK, xyzK) listed with different cells
key = {}
for e, r in enumerate(rows):
    key.setdefault((r[0], r[1], r[2], r[3], r[5], r[6]), []).append(e)
mult = np.array([len(v) for v in key.values()])
out["image_multiplicity_counts"] = {int(m): int((mult == m).sum()) for m in np.unique(mult)}
# for multi-image groups: are the values equal across images (value divided by multiplicity)?
eq, neq = 0, 0
for v in key.values():
    if len(v) > 1:
        vv = vals[v]
        if np.allclose(vv, vv[0], rtol=1e-10, atol=1e-14):
            eq += 1
        else:
            neq += 1
out["multi_image_groups_equal_values"] = [eq, neq]

# ---------------- map displacement vectors (mirror, ALAMODE fractional -> phono3py fractional)
def to_P_frac(dcart):
    fA = dcart @ inv_primA            # ALAMODE primitive fractional
    fA = fA * np.array([1, 1, -1])    # mirror z -> -z acts on the c fractional component only
    return fA                          # same fractional in phono3py's primitive lattice

fJ, fK = to_P_frac(dJ), to_P_frac(dK)
# primitive atom of pair1 (ALAMODE primitive index a -> basis position, mirrored) -> phono3py primitive atom
basisA = (cartA[[tran1[a] for a in (1, 2, 3, 4)]] @ inv_primA) * np.array([1, 1, -1])
basisP = scaP[p2s] * np.array([3, 3, 2])      # phono3py primitive fractional of p2s atoms
amap = []
for a in range(4):
    d = basisP - basisA[a]
    err = np.abs(d - np.round(d)).max(axis=1)
    amap.append(int(err.argmin()))
    assert err.min() < 2e-3, err
out["alamode_prim_to_phono3py_prim"] = amap
pa = np.array([amap[r - 1] for r in rows[:, 0]])
scP_frac_prim = scaP * np.array([3, 3, 2])     # all supercell atoms in primitive fractional units

def find_atom(fr):
    """index of the phono3py supercell atom at primitive-fractional position fr (mod 3,3,2)."""
    d = fr[:, None, :] - scP_frac_prim[None, :, :]
    d = d - np.round(d / np.array([3, 3, 2])) * np.array([3, 3, 2])
    err = np.abs(d).max(axis=2)
    j = err.argmin(axis=1)
    return j, err[np.arange(len(j)), j]

base = scP_frac_prim[p2s[pa]]
J, eJ = find_atom(base + fJ)
K, eK = find_atom(base + fK)
out["match_err_max_frac"] = float(max(eJ.max(), eK.max()))
spec_ok = all(elP[j] == elA[rows[e, 2] - 1] for e, j in enumerate(J[:5000]))
out["species_check_first5000"] = bool(spec_ok)
sgn = np.ones(3)
sgn[2] = -1
s = sgn[rows[:, 1] - 1] * sgn[rows[:, 3] - 1] * sgn[rows[:, 6] - 1]
F = np.zeros((4, 72, 72, 3, 3, 3))
np.add.at(F, (pa, J, K, rows[:, 1] - 1, rows[:, 3] - 1, rows[:, 6] - 1), s * vals)
offd = ~((rows[:, 2] == rows[:, 5]) & (rows[:, 3] == rows[:, 6]) & (rows[:, 4] == rows[:, 7]))
np.add.at(F, (pa[offd], K[offd], J[offd], rows[offd, 1] - 1, rows[offd, 6] - 1, rows[offd, 3] - 1), s[offd] * vals[offd])
np.save(outdir / "r3_fc3_alamode_raw_independent.npy", F)
Fc = np.load(rdir / "alamode" / "fc3_alamode_raw.npy")
out["independent_vs_campaign_raw_maxabs_rel"] = float(np.abs(F - Fc).max() / np.abs(Fc).max())
out["asr_rel_independent"] = float(np.abs(F.sum(axis=2)).max() / np.abs(F).max())
# a few spot values
f0 = np.load(rdir / "fc3_ref_compact.npy")
spots = []
for (p, j, k, a, b, c) in [(0, 0, 0, 2, 2, 2), (0, 0, 0, 0, 0, 2), (2, 2, 2, 2, 2, 2)]:
    spots.append({"index": [p, j, k, a, b, c], "alamode": float(F[p, j, k, a, b, c]), "reference": float(f0[p, j, k, a, b, c])})
out["spot_values_onsite"] = spots

# ---------------- consistent diameters of the phono3py cluster orbits
CO = np.load(rdir / "cluster_orbits.npz")
blk = CO["col_block"]
bdmax = CO["block_dmax"]
ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsc()
nb = int(blk.max()) + 1
rep = {}
for col in range(ccm.shape[1]):
    b = blk[col]
    if b in rep:
        continue
    el = ccm.indices[ccm.indptr[col]:ccm.indptr[col + 1]][0]
    t = el // 27
    rep[b] = (t // (72 * 72), (t // 72) % 72, t % 72)
shifts = np.array(list(itertools.product(range(-2, 3), repeat=3)), float) @ latP
cons, indep = np.zeros(nb), np.zeros(nb)
for b, (p, j, k) in rep.items():
    ri = cartP[p2s[p]]
    rj = cartP[j] + shifts
    rk = cartP[k] + shifts
    dij = np.linalg.norm(rj - ri, axis=1)
    dik = np.linalg.norm(rk - ri, axis=1)
    djk = np.linalg.norm(rk[None, :, :] - rj[:, None, :], axis=2)
    dm = np.maximum(np.maximum(dij[:, None], dik[None, :]), djk)
    cons[b] = dm.min()
    indep[b] = max(dij.min(), dik.min(), djk.min())
alias = cons > indep + 1e-6
A = np.load(rdir / "alamode" / "alamode_orbits.npz")
zero = A["alanorm"] <= 1e-9 * A["refnorm"]
out["orbits_consistent_diameter_exceeds_shell_label"] = int(alias.sum())
out["aliased_orbits"] = [{"orbit": int(b), "dmax_label": float(bdmax[b]), "indep_check": float(indep[b]),
                          "consistent_diameter": float(cons[b]), "delta": float(A["rel"][b]),
                          "alamode_zero": bool(zero[b]), "ref_norm": float(A["refnorm"][b])}
                         for b in np.nonzero(alias)[0] if bdmax[b] < 5.0]
edges = [0, 0.01, 2.0, 3.2, 4.0, 5.0, 6.0, 20.0]
out["aliased_by_label_shell"] = {f"{lo}-{hi}": [int(((bdmax >= lo) & (bdmax < hi) & alias).sum()),
                                                 int(((bdmax >= lo) & (bdmax < hi)).sum())]
                                 for lo, hi in zip(edges[:-1], edges[1:])}
out["zero_orbits_consistent_diameter_min"] = float(cons[zero].min())
out["nonzero_orbits_consistent_diameter_max"] = float(cons[~zero].max())
np.savez(outdir / "r3_orbit_diameters.npz", cons=cons, indep=indep, label=bdmax, zero=zero)
(outdir / "r3_mapping.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
