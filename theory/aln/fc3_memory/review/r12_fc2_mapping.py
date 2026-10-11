"""Review R12: mapping test on a second tensor order. ALAMODE harmonic FC2 (harmonic443.xml, 4x4x3) mapped with
the same conventions as the fc3 comparison (pair1 = primitive atom, ALM image ordering, polarity mirror z -> -z
with factor (-1)^(#z), Ry/bohr^2 -> eV/A^2) onto the phono3py fc2 (5x5x3 supercell, compact (4, 300, 3, 3)) and
compared pair by pair for pairs shorter than 4.5 A (unique shortest images in both supercells). A wrong polarity or
orientation convention would show up as sign errors on the nearest-neighbour blocks.
Research env. usage: python -B r12_fc2_mapping.py <harmonic.xml> <fc2.hdf5> <static_r01.npz> <outjson>
"""
import itertools
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np

BOHR, RY = 0.529177210903, 13.605693123
xmlp, fc2p, statp, outp = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
root = ET.parse(xmlp).getroot()
latA = np.array([[float(x) for x in root.find(f"Structure/LatticeVector/a{i}").text.split()] for i in (1, 2, 3)]) * BOHR
posA = np.array([[float(x) for x in p.text.split()] for p in root.find("Structure/Position")])
cartA = posA @ latA
tran1 = {int(mp.attrib["atom"]): int(mp.text) - 1 for mp in root.find("Symmetry/Translations") if mp.attrib["tran"] == "1"}
imgs = np.array([(0, 0, 0)] + [t for t in itertools.product((-1, 0, 1), repeat=3) if t != (0, 0, 0)], float)
primA = latA / np.array([4, 4, 3])[:, None]
invA = np.linalg.inv(primA)
S = np.load(statp)
primP = S["prim_cell"]
basisP = S["prim_scaled"]
with h5py.File(fc2p, "r") as h:
    fc2 = h["force_constants"][:]
    p2s2 = h["p2s_map"][:]
# phono3py 5x5x3 phonon supercell positions (primitive fractional): atom index = (cell, basis) in phonopy order
# rebuild by matching: generate all n in 5x5x3 and basis k; phonopy orders supercell atoms basis-major
SC = np.array([5, 5, 3])
cells = np.array([(i, j, k) for k in range(SC[2]) for j in range(SC[1]) for i in range(SC[0])])
# phonopy supercell ordering: for each primitive atom, all lattice points (i fastest)
scfrac = np.array([cells[c] + basisP[b] for b in range(4) for c in range(len(cells))])
assert len(scfrac) == 300
pmap = [0, 75, 150, 225]
assert list(p2s2) == pmap, p2s2
rows, vals = [], []
for fc in root.find("ForceConstants/HARMONIC"):
    a1 = fc.attrib["pair1"].split()
    a2 = fc.attrib["pair2"].split()
    rows.append([int(a1[0]), int(a1[1]), int(a2[0]), int(a2[1]), int(a2[2])])
    vals.append(float(fc.text))
rows = np.array(rows)
vals = np.array(vals) * RY / BOHR ** 2
ri = cartA[[tran1[a] for a in rows[:, 0]]]
rj = cartA[rows[:, 2] - 1] + imgs[rows[:, 4] - 1] @ latA
dvec = rj - ri
dist = np.linalg.norm(dvec, axis=1)
fr = (dvec @ invA) * np.array([1, 1, -1])           # mirrored, primitive fractional (same in phono3py lattice)
basisA = (cartA[[tran1[a] for a in (1, 2, 3, 4)]] @ invA) * np.array([1, 1, -1])
amap = []
for a in range(4):
    d = basisP - basisA[a]
    err = np.abs(d - np.round(d)).max(axis=1)
    amap.append(int(err.argmin()))
pa = np.array([amap[r - 1] for r in rows[:, 0]])
target = basisP[pa] + fr
d = target[:, None, :] - scfrac[None, :, :]
d = d - np.round(d / SC) * SC
err = np.abs(d).max(axis=2)
J = err.argmin(axis=1)
sel = dist < 4.5
sgn = np.array([1.0, 1.0, -1.0])
A2 = np.zeros((4, 300, 3, 3))
cnt = np.zeros((4, 300, 3, 3))
for e in np.nonzero(sel)[0]:
    s = sgn[rows[e, 1] - 1] * sgn[rows[e, 3] - 1]
    A2[pa[e], J[e], rows[e, 1] - 1, rows[e, 3] - 1] += s * vals[e]
    cnt[pa[e], J[e], rows[e, 1] - 1, rows[e, 3] - 1] += 1
mask = cnt > 0
out = {"n_pairs_lt_4.5A": int(sel.sum()), "match_err_max": float(err[np.arange(len(J)), J][sel].max()),
       "fc2_rel_diff_on_mapped_pairs": float(np.linalg.norm((A2 - fc2)[mask]) / np.linalg.norm(fc2[mask])),
       "onsite_blocks": {"alamode": A2[0, 0].tolist(), "phono3py": fc2[0, 0].tolist()}}
# nearest neighbour blocks (axial bond Al-N and basal)
blocks = []
for e in np.nonzero(sel & (dist > 1.0) & (dist < 2.0))[0][:4]:
    p, j = pa[e], J[e]
    blocks.append({"dist": float(dist[e]), "alamode": A2[p, j].round(4).tolist(), "phono3py": fc2[p, j].round(4).tolist()})
out["nn_blocks"] = blocks
# the same with the polarity map omitted (control): positions not mirrored
fr2 = dvec @ invA
basisA2 = cartA[[tran1[a] for a in (1, 2, 3, 4)]] @ invA
amap2 = []
for a in range(4):
    dd = basisP - basisA2[a]
    er = np.abs(dd - np.round(dd)).max(axis=1)
    amap2.append((int(er.argmin()), float(er.min())))
out["control_no_mirror_basis_match"] = amap2
outp.write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
