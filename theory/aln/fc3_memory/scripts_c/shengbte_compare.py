"""E28: third first-principles workflow — the ShengBTE / thirdorder.py fc3 of the Phonon Olympics benchmark (Quantum
ESPRESSO, PBEsol PAW) folded into the 3x3x2 supercell of this campaign and projected onto the same symmetric + ASR
basis, as alamode_compare.py does for ALAMODE.

Mapping (position based, as review/r3_mapping.py): each FORCE_CONSTANTS_3RD block gives atoms i, j, k of the ShengBTE
unit cell and the Cartesian cell vectors R2, R3 (A) of j and k; displacement vectors d = R + tau - tau_i are written
in ShengBTE primitive fractional coordinates, mirrored z -> -z (the ShengBTE structure has the opposite polarity:
axial Al-N bond along +c, phono3py along -c), read as phono3py primitive fractional coordinates and matched to the
nearest phono3py supercell atom (mod 3,3,2). Component sign: product of per-axis factors (+1, +1, -1). Images that
fold onto the same supercell triplet are summed (what a 3x3x2 finite-displacement calculation would see).
Inputs: FORCE_CONSTANTS_3RD and CONTROL (structure) at commit 0640f077, checksums in the dataset folder.
usage: python -B shengbte_compare.py <FORCE_CONSTANTS_3RD> <CONTROL> <rdir m553> <outdir> [mirror|nomirror]
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

fcp, ctrl, rdir, outdir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
mode = sys.argv[5] if len(sys.argv) > 5 else "mirror"
outdir.mkdir(parents=True, exist_ok=True)
out = {"mode": mode}
BOHR = 0.529177210903

# ---------------- ShengBTE structure from CONTROL
txt = ctrl.read_text()
lfac = float(re.search(r"lfactor\s*=\s*([0-9.eEdD+-]+)", txt).group(1)) * 10.0      # nm -> A
lat = np.array([[float(x) for x in re.search(rf"lattvec\(:,{i}\)\s*=\s*([^\n]+)", txt).group(1).split()[:3]]
                for i in (1, 2, 3)]) * lfac
pos = np.array([[float(x) for x in re.search(rf"positions\(:,{i}\)\s*=\s*([^\n]+)", txt).group(1).split()[:3]]
                for i in (1, 2, 3, 4)])
types = [int(x) for x in re.search(r"types\s*=\s*([0-9 ]+)", txt).group(1).split()]
elS = ["Al" if t == 1 else "N" for t in types]
cartS = pos @ lat
inv_lat = np.linalg.inv(lat)
out["shengbte_lattice_A"] = lat.tolist()

# ---------------- FORCE_CONSTANTS_3RD
lines = [l.split() for l in fcp.read_text().splitlines() if l.strip()]
nblk = int(lines[0][0])
p = 1
I, J3, K3, R2, R3, V = [], [], [], [], [], []
for _ in range(nblk):
    p += 1                                   # block index
    r2 = [float(x) for x in lines[p]]; p += 1
    r3 = [float(x) for x in lines[p]]; p += 1
    i, j, k = (int(x) for x in lines[p]); p += 1
    vals = np.zeros((3, 3, 3))
    for _ in range(27):
        a, b, c, v = lines[p]; p += 1
        vals[int(a) - 1, int(b) - 1, int(c) - 1] = float(v)
    I.append(i - 1); J3.append(j - 1); K3.append(k - 1); R2.append(r2); R3.append(r3); V.append(vals)
I, J3, K3, R2, R3, V = map(np.array, (I, J3, K3, R2, R3, V))
out["n_blocks"] = int(nblk)
dJ = R2 + cartS[J3] - cartS[I]
dK = R3 + cartS[K3] - cartS[I]
dJK = dK - dJ
diam = np.maximum.reduce([np.linalg.norm(dJ, axis=1), np.linalg.norm(dK, axis=1), np.linalg.norm(dJK, axis=1)])
nz = np.abs(V).reshape(len(V), -1).max(axis=1) > 0
out["cluster_diameter_max_A"] = float(diam[nz].max())
out["pair_distance_max_A"] = float(max(np.linalg.norm(dJ[nz], axis=1).max(), np.linalg.norm(dK[nz], axis=1).max()))

# ---------------- phono3py supercell
S = np.load(rdir / "static_r01.npz")
latP, scaP, p2s, numsP = S["sc_cell"], S["sc_scaled"], S["p2s"], S["sc_numbers"]
elP = ["Al" if n == 13 else "N" for n in numsP]
SC = np.array([3, 3, 2])
scP_prim = scaP * SC
flip = np.array([1.0, 1.0, -1.0]) if mode == "mirror" else np.ones(3)


def frac(d):
    return (d @ inv_lat) * flip


basisS = (cartS @ inv_lat) * flip
basisP = scP_prim[p2s]
amap = []
for a in range(4):
    d = basisP - basisS[a]
    err = np.abs(d - np.round(d)).max(axis=1)
    amap.append(int(err.argmin()))
    out.setdefault("basis_match_err", []).append(float(err.min()))
assert max(out["basis_match_err"]) < 2e-3, out["basis_match_err"]
assert all(elP[p2s[amap[a]]] == elS[a] for a in range(4)), "species mismatch in basis map"
out["shengbte_to_phono3py_prim"] = amap


def find_atom(fr):
    d = fr[:, None, :] - scP_prim[None, :, :]
    d = d - np.round(d / SC) * SC
    err = np.abs(d).max(axis=2)
    j = err.argmin(axis=1)
    return j, err[np.arange(len(j)), j]


pa = np.array([amap[i] for i in I])
base = scP_prim[p2s[pa]]
Jp, eJ = find_atom(base + frac(dJ))
Kp, eK = find_atom(base + frac(dK))
out["match_err_max_frac"] = float(max(eJ.max(), eK.max()))
assert out["match_err_max_frac"] < 5e-3, out["match_err_max_frac"]
assert all(elP[Jp[e]] == elS[J3[e]] and elP[Kp[e]] == elS[K3[e]] for e in range(len(I))), "species mismatch"
sg = flip
F = np.zeros((4, 72, 72, 3, 3, 3))
for e in range(len(I)):
    s3 = V[e] * sg[:, None, None] * sg[None, :, None] * sg[None, None, :]
    F[pa[e], Jp[e], Kp[e]] += s3
f0 = np.load(rdir / "fc3_ref_compact.npy")
np.save(outdir / "fc3_shengbte_raw.npy", F)
out["asr_rel_shengbte_folded"] = float(np.abs(F.sum(axis=2)).max() / np.abs(F).max())
out["perm_jk_rel"] = float(np.abs(F - F.transpose(0, 2, 1, 3, 5, 4)).max() / np.abs(F).max())

# ---------------- projection onto the symmetric (+ASR) basis and orbit statistics (as alamode_compare.py)
ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
B = np.load(rdir / "basis_B.npy")
c0 = np.load(rdir / "basis_c0.npy")
CO = np.load(rdir / "cluster_orbits.npz")
blk, refn, bdmax = CO["col_block"], CO["block_refnorm"], CO["block_dmax"]
fS = F.ravel()
eta = 18.0 * (ccm.T @ fS)
out["sym_resid_rel_noASR"] = float(np.linalg.norm(ccm @ eta - fS) / np.linalg.norm(fS))
cS = B.T @ eta
fS_asr = ccm @ (B @ cS)
out["sym_ASR_resid_rel"] = float(np.linalg.norm(fS_asr - fS) / np.linalg.norm(fS))
np.save(outdir / "fc3_shengbte_folded.npy", fS_asr.reshape(f0.shape))
eta0, etaS = B @ c0, B @ cS
nb_ = len(refn)
dn = np.sqrt(np.bincount(blk, weights=(etaS - eta0) ** 2, minlength=nb_))
an = np.sqrt(np.bincount(blk, weights=etaS ** 2, minlength=nb_))
rel = dn / refn
out["frobenius_rel_vs_phono3py"] = float(np.linalg.norm(cS - c0) / np.linalg.norm(c0))
out["norm_ratio_over_phono3py"] = float(np.linalg.norm(cS) / np.linalg.norm(c0))
A = np.load(rdir / "alamode" / "alamode_orbits.npz")
cA = A["cA"]
etaA = B @ cA
dA = np.sqrt(np.bincount(blk, weights=(etaS - etaA) ** 2, minlength=nb_))
out["frobenius_rel_vs_alamode"] = float(np.linalg.norm(cS - cA) / np.linalg.norm(c0))
out["n_orbits_zero_in_shengbte"] = int((an <= 1e-9 * refn).sum())
w2 = refn ** 2 / np.sum(refn ** 2)
shells = [(0, 0.01), (0.01, 2.0), (2.0, 3.2), (3.2, 4.0), (4.0, 5.0), (5.0, 6.0), (6.0, 20.0)]
out["shells"] = []
for lo, hi in shells:
    m = (bdmax >= lo) & (bdmax < hi)
    if m.any():
        out["shells"].append({"shell_A": [lo, hi], "n_orbits": int(m.sum()), "ref_norm2_frac": float(w2[m].sum()),
                              "rel_vs_phono3py": float(np.sqrt(np.sum(dn[m] ** 2) / np.sum(refn[m] ** 2))),
                              "rel_vs_alamode": float(np.sqrt(np.sum(dA[m] ** 2) / np.sum(refn[m] ** 2))),
                              "alamode_vs_phono3py": float(np.sqrt(np.sum((A["rel"][m] * refn[m]) ** 2) / np.sum(refn[m] ** 2))),
                              "n_zero_shengbte": int((an[m] <= 1e-9 * refn[m]).sum())})
np.savez(outdir / "shengbte_orbits.npz", rel=rel, refnorm=refn, sbtnorm=an, dmax=bdmax, cS=cS, rel_vs_alamode=dA / refn)
(outdir / "shengbte_compare.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: v for k, v in out.items() if k != "shells"}, indent=1))
for s in out["shells"]:
    print(s)
