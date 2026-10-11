"""Review R1: anatomy of the per-orbit ALAMODE-vs-reference differences delta_o (input of the calibrated class).

Questions: (a) on how many orbits is the folded ALAMODE tensor identically zero (delta_o = 1 by construction, a
cutoff property of one workflow rather than a measured disagreement)? (b) per-orbit cosine between eta_A and
eta_0 where both are nonzero (correlated physics or noise?); (c) consistency of the 'shell' label (largest of
three independent shortest-image distances) with a single geometric embedding of the triplet; (d) share of the
calibrated-class box (k = 1) that sits on cutoff-zero orbits.
Research env. usage: python -B r1_orbit_anatomy.py <rdir m553> <alamode_dir> <outjson>
"""
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

rdir, adir, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
CO = np.load(rdir / "cluster_orbits.npz")
blk, refn, bdmax, cdmax = CO["col_block"], CO["block_refnorm"], CO["block_dmax"], CO["col_dmax"]
A = np.load(adir / "alamode_orbits.npz")
rel, an, cA = A["rel"], A["alanorm"], A["cA"]
B = np.load(rdir / "basis_B.npy", mmap_mode="r")
c0 = np.load(rdir / "basis_c0.npy")
eta0 = B @ c0
etaA = B @ cA
nb = len(refn)
dot = np.bincount(blk, weights=eta0 * etaA, minlength=nb)
n0 = np.sqrt(np.bincount(blk, weights=eta0 ** 2, minlength=nb))
nA = np.sqrt(np.bincount(blk, weights=etaA ** 2, minlength=nb))
cos = np.where((n0 > 0) & (nA > 0), dot / np.maximum(n0 * nA, 1e-300), np.nan)
ncol = np.bincount(blk, minlength=nb)
zeroA = nA <= 1e-9 * n0
out = {"n_orbits": int(nb), "n_alamode_zero_orbits": int(zeroA.sum())}
edges = [0, 0.01, 2.0, 3.2, 4.0, 5.0, 6.0, 20.0]
w2 = refn ** 2 / np.sum(refn ** 2)
out["shells"] = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (bdmax >= lo) & (bdmax < hi)
    mz = m & zeroA
    mnz = m & ~zeroA
    c = cos[mnz]
    out["shells"].append({
        "shell_A": [lo, hi], "n_orbits": int(m.sum()), "n_columns": int(ncol[m].sum()),
        "n_alamode_zero": int(mz.sum()), "ref_norm2_frac_on_zero_orbits": float(w2[mz].sum()),
        "ref_norm2_frac_shell": float(w2[m].sum()),
        "cos_nonzero_median": float(np.median(c)) if len(c) else None,
        "cos_nonzero_quantiles_10_50_90": np.quantile(c, [0.1, 0.5, 0.9]).tolist() if len(c) else None,
        "frac_cos_negative": float(np.mean(c < 0)) if len(c) else None,
        "rel_nonzero_median": float(np.median(rel[mnz])) if mnz.any() else None,
        "normratio_nonzero_median": float(np.median(nA[mnz] / n0[mnz])) if mnz.any() else None})
# orbits below 3.2 A with delta >= 0.5
low = np.nonzero((bdmax < 3.2) & (rel >= 0.5))[0]
out["short_orbits_delta_ge_0.5"] = [{"orbit": int(o), "dmax": float(bdmax[o]), "refnorm": float(refn[o]),
                                     "delta": float(rel[o]), "alanorm": float(nA[o]), "cos": None if np.isnan(cos[o]) else float(cos[o]),
                                     "ncol": int(ncol[o])} for o in low]
big = np.argsort(-rel)[:10]
out["largest_delta_orbits"] = [{"orbit": int(o), "dmax": float(bdmax[o]), "refnorm": float(refn[o]),
                                "refnorm_rel_to_max": float(refn[o] / refn.max()), "delta": float(rel[o]),
                                "alanorm": float(nA[o]), "cos": None if np.isnan(cos[o]) else float(cos[o])} for o in big]
# calibrated box (k = 1): absolute radius per orbit = max(delta, 0.005) * refnorm; share of box "volume"
# (sum of squared radii) on cutoff-zero orbits and per shell
rad = np.maximum(rel, 0.005) * refn
out["k1_box_radius2_share_on_zero_orbits"] = float(np.sum(rad[zeroA] ** 2) / np.sum(rad ** 2))
out["k1_box_radius2_share_by_shell"] = {f"{lo}-{hi}": float(np.sum(rad[(bdmax >= lo) & (bdmax < hi)] ** 2) / np.sum(rad ** 2))
                                        for lo, hi in zip(edges[:-1], edges[1:])}
np.savez(outp.with_suffix(".npz"), cos=cos, zeroA=zeroA, n0=n0, nA=nA, rel=rel, bdmax=bdmax, refn=refn)
outp.write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
