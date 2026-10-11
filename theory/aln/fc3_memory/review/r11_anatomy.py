"""Review R11: anatomy of witness fc3 changes by (a) the campaign's shell label (largest of three independent
shortest-image distances) and (b) the consistent single-embedding cluster diameter (r3_mapping), plus the share of
the change carried by orbits beyond ALAMODE's cubic cutoff (consistent diameter > 5.68 A, where delta_o = 1 by
construction). Change measured on the rebuilt tensor (projection 18 ccm^T df).
usage: python -B r11_anatomy.py <rdir> <diam.npz> <outjson> name=fc3.npy ...
"""
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

rdir, diamf, outp = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
f0 = np.load(rdir / "fc3_ref_compact.npy")
CO = np.load(rdir / "cluster_orbits.npz")
blk, refn = CO["col_block"], CO["block_refnorm"]
D = np.load(diamf)
cons, label, zero = D["cons"], D["label"], D["zero"]
edges = [0, 0.01, 2.0, 3.2, 4.0, 5.0, 6.0, 20.0]
res = {}
for spec in sys.argv[4:]:
    name, path = spec.split("=", 1)
    df = np.load(path) - f0
    deta = 18.0 * (ccm.T @ df.ravel())
    dn2 = np.bincount(blk, weights=deta ** 2, minlength=len(refn))
    rel = np.sqrt(dn2) / refn
    rec = {"share_change2_beyond_alamode_cutoff": float(dn2[zero].sum() / dn2.sum()),
           "share_change2_aliased_orbits": float(dn2[cons > label + 1e-6].sum() / dn2.sum())}
    for nm, d in (("label", label), ("consistent", cons)):
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (d >= lo) & (d < hi)
            if m.any():
                rows.append({"shell": [lo, hi], "n": int(m.sum()), "rel_change_shell": float(np.sqrt(dn2[m].sum() / (refn[m] ** 2).sum())),
                             "share_of_change2": float(dn2[m].sum() / dn2.sum())})
        rec[f"by_{nm}"] = rows
    m23 = (label >= 2.0) & (label < 3.2)
    rec["label_2_3.2_shell_change2_share_from_aliased"] = float(dn2[m23 & (cons > label + 1e-6)].sum() / max(dn2[m23].sum(), 1e-300))
    res[name] = rec
    print(name, json.dumps({k: v for k, v in rec.items() if not k.startswith("by_")}))
    for r in rec["by_consistent"]:
        print("   consistent", r)
outp.write_text(json.dumps(res, indent=1))
