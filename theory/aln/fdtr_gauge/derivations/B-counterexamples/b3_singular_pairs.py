"""B3: where is det M2(theta; w1, w2) = 0?  Sign map on a grid of (w1, w2), w1 < w2, in 30-digit arithmetic.

Off-diagonal sign changes of det M2 locate the singular curves (non-generic pairs).  For each theta, report
the number of grid cells (w1 < w2) adjacent to a sign change, and one refined singular pair (root of det along
w2 at fixed w1) with the singular values there.
"""
import sys, os, json
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from b3_two_freq import M2_mp

mp.mp.dps = 30
out = []
for th in [(1.0, 2.0, 0.5, 0.1), (1.0, 0.5, 0.1, 0.5), (1.0, 1.0, 0.5, 0.1), (1.0, 2.0, 0.5, 0.0), (1.0, 10.0, 3.0, 1.0)]:
    lg = np.linspace(np.log(1e-2), np.log(1e2), 33)
    D = np.full((len(lg), len(lg)), np.nan)
    for i in range(len(lg)):
        for j in range(i + 1, len(lg)):
            D[i, j] = float(mp.sign(mp.det(M2_mp(th, mp.exp(lg[i]), mp.exp(lg[j])))))
    changes = 0
    first = None
    for i in range(len(lg)):
        for j in range(i + 2, len(lg)):
            if D[i, j - 1] * D[i, j] < 0:
                changes += 1
                if first is None:
                    first = (i, j)
    rec = {"theta": th, "n_sign_changes_along_w2": changes}
    if first is not None:
        i, j = first
        w1 = mp.exp(lg[i])
        f = lambda l2: mp.det(M2_mp(th, w1, mp.exp(l2)))
        rt = mp.findroot(f, (mp.mpf(lg[j - 1]), mp.mpf(lg[j])), solver="anderson")
        M = M2_mp(th, w1, mp.exp(rt))
        sv = sorted([abs(x) for x in mp.svd_r(M, compute_uv=False)], reverse=True)
        rec.update({"w1": float(w1), "w2_root": float(mp.exp(rt)), "s_over_s1_at_root": [float(x / sv[0]) for x in sv]})
    out.append(rec)
    print(rec)
with open(os.path.join(HERE, "data", "b3_singular_pairs.json"), "w") as fh:
    json.dump(out, fh, indent=1)
