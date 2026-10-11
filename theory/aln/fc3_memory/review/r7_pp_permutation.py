"""Review R7 (phono3py env): permutation symmetry of phono3py's interaction strengths on the fixed event list.
For every event p -> a + b the strength |Phi|^2 is computed from three perspectives:
  P: q0 = q_p, triplet (q_p, -q_a, -q_b), bands (j_p, j_a, j_b)  (the campaign's convention, p3_events.py)
  A: q0 = q_a, triplet (q_a,  q_b, -q_p), bands (j_a, j_b, j_p)  (row of mode a in phono3py's collision matrix)
  B: q0 = q_b, triplet (q_b,  q_a, -q_p), bands (j_b, j_a, j_p)  (row of mode b)
For a permutation-symmetric fc3 and a consistent Fourier interpolation the three agree. Degenerate bands make
single values gauge dependent, so the comparison is made on sums over degenerate blocks (gauge invariant).
usage: python -B r7_pp_permutation.py <events.npz> <outdir> name=fc3.npy|REF ...
"""
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("RAYON_NUM_THREADS", "2")
sys.dont_write_bytecode = True
sys.path.insert(0, r"C:\Users\Koussay\ResearchLab\orchestration\campaigns\20261005-fc3-memory\scripts")
import numpy as np  # noqa: E402

import aln_common as common  # noqa: E402
from p3_events import load_ph3  # noqa: E402


def pp_perspectives(ph3, EV, sigma=0.1, T=300.0):
    from phono3py.phonon3.collision_matrix import CollisionMatrix
    ph3.init_phph_interaction()
    ph3.run_phonon_solver()
    itr = ph3.phph_interaction
    bz = itr.bz_grid
    N = int(np.prod(bz.D_diag))
    nb = itr.get_phonons()[0].shape[1]
    grg2bzg = bz.grg2bzg
    Jmap = EV["Jmap"]
    p, a, b = EV["ev_p"], EV["ev_a"], EV["ev_b"]
    qp, jp = p // nb, p % nb
    qa, ja = a // nb, a % nb
    qb, jb = b // nb, b % nb
    m = len(p)
    req = {"P": (qp, Jmap[qa], jp, ja, jb), "A": (qa, qb, ja, jb, jp), "B": (qb, qa, jb, ja, jp)}
    out = {k: np.full(m, np.nan) for k in req}
    col = CollisionMatrix(itr, rot_grid_points=None, lang="Rust")
    for g0 in range(N):
        col.set_grid_point(int(grg2bzg[g0]))
        col.set_sigma(sigma, sigma_cutoff=None)
        col.run_integration_weights()
        col.run_interaction(is_full_pp=False)
        pp = col._pp_strength
        gp2tp, tp2s, swapped = col._get_gp2tp_map()
        for k, (q0, q1, j0, j1, j2) in req.items():
            idx = np.nonzero(q0 == g0)[0]
            if len(idx) == 0:
                continue
            g1 = q1[idx]
            ti = gp2tp[g1]
            assert (ti >= 0).all()
            sw = swapped[g1].astype(bool)
            b1 = np.where(sw, j2[idx], j1[idx])
            b2 = np.where(sw, j1[idx], j2[idx])
            out[k][idx] = pp[ti, j0[idx], b1, b2]
    return out


def main():
    evp, outdir = sys.argv[1], Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    EV = np.load(evp)
    common.verify_inputs()
    ph3 = load_ph3(True, [5, 5, 3], 0.1)
    fc3_ref = np.array(ph3.fc3)
    f = EV["freqs"]
    nb = f.shape[1]
    # degenerate-block keys: (q, degenerate set) for each of the three modes
    def dkey(mode):
        q, j = mode // nb, mode % nb
        fq = f[q]
        lo = np.array([np.min(np.nonzero(np.abs(fq - fq[jj]) < 1e-6)[0]) for jj in range(nb)])
        return q * nb + lo[j]
    kp = np.array([dkey(x) for x in EV["ev_p"]])
    ka = np.array([dkey(x) for x in EV["ev_a"]])
    kb = np.array([dkey(x) for x in EV["ev_b"]])
    lo_, hi_ = np.minimum(ka, kb), np.maximum(ka, kb)
    _, blk = np.unique(np.stack([kp, lo_, hi_], 1), axis=0, return_inverse=True)
    blk = blk.ravel()
    res = {}
    for spec in sys.argv[3:]:
        name, path = spec.split("=", 1)
        fc3 = fc3_ref if path == "REF" else np.load(path)
        ph3.fc3 = np.ascontiguousarray(fc3, dtype="double")
        t0 = time.time()
        o = pp_perspectives(ph3, EV)
        np.savez_compressed(outdir / f"pp3_{name}.npz", **o)
        S = {k: np.bincount(blk, weights=v) for k, v in o.items()}
        big = S["P"] > 1e-6 * S["P"].max()
        rec = {"seconds": time.time() - t0, "n_blocks": int(len(S["P"])), "n_blocks_nonzero": int(big.sum())}
        for k in ("A", "B"):
            r = S[k][big] / S["P"][big] - 1
            rec[f"block_rel_{k}_vs_P"] = {"max_abs": float(np.abs(r).max()), "q99": float(np.quantile(np.abs(r), 0.99)),
                                          "q50": float(np.quantile(np.abs(r), 0.5)),
                                          "n_over_1pct": int((np.abs(r) > 0.01).sum()), "n_over_10pct": int((np.abs(r) > 0.1).sum())}
        res[name] = rec
        print(name, json.dumps(rec), flush=True)
    (outdir / "r7_pp_permutation.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
