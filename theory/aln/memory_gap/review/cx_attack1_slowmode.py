"""Attack 1d: what is the hidden slow mode of each witness? Is it a physical quasi-invariant
(crystal momentum of a sub-population, i.e. a phonon-hydrodynamic mode) or an arbitrary combination?

For each witness: slowest eigenvectors of C on H; weight by frequency band and branch; overlap with
the crystal-momentum vectors P_c = D * q_c (q in Cartesian 1/Angstrom, BZ representative), restricted
to modes below a cutoff f_c; the momentum-relaxation rate of P_c restricted to f < f_c under the
witness vs the reference (Rayleigh quotient P^T C P / P^T P); fraction of normal vs umklapp events
among the events touching the slow-mode support that the witness suppresses.
Output: review/cx_attack1_slowmode.json
"""
import json
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP, FOURPI, their_ev_orb, load_witness, basis_H  # noqa: E402

M = AlN()
gref = M.g_ref()
ev_orb_theirs, _ = their_ev_orb()
Z = M.Z
addr = Z["addr"].astype(float)              # BZ-representative addresses (GR grid)
mesh = Z["mesh"].astype(float)
lat = Z["lattice"]                            # rows: real-space lattice vectors (Angstrom)
rec_lat = 2 * np.pi * np.linalg.inv(lat).T    # rows: reciprocal vectors
qidx = np.nonzero(M.valid)[0] // M.nb
band = np.nonzero(M.valid)[0] % M.nb
qcart = (addr[qidx] / mesh) @ rec_lat         # (n, 3)
Q = basis_H(M.ehat)


def proj_H(v):
    return v - M.ehat * (M.ehat @ v)


out = {}
for nm in ["full_c0_max_08_vertex", "full_c0_max_04_logn3.0", "full_c2_max_01_logn0.5", "reference"]:
    g = gref if nm == "reference" else load_witness(nm, ev_orb_theirs)
    C = M.C(g)
    CH = Q.T @ C @ Q
    w, U = np.linalg.eigh(0.5 * (CH + CH.T))
    rec = {"lam_lowest4": [float(x) for x in w[:4]]}
    # slowest two eigenvectors (often a degenerate pair)
    V = Q @ U[:, :2]
    wt = np.sum(V ** 2, axis=1)
    wt /= wt.sum()
    rec["slow_weight_below_5THz"] = float(wt[M.eps < 5].sum())
    rec["slow_weight_5_10THz"] = float(wt[(M.eps >= 5) & (M.eps < 10)].sum())
    rec["slow_weight_above_10THz"] = float(wt[M.eps >= 10].sum())
    rec["slow_weight_by_band"] = [float(wt[band == j].sum()) for j in range(M.nb)]
    # overlap with restricted crystal momentum
    ov = {}
    for fc in (4.5, 5.0, 6.0, 8.0, 10.0, 30.0):
        best = 0.0
        for c in range(3):
            P = proj_H(M.D * qcart[:, c] * (M.eps < fc))
            if np.linalg.norm(P) == 0:
                continue
            P /= np.linalg.norm(P)
            s = float(np.sum((V.T @ P) ** 2))          # squared overlap with the slow 2-space
            relax = float(P @ C @ P)                   # momentum relaxation rate (Rayleigh quotient)
            relax_ref = float(P @ M.C(gref) @ P)
            best = max(best, s)
            ov[f"fc{fc}_c{c}"] = {"overlap2_slow2space": s, "rate": relax, "rate_ref": relax_ref,
                                   "rate_ratio": relax / relax_ref}
        ov[f"fc{fc}_best"] = best
    rec["momentum"] = ov
    out[nm] = rec
    print(nm, json.dumps({k: v for k, v in rec.items() if k != "momentum"}), flush=True)
    print("   momentum best overlaps:", {k: round(v, 4) for k, v in ov.items() if k.endswith("best")}, flush=True)
(CAMP / "review" / "cx_attack1_slowmode.json").write_text(json.dumps(out, indent=1))
