"""Attack 1c: are the numerically-zero reference vertices group-theoretic selection rules?

Reconstruct the 24 q-point maps as integer matrices R (addr[map q] = R addr[q] mod mesh), split them
into unitary C6v operations (R_zz = +1, q_z preserved) and antiunitary ones (time reversal x R, q_z ->
-q_z). A selection-rule zero of Phi(q_p j_p; q_a j_a; q_b j_b) needs a nontrivial unitary operation
fixing all three wave vectors (a common little group); for such events the vertex vanishes when the
product of the three phonons' characters does not contain the identity.
Output: review/cx_selection_rules.json
"""
import itertools, json, sys
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:\ResearchLab\orchestration\campaigns\20261004-memory-gap\review")
from cx_common import AlN, CAMP

M = AlN(); gref = M.g_ref()
forb = gref < 1e-12 * gref.max()
Z = M.Z
addr = Z["addr"]; mesh = Z["mesh"]; rot = Z["rot_maps"]
Nq = len(addr)
# integer matrices
mats = []
cands = [np.array(v) for v in itertools.product((-1, 0, 1), repeat=3)]
for k in range(len(rot)):
    found = None
    for rows in itertools.product(cands, repeat=3):
        R = np.array(rows)
        if abs(round(np.linalg.det(R))) != 1:
            continue
        img = (addr @ R.T) % mesh
        tgt = addr[rot[k]] % mesh
        if np.array_equal(img, tgt):
            found = R
            break
    assert found is not None, k
    mats.append(found)
mats = np.array(mats)
unitary = mats[:, 2, 2] == 1
print("unitary ops:", int(unitary.sum()), "antiunitary:", int((~unitary).sum()))
qidx = np.nonzero(M.valid)[0] // M.nb
qP, qA, qB = qidx[M.P], qidx[M.A], qidx[M.B]
fix = rot == np.arange(Nq)[None, :]            # fix[k, q]: op k fixes q
common_u = np.zeros(M.m, dtype=int)
common_all = np.zeros(M.m, dtype=int)
for k in range(len(rot)):
    f3 = fix[k, qP] & fix[k, qA] & fix[k, qB]
    common_all += f3
    if unitary[k]:
        common_u += f3
out = {"n_unitary": int(unitary.sum()),
       "forbidden_common_unitary_stabilizer_hist": {int(a): int(b) for a, b in zip(*np.unique(common_u[forb], return_counts=True))},
       "allowed_common_unitary_stabilizer_hist": {int(a): int(b) for a, b in zip(*np.unique(common_u[~forb], return_counts=True))},
       "frac_forbidden_with_nontrivial_common_unitary_stabilizer": float(np.mean(common_u[forb] > 1)),
       "frac_allowed_with_nontrivial_common_unitary_stabilizer": float(np.mean(common_u[~forb] > 1))}
# among events with a nontrivial common stabilizer, fraction forbidden
nt = common_u > 1
out["frac_forbidden_among_nontrivial_stabilizer_events"] = float(np.mean(forb[nt]))
out["n_events_nontrivial_stabilizer"] = int(nt.sum())
print(json.dumps(out, indent=1))
(CAMP / "review" / "cx_selection_rules.json").write_text(json.dumps(out, indent=1))
