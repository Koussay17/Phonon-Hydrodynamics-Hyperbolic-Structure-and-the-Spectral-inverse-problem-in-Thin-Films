"""Numerical audit, part 4: gauge sensitivity inside degenerate subspaces, and the e-component of b.

Inside a degenerate subspace (same q, |df| <= 1e-6 THz) the event rates refer to phono3py's eigenvector basis,
whereas phonopy's group velocities of degenerate states come from degenerate perturbation theory, direction by
direction; the relative basis of b and C inside such subspaces is therefore not fixed by the data. On even c-meshes
(9x9x6) the whole plane q_z = 1/2 consists of band pairs that stick together (6_3 screw axis plus time reversal) with
group velocities +-v_z. Test: K, N, tau for
  (a) b as stored;
  (b) b with the components permuted (order reversed) inside every degenerate group;
  (c) b with the components inside every degenerate group rotated by a random orthogonal matrix (5 draws). NOTE: (c)
      is NOT a gauge transformation (b_mu = d f <mu|v|mu> transforms as the diagonal of U^T V U, not as U b); it is
      kept only as a diagnostic of how strongly tau depends on the distribution of the current among degenerate
      partners. The legitimate ambiguity is (b), the assignment of the velocity eigenvalues to the partners.
  (d) the exchange symmetry of the operator itself: max |P C P - C| / max |C| for the reversal P inside every
      degenerate group (K, N are invariant under (b) exactly when P C P = C).
Also reported: the component of b along e and its contribution to K and N at the campaign's regularisation
gamma = tr C / n.
usage: python -B na_gauge.py m553|m775|m996 [out.json]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

import na_common as nc  # noqa: E402
from na_transfer import degeneracy_groups  # noqa: E402


def main():
    mk = sys.argv[1]
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "na_out" / f"gauge_{mk}.json"
    op = nc.Operator(nc.MESHES[mk][0])
    P = np.load(nc.MESHES[mk][1] + "/pp_r01_REF.npz")
    g = op.group_average(op.rates_from_pp(P["pp"], P["g0w"], float(P["conv"])))
    C = op.C(g)
    gam = np.trace(C) / op.n
    grp = degeneracy_groups(op)
    sizes = np.bincount(grp)
    deg_modes = np.nonzero(sizes[grp] > 1)[0]
    groups = [np.nonzero(grp == k)[0] for k in np.nonzero(sizes > 1)[0]]
    zface = np.abs(op.addr[op.q_of][:, 2]) * 2 == op.mesh[2]
    rng = np.random.default_rng(1)
    variants = {"stored": op.b.copy()}
    bp = op.b.copy()
    for ii in groups:
        bp[ii] = op.b[ii[::-1]]
    variants["permuted"] = bp
    for k in range(5):
        br = op.b.copy()
        for ii in groups:
            Qr, _ = np.linalg.qr(rng.standard_normal((len(ii), len(ii))))
            br[ii] = Qr @ op.b[ii]
        variants[f"rotated{k}"] = br
    perm = np.arange(op.n)
    for ii in groups:
        perm[ii] = ii[::-1]
    asym = float(np.max(np.abs(C[np.ix_(perm, perm)] - C)) / np.max(np.abs(C)))
    rows_deg = np.concatenate(groups) if groups else np.zeros(0, dtype=np.int64)
    asym_rows = float(np.max(np.abs(C[np.ix_(perm[rows_deg], perm)] - C[rows_deg])) / np.max(np.abs(C))) if len(rows_deg) else 0.0
    Creg = C + gam * np.outer(op.ehat, op.ehat)
    del C
    cf = sla.cho_factor(Creg, lower=True, overwrite_a=True, check_finite=False)
    R = {"mesh": op.mesh, "n_degenerate_modes": int(len(deg_modes)), "n_degenerate_groups": len(groups),
         "n_modes_on_qz_half_plane": int(zface.sum()), "exchange_asymmetry_PCP_minus_C_rel": asym,
         "exchange_asymmetry_degenerate_rows_rel": asym_rows, "variants": {}}
    for nm, B in variants.items():
        rec = {}
        for c in (0, 2):
            b = B[:, c]
            be = float(b @ op.ehat)
            bH = b - be * op.ehat
            x = sla.cho_solve(cf, bH, check_finite=False)
            K, N = float(bH @ x), float(x @ x)
            rec["xyz"[c]] = {"tau_ps": N / K / nc.FOURPI, "kappa": op.kf * K, "b_dot_ehat_rel": abs(be) / float(np.linalg.norm(b)),
                             "e_part_of_K_rel_at_gamma": be ** 2 / gam / K, "e_part_of_N_rel_at_gamma": (be / gam) ** 2 / N}
        R["variants"][nm] = rec
        print(nm, json.dumps({a: (round(rec[a]["tau_ps"], 6), round(rec[a]["kappa"], 4)) for a in rec}), flush=True)
    for a in "xz":
        t = [v[a]["tau_ps"] for v in R["variants"].values()]
        k = [v[a]["kappa"] for v in R["variants"].values()]
        R[f"spread_rel_tau_{a}"] = (max(t) - min(t)) / R["variants"]["stored"][a]["tau_ps"]
        R[f"spread_rel_kappa_{a}"] = (max(k) - min(k)) / R["variants"]["stored"][a]["kappa"]
    R["permuted_minus_stored_rel"] = {a: R["variants"]["permuted"][a]["tau_ps"] / R["variants"]["stored"][a]["tau_ps"] - 1 for a in "xz"}
    print({k: v for k, v in R.items() if k.startswith("spread") or k.startswith("n_") or k.startswith("exch") or k.startswith("perm")})
    out_path.write_text(json.dumps(R, indent=1))


if __name__ == "__main__":
    main()
