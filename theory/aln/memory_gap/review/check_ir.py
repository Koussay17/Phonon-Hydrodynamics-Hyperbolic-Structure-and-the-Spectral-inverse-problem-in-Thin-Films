r"""T7 (infrared criterion) checks (13-final-report.md, T7).

Part A - RTA mesh laws, d = 3, exact lattice sums (no fitting of the model):
  modes k in {-(N-1)/2..(N-1)/2}^3 \ {0}, q = 2 pi k / N, weight 1/N^3, b^2 = Omega_x^2 (classical limit
  b -> kT v with kT c = 1), r = |q|^alpha.  M_m = sum b^2 r^-m / N^3, tau = M_2 / M_1.
  Predictions (moment criterion): alpha < 3/2: converges; alpha = 3/2: ~ log N; 3/2 < alpha < 3: ~ N^(2 alpha - 3);
  alpha = 3: ~ N^3 / log N (the report states ~ N^alpha for alpha >= d); alpha > 3: ~ N^alpha.
Part B - is H4 (|x_mu| >= c |b_mu| / r_mu in the infrared) a property of the class or of the operator?
  Ratio rho_mu = x_mu r_mu / b_mu for the lowest-frequency modes (f < 5 THz) of the AlN 5x5x3 event operator at the
  reference rates and at the campaign's own extremal witnesses (same lifetimes, same K_xx, K_zz).
Part C - the one-line inequality tau >= |x| / |b| (from K = b.x <= |b||x|), which makes tau diverge whenever
  N = |x|^2 diverges with |b| bounded, without any upper bound on K.
Run: python -B check_ir.py   (writes check_ir.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

import numpy as np  # noqa: E402


def rta_sums(N, alphas):
    h = (N - 1) // 2
    k = np.arange(-h, h + 1)
    out = {a: [0.0, 0.0] for a in alphas}
    for kx in k:
        KY, KZ = np.meshgrid(k, k, indexing="ij")
        q2 = (kx ** 2 + KY ** 2 + KZ ** 2).astype(float)
        mask = q2 > 0
        q2 = q2[mask] * (2 * np.pi / N) ** 2
        ox2 = (kx ** 2) * (2 * np.pi / N) ** 2 / q2
        for a in alphas:
            r = q2 ** (a / 2)
            out[a][0] += np.sum(ox2 / r)
            out[a][1] += np.sum(ox2 / r ** 2)
    return {a: (v[0] / N ** 3, v[1] / N ** 3) for a, v in out.items()}


def part_a():
    alphas = [1.0, 1.5, 2.0, 3.0, 3.5]
    Ns = [21, 41, 81, 161]
    taus = {a: [] for a in alphas}
    for N in Ns:
        s = rta_sums(N, alphas)
        for a in alphas:
            taus[a].append(s[a][1] / s[a][0])
    res = {}
    for a in alphas:
        t = np.array(taus[a])
        loc = np.diff(np.log(t)) / np.diff(np.log(Ns))
        res[str(a)] = {"tau": t.tolist(), "local_exponents": loc.tolist(),
                       "tau_over_logN": (t / np.log(Ns)).tolist(),
                       "tau_logN_over_N3": (t * np.log(Ns) / np.array(Ns, float) ** 3).tolist()}
    return {"Ns": Ns, "alpha": res}


def part_b():
    from aln_geometry import load
    from symmetry import orbits
    ev = HERE.parent / "results" / "aln" / "events_m553_s0.1.npz"
    geom = load(str(ev))
    tie = orbits(geom, geom.labels["maps"])
    sizes = np.bincount(tie["ev_orb"])
    gref = (np.bincount(tie["ev_orb"], weights=geom.gphys) / sizes)[tie["ev_orb"]]
    r = geom.W @ gref
    wit = HERE.parent / "results" / "aln" / "scan_m553_s0.1"
    cases = {"reference": gref}
    for name in ("full_c0_min_02_logn1.0", "full_c0_max_08_vertex", "full_c2_min_06_vertex"):
        gam = np.load(wit / f"{name}.npz")["gamma"]
        cases[name] = gam[tie["ev_orb"]]
    f = geom.eps
    out = {}
    for name, g in cases.items():
        rr = geom.W @ g
        col = 2 if "c2" in name else 0
        R = geom.response(g, col)
        x, b = R["x"], geom.b[:, col]
        # lowest modes of this mesh: 3.70-4.35 THz (16 modes); the 5x5x3 mesh has nothing below 3.7 THz
        sel = (f < 5.0) & (np.abs(b) > 1e-3 * np.abs(b).max())
        rho = x[sel] * rr[sel] / b[sel]
        out[name] = {"col": col, "max_rel_diag_dev_from_ref": float(np.max(np.abs(rr - r) / r)),
                     "n_low_modes": int(sel.sum()), "tau_ps": R["tau"] / (4 * np.pi),
                     "rho_quantiles_10_50_90": np.quantile(rho, [0.1, 0.5, 0.9]).tolist(),
                     "fraction_rho_ge_0.5": float(np.mean(rho >= 0.5)),
                     "fraction_rho_le_0.1": float(np.mean(np.abs(rho) <= 0.1)),
                     "share_N_below_5THz": float(np.sum(x[f < 5.0] ** 2) / np.sum(x ** 2)),
                     "tau_ge_|x|/|b|": bool(R["tau"] >= np.linalg.norm(x) / np.linalg.norm(b)),
                     "tau_over_(|x|/|b|)": float(R["tau"] / (np.linalg.norm(x) / np.linalg.norm(b)))}
    return out


def main():
    rep = {"A_rta_mesh_laws_d3": part_a()}
    rep["B_H4_diagnostic_AlN"] = part_b()
    (HERE / "check_ir.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
