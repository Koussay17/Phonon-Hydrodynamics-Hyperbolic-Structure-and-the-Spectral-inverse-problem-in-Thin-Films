"""Task 3 (outer side): rigorous memory ceiling for norm-bounded symmetric fc3 perturbations.

For fc3 = fc3_0 + ccm B c with |c| <= rho |c0| (relative Frobenius size <= rho):
    |Phi_m(c)| >= max(0, |Phi0_m| - ||G_m|| rho |c0|) =: phi_m,   ||G_m|| = operator norm of the real-linear
    map c -> Phi_m (largest singular value of the 2 x n real matrix [Re G_m; Im G_m]),
so every (block-averaged) orbit rate obeys g_O >= g_low_O(rho) = W_O uc mean_m phi_m^2. Then C(g) - C(g_low)
is PSD, lambda_min(C(g)|H) >= lambda_min(C(g_low)|H) and, since C^+ b lies in H,
    tau_mem(g) <= 1 / lambda_min(C(g_low(rho))|H)        (P7 argument with per-event floors).
Valid for every symmetric fc3 in the ball, with or without lifetime / kappa data, both directions.
Also reported: the per-orbit 'fragility' rho*_O = min over the block of |Phi0_m| / (||G_m|| |c0|) (the smallest
relative fc3 change that can silence the orbit), and how many orbits can be silenced at each rho.
usage: python -B t3_ceiling.py <events.npz> <resultdir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI  # noqa: E402


def main():
    ev_path, rdir = sys.argv[1], Path(sys.argv[2])
    M = FC3Class(ev_path, rdir / "pp_r01_REF.npz")
    Z = np.load(rdir / "linearization.npz")
    Phi0, opn, row_orb, blk = Z["Phi0"], Z["opn"], Z["row_orb"], Z["blk"]
    W_orb, allowed, uc = Z["W_orb"], Z["allowed"], float(Z["uc"])
    c0n = float(np.linalg.norm(Z["c0"]))
    a0 = np.abs(Phi0)
    frag_row = a0 / (opn * c0n)
    frag_orb = np.full(M.n_orb, np.inf)
    np.minimum.at(frag_orb, row_orb, frag_row)
    # rate share of each orbit in the reference lifetimes
    g_blk = Z["g_blk_orb"]
    out = {"fragility_orbit_quantiles(0,0.01,0.05,0.1,0.25,0.5)":
           np.quantile(frag_orb[allowed], [0, 0.01, 0.05, 0.1, 0.25, 0.5]).tolist(), "rho": {}}
    base = M.metrics(g_blk[M.orb], cols=(0, 2))
    out["base"] = {"tau_x": base["tau_x"], "tau_z": base["tau_z"], "lam_min_H": base["lam_min_H"],
                   "ceiling_F1_ps": base["ceiling_F1_ps"], "tauCS_x": base["tauCS_x"], "tauCS_z": base["tauCS_z"]}
    for rho in [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
        phi = np.maximum(0.0, a0 - opn * rho * c0n)
        glow = np.zeros(M.n_orb)
        glow += np.bincount(row_orb, weights=W_orb[row_orb] * uc * phi ** 2 / blk[row_orb], minlength=M.n_orb)
        m = M.metrics(glow[M.orb], cols=())
        n_sil = int(np.sum(frag_orb[allowed] <= rho))
        # factor-F equivalent: largest per-orbit suppression allowed by the floors
        with np.errstate(divide="ignore"):
            Fmax = float(np.max(np.where(allowed, g_blk / np.maximum(glow, 1e-300), 1.0)))
        out["rho"][f"{rho:g}"] = {"lam_min_H_low": m["lam_min_H"],
                                  "tau_ceiling_ps": (1 / m["lam_min_H"] / FOURPI) if m["lam_min_H"] > 0 else float("inf"),
                                  "n_orbits_silenceable": n_sil, "max_orbit_suppression_factor": Fmax}
        print(rho, json.dumps(out["rho"][f"{rho:g}"]), flush=True)
    (rdir / "t3_ceiling.json").write_text(json.dumps(out, indent=1))
    np.save(rdir / "fragility_orbit.npy", frag_orb)


if __name__ == "__main__":
    main()
