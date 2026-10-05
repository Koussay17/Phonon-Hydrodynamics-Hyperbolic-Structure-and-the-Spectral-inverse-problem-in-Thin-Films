"""Transfer of an approximately feasible witness to the EXACT data (implicit-function / Newton step).

A witness with relative constraint residual ~1e-14 is exactly feasible only for perturbed data. For the
inner-interval claim "tau_max(exact AlN data) >= 6.9e3 ps" one needs an exactly feasible point nearby.
In log-rate variables u = log gamma (positivity automatic), the constraints are
  F(u) = ( (W_gamma e^u - r)/r ,  (K_xx(u) - K0_xx)/K0_xx , (K_zz(u) - K0_zz)/K0_zz ) = 0.
If J = dF/du has full row rank with smallest singular value s_min, the minimum-norm Newton correction has
|du|_2 <= |F|_2 / s_min, and tau changes by about |grad_u tau| |du|.  (A full Kantorovich bound also needs a
Lipschitz constant of J on the ball; with |du| ~ 1e-9 the second-order term is far below the first-order one.)
Run: python -B check_witness_ift.py  (writes check_witness_ift.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts"))

import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402

from aln_geometry import load  # noqa: E402
from symmetry import orbits  # noqa: E402


def main():
    geom = load(str(HERE.parent / "results" / "aln" / "events_m553_s0.1.npz"))
    tie = orbits(geom, geom.labels["maps"])
    orb = tie["ev_orb"]
    sizes = np.bincount(orb)
    gref = (np.bincount(orb, weights=geom.gphys) / sizes)[orb]
    r = geom.W @ gref
    reps = tie["reps"]
    K0 = {c: float(geom.b[:, c] @ geom.response(gref, c)["x"]) for c in (0, 2)}
    out = {}
    for name in ("full_c0_max_08_vertex", "full_c0_min_02_logn1.0", "full_c0_max_04_logn3.0"):
        gam = np.load(HERE.parent / "results" / "aln" / "scan_m553_s0.1" / f"{name}.npz")["gamma"]
        g = gam[orb]
        A = geom.A.tocsc()
        res = {}
        Fd = (geom.W @ g - r)[reps] / r[reps]
        xs = {c: geom.response(g, c, want_grad=True) for c in (0, 2)}
        FK = np.array([(xs[c]["K"] - K0[c]) / K0[c] for c in (0, 2)])
        Fv = np.concatenate([Fd, FK])
        # Jacobian in log-orbit variables
        Wg = np.zeros((len(reps), len(gam)))
        Wd = geom.W.tocsr()[reps].toarray()                    # (117, m)
        for o_row in range(len(reps)):
            Wg[o_row] = np.bincount(orb, weights=Wd[o_row], minlength=len(gam))
        Jd = Wg * gam[None, :] / r[reps][:, None]
        JK = []
        for c in (0, 2):
            p = A.T @ xs[c]["x"]
            dK_dgam = np.bincount(orb, weights=-(p * p), minlength=len(gam))
            JK.append(dK_dgam * gam / K0[c])
        J = np.vstack([Jd, np.array(JK)])
        s = sla.svdvals(J)
        smin = float(s[-1])
        du_bound = float(np.linalg.norm(Fv) / smin)
        # tau gradient in log variables for the objective column 0
        p0 = A.T @ xs[0]["x"]
        q0 = A.T @ xs[0]["u"]
        dN = np.bincount(orb, weights=-2.0 * p0 * q0, minlength=len(gam)) * gam
        dK = np.bincount(orb, weights=-(p0 * p0), minlength=len(gam)) * gam
        tau = xs[0]["tau"]
        dtau = (dN - tau * dK) / xs[0]["K"]
        res.update({"tau_ps": tau / (4 * np.pi), "residual_norm": float(np.linalg.norm(Fv)),
                    "residual_inf": float(np.abs(Fv).max()), "J_rank_rows": int(J.shape[0]),
                    "J_sigma_min": smin, "J_sigma_max": float(s[0]),
                    "newton_step_bound_|du|_2": du_bound,
                    "first_order_rel_change_tau_bound": float(np.linalg.norm(dtau) * du_bound / tau),
                    "min_rate_over_ref_orbit_rate": float(np.min(gam / (np.bincount(orb, weights=gref) / sizes)))})
        out[name] = res
    (HERE / "check_witness_ift.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
