"""Referee check of the infrared structure of tau_mem = N/K on the three meshes (reference fc3 only).

For each mesh: per-mode contributions N_mu = x_mu^2, K_mu = b_mu x_mu (x = C^+ b, full conserving operator),
low-frequency power laws of r = diag C, of b_mu and of N_mu, smallest |q| star, and the fit
tau(mesh) = a + s / f_min across meshes.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from referee_check import Mesh, ROOT, FOURPI  # noqa: E402

out = {}
for name, ev, pp in (("m553", ROOT / "m553" / "events_m553_s0.1.npz", ROOT / "m553" / "pp_r01_REF.npz"),
                     ("m775", ROOT / "m775x" / "events_m775_s0.1.npz", ROOT / "m775x" / "pp_r01_REF.npz"),
                     ("m996", ROOT / "m996" / "events_m996_s0.1.npz", ROOT / "m996" / "pp_r01_REF.npz")):
    M = Mesh(ev)
    Z = np.load(ev)
    g, _ = M.rates(pp)
    m = M.metrics(g, want_x=True)
    f = M.eps
    r = M.degavg(m["r"])
    rec = {"f_min": float(f.min()), "kT_THz": float(M.kT), "tau_x": m["tau_x"], "tau_z": m["tau_z"],
           "kappa_x": m["kappa_x"], "kappa_z": m["kappa_z"]}
    low = f < 5.0
    A = np.vstack([np.log(f[low]), np.ones(low.sum())]).T
    rec["fit_below_5THz"] = {"n_modes": int(low.sum())}
    # power law of the degeneracy-averaged diagonal rate r = diag C below 5 THz (meaningful only when the window
    # spans a factor > 2 in frequency, i.e. on 7x7x5 and 9x9x6)
    ok = low & (r > 0)
    Aa = np.vstack([np.log(f[ok]), np.ones(ok.sum())]).T
    rec["fit_below_5THz"]["r_exponent"] = float(np.linalg.lstsq(Aa, np.log(r[ok]), rcond=None)[0][0])
    rec["fit_below_5THz"]["f_window"] = [float(f[ok].min()), float(f[ok].max())]
    for nm in ("x", "z"):
        x = m[f"x_{nm}"]
        Nmu = x * x
        # share of N from the modes in the lowest-frequency 'q-star': modes with f < 1.15 f_min (acoustic, smallest |q|)
        star = f < 1.15 * f.min()
        rec[f"N_share_modes_below_1.15fmin_{nm}"] = float(Nmu[star].sum() / Nmu.sum())
        rec[f"n_modes_below_1.15fmin"] = int(star.sum())
    out[name] = rec
    print(name, json.dumps(rec), flush=True)

# tau versus 1/f_min across meshes
fm = np.array([out[k]["f_min"] for k in ("m553", "m775", "m996")])
for nm in ("x", "z"):
    t = np.array([out[k][f"tau_{nm}"] for k in ("m553", "m775", "m996")])
    A = np.vstack([1 / fm, np.ones(3)]).T
    (s, a), res, *_ = np.linalg.lstsq(A, t, rcond=None)
    pred = A @ np.array([s, a])
    out[f"fit_tau_{nm}_vs_inv_fmin"] = {"slope_ps_THz": float(s), "intercept_ps": float(a),
                                        "max_rel_residual": float(np.max(np.abs(pred / t - 1))),
                                        "pairwise_slopes": [float((t[1] - t[0]) / (1 / fm[1] - 1 / fm[0])),
                                                            float((t[2] - t[1]) / (1 / fm[2] - 1 / fm[1]))]}
    k = np.array([out[kk][f"kappa_{nm}"] for kk in ("m553", "m775", "m996")])
    out[f"kappa_{nm}_by_mesh"] = k.tolist()
print(json.dumps({k: v for k, v in out.items() if k.startswith("fit_tau") or k.startswith("kappa_")}, indent=1))
Path(Path(__file__).resolve().parent / "referee_ir.json").write_text(json.dumps(out, indent=1))
