"""Checks of the optimisation model (task 3 / task 5):
 (1) analytic d(log tau)/dy, dK/dy, dr/dy versus central finite differences along random directions;
 (2) exactness of the linear amplitude model: Phi0 + G c versus the independent engine evaluated on the
     explicit perturbed fc3 tensor (fc3_0 + ccm B c) for a random symmetric direction;
 (3) linear model versus phono3py's Rust kernel for the saved variant fc3_iso_0.10_s0.npy (all events,
     orbit-averaged phono3py-route rates versus block-averaged model rates: kappa, tau).
usage: python -B check_model.py <events.npz> <resultdir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.dont_write_bytecode = True
from optimize_fc3 import FC3Model  # noqa: E402
from amp_engine import AmpEngine  # noqa: E402
from fc3class import FOURPI  # noqa: E402


def main():
    ev_path, rdir = sys.argv[1], Path(sys.argv[2])
    model = FC3Model(ev_path, rdir)
    rng = np.random.default_rng(7)
    n = model.n_var
    out = {"fd": []}
    y = rng.standard_normal(n); y *= 0.02 / np.linalg.norm(y)
    ev = model.evaluate(y, grad=True)
    for k in range(3):
        d = rng.standard_normal(n); d /= np.linalg.norm(d)
        h = 1e-5
        ep, em = model.evaluate(y + h * d), model.evaluate(y - h * d)
        rec = {}
        for c in (0, 2):
            fd = (np.log(ep["tau"][c]) - np.log(em["tau"][c])) / (2 * h)
            Dorb = ev["dN"][c] / ev["N"][c] - ev["dK"][c] / ev["K"][c]
            an = model.chain(Dorb, ev["Phi"]) @ d
            rec[f"dlogtau_{c}"] = [float(fd), float(an)]
            fdK = (ep["K"][c] - em["K"][c]) / (2 * h)
            anK = model.chain(ev["dK"][c], ev["Phi"]) @ d
            rec[f"dK_{c}"] = [float(fdK), float(anK)]
        wm = rng.standard_normal(len(model.r0))
        fdr = (wm @ ep["r"] - wm @ em["r"]) / (2 * h)
        anr = model.chain(model.dr_dorb(wm), ev["Phi"]) @ d
        rec["dr"] = [float(fdr), float(anr)]
        out["fd"].append(rec)
        print(json.dumps(rec), flush=True)
    # (2) exactness of the linear model against the engine on an explicit fc3
    Z = np.load(rdir / "linearization.npz")
    ccm = sp.load_npz(rdir / "basis_ccm.npz").tocsr()
    B = np.load(rdir / "basis_B.npy")
    f0 = np.load(rdir / "fc3_ref_compact.npy")
    eng = AmpEngine(rdir / "static_r01.npz", r0avg=True, unit_conversion=float(Z["uc"]))
    yv = rng.standard_normal(n); yv *= 0.1 / np.linalg.norm(yv)
    fc3 = f0 + (ccm @ (B @ (yv * model.c0n))).reshape(f0.shape)
    phi_eng = eng.phi_events(fc3, Z["rows_tp"], Z["rows_b"])
    phi_lin = model.Phi0 + model.Gr @ yv + 1j * (model.Gi @ yv)
    out["linear_model_vs_engine_relmax"] = float(np.max(np.abs(phi_eng - phi_lin)) / np.max(np.abs(phi_eng)))
    out["linear_model_vs_engine_rel_median"] = float(np.median(np.abs(phi_eng - phi_lin) / np.abs(phi_eng)))
    print(json.dumps({k: v for k, v in out.items() if k != "fd"}), flush=True)
    # (3) versus phono3py Rust (saved variant)
    fv = rdir / "fc3_iso_0.10_s0.npy"
    V = np.load(rdir / "variants_iso.npz")
    names = list(V["names"])
    if fv.exists() and "iso_0.10_s0" in names:
        fc3v = np.load(fv)
        eta = 18.0 * (ccm.T @ (fc3v - f0).ravel())
        cv = B.T @ eta
        yv2 = cv / model.c0n
        ev2 = model.evaluate(yv2)
        pp = V["pp"][names.index("iso_0.10_s0")]
        M = model.M
        m_orb = M.metrics(M.rates(pp))
        # block-averaged phono3py-route rates are not available; report both
        out["iso_0.10_s0"] = {"model_tau_x": ev2["tau"][0] / FOURPI, "model_tau_z": ev2["tau"][2] / FOURPI,
                              "model_kappa_x": M.kf * ev2["K"][0], "model_kappa_z": M.kf * ev2["K"][2],
                              "p3_orbitavg_tau_x": m_orb["tau_x"], "p3_orbitavg_tau_z": m_orb["tau_z"],
                              "p3_orbitavg_kappa_x": m_orb["kappa_x"], "p3_orbitavg_kappa_z": m_orb["kappa_z"],
                              "eps_from_c": float(np.linalg.norm(yv2))}
        print(json.dumps(out["iso_0.10_s0"]), flush=True)
    (rdir / "check_model.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
