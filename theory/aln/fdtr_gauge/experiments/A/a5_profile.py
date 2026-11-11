"""A5: profile-likelihood checks of the Fisher numbers (N1, N2, N3), 500 nm Cattaneo, finite beams.

(a) Noiseless ("Asimov") profile with the twin: data = twin at the base point; for fixed ln tau_R = ln tau_0 + delta
    the nuisances of the set are refitted (instrument g, phi0, t_d by exact linear projection, thermal ones by
    Levenberg-Marquardt); Delta chi2_min(delta) is compared with the Fisher prediction (delta / sigma_F)^2.
    Includes the Fourier film (tau_R = 0) as the extreme alternative.
(b) Noisy synthetic data (twin `synthetic_measurement`, seeds 20261111 + i): joint maximum-likelihood fit of
    ln tau_R and the nuisances with the independent model (indep.py, exact autodiff Jacobian); the scatter of the
    estimates over seeds is compared with sigma_F. A few seeds are refitted with the twin as a check.
Usage: python a5_profile.py -> out/a5_profile.json
"""
from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import least_squares

from common import (FREQ, OUT, SETS_CAT, Model, model_for, jacobian, weighted, fisher, jsonable, WP, WR, tw)
from indep import Indep, base_from_model, PNAMES

T0 = time.perf_counter()
L = 500
m0 = model_for(L)
INSTR = ("g", "phi0", "td_ns")


def design(names, f, sa, sp):
    cols = []
    for n in names:
        if n == "g":
            cols.append(np.concatenate([np.ones(f.size) / sa, np.zeros(f.size)]))
        elif n == "phi0":
            cols.append(np.concatenate([np.zeros(f.size), np.ones(f.size) / sp]))
        elif n == "td_ns":
            cols.append(np.concatenate([np.zeros(f.size), -2 * np.pi * f * 1e-9 / sp]))
    return np.stack(cols, axis=1) if cols else None


def projector(A):
    if A is None:
        return None
    Q, _ = np.linalg.qr(A)
    return Q


def whiten(d, sa, sp):
    return np.concatenate([d.real / sa, d.imag / sp])


def proj_out(r, Q):
    return r if Q is None else r - Q @ (Q.T @ r)


# ---------------------------------------------------------------- (a) Asimov profile with the twin
def asimov_profile(setname, sig_phi_deg, deltas, include_fourier=True):
    names = SETS_CAT[setname]
    th_names = [n for n in names if n not in INSTR and n != "ln_tauR"]
    sa, sp = 0.01, np.deg2rad(sig_phi_deg)
    Q = projector(design([n for n in names if n in INSTR], FREQ, sa, sp))
    lnHd = m0.lnH(FREQ)

    def resid(x, delta, model=m0, fourier=False):
        th = dict(zip(th_names, x))
        if fourier:
            P = model.phys(th)
            P["tauR"] = 0.0
            h = model.thermal(FREQ, P)
            lnm = np.log(np.abs(h)) + 1j * np.unwrap(np.angle(h))
        else:
            th["ln_tauR"] = delta
            lnm = model.lnH(FREQ, th)
        d = lnHd - lnm
        d = np.log(np.abs(np.exp(d))) + 1j * np.angle(np.exp(d))      # principal branch of ln(H_d/H_m)
        return proj_out(whiten(d, sa, sp), Q)

    out = []
    x = np.zeros(len(th_names))
    for delta in deltas:
        if th_names:
            r = least_squares(resid, x, args=(delta,), method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                              diff_step=1e-4)
            chi2, xs = float(2 * r.cost), r.x.tolist()
        else:
            rr = resid(x, delta)
            chi2, xs = float(rr @ rr), []
        out.append({"delta_ln_tau": float(delta), "dchi2": chi2, "nuisance": dict(zip(th_names, xs))})
    if include_fourier:
        if th_names:
            r = least_squares(resid, np.zeros(len(th_names)), args=(0.0, m0, True), method="lm", xtol=1e-15,
                              ftol=1e-15, gtol=1e-15, diff_step=1e-4)
            chi2, xs = float(2 * r.cost), r.x.tolist()
        else:
            rr = resid(x, 0.0, m0, True)
            chi2, xs = float(rr @ rr), []
        out.append({"delta_ln_tau": "-inf (Fourier film, tau_R = 0)", "dchi2": chi2,
                    "nuisance": dict(zip(th_names, xs))})
    return out


J = jacobian(m0, FREQ, SETS_CAT["N5"], h=2e-3, scheme="richardson")
names5 = SETS_CAT["N5"]
res = {"L_nm": L, "tauR_ps": m0.tauR * 1e12, "asimov": {}, "seeds": {}}
ws = np.array([-0.5, -0.2, -0.08, -0.052, -0.028, 0.028, 0.052, 0.08, 0.2, 0.5, 1.0, 3.0])
deltas = np.log1p(ws)
for N in ("N0", "N1", "N2", "N3"):
    for s in (0.01, 0.1):
        idx = [names5.index(n) for n in SETS_CAT[N]]
        sF = fisher(weighted(J[:, idx], sig_phi_deg=s))["sigma"]
        prof = asimov_profile(N, s, deltas)
        for p, w in zip(prof, list(ws) + [-1.0]):
            p["w"] = float(w)
            dl = p["delta_ln_tau"]
            p["fisher_pred"] = float((dl / sF) ** 2) if not isinstance(dl, str) else None
        res["asimov"][f"{N}|{s}"] = {"sigma_F": sF, "profile": prof}
        print("asimov", N, s, "sigma_F", round(sF, 4),
              [(round(p["w"], 3), round(p["dchi2"], 6), (round(p["fisher_pred"], 6) if p["fisher_pred"] else None))
               for p in prof], f"{time.perf_counter() - T0:.0f}s", flush=True)
        (OUT / "a5_profile.json").write_text(json.dumps(jsonable(res), indent=1))

# ---------------------------------------------------------------- (b) noisy data, many seeds (independent model)
ind = Indep(base_from_model(m0), FREQ)
st = m0.stack(m0.phys())


def fit_indep(lnHd, setname, sa, sp, x0=None):
    names = SETS_CAT[setname]
    th_names = [n for n in names if n not in INSTR]
    idx = [PNAMES.index(n) for n in th_names]
    Q = projector(design([n for n in names if n in INSTR], FREQ, sa, sp))

    def full(x):
        th = np.zeros(len(PNAMES))
        th[idx] = x
        return th

    def resid(x):
        d = lnHd - ind.lnH(full(x))
        d = np.log(np.abs(np.exp(d))) + 1j * np.angle(np.exp(d))
        return proj_out(whiten(d, sa, sp), Q)

    def jac(x):
        Jc = ind.jacobian(th_names, full(x))
        Jw = -np.vstack([Jc.real / sa, Jc.imag / sp])
        return Jw if Q is None else Jw - Q @ (Q.T @ Jw)

    x0 = np.zeros(len(th_names)) if x0 is None else x0
    lo = np.array([-12.0 if n == "ln_tauR" else -3.0 for n in th_names])
    hi = np.array([4.0 if n == "ln_tauR" else 3.0 for n in th_names])
    r = least_squares(resid, x0, jac=jac, method="trf", bounds=(lo, hi), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return r.x[th_names.index("ln_tauR")], float(2 * r.cost), r


NSEED = 200
for N in ("N0", "N1", "N2"):
    for s in (0.01, 0.1):
        idx = [names5.index(n) for n in SETS_CAT[N]]
        sF = fisher(weighted(J[:, idx], sig_phi_deg=s))["sigma"]
        est, chi2s = [], []
        for i in range(NSEED):
            dat = tw.synthetic_measurement(FREQ, st, WP, WR, sigma_log_amplitude=0.01, sigma_phase_deg=s,
                                           seed=20261111 + i, passivity="ignore")
            lnHd = np.log(np.abs(dat["observed"])) + 1j * np.unwrap(np.angle(dat["observed"]))
            e, c2, _ = fit_indep(lnHd, N, 0.01, np.deg2rad(s))
            est.append(float(e))
            chi2s.append(c2)
        est = np.array(est)
        nfree = 2 * FREQ.size - len(SETS_CAT[N])
        qs = np.quantile(est, [0.025, 0.16, 0.5, 0.84, 0.975])
        res["seeds"][f"{N}|{s}"] = {
            "sigma_F": sF, "n_seeds": NSEED, "mean_hat": float(est.mean()), "std_hat": float(est.std(ddof=1)),
            "std_over_sigmaF": float(est.std(ddof=1) / sF), "quantiles_2.5_16_50_84_97.5": qs.tolist(),
            "robust_sigma_half_16_84": float((qs[3] - qs[1]) / 2),
            "frac_at_bounds": float(np.mean((est <= -11.999) | (est >= 3.999))),
            "coverage_1sigmaF": float(np.mean(np.abs(est) < sF)), "coverage_2sigmaF": float(np.mean(np.abs(est) < 2 * sF)),
            "mean_chi2": float(np.mean(chi2s)), "dof": nfree, "estimates": est.tolist()}
        print("seeds", N, s, "sigma_F", round(sF, 4), "std", round(est.std(ddof=1), 4), "mean", round(est.mean(), 4),
              "cov1", res["seeds"][f"{N}|{s}"]["coverage_1sigmaF"], "chi2/dof", round(np.mean(chi2s) / nfree, 3),
              f"{time.perf_counter() - T0:.0f}s", flush=True)
        (OUT / "a5_profile.json").write_text(json.dumps(jsonable(res), indent=1))

# twin refit of 3 seeds (N2, 0.01 deg) as a check of the independent-model fits
chk = []
names = SETS_CAT["N2"]
th_names = [n for n in names if n not in INSTR]
sa, sp = 0.01, np.deg2rad(0.01)
Q = projector(design([n for n in names if n in INSTR], FREQ, sa, sp))
for i in range(3):
    dat = tw.synthetic_measurement(FREQ, st, WP, WR, sigma_log_amplitude=0.01, sigma_phase_deg=0.01,
                                   seed=20261111 + i, passivity="ignore")
    lnHd = np.log(np.abs(dat["observed"])) + 1j * np.unwrap(np.angle(dat["observed"]))

    def resid(x):
        d = lnHd - m0.lnH(FREQ, dict(zip(th_names, x)))
        d = np.log(np.abs(np.exp(d))) + 1j * np.angle(np.exp(d))
        return proj_out(whiten(d, sa, sp), Q)
    r = least_squares(resid, np.zeros(len(th_names)), method="lm", xtol=1e-14, ftol=1e-14, gtol=1e-14,
                      diff_step=1e-4)
    chk.append({"seed": 20261111 + i, "twin_hat": float(r.x[0]), "indep_hat": res["seeds"]["N2|0.01"]["estimates"][i],
                "twin_chi2": float(2 * r.cost)})
    print("twin refit", chk[-1], flush=True)
res["twin_refit_N2_0.01"] = chk
res["seconds"] = time.perf_counter() - T0
(OUT / "a5_profile.json").write_text(json.dumps(jsonable(res), indent=1))
print("done", res["seconds"])
