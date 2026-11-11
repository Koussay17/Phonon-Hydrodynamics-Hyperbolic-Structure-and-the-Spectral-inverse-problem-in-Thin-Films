"""Shared set-up of workstream A (FDTR detectability of the film memory-time gap).

Forward model: src/twin_fdtr.py (imported read-only). Stack Au(80 nm) / R1 / AlN(L) / R2 / sapphire, values of
note 16 and campaign assumption A5; AlN layer Cattaneo (tau_R = tau_mem,z(L), ell = 0) or Guyer-Krumhansl
(tau_ell nuisance, GKBoundary default: slip C_s = 1, tangential 'slip', temperature 'local').
Observable: H_obs = instrument(gaussian_response) (finite beams) or the 1D large-beam limit
2 Z(k=0, p) / (pi (wp^2 + wr^2)) of the same observable.

Parameters (all positive physical ones in logarithms):
  ln_tauR ln_taul ln_lamz ln_lamr ln_C (AlN) ln_R1 ln_R2 ln_lamAu ln_CAu ln_dAu ln_dAlN ln_lamS ln_CS
  ln_wscale (common factor on both radii) ln_wp ln_wr ; instrument g (log gain), phi0 (rad), td_ns (ns).
The instrument columns of d ln H_obs / d theta are exact: g -> 1, phi0 -> i, td_ns -> -i 2 pi f 1e-9.
The thermal columns are central differences in the log parameter with a FIXED quadrature rule (step, x_low), so
that the rule does not move between the perturbed evaluations.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "3")

import numpy as np

sys.path.insert(0, r"C:/Users/Koussay/these/src")
import twin_fdtr as tw  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
FIG = HERE.parent.parent / "figures"
FIG.mkdir(exist_ok=True)

# note 16 / assumption A5
AU_LAM, AU_C, AU_D = 150.0, 2.49e6, 80e-9
SAP_LAM, SAP_C = 35.0, 3.03e6
R1, R2 = 1e-8, 2e-8
C_ALN = 2.41e6
WP, WR = 10e-6, 8e-6
ALPHA = 1.0 / 3.0

FREQ = np.geomspace(1e4, 2e8, 80)
BANDS = {"10k-200M": (1e4, 2e8), "1M-200M": (1e6, 2e8), "10k-50M": (1e4, 5e7)}
SIG_A = 0.01                                  # on ln|H|
SIG_PHI_DEG = (0.01, 0.03, 0.1)
W_ADMISSIBLE = (0.028, 0.037, 0.052, 0.08)    # half-widths of 00-question.md
W_BY_L = {30: 0.028, 100: 0.037, 300: 0.052}  # orb 0.1, three data meshes (100 nm calibrated classes: 0.08)

STEP_Q, XLOW_Q = 0.04, -22.0                  # fixed trapezoidal rule in x = ln(sqrt(eta) k) for production

INSTR = ("g", "phi0", "td_ns")


def load_base():
    d = json.loads((OUT / "a1_base_values.json").read_text())
    return d["base"]


SETS_CAT = {
    "N0": ["ln_tauR"],
    "N1": ["ln_tauR", "g", "phi0", "td_ns"],
    "N2": ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2"],
    "N3": ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C"],
    "N4": ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C", "ln_lamAu", "ln_CAu"],
    "N5": ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C", "ln_lamAu", "ln_CAu",
           "ln_wscale"],
}
SETS_GK = {k: list(v) for k, v in SETS_CAT.items()}
for _k in ("N4", "N5"):
    SETS_GK[_k] = SETS_GK[_k][:8] + ["ln_taul"] + SETS_GK[_k][8:]
# sensitivity sets (not declared in advance; reported as such)
EXTRA = {
    "N1-td": ["ln_tauR", "g", "phi0"],                                  # N1 without the delay
    "N2-td": ["ln_tauR", "g", "phi0", "ln_R1", "ln_R2"],                # N2 without the delay
    "N3+lamr": SETS_CAT["N3"] + ["ln_lamr"],
    "N5+thick": SETS_CAT["N5"] + ["ln_dAu", "ln_dAlN"],
    "N4+radii(wp,wr)": SETS_CAT["N4"] + ["ln_wp", "ln_wr"],              # exactly degenerate: H depends on eta only
}


@dataclass
class Model:
    """Base physical parameters and the forward map theta -> ln H_obs."""
    L_nm: float
    tauR: float
    lamz: float
    lamr: float
    law: str = "cattaneo"           # "cattaneo" or "gk"
    taul: float = 0.0
    geometry: str = "beam"          # "beam" (Gaussian average) or "1d" (k = 0 limit, same normalisation)
    wp: float = WP
    wr: float = WR
    step: float = STEP_Q
    x_low: float = XLOW_Q
    quad: str = "trapezoid"         # or "adaptive"
    epsrel: float = 1e-12
    extra: dict = field(default_factory=dict)

    def phys(self, theta: dict | None = None):
        """Physical parameters after applying log/linear perturbations theta (dict name -> increment)."""
        t = dict(theta or {})
        e = lambda n: float(np.exp(t.get(n, 0.0)))   # noqa: E731
        P = {
            "tauR": self.tauR * e("ln_tauR"), "taul": self.taul * e("ln_taul"),
            "lamz": self.lamz * e("ln_lamz"), "C": C_ALN * e("ln_C"),
            "R1": R1 * e("ln_R1"), "R2": R2 * e("ln_R2"),
            "lamAu": AU_LAM * e("ln_lamAu"), "CAu": AU_C * e("ln_CAu"), "dAu": AU_D * e("ln_dAu"),
            "dAlN": self.L_nm * 1e-9 * e("ln_dAlN"), "lamS": SAP_LAM * e("ln_lamS"), "CS": SAP_C * e("ln_CS"),
            "wp": self.wp * e("ln_wscale") * e("ln_wp"), "wr": self.wr * e("ln_wscale") * e("ln_wr"),
            "g": t.get("g", 0.0), "phi0": t.get("phi0", 0.0), "td_ns": t.get("td_ns", 0.0),
        }
        if self.law == "gk":
            P["lamr"] = P["lamz"]                    # GK layers of the twin are isotropic
        else:
            P["lamr"] = self.lamr * e("ln_lamr") * (e("ln_lamz") if self.extra.get("lamr_follows_lamz") else 1.0)
        return P

    def stack(self, P):
        if self.law == "gk":
            ell = float(np.sqrt(P["taul"] * P["lamz"] / ((1 + ALPHA) * P["C"])))
            aln = tw.Material(P["lamz"], P["C"], tau_r=P["tauR"], ell=ell, alpha=ALPHA)
            bc = tw.GKBoundary(slip=self.extra.get("slip", 1.0), temperature=self.extra.get("temperature", "local"))
        else:
            aln = tw.Material(P["lamz"], P["C"], lam_r=P["lamr"], tau_r=P["tauR"])
            bc = None
        return tw.au_aln_sapphire(aln, aln_thickness=P["dAlN"], au_thickness=P["dAu"],
                                  au=tw.Material(P["lamAu"], P["CAu"]),
                                  sapphire=tw.Material(P["lamS"], P["CS"]),
                                  r_au_aln=P["R1"], r_aln_sapphire=P["R2"], boundary=bc)

    def thermal(self, f, P):
        st = self.stack(P)
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            if self.geometry == "1d":
                h = tw.kernel(f, 0.0, st) * 2.0 / (np.pi * (P["wp"] ** 2 + P["wr"] ** 2))
            elif self.quad == "adaptive":
                h = tw.gaussian_response_adaptive(f, st, P["wp"], P["wr"], epsrel=self.epsrel)
            else:
                h = tw.gaussian_response(f, st, P["wp"], P["wr"], step=self.step, x_low=self.x_low)
        self._warnings = [str(w.message)[:120] for w in wl]
        return h

    def lnH(self, f, theta=None):
        """ln H_obs (complex; imaginary part = phase in rad, continuous across f)."""
        P = self.phys(theta)
        h = self.thermal(f, P)
        ln = np.log(np.abs(h)) + 1j * np.unwrap(np.angle(h))
        return ln + P["g"] + 1j * (P["phi0"] - 2 * np.pi * f * P["td_ns"] * 1e-9)


def jacobian(model: Model, f, names, h=1e-3, scheme="central"):
    """Complex Jacobian d ln H_obs / d theta, shape (len(f), len(names)).

    scheme: 'central' (O(h^2)) or 'richardson' (central at h and h/2, O(h^4)).
    Thermal log-parameters by finite differences of ln(H(+)/H(-)) (no branch issue); instrument columns exact.
    """
    f = np.asarray(f, float)
    cols = []
    for n in names:
        if n == "g":
            cols.append(np.ones(f.size, complex))
        elif n == "phi0":
            cols.append(np.full(f.size, 1j))
        elif n == "td_ns":
            cols.append(-1j * 2 * np.pi * f * 1e-9)
        else:
            def cd(hh):
                hp = model.thermal(f, model.phys({n: +hh}))
                hm = model.thermal(f, model.phys({n: -hh}))
                return np.log(hp / hm) / (2 * hh)
            if scheme == "central":
                cols.append(cd(h))
            elif scheme == "richardson":
                d1, d2 = cd(h), cd(h / 2)
                cols.append((4 * d2 - d1) / 3)
            else:
                raise ValueError(scheme)
    return np.stack(cols, axis=1)


def weighted(J, sig_a=SIG_A, sig_phi_deg=0.1, mask=None):
    """Whitened real Jacobian [Re J / sig_a ; Im J / sig_phi] for the selected frequencies."""
    if mask is not None:
        J = J[mask]
    return np.vstack([J.real / sig_a, J.imag / np.deg2rad(sig_phi_deg)])


def fisher(Jw, target=0, cond_max=1e10):
    """Marginal standard deviation of parameter `target` from the SVD of the column-normalised whitened Jacobian.

    Returns sigma (None if flagged), singular values (normalised), condition number, the right singular vector
    of the smallest singular value, and the conditional sigma (all nuisances fixed).
    """
    nrm = np.linalg.norm(Jw, axis=0)
    if np.any(nrm == 0):
        return {"flag": "zero column", "sigma": None}
    Jn = Jw / nrm
    U, S, Vt = np.linalg.svd(Jn, full_matrices=False)
    cond = float(S[0] / S[-1]) if S[-1] > 0 else np.inf
    out = {"sv": S.tolist(), "cond": cond, "v_min": Vt[-1].tolist(),
           "sigma_conditional": float(1.0 / nrm[target])}
    if cond > cond_max:
        out["flag"] = f"numerically unidentifiable (cond {cond:.3g} > {cond_max:.0e})"
        out["sigma"] = None
        return out
    var = np.sum((Vt[:, target] / S) ** 2) / nrm[target] ** 2
    out["sigma"] = float(np.sqrt(var))
    out["flag"] = ""
    # variance inflation and the parameter most correlated with the target
    cov = (Vt.T / S ** 2) @ Vt / np.outer(nrm, nrm)
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    out["corr_target"] = corr[target].tolist()
    return out


def band_mask(f, band):
    lo, hi = BANDS[band]
    return (f >= lo * (1 - 1e-12)) & (f <= hi * (1 + 1e-12))


def model_for(L, law="cattaneo", taul_factor=0.0, geometry="beam", base=None, lamr_mode="rta", **kw):
    """Model at thickness L (nm) with the adopted base values (event operator; 500 nm by the ratio method).

    lamr_mode: 'rta' (in-plane film value kappa_x(L), RTA 64x64x35), 'iso' (lam_r = lam_z), 'event12' (12x12x8).
    """
    base = base or load_base()
    b = base[str(int(L))]
    lamz, tau = b["lam_z"], b["tau_R_ps"] * 1e-12
    if lamr_mode == "rta":
        lamr = base[str(int(L)) + "_rta"]["lam_r"]
    elif lamr_mode == "iso":
        lamr = lamz
    elif lamr_mode == "event12":
        lamr = b["lam_r"]
    elif lamr_mode == "xextrap":
        lamr = base[str(int(L)) + "_xextrap"]["lam_r"]
    else:
        raise ValueError(lamr_mode)
    return Model(L_nm=float(L), tauR=tau, lamz=lamz, lamr=lamr, law=law, taul=taul_factor * tau,
                 geometry=geometry, **kw)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    return o
