"""AlN event geometry (from aln_events.py export) as an EventGeometry, plus validation.

Modes: full-grid (q, j) with frequency above phono3py's cutoff (the three acoustic Gamma modes are
excluded). Energies in THz (eps = f), kT = k_B T / THzToEv (THz). Velocities in THz*Angstrom.
Event rates g in THz (phono3py collision units); physical rates are 4*pi*g (1/ps), so
tau_mem[ps] = (N/K)/(4 pi) with N, K computed from the THz-valued operator.
kappa_aa [W/(m K)] = unit_to_WmK/volume * THzToEv^2 * K_aa / (2 N_q k_B T^2).
"""
from __future__ import annotations

import json
import sys

import numpy as np

from memgap import EventGeometry


def load(path, nsig_max=None):
    Z = np.load(path)
    f = Z["freqs"]
    Nq, nb = f.shape
    cutoff = float(Z["cutoff"])
    valid = (f > cutoff).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(valid.sum())
    eps = f.ravel()[valid]
    T = float(Z["T"])
    kT = float(Z["KB"]) * T / float(Z["THzToEv"])
    vel = Z["gv"].reshape(Nq * nb, 3)[valid]
    p, a, b = Z["ev_p"], Z["ev_a"], Z["ev_b"]
    g = Z["ev_g"]
    dl = Z["ev_dl"]
    keep = np.ones(len(p), dtype=bool)
    if nsig_max is not None:
        keep &= np.abs(dl) <= nsig_max * float(Z["sigma"])
    p, a, b, g, dl = p[keep], a[keep], b[keep], g[keep], dl[keep]
    umk = Z["ev_umk"][keep]
    P, A, B = newidx[p], newidx[a], newidx[b]
    assert (P >= 0).all() and (A >= 0).all() and (B >= 0).all()
    rep = A == B
    fp, fa, fb = eps[P], eps[A], eps[B]
    lp = np.sqrt((fa + fb) / fp)
    ld = 1.0 / lp
    coef = np.zeros((len(P), 3))
    coef[:, 0] = -lp
    coef[:, 1] = np.where(rep, 2 * ld, ld)
    coef[:, 2] = np.where(rep, 0.0, ld)
    ev = np.stack([P, A, np.where(rep, -1, B)], axis=1)
    # mode-level symmetry maps (rotations incl. time reversal), restricted to valid modes
    rot = Z["rot_maps"]                      # (nrot, Nq) grid maps
    maps = []
    for R in rot:
        full = (R[:, None] * nb + np.arange(nb)[None, :]).ravel()
        mp = newidx[full][valid]
        assert (mp >= 0).all()
        maps.append(mp)
    geom = EventGeometry(eps=eps, kT=kT, vel=vel, events=ev, coef=coef, kind=umk.astype(np.int64), gphys=g,
                         labels={"mesh": Z["mesh"].tolist(), "sigma": float(Z["sigma"]), "T": T,
                                 "Nq": Nq, "nb": nb, "valid": valid, "newidx": newidx, "maps": maps,
                                 "delta": dl, "unit_to_WmK": float(Z["unit_to_WmK"]),
                                 "volume": float(Z["volume"]), "THzToEv": float(Z["THzToEv"]),
                                 "KB": float(Z["KB"]), "gam_raw": Z["gam_raw"].ravel()[valid]})
    return geom


def kappa_factor(geom):
    L = geom.labels
    return L["unit_to_WmK"] / L["volume"] * L["THzToEv"] ** 2 / (2 * L["Nq"] * L["KB"] * L["T"] ** 2)


def tau_ps(tau_thz_units):
    return tau_thz_units / (4 * np.pi)


if __name__ == "__main__":
    path = sys.argv[1]
    geom = load(path)
    g = geom.gphys
    r = geom.W @ g
    gam = geom.labels["gam_raw"]
    rel = (r - gam) / gam
    out = {"n": geom.n, "m": geom.m, "umklapp_fraction": float(geom.kind.mean()),
           "repeated": int((geom.events[:, 2] < 0).sum()),
           "diag_vs_phono3py_gamma_rel": {"median": float(np.median(np.abs(rel))), "max": float(np.abs(rel).max()),
                                           "mean_signed": float(rel.mean())}}
    C = geom.C(g)
    out["energy_residual"] = float(np.linalg.norm(C @ geom.e) / np.linalg.norm(np.diag(C) * geom.e))
    kf = kappa_factor(geom)
    for col, name in ((0, "x"), (2, "z")):
        R = geom.response(g, col)
        b = geom.b[:, col]
        KR = float(np.sum(b * b / r))
        NR = float(np.sum(b * b / r ** 2))
        Kg = float(np.sum(b * b / gam))
        out[name] = {"K": R["K"], "kappa_event_WmK": kf * R["K"], "kappa_RTA_diag_WmK": kf * KR,
                     "kappa_RTA_phono3py_gamma_WmK": kf * Kg,
                     "tau_mem_ps": tau_ps(R["tau"]), "tau_RTA_ps": tau_ps(NR / KR)}
    print(json.dumps(out, indent=1))
