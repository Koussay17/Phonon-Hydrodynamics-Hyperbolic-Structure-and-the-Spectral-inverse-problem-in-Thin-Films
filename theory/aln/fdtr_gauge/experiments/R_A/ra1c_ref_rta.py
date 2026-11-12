"""R-A1c: RTA kinetic times of the REF tensor on the event meshes (diagonal of the event operator), for the mesh
correction of tau_l (A-2): the full-operator tau_l is available on 7x7x5 and 9x9x6 only, the RTA tau_l also on the
converged phono3py mesh 64x64x35 (ra2).

Diagonal r_i = sum_events g a_i^2 (= diag C(g)), g = W pp of the REF amplitudes WITHOUT orbit averaging (the reference
rates are symmetric; checked on 7x7x5 against the orbit-tied values of ra1b_dense_check.json). Energy-conserving RTA
operator R - (R e)(R e)^T/(e^T R e), R = 4 pi (r + d), d = c |v_z| / (4 pi L 10) (film_event units):
  tau_R,RTA = sum x^2 / K, x = b / R;  tau_l,RTA = |e|^2 sum (Q V x)_i^2 / R_i / K^2.
Light: no dense matrix, no orbit computation. Output: out/ra1c_ref_rta.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
FOURPI = 4 * np.pi
MESH = {"m775": ("results_c/m775x/events_m775_s0.1.npz", "results_c/m775x/film2_pp"),
        "m996": ("results_c/m996/events_m996_s0.1.npz", "results_c/m996/film2_pp")}
sys.path.insert(0, str(HERE))
from ra1_film_event_ext import LS_NM, lkey  # noqa: E402


def ref_diag(mesh):
    Z = np.load(FC3 / MESH[mesh][0])
    f = Z["freqs"]
    Nq, nb = f.shape
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(valid.sum())
    eps = f.ravel()[valid]
    kT = float(Z["KB"]) * float(Z["T"]) / float(Z["THzToEv"])
    vel = Z["gv"].reshape(Nq * nb, 3)[valid]
    P, A, B = newidx[Z["ev_p"]], newidx[Z["ev_a"]], newidx[Z["ev_b"]]
    rep = A == B
    fp, fa, fb = eps[P], eps[A], eps[B]
    lp = np.sqrt((fa + fb) / fp)
    ld = 1.0 / lp
    x = eps / kT
    nocc = 1.0 / np.expm1(x)
    dsc = np.sqrt(nocc * (1 + nocc))
    s = np.sinh(x / 2)
    pp = np.load(FC3 / MESH[mesh][1] / "pp_r01_REF.npz")
    conv, g0w = float(pp["conv"]), pp["g0w"]
    Wv = np.where(rep, conv * g0w / (8 * s[P] * s[A] * s[A]), conv * g0w / (4 * s[P] * s[A] * s[np.where(rep, A, B)]))
    g = Wv * pp["pp"]
    n = eps.size
    r = np.bincount(P, weights=g * (lp / dsc[P]) ** 2, minlength=n)
    r += np.bincount(A, weights=g * (np.where(rep, 2 * ld, ld) / dsc[A]) ** 2, minlength=n)
    r += np.bincount(B[~rep], weights=g[~rep] * (ld[~rep] / dsc[B[~rep]]) ** 2, minlength=n)
    e = dsc * eps
    return r, e, vel[:, 2]


res = {}
for mesh in ("m775", "m996"):
    r, e, vz = ref_diag(mesh)
    ehat = e / np.linalg.norm(e)
    b = e * vz
    row = {}
    for L in LS_NM:
        R = FOURPI * (r + np.abs(vz) / (FOURPI * L * 10.0))
        x = b / R
        K = float(b @ x)
        u = vz * x
        u = u - ehat * float(ehat @ u)
        row[lkey(L)] = {"tau_RTA": float(x @ x) / K, "taul_RTA": float(e @ e) * float(np.sum(u * u / R)) / K ** 2}
    res[mesh] = row
    print(mesh, {k: (round(v["tau_RTA"], 3), round(v["taul_RTA"], 2)) for k, v in row.items()
                 if k in ("30", "100", "300", "500", "1000")}, flush=True)
chk = json.loads((OUT / "ra1b_dense_check.json").read_text())
res["check_vs_orbit_tied_m775"] = max(max(abs(res["m775"][k.split("|")[1]]["tau_RTA"] / v["tau_RTA_ec"] - 1),
                                          abs(res["m775"][k.split("|")[1]]["taul_RTA"] / v["taul_RTA_ec"] - 1))
                                      for k, v in chk.items() if k.startswith("REF|"))
print("check vs orbit-tied (7x7x5):", res["check_vs_orbit_tied_m775"])
(OUT / "ra1c_ref_rta.json").write_text(json.dumps(res, indent=1))
