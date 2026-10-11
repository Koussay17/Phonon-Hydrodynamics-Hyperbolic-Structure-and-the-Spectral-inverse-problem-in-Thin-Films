"""Review R5: independent implementation of the event operator from the written definitions (02-assumptions A, B),
no campaign code imported: own orbit tying (rot_maps of the event file), own assembly (dense), own eigensolver.
Variants of the discretisation convention (same fc3, same events, same pp):
  'cons' : conserving rescaled stoichiometry (-lp, ld, ld), lp = sqrt((fa + fb)/fp)  (campaign convention)
  'raw'  : raw stoichiometry (-1, 1, 1) (repeated daughter (-1, 2)); C e != 0, inverse on the odd sector
Orbit tying: average of per-event rates over the 24-operation event orbits (memory-gap convention).
Outputs kappa_xx, kappa_zz (W/mK), tau_x, tau_z (ps), r = diag C, for REF and each pp file; ratios to REF.
Research env. usage: python -B r5_event_operator.py <events.npz> <pp_ref.npz> <outjson> name=pp.npz ...
"""
import json
import sys
from pathlib import Path

import numpy as np

evp, ppref, outp = sys.argv[1], sys.argv[2], Path(sys.argv[3])
specs = sys.argv[4:]
Z = np.load(evp)
f = Z["freqs"]
Nq, nb = f.shape
valid = (f > float(Z["cutoff"])).ravel()
idx = -np.ones(Nq * nb, int)
idx[valid] = np.arange(valid.sum())
n = int(valid.sum())
eps = f.ravel()[valid]
T = float(Z["T"])
kT = float(Z["KB"]) * T / float(Z["THzToEv"])
nocc = 1.0 / np.expm1(eps / kT)
d = np.sqrt(nocc * (nocc + 1))
vel = Z["gv"].reshape(Nq * nb, 3)[valid]
P, A, Bm = idx[Z["ev_p"]], idx[Z["ev_a"]], idx[Z["ev_b"]]
assert (P >= 0).all() and (A >= 0).all() and (Bm >= 0).all()
rep = A == Bm
m = len(P)
s = np.sinh(eps / (2 * kT))

# orbits of events under the 24 operations (mode map (q, j) -> (R q, j))
rot = Z["rot_maps"]
full_maps = [idx[(R[:, None] * nb + np.arange(nb)[None, :]).ravel()][valid] for R in rot]
key = {}
canon = lambda p, a, b: (p, min(a, b), max(a, b))
for i in range(m):
    key[canon(P[i], A[i], Bm[i])] = i
parent = np.arange(m)


def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


for mp in full_maps:
    for i in range(m):
        j = key[canon(mp[P[i]], mp[A[i]], mp[Bm[i]])]
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
roots = np.array([find(i) for i in range(m)])
_, orb = np.unique(roots, return_inverse=True)
osz = np.bincount(orb)

e_vec = d * eps
b_all = (d * eps)[:, None] * vel
kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (2 * Nq * float(Z["KB"]) * T ** 2)


def coefs(conv):
    fp, fa, fb = eps[P], eps[A], eps[Bm]
    if conv == "cons":
        lp = np.sqrt((fa + fb) / fp)
        ld = 1 / lp
    else:
        lp = np.ones(m)
        ld = np.ones(m)
    c = np.zeros((m, 3))
    c[:, 0] = -lp
    c[:, 1] = np.where(rep, 2 * ld, ld)
    c[:, 2] = np.where(rep, 0.0, ld)
    return c


def operator(g, conv):
    c = coefs(conv)
    C = np.zeros((n, n))
    ids = [P, A, Bm]
    for k in range(3):
        for l in range(3):
            w = g * c[:, k] * c[:, l] / (d[ids[k]] * d[ids[l]])
            if k == 2 or l == 2:
                w = np.where(rep, 0.0, w)
            np.add.at(C, (ids[k], ids[l]), w)
    return C


def metrics(g, conv):
    C = operator(g, conv)
    w, V = np.linalg.eigh(C)
    eh = e_vec / np.linalg.norm(e_vec)
    out = {"r": np.diag(C).copy(), "energy_residual": float(np.linalg.norm(C @ eh) / np.linalg.norm(np.diag(C)))}
    if conv == "cons":
        keep = np.ones(n, bool)
        keep[np.argmax(np.abs(V.T @ eh))] = False
    else:
        keep = np.ones(n, bool)
    for col, nm in ((0, "x"), (2, "z")):
        bc = V[:, keep].T @ b_all[:, col]
        K = float(np.sum(bc ** 2 / w[keep]))
        N = float(np.sum(bc ** 2 / w[keep] ** 2))
        out[f"kappa_{nm}"] = kf * K
        out[f"tau_{nm}"] = N / K / (4 * np.pi)
    out["lam_min"] = float(np.sort(w[keep])[0])
    return out


def rates(ppfile):
    Pz = np.load(ppfile)
    conv = float(Pz["conv"])
    g0w, pp = Pz["g0w"], Pz["pp"]
    W = np.where(rep, conv * g0w / (8 * s[P] * s[A] ** 2), conv * g0w / (4 * s[P] * s[A] * s[Bm]))
    g = W * pp
    return (np.bincount(orb, weights=g) / osz)[orb]


res = {"n_modes": n, "n_events": m, "n_orbits": int(orb.max() + 1)}
g_ref = rates(ppref)
res["ref_rate_check_vs_event_file"] = float(np.max(np.abs(g_ref - (np.bincount(orb, weights=Z["ev_g"]) / osz)[orb])) / g_ref.max())
for conv in ("cons", "raw"):
    m0 = metrics(g_ref, conv)
    res[f"REF_{conv}"] = {k: v for k, v in m0.items() if k != "r"}
    for spec in specs:
        name, path = spec.split("=", 1)
        mm = metrics(rates(path), conv)
        lr = mm["r"] / m0["r"] - 1
        res[f"{name}_{conv}"] = {"kappa_ratio": [mm["kappa_x"] / m0["kappa_x"], mm["kappa_z"] / m0["kappa_z"]],
                                 "tau_ratio": [mm["tau_x"] / m0["tau_x"], mm["tau_z"] / m0["tau_z"]],
                                 "tau_ps": [mm["tau_x"], mm["tau_z"]],
                                 "lifetime_raw_minmax": [float(lr.min()), float(lr.max())],
                                 "energy_residual": mm["energy_residual"]}
        print(name, conv, json.dumps(res[f"{name}_{conv}"]), flush=True)
outp.write_text(json.dumps(res, indent=1))
print(json.dumps({k: v for k, v in res.items() if k.startswith("REF") or k.startswith("n_") or k.startswith("ref")}, indent=1))
