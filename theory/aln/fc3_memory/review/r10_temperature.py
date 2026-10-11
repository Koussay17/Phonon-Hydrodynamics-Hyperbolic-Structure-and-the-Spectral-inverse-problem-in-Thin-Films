"""Review R10: temperature transfer of the witnesses (constraints were imposed at 300 K only). Event operator
(conserving convention, orbit-averaged rates, independent implementation as r5) rebuilt at temperature T with
the same events, Gaussian weights and |Phi|^2 (T independent): W(T) = conv g0w / (4 s_p s_a s_b), s = sinh(x/2),
entropy scales, current and kappa prefactor at T. Ratios witness / reference of kappa_xx, kappa_zz, tau_x, tau_z.
usage: python -B r10_temperature.py <events.npz> <pp_ref.npz> <outjson> <T1,T2,...> name=pp.npz ...
"""
import json
import sys
from pathlib import Path

import numpy as np

evp, ppref, outp, Ts = sys.argv[1], sys.argv[2], Path(sys.argv[3]), [float(t) for t in sys.argv[4].split(",")]
specs = sys.argv[5:]
Z = np.load(evp)
f = Z["freqs"]
Nq, nb = f.shape
valid = (f > float(Z["cutoff"])).ravel()
idx = -np.ones(Nq * nb, int)
idx[valid] = np.arange(valid.sum())
n = int(valid.sum())
eps = f.ravel()[valid]
vel = Z["gv"].reshape(Nq * nb, 3)[valid]
P, A, Bm = idx[Z["ev_p"]], idx[Z["ev_a"]], idx[Z["ev_b"]]
rep = A == Bm
m = len(P)
rot = Z["rot_maps"]
maps = [idx[(R[:, None] * nb + np.arange(nb)[None, :]).ravel()][valid] for R in rot]
key = {(P[i], min(A[i], Bm[i]), max(A[i], Bm[i])): i for i in range(m)}
parent = np.arange(m)


def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


for mp in maps:
    for i in range(m):
        a_, b_ = mp[A[i]], mp[Bm[i]]
        j = key[(mp[P[i]], min(a_, b_), max(a_, b_))]
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
_, orb = np.unique(np.array([find(i) for i in range(m)]), return_inverse=True)
osz = np.bincount(orb)
lp = np.sqrt((eps[A] + eps[Bm]) / eps[P])
ld = 1 / lp
coef = np.stack([-lp, np.where(rep, 2 * ld, ld), np.where(rep, 0.0, ld)], 1)


def at_T(T):
    kT = float(Z["KB"]) * T / float(Z["THzToEv"])
    nocc = 1 / np.expm1(eps / kT)
    d = np.sqrt(nocc * (nocc + 1))
    s = np.sinh(eps / (2 * kT))
    kf = float(Z["unit_to_WmK"]) / float(Z["volume"]) * float(Z["THzToEv"]) ** 2 / (2 * Nq * float(Z["KB"]) * T ** 2)
    return d, s, kf, d * eps, (d * eps)[:, None] * vel


def metrics(ppfile, T):
    d, s, kf, e, b = at_T(T)
    Pz = np.load(ppfile)
    conv = float(Pz["conv"])
    W = np.where(rep, conv * Pz["g0w"] / (8 * s[P] * s[A] ** 2), conv * Pz["g0w"] / (4 * s[P] * s[A] * s[Bm]))
    g = W * Pz["pp"]
    g = (np.bincount(orb, weights=g) / osz)[orb]
    C = np.zeros((n, n))
    ids = [P, A, Bm]
    for k in range(3):
        for l in range(3):
            w = g * coef[:, k] * coef[:, l] / (d[ids[k]] * d[ids[l]])
            if k == 2 or l == 2:
                w = np.where(rep, 0.0, w)
            np.add.at(C, (ids[k], ids[l]), w)
    wv, V = np.linalg.eigh(C)
    eh = e / np.linalg.norm(e)
    keep = np.ones(n, bool)
    keep[np.argmax(np.abs(V.T @ eh))] = False
    out = {}
    for col, nm in ((0, "x"), (2, "z")):
        bc = V[:, keep].T @ b[:, col]
        K = float(np.sum(bc ** 2 / wv[keep]))
        N = float(np.sum(bc ** 2 / wv[keep] ** 2))
        out[f"kappa_{nm}"], out[f"tau_{nm}"] = kf * K, N / K / (4 * np.pi)
    return out


res = {"T": Ts}
for T in Ts:
    r0 = metrics(ppref, T)
    res[f"REF_{T:g}"] = r0
    for spec in specs:
        name, path = spec.split("=", 1)
        r = metrics(path, T)
        res.setdefault(name, {})[f"{T:g}"] = {"kappa_ratio": [r["kappa_x"] / r0["kappa_x"], r["kappa_z"] / r0["kappa_z"]],
                                             "tau_ratio": [r["tau_x"] / r0["tau_x"], r["tau_z"] / r0["tau_z"]]}
    print(T, "REF", r0, flush=True)
outp.write_text(json.dumps(res, indent=1))
for spec in specs:
    name = spec.split("=", 1)[0]
    print(name, " ".join(f"T{t}: k {v['kappa_ratio'][0]:.4f}/{v['kappa_ratio'][1]:.4f} tau {v['tau_ratio'][0]:.3f}/{v['tau_ratio'][1]:.3f}"
                         for t, v in res[name].items()))
