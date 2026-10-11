"""Numerical audit, part 7: infrared structure of the reference operator (S8).

For each event file: orbit-tied reference rates, r = diag C (mode relaxation rate in the operator's normalisation),
b = D (f v). Acoustic modes near Gamma (lowest three bands, f < f_max) are fitted by r = A f^p (per mesh and pooled
over meshes); the RTA analogue N_RTA = sum b^2 / r^2 is split into the contribution of the Gamma-adjacent shells.
Heuristic used in the audit: in three dimensions, with b ~ kT v (equipartition) and r ~ f^p ~ q^p near Gamma, the
continuum RTA integral N ~ int q^2 dq q^(-2p) diverges at small q when p >= 3/2 (as q_min^(3 - 2p)), whereas
K ~ int q^2 dq q^(-p) converges for p < 3. A mesh then cuts the divergence at q_min ~ 1/(mesh division).
usage: python -B na_ir.py out.json events1.npz [events2.npz ...]
"""
import json
import sys

import numpy as np


def analyse(path, fmax=4.0):
    Z = np.load(path)
    f = Z["freqs"]
    Nq, nb = f.shape
    valid = (f > float(Z["cutoff"])).ravel()
    newidx = -np.ones(Nq * nb, dtype=np.int64)
    newidx[valid] = np.arange(int(valid.sum()))
    fv = f.ravel()[valid]
    full = np.nonzero(valid)[0]
    band = full % nb
    q_of = full // nb
    T = float(Z["T"])
    kT = float(Z["KB"]) * T / float(Z["THzToEv"])
    nocc = 1 / np.expm1(fv / kT)
    d = np.sqrt(nocc * (1 + nocc))
    v = Z["gv"].reshape(Nq * nb, 3)[valid]
    b = (d * fv)[:, None] * v
    ev_p, ev_a, ev_b, g = Z["ev_p"], Z["ev_a"], Z["ev_b"], Z["ev_g"]
    m = len(g)
    M = Nq * nb
    key = (ev_p * M + ev_a) * M + ev_b
    order = np.argsort(key)
    skey = key[order]
    acc = np.zeros(m)
    for Rm in Z["rot_maps"]:
        p2 = Rm[ev_p // nb] * nb + ev_p % nb
        a2 = Rm[ev_a // nb] * nb + ev_a % nb
        b2 = Rm[ev_b // nb] * nb + ev_b % nb
        k2 = (p2 * M + np.minimum(a2, b2)) * M + np.maximum(a2, b2)
        acc += g[order[np.minimum(np.searchsorted(skey, k2), m - 1)]]
    g = acc / len(Z["rot_maps"])
    P, A, B = newidx[ev_p], newidx[ev_a], newidx[ev_b]
    rep = A == B
    lp = np.sqrt((fv[A] + fv[B]) / fv[P])
    r = (np.bincount(P, weights=g * (lp / d[P]) ** 2, minlength=len(fv))
         + np.bincount(A, weights=g * (np.where(rep, 2 / lp, 1 / lp) / d[A]) ** 2, minlength=len(fv))
         + np.bincount(B, weights=g * np.where(rep, 0.0, (1 / lp / d[B]) ** 2), minlength=len(fv)))
    ac = (band < 3) & (fv < fmax)
    out = {"mesh": Z["mesh"].tolist(), "n_acoustic_below_fmax": int(ac.sum()), "fmax": fmax}
    if ac.sum() >= 3:
        p, c = np.polyfit(np.log(fv[ac]), np.log(r[ac]), 1)
        out["fit_r_vs_f_acoustic"] = {"p": float(p), "A": float(np.exp(c))}
    out["lowest_acoustic"] = sorted([(float(fv[i]), float(r[i]), int(band[i]), Z["addr"][q_of[i]].tolist())
                                     for i in np.nonzero(ac)[0]])[:12]
    # RTA analogue and shell split (shell = distinct |q| of the Gamma-adjacent points, in mesh units)
    addr = Z["addr"][q_of] / Z["mesh"]
    rec = np.linalg.inv(Z["lattice"]).T * 2 * np.pi
    qn = np.linalg.norm(addr @ rec, axis=1)
    for c_, a in ((0, "x"), (2, "z")):
        w = b[:, c_] ** 2 / r ** 2
        NR, KR = float(w.sum()), float(np.sum(b[:, c_] ** 2 / r))
        qs = np.unique(np.round(qn[qn > 0], 6))[:3]
        out[f"RTA_{a}"] = {"tau_ps": NR / KR / (4 * np.pi), "N_over_Nq": NR / Nq, "K_over_Nq": KR / Nq,
                           "N_share_q_below_0.5": float(w[qn < 0.5].sum() / NR),
                           "N_share_first_shells": [float(w[np.isclose(qn, q0, atol=1e-6)].sum() / NR) for q0 in qs],
                           "first_shell_q": [float(t) for t in qs]}
    return out


def main():
    res = {}
    for pth in sys.argv[2:]:
        o = analyse(pth)
        res["x".join(map(str, o["mesh"]))] = o
        print(o["mesh"], o.get("fit_r_vs_f_acoustic"), {a: (round(o[f"RTA_{a}"]["tau_ps"], 3), round(o[f"RTA_{a}"]["N_share_q_below_0.5"], 3))
                                                        for a in "xz"}, flush=True)
    # pooled fit over meshes of the Gamma-adjacent acoustic modes
    fs, rs = [], []
    for o in res.values():
        for f_, r_, _, _ in o["lowest_acoustic"]:
            fs.append(f_)
            rs.append(r_)
    if len(fs) > 3:
        p, c = np.polyfit(np.log(fs), np.log(rs), 1)
        res["pooled_fit_lowest_acoustic"] = {"p": float(p), "A": float(np.exp(c)), "n": len(fs)}
        print("pooled", res["pooled_fit_lowest_acoustic"])
    with open(sys.argv[1], "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
