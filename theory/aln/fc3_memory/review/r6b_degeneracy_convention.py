"""Review R6b: does phono3py's degeneracy treatment explain the witness kappa discrepancy?
On the independent event operator (r5, conserving convention, orbit-averaged rates) apply phono3py-style
averaging: off-diagonal (in-scattering) part averaged over degenerate subspaces in rows and columns, diagonal
replaced by its degeneracy average; compare kappa ratios with the plain operator and with phono3py's LBTE.
usage: python -B r6b_degeneracy_convention.py <events.npz> <pp_ref.npz> <outjson> name=pp.npz ...
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

evp, ppref, outp = sys.argv[1], sys.argv[2], Path(sys.argv[3])
specs = sys.argv[4:]
saved = list(sys.argv)
sys.argv = ["r5", evp, ppref, str(outp.with_name("r6b_r5_dummy.json"))]
sp_ = importlib.util.spec_from_file_location("r5", str(Path(__file__).with_name("r5_event_operator.py")))
r5 = importlib.util.module_from_spec(sp_)
sp_.loader.exec_module(r5)
sys.argv = saved
n, nb = r5.n, r5.Z["freqs"].shape[1]
qv = np.nonzero(r5.valid)[0] // nb
f = r5.eps
# degenerate sets per q (tolerance 1e-6 THz, as the campaign)
grp = -np.ones(n, int)
g = 0
for q in np.unique(qv):
    ii = np.nonzero(qv == q)[0]
    ii = ii[np.argsort(f[ii])]
    st = 0
    for k in range(1, len(ii) + 1):
        if k == len(ii) or f[ii[k]] - f[ii[k - 1]] > 1e-6:
            grp[ii[st:k]] = g
            g += 1
            st = k
Pm = np.zeros((n, n))
for gg in range(g):
    ii = np.nonzero(grp == gg)[0]
    Pm[np.ix_(ii, ii)] = 1.0 / len(ii)


def kappa_tau(C):
    w, V = np.linalg.eigh(C)
    eh = r5.e_vec / np.linalg.norm(r5.e_vec)
    keep = np.ones(n, bool)
    keep[np.argmax(np.abs(V.T @ eh))] = False
    out = []
    for col in (0, 2):
        bc = V[:, keep].T @ r5.b_all[:, col]
        K = float(np.sum(bc ** 2 / w[keep]))
        N = float(np.sum(bc ** 2 / w[keep] ** 2))
        out += [r5.kf * K, N / K / (4 * np.pi)]
    return out


def p3style(C):
    D = np.diag(np.diag(C))
    Off = C - D
    return Pm @ Off @ Pm + np.diag(Pm @ np.diag(C))


res = {}
C0 = r5.operator(r5.rates(ppref), "cons")
base_plain, base_p3 = kappa_tau(C0), kappa_tau(p3style(C0))
res["REF_plain"], res["REF_p3style"] = base_plain, base_p3
for spec in specs:
    name, path = spec.split("=", 1)
    C = r5.operator(r5.rates(path), "cons")
    a, b = kappa_tau(C), kappa_tau(p3style(C))
    res[name] = {"plain_kappa_ratio": [a[0] / base_plain[0], a[2] / base_plain[2]],
                 "plain_tau_ratio": [a[1] / base_plain[1], a[3] / base_plain[3]],
                 "p3style_kappa_ratio": [b[0] / base_p3[0], b[2] / base_p3[2]],
                 "p3style_tau_ratio": [b[1] / base_p3[1], b[3] / base_p3[3]]}
    print(name, json.dumps(res[name]), flush=True)
outp.write_text(json.dumps(res, indent=1))
print("REF plain", base_plain, "p3style", base_p3)
