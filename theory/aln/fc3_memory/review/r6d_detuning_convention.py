"""Review R6d: detuned-event (finite sigma) weighting as the origin of the witness kappa discrepancy.
Event operator variant 'p3form' built from the same events, pp and Gaussian weights, with phono3py's Bose-factor
forms for detuned events (identical to the campaign's 'cons' operator when an event conserves energy exactly):
  off-diagonal (p,a), (p,b): -conv g0w pp / s_b, -conv g0w pp / s_a      (same as 'cons')
  off-diagonal (a,b):        +conv g0w pp / s_p                          ('cons': times ld^2)
  diag p: conv g0w pp sinh((x_a + x_b)/2) / (s_a s_b)                    ('cons': lp^2 s_p / (s_a s_b))
  diag a: conv g0w pp sinh((x_p - x_b)/2) / (s_p s_b), diag b analogous ('cons': ld^2 s_a / (s_p s_b))
  repeated daughter p -> a + a: off (p,a) -conv g0w pp / s_a; diag p conv g0w pp sinh(x_a) / (2 s_a^2);
                                diag a 2 conv g0w pp sinh((x_p - x_a)/2) / (s_p s_a)
Rates orbit-averaged as in r5 (W pp averaged over event orbits before the Bose forms are applied: the orbit
average is taken on pp g0w, which is what the orbit tie acts on). Inverse on the odd sector (pinv).
usage: python -B r6d_detuning_convention.py <events.npz> <pp_ref.npz> <outjson> name=pp.npz ...
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

evp, ppref, outp = sys.argv[1], sys.argv[2], Path(sys.argv[3])
specs = sys.argv[4:]
saved = list(sys.argv)
sys.argv = ["r5", evp, ppref, str(outp.with_name("r6d_r5_dummy.json"))]
sp_ = importlib.util.spec_from_file_location("r5", str(Path(__file__).with_name("r5_event_operator.py")))
r5 = importlib.util.module_from_spec(sp_)
sp_.loader.exec_module(r5)
sys.argv = saved
n, P, A, Bm, rep, s, eps, kT = r5.n, r5.P, r5.A, r5.Bm, r5.rep, r5.s, r5.eps, r5.kT
x = eps / kT
Z = r5.Z
nb = Z["freqs"].shape[1]
Jv = r5.idx[(Z["Jmap"][:, None] * nb + np.arange(nb)[None, :]).ravel()][r5.valid]
Podd = 0.5 * (np.eye(n) - np.eye(n)[Jv])


def amp(ppfile):
    Pz = np.load(ppfile)
    v = float(Pz["conv"]) * Pz["g0w"] * Pz["pp"]
    return (np.bincount(r5.orb, weights=v) / r5.osz)[r5.orb]


def op_p3form(v):
    C = np.zeros((n, n))
    sp, sa, sb = s[P], s[A], s[Bm]
    xp, xa, xb = x[P], x[A], x[Bm]
    d = ~rep
    np.add.at(C, (P[d], A[d]), -v[d] / sb[d])
    np.add.at(C, (A[d], P[d]), -v[d] / sb[d])
    np.add.at(C, (P[d], Bm[d]), -v[d] / sa[d])
    np.add.at(C, (Bm[d], P[d]), -v[d] / sa[d])
    np.add.at(C, (A[d], Bm[d]), v[d] / sp[d])
    np.add.at(C, (Bm[d], A[d]), v[d] / sp[d])
    np.add.at(C, (P[d], P[d]), v[d] * np.sinh((xa[d] + xb[d]) / 2) / (sa[d] * sb[d]))
    np.add.at(C, (A[d], A[d]), v[d] * np.sinh((xp[d] - xb[d]) / 2) / (sp[d] * sb[d]))
    np.add.at(C, (Bm[d], Bm[d]), v[d] * np.sinh((xp[d] - xa[d]) / 2) / (sp[d] * sa[d]))
    r = rep
    np.add.at(C, (P[r], A[r]), -v[r] / sa[r])
    np.add.at(C, (A[r], P[r]), -v[r] / sa[r])
    np.add.at(C, (P[r], P[r]), v[r] * np.sinh(xa[r]) / (2 * sa[r] ** 2))
    np.add.at(C, (A[r], A[r]), 2 * v[r] * np.sinh((xp[r] - xa[r]) / 2) / (sp[r] * sa[r]))
    return C


def kt(C):
    Co = Podd @ C @ Podd
    w, V = np.linalg.eigh(Co)
    keep = np.abs(w) > 1e-12 * np.abs(w).max()
    out = []
    for col in (0, 2):
        bc = V[:, keep].T @ r5.b_all[:, col]
        K = float(np.sum(bc ** 2 / w[keep]))
        N = float(np.sum(bc ** 2 / w[keep] ** 2))
        out += [r5.kf * K, N / K / (4 * np.pi)]
    return out


res = {}
v0 = amp(ppref)
# sanity: with the 'cons' convention, rebuild from v reproduces r5 (rates = v / (4 s s s) etc.)
b0 = kt(op_p3form(v0))
res["REF_p3form"] = {"kappa_x": b0[0], "tau_x": b0[1], "kappa_z": b0[2], "tau_z": b0[3]}
print("REF p3form", b0, flush=True)
for spec in specs:
    name, path = spec.split("=", 1)
    b = kt(op_p3form(amp(path)))
    res[name] = {"p3form_kappa_ratio": [b[0] / b0[0], b[2] / b0[2]], "p3form_tau_ratio": [b[1] / b0[1], b[3] / b0[3]]}
    print(name, json.dumps(res[name]), flush=True)
outp.write_text(json.dumps(res, indent=1))
