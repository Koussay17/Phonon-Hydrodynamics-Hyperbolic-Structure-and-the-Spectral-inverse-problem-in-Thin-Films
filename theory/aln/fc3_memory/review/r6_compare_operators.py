"""Review R6 (part 2, research env): entry-by-entry comparison of phono3py's symmetrised reducible collision matrix
Omega with the event operator C (independent implementation r5, conserving convention, orbit-averaged rates) on
the odd sector (the sector of the heat current), for the reference and witnesses.
usage: python -B r6_compare_operators.py <events.npz> <pp_ref.npz> <omega_dir> <outjson> name=pp.npz ...
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.argv_saved = list(sys.argv)
evp, ppref, odir, outp = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
specs = sys.argv[5:]
# reuse the independent event operator of r5 (module-level code needs argv; build a small shim)
sys.argv = ["r5", evp, ppref, str(outp.with_name("r6_r5_dummy.json"))]
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location("r5", str(Path(__file__).with_name("r5_event_operator.py")))
r5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r5)
sys.argv = sys.argv_saved

Z = np.load(evp)
Nq, nb = Z["freqs"].shape
valid = r5.valid
Jm = Z["Jmap"]
Jfull = (Jm[:, None] * nb + np.arange(nb)[None, :]).ravel()
newidx = r5.idx
Jv = newidx[Jfull][valid]             # J on valid modes
assert (Jv >= 0).all()
n = r5.n
Podd = 0.5 * (np.eye(n) - np.eye(n)[Jv])


def load_omega(name):
    O = np.load(odir / f"omega_{name}.npz")
    Om = O["Omega"][np.ix_(valid, valid)]
    X = O["X"][valid]
    return Om, X, O["kappa"]


res = {}
Om0, X0, k0 = load_omega("REF")
C0 = r5.operator(r5.rates(ppref), "cons")
b = r5.b_all
sx = float(X0[:, 0] @ b[:, 0] / (b[:, 0] @ b[:, 0]))
res["X_over_b_scale"] = sx
res["X_prop_b_resid"] = float(np.linalg.norm(X0 - sx * b) / np.linalg.norm(X0))
A0, B0 = Podd @ Om0 @ Podd, Podd @ C0 @ Podd
s = float(np.sum(A0 * B0) / np.sum(B0 * B0))
res["Omega_over_C_scale_odd"] = s
res["REF_odd_rel_diff"] = float(np.linalg.norm(A0 - s * B0) / np.linalg.norm(A0))
x0 = np.linalg.pinv(B0, rcond=1e-12) @ b[:, 0]
res["REF_resid_Omega_on_Cplus_b_x"] = float(np.linalg.norm(Podd @ (Om0 @ x0) - s * b[:, 0]) / np.linalg.norm(s * b[:, 0]))
for spec_ in specs:
    name, path = spec_.split("=", 1)
    Om, X, kk = load_omega(name)
    C = r5.operator(r5.rates(path), "cons")
    A, Bw = Podd @ Om @ Podd, Podd @ C @ Podd
    dA, dB = A - A0, s * (Bw - B0)
    rec = {"odd_rel_diff": float(np.linalg.norm(A - s * Bw) / np.linalg.norm(A)),
           "delta_rel_diff": float(np.linalg.norm(dA - dB) / np.linalg.norm(dA)),
           "delta_norm_ratio": float(np.linalg.norm(dB) / np.linalg.norm(dA))}
    # where do the changes differ? row-wise contribution
    E = dA - dB
    rown = np.sqrt(np.sum(E ** 2, axis=1))
    top = np.argsort(-rown)[:12]
    fval = r5.eps
    rec["top_rows"] = [{"mode": int(i), "q": int(np.nonzero(valid)[0][i] // nb), "band": int(np.nonzero(valid)[0][i] % nb),
                        "freq": float(fval[i]), "row_err": float(rown[i]),
                        "row_err_rel_to_rowdelta": float(rown[i] / max(np.sqrt(np.sum(dA[i] ** 2)), 1e-300)),
                        "diag_A": float(A[i, i]), "diag_sB": float(s * Bw[i, i]), "diag_A0": float(A0[i, i])} for i in top]
    rec["frac_err_in_top12_rows"] = float(np.sum(rown[top] ** 2) / np.sum(rown ** 2))
    # diagonal versus off-diagonal share of the discrepancy
    dg = np.diag(E)
    rec["diag_share_of_err2"] = float(np.sum(dg ** 2) / np.sum(E ** 2))
    # kappa from each operator restricted to odd sector
    for col, nm in ((0, "x"), (2, "z")):
        xA = np.linalg.pinv(A, rcond=1e-12) @ b[:, col]
        xB = np.linalg.pinv(s * Bw, rcond=1e-12) @ b[:, col]
        xA0 = np.linalg.pinv(A0, rcond=1e-12) @ b[:, col]
        xB0 = np.linalg.pinv(s * B0, rcond=1e-12) @ b[:, col]
        rec[f"kappa_ratio_{nm}_Omega"] = float((b[:, col] @ xA) / (b[:, col] @ xA0))
        rec[f"kappa_ratio_{nm}_C"] = float((b[:, col] @ xB) / (b[:, col] @ xB0))
    res[name] = rec
    print(name, json.dumps({k: v for k, v in rec.items() if k != "top_rows"}), flush=True)
outp.write_text(json.dumps(res, indent=1))
print(json.dumps({k: v for k, v in res.items() if not isinstance(v, dict)}, indent=1))
