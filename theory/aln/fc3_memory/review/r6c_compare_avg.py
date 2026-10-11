"""Review R6c: as R6 but with the degeneracy-averaged event operator (phono3py-style averaging, kappa-neutral),
so that the gauge split of degenerate partners does not dominate; locates the kappa-relevant discrepancy
between phono3py's Omega and the event operator, and splits it into diagonal / off-diagonal parts and by
the detuning of the events touching each mode.
usage: python -B r6c_compare_avg.py <events.npz> <pp_ref.npz> <omega_dir> <outjson> name=pp.npz ...
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np
evp, ppref, odir, outp = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
specs = sys.argv[5:]
saved = list(sys.argv)
sys.argv = ["r6b", evp, ppref, str(outp.with_name("r6c_dummy.json"))]
sp_ = importlib.util.spec_from_file_location("r6b", str(Path(__file__).with_name("r6b_degeneracy_convention.py")))
r6b = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(r6b)
sys.argv = saved
r5 = r6b.r5
valid, n = r5.valid, r5.n
Jm = r5.Z["Jmap"]; nb = r5.Z["freqs"].shape[1]
Jv = r5.idx[(Jm[:, None] * nb + np.arange(nb)[None, :]).ravel()][valid]
Podd = 0.5 * (np.eye(n) - np.eye(n)[Jv])
def om(name):
    O = np.load(odir / f"omega_{name}.npz"); return O["Omega"][np.ix_(valid, valid)]
A0 = Podd @ om("REF") @ Podd
C0 = r6b.p3style(r5.operator(r5.rates(ppref), "cons"))
B0 = Podd @ C0 @ Podd
s = float(np.sum(A0 * B0) / np.sum(B0 * B0))
res = {"scale": s, "REF_odd_rel_diff": float(np.linalg.norm(A0 - s * B0) / np.linalg.norm(A0))}
E0 = A0 - s * B0
res["REF_diag_share"] = float(np.sum(np.diag(E0) ** 2) / np.sum(E0 ** 2))
# per-mode max detuning share: fraction of diag C from events with |Delta| > 2 sigma
dl = np.abs(r5.Z["ev_dl"]); sig = float(r5.Z["sigma"])
for spec in specs:
    name, path = spec.split("=", 1)
    g = r5.rates(path)
    C = r6b.p3style(r5.operator(g, "cons"))
    A, B = Podd @ om(name) @ Podd, Podd @ C @ Podd
    dA, dB = A - A0, s * (B - B0)
    E = dA - dB
    rown = np.sqrt(np.sum(E ** 2, axis=1))
    top = np.argsort(-rown)[:10]
    rec = {"delta_rel_diff": float(np.linalg.norm(E) / np.linalg.norm(dA)),
           "diag_share": float(np.sum(np.diag(E) ** 2) / np.sum(E ** 2)),
           "frac_top10_rows": float(np.sum(rown[top] ** 2) / np.sum(rown ** 2)),
           "top": [{"mode": int(i), "freq": float(r5.eps[i]), "rowerr_rel": float(rown[i] / max(np.linalg.norm(dA[i]), 1e-300)),
                    "dA_diag": float(dA[i, i]), "dB_diag": float(dB[i, i]), "A0_diag": float(A0[i, i])} for i in top]}
    # kappa sensitivity of the discrepancy: first-order change of K_x from E: dK = -x^T E x, x = A^-1 b
    for col, nm in ((0, "x"), (2, "z")):
        x = np.linalg.pinv(A, rcond=1e-12) @ r5.b_all[:, col]
        K = r5.b_all[:, col] @ x
        Ed = np.diag(np.diag(E)); Eo = E - Ed
        rec[f"dK_rel_{nm}_from_diag"] = float(x @ Ed @ x / K)
        rec[f"dK_rel_{nm}_from_offdiag"] = float(x @ Eo @ x / K)
    res[name] = rec
    print(name, json.dumps({k: v for k, v in rec.items() if k != "top"}))
    for t in rec["top"][:6]: print("   ", t)
outp.write_text(json.dumps(res, indent=1))
print({k: v for k, v in res.items() if not isinstance(v, dict)})
