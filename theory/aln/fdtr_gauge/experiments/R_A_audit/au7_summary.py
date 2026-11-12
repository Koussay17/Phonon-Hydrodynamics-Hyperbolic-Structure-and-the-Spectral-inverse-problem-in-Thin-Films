"""AU7: combined margins for the resolvability statement of R-A6 (500 and 300 nm, note-16 and favourable stacks).

m = sigma(delta) / (visible half-range of delta), with the visible half-range = (campaign half-range) x f_mesh x T,
  f_mesh : geometric extrapolation of the witness excursions to the next mesh (au4; 1 = campaign reading),
  T      : transfer of a witness-like change of the static delta into the delta returned by a first-order GK analysis
           of the data: T = 1 (campaign model), T_k (exact all-k regularised RTA kernel, au5b), T_s (kinetic slab with
           black walls, au6, finest binning available);
sigma(delta) from au5 (first-order GK film) for N1/N2/N3/N2lam, and for the slab reading also at the slab-fitted base
point (tau_l = tau_R - delta_fit).
Output: out/au7_summary.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import au5_fisher1d as a5  # noqa: E402

AU4 = json.loads((OUT / "au4_delta_range.json").read_text())
AU5 = json.loads((OUT / "au5_fisher1d.json").read_text())
AU5B = json.loads((OUT / "au5b_ktrunc.json").read_text())
nb_files = sorted((p for p in OUT.glob("au6_slab_*.json") if p.stem.split("_")[-1].isdigit()),
                  key=lambda p: int(p.stem.split("_")[-1]))
AU6 = {int(p.stem.split("_")[-1]): json.loads(p.read_text()) for p in nb_files}
NBMAX = max(k for k, v in AU6.items() if "500" in v and "transfer" in json.dumps(v["500"]))
red = {p.stem: json.loads(p.read_text()) for p in OUT.glob("au6_slab_*_eta0.json")}
res = {"slab_nbins_used": NBMAX, "slab_convergence": {}, "margins": {}}
for nb, d in list(AU6.items()) + [(k.replace("au6_slab_", ""), v) for k, v in red.items()]:
    for L in ("300", "500"):
        if L not in d:
            continue
        for port in ("wall", "energy"):
            r = d[L].get(f"{port}|note16", {})
            for N in ("N2lam|0.01", "N3|0.01"):
                if N in r:
                    res["slab_convergence"][f"{nb}|{L}|{port}|{N}"] = {
                        "delta_fit_ps": r[N]["delta_fit_ps"], "delta_M_ps": r[N]["delta_M_ps"],
                        "ratio": r[N]["delta_fit_ps"] / r[N]["delta_M_ps"],
                        "transfer": [r[N]["transfer"][x]["transfer"] for x in r[N]["transfer"]],
                        "dchi2_gk_taulM": r[N]["dchi2_gk_taulM"], "dchi2_fourier": r[N]["dchi2_fourier"]}
for L in (300, 500):
    Lk = f"{float(L):.6g}"
    t4 = AU4["table"][Lk]["all"]
    hr = t4["base_point"]["halfrange_delta_ps_campaign_rule"]
    fmesh = t4["x12_geometric"]["w_delta"] / t4["robust"]["w_delta"]
    b = a5.A1[str(L)]
    for Rtag, (R1, R2) in {"note16": (1e-8, 2e-8), "favourable": (2e-9, 2e-9)}.items():
        kb = AU5B[str(L)][Rtag]
        sl = AU6[NBMAX][str(L)][f"wall|{Rtag}"]
        out = {}
        for N in ("N1", "N2", "N3", "N2lam"):
            for sp in ("0.01", "0.1"):
                key = f"{N}|{sp}"
                sig = AU5["fdtr"][f"{L}|{Rtag}|kinetic"][key]["sigma_delta_ps_via_tauR"]
                m0 = sig / hr
                Tk = None
                if key in kb:
                    Tk = float(np.mean([kb[key]["transfer"][x]["transfer_exact_over_static"] for x in ("-0.1", "0.1")]))
                Ts, sig_s, dfit = None, None, None
                if key in sl:
                    Ts = float(np.mean([sl[key]["transfer"][x]["transfer"] for x in ("-0.1", "0.1")]))
                    dfit = sl[key]["delta_fit_ps"]
                    # sigma(delta) at the slab-fitted base point
                    base = {"tauR": b["tau_R_ps"] * 1e-12, "taul": (b["tau_R_ps"] - dfit) * 1e-12, "lamz": b["lam_z"],
                            "C": a5.C_ALN, "R1": R1, "R2": R2, "d": L * 1e-9}
                    J = a5.jac(base)
                    nuis = {**a5.SETS, **a5.SETS_EXTRA}[N]
                    sig_s = a5.sigma_dir(J, (1, 0), nuis, float(sp))
                out[key] = {"sigma_delta_ps": sig, "halfrange_campaign_ps": hr, "m_campaign": m0,
                            "f_mesh": fmesh, "m_mesh": m0 / fmesh,
                            "T_ktrunc": Tk, "m_mesh_ktrunc": (m0 / fmesh / Tk) if Tk else None,
                            "T_slab": Ts, "delta_fit_slab_ps": dfit, "sigma_delta_at_slab_base_ps": sig_s,
                            "m_slab": (sig_s / (hr * Ts)) if (Ts and sig_s) else None,
                            "m_slab_mesh": (sig_s / (hr * Ts * fmesh)) if (Ts and sig_s) else None}
        res["margins"][f"{L}|{Rtag}"] = out
        print(L, Rtag, {k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()
                            if kk in ("m_campaign", "m_mesh", "m_mesh_ktrunc", "m_slab", "m_slab_mesh", "T_ktrunc", "T_slab")}
                        for k, v in out.items()}, flush=True)
(OUT / "au7_summary.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res["slab_convergence"], indent=0))
