"""SA1: independent verification of the 7x7x5 full-operator slab two-ports (REF; one range-setting witness).

Operator: own assembly (sa_core.operator) from the raw event list and stored amplitudes; checks: C e = 0, symmetry,
bulk K and delta against the audited values (rs1_validate B5: K = 272.33337519 W/mK, delta_bulk = -1151.5699 ps).
Slab: sa_core independent solver (no parity, no Schur complement, no eigendecomposition), graded z mesh
(hmin 0.02 nm, gamma 0.08, hmax 4 nm), fourth-order (cubic-source) scheme at bisection levels 0 and 1, Richardson
R01 = (16 Y1 - Y0)/15; second-order (linear-source) scheme at level 1 as an independent discretisation; compared with the stored modal two-ports rs2_m775_film2_full_all_REF.npz (black walls).
RTA reduction: the same solver on the energy-conserving RTA operator of the diagonal (K' = 0: inner solve only)
against rs2_m775_film2_rta_all_REF.npz.
usage: python -B sa1_verify775.py [stage ...]   stages: op, rta, full, wit   Output: out/sa1_verify775.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sa_core as sc  # noqa: E402
import numpy as np  # noqa: E402

RS = HERE.parent / "R_slab" / "out"
RES = HERE / "out" / "sa1_verify775.json"
T0 = time.perf_counter()
MESH = {"hmin": 0.02e-9, "gamma": 0.08, "hmax": 4e-9}
FIDX = [0, 20, 40, 60, 70, 79]          # indices into the 80 FDTR frequencies (stored 'all' index = 3 + i)


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


def cplx(Y):
    return [[[float(x.real), float(x.imag)] for x in row] for row in Y]


def save(res):
    RES.write_text(json.dumps(res, indent=1, default=float))


def compare(Yl, Thl, Yref, Tref, order=4):
    nrm, nrt = np.max(np.abs(Yref)), np.max(np.abs(Tref))
    q = 2.0 ** order
    R = (q * Yl[1] - Yl[0]) / (q - 1)
    RT = (q * Thl[1] - Thl[0]) / (q - 1)
    num = np.max(np.abs(Yl[1] - Yl[0]))
    out = {"err_l0": float(np.max(np.abs(Yl[0] - Yref)) / nrm), "err_l1": float(np.max(np.abs(Yl[1] - Yref)) / nrm),
           "err_rich": float(np.max(np.abs(R - Yref)) / nrm),
           "err_rich_Th": float(np.max(np.abs(RT - Tref)) / nrt),
           "err_rich_ImY11_rel_to_ImY11": float(abs((R[0, 0] - Yref[0, 0]).imag) / abs(Yref[0, 0].imag)),
           "err_rich_ImY12_rel_to_ImY12": float(abs((R[0, 1] - Yref[0, 1]).imag) / abs(Yref[0, 1].imag)),
           "change_l0_l1": float(num / nrm),
           "ratio_err_l0_over_l1": float(np.max(np.abs(Yl[0] - Yref)) / np.max(np.abs(Yl[1] - Yref)))}
    return out, R, RT


def main():
    stages = sys.argv[1:] or ["op", "rta", "full", "wit"]
    res = json.loads(RES.read_text()) if RES.exists() else {}
    G = sc.load_modes("m775")
    e, v = G["e"], G["vz"]
    log("modes", G["n"], "orbits", G["n_orb"], "missing images", G["missing_images"])
    C = sc.operator(G, "film2", "REF")
    if "op" in stages:
        kb = sc.kinetic_bulk(C, e, v)
        res["operator"] = {"n_modes": G["n"], "n_orbits": G["n_orb"], "missing_images": G["missing_images"],
                           "energy_residual": float(np.linalg.norm(C @ e) / np.linalg.norm(np.diag(C) * e)),
                           "asym": float(np.abs(C - C.T).max() / np.abs(C).max()),
                           "K_WmK": kb["K"], "delta_bulk_ps": kb["delta"] * 1e12, "tauR_bulk_ps": kb["tauR"] * 1e12,
                           "K_ref_rs1": 272.3333751891106, "delta_ref_rs1_ps": -1151.5698615425028,
                           "rel_K": abs(kb["K"] / 272.3333751891106 - 1),
                           "rel_delta": abs(kb["delta"] * 1e12 / -1151.5698615425028 - 1),
                           "C_tot": G["C_tot"]}
        log("operator", res["operator"])
        save(res)
    Zf = np.load(RS / "rs2_m775_film2_full_all_REF.npz")
    fa, da = Zf["f"], Zf["d"]
    Yf, Tf = Zf["Y"][0], Zf["Th"][0]
    if "rta" in stages:
        Zr = np.load(RS / "rs2_m775_film2_rta_all_REF.npz")
        Yr, Tr = Zr["Y"], Zr["Th"]
        Yr, Tr = (Yr[0], Tr[0]) if Yr.ndim == 5 else (Yr, Tr)
        D = np.diag(C).copy()
        De = D * e
        Z0 = np.zeros_like(C)
        rows = []
        for i in FIDX:
            k = 3 + i
            p = 2j * np.pi * fa[k]
            for di, d in enumerate(da):
                Yl, Thl = [], []
                for lev in (0, 1):
                    z = sc.make_mesh(d, level=lev, **MESH)
                    # the inner solver with D = diag(C) and K' = 0 is exactly the energy-conserving RTA operator CR
                    Y, Th, dg, _ = sc.twoport(None, v, e, p, z, rta_only=True, order=4, Kp=Z0, D=D)
                    Yl.append(Y)
                    Thl.append(Th)
                cmp_, R, RT = compare(Yl, Thl, Yr[di, k], Tr[di, k], order=4)
                row = {"f": float(fa[k]), "d_nm": float(d * 1e9), "N": [int(sc.make_mesh(d, level=l, **MESH).size - 1) for l in (0, 1)],
                       **cmp_, "Y_rich": cplx(R), "Y_modal": cplx(Yr[di, k])}
                rows.append(row)
                log("rta", {kk: vv for kk, vv in row.items() if not kk.startswith("Y_")})
                res["rta_vs_modal_v2"] = rows
                save(res)
    if "full" in stages:
        rows = res.get("full_vs_modal", [])
        done = {(r["f"], r["d_nm"]) for r in rows}
        for i in FIDX:
            k = 3 + i
            p = 2j * np.pi * fa[k]
            for di, d in enumerate(da):
                if (float(fa[k]), float(d * 1e9)) in done:
                    continue
                Yl, Thl, dgs = [], [], []
                x0s = None
                zc = None
                for lev in (0, 1):
                    z = sc.make_mesh(d, level=lev, **MESH)
                    if x0s is not None:
                        x0s = [sc.prolong(f_, zc, z) for f_ in x0s]
                    Y, Th, dg, fs = sc.twoport(C, v, e, p, z, rtol=2e-12, restart=30, maxiter=8, x0s=x0s, order=4)
                    x0s, zc = fs, z
                    Yl.append(Y)
                    Thl.append(Th)
                    dgs.append({"N": dg["N"], "ebal": dg["ebal"], "iters": [g["iters"] for g in dg["gmres"]],
                                "relres": [g["true_relres"] for g in dg["gmres"]], "t": dg["t_total"],
                                "reciprocity": dg["reciprocity"], "symmetry": dg["symmetry"]})
                    del fs
                cmp_, R, RT = compare(Yl, Thl, Yf[di, k], Tf[di, k], order=4)
                # independent discretisation: second-order scheme at level 1 (warm start from the order-4 solution)
                Y2, Th2, dg2, fs2 = sc.twoport(C, v, e, p, z, rtol=2e-12, restart=30, maxiter=8, x0s=x0s, order=2)
                del fs2, x0s
                cmp_["err_order2_l1"] = float(np.max(np.abs(Y2 - Yf[di, k])) / np.max(np.abs(Yf[di, k])))
                cmp_["order2_l1_iters"] = [g["iters"] for g in dg2["gmres"]]
                row = {"f": float(fa[k]), "d_nm": float(d * 1e9), **cmp_, "levels": dgs, "Y_rich": cplx(R),
                       "Th_rich": cplx(RT), "Y_l0": cplx(Yl[0]), "Y_l1": cplx(Yl[1]), "Y_modal": cplx(Yf[di, k])}
                rows.append(row)
                log("full", {kk: vv for kk, vv in row.items() if not kk.startswith(("Y_", "Th_"))})
                res["full_vs_modal"] = rows
                save(res)
    if "wit" in stages:
        nm = "cal3w_m3L100min"
        Cw = sc.operator(G, "film2", nm)
        Zw = np.load(RS / f"rs2_m775_film2_full_sub4_{nm}.npz")
        fw, dw, Yw, Tw = Zw["f"], Zw["d"], Zw["Y"][0], Zw["Th"][0]
        rows = res.get("witness_vs_modal", [])
        for k in (fw.size // 2, fw.size - 1):
            p = 2j * np.pi * fw[k]
            di = 1
            d = dw[di]
            Yl, Thl = [], []
            x0s, zc = None, None
            for lev in (0, 1):
                z = sc.make_mesh(d, level=lev, **MESH)
                if x0s is not None:
                    x0s = [sc.prolong(f_, zc, z) for f_ in x0s]
                Y, Th, dg, fs = sc.twoport(Cw, v, e, p, z, rtol=2e-12, restart=30, maxiter=8, x0s=x0s, order=4)
                x0s, zc = fs, z
                Yl.append(Y)
                Thl.append(Th)
            cmp_, R, RT = compare(Yl, Thl, Yw[di, k], Tw[di, k], order=4)
            row = {"tensor": nm, "f": float(fw[k]), "d_nm": float(d * 1e9), **cmp_, "Y_rich": cplx(R), "Y_modal": cplx(Yw[di, k])}
            rows.append(row)
            log("wit", {kk: vv for kk, vv in row.items() if not kk.startswith("Y_")})
            res["witness_vs_modal"] = rows
            save(res)
    res["seconds"] = time.perf_counter() - T0
    save(res)


if __name__ == "__main__":
    main()
