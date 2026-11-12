"""RS1: validation of the full-operator slab solver (rs_core.py).

Stages (argv[1], comma list; default all):
  A  synthetic: (A1) one pair, grey two-stream with energy conservation, against the closed form of au6_test;
     (A2) random parity-symmetric conserving operator (not RTA: dense even and odd sectors, local modes), paired
     solver vs unpaired first-order solver vs adding-doubling vs Chebyshev collocation of the full first-order system
     (N = 160, 240), at 1e4 ... 3e9 Hz.
  B1 RTA reduction: Chat -> energy-conserving RTA on diag(Chat) (REF, 7x7x5, unbinned), paired solver vs the audit's
     RTA slab solver au6_slab.slab_twoport with one class per mode pair (independent code), 300/500 nm.
  B2 full operator REF 7x7x5: paired vs unpaired (no parity, eigenvalues lam instead of lam^2).
  B3 full operator REF 7x7x5: paired vs adding-doubling (no eigendecomposition); Taylor order and number of doublings.
  B4 real p: eigenvalues lam of the paired route (geev on M, unpolished) vs the well-conditioned SVD route.
  B5 bulk limits: (a) diffusive eigenvalue mu_d(p) = (pC/K)(1 + delta p + ...) of the bulk operator against the KKT
     values K, tau_R - tau_l; (b) the same for the regularised operator C + QDQ against au3 (E-QDQ); (c) static
     conductance G(d) of the slab for d = 1 ... 100 um: 1/G = d/K + 2 R_w, slope against K_bulk.
Output: out/rs1_validate.json (+ log)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rs_core as rc  # noqa: E402

AUD = HERE.parent / "R_A_audit"
RES = OUT / "rs1_validate.json"
T0 = time.perf_counter()


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


def save(res):
    RES.write_text(json.dumps(res, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x)))


def relmax(Y1, Y2):
    return float(np.max(np.abs(np.asarray(Y1) - np.asarray(Y2))) / np.max(np.abs(Y2)))


# ------------------------------------------------------------------------------------------------ synthetic
def random_parity_operator(n_p, n_0, seed=1, rate_lo=1e10, rate_hi=3e12, v_lo=300.0, v_hi=8000.0):
    rng = np.random.default_rng(seed)
    ep = rng.uniform(0.5, 1.5, n_p) * 1e3
    e0 = rng.uniform(0.5, 1.5, n_0) * 1e3
    v = np.exp(rng.uniform(np.log(v_lo), np.log(v_hi), n_p))
    # even sector (n_p + n_0) with kernel e_E = (sqrt2 ep, e0); odd sector positive definite
    nE = n_p + n_0
    eE = np.concatenate([np.sqrt(2) * ep, e0])
    eEh = eE / np.linalg.norm(eE)
    rE = np.exp(rng.uniform(np.log(rate_lo), np.log(rate_hi), nE))
    AE = rng.standard_normal((nE, nE)) * np.sqrt(rE)[:, None]
    CE = np.diag(rE) + 0.5 * (AE @ AE.T) / nE
    Q = np.eye(nE) - np.outer(eEh, eEh)
    CE = Q @ CE @ Q
    rO = np.exp(rng.uniform(np.log(rate_lo), np.log(rate_hi), n_p))
    AO = rng.standard_normal((n_p, n_p)) * np.sqrt(rO)[:, None]
    CO = np.diag(rO) + 0.5 * (AO @ AO.T) / n_p
    Ce, Cp0, C00 = CE[:n_p, :n_p], CE[:n_p, n_p:] / np.sqrt(2), CE[n_p:, n_p:]
    Cpp, Cpm = 0.5 * (Ce + CO), 0.5 * (Ce - CO)
    Bk = rc.blocks_from_parts(ep, e0, v, Cpp, Cpm, Cp0, C00)
    # full matrix in (+, -, 0) ordering for the unpaired and collocation routes
    C = np.block([[Cpp, Cpm, Cp0], [Cpm, Cpp, Cp0], [Cp0.T, Cp0.T, C00]])
    vz = np.concatenate([v, -v, np.zeros(n_0)])
    e = np.concatenate([ep, ep, e0])
    G = {"vz": vz, "e": e, "C_tot": float(e @ e)}
    return Bk, C, G


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def colloc_full(C, G, d, p, N):
    """Chebyshev collocation of V f' + (p + C) f = 0 on [0, d] for all modes; black walls. Independent of rs_core."""
    vz, e = G["vz"], G["e"]
    nm = vz.size
    D, x = cheb(N)
    Dz = -2.0 / d * D                     # z = (1 - x) d / 2
    n1 = N + 1
    A = np.kron(np.diag(vz).astype(complex), Dz) + np.kron(C + p * np.eye(nm), np.eye(n1))
    out = []
    for TA, TB in ((1.0, 0.0), (0.0, 1.0)):
        Ab = A.copy()
        rb = np.zeros(nm * n1, complex)
        for i in range(nm):
            if vz[i] > 0:          # inflow at z = 0 (x = 1, node 0)
                r = i * n1 + 0
                Ab[r, :] = 0
                Ab[r, r] = 1
                rb[r] = TA * e[i]
            elif vz[i] < 0:        # inflow at z = d (x = -1, node N)
                r = i * n1 + N
                Ab[r, :] = 0
                Ab[r, r] = 1
                rb[r] = TB * e[i]
        sol = np.linalg.solve(Ab, rb).reshape(nm, n1)
        q0 = np.sum(vz * e * sol[:, 0])
        qd = np.sum(vz * e * sol[:, N])
        th0 = np.sum(e * sol[:, 0]) / G["C_tot"]
        thd = np.sum(e * sol[:, N]) / G["C_tot"]
        out.append((q0, qd, th0, thd))
    Y = np.array([[out[0][0], out[1][0]], [out[0][1], out[1][1]]])
    Th = np.array([[out[0][2], out[1][2]], [out[0][3], out[1][3]]])
    return Y, Th


def stage_A(res):
    r = {}
    # A1 two-stream
    c, v, tau, d = 1.3e6, 4000.0, 80e-12, 400e-9
    ep = np.array([np.sqrt(c)])
    Cpp = np.array([[0.5 / tau]])
    Cpm = np.array([[-0.5 / tau]])
    Bk = rc.blocks_from_parts(ep, np.zeros(0), np.array([v]), Cpp, Cpm, np.zeros((1, 0)), np.zeros((0, 0)))
    rows = []
    for f in (1e3, 1e6, 1e8, 1e9, 1e10):
        p = 2j * np.pi * f
        (Y, Th, dg), = rc.slab_paired(Bk, p, (d,))[0]
        g = p + 1 / tau
        kr = np.sqrt(p / g)
        lam = np.sqrt(g * p) / v
        ch, sh = np.cosh(lam * d), np.sinh(lam * d)
        q0 = 2 * c * v * kr * (ch + kr * sh) / ((1 + kr ** 2) * sh + 2 * kr * ch)
        Yd, Thd, _ = rc.slab_doubling(Bk, p, d, K=30, order=10)
        rows.append({"f": f, "paired_rel_err_Y11": float(abs(Y[0, 0] / q0 - 1)),
                     "doubling_rel_err_Y11": float(abs(Yd[0, 0] / q0 - 1)), "energy_balance": dg["energy_balance_rel"]})
    r["A1_two_stream"] = rows
    log("A1", rows)
    # A2 random parity-symmetric conserving operator
    rows = []
    Bk, C, G = random_parity_operator(6, 3, seed=7, rate_lo=1e10, rate_hi=3e11, v_lo=1000.0, v_hi=8000.0)
    d = 400e-9
    for f in (1e4, 1e7, 2e8, 3e9):
        p = 2j * np.pi * f
        (Yp, Thp, dg), = rc.slab_paired(Bk, p, (d,))[0]
        (Yu, Thu), = rc.slab_unpaired(G, C, p, (d,))
        Yd, Thd, di = rc.slab_doubling(Bk, p, d, K=None, order=10)
        Yc, Thc = colloc_full(C, G, d, p, 160)
        Yc2, Thc2 = colloc_full(C, G, d, p, 240)
        rows.append({"f": f, "paired_vs_unpaired_Y": relmax(Yp, Yu), "paired_vs_unpaired_Th": relmax(Thp, Thu),
                     "paired_vs_doubling_Y": relmax(Yp, Yd), "paired_vs_doubling_Th": relmax(Thp, Thd),
                     "paired_vs_colloc240_Y": relmax(Yp, Yc2), "paired_vs_colloc240_Th": relmax(Thp, Thc2),
                     "colloc160_vs_240_Y": relmax(Yc, Yc2), "energy_balance": dg["energy_balance_rel"],
                     "doubling_norm_hG": di["norm_hG"], "doubling_parity_R": di["R_top_vs_bottom"]})
        log("A2", rows[-1])
    r["A2_random_operator_6pairs_3local"] = rows
    # A2b stiffer operator (rates up to 1e13, v down to 30 m/s): paired vs unpaired vs doubling only
    rows = []
    Bk, C, G = random_parity_operator(40, 12, seed=11, rate_lo=1e9, rate_hi=1e13, v_lo=30.0, v_hi=1e4)
    for f in (1e4, 1e6, 2e8):
        p = 2j * np.pi * f
        (Yp, Thp, dg), = rc.slab_paired(Bk, p, (d,))[0]
        (Yu, Thu), = rc.slab_unpaired(G, C, p, (d,))
        Yd, Thd, di = rc.slab_doubling(Bk, p, d, K=None, order=10)
        rows.append({"f": f, "paired_vs_unpaired_Y": relmax(Yp, Yu), "paired_vs_doubling_Y": relmax(Yp, Yd),
                     "paired_vs_unpaired_Th": relmax(Thp, Thu), "paired_vs_doubling_Th": relmax(Thp, Thd),
                     "energy_balance": dg["energy_balance_rel"], "doubling_norm_hG": di["norm_hG"], "doubling_K": di["K"]})
        log("A2b", rows[-1])
    r["A2b_stiff_random_operator_40pairs_12local"] = rows
    res["A"] = r


def stage_A3(res):
    """Partially transmitting walls: paired (modal BC) vs doubling (walls as a separate linear system)."""
    rows = []
    for tag, args in (("6p3l", dict(n_p=6, n_0=3, seed=7, rate_lo=1e10, rate_hi=3e11, v_lo=1000.0, v_hi=8000.0)),
                      ("40p12l_stiff", dict(n_p=40, n_0=12, seed=11, rate_lo=1e9, rate_hi=1e13, v_lo=30.0, v_hi=1e4))):
        Bk, C, G = random_parity_operator(**args)
        d = 400e-9
        walls = [("black",), ("specular", 0.5), ("specular", 0.9), ("diffuse", 0.5), ("diffuse", 0.1)]
        for f in (1e4, 2e8):
            p = 2j * np.pi * f
            outs, info = rc.slab_paired(Bk, p, (d,), walls=walls)
            for w, o in zip(walls, outs):
                Y, Th, dg = o[0]
                Yd, Thd, di = rc.slab_doubling(Bk, p, d, K=None, order=10, wall=w)
                rows.append({"op": tag, "f": f, "wall": list(w), "rel_Y": relmax(Y, Yd), "rel_Th": relmax(Th, Thd),
                             "energy_balance": dg["energy_balance_rel"], "bc_residual": dg["bc_residual_rel"],
                             "G": float(Y[0, 0].real)})
                log("A3", rows[-1])
    res["A3_partial_walls"] = rows


# ------------------------------------------------------------------------------------------------ AlN REF
_cache = {}


def ref775():
    if "G" not in _cache:
        G = rc.geometry("m775")
        C = rc.operator(G, "film2", "REF")
        _cache.update({"G": G, "C": C, "Bk": rc.blocks(G, C)})
    return _cache["G"], _cache["C"], _cache["Bk"]


def stage_B1(res):
    sys.path.insert(0, str(AUD))
    import au6_slab as a6
    G, C, _ = ref775()
    Cr = rc.rta_operator(G, C)
    Bk = rc.blocks(G, Cr)
    r = np.diag(C)
    ip, i0 = G["ip"], G["i0"]
    cl = {"c": G["cap"][ip], "tau": 1.0 / r[ip], "v": G["vz"][ip], "c0": G["cap"][i0], "tau0": 1.0 / r[i0]}
    rows = []
    for f in (1e2, 1e4, 1e6, 3e7, 2e8):
        p = 2j * np.pi * f
        outs, info = rc.slab_paired(Bk, p, (300e-9, 500e-9))
        for (Y, Th, dg), d in zip(outs, (300e-9, 500e-9)):
            Ya, Tha, dga = a6.slab_twoport(cl, d, p)
            rows.append({"f": f, "d_nm": d * 1e9, "rel_Y": relmax(Y, Ya), "rel_Th": relmax(Th, Tha),
                         "energy_balance_full": dg["energy_balance_rel"], "energy_balance_au6": dga["energy_balance_rel"],
                         "Y11": [Y[0, 0].real, Y[0, 0].imag]})
            log("B1", rows[-1])
    res["B1_rta_vs_au6"] = rows


def stage_B2(res):
    G, C, Bk = ref775()
    rows = []
    for f in (1e4, 1e6, 5e7, 2e8):
        p = 2j * np.pi * f
        outs, info = rc.slab_paired(Bk, p, (300e-9, 500e-9))
        outu = rc.slab_unpaired(G, C, p, (300e-9, 500e-9))
        for (Y, Th, dg), (Yu, Thu), d in zip(outs, outu, (300, 500)):
            rows.append({"f": f, "d_nm": d, "rel_Y": relmax(Y, Yu), "rel_Th": relmax(Th, Thu),
                         "rel_Y11_minus_static_part": float(abs((Y[0, 0] - Yu[0, 0]) / (Y[0, 0].imag))),
                         "unpaired_reciprocity": float(abs(Yu[0, 1] + Yu[1, 0]) / abs(Yu[0, 1])),
                         "unpaired_symmetry": float(abs(Yu[0, 0] + Yu[1, 1]) / abs(Yu[0, 0])),
                         "eig_backward_max": info["eig_backward_max"], "polish_shift_rel": info["polish_shift_rel"]})
            log("B2", rows[-1])
        save(res)
    res["B2_paired_vs_unpaired"] = rows


def stage_B3(res):
    G, C, Bk = ref775()
    rows = []
    cases = [(1e6, 500e-9, 23, 10), (1e6, 500e-9, 23, 6), (1e6, 500e-9, 23, 4), (1e6, 500e-9, 20, 10),
             (1e6, 500e-9, 26, 10), (2e8, 500e-9, 23, 10), (1e4, 300e-9, 23, 10), (5e7, 300e-9, 23, 10)]
    ref = {}
    for f, d, K, order in cases:
        p = 2j * np.pi * f
        if (f, d) not in ref:
            ref[(f, d)] = rc.slab_paired(Bk, p, (d,))[0][0]
        Y, Th, dg = ref[(f, d)]
        Yd, Thd, di = rc.slab_doubling(Bk, p, d, K=K, order=order)
        rows.append({"f": f, "d_nm": d * 1e9, "K": K, "order": order, "rel_Y": relmax(Yd, Y), "rel_Th": relmax(Thd, Th),
                     "rel_ImY11": float(abs((Yd[0, 0] - Y[0, 0]).imag / Y[0, 0].imag)), **di})
        log("B3", rows[-1])
        res["B3_paired_vs_doubling"] = rows
        save(res)


def stage_B4(res):
    G, C, Bk = ref775()
    rows = []
    for p in (2 * np.pi * 1e4, 2 * np.pi * 1e6, 2 * np.pi * 2e8):
        ls = rc.slab_eigs_real_svd(Bk, p)
        lg = rc.slab_eigs_real_paired(Bk, p)
        rel = np.abs(lg - ls) / ls
        rows.append({"p": p, "lam_min_svd": float(ls[0]), "lam_min_geev": float(lg[0]), "rel_err_lam_min": float(rel[0]),
                     "rel_err_lam_2nd": float(rel[1]), "rel_err_max": float(rel.max()),
                     "rel_err_median": float(np.median(rel)), "lam_max": float(ls[-1]),
                     "spread_lam": float(ls[-1] / ls[0])})
        log("B4", rows[-1])
    res["B4_real_p_svd_vs_geev"] = rows


def stage_B5(res):
    G, C, Bk = ref775()
    au3 = json.loads((AUD / "out" / "au3_event775.json").read_text())["tensors"]["REF"]
    r = {}

    def mu_diff(Bk_, p):
        lam = rc.slab_eigs_real_svd(Bk_, p)
        return lam[0] ** 2

    for tag, Bk_, kt in (("bulk", Bk, rc.kinetic_times(G, C)),
                         ("reg300", rc.blocks(G, rc.regulated(G, C, 300e-9)), rc.kinetic_times(G, rc.regulated(G, C, 300e-9))),
                         ("reg500", rc.blocks(G, rc.regulated(G, C, 500e-9)), rc.kinetic_times(G, rc.regulated(G, C, 500e-9)))):
        ps = np.array([2.5e5, 5e5, 1e6, 2e6]) * (1.0 if tag == "bulk" else 10.0)
        y = np.array([mu_diff(Bk_, p) / p for p in ps])      # = (C/K)(1 + delta p + c2 p^2 + ...)
        # polynomial fit in p (cubic, exact through 4 points)
        cf = np.polyfit(ps, y, 3)
        a0, a1 = cf[-1], cf[-2]
        K_fit = G["C_tot"] / a0
        delta_fit = a1 / a0
        ent = {"K_eig": K_fit, "K_kkt": kt["K"], "rel_K": abs(K_fit / kt["K"] - 1),
               "delta_eig_ps": delta_fit * 1e12, "delta_kkt_ps": kt["delta"] * 1e12,
               "rel_delta": abs(delta_fit / kt["delta"] - 1), "p_values": ps.tolist()}
        if tag != "bulk":
            L = tag[3:]
            ent["delta_au3_ps"] = au3[L]["E-QDQ"]["tau_ps"] - au3[L]["E-QDQ"]["taul_ps"]
        r[tag] = ent
        log("B5", tag, ent)
    # static conductance vs thickness (bulk operator, no regulator)
    rows = []
    for d in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
        (Y, Th, dg), = rc.slab_paired(Bk, 2j * np.pi * 1.0, (d,))[0]
        rows.append({"d": d, "G": float(Y[0, 0].real), "invG": float(1 / Y[0, 0].real), "Gd_over_K": float(Y[0, 0].real * d / r["bulk"]["K_kkt"]),
                     "energy_balance": dg["energy_balance_rel"]})
        log("B5c", rows[-1])
    ds = np.array([x["d"] for x in rows])
    iG = np.array([x["invG"] for x in rows])
    sl, ic = np.polyfit(ds[-3:], iG[-3:], 1)
    r["large_d"] = {"rows": rows, "K_from_slope": 1 / sl, "rel_K_slope": abs(1 / sl / r["bulk"]["K_kkt"] - 1),
                    "2Rw_intercept": ic}
    log("B5c fit", r["large_d"]["K_from_slope"], r["large_d"]["rel_K_slope"], ic)
    res["B5_bulk_limits"] = r


def main():
    stages = sys.argv[1].split(",") if len(sys.argv) > 1 else ["A", "B1", "B2", "B3", "B4", "B5"]
    res = json.loads(RES.read_text()) if RES.exists() else {}
    for st in stages:
        {"A": stage_A, "A3": stage_A3, "B1": stage_B1, "B2": stage_B2, "B3": stage_B3, "B4": stage_B4, "B5": stage_B5}[st](res)
        res["seconds_" + st] = time.perf_counter() - T0
        save(res)


if __name__ == "__main__":
    main()
