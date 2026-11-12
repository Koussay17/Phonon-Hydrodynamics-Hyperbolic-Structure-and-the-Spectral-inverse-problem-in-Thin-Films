"""R-A1: film-regularised memory time of the conserving event operator at additional lengths (500 nm, and the
lengths L/c that represent a regulator normalisation c != 1), for REF and the twelve three-mesh optima of E33,
on 7x7x5 and 9x9x6; plus the quantities of objective A-2 (nonlocal time tau_l, normal/umklapp split).

Copy of fc3-memory `scripts_c/film_event.py` (E30; low-memory assembly of E36), changes:
  * length list LS_NM extended (500 nm added; L/c for c in {1/2, 2/3, 3/4, 4/3, 3/2, 2} and L in {30, 100, 300, 500});
    the regulator c |v_z| / L depends on L/c only, so tau(L; c) = tau(L/c; 1) exactly;
  * cross-plane column only (z);
  * the regulator is applied in its energy-conserving form Q D Q (Q = 1 - ehat ehat^T) plus gamma ehat ehat^T,
    gamma = tr C / n. On the odd sector (b and every C^-1 b) Q D Q = D, so tau_z and kappa_z are those of
    film_event.py (C + D); this is checked against the stored E37 values at 30/100/300/1000 nm. The even-sector
    solve needed for tau_l uses the same factor (energy conserved, temperature mode separated);
  * A-2 outputs per L: tau_l (kinetic, ra_kinetic.py), delta = tau_z - tau_l, Rayleigh quotients of the operator
    and of its normal part on u = Q V x (tau_c,RQ, tau_N,RQ), gray normal times from the normal/umklapp
    classification of the events (geom.kind = ev_umk; BZ-representative test of aln_events.py):
      G1  1/tau_N = 4 pi sum_i w_i r^N_i / sum_i w_i,  w_i = b_i^2 / (r_i + d_i)   (kappa_z weights),
      G2  1/tau_N = 4 pi sum_i e_i^2 r^N_i / sum_i e_i^2                          (capacity weights),
    with r^N = diag of the normal-event part of C(g), and the kappa-weighted normal fraction of the scattering rate.
Imports FC3Class from the campaign scripts (read-only). Threads: 2. Memory: one n^2 buffer per L as in E36, the
FastOps arrays not used here are released after construction.
usage: python -B ra1_film_event_ext.py m775|m996 [name ...]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import scipy.linalg as sla  # noqa: E402
import psutil  # noqa: E402

FC3 = Path("D:/ResearchLab/orchestration/campaigns/20261005-fc3-memory")
sys.path.insert(0, str(FC3 / "scripts_c"))
sys.dont_write_bytecode = True
from fc3class import FC3Class, FOURPI  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
FILM3 = Path("C:/Users/Koussay/these/theory/aln/fc3_memory/results_c/film3")

MESH = {"m775": ("results_c/m775x/events_m775_s0.1.npz", "results_c/m775x/film2_pp"),
        "m996": ("results_c/m996/events_m996_s0.1.npz", "results_c/m996/film2_pp")}
WITNESSES = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
            [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]
CS = (0.5, 2 / 3, 0.75, 1.0, 4 / 3, 1.5, 2.0)
LS_BASE = (30, 100, 300, 500)
LS_NM = sorted({round(L / c, 6) for L in LS_BASE for c in CS} | {1000.0})
CHECK_L = ("30", "100", "300", "1000")


def _C_padded(M, g):
    """C(g) on the (n+1)^2 bincount buffer with the dummy slot decoupled (film_event.py, E36)."""
    ops = M.ops
    n1 = ops.n1
    Cf = np.bincount(ops.pair_flat.ravel(), weights=(g[:, None] * ops.pair_val).ravel(), minlength=n1 * n1)
    Cf = Cf.reshape(n1, n1)
    Cf[-1, :] = 0.0
    Cf[:, -1] = 0.0
    Cf[-1, -1] = 1.0
    return Cf


def lkey(L):
    return f"{L:.6g}"


def film_metrics(M, g, kindN, Ls=LS_NM, kind_lo=None, kind_hi=None):
    geom, ops = M.geom, M.ops
    n = geom.n
    vz = geom.vel[:, 2]
    avz = np.abs(vz)
    b = geom.b[:, 2]
    e = geom.e
    ehat = geom.ehat
    e2 = float(e @ e)
    r = ops.diag_full(g)
    rN = ops.diag_full(g * kindN)
    rN_lo = ops.diag_full(g * kind_lo) if kind_lo is not None else None
    rN_hi = ops.diag_full(g * kind_hi) if kind_hi is not None else None
    trC = float(r.sum())
    gam = trC / n
    out = {}
    for L in Ls:
        d = avz / (FOURPI * L * 10.0)              # nm -> Angstrom; THz linewidth units (film_event.py)
        dt = d * ehat                               # D ehat
        c0 = float(ehat @ dt) + gam
        CL = _C_padded(M, g)
        di = np.arange(n)
        CL[di, di] += d
        for i0 in range(0, n, 1024):                # C + Q D Q + gam ehat ehat^T by row blocks (no n^2 temporary)
            i1 = min(i0 + 1024, n)
            CL[i0:i1, :n] += (-ehat[i0:i1, None] * dt[None, :] - dt[i0:i1, None] * ehat[None, :]
                              + c0 * ehat[i0:i1, None] * ehat[None, :])
        cf = sla.cho_factor(CL, lower=True, overwrite_a=True, check_finite=False)
        del CL
        x = sla.cho_solve(cf, np.append(b, 0.0), check_finite=False)[:n]
        K, N = float(b @ x), float(x @ x)
        u = vz * x
        u = u - ehat * float(ehat @ u)
        y = sla.cho_solve(cf, np.append(u, 0.0), check_finite=False)[:n]
        del cf
        Mu = float(u @ y)
        tau = N / K / FOURPI
        taul = e2 * Mu / (FOURPI * K * K)
        # Rayleigh quotients on u (full operator incl. regulator; normal part only)
        pa = ops.ATx(u)
        quN = float(np.sum(g * kindN * pa * pa))
        quT = float(np.sum(g * pa * pa))
        uu = float(u @ u)
        quD = float(np.sum(d * u * u))
        rL = r + d
        wk = b * b / rL
        rec = {"tau_z": tau, "kappa_z": M.kf * K, "tauRTA_z": float(np.sum(b * b / rL ** 2) / np.sum(wk)) / FOURPI,
               "taul_kin": taul, "delta_kin": tau - taul, "L2_A2": Mu / K / FOURPI ** 2,
               "tauc_RQ": uu / (FOURPI * (quT + quD)), "tauN_RQ": uu / (FOURPI * quN) if quN > 0 else None,
               "tauN_G1": float(np.sum(wk)) / (FOURPI * float(np.sum(wk * rN))),
               "tauN_G2": e2 / (FOURPI * float(np.sum(e * e * rN))),
               "fracN_kappa_weighted": float(np.sum(wk * rN) / np.sum(wk * rL)),
               "tauN_G1_orbitN_only": (float(np.sum(wk)) / (FOURPI * float(np.sum(wk * rN_lo)))) if rN_lo is not None else None,
               "tauN_G1_orbitN_any": (float(np.sum(wk)) / (FOURPI * float(np.sum(wk * rN_hi)))) if rN_hi is not None else None,
               "tauN_G2_orbitN_only": (e2 / (FOURPI * float(np.sum(e * e * rN_lo)))) if rN_lo is not None else None,
               "tauN_G2_orbitN_any": (e2 / (FOURPI * float(np.sum(e * e * rN_hi)))) if rN_hi is not None else None,
               "ehat_dot_x_rel": float(abs(ehat @ x) / np.linalg.norm(x)),
               "ehat_dot_y_rel": float(abs(ehat @ y) / np.linalg.norm(y))}
        out[lkey(L)] = rec
    return out


def _memwatch(proc, stop, period=20.0):
    import threading
    def run():
        t0 = time.perf_counter()
        while not stop.wait(period):
            mi = proc.memory_info()
            print(f"  [mem] t={time.perf_counter() - t0:6.0f}s rss={mi.rss / 2**30:.2f} GB "
                  f"peak={getattr(mi, 'peak_wset', mi.rss) / 2**30:.2f} GB", flush=True)
    th = threading.Thread(target=run, daemon=True)
    th.start()
    return th


def main():
    import threading
    mesh = sys.argv[1]
    names = sys.argv[2:] or WITNESSES
    evp, ppd = (FC3 / MESH[mesh][0]), (FC3 / MESH[mesh][1])
    proc = psutil.Process()
    stop = threading.Event()
    _memwatch(proc, stop)
    t0 = time.perf_counter()
    print("setup start", mesh, flush=True)
    M = FC3Class(str(evp), ppd / "pp_r01_REF.npz", allow_missing=True)
    # release arrays that the film metrics do not use (memory guard)
    for attr in ("wpair_flat", "wpair_val", "rp", "Wg_dense"):
        if hasattr(M.ops, attr):
            setattr(M.ops, attr, None)
    for attr in ("A", "W"):
        if hasattr(M.geom, attr):
            setattr(M.geom, attr, None)
    kindN = (M.geom.kind == 0).astype(float)
    geom = M.geom
    vz, b, e, ehat = geom.vel[:, 2], geom.b[:, 2], geom.e, geom.ehat
    QVb = vz * b - ehat * float(ehat @ (vz * b))
    # orbits that mix normal and umklapp events (BZ-boundary representatives): count and REF rate share
    kmin = np.full(M.n_orb, 2); kmax = np.full(M.n_orb, -1)
    np.minimum.at(kmin, M.orb, M.geom.kind); np.maximum.at(kmax, M.orb, M.geom.kind)
    mixed = kmin != kmax
    kind_lo = (kmax[M.orb] == 0).astype(float)          # normal only if every event of the orbit is classified normal
    kind_hi = (kmin[M.orb] == 0).astype(float)          # normal if any event of the orbit is classified normal
    gref = M.rates(np.load(ppd / "pp_r01_REF.npz")["pp"])
    mixed_share = float(gref[mixed[M.orb]].sum() / gref.sum())
    res = {"script": "ra1_film_event_ext.py", "events": str(evp), "pp_dir": str(ppd), "mesh": geom.labels["mesh"],
           "n_modes": int(geom.n), "n_events": int(geom.m), "n_normal": int(kindN.sum()),
           "L_nm": [lkey(L) for L in LS_NM], "c_values": list(CS),
           "n_missing_symmetry_images": M.tie["n_missing_images"],
           "gamma_g": float(e @ e) * float(QVb @ QVb) / float(b @ b) ** 2,
           "n_event_orbits": int(M.n_orb), "n_mixed_NU_orbits": int(mixed.sum()), "REF_rate_share_mixed": mixed_share,
           "REF_rate_share_normal": float(gref[M.geom.kind == 0].sum() / gref.sum()),
           "setup_seconds": time.perf_counter() - t0, "rss_gb_after_setup": proc.memory_info().rss / 2 ** 30}
    print(mesh, "n", geom.n, "m", geom.m, "normal", int(kindN.sum()), "gamma_g", round(res["gamma_g"], 5),
          f"setup {res['setup_seconds']:.0f}s rss {res['rss_gb_after_setup']:.2f} GB", flush=True)
    outp = OUT / f"ra1_film_{mesh}.json"
    for nm in ["REF"] + [x for x in names if x != "REF"]:
        t1 = time.perf_counter()
        g = M.rates(np.load(ppd / f"pp_r01_{nm}.npz")["pp"])
        res[nm] = film_metrics(M, g, kindN, kind_lo=kind_lo, kind_hi=kind_hi)
        res[nm + "_seconds"] = time.perf_counter() - t1
        mi = proc.memory_info()
        res[nm + "_peak_wset_gb"] = getattr(mi, "peak_wset", mi.rss) / 2 ** 30
        k5 = lkey(500.0)
        print(nm, f"tau_z(500) {res[nm][k5]['tau_z']:.4f} taul {res[nm][k5]['taul_kin']:.4f} "
                  f"delta {res[nm][k5]['delta_kin']:.4f} tauN_G1 {res[nm][k5]['tauN_G1']:.3f} "
                  f"{res[nm + '_seconds']:.0f}s peak {res[nm + '_peak_wset_gb']:.2f} GB", flush=True)
        tmp = outp.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(res, indent=1))
        tmp.replace(outp)
    # validation against the stored E37 film metrics (film_event.py, C + D) at 30/100/300/1000 nm
    chk = {}
    for nm in ["REF"] + names:
        src = FILM3 / (f"{names[0] if nm == 'REF' else nm}_{mesh}.json")
        if not src.exists():
            continue
        st = json.loads(src.read_text())
        blk = st["REF"] if nm == "REF" else st[nm]
        for L in CHECK_L:
            for q in ("tau_z", "kappa_z", "tauRTA_z"):
                chk[f"{nm}|{L}|{q}"] = abs(res[nm][lkey(float(L))][q] / blk[L][q] - 1)
    res["check_vs_E37_max_rel"] = max(chk.values()) if chk else None
    res["check_vs_E37"] = chk
    res["seconds"] = time.perf_counter() - t0
    tmp = outp.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(res, indent=1))
    tmp.replace(outp)
    print("check vs E37 max rel", res["check_vs_E37_max_rel"], "total", round(res["seconds"]), "s", flush=True)
    stop.set()


if __name__ == "__main__":
    main()
