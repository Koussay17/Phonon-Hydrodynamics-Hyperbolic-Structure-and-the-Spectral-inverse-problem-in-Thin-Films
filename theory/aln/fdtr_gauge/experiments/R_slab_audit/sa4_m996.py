"""SA4: full-operator slab on the 9x9x6 event mesh with the independent solver (no parity reduction needed: on 9x9x6
the zone-boundary plane q_z = 1/2 is mapped onto itself by q_z -> -q_z and carries moving modes, which the paired
reduction of rs_core cannot treat; the present solver treats every mode with its own v_z).

Stages:
  op     REF operator (own assembly), C e = 0, bulk K, delta; regularised (C + QDQ, |v_z|/L) kinetic times at 300/500 nm.
  slab   REF two-ports at the support frequencies (1 kHz + every 6th FDTR frequency + 200 MHz), 500 nm (and 300 nm if
         'slab300' is given), fourth-order scheme, level-0 graded mesh; a coarser second mesh at 10 kHz and 200 MHz as
         a discretisation-error estimate. Incremental npz output.
  fit    AAA interpolation to the 80 FDTR frequencies, GK fits (N2lam, N3; 0.01/0.1 deg; wall port, note-16), Fisher
         sigma at the fitted base point (sa_fdtr).
  proxy  c calibrated on REF so that delta_reg(L/c) = delta_fit (regularised C + QDQ), then delta_reg(L/c) of the 12
         campaign witnesses on 9x9x6 (predicted visible excursions).
  wit    slab of one witness (argv: name) at 500 nm, same frequencies; fit with the REF base.
usage: python -B sa4_m996.py <stage> [...]     Output: out/sa4_m996.json, out/sa4_m996_<tensor>_<d>.npz
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sa_core as sc  # noqa: E402
import sa_fdtr as fd  # noqa: E402
import numpy as np  # noqa: E402
from scipy.interpolate import AAA  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

RES = HERE / "out" / "sa4_m996.json"
T0 = time.perf_counter()
MESH = {"hmin": 0.02e-9, "gamma": 0.08, "hmax": 4e-9}
MESH_COARSE = {"hmin": 0.02e-9, "gamma": 0.12, "hmax": 6e-9}
F = fd.FREQ
FSUP = np.concatenate([[1e3], F[np.unique(np.r_[np.arange(0, 80, 6), 79])]])
WIT = [f"orb0.1_m3L{L}{s}" for L in (30, 100, 300) for s in ("max", "min")] + \
      [f"{c}_m3L100{s}" for c in ("cal1nz", "cal3wnz", "cal3w") for s in ("max", "min")]


def log(*a):
    print(*a, f"[{time.perf_counter() - T0:.0f}s]", flush=True)


# ------------------------------------------------------------------ memory guard (shared machine: abort, never compete)
import ctypes  # noqa: E402

FREE_MIN_GB = 1.5      # abort a heavy step if free physical memory is below this
OWN_MAX_GB = 1.9       # abort if this process's private memory exceeds this


class _MS(ctypes.Structure):
    _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]


class _PMC(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t), ("PrivateUsage", ctypes.c_size_t)]


def mem_state():
    ms = _MS()
    ms.dwLength = ctypes.sizeof(ms)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
    pmc = _PMC()
    pmc.cb = ctypes.sizeof(pmc)
    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    ps = ctypes.windll.psapi
    ps.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_ulong]
    ps.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
    return ms.ullAvailPhys / 2 ** 30, pmc.PrivateUsage / 2 ** 30, pmc.PeakWorkingSetSize / 2 ** 30


def mem_guard(tag):
    free, own, peak = mem_state()
    log(f"memory before {tag}: free {free:.2f} GB, own private {own:.2f} GB, own peak WS {peak:.2f} GB")
    if free < FREE_MIN_GB or own > OWN_MAX_GB:
        log(f"ABORT (memory guard) before {tag}: free {free:.2f} GB < {FREE_MIN_GB} or own {own:.2f} GB > {OWN_MAX_GB}")
        sys.exit(3)


def load_res():
    return json.loads(RES.read_text()) if RES.exists() else {}


def save(res):
    RES.write_text(json.dumps(res, indent=1, default=float))


def regulated_QDQ(C, e, vz, L):
    d = np.abs(vz) / L
    eh = e / np.linalg.norm(e)
    h = d * eh
    Ch = C.copy()
    Ch[np.diag_indices(e.size)] += d
    Ch -= np.outer(eh, h)
    Ch -= np.outer(h, eh)
    Ch += (eh @ h) * np.outer(eh, eh)
    return Ch


def reg_times(C, e, vz, L):
    mem_guard(f"bordered solve L = {L * 1e9:.1f} nm")
    out = sc.kinetic_bulk(C, e, vz, L=L)
    gc.collect()
    return out


def run_slab(G, name, d, freqs, tag):
    e, v = G["e"], G["vz"]
    mem_guard(f"operator assembly {name}")
    C = sc.operator(G, "film2", name)
    Kp, D = sc.build_Kp(C, e, inplace=True)
    del C
    for k in ("P", "A", "B", "rep", "aP", "aA", "aB", "orb"):     # event arrays no longer needed
        G.pop(k, None)
    gc.collect()
    path = HERE / "out" / f"sa4_m996_{name}_{d * 1e9:.0f}.npz"
    have = {}
    if path.exists():
        Z = np.load(path)
        have = {float(f): (Z["Y"][i], Z["Th"][i]) for i, f in enumerate(Z["f"])}
    Ys, Ths, fs, diags = [], [], [], []
    for f in freqs:
        if float(f) in have:
            Y, Th = have[float(f)]
        else:
            mem_guard(f"slab {name} {d * 1e9:.0f} nm {f:.4g} Hz")
            z = sc.make_mesh(d, level=0, **MESH)
            Y, Th, dg, sol = sc.twoport(None, v, e, 2j * np.pi * f, z, rtol=2e-12, restart=15, maxiter=20, order=4,
                                        Kp=Kp, D=D)
            del sol
            gc.collect()
            diags.append({"f": float(f), "N": dg["N"], "ebal": dg["ebal"], "iters": [g["iters"] for g in dg["gmres"]],
                          "relres": [g["true_relres"] for g in dg["gmres"]], "t": dg["t_total"],
                          "recip": dg["reciprocity"], "symm": dg["symmetry"]})
            log(name, d * 1e9, f"{f:.4g}", diags[-1])
        Ys.append(Y)
        Ths.append(Th)
        fs.append(float(f))
        np.savez(path, f=np.array(fs), Y=np.array(Ys), Th=np.array(Ths))
    # discretisation-error estimate with a coarser mesh at the band edges
    err = []
    for f in (F[0], F[-1]):
        mem_guard(f"coarse-mesh check {f:.4g} Hz")
        z = sc.make_mesh(d, level=0, **MESH_COARSE)
        Yc, Thc, dg, sol = sc.twoport(None, v, e, 2j * np.pi * f, z, rtol=2e-12, restart=15, maxiter=20, order=4,
                                      Kp=Kp, D=D)
        del sol
        i = fs.index(float(f))
        err.append({"f": float(f), "N_coarse": dg["N"],
                    "rel_diff_Y": float(np.max(np.abs(Yc - Ys[i])) / np.max(np.abs(Ys[i]))),
                    "rel_diff_ImY11": float(abs((Yc[0, 0] - Ys[i][0, 0]).imag) / abs(Ys[i][0, 0].imag))})
        log("coarse-mesh check", err[-1])
    del Kp
    gc.collect()
    return diags, err


def at_fdtr(fs, Y):
    out = np.zeros((F.size, 2, 2), complex)
    for a in range(2):
        for b in range(2):
            out[:, a, b] = AAA(2j * np.pi * np.asarray(fs), Y[:, a, b], rtol=1e-14)(2j * np.pi * F)
    return out


def main():
    stage = sys.argv[1]
    res = load_res()
    G = sc.load_modes("m996")
    for k in ("P", "A", "B") if stage not in ("op", "proxy") else ():
        pass
    e, v = G["e"], G["vz"]
    log("m996 modes", G["n"], "orbits", G["n_orb"], "missing", G["missing_images"],
        "moving", int(np.sum(np.abs(v) >= 1e-12)), "local", int(np.sum(np.abs(v) < 1e-12)))
    if stage == "op":
        C = sc.operator(G, "film2", "REF")
        ent = {"n": G["n"], "n_orb": G["n_orb"], "missing_images": G["missing_images"], "C_tot": G["C_tot"],
               "energy_residual": float(np.linalg.norm(C @ e) / np.linalg.norm(np.diag(C) * e))}
        ent["bulk"] = sc.kinetic_bulk(C, e, v)
        for L in (300e-9, 500e-9):
            ent[f"reg_{L * 1e9:.0f}"] = reg_times(C, e, v, L)
        res["op_REF"] = ent
        save(res)
        log("op", ent)
    elif stage in ("slab", "slab300"):
        d = 500e-9 if stage == "slab" else 300e-9
        diags, err = run_slab(G, "REF", d, FSUP, stage)
        res.setdefault("slab_REF", {})[f"{d * 1e9:.0f}"] = {"diags": diags, "coarse_mesh_check": err}
        save(res)
    elif stage == "wit":
        name = sys.argv[2]
        d = float(sys.argv[3]) * 1e-9 if len(sys.argv) > 3 else 500e-9
        diags, err = run_slab(G, name, d, FSUP, "wit")
        res.setdefault("slab_wit", {})[f"{name}|{d * 1e9:.0f}"] = {"diags": diags, "coarse_mesh_check": err}
        save(res)
    elif stage == "fit":
        out = res.get("fit", {})
        op = res["op_REF"]
        for path in sorted((HERE / "out").glob("sa4_m996_*_*.npz")):
            name, dnm = path.stem.split("_", 2)[2].rsplit("_", 1)
            Z = np.load(path)
            if Z["f"].size < FSUP.size:
                continue
            d = float(dnm) * 1e-9
            kr = op[f"reg_{dnm}"]
            base = np.array([kr["tauR"], kr["taul"], np.log(kr["K"]), np.log(G["C_tot"]), np.log(1e-8), np.log(2e-8)])
            Yraw = np.asarray(Z["Y"])
            # reciprocity- and mirror-symmetrised two-port (the q_z = 1/2 plane breaks both at 1e-6 ... 1e-3)
            Ysym = Yraw.copy()
            Ysym[:, 0, 0] = 0.5 * (Yraw[:, 0, 0] - Yraw[:, 1, 1])
            Ysym[:, 1, 1] = -Ysym[:, 0, 0]
            Ysym[:, 0, 1] = 0.5 * (Yraw[:, 0, 1] - Yraw[:, 1, 0])
            Ysym[:, 1, 0] = -Ysym[:, 0, 1]
            ent = {"G0_1kHz": float(Yraw[0][0, 0].real), "lam_eff_1kHz": float(Yraw[0][0, 0].real * d),
                   "max_recip_rel": float(np.max(np.abs(Yraw[:, 0, 1] + Yraw[:, 1, 0]) / np.abs(Yraw[:, 0, 0]))),
                   "max_symm_rel": float(np.max(np.abs(Yraw[:, 0, 0] + Yraw[:, 1, 1]) / np.abs(Yraw[:, 0, 0])))}
            for vtag, Yv in (("", Yraw), ("sym|", Ysym)):
                Yf = at_fdtr(Z["f"], Yv)
                y = np.log(fd.Z_slab(Yf, 1e-8, 2e-8, 2j * np.pi * F, port="wall"))
                for N in ("N2lam", "N3"):
                    for sp in (0.01, 0.1):
                        r = fd.fit_delta(y, d, base, sp, N)
                        ent[f"{vtag}{N}|{sp}"] = {"delta_ps": r["delta_ps"], "chi2": r["chi2"], "par": r["par"],
                                                  "n_at_min": r["n_at_min"], "spread_ps": r["spread_ps"]}
                        if name == "REF" and sp == 0.01 and vtag == "":
                            J = fd.cauchy_jac(np.array(r["par"]), d)
                            ent[f"sigma|{N}"] = {f"{NN}|{s}": fd.sigma_direction(J, 0, fd.SETCOLS[NN], s) * 1e12
                                                 for NN in ("N1", "N2", "N2lam", "N3") for s in (0.01, 0.03, 0.1)}
            out[f"{name}|{dnm}"] = ent
            log("fit", name, dnm, {k: round(v_["delta_ps"], 3) for k, v_ in ent.items() if isinstance(v_, dict) and "delta_ps" in v_})
        res["fit"] = out
        save(res)
    elif stage == "proxy":
        out = {}
        mem_guard("proxy REF operator")
        C = sc.operator(G, "film2", "REF")
        for dnm in ("300", "500"):
            if f"REF|{dnm}" not in res.get("fit", {}):
                continue
            d = float(dnm) * 1e-9
            for N in ("N2lam", "N3"):
                target = res["fit"][f"REF|{dnm}"][f"{N}|0.01"]["delta_ps"]
                c = brentq(lambda c_: reg_times(C, e, v, d / c_)["delta"] * 1e12 - target, 1.0, 8.0, xtol=1e-5)
                out[f"c|{dnm}|{N}"] = c
                log("proxy c", dnm, N, c)
        del C
        gc.collect()
        prox = {}
        for nm in ["REF"] + WIT:
            mem_guard(f"proxy operator {nm}")
            C = sc.operator(G, "film2", nm)
            prox[nm] = {}
            for key, c in out.items():
                _, dnm, N = key.split("|")
                prox[nm][f"{dnm}|{N}"] = reg_times(C, e, v, float(dnm) * 1e-9 / c)["delta"] * 1e12
            del C
            gc.collect()
            log("proxy", nm, {k: round(x, 3) for k, x in prox[nm].items()})
            res["proxy"] = {"c": out, "values": prox}
            save(res)
        summ = {}
        for key in prox["REF"]:
            ex = {w: prox[w][key] - prox["REF"][key] for w in WIT}
            up, low = max(list(ex.values()) + [0.0]), min(list(ex.values()) + [0.0])
            summ[key] = {"half_range_ps": 0.5 * (up - low), "up": up, "low": low, "argmax": max(ex, key=ex.get),
                         "argmin": min(ex, key=ex.get), "excursions": ex}
        res["proxy"]["ranges"] = summ
        save(res)
        log("proxy ranges", {k: round(v_["half_range_ps"], 3) for k, v_ in summ.items()})


if __name__ == "__main__":
    main()
