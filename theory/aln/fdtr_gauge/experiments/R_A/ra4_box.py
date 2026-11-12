"""R-A4: sensitivity box of the workstream-A verdict (objection A-1, minor A-m1, A-m3).

sigma(ln tau_R) / w(L) over a declared box, finite beams 10/8 um, Cattaneo AlN layer (A4 set-up, 80 log-spaced
frequencies 10 kHz - 200 MHz, sigma(ln|H|) = 1 %):
  R1 = R_Au/AlN, R2 = R_AlN/sapphire in {2, 3, 5, 7, 10, 14, 20} x 1e-9 m^2K/W (7 x 7),
  regulator normalisation c in {1/2, 2/3, 3/4, 1, 4/3, 3/2, 2}: base values at (L, c) from ra2 (event operator at
  c = 1 times the RTA 64x64x35 c-dependence), width w(L; c) = w(L/c) from ra3 (inner bounds),
  lam_z variants: 'event' (A1 base), 'rta' (tau_R and lam_z both from RTA 64x64x35), 'lam-20%', 'lam+20%',
  L in {30, 100, 300, 500} nm, sigma_phi in {0.01, 0.03, 0.1} deg, nuisance sets N1, N2, N3 (A4 definitions).
Forward model: the independent quadrupole model of branch A (experiments/A/indep.py: surface_impedance, Gauss-
Legendre Hankel rule), base values passed as traced arguments, exact forward-mode Jacobian (JAX, float64).
Marginal sigma by orthogonal projection of the whitened tau_R column on the nuisance span (Householder QR), not the
branch's SVD routine. Second implementation: the twin (src/twin_fdtr.py via experiments/A/common.py, Richardson
central differences, h = 2e-3) at the 32 corners R1, R2 in {2e-9, 2e-8}, c in {1/2, 2}, all L ('event' variant),
and the A4 base point (reproduction of A4).
Information measures (A-m1): ra_info.info(sigma / w) -> bits, prior variance removed.
Au-thickness sub-box (referee A-1 (iii)): d_Au in {20, 50, 80, 100} nm at c = 1, 'event', full (R1, R2) grid.
Outputs: out/ra4_box.npz (full array), out/ra4_box.json (summaries, checks).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ["XLA_FLAGS"] = "--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=2"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "2"
import numpy as np  # noqa: E402
import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
ADIR = HERE.parent / "A"
sys.path.insert(0, str(ADIR))
sys.path.insert(1, str(HERE))
from indep import surface_impedance, gl_rule  # noqa: E402
import common  # noqa: E402  (twin, A4 constants; read-only)
from ra_info import info  # noqa: E402

FREQ = common.FREQ
SIG_A = common.SIG_A
R_GRID = np.array([2e-9, 3e-9, 5e-9, 7e-9, 1e-8, 1.4e-8, 2e-8])
CS = (0.5, 2 / 3, 0.75, 1.0, 4 / 3, 1.5, 2.0)
LS = (30, 100, 300, 500)
VARS = ("event", "rta", "lam-20%", "lam+20%")
SIGPHI = (0.01, 0.03, 0.1)
SETS = {"N1": [0, 5, 6, 7], "N2": [0, 5, 6, 7, 3, 4], "N3": [0, 5, 6, 7, 3, 4, 1, 2]}
THN = ("ln_tauR", "ln_lamz", "ln_C", "ln_R1", "ln_R2", "g", "phi0", "td_ns")
U, WU = gl_rule()
U, WU = jnp.asarray(U), jnp.asarray(WU)
F = jnp.asarray(FREQ)


def lnH_real(th, bv):
    """[Re, Im] of ln H_obs; bv = (tauR, lamz, lamr, C, R1, R2, dAlN, lamAu, CAu, dAu, lamS, CS, wp, wr)."""
    e = jnp.exp
    P = {"tauR": bv[0] * e(th[0]), "taur_ip": bv[0] * e(th[0]), "lamz": bv[1] * e(th[1]), "lamr": bv[2],
         "C": bv[3] * e(th[2]), "R1": bv[4] * e(th[3]), "R2": bv[5] * e(th[4]), "dAlN": bv[6], "lamAu": bv[7],
         "CAu": bv[8], "dAu": bv[9], "lamS": bv[10], "CS": bv[11]}
    p = 2j * jnp.pi * F
    eta = (bv[12] ** 2 + bv[13] ** 2) / 8.0
    k2 = (U / jnp.sqrt(eta)) ** 2
    Z = surface_impedance(p[:, None], k2[None, :], P)
    H = (Z * (U * jnp.exp(-U * U) * WU)[None, :]).sum(axis=1) / (2.0 * jnp.pi * eta)
    z = jnp.log(H) + th[5] + 1j * (th[6] - 2 * jnp.pi * F * th[7] * 1e-9)
    return jnp.concatenate([z.real, z.imag])


jac_batch = jax.jit(jax.vmap(jax.jacfwd(lnH_real), in_axes=(None, 0)))
TH0 = jnp.zeros(8)


def sigma_marg(J, cols, sphi):
    """sigma of column 0 of cols, whitened rows [Re/sig_A ; Im/sig_phi(rad)], nuisances = cols[1:]."""
    nf = FREQ.size
    Jw = np.vstack([J[:nf, cols] / SIG_A, J[nf:, cols] / np.deg2rad(sphi)])
    jt, Jn = Jw[:, 0], Jw[:, 1:]
    sc = np.linalg.norm(Jn, axis=0)
    Q, _ = np.linalg.qr(Jn / sc)
    res = jt - Q @ (Q.T @ jt)
    res = res - Q @ (Q.T @ res)                      # re-orthogonalisation
    return 1.0 / np.linalg.norm(res)


def base_vec(L, c, var, R1, R2, bb):
    b = bb[f"{L}|{c:.6g}"]
    tau, lamz = b["tau_R_ps"] * 1e-12, b["lam_z"]
    if var == "rta":
        tau, lamz = b["tau_R_ps_rta"] * 1e-12, b["lam_z_rta"]
    elif var == "lam-20%":
        lamz = 0.8 * lamz
    elif var == "lam+20%":
        lamz = 1.2 * lamz
    return [tau, lamz, b["lam_r"], common.C_ALN, R1, R2, L * 1e-9, common.AU_LAM, common.AU_C, common.AU_D,
            common.SAP_LAM, common.SAP_C, common.WP, common.WR]


def main():
    t0 = time.perf_counter()
    bb = json.loads((OUT / "ra2_rta_film_c.json").read_text())["base_box"]
    wd = json.loads((OUT / "ra3_widths.json").read_text())["per_L"]
    key = lambda L: f"{L:.6g}"  # noqa: E731
    ia_, ib_ = int(np.argmin(abs(R_GRID - 1e-8))), int(np.argmin(abs(R_GRID - 2e-8)))
    shape = (len(LS), len(CS), len(VARS), len(R_GRID), len(R_GRID), len(SETS), len(SIGPHI))
    sig = np.full(shape, np.nan)
    W = {cl: np.zeros((len(LS), len(CS))) for cl in ("orb0.1", "all")}
    W12 = {cl: np.zeros((len(LS), len(CS))) for cl in ("orb0.1", "all")}
    for i, L in enumerate(LS):
        for j, c in enumerate(CS):
            row = wd[key(round(L / c, 6))]
            for cl in W:
                W[cl][i, j] = row[cl]["w_half"]
                W12[cl][i, j] = row[cl]["w_half_m12128"]
    for i, L in enumerate(LS):
        for j, c in enumerate(CS):
            for k, var in enumerate(VARS):
                for a, R1 in enumerate(R_GRID):                             # batches of 7 (memory)
                    bvs = np.array([base_vec(L, c, var, R1, R2, bb) for R2 in R_GRID])
                    Jb = np.asarray(jac_batch(TH0, jnp.asarray(bvs)))       # (7, 160, 8)
                    for bidx, J in enumerate(Jb):
                        for s, cols in enumerate(SETS.values()):
                            for t, sp in enumerate(SIGPHI):
                                sig[i, j, k, a, bidx, s, t] = sigma_marg(J, cols, sp)
        import psutil
        mi = psutil.Process().memory_info()
        print(f"L={L} done {time.perf_counter() - t0:.0f}s, rss {mi.rss / 2**30:.2f} GB, "
              f"peak {getattr(mi, 'peak_wset', mi.rss) / 2**30:.2f} GB", flush=True)
    out = {"axes": {"L_nm": LS, "c": CS, "lam_variant": VARS, "R1": R_GRID.tolist(), "R2": R_GRID.tolist(),
                    "set": list(SETS), "sigma_phi_deg": SIGPHI}, "checks": {}}
    # ---------------- Au-thickness sub-box (referee A-1 (iii)): c = 1, 'event', full R grid, d_Au in DAU
    DAU = (20e-9, 50e-9, 80e-9, 100e-9)
    sig_au = np.full((len(LS), len(DAU), len(R_GRID), len(R_GRID), len(SETS), len(SIGPHI)), np.nan)
    for i, L in enumerate(LS):
        for j, dau in enumerate(DAU):
            for a, R1 in enumerate(R_GRID):
                bvs = np.array([base_vec(L, 1.0, "event", R1, R2, bb) for R2 in R_GRID])
                bvs[:, 9] = dau
                Jb = np.asarray(jac_batch(TH0, jnp.asarray(bvs)))
                for bidx, J in enumerate(Jb):
                    for s, cols in enumerate(SETS.values()):
                        for t, sp in enumerate(SIGPHI):
                            sig_au[i, j, a, bidx, s, t] = sigma_marg(J, cols, sp)
    au = {}
    for i, L in enumerate(LS):
        w0 = W["all"][i, CS.index(1.0)]
        wo = W["orb0.1"][i, CS.index(1.0)]
        for s, N in enumerate(SETS):
            for t, sp in enumerate(SIGPHI):
                blk = sig_au[i, :, :, :, s, t]
                idx = np.unravel_index(np.nanargmin(blk), blk.shape)
                au[f"{L}|{N}|{sp}"] = {"min_sigma": float(blk[idx]), "at": {"dAu_nm": DAU[idx[0]] * 1e9,
                                       "R1": float(R_GRID[idx[1]]), "R2": float(R_GRID[idx[2]])},
                                       "min_ratio_all": float(blk[idx] / w0), "min_ratio_orb": float(blk[idx] / wo),
                                       "sigma_80nm_note16": float(sig_au[i, DAU.index(80e-9), ia_, ib_, s, t])}
    out["au_subbox"] = au
    np.save(OUT / "ra4_au_subbox.npy", sig_au)
    print("Au sub-box done", flush=True)
    # ---------------- reproduction of A4 (event, c = 1, R1 = 1e-8, R2 = 2e-8)
    a4 = {(500, "N2", 0.01): 0.340379, (500, "N1", 0.01): 0.205122, (500, "N3", 0.01): 2.25388,
          (500, "N2", 0.1): 2.73316, (300, "N2", 0.01): 0.662085, (300, "N2", 0.1): 5.36486,
          (100, "N2", 0.1): 47.4065, (30, "N2", 0.1): 1259.692}
    jc, ia, ib = CS.index(1.0), int(np.argmin(abs(R_GRID - 1e-8))), int(np.argmin(abs(R_GRID - 2e-8)))
    out["checks"]["A4_reproduction_rel"] = {
        f"{L}|{N}|{sp}": float(sig[LS.index(L), jc, 0, ia, ib, list(SETS).index(N), SIGPHI.index(sp)] / v - 1)
        for (L, N, sp), v in a4.items()}
    # ---------------- twin cross-check at the corners
    names = ["ln_tauR", "g", "phi0", "td_ns", "ln_R1", "ln_R2", "ln_lamz", "ln_C"]
    tw_cols = {"N1": [0, 1, 2, 3], "N2": [0, 1, 2, 3, 4, 5], "N3": [0, 1, 2, 3, 4, 5, 6, 7]}
    worst = 0.0
    corner = {}
    for L in LS:
        for c in (0.5, 2.0):
            for R1 in (2e-9, 2e-8):
                for R2 in (2e-9, 2e-8):
                    b = bb[f"{L}|{c:.6g}"]
                    common.R1, common.R2 = R1, R2
                    m = common.Model(L_nm=float(L), tauR=b["tau_R_ps"] * 1e-12, lamz=b["lam_z"], lamr=b["lam_r"])
                    Jt = common.jacobian(m, FREQ, names, h=2e-3, scheme="richardson")
                    Jtr = np.vstack([Jt.real, Jt.imag])
                    i, j = LS.index(L), CS.index(c)
                    a, bidx = list(R_GRID).index(R1), list(R_GRID).index(R2)
                    for s, N in enumerate(SETS):
                        for t, sp in enumerate(SIGPHI):
                            st = sigma_marg(Jtr, tw_cols[N], sp)
                            sj = sig[i, j, 0, a, bidx, s, t]
                            d = abs(st / sj - 1)
                            worst = max(worst, d)
                            corner[f"{L}|{c:.3g}|{R1:.0e}|{R2:.0e}|{N}|{sp}"] = [float(sj), float(st), float(d)]
    common.R1, common.R2 = 1e-8, 2e-8
    out["checks"]["twin_vs_indep_corners_max_rel"] = worst
    out["checks"]["twin_vs_indep_corners"] = corner
    print("twin corner check max rel", worst, flush=True)
    # ---------------- ratios, information, summaries
    summ = {}
    for cl, Wc in W.items():
        for basis, Wuse in (("robust", Wc), ("m12128", W12[cl])):
            ratio = sig / Wuse[:, :, None, None, None, None, None]
            tag = f"{cl}|{basis}"
            for i, L in enumerate(LS):
                for s, N in enumerate(SETS):
                    for t, sp in enumerate(SIGPHI):
                        blk = ratio[i, :, :, :, :, s, t]
                        idx = np.unravel_index(np.nanargmin(blk), blk.shape)
                        rmin = float(blk[idx])
                        bits, vr = (info(rmin)["bits"], info(rmin)["var_removed"])
                        flip = blk <= 1.0
                        flip2 = blk <= 0.5
                        reg = {}
                        if flip.any():
                            jj, kk, aa, bb_ = np.nonzero(flip)
                            reg = {"c": sorted({round(CS[x], 4) for x in jj}),
                                   "var": sorted({VARS[x] for x in kk}),
                                   "R1_max": float(R_GRID[aa].max()), "R2": sorted({float(R_GRID[x]) for x in bb_})}
                        base = float(ratio[i, CS.index(1.0), 0, ia, ib, s, t])
                        summ[f"{tag}|{L}|{N}|{sp}"] = {
                            "min_ratio": rmin, "at": {"c": CS[idx[0]], "var": VARS[idx[1]], "R1": float(R_GRID[idx[2]]),
                                                      "R2": float(R_GRID[idx[3]])},
                            "max_ratio": float(np.nanmax(blk)), "base_ratio": base,
                            "bits_at_min": bits, "var_removed_at_min": vr,
                            "bits_at_base": info(base)["bits"], "var_removed_at_base": info(base)["var_removed"],
                            "frac_flip_1sigma": float(flip.mean()), "frac_flip_2sigma": float(flip2.mean()),
                            "flip_region": reg}
    out["summary"] = summ
    out["widths"] = {cl: W[cl].tolist() for cl in W}
    out["widths_m12128"] = {cl: W12[cl].tolist() for cl in W12}
    out["seconds"] = time.perf_counter() - t0
    np.savez_compressed(OUT / "ra4_box.npz", sigma=sig, w_orb=W["orb0.1"], w_all=W["all"],
                        w12_orb=W12["orb0.1"], w12_all=W12["all"], L=np.array(LS), c=np.array(CS),
                        R=R_GRID, sigphi=np.array(SIGPHI))
    (OUT / "ra4_box.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["checks"]["A4_reproduction_rel"], indent=1))
    for k, v in summ.items():
        if k.startswith("orb0.1|robust") or k.startswith("all|robust"):
            print(k, f"min {v['min_ratio']:.3g} at {v['at']}  base {v['base_ratio']:.3g}  max {v['max_ratio']:.3g} "
                     f"flip1 {v['frac_flip_1sigma']:.3f} bits@min {v['bits_at_min']:.3g}")
    print("total", round(out["seconds"]), "s")


if __name__ == "__main__":
    main()
