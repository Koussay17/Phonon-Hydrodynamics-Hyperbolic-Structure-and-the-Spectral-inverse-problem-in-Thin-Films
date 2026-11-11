"""Practical (band-limited) conditioning: singular values of the 4-parameter problem in a finite band.

Band [omega_b / 2e4, omega_b] (the 4.3 decades of 10 kHz - 200 MHz), slid across the film scales.
Metric: RMS-normalised singular values  s_hat_k = s_k / sqrt(L),  L = ln(2e4),  of the relative
Jacobian in log parameters with weights d ln(omega): a unit step of ln(theta) along the k-th right
singular vector changes Z by s_hat_k in RMS relative value over the band.  The standard deviation
of that combination at RMS relative noise sigma per frequency (n frequencies) is ~ sigma/(s_hat_k sqrt n).
Asymptotic laws fitted: thin film (omega_b xi1^2 -> 0), thick film (omega_a xi1^2 -> inf),
low-frequency non-Fourier (omega_b tau -> 0 at fixed omega_b xi1^2).
"""
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zmodel as zm

RATIO = 2e4
L = np.log(RATIO)


def band_svd(theta, wb, R=None):
    xi1, beta, tR, tl = theta
    w, wts, ok = zm.resolved_grid(xi1, beta, tR, tl, 1.0, wmin=wb / RATIO, wmax=wb, nmax=60000)
    M = zm.real_rel_jac_np(w, xi1, beta, tR, tl, 1.0, 0.0, "code", False, "log", wts)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    return S / np.sqrt(L), Vt


def fit(x, y):
    k, c = np.polyfit(np.log(x), np.log(y), 1)
    return float(k), float(np.exp(c))


if __name__ == "__main__":
    out = {}
    shapes = {"gk_tR>tl": (2.0, 0.5, 0.1), "gk_tl>tR": (2.0, 0.1, 0.5), "cattaneo": (2.0, 0.5, 0.0),
              "gk0": (2.0, 0.0, 0.5), "beta0.5_gk": (0.5, 0.5, 0.1)}
    # 1. slide the band across the film scale (times fixed relative to xi1^2)
    wbs = 10.0 ** np.arange(-5, 6.01, 0.5)
    for name, (beta, tR, tl) in shapes.items():
        rows = []
        for wb in wbs:
            S, Vt = band_svd((1.0, beta, tR, tl), wb)
            rows.append([wb, *S.tolist(), *Vt[-1].tolist()])
        rows = np.array(rows)
        out[f"slide_{name}"] = rows.tolist()
        print(f"\n[{name}] beta={beta} tR={tR} tl={tl}: band top omega_b xi1^2, s_hat_1..4, null vector")
        for r in rows[::2]:
            print("  wb=%8.1e  s=%s  v4=%s" % (r[0], np.array2string(r[1:5], precision=2),
                                                np.array2string(r[5:], precision=2)))
        thin = rows[:, 0] <= 1e-3
        if thin.sum() >= 3:
            out[f"thin_fit_{name}"] = {f"s{k+1}": fit(rows[thin, 0], rows[thin, 1 + k]) for k in range(4)}
            print("  thin-film fits s_k ~ c (omega_b xi1^2)^k:",
                  {kk: "%.3f" % v[0] for kk, v in out[f"thin_fit_{name}"].items()})
    # 2. low-frequency non-Fourier regime: fixed band position, times -> 0
    print("\n[low omega tau] band top omega_b xi1^2 = 2, times scaled by eps, tR/tl = 5, beta = 2")
    rows = []
    for eps in 10.0 ** np.arange(0, -4.01, -0.5):
        S, Vt = band_svd((1.0, 2.0, 0.5 * eps, 0.1 * eps), 2.0)
        rows.append([eps, *S.tolist()])
        print("  omega_b tR = %.1e: s=%s" % (2.0 * 0.5 * eps, np.array2string(S, precision=3)))
    rows = np.array(rows)
    out["lowtau_rows"] = rows.tolist()
    out["lowtau_fit"] = {f"s{k+1}": fit(rows[3:, 0], rows[3:, 1 + k]) for k in range(4)}
    print("  fits s_k ~ eps^k:", {kk: "%.3f" % v[0] for kk, v in out["lowtau_fit"].items()})
    # 3. thick film: band bottom omega_a xi1^2 -> inf, s4 vs exp(-2 Re A(omega_a))
    print("\n[thick] beta=2, tR=0.5, tl=0.1 (times in units of xi1^2), band bottom omega_a xi1^2")
    rows = []
    for wa in 10.0 ** np.arange(-1, 2.6, 0.25):
        S, Vt = band_svd((1.0, 2.0, 0.5, 0.1), wa * RATIO)
        A = zm.pieces_np(np.array([wa]), 1.0, 2.0, 0.5, 0.1)["A"][0]
        rows.append([wa, *S.tolist(), float(np.exp(-2 * A.real))])
        print("  omega_a=%.2e: s=%s  exp(-2ReA(omega_a))=%.2e" % (wa, np.array2string(S, precision=2), rows[-1][-1]))
    rows = np.array(rows)
    out["thick_rows"] = rows.tolist()
    m = rows[:, -1] < 1e-3
    if m.sum() >= 3:
        k, c = fit(rows[m, -1], rows[m, 4])
        out["thick_fit_s4_vs_E"] = [k, c]
        print("  fit s4 ~ c exp(-2 Re A)^k: k = %.3f" % k)
    with open(os.path.join(HERE, "data", "practical_bands.json"), "w") as fh:
        json.dump(out, fh, indent=1)
