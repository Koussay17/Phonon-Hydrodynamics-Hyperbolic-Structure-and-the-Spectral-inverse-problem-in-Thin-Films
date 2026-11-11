"""Figures of workstream A from the saved JSON (no recomputation).

figures/A3_dphi.png   |Delta phi(f)| at the admissible half-width w(L) for 30/100/300/500 nm against the per-point
                      phase noise; signed Delta phi(f) at 500 nm for w = +-2.8, 5.2, 8 %.
figures/A4_wstar.png  w*_1sigma = sigma(ln tau_R) per nuisance set and thickness (Cattaneo, beams 10/8 um, full band)
                      at sigma_phi = 0.01 and 0.1 deg, against the admissible half-widths.
figures/A5_profile.png noiseless profile Delta chi2(w) against the Fisher parabola (500 nm), if a5 output exists.
Palette: validated categorical slots 1-4 (light surface), direct labels as the contrast relief.
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from common import OUT, FIG, W_BY_L  # noqa: E402

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd5"
COL = {30: "#2a78d6", 100: "#eb6834", 300: "#1baf7a", 500: "#eda100"}
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "font.size": 9, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "lines.linewidth": 1.6, "legend.frameon": False})


def fig_dphi():
    d = json.loads((OUT / "a3_dphi.json").read_text())
    f = np.array(d["freq_Hz"])
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    wL = {30: 0.028, 100: 0.037, 300: 0.052, 500: 0.052}
    for L in (30, 100, 300, 500):
        c = d["cases"][f"L{L}_cat_beam"]["w"][f"{wL[L]:+.3f}"]
        y = np.abs(np.array(c["dphi_deg"]))
        ax[0].loglog(f, y, color=COL[L], label=f"{L} nm, w = {100 * wL[L]:.1f} %")
        i = int(np.argmax(y))
        off = {30: (4, 4), 100: (4, 4), 300: (6, -12), 500: (-46, 6)}[L]
        ax[0].annotate(f"{L} nm", (f[i], y[i]), textcoords="offset points", xytext=off, color=INK, fontsize=8)
    for s, ls in ((0.01, "--"), (0.03, "-."), (0.1, ":")):
        ax[0].axhline(s, color=INK2, ls=ls, lw=1.0)
        ax[0].text(1.2e4, s * 1.15, f"phase noise {s}°", color=INK2, fontsize=7.5)
    ax[0].set_ylim(1e-7, 0.3)
    ax[0].set_xlabel("frequency (Hz)")
    ax[0].set_ylabel("|Δφ| (deg)")
    ax[0].set_title("(a) Cattaneo AlN, beams 10/8 µm, admissible w(L)", fontsize=9, loc="left")
    ax[0].legend(fontsize=7.5, loc="lower right")
    case = d["cases"]["L500_cat_beam"]["w"]
    for w, ls in ((0.028, ":"), (0.052, "--"), (0.08, "-")):
        for sgn in (+1, -1):
            y = np.array(case[f"{sgn * w:+.3f}"]["dphi_deg"]) * 1e3
            ax[1].semilogx(f, y, color=COL[500] if sgn > 0 else COL[30], ls=ls,
                           label=f"w = {sgn * 100 * w:+.1f} %")
    ax[1].set_xlabel("frequency (Hz)")
    ax[1].set_ylabel("Δφ (millidegree)")
    ax[1].set_title("(b) 500 nm, τ_R = 25.0 ps: signed Δφ", fontsize=9, loc="left")
    ax[1].legend(fontsize=7.5, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "A3_dphi.png", dpi=170)
    plt.close(fig)


def fig_wstar():
    d = json.loads((OUT / "a4_fisher.json").read_text())
    sets = ["N0", "N1", "N2", "N3", "N4", "N5"]
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
    for k, s in enumerate(("0.01", "0.1")):
        for j, L in enumerate((30, 100, 300, 500)):
            c = d["cases"][f"L{L}_cat_aniso_beam"]["sets"]
            y = [c[N]["10k-200M"][s]["sigma"] for N in sets]
            x = np.arange(len(sets)) + (j - 1.5) * 0.08
            ax[k].semilogy(x, y, "o-", color=COL[L], ms=4.5, label=f"{L} nm")
            ax[k].annotate(f"{L} nm", (x[-1], y[-1]), textcoords="offset points", xytext=(5, 0), color=INK,
                           fontsize=7.5, va="center")
        for L, w in ((30, 0.028), (100, 0.037), (300, 0.052)):
            ax[k].axhline(w, color=COL[L], ls="--", lw=1.0)
        ax[k].axhline(0.08, color=INK2, ls=":", lw=1.0)
        ax[k].text(-0.3, 0.085, "8 % (100 nm, calibrated classes)", fontsize=7, color=INK2)
        ax[k].text(-0.3, 0.020, "dashed: admissible w(L)", fontsize=7, color=INK2)
        ax[k].set_xticks(range(len(sets)), sets)
        ax[k].set_xlim(-0.4, len(sets) - 0.3)
        ax[k].set_title(f"({'ab'[k]}) σ_φ = {s}°, σ_ln|H| = 1 %, 10 kHz–200 MHz", fontsize=9, loc="left")
        ax[k].set_xlabel("nuisance set")
    ax[0].set_ylabel("w*₁σ = σ(ln τ_R)")
    ax[0].legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "A4_wstar.png", dpi=170)
    plt.close(fig)


def fig_profile():
    p = OUT / "a5_profile.json"
    if not p.exists():
        return
    d = json.loads(p.read_text())
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=False)
    cols = {"N0": "#2a78d6", "N1": "#eb6834", "N2": "#1baf7a", "N3": "#eda100"}
    for k, s in enumerate(("0.01", "0.1")):
        for N in ("N0", "N1", "N2", "N3"):
            key = f"{N}|{s}"
            if key not in d["asimov"]:
                continue
            pr = d["asimov"][key]["profile"]
            w = np.array([q["w"] for q in pr])                     # relative change of tau_R; -1 = Fourier film
            y = np.array([q["dchi2"] for q in pr])
            o = np.argsort(w)
            ax[k].plot(w[o], y[o], "o", color=cols[N], ms=4.5, label=f"{N}: profile (dots), (w/σ_F)² (line)")
            xx = np.linspace(-1.0, 3.0, 400)
            ax[k].plot(xx, (xx / d["asimov"][key]["sigma_F"]) ** 2, color=cols[N], lw=1.0)
        ax[k].set_yscale("log")
        ax[k].set_ylim(1e-6, 1e3)
        ax[k].axhline(1, color=INK2, ls="--", lw=1.0)
        ax[k].axhline(4, color=INK2, ls=":", lw=1.0)
        ax[k].text(2.2, 1.25, "1σ", color=INK2, fontsize=7.5)
        ax[k].text(2.2, 5.0, "2σ", color=INK2, fontsize=7.5)
        for wa in (0.028, 0.052, 0.08):
            ax[k].axvspan(-wa, wa, color="#e9e8e3", lw=0, zorder=0)
        ax[k].set_xlabel("w = τ_R / τ_R,0 − 1   (w = −1: Fourier film)")
        ax[k].set_title(f"({'ab'[k]}) 500 nm, noiseless profile, σ_φ = {s}°", fontsize=9, loc="left")
    ax[0].set_ylabel("Δχ²_min")
    ax[0].legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIG / "A5_profile.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    fig_dphi()
    if (OUT / "a4_fisher.json").exists():
        fig_wstar()
    fig_profile()
    print("figures written to", FIG)
