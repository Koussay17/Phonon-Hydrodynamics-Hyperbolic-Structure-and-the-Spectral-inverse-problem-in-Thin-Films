"""Central figure of Part IV: three related points at the same abscissa.

Everything on this figure is analytic and can be checked by hand.

Upper panel, two phases.

    arg(E) = -arctan(1 / (omega tau))
        Argument of the energy of the Schrodinger analogy, E = tau omega^2
        - i omega. It runs from -90 degrees in the diffusive limit to 0 in the
        wave limit, and passes through -45 degrees, the midpoint of its
        excursion, exactly at omega tau = 1.

    phi = -45 + (1/2) arctan(omega tau)   in degrees
        Measurable phase of the front-face response of a semi-infinite
        Cattaneo medium in the one-dimensional regime. It runs from -45 to 0
        degrees and passes through -22.5, again the midpoint of its excursion,
        exactly at omega tau = 1. Derived from Camacho de la Rosa,
        Esquivel-Sirvent and Becerril (2025), equation 15, restricted to the
        half-space.

Lower panel, the uncertainty this phase allows on the relaxation time, for a
single frequency and a phase resolution of one hundredth of a degree:

    sigma(ln tau) = sigma_phi / |d phi / d ln tau|,
    d phi / d ln tau = (1/2) omega tau / (1 + (omega tau)^2)

The sensitivity is maximal at omega tau = 1, so the uncertainty is minimal
there. It decays as omega tau below the threshold and as its inverse above:
raising the frequency indefinitely does not help.

The vertical markers show where the ceiling of a frequency-domain
thermoreflectance bench, 200 MHz, falls on this axis for two relaxation times.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Analytic curves
# --------------------------------------------------------------------------

wt = np.logspace(-4, 4, 2000)                 # omega tau

arg_E = -np.degrees(np.arctan(1.0 / wt))      # spectral analogy
phase = -45.0 + 0.5 * np.degrees(np.arctan(wt))   # measurable phase

S = 0.5 * wt / (1.0 + wt ** 2)                # d phi / d ln tau, in radians
SIGMA_PHI_DEG = 0.01                          # lock-in phase resolution
sigma = np.radians(SIGMA_PHI_DEG) / S         # relative uncertainty on tau

# Asymptotes of the sensitivity, drawn only where they apply
low = wt < 1.0
high = wt > 1.0
sigma_low = np.radians(SIGMA_PHI_DEG) * 2.0 / wt
sigma_high = np.radians(SIGMA_PHI_DEG) * 2.0 * wt

# --------------------------------------------------------------------------
# Where a 200 MHz bench reaches, for two relaxation times
# --------------------------------------------------------------------------

F_MAX = 2.0e8                                  # hertz, ceiling of a FDTR bench
OMEGA_MAX = 2.0 * np.pi * F_MAX

CASES = [
    (1.0e-8, r"exemple, $\tau = 10$ ns", "#1f6f8b"),
    (1.0e-11, r"exemple, $\tau = 10$ ps", "#b03030"),
]

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.4), sharex=True)

# ---- upper panel ---------------------------------------------------------
ax1.semilogx(wt, arg_E, lw=2, color="#1f6f8b",
             label=r"$\arg(E)$, analogie spectrale")
ax1.semilogx(wt, phase, lw=2, color="#b03030",
             label=r"$\varphi$, phase mesurable")

ax1.axvline(1.0, ls="--", lw=1, color="0.35")
ax1.axhline(-45.0, ls=":", lw=0.8, color="#1f6f8b")
ax1.axhline(-22.5, ls=":", lw=0.8, color="#b03030")

ax1.plot([1.0], [-45.0], "o", ms=7, color="#1f6f8b")
ax1.plot([1.0], [-22.5], "o", ms=7, color="#b03030")

ax1.annotate(r"$-45^\circ$", (1.0, -45.0), textcoords="offset points",
             xytext=(10, -14), fontsize=9, color="#1f6f8b")
ax1.annotate(r"$-22{,}5^\circ$", (1.0, -22.5), textcoords="offset points",
             xytext=(10, 8), fontsize=9, color="#b03030")

ax1.set_ylabel("phase [degrés]")
ax1.set_ylim(-95, 8)
ax1.grid(alpha=0.3, which="both")
ax1.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax1.set_title("Les deux phases franchissent leur médiane en "
              r"$\omega\tau = 1$", fontsize=11)

# ---- lower panel ---------------------------------------------------------
ax2.loglog(wt, sigma, lw=2, color="#5b4a8a")
ax2.loglog(wt[low], sigma_low[low], ls=":", lw=1.2, color="0.45")
ax2.loglog(wt[high], sigma_high[high], ls=":", lw=1.2, color="0.45")

ax2.axvline(1.0, ls="--", lw=1, color="0.35")
best = np.radians(SIGMA_PHI_DEG) / 0.25
ax2.plot([1.0], [best], "o", ms=7, color="#5b4a8a")
ax2.annotate(f"minimum, {100 * best:.2f} %", (1.0, best),
             textcoords="offset points", xytext=(12, -4), fontsize=9,
             color="#5b4a8a")

ax2.annotate(r"pente $-1$", (3e-3, np.radians(SIGMA_PHI_DEG) * 2.0 / 3e-3),
             textcoords="offset points", xytext=(6, 6), fontsize=9, color="0.4")
ax2.annotate(r"pente $+1$", (1e3, np.radians(SIGMA_PHI_DEG) * 2.0 * 1e3),
             textcoords="offset points", xytext=(-58, 4), fontsize=9, color="0.4")

for tau, label, colour in CASES:
    x = OMEGA_MAX * tau
    ax2.axvline(x, lw=1.4, color=colour, alpha=0.85)
    ax2.annotate(label, (x, 30.0), rotation=90, fontsize=8.5, color=colour,
                 ha="right", va="top",
                 textcoords="offset points", xytext=(-5, 0))

ax2.set_xlabel(r"$\omega\tau$")
ax2.set_ylabel(r"$\sigma(\tau)\,/\,\tau$")
ax2.set_ylim(1e-4, 1e2)
ax2.grid(alpha=0.3, which="both")
ax2.set_title("Incertitude sur le temps de relaxation, résolution de phase "
              f"de {SIGMA_PHI_DEG:g}" + r"$^\circ$", fontsize=11)

plt.tight_layout()
plt.savefig(Path(__file__).resolve().parents[1] / "figures/02_omega_tau_threshold.png", dpi=180)
plt.show()

# --------------------------------------------------------------------------
# Numerical check of the three related values
# --------------------------------------------------------------------------

print("Contrôles à omega tau = 1")
print(f"  arg(E)                 = {-np.degrees(np.arctan(1.0)):+8.3f} deg")
print(f"  phase mesurable        = {-45.0 + 0.5 * np.degrees(np.arctan(1.0)):+8.3f} deg")
print(f"  d phi / d ln tau       = {0.25:8.4f} rad = "
      f"{np.degrees(0.25):.2f} deg par facteur e")
print(f"  sigma(tau) / tau       = {100 * best:8.3f} %")
print()
print("Plafond d'un banc à 200 MHz")
for tau, label, _ in CASES:
    print(f"  {label:<26} omega_max tau = {OMEGA_MAX * tau:9.3e}")
print(f"  temps au centre à f_max     = {1.0 / OMEGA_MAX:9.3e} s")
