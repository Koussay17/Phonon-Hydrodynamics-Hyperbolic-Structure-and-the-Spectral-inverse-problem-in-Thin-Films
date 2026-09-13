"""Experiment design map: what a measurement must satisfy to be informative.

The other scripts answer "here are the data, what are the parameters". This
one answers the question that comes before: "here is a sample, what must be
measured for the answer to mean anything".

Two design variables are scanned.

    Effusivity contrast, substrate over film
        A choice of substrate, not a property of the film. The thermal wave
        only returns a signal from the interface when the two effusivities
        differ; at unit contrast the interface is strictly invisible and the
        transit time cannot be recovered at all.

    Position of the measurement band, in decades relative to the
    characteristic frequency of the film
        Far below, the response is that of the substrate. Far above, the wave
        is confined to the film and only the effusivity survives, conductivity
        and heat capacity becoming degenerate.

The output is the worst relative uncertainty on the pair (conductivity, heat
capacity), obtained from the Fisher information, with the film thickness taken
as known. Without an independent thickness the problem is not identifiable at
all, which is the result of Krapez and Rigollet (2017).

Normalisation
-------------
The uncertainty is proportional to the noise level, so the map is computed at
a reference noise and scales linearly: doubling the noise doubles every number
on it. The reference is one per cent on amplitude and a tenth of a degree on
phase, a plausible lock-in pairing. Both axes are dimensionless, so the map
holds for any film thickness and any diffusivity.
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import forward_model as fm  # noqa: E402
import inversion as inv  # noqa: E402
import quadrupoles as q  # noqa: E402

# --------------------------------------------------------------------------
# Fixed film, reference noise
# --------------------------------------------------------------------------

FILM_LAM = 60.0
FILM_RHO_C = 2.41e6
THICKNESS = 500e-9
SUB_RHO_C = 3.03e6

SIGMA_REL = 0.01          # relative, on amplitude
SIGMA_PHASE = 0.1         # degrees, on phase
N_POINTS = 60             # measurement points in the band
BAND_WIDTH = 1.0          # decades

NAMES = ["film_lam", "film_rho_c"]

B_FILM = q.effusivity(FILM_LAM, FILM_RHO_C)


def sample_for(contrast):
    """Film fixed, substrate chosen to realise the requested contrast."""
    b_sub = contrast * B_FILM
    return fm.Sample(film_lam=FILM_LAM, film_rho_c=FILM_RHO_C,
                     thickness=THICKNESS,
                     sub_lam=b_sub ** 2 / SUB_RHO_C, sub_rho_c=SUB_RHO_C)


def worst_sigma(contrast, decades):
    """Worst relative uncertainty of the pair, or NaN where unattainable."""
    s = sample_for(contrast)
    centre = np.log10(s.characteristic_frequency) + decades
    freq = np.logspace(centre - BAND_WIDTH / 2, centre + BAND_WIDTH / 2, N_POINTS)
    try:
        _, cov = inv.fisher_analysis(freq, s, NAMES, SIGMA_REL, SIGMA_PHASE)
    except (OverflowError, np.linalg.LinAlgError, ValueError):
        return np.nan
    d = np.diag(cov)
    if np.any(d <= 0.0):
        return np.nan
    return float(np.sqrt(d).max())


# --------------------------------------------------------------------------
# Scan
# --------------------------------------------------------------------------

contrasts = np.logspace(np.log10(0.3), np.log10(3.0), 61)
positions = np.linspace(-3.0, 2.0, 61)   # au-delà, la matrice de transfert déborde

grid = np.empty((contrasts.size, positions.size))
for i, c in enumerate(contrasts):
    for j, d in enumerate(positions):
        grid[i, j] = worst_sigma(c, d)

# --------------------------------------------------------------------------
# Readouts
# --------------------------------------------------------------------------

finite = np.isfinite(grid)
best_flat = np.nanargmin(np.where(finite, grid, np.inf))
bi, bj = np.unravel_index(best_flat, grid.shape)

print("Carte de conception, bruit de référence 1 % et 0,1 degré")
print(f"  points par bande      : {N_POINTS}")
print(f"  largeur de bande      : {BAND_WIDTH:.1f} décade")
print()
print(f"  meilleur point du balayage : {100 * grid[bi, bj]:.3f} %")
print(f"    contraste                : {contrasts[bi]:.3f}")
print(f"    position de bande        : {positions[bj]:+.2f} décade")
print()
best_positions = [positions[int(np.nanargmin(np.where(np.isfinite(r), r, np.inf)))]
                  for r in grid if np.any(np.isfinite(r))]
print(f"  position optimale de bande, sur tous les contrastes :")
print(f"    médiane {np.median(best_positions):+.2f} décade, "
      f"étendue {min(best_positions):+.2f} à {max(best_positions):+.2f}")
print()
print("Incertitude au meilleur placement de bande, selon le contraste :")
print(f"  {'contraste':>10} {'meilleure sigma':>17} {'position':>10}")
for c_target in (0.33, 0.5, 0.8, 0.95, 1.0, 1.05, 1.25, 2.0, 3.0):
    i = int(np.argmin(np.abs(contrasts - c_target)))
    row = grid[i]
    if not np.any(np.isfinite(row)):
        continue
    j = int(np.nanargmin(np.where(np.isfinite(row), row, np.inf)))
    print(f"  {contrasts[i]:10.3f} {100 * row[j]:16.3f} % {positions[j]:+9.2f}")

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

fig = plt.figure(figsize=(8.0, 8.6))
gs = fig.add_gridspec(3, 2, height_ratios=[1.55, 1.0, 1.0], hspace=0.42,
                      wspace=0.30)

ax = fig.add_subplot(gs[0, :])
masked = np.ma.masked_invalid(grid)
mesh = ax.pcolormesh(positions, contrasts, masked,
                     norm=LogNorm(vmin=2e-3, vmax=1e1),
                     cmap="viridis_r", shading="auto")
cs = ax.contour(positions, contrasts, masked,
                levels=[0.01, 0.05, 0.2, 1.0], colors="white",
                linewidths=0.9)
ax.clabel(cs, fmt=lambda v: f"{100 * v:g} %", fontsize=8)

ax.axhline(1.0, color="#b03030", lw=1.6)
ax.annotate("interface aveugle", (-2.85, 1.0), color="#b03030", fontsize=9,
            ha="left", va="bottom")
ax.plot([positions[bj]], [contrasts[bi]], "o", ms=7,
        mfc="none", mec="white", mew=1.8)

ax.set_yscale("log")
ax.set_xlabel("position de la bande, en décades depuis la fréquence caractéristique")
ax.set_ylabel("contraste d'effusivité")
ax.set_title("Incertitude sur le couple conductivité et capacité,\n"
             "épaisseur supposée connue", fontsize=11)
fig.colorbar(mesh, ax=ax, label=r"$\sigma$ relative, pire des deux")

# ---- cut at the optimal band position ------------------------------------
ax1 = fig.add_subplot(gs[1, :])
col = grid[:, bj]
ax1.loglog(contrasts, col, lw=2, color="#5b4a8a")
ax1.axvline(1.0, color="#b03030", lw=1.4)
ax1.annotate("interface aveugle", (1.0, col.max() * 0.7), rotation=90,
             fontsize=8.5, color="#b03030", ha="right", va="top")
ax1.set_xlabel("contraste d'effusivité")
ax1.set_ylabel(r"$\sigma$ relative")
ax1.grid(alpha=0.3, which="both")
ax1.set_title(f"Coupe à la bande optimale, {positions[bj]:+.1f} décade",
              fontsize=10)

# ---- cut at a realistic contrast ------------------------------------------
ax2 = fig.add_subplot(gs[2, :])
styles = ((0.35, "#1f6f8b", "-"), (0.95, "#b03030", "-"), (2.8, "#2d7a4f", "--"))
for c_target, colour, ls in styles:
    i = int(np.argmin(np.abs(contrasts - c_target)))
    ax2.semilogy(positions, grid[i], lw=2, color=colour, ls=ls,
                 label=f"contraste {contrasts[i]:.2f}")
ax2.axvline(0.0, color="0.35", ls="--", lw=1)
ax2.set_xlabel("position de la bande, en décades")
ax2.set_ylabel(r"$\sigma$ relative")
ax2.grid(alpha=0.3, which="both")
ax2.legend(fontsize=8.5, loc="upper left")
ax2.set_title("Coupes à contraste fixé", fontsize=10)

plt.savefig("../figures/03_experiment_design_map.png", dpi=180,
            bbox_inches="tight")
plt.show()
