"""Phonon transport regimes, and where the blind condition falls.

The regimes are fixed by the hierarchy of collision times. Guyer's criterion
for the hydrodynamic window reads

    tau_N  <<  tau_B  <<  tau_R ,

with tau_N the normal scattering time, tau_R the resistive one and tau_B the
time associated with boundary scattering. Two dimensionless numbers follow:

    x = tau_R / tau_N      how far resistive scattering is outrun by normal
    y = tau_B / tau_N      sample thickness in units of the normal mean free
                           path, since tau_B = d / v and tau_N = Lambda_N / v

In these variables the window becomes 1 << y << x, a wedge whose apex sits at
(1, 1). Three regions partition the plane:

    y < 1        boundary scattering outruns normal processes: ballistic
    y > x        resistive scattering outruns the boundary: diffusive
    1 < y < x    hydrodynamic

The map is material-independent. A given material traces a curve across it as
temperature varies, since tau_N and tau_R follow different laws; the function
`trajectory` converts a Callaway parameter set into that curve.

The blind condition established elsewhere, tau_R = tau_l with
tau_l = 9 tau_N / 5, reads x = 1.8 here: a vertical line passing just to the
right of the apex. At that abscissa the hydrodynamic window is less than half
a decade wide.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

BLIND_X = 1.8          # tau_R / tau_N at which the medium mimics Fourier

# --------------------------------------------------------------------------
# AlN at 300 K
# --------------------------------------------------------------------------
#
# Relaxation times read off figure 7 of Ma, Li and Luo, Phys. Rev. B 90,
# 035203 (2014), for wurtzite AlN at 300 K, over the range 2 to 10 THz.
# The readings are made by eye on a logarithmic plot and carry a factor of
# about two; they validate themselves against the published scaling laws,
# the fitted exponents of tau_N coming out at -2.00 for LA and -1.00 for TA,
# exactly the values the paper states.
#
# Sound velocity taken as 6000 m/s, the transverse value.

ALN_X_RANGE = (10.0, 200.0)            # tau_R / tau_N across the acoustic spectrum
ALN_TAU_N = (4.0e-10, 1.0e-8)          # seconds, from 10 THz down to 2 THz
ALN_VELOCITY = 6000.0                  # m/s

ALN_THICKNESSES = [
    (500e-9, "500 nm", "#b03030"),
    (5e-6, "5 µm", "#8a6a20"),
    (20e-6, "20 µm", "#1f6f8b"),
    (500e-6, "500 µm", "#20591f"),
]


def aln_y_range(thickness, velocity=ALN_VELOCITY, tau_n=ALN_TAU_N):
    """Range of y = tau_B / tau_N spanned by the acoustic spectrum."""
    tau_b = thickness / velocity
    return tau_b / tau_n[1], tau_b / tau_n[0]


def regime(x, y):
    """Name of the regime at a point of the map."""
    if y < 1.0:
        return "ballistic"
    if y > x:
        return "diffusive"
    return "hydrodynamic"


def window_width_decades(x):
    """Width of the hydrodynamic window, in decades, at a given x."""
    return np.log10(np.maximum(x, 1.0))


def trajectory(temperatures, tau_n, tau_r, thickness, velocity):
    """Curve traced by a material across the map as temperature varies.

    Parameters
    ----------
    temperatures : array-like
        Temperatures, in kelvin.
    tau_n, tau_r : callable
        Functions of temperature returning the normal and resistive times, in
        seconds. Obtained by fitting a Callaway model to conductivity data, or
        from first principles. No reliable values for AlN are used here.
    thickness : float
        Film thickness, in metres.
    velocity : float
        Mean phonon group velocity, in metres per second.

    Returns
    -------
    x, y : ndarray
        Coordinates on the map.
    """
    t = np.asarray(temperatures, dtype=float)
    tn = np.asarray([tau_n(ti) for ti in t], dtype=float)
    tr = np.asarray([tau_r(ti) for ti in t], dtype=float)
    tb = thickness / velocity
    return tr / tn, tb / tn


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

X_MIN, X_MAX = 1e-1, 1e4
Y_MIN, Y_MAX = 1e-3, 1e4

fig, ax = plt.subplots(figsize=(7.4, 6.4))

# ---- ballistic band -------------------------------------------------------
ax.fill_between([X_MIN, X_MAX], Y_MIN, 1.0, color="#c8d8e4", alpha=0.75, lw=0)

# ---- diffusive region, above the diagonal --------------------------------
xs = np.logspace(np.log10(X_MIN), np.log10(X_MAX), 400)
ax.fill_between(xs, np.maximum(xs, 1.0), Y_MAX, color="#e8dcc8", alpha=0.75, lw=0)

# ---- hydrodynamic wedge ---------------------------------------------------
wedge_x = np.logspace(0.0, np.log10(X_MAX), 200)
ax.fill_between(wedge_x, 1.0, wedge_x, color="#c9e0c9", alpha=0.85, lw=0)

# ---- boundaries -----------------------------------------------------------
ax.plot([X_MIN, X_MAX], [1.0, 1.0], color="0.3", lw=1.3)
ax.plot(xs, xs, color="0.3", lw=1.3)

# ---- the blind line -------------------------------------------------------
ax.axvline(BLIND_X, color="#b03030", lw=2.0)
ax.annotate(r"$\tau_R = \tau_\ell = 1{,}8\,\tau_N$" "\n" "réponse de Fourier exacte",
            (BLIND_X, 4e-3), color="#b03030", fontsize=9.5,
            ha="left", va="bottom",
            textcoords="offset points", xytext=(8, 0))

# ---- apex -----------------------------------------------------------------
ax.plot([1.0], [1.0], "o", ms=7, mfc="white", mec="0.2", mew=1.6, zorder=5)
ax.annotate("sommet du coin", (1.0, 1.0), fontsize=9, color="0.2",
            ha="right", va="top", textcoords="offset points", xytext=(-8, -6))

# ---- region labels --------------------------------------------------------
ax.text(1.5e3, 12.0, "hydrodynamique", fontsize=12, color="#20591f",
        ha="center", va="center", rotation=32)
ax.text(0.3, 2e3, "diffusif", fontsize=12, color="#8a6a20",
        ha="center", va="center")
ax.text(3e3, 3e-2, "balistique", fontsize=12, color="#2a5a80",
        ha="center", va="center")

# ---- AlN bands ------------------------------------------------------------
for thickness, label, colour in ALN_THICKNESSES:
    y_lo, y_hi = aln_y_range(thickness)
    ax.add_patch(Polygon(
        [(ALN_X_RANGE[0], y_lo), (ALN_X_RANGE[1], y_lo),
         (ALN_X_RANGE[1], y_hi), (ALN_X_RANGE[0], y_hi)],
        closed=True, facecolor="none", edgecolor=colour, lw=1.8, zorder=6))
    ax.annotate(f"AlN, {label}", (ALN_X_RANGE[1], y_hi), color=colour,
                fontsize=9.5, ha="left", va="center",
                textcoords="offset points", xytext=(8, 0), zorder=7)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(X_MIN, 1e5)
ax.set_ylim(Y_MIN, Y_MAX)
ax.set_xlabel(r"$x = \tau_R / \tau_N$   —   les processus normaux dominent vers la droite")
ax.set_ylabel(r"$y = \tau_B / \tau_N$   —   épaisseur en libres parcours normaux")
ax.set_title("Régimes de transport phononique, et position de la cécité",
             fontsize=12)
ax.grid(alpha=0.25, which="both")

plt.tight_layout()
plt.savefig("../figures/05_regime_map.png", dpi=180)
plt.show()

# --------------------------------------------------------------------------
# Readouts
# --------------------------------------------------------------------------

print("Carte des régimes, variables réduites")
print(f"  sommet du coin hydrodynamique : x = 1, y = 1")
print(f"  droite de cécité              : x = {BLIND_X}")
print()
print(f"  {'x = tau_R/tau_N':>18} {'largeur de la fenêtre':>24}")
for x in (1.2, BLIND_X, 5.0, 10.0, 100.0, 1000.0):
    print(f"  {x:18.1f} {window_width_decades(x):19.2f} décade(s)")
print()
print("À la droite de cécité, la fenêtre hydrodynamique mesure "
      f"{window_width_decades(BLIND_X):.2f} décade.")
print("La condition de cécité est donc franchie au moment précis où la fenêtre s'ouvre.")
print()
print("Contrôle des régions :")
for x, y in ((10.0, 0.5), (2.0, 100.0), (100.0, 10.0)):
    print(f"  x = {x:6.1f}, y = {y:6.1f}  ->  {regime(x, y)}")
print()
print("AlN à 300 K, d'après la figure 7 de Ma, Li et Luo")
print(f"  x = tau_R/tau_N sur les branches acoustiques : "
      f"{ALN_X_RANGE[0]:.0f} à {ALN_X_RANGE[1]:.0f}")
print(f"  facteur de sécurité vis-à-vis de la cécité   : "
      f"{ALN_X_RANGE[0] / BLIND_X:.0f} au pire")
print(f"  libre parcours normal, v tau_N               : "
      f"{ALN_VELOCITY * ALN_TAU_N[0] * 1e6:.1f} à "
      f"{ALN_VELOCITY * ALN_TAU_N[1] * 1e6:.0f} µm")
print()
print(f"  {'épaisseur':>10} {'y = tau_B/tau_N':>24} {'régime':>16}")
for thickness, label, _ in ALN_THICKNESSES:
    y_lo, y_hi = aln_y_range(thickness)
    if y_hi < 1.0:
        reg = "balistique"
    elif y_lo > ALN_X_RANGE[1]:
        reg = "diffusif"
    elif y_lo > 1.0:
        reg = "hydrodynamique"
    else:
        reg = "mixte"
    print(f"  {label:>10} {y_lo:9.3f} à {y_hi:9.3f} {reg:>16}")
