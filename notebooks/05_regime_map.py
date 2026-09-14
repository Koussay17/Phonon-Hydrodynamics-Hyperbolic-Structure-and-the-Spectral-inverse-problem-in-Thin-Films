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
Y_MIN, Y_MAX = 1e-1, 1e4

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
            (BLIND_X, 3.2e3), color="#b03030", fontsize=9.5,
            ha="left", va="top",
            textcoords="offset points", xytext=(8, 0))

# ---- apex -----------------------------------------------------------------
ax.plot([1.0], [1.0], "o", ms=7, mfc="white", mec="0.2", mew=1.6, zorder=5)
ax.annotate("sommet du coin", (1.0, 1.0), fontsize=9, color="0.2",
            ha="right", va="top", textcoords="offset points", xytext=(-8, -6))

# ---- region labels --------------------------------------------------------
ax.text(3e2, 3.0, "hydrodynamique", fontsize=12, color="#20591f",
        ha="center", va="center", rotation=32)
ax.text(0.28, 1.2e3, "diffusif", fontsize=12, color="#8a6a20",
        ha="center", va="center")
ax.text(3e2, 0.3, "balistique", fontsize=12, color="#2a5a80",
        ha="center", va="center")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(X_MIN, X_MAX)
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
