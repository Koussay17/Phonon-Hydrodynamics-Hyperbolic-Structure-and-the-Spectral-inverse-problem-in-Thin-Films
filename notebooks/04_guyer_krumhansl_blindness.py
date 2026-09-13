"""Identifiability of the two relaxation times of Guyer--Krumhansl.

In the one-dimensional case the two nonlocal terms combine, and eliminating
the flux through the energy balance leaves a frequency-dependent effective
conductivity. The front-face response then depends on the film through four
quantities only, and on no length other than the transit time:

    sigma e       = xi_1 sqrt[ p (1 + tau_R p) / (1 + tau_l p) ]
    lambda_eff s  = b    sqrt[ p (1 + tau_l p) / (1 + tau_R p) ]

    xi_1 = e / sqrt(a),   b = sqrt(lambda rho c),
    tau_R = resistive relaxation time,
    tau_l = 3 ell^2 / a, the time to diffuse over the nonlocal length.

Two consequences follow.

The nonlocal length never appears alone: it enters only through tau_l, which
is a time. The scale invariance that survives Fourier and Cattaneo therefore
survives Guyer--Krumhansl as well.

The two times exchange their roles between the two combinations. When they
coincide the square roots cancel and the response is exactly that of a Fourier
medium, whatever their common value. Since tau_l = 9 tau_N / 5, this blindness
occurs at tau_R = 1.8 tau_N, a ratio of collision times that varies with
temperature.

This script maps how far from that diagonal a measurement must sit for the two
times to be separable.
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import quadrupoles as q  # noqa: E402

# --------------------------------------------------------------------------
# Film and bench
# --------------------------------------------------------------------------

FILM_LAM = 60.0
FILM_RHO_C = 2.41e6
THICKNESS = 500e-9
B_SUB = 10298.0

SIGMA_REL = 0.01          # relative, on amplitude
SIGMA_PHASE = 0.1         # degrees, on phase
N_POINTS = 80

A_FILM = FILM_LAM / FILM_RHO_C
XI1 = THICKNESS / np.sqrt(A_FILM)
B_FILM = q.effusivity(FILM_LAM, FILM_RHO_C)
F_C = 1.0 / (2.0 * np.pi * XI1 ** 2)

NAMES = ("xi1", "b", "tau_R", "tau_l")


def response(freq, xi1, b, tau_r, tau_l):
    p = 2j * np.pi * np.asarray(freq, dtype=float)
    x = xi1 * np.sqrt(p * (1.0 + tau_r * p) / (1.0 + tau_l * p))
    ls = b * np.sqrt(p * (1.0 + tau_l * p) / (1.0 + tau_r * p))
    a_ = d_ = np.cosh(x)
    b_ = np.sinh(x) / ls
    c_ = ls * np.sinh(x)
    z = 1.0 / (B_SUB * np.sqrt(p))
    return (a_ * z + b_) / (c_ * z + d_)


def observables(freq, theta):
    """Log amplitude and phase, each divided by its noise level."""
    return np.concatenate([np.log(np.abs(theta)) / SIGMA_REL,
                           np.degrees(np.angle(theta)) / SIGMA_PHASE])


def covariance(freq, params, step=1e-5):
    """Covariance in logarithmic coordinates, or None if unusable."""
    log0 = np.log(np.array(params, dtype=float))

    def model(lv):
        return observables(freq, response(freq, *np.exp(lv)))

    base = model(log0)
    if not np.all(np.isfinite(base)):
        return None

    jac = np.empty((base.size, log0.size))
    for i in range(log0.size):
        up = log0.copy()
        up[i] += step
        col = model(up)
        if not np.all(np.isfinite(col)):
            return None
        jac[:, i] = (col - base) / step

    fisher = jac.T @ jac
    try:
        cov = np.linalg.inv(fisher)
    except np.linalg.LinAlgError:
        return None
    if np.any(np.diag(cov) <= 0.0):
        return None
    return cov


def band_for(tau_r, tau_l, decades=4.0, x_max=400.0):
    """Fixed band starting at the characteristic frequency of the film.

    The upper limit is capped so that the argument of the hyperbolic functions
    stays within reach of double precision. At high frequency the two factors
    partly cancel and the argument tends to xi_1 sqrt(tau_R / tau_l) sqrt(p),
    which sets the cap.

    A band whose limits jump with the parameters introduces spurious structure
    in the scan; the band is therefore fixed except for that cap.
    """
    omega_cap = (x_max / (XI1 * np.sqrt(tau_r / tau_l))) ** 2
    f_top = min(F_C * 10.0 ** decades, omega_cap / (2.0 * np.pi))
    if f_top <= F_C * 10.0:
        return None
    return np.logspace(np.log10(F_C), np.log10(f_top), N_POINTS)


def sigma_times(tau_r, tau_l):
    """Worst relative uncertainty of the two times."""
    freq = band_for(tau_r, tau_l)
    if freq is None:
        return np.nan
    try:
        cov = covariance(freq, (XI1, B_FILM, tau_r, tau_l))
    except (OverflowError, FloatingPointError):
        return np.nan
    if cov is None:
        return np.nan
    return float(np.sqrt(np.diag(cov))[2:].max())


# --------------------------------------------------------------------------
# Scan
# --------------------------------------------------------------------------

TAU_R = 1.0e-10
ratios = np.logspace(-2, 2, 101)          # tau_l / tau_R
sigma = np.array([sigma_times(TAU_R, TAU_R * r) for r in ratios])

print("Identifiabilité des deux temps de Guyer--Krumhansl")
print(f"  film            : {THICKNESS * 1e9:.0f} nm, "
      f"fréquence caractéristique {F_C:.3e} Hz")
print(f"  tau_R fixé à    : {TAU_R:.1e} s")
print(f"  bruit           : {100 * SIGMA_REL:g} % et {SIGMA_PHASE:g} degré")
print()
print(f"  {'tau_l / tau_R':>14} {'sigma des temps':>17}")
for r_target in (0.01, 0.1, 0.5, 0.9, 0.99, 1.01, 1.1, 2.0, 10.0, 100.0):
    i = int(np.argmin(np.abs(ratios - r_target)))
    v = sigma[i]
    txt = "non calculable" if not np.isfinite(v) else f"{100 * v:15.2f} %"
    print(f"  {ratios[i]:14.3f} {txt:>17}")

finite = np.isfinite(sigma)
if np.any(finite):
    j = int(np.nanargmin(np.where(finite, sigma, np.inf)))
    print()
    print(f"  meilleur rapport : {ratios[j]:.3f}, soit {100 * sigma[j]:.3f} %")

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.0))

ax1.loglog(ratios, sigma, lw=2, color="#5b4a8a")
ax1.axvline(1.0, color="#b03030", lw=1.6)
ax1.annotate(r"$\tau_R = \tau_\ell$ : réponse de Fourier exacte",
             (1.0, np.nanmax(sigma[np.isfinite(sigma)]) * 0.6), rotation=90,
             fontsize=9, color="#b03030", ha="right", va="top")
ax1.set_xlabel(r"$\tau_\ell\,/\,\tau_R$")
ax1.set_ylabel(r"$\sigma$ relative, pire des deux temps")
ax1.grid(alpha=0.3, which="both")
ax1.set_title("Divergence à la cécité", fontsize=11)

# ---- two-dimensional map --------------------------------------------------
tr = np.logspace(-12, -8, 41)
tl = np.logspace(-12, -8, 41)
grid = np.empty((tr.size, tl.size))
for i, a_ in enumerate(tr):
    for j, b_ in enumerate(tl):
        grid[i, j] = sigma_times(a_, b_)

masked = np.ma.masked_invalid(grid)
mesh = ax2.pcolormesh(tl, tr, masked, norm=LogNorm(vmin=1e-3, vmax=1e1),
                      cmap="viridis_r", shading="auto")
ax2.plot(tl, tl, color="#b03030", lw=1.6)
ax2.annotate("diagonale aveugle", (2e-11, 1.1e-11), color="#b03030",
             fontsize=9, rotation=39)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r"$\tau_\ell$  [s]")
ax2.set_ylabel(r"$\tau_R$  [s]")
ax2.set_title("Incertitude sur les deux temps", fontsize=11)
fig.colorbar(mesh, ax=ax2, label=r"$\sigma$ relative")

plt.tight_layout()
plt.savefig("../figures/04_guyer_krumhansl_blindness.png", dpi=180)
plt.show()
