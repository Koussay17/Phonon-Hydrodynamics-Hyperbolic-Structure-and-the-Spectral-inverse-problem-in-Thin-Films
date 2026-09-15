# Phonon Hydrodynamics, Hyperbolic Structure, and the Spectral Inverse Problem in Thin Films

**From Kinetic Closure to the Limits of Thermal Depth Profiling**

> **Status — work in progress.** Research internship, September 2026 to January 2027, CRTEn.
> Nothing here is peer reviewed. Results, figures and claims may change without notice.

---

## Scope

This repository accompanies a study of a single question: **when can thermal transport parameters
actually be recovered from a surface measurement on a thin film, and when can they not?**

The physical setting is aluminium nitride thin films, where the phonon mean free path and the phonon
relaxation time are not negligible compared with the film thickness and the modulation period. In
that regime Fourier's law is no longer the obvious starting point, and the inverse problem changes
character.

The work has four strands:

- **Part I — kinetics.** Deriving the Guyer–Krumhansl equation from the Boltzmann transport equation
  through a moment hierarchy and a closure, using the Callaway separation between normal and umklapp
  scattering.
- **Part II — admissibility.** Which non-Fourier models are compatible with a convex entropy and with
  finite characteristic speeds, and what that imposes on the relaxation time and the nonlocal length.
- **Part III — regimes.** A transport regime map for AlN as a function of temperature and film
  thickness, anchored to published measurements rather than to estimates.
- **Part IV — the inverse problem.** After a Liouville transformation the spatial operator takes the
  form of a stationary Schrödinger operator. The question is what the position of the associated
  spectral parameter in the complex plane implies for the conditioning of parameter recovery.

A single quantity runs through all four: the **relaxation time τ**.

---

## Current state

| Part | Status |
|---|---|
| I — kinetic closure | Course notes and an exercise sheet. The derivation for the manuscript is not written. |
| II — admissibility | Note written, results tested. Entropy production and propagation speeds. |
| III — AlN regime map | Note written, figure produced, AlN placed from published relaxation times. Temperature trajectory still missing. |
| IV — spectral inverse problem | Four results established and tested. |
| Laboratory deliverable | Complete except for real data. Forward model in both regimes, inversion, uncertainties validated by Monte Carlo, identifiability analysis, experiment design map. |

### Results established in Part IV

**The scale invariance of Krapez and Rigollet survives every non-Fourier law examined.** Neither a
finite relaxation time nor a nonlocal term restores the identifiability lost under Fourier. The
reason is structural: the response depends on the layer through `ξ₁`, `b`, and two *times*, and no
length appears separately.

**The relaxation times enter as genuinely independent parameters.** Estimating them degrades neither
the effusivity nor the transit time.

**Their identifiability obeys an inverse scaling law** whose threshold sits at `ω τ = 1`, the same
value at which the energy of the Schrödinger analogy crosses the midpoint of its excursion. Three
independent routes give the same number.

**A medium with `τ_R = τ_ℓ` is thermally indistinguishable from a Fourier medium**, whatever the
common value. In microscopic terms this reads `τ_R = 1.8 τ_N`, a ratio of collision times that varies
with temperature. On the regime map that condition falls just at the apex of the hydrodynamic wedge:
a material crosses it exactly as the window opens.

### Results established in Part II

**Thermodynamic admissibility and finite propagation speed are independent criteria.** All three laws
admit a convex entropy with non-negative production; only Cattaneo propagates at finite speed. The
nonlocal term of Guyer–Krumhansl is diffusive in the flux and restores the infinite speed that
Cattaneo had removed. The hierarchy Fourier, Cattaneo, Guyer–Krumhansl is not a monotone refinement.

### Results established in Part III

**A submicron AlN film sits two decades below the hydrodynamic window** and is transitional rather
than ballistic: the total mean free path is near 67 nm against a thickness of 500 nm. A conductivity
extracted there is an apparent, thickness-dependent value, not an intrinsic property.

The normal-process mean free path, 2.4 to 60 micrometres, and the total one, 67 nm, answer different
questions and must not be confused.

### Instrumental thresholds, for the laboratory

| Quantity | Formula | Value for a 500 nm AlN film |
|---|---|---|
| Characteristic frequency | `1/(2π ξ₁²)` | 16 MHz — below it the film is invisible |
| Blind interface | `b_film = b_substrate` | 44 W m⁻¹K⁻¹ on sapphire. Measured films sit at 321, so the contrast is 0.46: **not a practical risk for this system** |
| Relaxation threshold | `ω_max τ ≥ 1` | `τ ≥ 8×10⁻¹⁰ s` for a 200 MHz bench. AlN phonon times are two to three decades below: **out of reach of FDTR** |
| Blind diagonal | `τ_R = τ_ℓ` | intrinsic to the material, cannot be worked around. AlN sits at `τ_R/τ_N = 10` to `200` at 300 K, far from it |
| Knudsen number | `3λ/(ρc v d)` | `0.13` — transitional, so the extracted conductivity is an apparent value |

The first question to settle with the laboratory is the **film thickness**: it decides what the
deliverable can claim, ahead of the measurement bandwidth.

---

## Repository layout

| Folder | Contents |
|---|---|
| `notes/` | Reading notes, one file per paper. Conventions and verified thresholds. Part I and Part IV notes. |
| `theory/` | Derivations in LaTeX, one file per step of the chain. |
| `src/` | Quadrupole assembly, numerical Laplace inversion, forward model, inversion, Bayesian sampling. |
| `notebooks/` | Figure-producing scripts, numbered in the order they were written. |
| `tests/` | Unit tests. Every routine must pass an analytical benchmark before it touches real data. |
| `data/` | **Never committed.** |
| `figures/` | Generated figures, regenerable from `notebooks/`. |
| `paper/` | Manuscript sources. |

`notes/conventions.md` is the reference for every convention, threshold and verification used in
`src/`. Any external expression must be converted to those conventions before use.

---

## Testing

```
pip install -r requirements.txt
python -m pytest tests/ -v
```

Two classes of test are kept deliberately distinct.

**Internal consistency** — unimodular determinant, homogeneous limit, composition, layer splitting.
Necessary, and insufficient: an error consistent with itself passes all of them. One did, for a
hundred and thirteen tests.

**External validation** — comparison with a closed form taken from outside the codebase. Every
constitutive law implemented carries at least one. The list is in `notes/conventions.md`, section 14.

Three errors were caught during development, and each was caught the same way: by confronting the
work with something outside it. A wrong flux coefficient, by a published closed form; an overstated
risk, by measured data; a confusion between two mean free paths, by a definitional check. None was
caught by internal consistency.

---

## Data policy

**Laboratory data is not published in this repository, and never will be.** The `data/` folder is
excluded by `.gitignore`, which was added before the first commit. Notebooks requiring real
measurements will fail on a fresh clone; that is intended. Synthetic data generators are provided in
`src/` so that every theoretical result can be reproduced without laboratory access.

Papers are not committed either. Bibliographic records live in a Zotero library; this repository
holds citations, not copies.

---

## Prior work this builds on

Listed so that the boundary between what is inherited and what is new is visible from the outset.

**J.-C. Krapez, F. Rigollet**, *Comment on "Photothermal radiometry parametric identifiability theory
for reliable and unique nondestructive coating thickness and thermophysical measurements"*,
arXiv:1708.07362 (2017). Establishes that thickness, diffusivity and conductivity of a coating are
structurally correlated under front-face measurement, and identifies the transit time and the
effusivity as the two quantities the measurement determines. **The degeneracy this work extends.**

**J.-C. Krapez**, *Comment on "Simultaneous density and thermal conductivity depth profile
reconstructions…"*, J. Appl. Phys. **134**, 056101 (2023). The graded counterpart: the same
invariance, of functional dimension. Contains a section on inverse crime.

**J.-C. Krapez**, *Linear, trigonometric and hyperbolic profiles of thermal effusivity in the
Liouville space and related quadrupoles*, Int. J. Therm. Sci. **136**, 182–199 (2018). The Liouville
transformation of the heat equation and the identification of the Schrödinger potential. Explicitly
declines the spectral reading of the resulting operator.

**A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, *Relaxation times of non-Fourier
materials using frequency-domain thermoreflectance*, J. Appl. Phys. **137**, 155103 (2025). Provides
the closed form used as external validation of the Cattaneo response, and the instrumental bandwidth
that fixes the measurability threshold.

**J. Ma, W. Li, X. Luo**, *Examining the Callaway model for lattice thermal conductivity*, Phys. Rev.
B **90**, 035203 (2014). Source of the frequency scaling laws and of the relaxation times used to
place AlN on the regime map. Establishes that the Callaway model corrects the relaxation-time
approximation by only 0.7 per cent in the cross-plane direction, where 12.3 per cent is needed.

**Y. Cheng et al.**, *Experimental observation of high intrinsic thermal conductivity of AlN*, Phys.
Rev. Materials **4**, 044602 (2020). Measured conductivities, including 321 W m⁻¹K⁻¹ for MOCVD films
on sapphire.

**P. Chen, I. M. Gamba, Q. Li, L. Wang**, *Reconstruction of heat relaxation index in phonon transport
equation*, arXiv:2502.19533 (2025), accepted in SIAM J. Appl. Math. Numerical reconstruction of the
relaxation time; leaves the sensitivity question open.

**M. S. B. Hoque et al.**, *Experimental observation of ballistic to diffusive transition in phonon
thermal transport of AlN thin films*, Appl. Phys. Lett. **125**, 262201 (2024).

**D. Maillet, S. André, J.-C. Batsale, A. Degiovanni, C. Moyne**, *Thermal Quadrupoles*, Wiley (2000).
Reference text for the transfer matrix formalism used throughout `src/`.

**R. A. Guyer, J. A. Krumhansl**, Phys. Rev. **148**, 766 and 778 (1966); **J. Callaway**, Phys. Rev.
**113**, 1046 (1959). The kinetic foundation of Part I.

---

## Licence

Code in `src/`, `tests/` and `notebooks/` is released under the MIT licence. Manuscript text and
figures are not covered by it. Laboratory data is not distributed.

---

## Contact

Koussay Mansouri — internship supervised at CRTEn.
Issues and corrections are welcome, including on the derivations.
