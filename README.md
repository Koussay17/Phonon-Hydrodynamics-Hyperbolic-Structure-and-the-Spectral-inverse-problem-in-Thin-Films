# Phonon Hydrodynamics, Hyperbolic Structure, and the Spectral Inverse Problem in Thin Films

**From Kinetic Closure to the Limits of Thermal Depth Profiling**

> Research in progress — CRTEn internship, September 2026–January 2027.
> Corrected baseline: 15 September 2026. No experimental validation or general novelty claim.

## Scope

When can a surface thermal measurement determine thin-film properties separately?
This project studies Fourier, Cattaneo and Guyer–Krumhansl (GK) models,
their kinetics, entropy, transport regimes and inverse problems.

## Current state

| Part | Status |
|---|---|
| I — kinetic closure | Course notes/exercises; complete closure remains unwritten. |
| II — admissibility | Near-equilibrium entropy and homogeneous propagation checks; GK infinite speed requires a nonzero nonlocal term. |
| III — AlN regimes | Paired approximate mode readings at 300 K; no full-spectrum classification or temperature trajectory. |
| IV — inverse problem | Homogeneous scale invariance, conditional sensitivities and reproduction of Fourier resonance. |
| Laboratory tools | Synthetic tools available; sample geometry, calibration, nuisance parameters and real-data validation remain required. |

## Corrected conclusions

- Homogeneous 1D scale invariance holds for free thermal parameters through
  `b, ξ₁, τ_R, τ_ℓ`. Independent thickness or constitutive constraints can restrict it.
  Graded non-Fourier layers are not implemented.
- Relaxation parameters can correlate with other fitted properties. One fixed-band
  example raises conductivity/capacity standard errors by factors 3.13 and 1.61.
- `ωτ=1` maximizes single-frequency phase sensitivity, not detectability.
  Spectral angle and half-space phase are algebraically related.
- `τ_R=τ_ℓ` gives **Fourier resonance**, established prior art, under the
  homogeneous source-free, zero-initial-perturbation boundary-value problem.
- The formal kinetic ratio `τ_R/τ_N=1.8` lacks the strong separation required
  for hydrodynamic closure and does not prove a temperature-dependent crossing.
- `3λ/(Cv)` is a **grey conductivity-derived length**, not an all-collision
  mean free path. Comparing apparent film conductivity with bulk is meaningful.

## Numerical scenarios, kept separate

For d=500 nm, C=2.41×10⁶ J·m⁻³·K⁻¹, v=6000 m/s and the illustrative sapphire substrate:

| Quantity | λ=60 W·m⁻¹·K⁻¹ | λ=321 W·m⁻¹·K⁻¹ |
|---|---|---|
| Characteristic frequency a/(2πd²) | 15.85 MHz | 84.80 MHz |
| Fourier perfect-contact reflection Γ | 0.077 | 0.460 |
| Grey transport length | 12.45 nm | 66.60 nm |
| Grey Knudsen number | 0.0249 | 0.1332 |

The default 60 is illustrative. Cheng's 321 reference comes from 18 and 22.5 µm
films, not the hypothetical 500 nm film.
Equal effusivity at λ≈44 makes the ideal Fourier interface invisible;
contacts and differing non-Fourier times require a fuller analysis.

A 200 MHz upper frequency reaches ωτ=1 at 0.796 ns. This does not prove shorter
times inaccessible: an ideal half-space with τ=10 ps has a 0.360° phase shift there.
Real detectability depends on noise, calibration and nuisance parameters.

## Repository layout

| Folder | Contents |
|---|---|
| `notes/` | Authored LaTeX/PDF notes, conventions and third-party reference PDFs. |
| `src/` | Transfer matrices, stable homogeneous impedance responses, Stehfest inversion, fitting and Bayesian sampling. |
| `notebooks/` | Four figure scripts and one notebook. |
| `tests/` | Analytical benchmarks, regressions and synthetic estimator checks. |
| `figures/` | Generated figures. |
| `theory/`, `paper/` | Placeholders; no completed manuscript is implied. |
| `data/` | README; laboratory files excluded by .gitignore. |
| `scripts/` | Figure and PDF rebuild command. |

## Install, test and rebuild

```powershell
python -m pip install -r requirements-dev.txt
python -B -m pytest tests/ -q -p no:cacheprovider
python scripts/rebuild.py
```

PDF builds require XeLaTeX, pdfLaTeX, declared LaTeX packages and DejaVu Sans.
Each note retains its existing fonts, margins and heading style; Camacho retains
its separate Latin Modern layout. Only PDFs with matching authored TeX are rebuilt.
See [conventions](notes/conventions.md) and [correction log](CORRECTIONS.md).

`fisher_analysis` raises `NonIdentifiableError` for unsupported separate uncertainties.
A fit's `success` reports optimizer convergence; inspect covariance too.
Graded variable diffusivity requires an explicit `graded_xi1`.
Bayesian priors are bounded log-uniform; seed controls the full chain, and
autocorrelation reliability is reported separately.

## Data and references

No laboratory data are included. Synthetic tests do not establish model adequacy
for a real FDTR bench. Several third-party reference PDFs already exist in `notes/`;
they are distinct from the authored notes and remain intact.

## Work before extending the conclusions

Complete the kinetic closure, justify spectral averaging and a temperature trajectory,
compare remaining prior work, then specify the actual sample, beams, interfaces,
calibration, frequencies and noise for experimental inversion.

---

## Prior work this builds on

Listed so that the boundary between what is inherited and what is new is visible from the outset.

**J.-C. Krapez, F. Rigollet**, *Comment on "Photothermal radiometry parametric identifiability theory
for reliable and unique nondestructive coating thickness and thermophysical measurements"*,
arXiv:1708.07362 (2017). Establishes that thickness, diffusivity and conductivity of a coating are
structurally correlated under front-face measurement, and identifies the transit time and the
effusivity as the two quantities the measurement determines. The homogeneous non-Fourier code reproduces the analogous invariance.

**J.-C. Krapez**, *Comment on "Simultaneous density and thermal conductivity depth profile
reconstructions…"*, J. Appl. Phys. **134**, 056101 (2023). The graded counterpart: the same
invariance, of functional dimension. Contains a section on inverse crime.

**J.-C. Krapez**, *Linear, trigonometric and hyperbolic profiles of thermal effusivity in the
Liouville space and related quadrupoles*, Int. J. Therm. Sci. **136**, 182–199 (2019; online 2018). The Liouville
transformation of the heat equation and the identification of the Schrödinger potential. Explicitly
declines the spectral reading of the resulting operator.

**A. Camacho de la Rosa, R. Esquivel-Sirvent, D. Becerril**, *Relaxation times of non-Fourier
materials using frequency-domain thermoreflectance*, J. Appl. Phys. **137**, 155103 (2025). Provides
the closed form used as external validation of the Cattaneo response, and the instrumental bandwidth
used for conditional sensitivity analysis; bandwidth alone is not a detectability threshold.

**J. Ma, W. Li, X. Luo**, *Examining the Callaway model for lattice thermal conductivity*, Phys. Rev.
B **90**, 035203 (2014). Source of the frequency scaling laws and of the relaxation times used to
place AlN on the regime map. Establishes that the Callaway model corrects the relaxation-time
approximation by only 0.7 per cent in the cross-plane direction, where 12.3 per cent is needed.

**Z. Cheng et al.**, *Experimental observation of high intrinsic thermal conductivity of AlN*, Phys.
Rev. Materials **4**, 044602 (2020). Measured conductivities, including 321 W m⁻¹K⁻¹ for MOCVD films
on sapphire, 18 and 22.5 µm thick.

**P. Chen, I. M. Gamba, Q. Li, L. Wang**, *Reconstruction of heat relaxation index in phonon transport
equation*, arXiv:2502.19533 (2025), preprint. Numerical reconstruction of the
relaxation time; leaves the sensitivity question open.

**M. S. B. Hoque et al.**, *Experimental observation of ballistic to diffusive transition in phonon
thermal transport of AlN thin films*, Appl. Phys. Lett. **125**, 262201 (2024).

**D. Maillet, S. André, J.-C. Batsale, A. Degiovanni, C. Moyne**, *Thermal Quadrupoles*, Wiley (2000).
Reference text for the transfer matrix formalism used throughout `src/`.

**R. A. Guyer, J. A. Krumhansl**, Phys. Rev. **148**, 766 and 778 (1966); **J. Callaway**, Phys. Rev.
**113**, 1046 (1959). The kinetic foundation of Part I.

---

**R. Kovács**, [Analytic solution of Guyer–Krumhansl equation for laser flash experiments](https://arxiv.org/abs/1804.05225) (2018). Fourier resonance is established prior art.

**M. G. Hennessy, T. G. Myers**, [Guyer–Krumhansl Heat Conduction in Thermoreflectance Experiments](https://doi.org/10.1007/978-3-030-64272-3_2) (2021). Complete comparison remains pending.

**G. Lebon, P. C. Dauby**, [Phys. Rev. A 42, 4710](https://doi.org/10.1103/PhysRevA.42.4710) (1990). Kinetic coefficients and assumptions.

## Licence

Code in `src/`, `tests/` and `notebooks/` is released under the MIT licence. Manuscript text and
figures are not covered by it. Laboratory data is not distributed.

---

## Contact

Koussay Mansouri — internship supervised at CRTEn.
Issues and corrections are welcome, including on the derivations.
