# Phonon Hydrodynamics, Hyperbolic Structure, and the Spectral Inverse Problem in Thin Films

**From Kinetic Closure to the Limits of Thermal Depth Profiling**

> Research in progress — CRTEn internship, September 2026–January 2027.
> Corrected baseline and research continuation: 16 September 2026. No experimental validation or general novelty claim.

## Scope

When can a surface thermal measurement determine thin-film properties separately?
This project studies Fourier, Cattaneo and Guyer–Krumhansl (GK) models,
their kinetics, entropy, transport regimes and inverse problems.

## Current state

| Part | Status |
|---|---|
| I — kinetic closure | Linear grey conserving closure derived and checked (note 14); spectral data and conserving projection framework added (note 17); operator identifiability and event-cone limits audited (note 18); quantitative AlN GK closure remains open. |
| II — admissibility | Near-equilibrium entropy and homogeneous propagation checks; GK infinite speed requires a nonzero nonlocal term. |
| III — AlN regimes | Twelve-branch RTA spectrum at 300 K and published bulk conductivity trajectory; no validated N/U hydrodynamic classification. |
| IV — inverse problem | Homogeneous scale invariance, conditional sensitivities and reproduction of Fourier resonance. |
| Laboratory tools | Axisymmetric Fourier FDTR model and synthetic nuisance-parameter study (note 16); actual sample/calibration/data required. |

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
- Under the historical GK convention, the formal kinetic ratio `τ_R/τ_N=1.8` lacks the strong separation required
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
| `notebooks/` | Four figure scripts and one notebook; the experimental figure is produced by scripts/analyse_experiment.py. |
| `tests/` | Analytical benchmarks, regressions and synthetic estimator checks. |
| `figures/` | Generated figures. |
| `theory/` | Reproducible synthetic results, source-access record and experimental input checklist. |
| `paper/` | Placeholder; no completed manuscript is implied. |
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

The grey kinetic derivation is complete in note 14. It gives a longitudinal
coefficient 4 ell²/3, compared with 3 ell² for historical GK. Both use the
existing solver through nonlocal_time=L²/a; interpreting that time requires
an explicit closure convention. The conserving leading-order conversion is
tau_l=4 tau_N/5, rather than the historical 9 tau_N/5.

Note 15 compares the accessible primary sources and checks the Hennessy–Myers
1D limit and an independently constructed Darboux profile. The full Krapez
IJHMT 2016 text and Camacho supplement remain unavailable; their detailed
comparison is explicitly incomplete.

Note 16 supplies an axisymmetric Fourier FDTR baseline and an entirely synthetic
design study, including a transducer, two interfaces, Gaussian beams and
calibration parameters. At the illustrative conductivity of 60 W/(m K),
the local uncertainty grows from 0.58% with two free thermal parameters to
7.33% with ten free parameters. These are conditional design calculations,
not measured AlN properties.

Next inputs needed: actual sample, beams, interfaces, calibration, frequencies
and noise. Note 17 adds spectral AlN sums and a published bulk temperature
trajectory; separate normal/resistive rates and a validated dynamic closure
remain missing.

### New deliverables

- [14 — kinetic derivation](notes/14_Derivation_fermeture_cinetique.pdf)
- [15 — comparison with prior work](notes/15_Comparaison_travaux_anterieurs.pdf)
- [16 — experimental preparation](notes/16_Analyse_experimentale_FDTR.pdf)
- [17 — spectral AlN, temperature and boundaries](notes/17_AlN_spectral_temperature.pdf)
- [18 — audited spectral closure and physical limits](notes/18_Audit_fermeture_spectrale.pdf)
- [Research synthesis, reviews and reproducible experiments](theory/aln/closure/13-final-report.md)
- [Spectral data, attribution and reproducibility](theory/aln/README.md)
- [Synthetic results](theory/experimental_results.json)
- [Experimental input checklist](theory/experimental_inputs.template.json)
- [Source access record](theory/source_access.md)

Reproduce the synthetic study with: python -X utf8 -B scripts/analyse_experiment.py.
It uses no laboratory data. The full rebuild command also regenerates this study and the AlN spectral figures.
Reproduce the spectral study offline with: python -X utf8 -B scripts/analyse_aln_spectrum.py.

The new AlN inputs are published calculations, not laboratory measurements.
At 300 K, the source RTA capacity and conductivities are reproduced; this
checks unit conversion and quadrature weights, not mesh convergence. The
in-plane surface calculation is not a cross-plane FDTR prediction.
The full test suite now contains 205 passing tests (17 September 2026).

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

**M. G. Hennessy, T. G. Myers**, [Guyer–Krumhansl Heat Conduction in Thermoreflectance Experiments](https://doi.org/10.1007/978-3-030-64272-3_2) (2021). Accepted manuscript compared in note 15; its 1D phase agrees after convention conversion.

**G. Lebon, P. C. Dauby**, [Phys. Rev. A 42, 4710](https://doi.org/10.1103/PhysRevA.42.4710) (1990). Kinetic coefficients and assumptions.

**L. Sendra et al.**, [Phys. Rev. B 106, 155301 (2022)](https://doi.org/10.1103/PhysRevB.106.155301). Conserving longitudinal coefficient and general dispersion.

## Licence

Code in `src/`, `tests/` and `notebooks/` is released under the MIT licence. Manuscript text and
figures are not covered by it. The converted Rao data retain CC BY 4.0;
see theory/aln/README.md for attribution. Laboratory data is not distributed.

---

## Contact

Koussay Mansouri — internship supervised at CRTEn.
Issues and corrections are welcome, including on the derivations.

## Audited spectral continuation — 17 September 2026

Four independent mathematical branches, a second-generation experiment, four independent hostile reviews, and attacks on four alternative formulations are archived under theory/aln/closure/. The final report supersedes the frozen first-pass reports where explicit errata apply.

A broad positive-matrix counterexample does **not** establish microscopic AlN nonidentifiability. In a fixed finite event cone, bounded DC response bounds the first memory moment; the bound is not uniform under mesh refinement. The supplied AlN data still lack a validated physical collision action and N/U temperature trajectory.

The bulk RTA diagnostic records provenance and evaluates pole errors without low-frequency cancellation. Its 200 MHz memory-pole errors are 6.49% basal and 14.73% along c relative to total response; these are discrete bulk-model errors, not FDTR measurement errors. Static agreement does not establish continuum moment convergence.

The inspected phono3py v4.5.0 reducible conductivity matrix is equivalent to the physical operator on odd populations; its shape alone does not certify even-sector energy/viscosity dynamics. Source and parity checks record exact scope. No large first-principles rerun or experimental validation was performed.
