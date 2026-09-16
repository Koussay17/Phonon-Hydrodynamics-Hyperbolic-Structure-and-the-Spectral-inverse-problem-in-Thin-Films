# Corrected baseline — 15–16 September 2026

Prepared locally from commit e80fd6c442e9ff810ffb15ba1e90ae5f9136284a.
The corrected baseline was committed as `8ae5b50` and pushed to GitHub `main`
on 16 September 2026. The two reference PDFs added remotely in `642e6a5`
were incorporated before that push.

## What changed

All 14 authored LaTeX notes were revised and their PDFs rebuilt.
Existing document classes, fonts, font scales, page sizes, margins, heading
styles and table styles were retained. Long paths and references gained
line-break opportunities. Text corrections necessarily change pagination.
The Camacho reading note retains its distinct original Latin Modern design.
The five generated figures and notebook outputs were refreshed.
Published third-party reference PDFs were not edited.

### Scientific corrections

- Fourier resonance is attributed to prior work, including Kovács (2018);
  no claim of discovery or universal indiscernibility is made.
- Homogeneous scale invariance is conditional on free thermal parameters;
  independently imposed constitutive information can restrict it.
- The relaxation crossover is distinguished from detectability. Noise,
  calibration, parameter correlations and a fixed measurement band matter.
- The 60 and 321 conductivity scenarios are separated. The latter comes from
  18 and 22.5 micrometre films and is not assigned to a 500 nm sample.
- The conductivity-derived length is explicitly a grey estimate, not a
  full-spectrum or all-collision mean free path.
- The regime map uses paired approximate modal readings. The formal ratio
  1.8 is outside the strong separation needed for the hydrodynamic closure.
- Entropy production is quadratic near equilibrium; physical entropy is
  concave and its negative is the convex mathematical entropy.
- The Mandelis reply is summarized; successful inversion is distinguished
  from a uniqueness proof. Bibliographic dates and identified initials corrected.

### Numerical corrections

- DC transfer-matrix limits, small-argument evaluation and real-part overflow checks.
- Stable homogeneous impedance recursion instead of exponentially large matrices.
- Explicit graded Liouville thickness when variable diffusivity is specified;
  incomplete graded inputs and nonphysical sample parameters are rejected.
- SVD-based covariance with column-normalized rank detection; exact scale gauge
  raises NonIdentifiableError in design analysis and gives infinite fit covariance.
- Dimensionless fit condition number and wrapped phase residuals.
- Full emcee seed reproducibility, prior validation and autocorrelation reliability.
- Stehfest limitations describe benchmarks, not a universal usable time window.

## Validation

- 160 tests passed in the final numerical run, without warnings.
- Independent ODE integration checks homogeneous Fourier/Cattaneo/GK matrices
  and both graded Fourier forms, including the DC limit.
- Regressions cover structural degeneracy, ideal sub-crossover phase sensitivity,
  low-frequency film sensitivity, imperfect-contact contrast and RNG reproducibility.
- All five figures regenerated from source.
- All authored PDFs compiled with their original engines. Page and font checks
  and visual inspection of rendered page sheets were performed.
- No broken Python requirements were reported by pip check.

The detailed reproducibility commands are in README.md.
The build script stages PDFs and only replaces requested repository PDFs after
their builds succeed. Local logs and review renders are excluded in .build/.

## Continuation — 16 September 2026

The continuation adds notes 14–16, retaining the existing note typography,
page size, margins and heading styles. Existing notes and conventions were
updated to agree with the new derivation. Pagination can change with content.

### Derivation and literature comparisons

- Note 14 completes the conserving linear grey kinetic closure under its
  stated isotropic, constant-relaxation-time assumptions. It distinguishes
  the historical longitudinal coefficient 3 ell² from the conserving
  coefficient 4 ell²/3. The existing one-dimensional solver uses the
  longitudinal nonlocal time directly and supports either convention.
- Note 15 compares the accessible Hennessy/Myers and Beardo manuscripts,
  documents a factor-three discrepancy in the Hennessy/Myers accepted
  manuscript, and independently checks a Darboux graded-profile example.
- Source access and its limits are recorded in theory/source_access.md.
  Abstracts and indexed excerpts are not presented as full-paper reviews.

### Experimental preparation

- Note 16 and src/fdtr.py provide an axisymmetric Fourier multilayer baseline
  with Gaussian pump/probe averaging and interface resistances.
- scripts/analyse_experiment.py generates the synthetic design study,
  theory/experimental_results.json and figures/06_experimental_design.png.
  Calibration uncertainty and correlations are included explicitly.
- The two hypothetical conductivity scenarios are not measurements of a
  real sample. The synthetic multistart fit is a reproducibility check.
- theory/experimental_inputs.template.json lists missing experimental
  inputs; it is a checklist, not a runnable configuration or a dataset.

### Continuation validation

- 182 tests passed without warnings in the full numerical run.
- Added independent angular-moment, published-formula, graded-profile ODE,
  Gaussian half-space and anisotropic radial-ODE checks.
- Repeated synthetic study runs gave identical outputs; three fit starts
  converged successfully without reaching bounds.
- All new and modified notes compiled with zero reported build warnings.
  Rendered page sheets of notes 14–16 were visually inspected.

## Research that remains open

The grey derivation is complete within its assumptions. A spectral AlN
dynamic closure and the normal/resistive temperature trajectory remain
open. Note 17 now supplies observable-specific spectral averages, a 300 K
RTA spectrum, harmonic capacity and a published bulk conductivity trajectory. The full Krapez 2016 sources and Camacho supplementary material are
still missing; the literature comparison explicitly identifies these gaps.
Real experimental analysis requires sample properties, geometry, calibration
and measured data. The new Fourier baseline is not a three-dimensional GK
model and does not reproduce every boundary condition in Beardo's model.

## Updating GitHub

The corrected files are available on GitHub `main` in commit `8ae5b50`.
The local update archive contains the 62 changed/new project files from the
correction pass, preserving their relative paths; it is not a full repository clone.
The original repository was backed up separately before editing.

The continuation described above is included in this revision. The earlier
update archive covers only the initial correction pass.

## Spectral AlN continuation — 16 September 2026

- Note 17 preserves the original note layout and documents mode weights,
  two distinct lifetime averages, conserving collision projections and
  free-flight versus viscous wall times.
- Published CC BY AlN data: 12 branches, 793 irreducible points, 24^3 grid,
  300 K; archive checksum verified and SI conversion reproducible.
- Published bulk temperature table is kept separate from the spectral source.
- Steady in-plane boundary suppression is checked against an independent
  transport ODE; no cross-plane FDTR or hydrodynamic prediction is claimed.
- 200 tests passed without warnings in the full numerical run.
- Source-rate separation, collision matrix and q-grid convergence remain open.

Notes 14–16 were pushed in commit 951d254. This revision includes the
spectral AlN continuation described above.
