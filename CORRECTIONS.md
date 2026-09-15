# Corrected baseline — 15–16 September 2026

Prepared locally from commit e80fd6c442e9ff810ffb15ba1e90ae5f9136284a.
No commit or push was performed.

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

## Research that remains open

This update does not supply missing laboratory data or finish an unwritten derivation.
Remaining work includes the complete Part I closure, spectral averaging and a
temperature trajectory, experimental geometry and calibration, source supplements,
and a complete comparison with Hennessy/Myers and Krapez 2016.
The historical Camacho reading note explicitly lists unconsulted supplementary
material. These limitations are preserved rather than represented as resolved.

## Updating GitHub

Review the local diff and rebuilt PDFs, then commit and push through your usual
Git client. An update archive accompanies this local version and contains only
changed/new project files, preserving their relative paths.
The original repository was backed up separately before editing.
