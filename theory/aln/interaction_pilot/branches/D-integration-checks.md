# D — Bounded numerical integration checks

**Independent first pass complete, 24 September 2026.** Design only: no material calculation or heavy run performed; no peer/PI results read. At initial inspection the target environment contained phonors but phono3py was not yet present. This report does not certify the installed 4.5.0 implementation.

## 1. Smallest useful experiment

Use one explicitly labelled coarse, full reciprocal grid and one non-Gamma target q with its distinct inversion partner. A Gamma-centred 3x3x3 grid permits such a pair; on a 2x2x2 grid every point is self-inverse modulo reciprocal lattice vectors, so it cannot test a generic two-member reversal orbit. Gamma can be a separate harmonic diagnostic. Limit cubic output to these targets and an explicit band selection; record whether pp is complete or masked/compressed.

Keep the structure, force constants, symmetry/NAC settings, cutoff, grid, backend, target modes and normalization fixed. Reuse the same interactions for integration comparisons. Save hashes, actual versions, array axes/units, chosen bands, reciprocal maps and multiplicities. Identical output after reload verifies serialization only.

## 2. Acceptance checks before interpreting a rate

| Object | Bounded check and failure interpretation |
|---|---|
| Phonons | Record negative/imaginary, zero and positive frequencies separately, by q and branch. Do not replace negatives by absolute values. A tiny Gamma acoustic residual may be numerical; classify it using dynamical-matrix eigenpair residuals/scale and stated acoustic-sum-rule treatment, not a universal frequency threshold. |
| Eigenvectors | Check normalization and orthogonality; if the dynamical matrix is accessible, check its Hermitian defect and scaled eigenpair residuals. Near degeneracy, compare subspaces or summed interactions: individual band eigenvectors and pp entries can change under a basis rotation. |
| Zero/soft modes | List all excluded modes and every cutoff. Bose occupation at exactly zero frequency is undefined. A zero pp/gamma entry may mean skipped or masked work; it is not evidence that the physical coupling vanishes. Report sensitivity to the cutoff only as sensitivity, not convergence. |
| Squared interactions | Require finite, nonnegative pp where the dataset is actually squared amplitude. Record minimum, maximum, zero/masked fractions and any negative values before clipping. Use both absolute scale and relative discrepancy; relative error at a true zero is meaningless. Squared amplitudes alone are not rates. |
| Quadrature/counting | Check periodic triplet closure using integer reciprocal-grid maps. If weights count partner-grid points for one fixed target, their sum must equal the number of those points; confirm that convention before applying this identity. Do not count Brillouin-zone boundary images twice. This does not establish the prefactor's volume/grid normalization. |
| Rate-like outputs | Distinguish squared amplitudes, integration weights, occupation factors and final damping. Individual algebraic integration terms may be signed; do not impose pp positivity on every intermediate. Record any negative final damping and its scale rather than silently clipping. |

Public documentation identifies frequency as THz in cycles, and gamma as a half-linewidth with an angular-frequency conversion needed for lifetime. These are cross-checks, not proof of this pinned export's units or counting convention. [Official output documentation](https://phonopy.github.io/phono3py/input-output-files.html) currently displays version 4.4.0; the actual 4.5.0 source/metadata must settle conversion to s^-1.

## 3. Energy-mismatch evidence to save

For every actual signed channel used by the backend, compute its mismatch, e.g. Delta_f=f0-f1-f2 or f0+f1-f2. Retain signs, units and channel labels. Evaluate the three-term sum with compensated summation and high precision on a few extreme cases; ordinary subtraction can erase a thermally significant small term between two large frequencies.

At one declared temperature, e.g. 300 K, also record theta=h*10^12*Delta_f/(k_B*T) when f is in cycles/ps. Save count-weighted and multiplicity-times-pp-weighted histograms of signed Delta_f, |Delta_f| and |theta|; report exact-zero counts, quantiles, and cumulative fractions within each proposed width. Keep invalid-frequency channels separately. Neither an unweighted histogram nor a median mismatch bounds the dominant rate contribution.

For a literal detuned Bose reaction, the forward/reverse equilibrium-factor ratio is exp(-theta), so relative-energy agreement alone does not control detailed balance. This diagnostic does not turn a broadened integration term into an exactly resonant, conserving population event.

## 4. Broadening: what a fixed tiny grid can demonstrate

If a Gaussian option is available, use three explicitly recorded widths, such as sigma/2, sigma and 2*sigma, with the actual kernel convention. Choose a starting width from the resolved mismatch distribution and label it illustrative. Compare target/band sums, active support and dominant-triplet concentration; retain relative changes and absolute differences. A tetrahedron result, if already available within the same small budget, is a different discretization check, not an exact reference. [Official integration options](https://phonopy.github.io/phono3py/command-options.html)

**Order-of-limits obstruction.** For a normalized Gaussian delta_sigma(Delta), a finite grid with every mismatch nonzero has sum -> 0 as sigma -> 0. Exactly zero mismatches instead contribute proportional to 1/sigma. Thus a narrow-width plateau, collapse or divergence on this grid is not a continuum rate limit. Tetrahedron interpolation also needs mesh refinement; absence of sigma does not remove discretization error.

A subsequent convergence study would first refine the grid at fixed physical width, then reduce width while resolving the resonance-normal energy spacing. For regular resonance surfaces the desired hierarchy is local mesh-induced energy spacing << sigma << variation scale of the integrand. Degenerate surfaces or van Hove structure require separate treatment. Two tiny grids alone cannot establish an asymptotic order or justify Richardson extrapolation. This pilot does not perform that study.

## 5. Permitted conclusion and explicit status

**Can certify, if checks pass:** reproducible execution/serialization, declared array/map integrity, sampled positivity and harmonic residuals, selected counting identities, and sensitivity of specified coarse sums to the chosen integration settings.

**Cannot certify:** a converged lifetime, finite continuum memory, complete physical event operator, unique energy-only nullspace, broadening-independent rates, experimental agreement, or equivalence of two backends without matched conventions and independent numerical evidence. Raw pp and gamma_detail are insufficient for those claims.

**Numerical status now:** material results and source-to-rate conversion are **inconclusive/not yet assessed**. The designed fixed-grid resonance integration is potentially **under-resolved**; precision, cutoff and width sensitivities are separate from mesh convergence. No temporal integrator, PDE boundary, CFL limit, iterative solver tolerance or spatial propagation scheme is involved. A successful coarse pilot remains an export and consistency test.
