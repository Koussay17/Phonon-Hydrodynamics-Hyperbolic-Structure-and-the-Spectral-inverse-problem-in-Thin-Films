# Physical event validation ? final scientific synthesis
Date: 24 September 2026. Bounded validation step complete; material import open.

## Question and strongest conclusion
Can an exactly resonant reversible-event reference validate the full population
action needed before importing AlN data? Under the stated finite-mode Bose
assumptions, yes for the reference algebra and the tested numerical cases.
It cannot establish physical event selection, absolute rates or continuum
energy integration.

Four genuinely independent approaches, a controlled second generation and four
isolated hostile reviews were completed. The latter found a real admission bug:
ordered subtraction reported zero for exact stored-input thermal detunings
+1 and -1. Independent nonlinear Jacobians differed by about 44.70% and 59.33%.
The frozen implementation and reproductions remain archived. The live code uses
compensated summation, records approximate detunings, and rejects both cases.

## Evidence categories
- **PROVED:** conditional finite algebra accepted by the proof audit: the
  rank-one linearization and the existence of an odd-invisible finite event
  counting error. The physical kinetic law itself remains an assumption.
- **FORMALLY VERIFIED:** none.
- **DERIVED UNDER ASSUMPTIONS:** conservation, entropy coordinates, reciprocal
  parity structure and detuned thermal energy production.
- **NUMERICALLY DEMONSTRATED:** independent differentiation, conditioning
  failures, finite streaming example, 13-check integration comparison and
  corrected-code regressions. See validation.json for final execution totals.
- **LITERATURE SUPPORTED:** signed-event decomposition, counting conventions,
  entropy normalization and the distinction introduced by broadening.
- **SUPPORTED BUT UNPROVED:** usefulness of this reference in a future actual
  material exporter; no source-to-event normalization has been certified.
- **CONJECTURED:** a computationally affordable conserving AlN event action
  may be obtained from the available force constants; not demonstrated here.
- **DISPROVED:** the frozen tolerance contract; odd-only certification as a
  sufficient general counting test; zero equilibrium residual as sufficient
  numerical evidence.
- **UNKNOWN:** real-AlN event normalization, material selection rules,
  continuum significance of the self-reciprocal toy, hydrodynamic window,
  temperature trajectory and experimental validation.

## Scope and remaining risks
Nonzero admitted roundoff mismatch gives a PSD approximation, not the exact
Jacobian at a nonstationary Bose state. Tiny mismatch does not bound inverse
response error without additional conditioning information. Very cold-mode
entropy scales can underflow; a finite limiting matrix does not restore the
lost invertible coordinate map. The API rejects repeated indices and unequal
weights. No physical momentum map is inferred from an energy triple.
A self-reciprocal counting error can be important in a finite test while its
weight disappears in a continuum limit. No macroscopic AlN error is inferred.

## Deliverables
Corrected src/collision_events.py and tests; note 19 with the prior note's
preamble; first passes, reviews, failed candidate, scripts and data; manuscript
scope plan. Five pinned public force-constant inputs were downloaded to D,
then rehashed after reconnection. Those inputs are retained locally, not
redistributed here. Their presence is not a BTE calculation.

## Reproduction and next decision
See 07-experiments.md. Historical review scripts intentionally reproduce the
old bug; current tests exercise the corrected source. Post-fix integration
replay is a PI regression, not a new independent derivation.
The next substantive milestone is a tiny, version-pinned interaction export
with a fully traced unit/counting contract and an energy-integration strategy.
Repeating abstract matrix examples cannot supply that missing evidence.

**What is the strongest remaining reason our conclusion could be wrong?**
An actual export may use different event counting, coordinates or energy
integration from the stipulated reference. Passing these finite tests would
then fail to certify its material interpretation. Those extensions remain
explicitly unvalidated, and no novelty or experimental claim is made.

## Final execution checkpoint
233 tests passed in 189.40 s; 28 event tests and 13 post-fix comparison checks pass. Note 19 has 3 pages, zero reported warnings and an identical preamble to note 18. Visual inspection remains uncompleted because the image tool failed. All five restored input files match their original hashes.
