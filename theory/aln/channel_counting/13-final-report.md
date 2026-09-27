# Final report: bounded channel-counting and closure audit

## Question and outcome
How does a complex all-incoming cubic tensor map to a reversible population event, and which additional assumptions enter?
The conditional normalization is now explicit and independently audited. This completes the bounded convention investigation. It does not complete the AlN material operator, temperature trajectory, film analysis or scientific paper.

## Evidence and methods
Four isolated first passes used Hamiltonian expansion, permutation-orbit counting, sparse Fock matrices, and source/dimensional analysis. Independent literature inspection preceded controlled cross-examination. Second-generation work derived compatible sewing and tested moment closure. Four hostile reviews were completed before PI revision.
Source pp contains its Hamiltonian factorial and mesh normalization. Physical channel coefficients, Fock amplitudes, spectral rate densities, mean flux coefficients and equilibrium event weights are distinct.
Fresh full replay passed 13 top-level checks. The independent proof audit added 30 exact symbolic checks. Numerical review independently checks gauge covariance, interior ladder paths and rational tails; source normalization is outside that review's scope.

## PROVED
Within explicitly declared finite models: 6V/3V counting; frequency-compatible sewing identity; resonant positive entropy factorization; geometric-family noninvariance for the assumed Fock jump process. Definitions/domains and proofs are in branch reports; red-proof.md audits them. The simple stationary sector has exact flux zero versus geometric mean flux one.

## FORMALLY VERIFIED
None. No proof assistant was used.

## DERIVED UNDER ASSUMPTIONS
The population coefficient kappa=36 pp D_f/[hbar_ps²(1+delta_ab)] requires correctly normalized/sewn amplitudes, a justified kinetic spectral measure and statistical closure. Repeated population/source-comparator ratios (1,2) do not define a physical lifetime correction.

## NUMERICALLY DEMONSTRATED
One selected AlN tensor has permutation residual <=2.4e-15 and oriented conjugation residual about 1.50e-14 globally. Retained entrywise residuals can reach about 1.63e-12.
For 1472 entries above the stated pp floor, first-leg detuning has minimum .0064022820850526685 THz and none falls below 1e-10 THz. This is neither a whole-mesh classification nor absence of bulk scattering.
Finite Fock checks and exact-rational moment remainders agree with derivations; no timestep extrapolation is required for the instantaneous obstruction.

## LITERATURE SUPPORTED
Ordered cubic/self-energy conventions, occupation hierarchies, molecular-chaos/quasifree closure and the need to distinguish spectral and population responses. Primary references and access/coverage limits appear in L/L2. No novelty is claimed.

## SUPPORTED BUT UNPROVED
The reviewer's additional normalized correlation plateau 1/18 has its own derivation and rational experiment but was not independently proof-audited. It remains outside the accepted headline result.

## CONJECTURED
No new hypothesis is promoted to an accepted physical conclusion.

## DISPROVED
Mean occupation alone determines repeated-channel flux: explicit same-means laws give different fluxes.
The assumed geometric product family is exactly invariant away from equilibrium: exact tangency defect disproves it.
A universal identification of the full population diagonal with the source comparator: conditional repeated daughter gives ratio two.
None of these disproves a controlled continuum phonon kinetic limit.

## UNKNOWN
A valid microscopic-to-kinetic reduction for this AlN application; converged energy integration and coincident-channel weight; complete material collision action; material closure error; physical transport regime; experimental response.

## Review resolutions and limits
12-peer-review.md records every objection and its resolution or explicit non-resolution. Note 22 now fixes dimensionful wording and Planck units, adds closure literature, separates physical from algebraic validation, and records weak-channel precision and closure limits. Its preamble is unchanged from note 21; the PDF builds with zero warnings. Visual PDF inspection has not been performed.

## Reproduction and remaining work
See README.md and scripts/reproduce_channel_counting.py. Finite checks are portable. Optional source/material tests need pinned packages and retained hashed pilot inputs; raw data stay on D. Fresh local replay is established, clean external-machine replay is not.
Next: define a conserving resonance integration and kinetic regime before constructing a material operator. The manuscript still needs a focused original contribution, complete claim mapping and independent final manuscript review.

## Strongest remaining reason the conclusion could be wrong
A correct finite algebraic audit can still be attached to the wrong physical kinetic model. None of the roundoff tests certifies the material approximation, continuum limit or experimental validity.
