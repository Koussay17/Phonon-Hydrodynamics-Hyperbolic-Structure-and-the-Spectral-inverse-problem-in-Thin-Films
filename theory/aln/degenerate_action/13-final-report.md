# Degenerate collision action: audited checkpoint
Completed 26 September 2026. This report supersedes broad readings of the frozen first-pass candidate; historical reports remain unchanged.

## Question and outcome
Can squared interactions or invariant block sums supply an autonomous, basis-independent collision action? **Not in general.** Explicit finite counterexamples separate three questions: phase information, invariant observables, and dynamical closure. A supplied fixed-basis incoherent population model may still be built from squared interactions plus correct counting. No failure of that approximation for real AlN is established here.

The practical material step is positive but narrower: complex amplitudes were recovered for one selected AlN triplet from the retained force constants, and their squared magnitudes reproduce the previous native export. This is contraction validation, not a completed collision importer.

## Mathematical content
For finite dot x=-Lx, block totals N=Bx close for every initial perturbation iff BL=GB. In the stipulated reversible event model L=MW^-1, W is positive Bose variance and
Sigma=BWB^T, H=WB^T Sigma^-1, G=BLH.
The closure residual is BL(I-HB). G alone is an instantaneous compression. Within an exactly degenerate block, H divides total population by dimension, and Sigma contains dimension times per-mode variance.

C's exactly resonant four-mode model compares positive prefactors (2,2) and (3,1). The total coupling strength and compressed G agree, but only the equal-rate case closes. Equal initial occupations do not repair the other case: hidden imbalance is generated. These are synthetic scalar events, not AlN collisions.

D's matrices [[1,1],[1,1]] and [[1,1],[1,-1]] have identical squared entries and Gram spectra{0,4} versus{2,2}. This proves non-injectivity in general, not a universal claim for all tensor dimensions. Diagonal population data cannot undergo arbitrary unitary changes without missing coherences.

For a supplied symmetric positive entropy generator partitioned as [[A,K],[K*,D]], exact elimination gives
a'=-Aa+integral K exp[-D(t-s)]K* a(s) ds-K exp(-Dt)h(0).
Its Schur operator is zI+A-K(zI+D)^-1K*, Re z>0. This is known projection theory. Positivity decouples hidden zero modes, but conserved full-system modes can still prevent a DC inverse. The11-check second-generation experiment includes80-digit arithmetic and shows conserved-pole normalization can hide a 25% relaxing-channel compression error.

## Material evidence and discovered failure
The one-triplet amplitude calculation uses the same reciprocal force tensor and an independent normal-coordinate contraction. Globally scaled squared-interaction errors are below 8e-16. Complete-basis transport gives 2.15e-15; restricting to numerical degenerate blocks gives 2.46e-14 relative Frobenius error. Channelwise relative errors can be 2.28e-11. Neither global metric proves accuracy of all weak channels.

The numerical reviewer found a real generalization bug: nearby-frequency chains could produce overlapping groups and a nonunitary transform. It does not affect this saved triplet, whose groups are disjoint over three tested tolerances. The portable script now rejects ambiguous chains, compares partitions and checks unitarity. Eight regression checks and a fresh amplitude replay preserve the reported results. The broken algorithm and counterexample are retained.

Physical exact degeneracy, decay/absorption orientation, reciprocal conjugation and repeated-leg prefactors remain unvalidated. Complex all-incoming amplitudes are not yet canonical events or an irreversible generator.

## Evidence categories
PROVED: scoped finite counterexamples, closure criterion and exact elimination identities, examined by independent proof audit.
FORMALLY VERIFIED: none.
DERIVED UNDER ASSUMPTIONS: reversible-event compression and separately stated matrix/bath examples.
NUMERICALLY DEMONSTRATED: explicit toy checks, repository-action comparison, selected material contraction and conservative clustering replay.
LITERATURE SUPPORTED: exact lumping, phase-sensitive matrix kinetics and projection memory.
SUPPORTED BUT UNPROVED: no additional material theorem; a richer representation is a justified direction, not a validated AlN model.
CONJECTURED: which reduced variables will suffice for the intended material observable remains a research hypothesis.
DISPROVED: general squared-data sufficiency; invariant-sum-implies-closure in the stated class; general correctness of the old clustering routine.
UNKNOWN: complete physical collision action, material symmetry/secular regime, transport/temperature convergence and experimental predictions.

## Ten-question checkpoint
1. Learned: norms, phases, observables and closed dynamics are different data requirements.
2. Strongest: audited finite obstruction and closure proofs; bounded material contraction.
3. Weakened: direct promotion of averaged interactions to autonomous block kinetics.
4. Dominant assumptions: chosen kinetic approximation, exact versus clustered degeneracy, channel counting and retained state space.
5. Disagreements: sign conventions and population versus matrix models resolved by explicit domains; material model choice remains open.
6. Best next experiment: orient and normalize one physical channel, then test its conservation/covariance before expanding.
7. Repeated failures: do not re-test scalar linewidth agreement as a substitute for operator validation.
8. Literature: all general principles are established; no novelty claimed.
9. Leading hypothesis: phase-preserving extraction is useful, but does not choose the kinetic law.
10. Resources: prioritize channel conventions and model definition; postpone large grids and experimental claims.

## Strongest remaining reason our conclusion could be wrong
The finite proofs need not control the physical approximation relevant to AlN: crystal symmetry, complete channel sets or dephasing may restore a useful closure absent in toy examples. The material check also shares upstream conventions and covers one triplet only. These limits prevent extrapolation to a real collision operator.

## Completion boundary
The bounded sufficiency question is resolved by explicit counterexamples and audited criteria; the one-triplet phase check is complete. The full research programme and scientific paper are not complete. Next work is the importer contract in branches/second-importer-contract.md, starting with explicit channel orientation/counting and a declared kinetic approximation.
