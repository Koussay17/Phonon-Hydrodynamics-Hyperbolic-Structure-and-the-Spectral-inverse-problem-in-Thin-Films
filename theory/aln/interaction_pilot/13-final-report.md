# AlN interaction export pilot: research checkpoint, 25 September 2026

## Main conclusion
A real public AlN force-constant dataset now produces a reproducible small interaction export. A separately written NumPy accumulator reconstructs 12 half-linewidths at one external wavevector: maximum absolute discrepancy 6.94e-18 THz, or 2.96e-16 divided by the largest reference linewidth. This validates selected scalar accounting under the pinned conventions, not upstream microscopic normalization or a physical collision operator.

The original duplicate-frequency check FAILS: 6.8181e-8 THz exceeds 1e-10. Matrix controls localize it to the polar correction. Holding Lambda fixed while increasing its reciprocal cutoff reduces the tested gauge-aligned matrix defect from 1.35e-8 to about 5e-16. This supports a reciprocal-truncation explanation at that pair. No corrected full interaction export has been performed.

Changing Gaussian width 0.1 to 0.05 THz changes branch linewidth ratios to 0.463--1.663; changing to 0.2 gives 0.707--2.933. These are sensitivity observations, not material error bars.

## Evidence categories
- PROVED: no new theorem.
- FORMALLY VERIFIED: none.
- DERIVED UNDER ASSUMPTIONS: Gaussian channel formula and half-linewidth conversion in the pinned default conventions; independently audited.
- NUMERICALLY DEMONSTRATED: scalar reconstruction, integer orbit multiplicities, reciprocal comparison, width sensitivity, matrix and fixed-Lambda cutoff controls.
- LITERATURE SUPPORTED: established three-phonon linewidth framework and documented software conventions.
- SUPPORTED BUT UNPROVED: reciprocal truncation explains the local default-NAC discrepancy; full-pipeline propagation remains open.
- CONJECTURED: a correct conserving full-population importer may be constructed from suitably expanded material events; not implemented here.
- DISPROVED: the original export passes every declared validation check.
- UNKNOWN: canonical event normalization, complete corrected export, converged collision action, material temperature trajectory and film/experimental validity.

## Ten-question checkpoint
1. Learned: real arrays and conventions reproduce a selected scalar sum; a numerical defect has a controlled local remedy.
2. Strongest: direct reconstruction and exact integer orbit accounting for this pilot.
3. Weakened: any inference from scalar consistency to material or operator validity.
4. Dominant assumptions: shared upstream FC/phonons, default normalization/cutoff, coarse grid and Gaussian integration.
5. Unresolved: canonical counting, conserving broadening and propagation of corrected NAC.
6. Best experiment: complete corrected-NAC export, followed by mesh/width refinement.
7. Repetition: do not repeat the invalid NAC flag-only toggle or scalar checks as substitutes for operator validation.
8. Literature: source conventions and reciprocal construction clarified implementation; no new physical theory is asserted.
9. Leading hypothesis: local truncation is supported; material hydrodynamics remains unclassified.
10. Resources: move from more scalar consistency tests to corrected export and canonical-action design.

## Strongest remaining reason our conclusion could be wrong
A shared upstream convention or force-constant error can preserve scalar reconstruction and reciprocity while changing absolute rates. Even correct scalar rates do not identify off-diagonal collision action. The local cutoff remedy could leave other representatives, eigenvectors or interactions insufficiently converged.

## Scope of completion
This closes the bounded export/accounting checkpoint. It does not complete the scientific paper or the AlN research programme. Negative controls and original failures are preserved; no novelty or all-tests-passed claim is made.

## Subsequent cutoff-export checkpoint

The selected-point exports have now been recalculated at cutoff factors 1, 1.5 and 2 with fixed Lambda. Both enlarged-cutoff exports pass the unchanged inspection, including all 18 duplicate-frequency pairs (maximum 2.49e-14 THz), and independent linewidth reconstruction. The original baseline still fails.

Maximum branchwise relative linewidth change from baseline is 4.63e-7; between factors 1.5 and 2 it is 1.39e-13. Raw pp entries differ by 0.110 scaled by the maximum pp entry between the enlarged cutoffs, requiring a degeneracy-aware comparison. This is a selected-point, fixed-mesh numerical result, not a converged material operator. See export-cutoff-comparison.json and continuation report for subsequent evidence.

The subsequent independent degeneracy diagnostic supports basis redistribution: complete product-degenerate block sums agree to 6.29e-15 of their maximum. Complex phases and full operator equivalence remain unverified. See continuation-export-report.md and branches/review-export-cutoff.md.
