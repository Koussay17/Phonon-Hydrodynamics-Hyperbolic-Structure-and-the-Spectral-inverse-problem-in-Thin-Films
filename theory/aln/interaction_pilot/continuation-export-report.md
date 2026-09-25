# Corrected polar-cutoff export: 25 September 2026

## Question
Does the fixed-Lambda reciprocal-cutoff remedy persist through phonon and interaction export, and what changes in the selected linewidths?

## Controlled experiment
Same hashed FC2, FC3, structure and Born inputs; same software, mesh 3x3x3, external point 1, temperature 300 K and Gaussian width 0.1 THz. Fresh directories and fresh package objects for each cutoff factor 1, 1.5, 2. Lambda is fixed at 0.26615773936469767. All original data remain unchanged. An independent second execution produced identical recorded numerical comparisons.

## Results
- Default reproduces the original 6.8181e-8 THz duplicate-frequency failure.
- Both enlarged cutoffs reduce the maximum over all 18 duplicate pairs to 2.4869e-14 THz and pass the unchanged 1e-10 threshold.
- Both pass the original full inspection and independent pp-to-gamma reconstruction.
- Maximum branchwise relative gamma change from default: 4.6277e-7.
- Maximum branchwise relative gamma change between factors 1.5 and 2: 1.3850e-13.
- Individual pp values between those factors differ by 0.10971 scaled by the maximum pp entry. This is an observed array difference, not a relative uncertainty on a physical event rate.

## Scope and status
NUMERICALLY DEMONSTRATED: the stated fixed-grid selected-point checks and repeated-execution comparisons.
DERIVED UNDER ASSUMPTIONS / LITERATURE SUPPORTED: conventions inherited from the audited pilot.
PROVED / FORMALLY VERIFIED: no new theorem or formal verification.
SUPPORTED BUT UNPROVED: scalar stability alone is insufficient for full operator import.
CONJECTURED: a degeneracy-aware interaction representation will clarify the entrywise pp differences; independent analysis is recorded separately.
DISPROVED: the original default export meets the duplicate-frequency tolerance.
UNKNOWN: complete canonical collision action, independent absolute normalization, mesh/width and temperature convergence, film and experimental predictions.

## Research checkpoint
1. Learned: local NAC remedy extends to every duplicate-frequency pair on this grid and stable selected linewidths.
2. Strongest: unchanged inspections and reconstruction pass, with exactly repeated recorded results.
3. Weakened: the assumption that stable phonons and rates imply stable individual pp entries.
4. Dominant assumptions: same force constants, finite grid, finite Gaussian width, chosen external point.
5. Unresolved: interpretation and mapping of degenerate-mode interactions.
6. Best next experiment: compare interaction block sums and eigenspace projectors before constructing event rows.
7. Failed approaches: no threshold relaxation; do not compare arbitrary degenerate bases as fixed mode identities.
8. Literature: no new literature-based claim is made in this continuation.
9. Leading hypothesis: basis ambiguity is a candidate, not a conclusion from linewidth stability.
10. Resources: prioritize representation/canonical counting over another repetition of scalar accumulation.

## Strongest remaining reason our conclusion could be wrong
Shared upstream errors could survive every scalar check. A stable diagonal response does not constrain all off-diagonal collision action, and a degenerate-mode basis may obstruct a naive scalar-population event representation.

## Independent degeneracy-aware follow-up

The numerical reviewer locates all 4,296 entries changing by more than 1e-10 times the global pp maximum on at least one degenerate leg. Complete product-degenerate block sums agree to 6.29e-15 globally scaled (L1 change 6.85e-15). For the worst entry, a two-dimensional exactly degenerate subspace has substantial off-diagonal eigenvector overlap but projector difference only 8.10e-15. This supports basis redistribution, with detailed metrics and clustering tolerances in the independent report.

This strengthens the basis explanation from conjecture to numerical support for these saved arrays. It does not reconstruct complex interaction phases or establish full operator equivalence. The most informative next task is a representation-aware canonical action; treating each arbitrary degenerate eigenvector as a uniquely identified scalar event mode is unjustified.
