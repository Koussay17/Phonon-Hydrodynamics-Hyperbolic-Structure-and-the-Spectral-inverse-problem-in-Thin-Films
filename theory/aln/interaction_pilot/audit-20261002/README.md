# Default-cutoff block comparison (audit, 2 October 2026)

Question left open by `review_export_cutoff.py` (which compared only factors 1.5 and 2.0): does the
original default-cutoff export (factor 1.0) differ from the corrected exports only by basis
redistribution inside degenerate blocks?

`review_export_cutoff_pair.py FA FB` is `review_export_cutoff.py` with the factor pair as arguments and
an added clustering tolerance of 1e-6 THz. It reads only the saved exports listed in
`../export-cutoff-comparison.json` (stored on D:). Run with the isolated research Python and h5py.

| Pair | Raw pp max change / max(pp) | Block change / max block (tol ≤ 1e-10 THz) | (tol 1e-8) | (tol 1e-6) | Blocks at 1e-6 |
|---|---|---|---|---|---|
| 1.5 vs 2.0 | 0.110 | 6.3e-15 | 6.3e-15 | 6.3e-15 | 11376 |
| 1.0 vs 1.5 | 0.343 | 0.239 | 0.217 | 3.7e-8 | 11376 |
| 1.0 vs 2.0 | 0.275 | 0.275 | 0.232 | 3.7e-8 | 11376 |

Interpretation. The first row reproduces the published result and validates the generalized script.
In the default export the polar-correction defect splits exact degeneracies by up to about 1e-7 THz,
so tolerances at or below 1e-8 THz break degenerate clusters (12528 or 11808 blocks instead of 11376)
and the apparent block changes are artefacts of that splitting. At 1e-6 THz the cluster partition is
identical to that of the corrected exports, and block sums agree to 3.7e-8 of the largest block
(L1 relative change 1.5e-8), the scale of the reciprocal-cutoff defect (matrix defect 1.35e-8).

Conclusion (NUMERICALLY DEMONSTRATED, this q-point and mesh only): the default export differs from the
corrected ones by basis redistribution within degenerate blocks plus a perturbation of order 4e-8 of
max(pp), consistent with the diagnosed polar-cutoff truncation. It does not certify material convergence.
