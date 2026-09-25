# Independent review: cutoff exports and degenerate interaction blocks

Status: COMPLETE bounded saved-array review, 2026-09-25. No new material calculation and no repository edits.

## Finding

The 10.97% globally scaled change in individual squared-interaction (`pp`) entries between cutoff factors 1.5 and 2 is strongly supported as redistribution among degenerate eigenvector labels. The evidence consists of matched eigenspaces, substantial rotations within them, and complete product-block sums agreeing near roundoff. This does not establish an exact transformation of the complex interaction amplitudes or equivalence of a full collision operator.

Scope: the saved selected-point exports on the 3 x 3 x 3 mesh, T = 300 K, Gaussian sigma = 0.1 THz, and fixed baseline Lambda. The two run metadata files record identical input-hash dictionaries; saved grid addresses, reduced triplets and weights are identical. This review independently recomputes differences from the HDF5 arrays. It is an independent reduction of the same material calculation, not an independent force-constant model or interaction implementation.

## Quantitative evidence

Let A denote factor 1.5 and B factor 2. All indices below are zero-based. The raw pp arrays have shape (8, 12, 12, 12).

| Comparison | Result |
|---|---:|
| max abs(pp_B - pp_A) / max(pp_A) | 0.1097100507683089 |
| sum abs(pp_B - pp_A) / sum(pp_A) | 0.08601756992745939 |
| max branch abs(gamma_B - gamma_A) / abs(gamma_A) | 1.384998383719554e-13 |
| max block-sum difference / largest A block | 6.287940507407432e-15 |
| L1 block-sum difference / total A block sum | 6.852431723382684e-15 |
| Maximum active degenerate-subspace projector difference, spectral norm | 3.889004043047765e-14 |

Blocks are complete Cartesian products of the frequency clusters on all three legs of each triplet. The diagnostic groups contiguous active frequencies using the mean of the two spectra and repeats the calculation at 1e-12, 1e-10 and 1e-8 THz clustering thresholds. All three give the same 11,376 blocks and the same reported errors. Modes below 1e-4 THz in both files are grouped as inactive; they are excluded from the active-subspace projector check. This grouping is not a demonstration of the physics of zero or imaginary modes.

Every one of the 4,296 raw entries changing by more than 1e-10 times the largest pp entry has at least one nonsingleton frequency cluster. None has exclusively nondegenerate legs. The largest absolute block-sum difference is 3.234276195082338e-22, between 2.123618343702346e-9 and 2.1236183437026695e-9 at triplet 7, bands (6, 8, 10).

**Relative-error qualification:** the maximum relative block error is 8.636051900578212e-10 when restricting to A blocks larger than 1e-12 times the largest A block. The global errors near 1e-14 must not be reported as relative accuracy of every small block; no uniform relative claim is made below that floor.

## Direct example at the largest raw discrepancy

The largest entry change occurs at pp[4, 7, 7, 7], with q rows (1, 13, 28). Its second leg belongs to q row 13, address (0, 0, 1), bands (6, 7). Those two frequencies are exactly equal in each saved double-precision spectrum, around 19.707804453791674 THz.

For the complete second-leg pair, with first and third bands fixed at 7:

| Cutoff factor | pp[4, 7, 6, 7] | pp[4, 7, 7, 7] | Sum |
|---|---:|---:|---:|
| 1.5 | 2.0307162358519215e-8 | 5.803041963085531e-9 | 2.6110204321604744e-8 |
| 2 | 2.5950227924155688e-8 | 1.5997639744905493e-10 | 2.6110204321604744e-8 |

For the corresponding two-column eigenvector matrices V_A and V_B, the largest off-diagonal entry of abs(V_A^* V_B)^2 is 0.18933831112406274, while norm(V_A V_A^* - V_B V_B^*, 2) is 8.096526753370513e-15. The minimum overlap singular value is 1.0000000000000002, consistent with unity at roundoff. Thus the basis rotates appreciably while the underlying subspace barely changes. Equal saved frequencies establish numerical degeneracy at this precision, not an independent proof of exact physical degeneracy.

## What the evidence does and does not establish

For a complex interaction tensor, independent unitary basis changes on three exactly degenerate legs preserve its squared Frobenius norm: the sum of abs(V_abc)^2 over a complete product of the three subspaces. The observed block sums, eigenvector rotations, stable projectors and threshold robustness provide substantial evidence for this explanation of the raw pp changes. Individual mode labels inside a degenerate eigenspace are ill-conditioned; a small matrix perturbation can select substantially different bases without changing the subspace.

The saved pp data contain squared magnitudes, not complex amplitudes. Consequently this review cannot apply the observed unitary rotations to the original amplitude tensor and verify the transformed tensor entry by entry. Block norms alone are necessary invariants, not sufficient evidence of exact tensor equivalence. The gauge explanation is strongly supported numerically, but exact coherent-amplitude equivalence remains unverified.

**Strongest remaining limitation of scalar agreement:** Bose and Gaussian factors are identical within an exactly degenerate partner block, so gamma sums can be insensitive to redistribution within that block. Block norms and scalar linewidths may agree while finer channel structure, an action on unequal populations inside the block, or off-diagonal coherences remain undetermined. A matrix of squared amplitudes cannot generally be rotated as if it were the complex amplitude tensor. Neither the linewidth agreement nor this block test certifies a complete population or coherence-resolved collision operator.

The parent comparison record reports all 18 duplicate-frequency pairs passing the unchanged 1e-10 THz threshold at both larger cutoffs, with maximum difference 2.4868995751603507e-14 THz. The original factor-1 artifact still fails, at 6.818098086114333e-8 THz; the fresh exports do not retroactively repair it. This review resolves the formerly untested pp/gamma comparison for the fresh exports only.

Numerical classification: selected invariant block sums and scalar gamma appear stabilized between factors 1.5 and 2, with globally scaled block differences near floating-point error. Raw labeled pp values remain basis sensitive. Full-operator equivalence is unresolved. Mesh refinement, Gaussian-width limits, force-constant errors and continuum/material convergence have not been tested here. No observed convergence order, rigorous cutoff tail bound or physical linewidth accuracy follows from this comparison.

## Saved evidence and reproduction

Campaign root:
`D:\ResearchLab\orchestration\campaigns\20260924-105051-aln-interaction-pilot`

- Diagnostic: `review_export_cutoff.py`.
- Persisted complete output: `review-export-cutoff.json` (131,755 bytes), including all threshold studies, rotating-cluster data, the two-entry example and reproduction metadata.
- Material export manifest: `export-cutoff-comparison.json`.
- Script SHA-256: `b74fbcbb9224d1efa9c6a2b7d1d5409f8e4e70f8905a8a98eb86f241fe8c98f2`.

Reproduce the base numerical output in PowerShell:

```powershell
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
& 'D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe' 'D:\ResearchLab\orchestration\campaigns\20260924-105051-aln-interaction-pilot\review_export_cutoff.py'
```

The script prints JSON. The saved capture additionally includes the displayed two-entry slice and provenance. It locates the comparison JSON relative to itself, but resolves HDF5 inputs and run metadata through the absolute run paths stored in that comparison JSON. Copying only the script/report/output is insufficient for an independent rerun: preserve the referenced inputs or remap those absolute paths. No broader calculation is required to reproduce this diagnostic.
