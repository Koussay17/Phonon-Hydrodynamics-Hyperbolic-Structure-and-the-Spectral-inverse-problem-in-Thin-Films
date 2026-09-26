# Independent red-team numerical audit of PI and C

**Status: COMPLETE bounded audit, 2026-09-26.** Targets: frozen `candidate-for-review.md`, PI complex-triplet scripts/JSON/NPZ, portable reconstruction source, and C lumpability source/results. No D self-review, no other red reports, no material recalculation, and no changes to the candidate implementation. Independent diagnostic: `experiments/R-numerical-audit.py`; persisted output: `experiments/R-numerical-audit.json`.

## Findings and severity

### 1. MEDIUM: frequency grouping can create overlapping blocks

Both PI reconstruction scripts choose all frequencies within 1e-10 THz of an unused seed. The `used` set excludes a seed already assigned, but does not exclude already assigned modes from a later group. Closeness to a seed is not transitive. Consequently the resulting groups need not form a partition, and overlapping polar factors need not assemble into a unitary matrix.

Reproduction in the independent diagnostic: frequencies begin with (1, 1+0.75e-10, 1+1.5e-10) THz, followed by nine separated positive values; the overlap is a unitary 3-point Fourier matrix embedded in the 12-dimensional identity. The source algorithm returns groups (0,1) and (1,2), assigning mode 1 twice.

| Quantity | Value |
|---|---:|
| Input overlap norm(U†U-I), Frobenius | 7.682120993550905e-16 |
| Assembled polar matrix norm(U_block†U_block-I), Frobenius | 0.8660254037844388 |
| Largest assembled singular value | 1.2247448713915892 |
| Smallest assembled singular value | 0.7071067811865477 |

A vector in the largest singular direction has its squared norm multiplied by 1.5, so this can corrupt the invariant the calculation is intended to test. This is an actual counterexample to the grouping routine's general reliability, not a claim that these frequencies occur in the selected AlN data.

**Effect on the frozen artifact:** no overlapping groups occur in the saved triplet. Every mode has coverage one. Thresholds 1e-12, 1e-10 and 1e-8 THz produce identical partitions and residuals. The only nonsingleton groups are second-leg bands (0,1), (2,3), (6,7), (8,9). Their maximum saved frequency span is 2.753353101070388e-14 THz in A and 1.243449787580175e-14 THz in B. The minimum gap between distinct groups on any leg is 0.07571928011952878 THz. All polar-block singular values lie between 0.9999999999999982 and 1.0000000000000013. Thus the observed artifact is not near this algorithmic failure regime. This observation neither proves physical exact degeneracy nor validates arbitrary near-degenerate clustering.

### 2. MEDIUM INTERPRETATION LIMIT: full-basis transport is an algebraic consistency check

The PI uses one reciprocal force-constant tensor for both reconstructions. For complete orthonormal eigenbases E_A and E_B, the overlap U=E_A†E_B obeys E_A U=E_B. Multilinearity then transports the contraction for *any* reciprocal tensor. Near-roundoff full-basis transport therefore does not independently validate that tensor, an all-incoming versus decay-channel convention, or absolute normalization.

The squared-amplitude comparison does check the new normal-coordinate contraction against native squared exports. It shares the reciprocal conversion implementation and its conventions, and squared values alone cannot select all phase conventions. The frozen candidate explicitly excludes canonical channel orientation, absolute physical normalization and full collision equivalence. This limitation is not a refutation of its narrower numerical statement; interpreting these residuals as independent validation of those excluded objects would be unsupported.

### 3. LOW: residual labels and precision must retain their denominators

Using sequential `tensordot` operations instead of the PI's optimized `einsum`, with frequency normalization restored before full-basis transport, gives:

| Saved-array check | Independently recomputed result |
|---|---:|
| Full complex transport, relative Frobenius | 2.113186121390093e-15 |
| Difference between independent and saved full predictions, relative Frobenius | 2.740904142404731e-16 |
| Block-only complex transport, relative Frobenius | 2.457103956627176e-14 |
| Max squared-amplitude error / max native pp, factor 1.5 | 7.726554968823502e-16 |
| Same metric, factor 2 | 6.375131771255460e-16 |

The PI's full transport value 2.147819025539619e-15 differs only at the scale expected from changing contraction order. The squared-amplitude absolute maxima are 1.985233470127266e-23 and 1.654361225106055e-23 in the stored convention.

These are global error metrics. For entries above 1e-12 times the respective largest reference entry, the largest relative squared-amplitude error reaches 6.231562438644929e-13, and the largest relative block-transport amplitude error is 2.275817566720910e-11. The latter also occurs above the stricter 1e-6 global-amplitude floor. A norm residual of 2.46e-14 does not mean every channel has that relative accuracy. No uniform relative claim is justified for arbitrarily small entries.

The field `off_block_overlap_norms` actually measures norm(U-U_block), including changes from polar normalization *inside* blocks. The pure off-block norms are (1.833010665904054e-14, 1.936190911343129e-14, 6.259549633751165e-14), slightly different from that field. This naming issue does not change the transport residual.

## Independent attack of C's propagation evidence

I reconstructed the four-mode reversible event matrices and evaluated their exponential by a quadratic polynomial in L using its two nonzero eigenvalues, rather than C's SciPy matrix exponential. Those eigenvalues are (5,9.4) for equal rates and (3.126917629116740,11.27308237088326) for unequal rates. The compressed exponential follows from G²=(47/5)G.

The independent unequal-rate closure defect has Frobenius norm 6.123724356957945, and the hidden-state block derivative is exactly represented as (-5,5,5). The second-derivative defect has Frobenius norm 13.112494041943352. At t=0.5 the relative block-propagator discrepancy is 0.04220725492119965, versus C's 0.04220725492119959. Across all four saved times, differences between the reported and independently calculated discrepancy metrics are at most 1.11e-16; reconstructed state components differ by at most 5.21e-18.

This attack did not produce a numerical counterexample to the declared finite C model or its failure of closure. Its nonzero algebraic defect is not a time-stepping or cancellation artifact. The displayed percentage is a Euclidean norm comparison in the chosen block coordinates; it is not a universal percentage error for every initial state, an entropy-weighted norm, or a material transport prediction. No such stronger conclusion was tested.

## Provenance and unresolved scope

The SHA-256 values of the frozen PI script, JSON and NPZ match `complex-triplet-provenance.json`. Both native phonon and pp files also match its recorded hashes. The saved frequencies equal the native slices, and the saved overlaps equal direct E_A†E_B recomputation. All NPZ arrays are finite; the minimum selected frequency is 3.700590101495183 THz, so this artifact does not exercise zero/imaginary-frequency normalization.

I did not independently rehash or regenerate the large FC2/FC3 inputs, or revalidate the installed reciprocal-transform implementation. The original and portable reconstruction scripts share that implementation; replaying the portable script is reproducibility evidence, not an independent reciprocal-space derivation. Both scripts are specialized to this selected 12-band, four-atom, 3 x 3 x 3 configuration. Their results establish no mesh, broadening, force-constant, infinite-cutoff or material convergence.

`second_memory_check.py` and `.json` became available while this bounded diagnostic ran; their contents were **not inspected or independently executed**. This report therefore does not certify that implementation or the sign/normalization of its memory kernel. No additional experiment is pending for this report.

**Numerical classification:** saved global contraction errors are near floating-point precision; smaller-channel relative accuracy is weaker. Selected polar blocks are well conditioned and insensitive to the three tested grouping tolerances. Generic grouping has the explicit nonunitary failure above. C's finite nonclosure is algebraically and numerically resolved. Physical exact degeneracy, canonical amplitude conventions, full material collision action and material convergence remain unresolved.

## Reproduction

Run `experiments/R-numerical-audit.py` with the isolated interpreter `D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe`. It writes only its own `R-numerical-audit.json`; no material solver is called. The script resolves candidate artifacts relative to itself and has an explicit absolute pilot-root dependency:

`D:\ResearchLab\orchestration\campaigns\20260924-105051-aln-interaction-pilot`

A moved archive needs the referenced native HDF5 inputs and a corresponding path adjustment. The diagnostic source hash is recorded in its JSON. This audit did not alter or harden the frozen candidate.
