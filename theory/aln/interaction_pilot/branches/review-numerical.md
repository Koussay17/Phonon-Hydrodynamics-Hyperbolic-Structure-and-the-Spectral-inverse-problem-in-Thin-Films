# Independent numerical review: pilot controls

**Complete bounded review, 24 September 2026. Status: the duplicate-frequency validation remains FAILED and its cause UNRESOLVED.** No tolerance was changed, no repository was edited, and no new material calculation was run. Evidence: the three supplied control/inspection JSON records, their scripts, and a separate comparison of the four saved harmonic NPZ arrays.

## 1. What the NAC controls establish

| Requested backend | Actual NAC loaded | Maximum BZ-duplicate difference, THz | Original 1e-10 test |
|---|---:|---:|---|
| Rust | No | 3.9079850467e-14 | Pass |
| C | No | 1.8651746814e-14 | Pass |
| Rust | Yes | 6.8180980861e-8 | **Fail** |
| C | Yes | 6.8180980861e-8 | **Fail** |

The original `harmonic-controls-invalid-toggle.json` is not a valid NAC-on/off experiment. Its alleged off cases match the on cases, and the retained failure record states that explicitly supplying BORN overrode the intended off setting. The corrected script conditionally supplies BORN and records `ph.nac_params is not None`; this is materially stronger evidence of which configuration was loaded.

The corrected controls establish a **NAC-associated occurrence of the duplicate discrepancy in this setup**. They do not identify a mathematical error in the physical nonanalytic correction, a particular reciprocal cutoff, a backend defect, or an acceptable error bound. NAC changes the dynamical matrix, so a change in numerical error when it is enabled is not by itself a diagnosis.

Disabling NAC is a control, not a validation of the intended polar calculation. Averaging duplicate frequencies or enlarging the tolerance would conceal this unresolved check.

## 2. The failed pair is not a Gamma acoustic sign problem

Independent inspection of the saved arrays locates the NAC-on maximum at:

- BZ rows **23 and 24**, addresses **(-1,-1,1)** and **(-1,2,1)** on mesh (3,3,3);
- equal reciprocal-grid class, differing by one reciprocal lattice vector;
- zero-based band **0**;
- Rust frequencies **8.008505264243189** and **8.008505196062208 THz**;
- difference **6.818098086114333e-8 THz**, approximately **682 times** the unchanged acceptance threshold and approximately 8.51e-9 relative.

The C configuration gives the same difference at the same pair and band. This is a non-Gamma frequency near 8 THz, where the square-root frequency map is not singular. The tiny acoustic signs at Gamma therefore do not explain this particular failed check.

After confirming identical address arrays, the maximum Rust/C frequency difference away from Gamma is only **2.8421709430e-14 THz** with NAC off and **3.4638958368e-14 THz** with NAC on. These data make a failure unique to one requested backend unlikely. They do **not** provide independent verification of the shared NAC construction, conventions, inputs or dependencies. Backend labels and agreement are insufficient to localize the common cause.

## 3. Gamma modes: numerically sensitive, physical stability unresolved

The saved-array comparison finds no negative frequencies away from Gamma. At Gamma:

- NAC off: both configurations have three negative signed frequencies; minima are -3.3333346811e-7 and -4.8410048616e-7 THz;
- NAC on: Rust has no negative frequencies and minimum +1.6520334449e-7 THz; C has one negative frequency and minimum -4.0762002107e-7 THz.

The largest NAC-on Rust/C difference, **5.7282336556e-7 THz**, is at this sensitive Gamma sector. Near a zero dynamical-matrix eigenvalue, the signed square-root frequency has poor relative conditioning; a tiny eigenvalue perturbation can change the displayed frequency sign. These observations are **consistent with** near-zero numerical sensitivity. They do not prove that the force constants satisfy exact acoustic sum rules or that a physical instability is absent.

The controls used `symmetrize_fc=False`. No dynamical-matrix residuals or independent acoustic-nullspace residuals were saved. Consequently neither the positive Rust minimum nor the negative C minimum settles harmonic physical stability. Their signed squared frequencies, matrix scale and eigenpair residuals are the relevant next evidence; the sign alone is inadequate.

## 4. False confidence available from the inspection script

1. **Orthonormality is not an eigenpair residual.** The reported 2.8866e-15 eigenvector orthogonality error checks a unitary basis. It does not check D V - V Lambda, periodic gauge equivalence of D, or correct assignment of a matrix to a q point.
2. **The quotient uses a first representative.** `inspect_pilot.py` groups addresses modulo the mesh and selects the first row. The integer involution and 27-point coverage can pass while frequencies on different representatives disagree. A consistent single-valued physical spectrum on that quotient has not passed its own test.
3. **Reconstruction is accounting evidence.** The 1.4813e-16 scaled gamma-detail reconstruction discrepancy and exact N+U sum check related outputs of the same calculation. They do not independently verify the broadening integral, interaction prefactor, event multiplicities or physical linewidth convention. Shared errors can survive these identities.
4. **Positivity does not establish completeness or convergence.** Finite nonnegative pp and plausible positive gamma do not resolve the duplicate-frequency failure or identify the full collision operator.
5. **One channel's mismatch histogram is limited.** The script examines f0-f1-f2 among nonzero pp, without all signed channels, Bose weights or actual integration contributions. Its 100% fraction above 1e-8 THz cannot be promoted to a statement about every physical decay/absorption contribution. The script's scope warning must remain attached.
6. **No width or mesh sequence exists in these records.** One 3x3x3 grid and sigma=0.1 THz do not establish a continuum rate or a width-independent result. The controls are harmonic; they do not measure how the duplicate issue changes pp or gamma.

The inspection writes its JSON and then asserts all checks. The false duplicate flag is an actual failed validation, not an informational warning to be omitted from the pilot outcome.

## 5. One decisive bounded follow-up

**At the worst pair only**, extract the actual 12x12 dynamical matrices at the two unwrapped reciprocal coordinates, with identical input/settings, separately for NAC on and off. Retain each NAC increment (D_on-D_off). Establish the installed code's Bloch-basis phase convention, then transform the q+G matrix into the q basis using the corresponding atomic-position phase matrix; do not compare raw matrix entries before this gauge transformation.

For these four small matrices:

- measure Hermitian defects and the gauge-aligned matrix differences;
- diagonalize each using the same independent Hermitian solver and record scaled eigenpair residuals;
- compare the resulting spectra with the saved grid spectra, using the existing frequency threshold unchanged.

This distinguishes a nonperiodic/shared NAC matrix construction or truncation from a grid/cache/basis association problem and from eigensolver error. If the matrix-level mismatch is already present in the NAC increment, frequency tolerance cannot explain it away. If gauge-equivalent matrices agree but the saved spectra do not, the failure lies downstream or in the stored matrix/point association. No dense BTE matrix or additional interaction mesh is needed.

**This follow-up has not been executed in this review.** The current cause, an error bound, and the effect on the exported interaction/rate data remain unresolved.

## 6. Artifacts and classification

New audit artifacts only:

- `review_numerical_arrays.py`: read-only comparison of saved NPZ arrays;
- `review-numerical-data.json`: offending pair, negative-mode counts and aligned backend differences;
- this report.

Reproduce the small inspection with the pinned environment's Python and `review_numerical_arrays.py`; it performs no phonon solve. Source scripts inspected directly were `inspect_pilot.py`, `harmonic_controls.py` and `pilot_run.py`. No peer review findings were used. A broad filename search displayed citations from completed first-pass branches; this review's numerical conclusions come from the listed actual records and saved arrays.

**Classification:** duplicate-frequency consistency is **failed/unresolved**; Gamma signs are **conditioning-sensitive and physically inconclusive**; array/accounting checks are passed only within their stated scope; continuum mesh and broadening convergence are **unestablished**. No evidence here justifies calling the intended NAC-on export fully numerically validated, blaming only Rust, treating the invalid toggle as a control, or loosening the threshold.

## 7. Addendum, 25 September 2026: executed matrix control reviewed

**The proposed bounded matrix test has now been executed by the PI and independently inspected here.** This addendum supersedes the earlier time-specific statements that the matrix test was unexecuted and that no width sequence was available. The original failed frequency threshold remains unchanged and **failed**. No new material runs were performed for this review; only the saved four 12x12 matrices were read and diagonalized.

### Matrix evidence localizes the discrepancy before diagonalization

I inspected `matrix_control.py`, `matrix-control.json`, and both `matrix-control-nac*.npz` files. The stored coordinates differ by (0,1,0), and the atomic-position phase matrix is unitary to 3.90e-16 in Frobenius norm. Independent matrix multiplication for gauge alignment gives:

| Configuration | Gauge-aligned relative matrix discrepancy | Scaled eigenpair residuals |
|---|---:|---:|
| NAC off | 3.14e-16 | 6.25e-16, 8.22e-16 |
| NAC on | 1.3539770654e-8 | 7.06e-16, 7.28e-16 |

Both saved matrices in each configuration are exactly Hermitian in their stored entries. The original script's off discrepancy is 3.21e-16; the small change to 3.14e-16 from a different multiplication ordering is at roundoff level. The NAC-on discrepancy is unchanged at meaningful precision. Subtracting the off matrices from the on matrices gives a gauge-periodicity defect of **1.3539770766e-8**, normalized by the on-matrix norm, in the NAC increment itself.

The script records `DynamicalMatrixGL` with method `gonze` for NAC on. Its direct matrix diagonalization reproduces the frequency discrepancy as **6.8180989743e-8 THz**, consistent with the saved-grid failure to approximately 9e-15 THz. Therefore the observed duplicate-frequency discrepancy is already present in the **assembled NAC-corrected dynamical matrices**. It is not generated by the eigenvalue solver, and this reproduction does not require a grid-cache or frequency-serialization error. The opposite gauge convention gives an O(0.48) discrepancy even without NAC, whereas the selected convention restores off-matrix agreement at roundoff; choosing the opposite gauge cannot resolve the issue.

**What this does not establish:** the deeper cause within construction, finite reciprocal sums/cutoffs, parameters, or periodic-coordinate handling is still unknown. The result does not show that the physical nonanalytic correction is wrong. Its effect on pp or gamma is not bounded by these matrix norms alone, especially near degeneracies or small integration denominators. Backend agreement still cannot validate a shared matrix defect. No tolerance relaxation is justified by this localization.

### Newly available rate controls: stronger arithmetic evidence, limited scope

`gamma-reconstruction.json` reports a separate Gaussian/Bose sum over the exported arrays, without calling phono3py's accumulation routine. Its maximum absolute discrepancy is **6.94e-18 THz**, and its maximum error divided by the largest exported gamma is **2.96e-16**. This is stronger evidence for the stipulated accumulation than the earlier sum of related output fields. It shares the exported pp, frequencies, integration convention, cutoff and conversion constants; it is not an independent physical model or material validation. An extra division by 27 gives scaled error 0.962963 relative to this export, so that alteration is numerically incompatible with the recorded convention.

`pair-checks.json` now supplies widths 0.05, 0.1 and 0.2 THz and reciprocal targets 1 and 2. Its largest reported reciprocity discrepancy, scaled by the maximum rate, is **1.97e-13**. This tests the sampled scalar linewidth relation, not equivalence of the entire interaction tensor or a collision operator.

Width sensitivity is substantial: the bandwise ratios gamma(0.05)/gamma(0.1) range **0.4633-1.6629**, and gamma(0.2)/gamma(0.1) range **0.7069-2.9333**. These are observed sensitivity ratios, not a physical uncertainty interval or a continuum error bound. On one mesh they cannot separate finite-width bias from unresolved resonance quadrature. No grid refinement or width-independent material rate has been demonstrated.

### Updated conclusion

The decisive small test is complete and localizes the duplicate discrepancy to the assembled NAC contribution rather than eigensolver error. The specified periodic-frequency validation remains **failed**; the deeper construction cause and impact on exported interactions/rates remain **unresolved**. Independent accumulation and reciprocal-width controls pass their limited arithmetic checks, while actual material convergence remains **inconclusive**. This review is complete without further experiment expansion.
