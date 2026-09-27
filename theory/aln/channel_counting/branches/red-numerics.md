# Independent numerical red-team audit — PI/C channel counting

**Status: COMPLETE bounded audit, 2026-09-27.** Reviewed the frozen candidate, note/synthesis, PI permutation and sewing scripts/results, C finite-Fock checks and C2 closure checks. No other red reports were read. D's source-normalization branch and C2's imported D comparator were excluded from approval; this report does not independently revalidate that source comparison. No material calculation or candidate edits were performed.

## Outcome and objections

No numerical counterexample to the narrowly stated PI/C claims was found in the checks below. This is a failed bounded falsification attempt, not certification of a full material kinetic operator. The surviving limits are substantive:

1. **Global residuals do not give uniform channel-relative precision.** Individual retained channels have relative differences around 1.6e-12 even when the full tensor norm residual is around 1e-14 or smaller.
2. **Sewing tests establish compatibility of the selected numerical representation.** Shared force constants and reciprocal-transform code do not independently establish physical normalization, exact degeneracy or kinetic closure. The arbitrary complex-gauge check is meaningful only with the stipulated compatible block transformations.
3. **Small finite-Fock matrix residuals do not estimate infinite-space truncation error.** Selected untruncated ladder paths can be checked exactly; operator-norm convergence of an unbounded cubic Hamiltonian does not follow.
4. **Omitted probability is not an observable-tail bound.** The actual moment and generator-expectation errors exceed omitted probability substantially. C/C2 report those errors separately; they must remain separate in any summary.

These objections constrain broader interpretations; they do not contradict the corresponding qualified frozen statements.

## 1. Independent saved-material checks

The material NPZ hash matches `material-provenance.json` and the channel record. Current PI/C script hashes match their respective result files. Direct HDF5 inspection confirms the triplet rows (1,13,28), exactly reversed addresses at (2,26,14), saved frequency slices, and all three raw sewing matrices E_q^T E_-q. Their recomputation differences are zero. Large force-constant files and the material reconstruction were not regenerated.

Using sequential tensor contractions independently of the PI's optimized einsum gives:

| Quantity | Independent result |
|---|---:|
| Largest six-permutation relative Frobenius residual | 2.387400817453143e-15 |
| Forward/inverse channel conjugation, relative Frobenius | 1.502654752085260e-14 |
| Same channel comparison, maximum absolute amplitude difference | 1.511936490468164e-17 |
| Largest permutation entry-relative difference above power floor | 1.634772135506966e-12 |
| Channel entry-relative difference above power floor | 1.576748365898744e-12 |

The entry mask is |reference amplitude|² > 1e-12 times its largest squared amplitude; it retains 1,472 entries. The global and entrywise denominators are different. No relative-accuracy claim is made for vanishing or excluded entries.

The block sewing matrices have singular values between 0.9999999999999998 and 1.0000000000000004. Their largest unitarity Frobenius residual is 7.122522901182137e-16, so there is no indicated ill-conditioning in the tested block maps. The frequency intertwining gate is an absolute THz tolerance; passing it is not a proof of exact physical degeneracy.

The first-leg detuning calculation independently reproduces 1,472 above-floor entries, minimum absolute detuning 0.0064022820850526685 THz, and zero entries below 1e-10 THz. It does not classify other orientations, excluded weak entries, or the continuum resonance surface. Replacing the declared tolerance by an arbitrary broadened admission rule would change the question.

## 2. Fresh gauge attack and negative controls

I independently changed the original and reversed bases using three fresh sets of complex unitary matrices confined to the recorded blocks. If those changes are A and B, the sewing transforms as S' = A^T S B, rather than A† S B. Amplitudes were transformed on every leg before recomputing channels.

Across these trials, channel-conjugation relative residuals are 1.5033e-14 to 1.5049e-14; forward covariance residuals are 4.84e-16 to 5.61e-16. Deliberately using A† instead of A^T produces relative errors 1.3128, 1.3934 and 1.4570. This detects a substantial arbitrary-gauge error, not merely a sign convention invisible in real eigenvectors.

The retained numerical blocks have tiny nonzero frequency differences. These algebraic gauge trials do not license arbitrary rotations between unequal-frequency scalar kinetic modes. The PI raw-contraction path explicitly removes and restores oscillator normalization for the full transformation; its simpler channel formulas rely on compatible block sewing. The frozen candidate preserves that distinction.

The sewing CLI hardcodes the selected physical row labels in its JSON. Those labels are correct for this audited NPZ; the script is a selected-case diagnostic, not a generic importer whose row metadata have been inferred from every possible input file. Similarly, provenance hashes bind these artifacts but do not replace complete environment/input validation for a new reconstruction.

## 3. Independent finite-Fock attack

I evaluated the original ordered Hamiltonian on occupation states by direct untruncated ladder-path enumeration, without sparse Kronecker matrices or a local occupation cutoff. All eight distinct-mode and four repeated-mode sample amplitudes equal C's stored sparse-matrix samples exactly in this calculation; reverse-conjugation differences are zero. Deliberately retaining duplicate repeated-index permutations gives squared-amplitude ratio 4.000000000000002, reproducing the negative control.

This independently checks that the selected paths are unaffected by the finite ladder boundary. It does not validate arbitrary states at or beyond that boundary. The sparse matrix sizes, 729 and 256, supply no infinite-Hilbert-space norm error estimate. Also, C's zero energy commutator after exact integer-energy projection is enforced by the projection construction; it is not evidence for irreversible energy-conserving dynamics. The tested incoming and reciprocal orbits are disjoint, so the experiment does not certify counting for an orbit equal to its own reciprocal orbit.

## 4. Rational tails and C2 closure

I independently obtained finite geometric moments from the memoryless-tail identity

T_k(r,M) = r^(M+1) sum_j binomial(k,j) (M+1)^(k-j) E[N^j],

using exact fractions. This uses neither C's direct partial sums nor C2's symbolic differentiation of the generating function. Every reported C thermal relative error is reproduced exactly after conversion to float. At ratio r=3/4 and cutoff 128, the annihilation-moment relative tail is 7.341367736380733e-14, although omitted probability is only 7.636651598176907e-17. This moment remains truncation-limited at that cutoff.

For C2, I expanded the jump action on the six observables independently and evaluated it using products of these closed partial moments. All three cutoffs' partial expectations and maximum errors match the saved rational values **exactly**. Maximum absolute errors are:

| Parent/daughter cutoffs | Maximum generator-expectation error |
|---|---:|
| (8,12) | 0.0030901356735375494 |
| (16,24) | 1.389412129683222e-7 |
| (32,48) | 1.159845764887653e-16 |

The final error is 8,555.5 times the omitted probability. Its value is an exact rational-tail remainder, not a floating-point solver tolerance. Destination states were not clipped, so this check is an expectation under a partial initial distribution, not a finite-state simulation with a boundary closure.

At parent/daughter means (1/3,1/2), independent polynomial moments give flux J=5/6, daughter factorial derivative 20/3, geometric tangent derivative 10/3, and non-tangency defect 10/3 = 2(2b+1)J, where b is the daughter mean. The obstruction is already a nonzero exact instantaneous derivative; no timestep or long-time numerical convergence is needed to exhibit it. This does not itself construct the full projected memory dynamics or a physical linewidth pole.

## Numerical classification and evidence

Selected tensor comparisons are globally near floating-point precision, with weaker entrywise relative precision. Sewing is well conditioned in this case. The finite ladder samples are boundary-free exact algebra checks; no infinite-operator convergence was established. Thermal moment truncations have explicitly verified tails at the tested ratios, while the C2 derivative obstruction is exact within its assumed Markov jump model. Material discretization, broadening, absolute physical conventions and full kinetic closure remain outside this audit.

Saved diagnostic and data: `experiments/R_pi_c_array_audit.py` and `experiments/R_pi_c_array_audit.json`. Run with the isolated `D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe`. The script resolves campaign files relative to itself and uses absolute material/pilot paths from `material-provenance.json`; a moved archive needs those inputs or remapped paths. It writes only its own JSON. No unfinished checks remain within this bounded scope.
