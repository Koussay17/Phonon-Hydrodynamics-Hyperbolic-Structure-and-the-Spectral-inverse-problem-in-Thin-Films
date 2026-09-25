# Independent hostile peer review — bounded export pilot

**Date:** 2026-09-25. **Scope:** `00-question.md`, assumptions, pilot-run/inspection, gamma-reconstruction, pair-checks, matrix-control, harmonic controls, C-unit-source-trace, and the associated short scripts. No numerical-review report/data was read. No new numerical runs, broad literature search, or repository edits were performed.

## Verdict

The narrow claim is supported by the saved calculation: the exported pp values, frequencies, multiplicities, specified Gaussian integration and Bose factors reconstruct the **12 branch linewidths at one external grid point** to a maximum absolute difference of 6.94e-18 THz. Inspection of `reconstruct_gamma.py` confirms that it writes the accumulation directly and imports physical constants, not the package's self-energy accumulation kernel.

This is **NUMERICALLY DEMONSTRATED implementation consistency for one stipulated coarse calculation**. It does not independently validate upstream interaction normalization, material rates, continuum integration, canonical event counting, or a full conserving collision operator. The failed duplicate-frequency check remains an unresolved qualification on the export. No fatal error in the narrow reconstruction claim was identified by this bounded review.

## Major objections

### 1. Independence ends at the exported interactions and shared physical model

The reconstruction reuses the same phonons, pp, triplet weights, constants, cutoff convention, temperature and Gaussian prescription that generated the reference. It independently checks the accumulation formula and the compatibility of those exported quantities. It cannot check their upstream absolute physical normalization.

Concrete alternative explanation: multiplying an input FC3 family by a common factor a scales pp and the corresponding linewidths by a². The independent accumulation would continue to agree exactly. The same common-mode failure applies to a mistaken FC3 unit or mass/volume convention shared by the exporter and reference. The source trace provides separate normalization evidence; the numerical match cannot substitute for that trace or an independent force-constant-to-interaction benchmark.

The `incorrect_extra_division_by_Nq` failure is not an additional independent discovery: once the reconstructed vector matches, dividing it by 27 necessarily creates the reported 26/27 discrepancy. It rejects that extra division for the current convention. It is not a numerical verification of mesh scaling across grids.

**Severity:** major if described as independent physical or absolute-normalization validation. The JSON's stated limitation is appropriate and must remain in summaries.

### 2. Duplicate BZ frequencies still fail a declared check

`pilot-inspection.json` records `BZ_duplicate_frequencies=false`: 6.818098086e-8 THz exceeds the declared 1e-10 THz threshold. `inspect_pilot.py` writes the result and then asserts all checks; that invocation therefore cannot be described as an all-passing inspection.

The later controls clarify but do not close the issue:

- Without NAC, gauge-aligned dynamical matrices differ by 3.21e-16 relatively, with frequency differences around 1e-14 THz.
- With NAC, the matrix discrepancy is 1.354e-8 and the frequency difference remains 6.8181e-8 THz.
- Eigenpair residuals remain around 1e-15, so accurate diagonalization of two different matrices is consistent with the discrepancy.
- C and Rust harmonic runs share the NAC-enabled frequency discrepancy. Shared behavior does not establish backend independence or validate the shared harmonic/NAC implementation.

The evidence supports an **NAC-associated periodicity discrepancy after the tested phase alignment**. It does not yet identify the ultimate cause as truncation, source implementation, input inconsistency, or another mechanism. The earlier invalid NAC toggle is correctly retained as a failed control and cannot be counted as evidence about NAC removal.

The error is tiny relative to sigma=0.1 THz, but that ratio alone does not bound errors in eigenvectors, squared couplings, small rates, future narrow quadratures, or canonical mode identification. The gamma reconstruction inherits the same exported frequencies and therefore cannot detect this defect.

**Severity:** unresolved export-quality issue; not a falsification of the scalar reconstruction. Its propagation into pp/gamma and the unique-mode mapping remains unestablished.

### 3. Width dependence is strong; no continuum rate is established

On the same 3x3x3 mesh, changing sigma from 0.1 to 0.05 THz gives branch ratios 0.463–1.663; changing it to 0.2 gives 0.707–2.933. Thus some computed rates change by approximately -54% or +193% relative to the chosen setting. This is direct evidence of substantial quadrature sensitivity on the tested grid.

The reciprocal-pair agreement is a useful symmetry check within that same pipeline. It is not an independent physical prediction or a grid-convergence study. No uncertainty interval, converged lifetime, transport coefficient, or hydrodynamic parameter follows from these width ratios. This review did not assess the separately listed second-grid branch and makes no claim about its results.

The acoustic near-zero frequencies also change sign between harmonic controls, while the reconstruction excludes frequencies below the hard-coded 1e-4 THz cutoff. The chosen cutoff can hide those problematic modes in this scalar test. It does not validate a physical treatment of translational modes or the infrared limit.

**Severity:** fatal to converged-material-rate claims; no contradiction with the explicitly coarse pilot scope.

### 4. Scalar reconstruction does not complete event reconstruction

The eight reduced triplets and weights summing to 27 establish accounting for the selected target's scalar sum. They do not establish a canonical full-population event list. Unresolved information includes orbit expansion and mode/branch permutations, phase/gauge and degeneracy handling, physical channel orientation, repeated-index factors, reverse-process counting, and the mapping from pp to a reversible event prefactor.

The source trace explicitly says the writer call does not supply `triplet_map`. `triplet_all`, total multiplicities, and scalar agreement are insufficient by themselves to certify the signed off-diagonal products of a collision operator. Two incorrect event expansions can have the same diagonal linewidths or total weighted sum.

The saved mismatch statistic checks only f0-f1-f2 for nonzero pp. Its 100% failure at the specified exact-resonance threshold cannot be extended to every energy channel or interpreted as a physical scattering-weight distribution. Conversely, nonzero Gaussian contributions cannot be relabelled exactly resonant events. No discrete energy invariant, nonlinear Bose Jacobian, or full-population detailed-balance identity is verified by the linewidth match.

**Severity:** central incomplete deliverable if called a physical-event importer; explicitly outside the currently supported export/serialization result.

## Minor objections and reproducibility limits

- The maximum 'scaled error' divides by the largest reference linewidth; it is not the largest branchwise relative error. The full vectors are saved, which limits ambiguity, but the metric should retain its definition.
- The reconstruction hard-codes `cutoff=1e-4`; this is a version/default assumption, not an independently exported setting read by that script. The selected test does not validate all cutoff or zero-mode branches.
- `is_full_pp` means a full interaction array for the selected external points, not the full crystal collision dataset. The kappa-named output file does not turn this calculation into conductivity validation.
- The run script inventories matching HDF5 files already present in its run directory. A reproducibility claim should be tied to the saved file hashes and script/settings, not merely successful execution in a reused directory.
- I inspected stored evidence and code, but did not rehash HDF5 payloads or rerun the pilot. The reported numerical tolerances are the campaign's saved results, not a second execution by this reviewer.

## Strongest parts

Pinned input/output hashes, explicit formula and unit-source trace, a separately written accumulator, retained failed controls, and honest coarse-grid scope. The narrow exported-array reconstruction survived this review.

## Fatal issues, if any

None found for reconstruction of the selected coarse-grid linewidth vector. A claim that all export checks pass would be false. A claim of validated material rates, canonical events, or a conserving AlN operator would exceed the evidence.

## Major issues

Shared upstream normalization, unresolved NAC-associated BZ duplication, large Gaussian-width sensitivity, and missing source-to-canonical-event reconstruction.

## Minor issues

Global rather than branchwise error normalization, hard-coded cutoff, full-pp terminology, and replay/provenance limits.

## Decisive tests

For claims beyond this pilot: independently benchmark force-constant-to-pp normalization; explain and bound the duplicate-representation discrepancy; verify canonical orbit/channel expansion against a full-population action; and establish a controlled mesh/energy-integration limit. No such new work was run in this review.

## Unresolved questions

Does the duplicate-representation discrepancy materially alter exported interactions? Can symmetry-expanded canonical events reproduce the required even and odd collision actions with correct counting? What conserving integration connects broadened mesh data to the prior exact-event reference? What mesh/width/cutoff sequence yields converged rates or inverse moments? The present scalar match does not answer these questions.
