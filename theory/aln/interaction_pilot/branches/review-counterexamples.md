# Counterexample review — interaction pilot

**Completed, 25 September 2026.** Independent bounded audit of the named JSON/scripts and campaign assumptions. No other review, new material run, or repository edit. Severity concerns extrapolation to physical validity.

## Observed limitations

**MEDIUM — actual failed check.** pilot-inspection.json reports BZ_duplicate_frequencies=false: equivalent grid representatives differ by 6.818098e-8 THz against a 1e-10 criterion. Its script's final all-check assertion therefore fails. With NAC, matrix-control.json retains a 1.354e-8 relative discrepancy after the correct gauge alignment; without NAC it is 3.205e-16. Eigenpair residuals remain about 7e-16. Accurate diagonalization and gauge alignment do not eliminate the tested periodicity defect. Cause and material impact remain unresolved.

**MEDIUM — observed width sensitivity.** Changing sigma from 0.1 to 0.05 THz gives rate ratios 0.463–1.663; changing it to 0.2 gives 0.707–2.933. Reciprocity passes while individual rates change almost threefold. This is not a failed reciprocity check, continuum convergence study, or physical uncertainty interval.

## Counterexamples to stronger inferences

**HIGH — single-orientation mismatch cannot classify all channels.** For (f0,f1,f2)=(1,2,3), f0-f1-f2=-4 but f0+f1-f2=0. reconstruct_gamma.py includes three Gaussian orientations. Furthermore, nonzero pp need not carry appreciable Gaussian/Bose rate weight. The reported raw-pp mismatch fraction certifies neither resonance nor its absence across channels.

**HIGH — symmetry-reduced sums do not identify physical event rows.** Eight weighted representatives covering 27 partners at one target verify selected scalar accounting. They do not establish all target modes, reciprocal orbit sizes, unordered daughter counts, branch/eigenvector symmetry maps, or event-level invariants. Different lifts can preserve weighted sums. Finite-width contributions also do not become exactly resonant events because gamma is positive: for detuning Delta, the unmodified Bose forward/reverse ratio is exp(-Delta). A conserving full-population lift was not tested.

**HIGH — common normalization errors survive.** Multiplying exported pp, gamma_detail, gamma_N, gamma_U and gamma by the same positive constant preserves the relative reconstruction, positivity, reciprocity and grid checks while changing every absolute lifetime. This is a hypothetical common upstream error, not an observed one. An inverse compensation between interaction normalization and the conversion prefactor is another algebraic ambiguity. The successful direct sum verifies its stipulated convention; rejecting an extra division by 27 does not independently certify microscopic prefactors or counting.

**Verdict:** the checks support their explicitly limited serialization/accounting scope. They do not certify a conserving physical collision operator, complete event counts, absolute microscopic normalization, or converged material rates. The failed periodicity check, observed width dependence, and hypothetical blind spots above are distinct. No additional experiment is claimed.

