# Second-generation frozen-operator import check

**Complete.** Computation executed on 2026-09-23; final report saved after the explicit resume on 2026-09-24. The script and results survived the interruption and were verified before this report was finalized. No further computation or scope expansion was needed.

Only permitted A/B/D/L conclusions, the independent C implementation, and the frozen source/manifest were read. No red report or red-team discussion was read. No repository implementation, frozen snapshot, or first-pass artifact was changed.

## 1. Outcome

The frozen DecayEventOperator agrees with the independently implemented seven-mode C model on full, even, and odd collision actions, nonlinear population updates, and finite-wavevector resolvents in the tested regime.

The odd-sector limitation survives import: duplicating the self-reciprocal event changes the full/even action while leaving the entire odd collision block unchanged. The constructor sums the supplied rows and does not canonicalize or deduplicate unordered events. Correct counting remains the caller's responsibility under its stated once-per-event contract.

A clearly detuned reciprocal pair is rejected. Genuine repeated indices are explicitly unsupported and rejected. **Material export remains unvalidated.**

## 2. Source provenance and execution

Tested snapshot: review-inputs/collision_events.py.

SHA-256, matching review-inputs/manifest.json:
~~~
8daa09d8f851da3b7823bdd5cc1ad86f98ced5892c814a380dcb58ce50f020c6
~~~

Independent source: experiments/C_parity_events.py, SHA-256:
~~~
de62b4b2906f8d74ed53b48999f5a598a8cf31f091dbd0a524d9b3e0335092b6
~~~

Both hashes were checked before and after the executed comparison and were unchanged. The JSON also records the new script's hash. Bytecode writing was disabled.

Owned artifacts:

- experiments/second_import_check.py
- experiments/second_import_results.json
- branches/second-generation-import-check.md

Executed command:
~~~powershell
py -B 'C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260917-232509-aln-physical-events\experiments\second_import_check.py'
~~~

Environment: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1. Exit code 0; **all 13 recorded checks passed**. The resumed verification confirmed both experiment files exist and all saved check flags are true.

## 3. Adapter and assumptions

C stores (daughter, daughter, parent). DecayEventOperator requires (parent, daughter, daughter):

| Independent tuple | Adapted tuple |
|---|---|
| (0,1,2) | (2,0,1) |
| (3,4,5) | (5,3,4) |
| (0,3,6) | (6,0,3) |

The final event is self-reciprocal, but its daughters 0 and 3 are distinct full mode indices; it is inside the supported contract.

Dimensionless energies are (1,2,3,1,2,3,2), with equal mode weights and prefactors (1,1,1). These prefactors are the bare Gamma of the stipulated nonlinear event bracket in an arbitrary test time scale. The implementation supplies the equilibrium Bose factor internally. No linewidth, lifetime, or material matrix element was used.

The independent nonlinear implementation uses absorption-minus-decay with daughter-positive stoichiometry. The frozen implementation uses decay-minus-absorption with the corresponding population updates. Reordering and sign accounting give the same nonlinear action; the comparison tests this explicitly rather than relying only on sign-invariant outer products.

## 4. Agreements and differences with A/B/D/L

| Topic | Reconciliation with C |
|---|---|
| Entropy coordinates | A/B/D/L and C use delta n divided by sqrt(n0(1+n0)); the resonant event is Gamma times one equilibrium Bose factor times a stoichiometric outer product. |
| Reaction counting | Forward and reverse reactions are already combined in the net bracket; an extra reverse-event row would duplicate the contribution. |
| Reciprocal counting | A/B/L and C distinguish a distinct reciprocal partner from a self-reciprocal unordered event. The latter is counted once. |
| Finite differences | D agrees that the net bracket is quadratic. Central directional differences have no algebraic second-order truncation term; shrinking the step primarily exposes cancellation here. |
| Detuning | A/B/D/L require a separate derivation for broadening. This API rejects the tested physical detuning. |
| Repeated daughters | A/B/D derive conditional results for an explicitly stipulated repeated-mode closure. The frozen implementation excludes that larger domain. Opposite-wavevector daughters in C are distinct indices and need no such extension. |
| Material coefficients | L does not identify modal linewidths with new event coefficients. This check uses only the synthetic Gamma convention. |

No unresolved mathematical disagreement was identified within the common tested domain. Repeated-index support is an explicit scope difference, not an agreement that the API implements those extended formulas. D's extreme-energy range tests and L's material normalization/export barriers remain outside this run.

## 5. Full operator and nonlinear comparisons

The frozen matrix was assembled by applying its public action method to every coordinate basis vector. The reference matrix came from the unchanged independent C implementation; no frozen factor entries were reused to construct it.

All matrix norms are spectral 2-norms.

| Comparison | Relative discrepancy |
|---|---:|
| Complete collision operator | \(1.62\times10^{-16}\) |
| Even collision block | \(1.72\times10^{-16}\) |
| Odd collision block | \(1.25\times10^{-16}\) |
| Complex-vector action | \(1.67\times10^{-16}\) |
| Complex LinearOperator matrix action | \(2.01\times10^{-16}\) |
| Energy-null residual, normalized by operator/vector norms | \(5.02\times10^{-17}\) |
| Momentum-null residual, normalized likewise | \(3.48\times10^{-17}\) |

The smallest computed eigenvalue is \(-3.45\times10^{-16}\), consistent with roundoff about an exact zero. The positive eigenvalues match the independent model.

Nonlinear updates were compared at positive perturbed populations in even, odd, and generic directions with entropy-coordinate step \(h=10^{-3}\). The largest relative RHS discrepancy was \(2.73\times10^{-12}\). Central derivatives of the frozen nonlinear RHS agreed with the independent entropy action to at worst \(9.85\times10^{-13}\).

The largest nonlinear relative discrepancy occurs for the tested odd direction, whose net update is small and involves direct-product cancellation. This run uses one moderate step; it is not a new precision-refinement study. Earlier independent C/D precision checks are not represented as newly executed here.

## 6. Full, even, and odd response

Compare all columns of
\[
 R(z,k)=(zI+C+ikV)^{-1},
\]
and separately compare its complete even-even and odd-odd blocks. Cross-parity discrepancies are also recorded.

The six cases are z = 0.25 with k = 0, 0.03, 0.1, 0.3, 1, and z = 0.25 + 0.4i with k = 0.3. The test therefore includes noncommuting streaming and complex-frequency response.

| Diagnostic over the six cases | Maximum |
|---|---:|
| Full inverse relative discrepancy | \(8.90\times10^{-16}\) |
| Even response-block relative discrepancy | \(7.77\times10^{-16}\) |
| Odd response-block relative discrepancy | \(4.44\times10^{-16}\) |
| Current response relative discrepancy | \(3.42\times10^{-16}\) |
| Energy response relative discrepancy | \(1.93\times10^{-16}\) |
| Normalized equation residual | \(7.00\times10^{-17}\) |
| Full-matrix condition number | 13.53 |

The equation residual is \(\|MR-I\|/(\|M\|\|R\|+1)\), with \(M=zI+C+ikV\). These comparisons are well conditioned. They do not certify near-pole accuracy, a continuum limit, or a material response.

## 7. Deliberate duplicate as a negative control

Append (6,3,0), the same self-reciprocal unordered event as (6,0,3), with daughters exchanged and the same prefactor. The constructor accepts all four rows and adds them.

The resulting operator agrees with the independent doubled-self-event model to \(1.65\times10^{-16}\). Relative to the correctly counted list:

| Diagnostic | Result |
|---|---:|
| Full operator relative change | 0.9627412865 |
| Even collision-block relative change | 0.9627412865 |
| Odd collision-block absolute change | \(3.79\times10^{-35}\), analytically zero |
| Relative action change on the tested even direction | 0.9866883026 |
| Energy-response relative change at z = 0.25, k = 1 | 0.0458004711 |
| Current-response relative change at z = 0.25, k = 1 | 0.0068167203 |

This reproduces C's counterexample through the frozen implementation. A homogeneous odd-only check would miss the incorrect event count; full/even checks detect it.

The implementation correctly computes the sum for the supplied rows. The invalid step is claiming a unique-event input while inserting the same event and weight twice. Different legitimate channel weights must be aggregated according to their declared convention; this result is not an instruction to silently discard every repeated tuple.

## 8. Rejections and failed approaches

Each constructor control raised ValueError:

1. **Wrong tuple order:** passing the original daughters-first C tuples directly failed with “events must conserve energy; broadening is not supported”.
2. **Detuned reciprocal pair:** add \(10^{-6}\) to both reciprocal parent energies at indices 2 and 5, preserving their reciprocal equality. With default roundoff tolerance \(10^{-12}\), construction failed with the same energy-conservation message. Tolerance-threshold behavior was not surveyed.
3. **Genuine repeated daughter:** energies (1,2) and event (1,0,0), exactly resonant, failed with “repeated indices require a separately validated counting convention”.

The distinct-index self-reciprocal event (6,0,3) succeeds.

Failed validation approaches exposed here are ignoring tuple orientation, treating physical detuning as the supported resonant model, attempting repeated-index import into this API, and using homogeneous odd agreement to certify a unique-event count. No script execution failed or required correction. The final report write alone hit a usage-limit rejection and was completed after the explicit resume.

## 9. Limits and next discriminating step

Numerically supported: for this exactly resonant, equal-weight, distinct-index seven-mode model, the frozen implementation reproduces the independent full-population entropy action and streaming response, including both parity sectors and the self-reciprocal event.

**Material export is unvalidated.** This work does not validate AlN matrix elements, physical absolute rates, event completeness, material grid/reciprocal maps, repeated-index combinatorics, broadened integration, extreme-energy numerical range, or continuum convergence.

The next importer acceptance test should explicitly account for unordered-event multiplicities and reciprocal fixed points before constructing the operator, and retain full/even checks alongside homogeneous odd response.

No repository edits were made. Independent red-team conclusions were neither read nor preempted.

