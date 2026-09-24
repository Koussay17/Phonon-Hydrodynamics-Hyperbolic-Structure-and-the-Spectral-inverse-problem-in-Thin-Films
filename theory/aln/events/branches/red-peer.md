# Independent hostile peer review — final bounded report

**Date:** 2026-09-23. **Scope:** frozen `review-target.md`, note 19, `src/collision_events.py`, its tests, C/D first-pass reports, L sources, and the C-drive download manifest. No other red report read, no cross-agent discussion, no repository edits. Only small read-only reproductions were run.

## Verdict

No fatal algebraic error was found in the stated **exactly resonant, distinct-mode, finite Bose-event theorem**. Reject any broader claim that this validates physical AlN events, a production importer, or macroscopic transport. The main concrete defect is a mismatch between an exact-Jacobian description and the API's admitted nonzero-detuning inputs. Its small default size does not turn the approximation into an identity. The strongest parity result remains a finite implementation counterexample; its material significance is unestablished.

## Major objections

### M1 — Accepted detuning breaks exact Jacobian equivalence

**Proposition:** every accepted input implements the stated nonlinear Bose flux linearized at equilibrium.

**Counterexample:** construct `DecayEventOperator([3+5e-13,1,2], [[0,1,2]], [1])` with the default tolerance. It is accepted. For the independent exact quadratic-flux gradient

    g = (1+n_j+n_k, n_i-n_k, n_i-n_j),
    C_J = -outer(s/d, g*d),   s=(-1,1,1),

this review measured:

| Diagnostic | Value |
|---|---:|
| Stored detuning | 5.000444502911705e-13 |
| Relative difference between R.T R and C_J | 3.2048286208645094e-13 |
| Relative asymmetry of C_J | 4.043173010189e-13 |
| Nonlinear Bose-state drift norm | 8.304795101819017e-14 |
| Norm of R applied to the energy vector x*d | 1.546540673302843e-13 |

With detuning `5e-9` and permitted `resonance_rtol=1e-8`, the matrix discrepancy is `3.2033468960218787e-9`, the asymmetry `4.04231477940278e-9`, and the drift `8.30192051297321e-10`.

These are expected consequences of detuning, not evidence of a large practical error at the default tolerance. They decisively limit the theorem/API correspondence: the admitted matrix is a positive surrogate when the supplied numbers are nonresonant. The test cannot determine whether a mismatch arose from rounding or physical detuning. Unweighted mismatch diagnostics do not bound a low-frequency inverse-response error when a gap is small.

**Status:** NUMERICALLY DEMONSTRATED; analytically explained by the note's own negative control. **Severity:** major contract/interpretation issue; not a refutation of the exact-resonance theorem. The frozen target partly acknowledges it, but the note's exact language must not be applied to all accepted arrays.

### M2 — The material import barrier is not just a unit conversion

The implementation receives no wavevectors, branches, reciprocal vectors, crystal selection rules, full-grid maps, or cubic-amplitude convention. Thus its PSD construction cannot distinguish an allowed AlN event from an energy-matched but momentum-forbidden or symmetry-forbidden triple. User-supplied Gamma deliberately absorbs the unresolved factors.

The more serious barrier is energy integration. Generic finite-grid frequencies need not lie exactly on the continuum resonance surfaces. Rejecting every detuned triple may remove real continuum scattering and create false invariants. Accepting broadened triples does not satisfy the exact nonlinear theorem. No source-to-event normalization or conserving energy-integration method has been validated here.

This distinction is present in the primary literature: Fugallo et al., Section IV after Eq. (30), state that finite Gaussian broadening makes detailed balance approximate and changes the equivalence of two collision definitions; their chosen rearrangement preserves symmetry/positivity. I independently inspected that passage in the [actual author manuscript](https://arxiv.org/html/1212.0470v2). It does not justify identifying an arbitrary broadened matrix with the unmodified Bose Jacobian.

**Status:** material-export validity UNSUPPORTED. **Severity:** fatal only to material validation/import-completion claims, which the note correctly disclaims. Synthetic differentiation and a large repository test count do not resolve this barrier.

### M3 — The finite parity example does not establish a persistent bulk error

The self-reciprocal event is invisible in the homogeneous odd block by construction. But an exactly resonant positive-energy decay has a unique highest-energy parent. If reversal maps the canonical event to itself, it must fix that parent; its wavevector therefore satisfies `q_parent = -q_parent modulo G`. Such parents occupy isolated time-reversal-invariant wavevectors, including Gamma.

Under regular continuum weights these special events can have vanishing weight in a bulk observable as the grid is refined. Singular rates, degeneracies, finite-size physics, or specific observables could change that conclusion; none is established here. The three/seven-mode examples do not settle this scaling. This is a reason to withhold an AlN error estimate, not a refutation of the finite counterexample.

Likewise, nonzero spatial wavevector is not a universal detector: the changed even subspace must be reached by the actual streaming operator and probed by the chosen source/observable. C demonstrates that condition for its assigned velocities. No AlN dispersion realizes those velocities/events in the evidence presented.

**Status:** finite toy differences NUMERICALLY DEMONSTRATED; macroscopic persistence and material magnitude UNSUPPORTED. **Severity:** major only if the example is promoted beyond its stated validation scope.

### M4 — Finite stored action does not ensure usable entropy coordinates

Read-only reproduction of the existing extreme case

    DecayEventOperator([3e16,1e16,2e16], [[0,1,2]], [2])

returns `entropy_scale = [0,0,0]`, while the collision diagonal is approximately `[2,0,0]`. The asymptotic decay action remains finite, but the stored population-to-entropy map is noninvertible and its stored energy vector vanishes. A zero conservation residual or population reconstruction through those arrays would be meaningless.

The frozen target already warns of underflow. That warning must remain attached to any claim of validated full-population coordinates or invariants in this regime. The high-precision comparison tests a limiting generator; it does not restore the binary64 coordinate map.

**Status:** NUMERICALLY DEMONSTRATED. **Severity:** medium numerical-contract limitation, not a failure of the limiting decay coefficient.

## Minor objections

1. The module docstring says wavevector reversal is a separate event without the self-reciprocal exception that motivates the strongest counterexample. The literal wording can instruct the very double counting the note diagnoses.
2. The nonresonant energy-production formula is strictly positive only for an active event, Gamma>0. With the allowed Gamma=0 it is zero.
3. The seven-mode toy has accidental invariants; its resolvent uses positive real z. Neither a finite DC conductivity nor a harmonic experimental response is validated by those numbers.
4. The campaign assumptions file still says pending. The actual exclusions reside in other files, so that ledger cannot yet be cited as the finalized assumptions record.
5. The current review did not rerun the complete 225-test suite or compile the PDF. Those totals remain reported campaign evidence. It independently inspected the 20-test file and reproduced the two contract limitations above.
6. The C-drive manifest records the pinned revision, file sizes and SHA-256 hashes for the downloaded AlN inputs. The D-drive payloads are currently unavailable in this review environment, so their present existence/content could not be independently rehashed. This does not contradict the historical download record.

## Prior art and evidence classification

The positive signed-event decomposition and Bose-coordinate congruence are established collision theory; the note makes no novelty claim. No broad prior-art search was needed for this bounded validation review, and no novelty conclusion follows. The inspected Fugallo reference supports the convention/smearing cautions attributed to it.

- Exact-resonance identity: complete conditional derivation; no contradiction found in this review. Final PROVED classification belongs to the independent proof audit. No formal verification performed.
- Reference implementation on the displayed cases: NUMERICALLY DEMONSTRATED, with the tolerance and underflow limits above.
- Odd-only checks can miss a finite self-reciprocal event miscount: conditional algebra plus NUMERICALLY DEMONSTRATED toy evidence.
- Actual AlN event rates, importer normalization, continuum conservation, hydrodynamic window, or experimental validity: UNSUPPORTED by this task.

## Strongest parts

The exact resonant theorem is explicit, the input exclusions are mostly clear, and the negative detuning control correctly separates PSD structure from physical conservation. The strongest finite counterexample survived this bounded attack.

## Fatal issues, if any

None found for the expressly restricted finite theorem. Calling the API an exact nonlinear Jacobian for every accepted input is false; calling the work a completed material importer or AlN validation would exceed the evidence.

## Major issues

Accepted-detuning approximation, unresolved physical event support/normalization/energy integration, unknown continuum significance of the parity toy, and underflow-induced loss of the coordinate map.

## Minor issues

Reciprocal-orbit wording, Gamma=0 strictness, incomplete assumptions ledger, and limits on independently checked test/download evidence.

## Decisive tests

Audit the exact theorem/API domain distinction; independently map one pinned exporter channel including all units and multiplicities; test its energy integration and invariants under physical mesh/broadening refinement; establish whether the self-reciprocal defect persists in a normalized observable; test invariant/reconstruction diagnostics when entropy scales underflow. These are required evidence, not extensions performed by this reviewer.

## Unresolved questions

What conserving discrete integration will connect the ideal reference to actual mesh data? Which omitted repeated-index or zero-mode channels matter? Does the self-reciprocal error persist at useful magnitude in AlN? What error norm and gap assumptions turn resonance tolerances into a response bound? None is resolved by the present evidence.
