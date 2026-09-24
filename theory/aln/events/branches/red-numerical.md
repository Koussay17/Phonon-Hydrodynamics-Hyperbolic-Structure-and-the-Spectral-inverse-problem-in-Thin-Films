# Independent red-team numerical audit

**Complete, 2026-09-23.** Frozen target: `review-target.md` and the three files under `review-inputs/`. No other red-team report was read or discussed. The repository and frozen candidate were not edited. Reproductions import the frozen implementation directly.

## Verdict

**One high-severity numerical contract failure was reproduced:** the absolute resonance check can erase a daughter energy, accept an event with exact thermal detuning -1, and expose a false zero detuning. The resulting symmetric sparse matrix differs by **59.325%** from the nonlinear population Jacobian for those accepted inputs. Swapping the daughters changes acceptance.

The remaining findings concern zero-tolerance behavior, inconsistent vector/batch validation, and limits of the claimed reference precision. They do not require expansive experiments or actual AlN data.

## 1. HIGH: absolute detuning is not reliably measured

**Proposition attacked:** accepted events satisfy both the relative energy check and the absolute thermal-detuning check to `resonance_rtol`, and the exposed `detailed_balance_log_residual` accurately reports that mismatch.

Frozen `collision_events.py`, lines 49-53, computes

```python
residual = x_parent/scale - x_daughter1/scale - x_daughter2/scale
detuning = (x_parent - x_daughter1) - x_daughter2
```

Small reproduction:

```python
op = DecayEventOperator([1e16, 1., 1e16], [[0, 1, 2]], [1.])
```

All three energies are finite, positive, and exactly representable binary64 values; the mode indices are distinct. The call **accepts** and reports:

- relative energy residual: `-1.1102230246251565e-16`;
- detailed-balance log residual: `0.0`;
- default tolerance: `1e-12`.

But the exact mismatch of the stored inputs is

    Delta = 1e16 - 1 - 1e16 = -1.

The first subtraction rounds away the unit daughter energy. The forward/reverse Bose-factor ratio is exp(-Delta)=e, and the normalized equilibrium net flux is e-1=1.718281828459045. This is not a thermal detuning bounded by the supplied tolerance.

The ordering-equivalent call

```python
DecayEventOperator([1e16, 1., 1e16], [[0, 2, 1]], [1.])
```

**rejects** with `ValueError`. The physical daughter permutation changes neither the event nor its exact detuning, but changes the numerical admission decision.

A second reproduction, `[1e5, 5e-12, 1e5]` with event `[0,1,2]`, also accepts and reports both residuals zero even though the exact stored-input detuning is approximately -5e-12, exceeding the default absolute tolerance.

### Consequence: accepted sparse action is not its nonlinear Jacobian

The first accepted input returns, to the shown digits,

    C_candidate =
      [ 1.581976706869327  0  -1.581976706869327 ]
      [ 0                  0   0                 ]
      [-1.581976706869327  0   1.581976706869327 ].

A 100-digit calculation of the actual population-polynomial derivative, without imposing equilibrium detailed balance, gives

    C_nonlinear =
      [ 1.581976706869326  0  -0.581976706869326 ]
      [ negligible         0   negligible        ]
      [-1.581976706869326  0   0.581976706869326 ].

The omitted entries are retained in the JSON at their tiny nonzero magnitude. Relative Frobenius discrepancy is **0.593250138083519352**. The mathematical exact-resonance identity is not the proposition falsified here: these inputs are nonresonant and should fail the advertised validation. The implementation admits them and then substitutes a different, symmetric PSD operator.

**Why existing evidence misses it:** the negative-control test at `test_collision_events.py:112` checks a large mismatch that survives its subtraction order. The high-precision tests use balanced energy ratios `(3,1,2)`, not a large daughter plus a thermally significant unresolved small one. Symmetry and PSD tests cannot detect this failure because `factor.T @ factor` has those properties by construction, even for mistakenly accepted events.

**What remains unestablished:** the frozen constructor cannot be relied on to enforce or faithfully report the absolute detailed-balance tolerance over its advertised positive finite energy range. The accompanying statement that absolute mismatch is controlled is too strong for this implementation.

## 2. LOW: an allowed zero tolerance rejects exactly resonant data

The constructor permits `resonance_rtol=0`, but

```python
DecayEventOperator([3., 1., 2.], [[0, 1, 2]], [1.], resonance_rtol=0)
```

rejects, although the exact stored-input detuning is zero. The division-based relative residual introduces roundoff before it is compared with zero. Swapping the daughters also rejects this particular input.

This is an edge-case contract mismatch. It does not falsify the resonance formula, but the allowed tolerance interval is not a guarantee that exactly resonant represented data are accepted at its lower endpoint.

## 3. LOW: sparse vector and batch paths disagree on finite-input validation

`action()` checks finiteness at `collision_events.py:88`; `as_linear_operator()` routes vector calls through that check, while `matmat()` at line 100 multiplies the sparse factor directly.

For an ordinary resonant three-mode operator:

```python
linear @ np.array([np.nan, 0., 0.])
```

raises `ValueError`, whereas

```python
linear @ np.array([[np.nan, 1.], [0., 0.], [0., 0.]])
```

returns a nonfinite `(3,2)` array without rejecting it. The batch adjoint path uses the same direct multiplication. Thus callers do not receive consistent finite-input validation when switching to multiple right-hand sides.

This is an API error-handling defect, not a demonstrated wrong result for valid finite vectors. The small reproduction does not claim a general failure of SciPy's `LinearOperator` protocol.

## 4. MODERATE evidence limitation: 100-digit reference does not imply 100-digit validation

`test_collision_events.py:35` constructs a nonlinear mpmath derivative at 100 digits, but line 46 casts the entire reference to binary64. Line 47 then compares with `rtol=2e-12, atol=1e-300`.

The supported claim is a binary64 comparison to a high-precision-derived reference on four selected energy scales. It does not establish:

- 100-digit accuracy of the sparse action;
- relative accuracy of entries below the absolute tolerance;
- preservation of unrepresentable occupation factors or entropy-coordinate invariants;
- uniform accuracy across energy ratios or admitted detunings;
- microscopic normalization or material event counting.

The nonlinear derivative is a separate arithmetic path from the sparse factor, but both implement the same stipulated Bose population law and the same user-supplied prefactor convention. It is not an independent derivation of that microscopic convention. The absolute-resonance counterexample above demonstrates a gap left by the sampled reference tests.

### `atol=0` does not test exact equality

At `test_collision_events.py:60`, the parity check passes `atol=0` but leaves the relative tolerance at its default. In the audit runtime,

```python
np.testing.assert_allclose(1.00000005*c, c, atol=0)
```

passes for the ordinary test matrix. Therefore that assertion cannot support bitwise or exact numerical equality. This finding concerns the strength of the assertion, not a newly observed parity defect in the supplied reciprocal pair.

## 5. Documented cancellation/underflow limits remain real

The frozen nonlinear `population_rhs()` explicitly restricts its intended use to moderate occupations. A small reproduction using x=(3e-18,1e-18,2e-18), positive relative population perturbations, and a centered directional difference returns an **entirely zero derivative**, while the sparse action is nonzero: relative disagreement is 1.0. The base-state RHS is also exactly zero in binary64.

This is a retained failure of the raw-product reference outside its stated moderate-input regime, not a new in-contract moderate-input bug. It proves that a vanishing equilibrium residual alone cannot certify that reference derivative. The earlier independent D precision and asymptotic results explain the cancellation and are retained separately.

Likewise, high-energy tests that permit underflow may validate a normwise limiting matrix while losing tiny entries and entropy scales. They do not establish a uniformly invertible population-to-entropy coordinate map. The frozen target acknowledges underflow; no stronger accuracy claim is supported by the existing tests.

## 6. Reproduction and numerical classification

Owned artifacts:

- `experiments/red_numerical_api.py`
- `experiments/red_numerical_api_results.json`
- this report

Run:

```powershell
py -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260917-232509-aln-physical-events\experiments\red_numerical_api.py
```

This script imports `review-inputs/collision_events.py` directly and performs only small three-mode tests. It does not import or modify the live repository. The audited frozen implementation SHA-256 is

`8daa09d8f851da3b7823bdd5cc1ad86f98ced5892c814a380dcb58ce50f020c6`.

Runtime metadata: Python 3.14.7, NumPy 2.5.3, mpmath 1.3.0. Full reproduction outputs and exception messages are saved in the JSON.

- **Consistency failure:** the accepted Delta=-1 event yields an operator different from its actual nonlinear entropy-coordinate Jacobian by 59.325%.
- **Floating-point cancellation:** the thermal-detuning diagnostic reports a false zero and depends on daughter ordering. Raw-product differentiation also fails at small energies as documented.
- **Conditioning:** the detuning is small relative to huge energies but not small in thermal units; fractional energy accuracy therefore does not control detailed balance. The exposed residual loses exactly the information the second check was meant to retain.
- **Stability:** PSD construction does not establish physical consistency of an erroneously admitted event. No claim of dynamic spectral instability is inferred from these tests.
- **Convergence:** no new mesh or timestep convergence study was performed. The bounded audit establishes counterexamples, not uniform error bounds.
- **Tolerance-limited evidence:** comparisons after binary64 conversion certify only the specified mixed relative/absolute tolerance at sampled inputs.
- **Inconclusive:** actual AlN event validity, absolute prefactors, grid/broadening convergence, and material response.

There is no time integrator, PDE spatial grid, finite-domain boundary problem, or iterative collision solver in these reproductions. CFL, temporal/spatial discretization error, boundary contamination, and nonlinear spectral aliasing do not apply. The full repository suite and PDF were not rerun; their reported success does not remove the reproduced counterexample.

## Required interpretation of this audit

The frozen candidate has a concrete failure in the validation condition that defines its physical regime. Existing pass counts, a small reported detuning, and a symmetric PSD sparse action are insufficient to certify that an accepted event lies in that regime. No code improvements or repairs were made during this independent review.
