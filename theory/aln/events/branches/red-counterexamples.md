# Red-team counterexample review: physical event reference

Date: 2026-09-24. **Status: COMPLETED bounded independent audit.**

Read review-target.md, the frozen code/tests/note, and branch C. No other red report was read. Repository and frozen files were not edited. The executable reproducer and JSON were saved successfully on 23 September; only the report update was interrupted by the usage-limit rejection. Their saved contents were verified on 24 September before this finalization.

## Verdict

**HIGH severity: the constructor's thermal-resonance guard admits a detuning of 1 while reporting zero.** Acceptance depends on daughter order. On the admitted input, its symmetric factor matrix differs by **44.6977%** from the high-precision entropy-coordinate derivative of its stated nonlinear Bose flux.

This falsifies the numerical acceptance contract in the explicitly tested extreme-energy regime. It does not falsify the analytical rank-one formula under genuinely exact resonance. The example is not a physical AlN-temperature claim.

## 1. Exact reproduction

Import DecayEventOperator from the frozen review-inputs/collision_events.py:

~~~python
x = [float(2**53 + 2), 1.0, float(2**53)]
op = DecayEventOperator(x, [[0, 1, 2]], [1.0])

print(op.relative_energy_residual)       # [1.11022302e-16]
print(op.detailed_balance_log_residual)  # [0.]

DecayEventOperator(x, [[0, 2, 1]], [1.0]) # raises ValueError
~~~

All three energy values are exactly representable binary64 numbers:

    (9007199254740994, 1, 9007199254740992).

Their exact stored-value detuning is

    Delta = x_parent - x_daughter1 - x_daughter2 = 1.

The nonlinear thermal forward/reverse Bose-factor ratio is exp(-1)=0.36787944117144233, rather than approximately one.

### Mechanism

The constructor evaluates

    detuning = (energies[:,0]-energies[:,1])-energies[:,2].

Subtracting 1 from 2**53+2 lands halfway between adjacent binary64 values and rounds to 2**53. The final subtraction returns zero. Swapping the daughters makes the first subtraction exactly 2 and the second exactly 1, causing rejection of the same unordered event.

The normalized residual is only 1.1102230246251565e-16 and therefore passes the default 1e-12 threshold. The absolute thermal check was intended to reject this regime, but its evaluation order defeats it. The admitted mismatch exceeds that absolute threshold by a factor of one trillion.

## 2. Independent derivative: material numerical error

The reference differentiates the quadratic population flux directly, independently of the factor formula:

    grad F = Gamma*(1+n_j+n_k, n_i-n_k, n_i-n_j).

At 100 decimal digits, evaluate Bose occupations and entropy scales at the exact stored input values, then form

    C_reference = -D^(-1) s (grad F) D.

On the active mode indices (0,2), its entries are

    [[ 1.5819767068693265, -1.5819767068693265],
     [-0.5819767068693265,  0.5819767068693265]].

The implementation instead constructs

    [[ 1.5819767068693267, -0.5819767068693265 ],
     [-0.5819767068693265,  0.21409726569788412]].

The remaining row/column underflow to zero on conversion of the final reference entries to binary64. Several displayed discrepancies are order one and are not lost tiny factors.

Executed discrepancy:

    ||C_factor-C_reference||_F / ||C_reference||_F
      = 0.44697673367510304.

The reference nonsymmetry norm is sqrt(2). This is expected away from detailed-balance equilibrium; the code incorrectly applies the equilibrium symmetric formula after admitting a nonresonant event.

### Negative control at the same large energy scale

Replacing the small daughter energy 1 by 2 gives exactly resonant stored values:

    [2**53+2, 2, 2**53].

The same independent derivative then agrees with the factor to relative Frobenius error 1.70423907900154e-16. Thus huge energies or occupation underflow alone do not explain the failed case.

The existing frozen test [3e15+100,1e15,2e15] does not exercise the relevant ordered-subtraction cancellation. The suite can pass while this admission defect persists.

## 3. Other executed attacks and negative results

These small controls were already executed before the interruption; no additional checks were added on resumption.

### Self-reciprocal counting

For modes (+q,-q,0), energies (1,1,2), and the event [2,0,1], listing [2,1,0] again duplicates the same unordered event. The independent script compared one listing with two:

| Check | Executed result |
|---|---:|
| Odd-action change on (1,-1,0)/sqrt(2) | 1.60e-17; analytically zero |
| Full collision-matrix change norm | 3.014871541977974 |
| Energy residuals, one/two listings | 2.69e-16 / 5.37e-16 |
| Current response at z=0.25, k=0.3, one listing | 3.5673008347861406 |
| Same response, duplicate listing | 3.593303644699138 |
| Relative finite-k response change | 0.00728921140023687 |

The candidate counting example survives this attack. Conservation, positivity, reciprocal parity and uniform odd-sector agreement do not diagnose that count. Finite-k exposure here is an existence result for the specified streaming and observable, not a universal detection guarantee.

The code intentionally accepts duplicate rows under caller-supplied counting conventions; this is not an undisclosed deduplication failure. **Minor documentation ambiguity:** its module-level sentence “Wave-vector reversal is a separate event” omits the self-reciprocal exception that the draft note correctly supplies.

### Sparse factor

A moderate-scale complex-vector action agreed with the small dense factor product to absolute norm 4.71e-16. No conjugation/sign defect was found in that bounded sparse-action test. It does not certify all floating-point regimes or resource scaling.

### Limits of this review

Repeated indices and unequal weights are explicitly excluded, so they supply no in-contract counterexample. Additional invariants are acknowledged by the draft; no global gap or finite-DC claim was established or attacked as if assumed. Microscopic AlN event normalization, physical broadening, material export, and boundary behavior remain unassessed.

No exhaustive input search or additional numerical test is claimed. These unassessed topics are outside this completed bounded review. Failure to find further counterexamples is not proof.

## 4. Actual executed provenance

- Reproducer: experiments/red_counter_checks.py.
- Saved results: experiments/red_counter_checks.json.
- Command: py -B followed by the full campaign path to that script.
- Exit code: 0; all assertions passed.
- Environment: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.3.0.
- Frozen imported code SHA-256:

      8daa09d8f851da3b7823bdd5cc1ad86f98ced5892c814a380dcb58ce50f020c6

  This matches review-inputs/manifest.json.

The high-precision result is an analytic derivative of the nonlinear flux, not a finite-difference or formal-verification claim. The September 23 write rejection affected the final report only: the script, result JSON, and successful execution already existed. No rejected operation was bypassed. This September 24 update closes the saved report with those verified results.

