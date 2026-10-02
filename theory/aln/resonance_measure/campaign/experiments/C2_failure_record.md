# C2 initial failed assertion: preserved evidence

The initial execution of C2_weak_form.py exited with code 1. It produced no results JSON because the exception occurred while constructing the result. The first failing case below is identified from the same numerical calculation in the final saved JSON; the numerical formulas, case order, nodes, weights and parameters were unchanged. No additional experiment was run to obtain these values.

## Original execution and traceback

Command:
~~~powershell
py -B 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments\C2_weak_form.py' --output-dir 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments'
~~~

The saved conversational tool output reported:
~~~text
Traceback (most recent call last):
  File "...C2_weak_form.py", line 411, in <module>
    main()
  File "...C2_weak_form.py", line 392, in main
    "two_node":two_node_checks(),"periodic":periodic_checks(),
  File "...C2_weak_form.py", line 348, in periodic_checks
    assert abs(energy_event-positive_square) < 1e-12*max(positive_square,1e-12)
AssertionError
~~~

Paths above are abbreviated; the command supplies the actual path. The tool reported process exit code 1.

## First case that violates the original predicate

Cases run in order rho=-0.01,+0.01,-0.001,+0.001,-0.0001,... .
The first violation is rho=-0.0001.

| Quantity | Recorded value |
|---|---:|
| Prescribed energy residual rho | -0.0001 |
| Maximum actual magnitude of b e | 0.00010000000000065512 |
| Observed event heating, -sum omega F Delta | 8.1378800375857028e-11 |
| Expected thermal positive-square expression, beta sum omega Lambda Delta^2 | 8.13788003763733e-11 |
| Absolute difference | 5.1627702450672954e-22 |
| Original allowed absolute difference | 8.13788003763733e-23 |
| Maximum difference between b(beta e) and beta(b e) | 4.4408920985006262e-16 |
| Final diagnostic allowance from the measured difference plus summation rounding | 5.174342788440831e-22 |

The original condition therefore failed, and its failure is retained. The small absolute discrepancy did not meet the original requested tolerance.

## What changed, and why

The original test implicitly treated a=b(beta e) and beta Delta=beta(b e) as numerically identical. They are identical algebraically but are different floating-point multiplication orders.

Subtracting the two heating expressions, using F=-Lambda a, gives the diagnostic term
sum omega Lambda Delta (a-beta Delta).
Its absolute value is bounded by
sum omega Lambda |Delta| |a-beta Delta|.

The revised script records that positive absolute sum and adds
32 machine_epsilon (|observed_heating|+|positive_square|)
for the subsequent arithmetic. The comparison now uses this explicit scale-aware allowance. The JSON records both heating expressions, their discrepancy, the allowance, and the affinity multiplication-order residual for every residual case.

No physical coefficient, quadrature, energy, case, root, or collision formula was changed. No case was deleted. The final rerun passed. This is a numerical-conditioning diagnosis supported by the measured residual, not a machine-verified floating-point error theorem.

Final artifacts: C2_weak_form.py and C2_weak_form.json.
Final executed script SHA-256:
933e7ae55004274b2c33c6fcc709879bc1f49bb1bd7b3a5cf425d69e7c136e5c.

