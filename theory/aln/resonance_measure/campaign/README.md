# Conserving resonance measure — research campaign
Current status: independent first and second rounds complete; four hostile reviews pending. No candidate theorem is accepted yet. Read REVIEW_CANDIDATE.md for the frozen claims and FIRST_SYNTHESIS.md for the comparison.

## Reproduce
From the repository root:
~~~powershell
python -B scripts/reproduce_resonance_measure.py --output-dir D:/ResearchLab/scratch/resonance-new-run
~~~
Choose a new output directory. Requires Python >=3.10, NumPy, SciPy, SymPy and mpmath. No material files or network required. The runner copies Python inputs only, so stale saved JSON cannot serve as fresh evidence.
Fresh complete replay passed 12 top-level checks, including all individual script assertions. reproduction-status.json records the scope. This is bounded analytic/synthetic validation, not material convergence.

## Contents
A: regular coarea and entropy example; B: finite entropy weak form; C: independent analytic resonance benchmark; D: mesh/width asymptotics; L: primary prior art.
A2/B2 state physical/mathematical gates and alternatives. C2 checks finite nonlinear equations and derivatives, with the initial failed residual assertion and correction preserved.
The PI interface map explains why native integration arrays are not a unique positive resonant-event list.
No production collision API was changed. Draft note23 retains the preceding notes' preamble and builds with zero warnings, but awaits independent review.

## Limits
Finite quadrature energy is not automatically continuum energy. Positive entropy production is not convergence. Smooth-root formulas do not apply blindly at critical/acoustic points. No global moment-domain theorem, material run, physical hydrodynamic map, experimental validation, formal verification or novelty claim.
