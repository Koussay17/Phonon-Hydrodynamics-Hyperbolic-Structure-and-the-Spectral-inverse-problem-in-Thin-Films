# Experiments
Independent first passes use different analytic references and preserved negative controls. Measure normalizations differ: compare definitions before numbers.
Root quadrature is checked against exact integrals; Gaussian mesh errors are separated from finite-width bias; critical scaling is not mistaken for a finite measure.
Reproduction entry point being assembled: scripts/reproduce_resonance_measure.py, new --output-dir required. It copies executable inputs only, never saved JSON outputs. Python, NumPy, SciPy, SymPy and mpmath required; no material files.
Fresh first replay: D:/ResearchLab/scratch/resonance-first-reproduction-20260928, 11 checks passed.
C2 adds the finite moment representation and derivative tests; its assertions and scope will be recorded when complete.
