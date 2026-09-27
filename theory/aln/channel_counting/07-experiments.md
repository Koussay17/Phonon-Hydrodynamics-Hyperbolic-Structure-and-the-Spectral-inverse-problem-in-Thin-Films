# Experiments and reproduction
C Fock experiment uses interior transitions in finite sparse oscillator matrices and explicit duplicate-permutation controls.
C2 uses symbolic polynomial generator moments and a separate rational initial-state summation with no clipped destination transitions.
D calls a pinned source kernel on synthetic frequencies; a formal common spectral density compares coefficients without assigning delta(0).
PI reconstructs all six permutations of the selected material tensor using hashed saved data; two routes check opposite-q conjugation.
PI oriented-channel test includes complex random unitaries and a deliberately wrong conjugation.
Large arrays remain on D. Paths and hashes are in experiments/material-provenance.json.

Portable entry point in the repository: scripts/reproduce_channel_counting.py.
Default: finite Fock, exact moments and synthetic sewing (NumPy/SciPy/SymPy).
--source-python: pinned phono3py/phonopy 4.5.0 source audit.
--pilot-root: additionally reconstruct the selected material tensor from the retained pilot; requires --source-python.
All output goes to a new --output-dir; existing directories are rejected.
