# Experiments and reproduction
C_lumpability_check.py: exact SymPy rational algebra and finite SciPy exponentials for a four-mode reversible event network. Files live in experiments/. Script writes results beside itself; run a copy in scratch to preserve the archive.
D-phase-information.py: input-free NumPy phase ambiguity and conditional kinetic examples; use --output with an absolute scratch JSON path.
experiments/reconstruct_complex_triplet.py: portable selected AlN amplitude contraction. Requires prior corrected-export workspace from the interaction_pilot campaign, pinned phono3py/phonopy4.5.0 and identical input/output hashes. Supply --pilot-root and a new --output-dir. Does not overwrite an existing output directory.
Original exploratory complex_triplet.py and its provenance are retained. Its hardcoded D paths are historical, not the portable entry point.
No timestep or grid refinement is implied by exact matrix exponentials or one-triplet tensor checks. The material calculations retain the previous3x3x3 grid and source assumptions.
