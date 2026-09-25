# Literature and source scope
Primary source and unit trace: branches/C-unit-source-trace.md, pinned phono3py commit 21fa8f3817fbcc603254656f525bb5aec113afb6 (4.5.0).
- Togo, Chaput and Tanaka, Phys. Rev. B 91, 094306 (2015): https://arxiv.org/abs/1501.00691
- Official exported linewidth convention: https://phonopy.github.io/phono3py/input-output-files.html#gamma
- Input force constants: PhononOlympics commit 0640f07735059be9717a7565c2a0f22dc0da7a17, Aluminum Nitride/phono3py/AlN_kappa_input_files; input hashes in pilot-run.json.
- Installed phonopy 4.5.0 harmonic/dynamical_matrix.py: Gonze reciprocal construction and G_cutoff/Lambda parameters. The source explicitly notes that its dipole term is not strictly periodic over G; this motivated the controlled numerical sweep.
This is a targeted implementation/prior-art check, not an exhaustive novelty search. No novelty is claimed.

Pinned harmonic implementation: https://github.com/phonopy/phonopy/blob/v4.5.0/phonopy/harmonic/dynamical_matrix.py (live source checked 25 September 2026). Online phono3py documentation currently displays 4.4.0; version-specific conclusions above use installed/pinned 4.5.0 source.
