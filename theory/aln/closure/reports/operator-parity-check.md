# Independent operator/parity check: phono3py v4.5.0

Date: 2026-09-17. Scope: commit `21fa8f3817fbcc603254656f525bb5aec113afb6`, primary equations and source inspection only. No L2 conclusion or other review was read. This is a bounded independent check; no material calculation is being run.

## Early saved verdict

**The reducible collision builder represents Chaput's conductivity-equivalent Omega, not the physical scalar-population generator Omega-prime on both parity sectors.** It has the correct intrinsic three-phonon action on odd distributions, subject to normalization, grid, and symmetry qualifications below. Its even sector cannot be used directly for energy conservation, viscosity cell problems, or finite-wave-vector kinetic evolution.

The requested path `c/reducible_collision_matrix.c` does not exist at this exact commit. The implementation is in `c/collision_matrix.c`, functions `col_get_reducible_collision_matrix` and `get_reducible_collision_matrix_at_gp`. GitHub identifies the specified commit as the merge preparing v4.5.0.

## Primary evidence already checked

1. Chaput, arXiv:1303.4062, PDF page 3, equation (2), defines the physical Omega-prime by linearizing the population collision term. It contains negative contributions with momentum constraint q-q'+q_b=G and positive contributions with q+q'+q_b=G. PDF page 4 explicitly uses the odd response f(-q)=-f(q) to replace Omega-prime by Omega. This is equality of action on the odd response, not equality of operators. Omega then has the sum of three positive delta channels. Source: https://arxiv.org/pdf/1303.4062
2. At the specified commit, `phono3py/phonon3/collision_matrix.py`, `_run_py_reducible_collision_matrix`, explicitly identifies its output as the second term of Chaput's Omega. Its arithmetic multiplies nonnegative interaction strength, inverse sinh, and `g[2]` and adds the result. It does not restore the negative channels of Omega-prime. Source: https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/collision_matrix.py
3. `phono3py/phonon3/triplets.py`, `get_triplets_integration_weights`, defines g0=delta(f0-f1-f2), g1=delta(f0+f1-f2), g2=delta(f0-f1+f2), and stores `g[2]=g0+g1+g2`. Triplets use q0+q1+q2=G. Both Gaussian and tetrahedron Python paths use this positive sum. Source: https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/phono3py/phonon3/triplets.py
4. `c/collision_matrix.c`, `col_get_reducible_collision_matrix`, advances the g pointer by two complete blocks, thus selecting Python g[2]. `get_reducible_collision_matrix_at_gp` adds strength times that weight times inverse sinh times the conversion factor. There is no parity-dependent sign change. Source: https://github.com/phonopy/phono3py/blob/21fa8f3817fbcc603254656f525bb5aec113afb6/c/collision_matrix.c

## Exact algebra behind the restriction

Let J reverse the mode wave vector, including the reciprocal branch correspondence, and let P_even=(I+J)/2 and P_odd=(I-J)/2. Separate the three-phonon formulas into the diagonal lifetime part D0, the two frequency-difference channels H with momentum constraint q+q'+q_b=G, and the sum-frequency channel Y with the same momentum constraint. Then, using reciprocal frequencies and interaction strengths,

    C_physical = Omega-prime = D0 - H J + Y,
    C_equivalent = Omega = D0 + H + Y,
    Omega - Omega-prime = H(I+J) = 2 H P_even.

Hence `(Omega-Omega-prime)P_odd=0`. There is no corresponding identity on even vectors. The step is a column reindexing and oddness substitution, not a similarity transformation of the full population operator.

For finite nonzero frequencies, the equilibrium temperature vector in entropy variables is proportional to e_q=epsilon_q/sinh(epsilon_q/(2 k_B T)), which is even and positive. The exact energy-conserving physical operator annihilates it. The conductivity-equivalent matrix has nonnegative entries for positive integration weights and positive lifetime diagonal, so it generally does not annihilate this vector. Rescaling the matrix, applying positive diagonal coordinate weights, or symmetrizing it cannot repair this structural difference.

## Normalization and export caveats (being finalized)

- Chaput uses f=sinh(epsilon/(2 k_B T)) delta_n. Standard entropy coordinates delta_n/sqrt(n0(1+n0)) equal 2f for uniform mode weights; this common factor does not change the generator or its parity restriction.
- The low-level CollisionMatrix class stores the lifetime diagonal separately from the other term. A raw row is therefore not yet the full assembled operator.
- The v4.5.0 collision-matrix kernel explicitly documents that its matrix is half the one in Chaput's paper; its main diagonal is gamma, and its conductivity code additionally handles ordinary-frequency versus angular-frequency units. Those factors must be traced before using an exported array as a time generator.
- Equal full-grid quadrature weights and reciprocal branch pairing are assumed by the simple J above. Irreducible Cartesian-vector matrices and full scalar mode matrices are different representations. Degeneracies, symmetry unfolding, finite integration widths, frequency cutoffs, and any diagonal isotope/boundary additions require separate checks.

## Practical decision

Do not certify the exported reducible array as the collision action needed by the AlN hydrodynamic campaign. It can provide the intrinsic odd-sector conductivity action after the documented normalization and assembly are verified. Energy is even; the viscous source V_i P_alpha is even because both velocity and crystal momentum are odd. Those quantities need the physical signed event operator or separately validated even-sector action. At nonzero spatial wave vector, streaming changes parity, so odd-sector equality alone does not certify a conserving nonlocal propagator.

Pending bounded checks: a single decay event plus its reciprocal partner, precise source line pointers, and the final unit factor. No peer conclusions have been consulted.