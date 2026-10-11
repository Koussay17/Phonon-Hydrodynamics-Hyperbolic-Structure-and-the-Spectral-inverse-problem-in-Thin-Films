# 02 — Definitions, assumptions and declared modelling choices

## A. Inherited setting (memory-gap campaign, unchanged)
* Entropy coordinates, event operator C(g) = sum_alpha g_alpha a_alpha a_alpha^T on H = e^perp (energy invariant
  e = D eps), current b = D (v eps), K = b^T C^+ b, N = |C^+ b|^2,
  tau_mem = N/K = b^T (C^+)^2 b / b^T C^+ b, tau_CS = K/|b|^2 (proved lower bound on tau_mem),
  factor-F ceiling tau <= F / lambda_min(C_ref|H) (P7 of the memory-gap report).
* Units: g in THz (phono3py collision units); tau[ps] = (N/K)/(4 pi); lambda_min ceiling 1/(4 pi lambda_min) ps.
* Event set: the memory-gap AlN export (aln_events.py), FIXED: 5x5x3 mesh, 897 modes, 44,406 events
  (p -> a + b decay channel of phono3py triplets, |w0 - w1 - w2| <= 4 sigma, sigma = 0.1 THz Gaussian, 300 K),
  2,627 event orbits under 6mm x time reversal (24 operations), 117 mode orbits; 7x7x5: 2,937 modes,
  481,899 events (review export). Conserving stoichiometry and rank-one event vectors as declared there.
* Lifetimes = r = diag C (memory-gap convention, not phono3py's self-energy linewidth).

## B. The fc3 rate class (this campaign)
* **Rates.** g_alpha = W_alpha pp_alpha(fc3), pp = |Phi|^2 * Interaction.unit_conversion (phono3py interaction
  strength), W_alpha = conv * g0w / (4 s_p s_a s_b) (distinct daughters) or conv * g0w / (8 s_p s_a^2) (repeated
  daughter), s = sinh(x/2), g0w = Gaussian delta weight. W is independent of fc3 (fc2, phonons, Bose factors,
  smearing fixed). Check: W * pp_phono3py reproduces the memory-gap rates to 1.3e-15 relative.
* **Unknown.** fc3 in the space S of third-order force constants of the 3x3x2 supercell that satisfy the space-group
  symmetry (P6_3mc), index-permutation symmetry and the acoustic sum rule. Basis: symfc 1.7.3 FCBasisSetO3 (no
  cutoff): 8,098 symmetric coordinates, 7,752 after the sum rule. fc3 = fc3_0 + ccm B c.
  The hash-pinned Olympics fc3 lies in S to 4.5e-15 (relative projection residual), so the reference is a member.
* **Size of a perturbation.** eps = ||fc3 - fc3_0||_F / ||fc3_0||_F in the compact (4 x 72 x 72 x 27) form
  (= |c| / |c0| in the orthonormal basis coordinates).
* **Phase convention.** phono3py make_r0_average = True (4.5.0 default, used by the memory-gap reference); False
  is variant (a).
* **Orbit tying.** Two conventions:
  (i) orbit average of per-event rates (memory-gap convention; gauge dependent for events touching degenerate
      bands: 128 of 2,292 allowed orbits differ from (ii) by up to a factor 18);
  (ii) block average: mean of |Phi|^2 over the degenerate block D(j0) x D(j1) x D(j2) at the representative
      triplet (gauge invariant; constant on orbits for symmetric fc3; degeneracy tolerance 1e-6 THz, no frequency
      gap between 1e-6 and 1e-2 THz on this mesh).
  Reference values: (i) kappa 214.300 / 201.967 W/(m K), tau 34.2530 / 27.4465 ps;
  (ii) kappa 214.300 / 201.904, tau 34.2530 / 27.4410 ps. Task 3 uses (ii) as base point; the phono3py-route
  variants use (i).
* **Selection-rule events.** 335 orbits with numerically zero reference rate (<= 1e-9 of the max) stay zero for
  every symmetric fc3 (their reference values, <= 4.4e-15 relative, are kept).

## C. Constrained problem (task 3)
* Constraints: every mode lifetime within +-dr (dr = 5 % or 10 %) of the base point (equivalently the 117
  mode-orbit representatives, rates being symmetric), kappa_xx and kappa_zz within +-1 %, and optionally
  eps <= rho (rho = 5, 10, 20, 50, 100 %, or unbounded).
* Internal margins: the optimiser works with 0.998 dr, 0.98 dK, (1 - 1e-4) rho, and only iterates satisfying the
  NOMINAL tolerances are reported. All values are attained by explicit fc3 tensors (inner bounds).

## D. Not included / scope
* fc2 (harmonic frequencies, eigenvectors, velocities) is fixed: the class varies anharmonicity only.
* No four-phonon, isotope or boundary scattering; one temperature (300 K) for the constraints.
* fc3 restricted to the 3x3x2 supercell range of the dataset (no longer-range terms).
* 5x5x3 (and 7x7x5) are coarse meshes; the statements concern these discretised operators, not converged AlN
  values (kappa_LBTE at 31x31x17 is 302 / 286 W/(m K)).
* The conserving event stoichiometry of the memory-gap campaign is kept as declared there.
