# 02 — Definitions, assumptions and declared modelling choices

## A. Inherited setting (memory-gap campaign, unchanged)
* Entropy coordinates; C(g) = sum_alpha g_alpha a_alpha a_alpha^T on H = e^perp (energy invariant e = D eps);
  b = D (v eps); K = b^T C^+ b; N = |C^+ b|^2; tau_mem = N/K; tau_CS = K/|b|^2 <= tau_mem (proved);
  factor-F ceiling tau <= F / lambda_min(C_ref|H) (P7 of the memory-gap report).
* Units: g in THz (phono3py collision units); tau[ps] = (N/K)/(4 pi); ceiling 1/(4 pi lambda_min) ps.
* Event set FIXED: memory-gap AlN export (aln_events.py): 5x5x3 mesh, 897 modes, 44,406 events (decay channel of
  phono3py triplets, |w0 - w1 - w2| <= 4 sigma, sigma = 0.1 THz, 300 K), 2,627 event orbits under 6mm x time
  reversal (24 operations), 117 mode orbits. 7x7x5: 2,937 modes, 481,899 events, 23,674 / 285 orbits.
  Conserving stoichiometry and rank-one event vectors as declared there.
* Lifetimes: r = diag C (memory-gap convention, not the self-energy linewidth).

## B. The fc3 rate class
* **Rates.** g_alpha = W_alpha pp_alpha(fc3), pp = |Phi|^2 * Interaction.unit_conversion (phono3py interaction
  strength); W = conv g0w / (4 s_p s_a s_b) (distinct daughters) or conv g0w / (8 s_p s_a^2) (repeated daughter),
  s = sinh(x/2), g0w the Gaussian delta weight. W does not depend on fc3 (fc2, phonons, Bose factors, smearing
  fixed). W * pp_phono3py reproduces the memory-gap rates to 1.3e-15.
* **Unknown.** fc3 in the space S of third-order force constants of the 3x3x2 supercell satisfying the space group
  (P6_3mc), index permutation symmetry and the acoustic sum rule: symfc 1.7.3 FCBasisSetO3 without cutoff, 8,098
  symmetric coordinates, 7,752 after the sum rule; fc3 = fc3_0 + ccm B c. The hash-pinned Olympics fc3 lies in S
  (relative projection residual 4.5e-15). The 8,098 columns form 456 symmetry-distinct cluster orbits; each column
  lies in one distance shell (largest pair distance of the atom triplet).
* **Sizes of a perturbation.** (i) relative Frobenius size eps = ||fc3 - fc3_0||_F / ||fc3_0||_F (compact form);
  (ii) per-cluster relative size: max over the 456 cluster orbits of ||d fc3_o|| / ||fc3_0,o||;
  (iii) calibrated (10 October): ||d fc3_o|| / ||fc3_0,o|| <= eps_o = max(k delta_o, 0.005), delta_o the relative
  difference between the ALAMODE and phono3py Olympics workflows in orbit o (E13), k = 1, 2. Norms are taken in the
  symmetric compressed coordinates eta = B c (sum rule included).
* **Phase convention.** make_r0_average = True (phono3py 4.5.0 default, memory-gap reference); False = variant (a).
* **Orbit tying.** (i) orbit average of per-event rates (memory-gap convention; gauge dependent for events touching
  degenerate bands); (ii) block average over the degenerate block D(j0) x D(j1) x D(j2) at the representative
  triplet (gauge invariant; constant on orbits for symmetric fc3; tolerance 1e-6 THz, no frequency gap between
  1e-6 and 1e-2 THz). References: (i) kappa 214.300 / 201.967, tau 34.2530 / 27.4465 ps; (ii) kappa
  214.300 / 201.904, tau 34.2530 / 27.4410 ps. Task 3 uses (ii) as base point; phono3py-route variants use (i).
  Lifetime comparisons in the orbit-average convention are made after averaging over degenerate partners.
* **Selection-rule events.** 335 orbits with numerically zero reference rate stay zero for every symmetric fc3.

## C. Constrained problem (task 3)
* Constraints: every mode-orbit lifetime within +-dr (dr = 10 %; 5 % not run), kappa_xx and kappa_zz within +-1 %,
  plus one of: no further bound; relative Frobenius size <= rho; restriction to the cutoff-3.5 A subspace (218
  parameters); per-cluster relative size <= eps; calibrated per-orbit bounds eps_o.
* Augmented Lagrangian with L-BFGS-B; internal margins 0.998 dr, 0.98 dK, (1 - 1e-4) rho, (1 - 1e-3) eps; only
  iterates feasible for the nominal tolerances are reported. All values are attained by explicit fc3 tensors:
  inner bounds on the attainable range. Lower ends are bounded below by tau_CS (1 - dK) (rigorous).

## D. Not included / scope
* fc2 fixed (harmonic frequencies, eigenvectors, velocities): the class varies anharmonicity only.
* No four-phonon, isotope, boundary scattering; constraints at 300 K only.
* fc3 restricted to the 3x3x2 supercell range of the dataset (the ALAMODE tensor is folded into it).
* The calibration uses ONE pair of workflows (phono3py and ALAMODE teams of the Phonon Olympics; the third, ShengBTE,
  was not mapped). delta_o is a single-sample spread, not a statistical uncertainty.
* 5x5x3 and 7x7x5 are coarse meshes; statements concern these discretised operators, not converged AlN values.
