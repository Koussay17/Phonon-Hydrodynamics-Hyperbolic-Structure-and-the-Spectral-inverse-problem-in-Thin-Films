import ResearchProofs.FdtrGauge

/-! Axiom report for the statements of `ResearchProofs.FdtrGauge` (run with `lake env lean`). -/

-- Theorem 1, Steps 4–6
#print axioms FdtrGauge.crossmul_iff_coeffs
#print axioms FdtrGauge.step5
#print axioms FdtrGauge.step5_iff
#print axioms FdtrGauge.crossmul_of_two_points
#print axioms FdtrGauge.step5_two_points
#print axioms FdtrGauge.step5_g
#print axioms FdtrGauge.step4
#print axioms FdtrGauge.mem_fibre_iff_of_not_D2
#print axioms FdtrGauge.thm1_of_steps
#print axioms FdtrGauge.hasDerivAt_g_tR
#print axioms FdtrGauge.hasDerivAt_g_tl
#print axioms FdtrGauge.g_tangent_resonance
-- Two-load kernel lemma and Theorem 2, Proof A
#print axioms FdtrGauge.twoLoad_minor01
#print axioms FdtrGauge.twoLoad_minor02
#print axioms FdtrGauge.twoLoad_minor12
#print axioms FdtrGauge.twoLoad_minors_eq_zero_iff
#print axioms FdtrGauge.twoLoad_ker_iff
#print axioms FdtrGauge.quadrupole_mem_ker
#print axioms FdtrGauge.ker_line
#print axioms FdtrGauge.step3_pointwise
#print axioms FdtrGauge.thm2_pointwise
#print axioms FdtrGauge.F_sub_inv_S
#print axioms FdtrGauge.F_eq_inv_S_iff
#print axioms FdtrGauge.thm2_endgame
#print axioms FdtrGauge.thm2_kernel_iff
#print axioms FdtrGauge.thm2_D2_kernel_iff
-- Moment coefficients and gradients
#print axioms FdtrGauge.Moments.c0_eq_proofC
#print axioms FdtrGauge.Moments.c1_eq_proofC
#print axioms FdtrGauge.Moments.c2_eq_proofC
#print axioms FdtrGauge.Moments.c3_eq_proofC
#print axioms FdtrGauge.Moments.grad_c0
#print axioms FdtrGauge.Moments.grad_c1
#print axioms FdtrGauge.Moments.grad_c2
#print axioms FdtrGauge.Moments.grad_c3
#print axioms FdtrGauge.Moments.grad_c4
#print axioms FdtrGauge.Moments.grad_c5
#print axioms FdtrGauge.Moments.grad_c6
-- Moment minors (Proof C, Theorem N2)
#print axioms FdtrGauge.Moments.G_0123
#print axioms FdtrGauge.Moments.G_0124
#print axioms FdtrGauge.Moments.G_0234
#print axioms FdtrGauge.Moments.G_0234_matched
#print axioms FdtrGauge.Moments.G_01_matched
#print axioms FdtrGauge.Moments.G_0125
#print axioms FdtrGauge.Moments.G_0134
#print axioms FdtrGauge.Moments.G_0126
#print axioms FdtrGauge.Moments.G_0135
#print axioms FdtrGauge.Moments.moment_minor_ne_zero
#print axioms FdtrGauge.Moments.moment_kernel_trivial
#print axioms FdtrGauge.Moments.resonance_direction
-- Theorem N2, step (3)
#print axioms FdtrGauge.LeadingOrder.t_0123
#print axioms FdtrGauge.LeadingOrder.t_0124
#print axioms FdtrGauge.LeadingOrder.t_0234
#print axioms FdtrGauge.LeadingOrder.t_0124_ne_zero
#print axioms FdtrGauge.LeadingOrder.t_0234_ne_zero
#print axioms FdtrGauge.LeadingOrder.index_sets_le_nine
#print axioms FdtrGauge.LeadingOrder.index_sum_ge_six
#print axioms FdtrGauge.LeadingOrder.index_sum_six
#print axioms FdtrGauge.LeadingOrder.index_sum_seven
#print axioms FdtrGauge.LeadingOrder.lowSets_contain_01
#print axioms FdtrGauge.LeadingOrder.leading_coeff_ne_zero
#print axioms FdtrGauge.LeadingOrder.matched_low_order
-- Riccati hierarchy
#print axioms FdtrGauge.Riccati.Zfilm_zero
#print axioms FdtrGauge.Riccati.hasDerivAt_Zfilm
#print axioms FdtrGauge.Riccati.hasDerivAt_Zfilm_model
#print axioms FdtrGauge.Riccati.gam_series
#print axioms FdtrGauge.Riccati.hierarchy
#print axioms FdtrGauge.Riccati.hierarchy_init
#print axioms FdtrGauge.Riccati.hierarchy_unique
