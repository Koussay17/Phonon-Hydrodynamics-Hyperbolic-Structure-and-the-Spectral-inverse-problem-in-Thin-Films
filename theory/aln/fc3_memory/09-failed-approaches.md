# 09 — Failed approaches, bugs and negative results

1. **phono3py's pure-Python real-to-reciprocal prototype disagrees with its production kernels.**
   `phono3py.phonon3.real_to_reciprocal.RealToReciprocal` (4.5.0) builds the pre-phase from
   `np.where(np.all(triplet != 0, axis=0), triplet.sum(axis=0), 0)`: a Cartesian component of G = q0 + q1 + q2 is
   dropped unless all three addresses are non-zero in it. Re-implementing that rule, 1,062 of 41,131 allowed events
   (38 of 1,613 triplets, all Umklapp with G = (0, +-5, 0)) disagreed with the Rust interaction strengths by up to a
   factor 440. With the full G in every component the engine reproduces the Rust kernel for every event (max
   relative difference 6.5e-12 with r0 averaging, 2.5e-11 without). Use the C/Rust path for Umklapp-sensitive work.

2. **Orbit averaging of per-event rates is not gauge invariant.** For events touching degenerate bands the
   memory-gap orbit average keeps a gauge-dependent split between degenerate partners: 128 of 2,292 allowed orbits
   differ from the gauge-invariant block average by up to a factor 18 (5x5x3; 733 orbits up to a factor 207 at
   7x7x5). Effect on the reference: kappa_zz -3e-4, tau_z -2e-4 (5x5x3); x unchanged to 1e-9. For lifetimes the
   effect is visible: for the tau_x = 481 ps witness one doubly degenerate mode (19.97 THz) moves by 58 % in the
   orbit-average convention while all degeneracy-averaged lifetimes move by at most 9.98 % (the tolerance); 44 %
   versus 9.98 % for the tau_z = 17,121 ps witness. Per-mode lifetimes of degenerate partners in the orbit-average
   convention are gauge artefacts (22 degenerate modes at 5x5x3; this also qualifies the memory-gap phrase
   "all 897 lifetimes fixed").

3. **fc3 storage:** the dataset fc3.hdf5 is full (72 x 72 x 72 x 3 x 3 x 3); the compact tensor is fc3[p2s]
   (first symfc projection attempt failed on the shape).

4. **Per-event outer ceiling is vacuous beyond a 0.5 % Frobenius change.** tau <= 1/lambda_min(C(g_low(rho)))
   with per-amplitude floors |Phi| >= |Phi0| - ||G_m|| rho |c0|: 81.5 ps (rho = 0), 481 ps (0.1 %), 3,287 ps
   (0.2 %), infinite from 0.5 % on (individual amplitudes can be silenced by tiny worst-case changes; fragility
   median 0.26 %, max 2.2 %). Rigorous but far from the attainable values.

5. **The Frobenius norm is a misleading plausibility measure.** The tau_x = 481 ps and tau_z = 17,121 ps witnesses
   have Frobenius size 9.7 % / 9.6 % but leave the on-site and nearest-neighbour clusters (99 % of ||fc3_0||^2)
   changed by only 1.6 - 2.3 %, change the 2 - 3.2 A shell by 34 - 35 % and rewrite the clusters beyond 3.5 A
   (0.2 % of ||fc3_0||^2) by factors 1.1 - 5.0 (87 % of the change). A small global norm hides order-one relative
   changes of small long-range constants; the physically meaningful classes are short-range subspaces and
   per-cluster relative bounds.

6. **Local optimiser, non-convex problem.** The unbounded full-space run for max tau_x stopped at 453.5 ps although
   the rho <= 0.1 run (a subset of the unbounded class) reached 481.4 ps; the 3.5 A run for max tau_x with
   rho = 0.05 ended without a nominally feasible iterate (25 outer iterations). All maxima are inner bounds; for
   nested classes the best value over all runs is quoted.

7. **Code bug (caught before use):** a patch inserted two methods in the middle of FC3Model.__init__ (D: copy),
   truncating the constructor; the first per-cluster launch failed at once and was relaunched after the repair.

8. **Per-cluster bounds in Frobenius coordinates are ill-conditioned.** Cluster reference norms span 0.058 - 1,166;
   in the coordinates y = c/|c0| the smallest clusters allow steps of order 1e-6 while the objective varies on
   order 1, and L-BFGS-B stopped after one iteration per outer step (first D: batch: no result in 60 min of
   4 x CPU; first C: batch: 'nit' = 1). Fixed by the change of variables u = M^(1/2) c, M = B^T D^-2 B
   (orbit_scaling.py), in which every per-cluster constraint is ||(W u)_o|| <= eps with W isometric
   (cond(M) = 4.0e8; isometry error 3e-13).

9. **Loss of drive D: (external) at about 16:50 on 5 October.** All D: processes died (per-cluster batch,
   7x7x5 local analysis); the D: phono3py environment, the hash-pinned inputs and the memory-gap event files became
   unavailable. Recovery: the replica interpreter `C:\Users\Koussay\ResearchLab\envs\aln-phono3py-4.5.0-replica`
   (same pinned versions) and a re-download of the pinned inputs from the Olympics repository (commit 0640f077),
   SHA-256 identical to the pins. The regenerated 5x5x3 event list, reference, engine check, basis and linearisation
   reproduce the D: numbers to every printed digit (kappa_xx 214.29993549566584, tau_x 34.252990795285065,
   G c0 = Phi0 to 1.29e-11). The 7x7x5 linearisation (D:, 23,313 rows, verified) is not recoverable without D:; the
   7x7x5 constrained analysis was not repeated (economy).

10. **PowerShell process-tree helper looped** (null children kept the frontier non-empty); stopped after the first
    poll; processes were then killed by explicit PID.

## Resumed 10 October 2026 (drive D: available again)

11. **First per-cluster batch (D:, unscaled coordinates) superseded.** The runs of 5 October on D:
    (results/m553/t3_orb/c*_orb*.json: per-cluster eps 0.2 -> tau_x <= 36.80 ps, eps 0.5 -> <= 46.15 ps) stopped
    after one L-BFGS-B iteration per outer step (item 8). They are feasible but far from local optima: the scaled
    run reaches 55.37 ps at eps 0.2. Kept for the record; not used.

12. **Sign bug in the first ALAMODE mapping (caught by an equivalence test).** The polarity map z -> -z multiplies
    each fc3 component by (-1)^(number of z indices). The first version wrote this as
    (-1) ** ((x1 == 3) + (xj == 3) + (xk == 3)) with NumPy integer indices: the comparisons are numpy.bool_, and
    numpy.bool_ + numpy.bool_ is a logical OR, so components with two z indices received -1 instead of +1. The
    symmetry-projection residual (5e-14) and the sum rule (7e-13) did NOT detect it: the wrongly signed tensor is
    still P6_3mc-symmetric, since R = diag(1,1,-1) normalises the point group. It was caught because the
    inversion map (r -> -r + c/2, factor -1 for every component), which must give the same tensor up to the 2_1
    screw, disagreed (Frobenius distance to the reference 16.5 % versus 2.4 %, nearest-neighbour basal orbits
    15 - 20 % versus 0.3 %). With the sign written as a product of per-axis factors both maps agree to 1.2e-16.
    Lesson: symmetry residuals do not validate improper-isometry conventions; use two maps that differ by a
    crystal symmetry.

13. **Calibrated class in unweighted coordinates stalled.** With per-orbit bounds eps_o = max(k delta_o, 0.005)
    (0.5 % ... 2,958 % at k = 2) imposed in the coordinates of item 8, the constraint curvature scales like
    1/eps_o^2 (range 3.5e7) and L-BFGS-B made zero steps from the second outer iteration on (cal2_c0_max: no
    nominally feasible iterate in 15 outer iterations; best infeasible tau_x 37.4 ps, lifetime violation 10.4 %).
    Batch stopped after 10 minutes (results/m553/t3_cal_superseded/). Fixed as in item 8 with eps-weighted
    reference norms (orbit_scaling.py --eps; cond(M) = 2.4e4 instead of 4.0e8), in which every constraint reads
    ||(W u)_o|| <= 1.

14. **The single-mesh witnesses of E15–E18 were sampling artefacts.** The calibrated maxima (2.35 / 2.15 at k = 1,
    3.72 / 5.95 at k = 2) and the per-cluster maxima were found with the data imposed on the 5x5x3 event set only.
    The same tensors on 7x7x5 give tau ratios 0.83 / 1.08 / 0.71 / 0.53 and violate kappa by up to 30 % and single
    lifetimes by up to +338 % (review section 4.2; reproduced with the campaign route, E19). The optimiser had
    tuned amplitudes at the triplets the 75-point mesh samples (48–71 orbit rates suppressed by > 10x). The final
    report of 10 October drew a physical conclusion ("tau_mem pinned only to a factor ~2") from these witnesses;
    it is withdrawn. Lesson: a witness for a statement about a material must be validated on a sampling not used
    in its construction.

15. **Two-mesh constraints do not cure it (E20–E21).** Imposing the data on 5x5x3 and 7x7x5 simultaneously still
    allows tau ratios 1.80 / 3.05 (calibrated) and 1.33 / 1.69 (per-cluster) on the objective mesh; on the unused
    9x9x6 mesh these fall to 1.18 / 0.97 / 1.18 / 1.14 with kappa violated by up to 8 %. The freedom moves to
    whatever the objective mesh samples. A finite list of meshes is not a substitute for a converged operator.

16. **Optimiser budget was binding in every run.** All inner L-BFGS-B solves (E15–E18, E20) stopped at maxiter; the
    review's continuations (E22) moved the per-cluster eps = 0.1 range from [0.84, 1.31] x [0.81, 1.37] to at least
    [0.83, 1.40] x [0.80, 1.44]. The multistart agreement of E16 (2.6 % / 0.5 %) was agreement between equally
    truncated runs, not evidence of a maximum.

17. **Cutoff zeros counted as inter-workflow differences.** On 266 of 456 cluster orbits the folded ALAMODE tensor
    is identically zero because of ALAMODE's cubic cutoff (largest cluster 5.68 A), so delta_o = 1 there measured
    the reference norm, not a disagreement; these orbits carried 45 % of the k = 1 box (review 1.2). E13's
    statement that cutoff zeros occur only beyond 6 A was inaccurate (126 cutoff-zero orbits lie below 6 A by the
    campaign's shell labels; the labels are not cluster diameters for 60 orbits).

18. **Over-generalised "artefact" conclusion (second report, 11 October).** The second report stated that every
    single-mesh witness loses or reverses its memory change on a finer mesh and that no witness changes tau by more
    than 1.18 x on an unused mesh. Its own E19 table contained the counter-example (E15 orb 0.1 max z: 1.370 on
    5x5x3, 1.347 on 7x7x5), the section-1 summary quoted only the calibrated rows, and the one surviving witness was
    not transferred to 9x9x6. The audit's transfer gives 1.518 (and 1.599 for the review continuation). Superseded
    by the third report. Lesson: state what every row of a table shows, and transfer the witness that contradicts
    the hypothesis first.

19. **A converged bulk tau_mem was presupposed.** Both earlier reports treated the mesh dependence of tau_mem as a
    discretisation error to be removed by a convergence study. The memory moment N = |C^+ b|^2 is carried by the
    Gamma-adjacent acoustic modes; it grows without a convergent regime through 12x12x8 (event operator) and doubles
    from 15x15x9 to 31x31x17 (RTA, tetrahedron) while kappa changes by 11 % (E23, E24); a long-wavelength linewidth
    law gamma ~ f^n with n >= 3/2 makes it infrared-divergent. The identifiability question was posed for a quantity
    whose existence in bulk is unresolved.
