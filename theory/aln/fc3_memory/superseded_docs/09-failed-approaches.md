# 09 — Failed approaches, bugs and negative results

1. **phono3py's pure-Python real-to-reciprocal prototype disagrees with its production kernels.**
   `phono3py.phonon3.real_to_reciprocal.RealToReciprocal` (4.5.0) builds the pre-phase from
   `np.where(np.all(triplet != 0, axis=0), triplet.sum(axis=0), 0)`, i.e. a Cartesian component of
   G = q0 + q1 + q2 is dropped unless all three addresses are non-zero in that component. Re-implementing that
   rule, 1,062 of 41,131 allowed events (38 of 1,613 triplets, all Umklapp with G = (0, +-5, 0)) disagreed with the
   Rust interaction strengths by up to a factor 440 (t5_engine_check, first run). Using the full G for every
   component reproduces the Rust kernel for every event (max relative difference 6.5e-12 with r0 averaging,
   2.5e-11 without). The production C/Rust path is the one used by the validated pipeline; the Python prototype
   (lang="Python") should not be used for Umklapp-sensitive work.

2. **Orbit averaging of per-event rates is not gauge invariant.** The memory-gap reference ties rates by averaging
   over the 24-operation orbit. For events touching degenerate bands this keeps a gauge-dependent split between
   degenerate partners: 128 of 2,292 allowed orbits differ from the gauge-invariant block average by up to a factor
   18. Effect on the reference: kappa_zz 201.967 -> 201.904 W/(m K) (-3e-4), tau_z 27.4465 -> 27.4410 ps (-2e-4);
   x unchanged to 1e-9. More visible for lifetimes: for the full-space witness (tau_x = 481 ps) the orbit-averaged
   phono3py-route lifetime of one doubly degenerate mode (19.97 THz) moves by 58 % while the degeneracy-averaged
   lifetimes move by at most 9.98 % (the imposed tolerance). Per-mode lifetimes of degenerate partners in the
   orbit-average convention are gauge artefacts; constraints must be imposed on degeneracy-averaged quantities.
   (This also applies to the memory-gap statement "all 897 lifetimes fixed": 22 modes are degenerate partners.)

3. **First symfc projection attempt failed (shape mismatch):** the dataset fc3.hdf5 is stored in full form
   (72 x 72 x 72 x 3 x 3 x 3), not compact; the compact tensor is fc3[p2s]. Fixed in probe_symfc.py.

4. **Per-event outer ceiling is vacuous beyond a 0.5 % fc3 change.** The rigorous bound
   tau <= 1/lambda_min(C(g_low(rho))) with per-amplitude floors |Phi| >= |Phi0| - ||G_m|| rho |c0| gives 81.5 ps
   (rho = 0), 481 ps (0.1 %), 3,287 ps (0.2 %), infinity from 0.5 % on: individual amplitudes can be silenced by
   tiny worst-case fc3 changes (fragility rho* = |Phi0|/(||G_m|| |c0|): median 0.26 %, max 2.2 %). The bound treats
   all events as simultaneously worst case and is far from the attainable values; it is reported, not used.

5. **The Frobenius norm is a misleading plausibility measure for fc3 changes.** The first full-space witness
   (tau_x = 481 ps at a 9.7 % Frobenius change) leaves the nearest-neighbour clusters (68 % of ||fc3_0||^2) almost
   unchanged (1.7 %) and rewrites the clusters beyond 3.5 A (0.2 % of ||fc3_0||^2) by factors 1.1 - 5.0, where 87 %
   of the change sits. A small global norm can hide order-one relative changes of the small long-range constants.
   Consequence: the physically meaningful classes are (i) short-range subspaces (cutoff) and (ii) per-cluster-orbit
   relative bounds; the Frobenius-bounded results are kept as the mathematical analogue only.

6. **Code bug (caught before use):** a patch inserted two methods in the middle of FC3Model.__init__ in
   optimize_fc3.py, truncating the constructor (AttributeError at the first evaluation). The first launch of the
   per-orbit batch failed immediately with no output; the file was repaired and the batch relaunched. Batches
   started before the patch had already imported the correct module.

7. **Augmented-Lagrangian convergence is slow near active lifetime constraints:** mu escalates to 1e7 before the
   violation falls below 1e-7 (about 15 min for one full-space run). Internal margins (0.998 dr, 0.98 dK) make the
   reported iterates satisfy the nominal tolerances exactly; only nominally feasible iterates are reported.
