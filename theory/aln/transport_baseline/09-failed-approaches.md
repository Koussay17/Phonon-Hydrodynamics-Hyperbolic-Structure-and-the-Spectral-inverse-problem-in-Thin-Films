# 09 — Failed approaches and negative results (kept on purpose)

1. **phono3py 4.5.0 defaults do not reproduce the Olympics phono3py numbers.** At 9x9x5,
   300 K, RTA = 231.18/212.42 W/(m K) vs the team's 224.93/206.23 (+2.8 %/+3.0 %).
   Cause isolated by single-switch runs (`runs/t1-rta-m995-*`): the default
   `make_r0_average=True` (fc3 real-to-reciprocal phase averaged over the three atoms,
   introduced after v2.x). With `make_r0_average=False` the difference falls to
   +8e-5/+1.3e-4. Finite-difference vs analytic group velocities changes kappa by 5e-9.
   Not an error of either code: the two conventions differ because fc3 is truncated by the
   3x3x2 supercell; the 2.8–3 % spread is a measure of that truncation.
2. **C backend + analytic NAC group velocities**: `NotImplementedError` (the C path has no
   Gonze-Lee derivative). The C backend must be run with `gv_delta_q=1e-5`
   (`runs/logs/FAILED-t1-rta-m995-lang-C-analyticNACgv.log`).
3. **C vs Rust kernels are not bit-equivalent.** kappa differs by 1.3e-4 (xx)/3.5e-4 (zz) at
   9x9x5 and 2e-5/1.2e-4 at 15x15x9. Localized (`scripts/diag_c_vs_rust.py`,
   `results/diag_c_vs_rust_m995.json`): integration weights and interaction strengths are
   identical given identical phonons (except q0 = Gamma, which does not enter kappa); the
   difference comes from polar phonon frequencies differing by <= 5.7e-7 THz between the two
   Gonze-Lee implementations (1.2e-7 THz even with a doubled reciprocal cutoff,
   `results/nac_backend_check_m995.json`). Swapping only the phonon set changes individual
   tetrahedron linewidths by up to 1.6 % on 9x9x5: the tetrahedron linewidth is
   ill-conditioned with respect to tiny frequency perturbations on coarse meshes.
4. **Python reference tetrahedron weights** (`lang="Python"`) differ from C/Rust by up to
   0.28 (absolute, g_max ~2–6) at some grid points; C and Rust agree exactly. Not pursued
   (the Python path is not used by production code).
5. **Plain Jacobi-PCG on the irreducible matrix diverged** (kappa_xx = 4.8e4 W/(m K) at 9x9x5,
   2000 iterations, residual 8e-5) because symmetry-forbidden directions (stabilizer
   subspaces not aligned with x, y, z) are exact null directions that a diagonal criterion
   does not detect, and finite-difference velocities have ~1e-9 components there. Fixed by
   projecting the right-hand side, residual and preconditioned residual onto
   range(P x S_i) (degeneracy averaging x stabilizer average): 19–20 iterations, kappa equal
   to phono3py's pseudo-inverse to 4e-13 (9x9x5) and 1.2e-12 (15x15x9).
6. **Dense phono3py LBTE at 31x31x17 is infeasible on this machine**: 7.2 GiB matrix
   (n = 31104) with ~3–4 GB free RAM and ~4–5 GB free commit; scipy `dsyev` took 183 s at
   n = 4860, i.e. ~13 h by n^3 scaling. Replaced by the out-of-core route (rows on disk +
   own assembly + PCG), validated against the dense solver.
7. **Mixing `tetra` and Gaussian widths in one phono3py run** silently reuses interaction
   strengths computed with tetrahedron zero-skipping for the Gaussian widths
   (`LBTECollisionSolver._run_interaction`); the driver refuses such runs.
8. **Environment interruptions**: the external drive D: disappeared during the campaign
   (02:40–11:45 on 3 Oct); work continued on a C: staging copy with a replica interpreter
   (bit-identical kappa on the 9x9x5 test) and hash-verified inputs, then everything was
   copied back (cmp-verified). The Windows memory query initially returned `{}` (wrong API
   binding); fixed (`K32GetProcessMemoryInfo`).
9. **Energy is not a null vector of phono3py's matrix** (expected, structural): see
   07-experiments.md §5. Not a numerical failure, but it means the exported matrix cannot be
   used as an even-sector (energy, viscosity) operator.
