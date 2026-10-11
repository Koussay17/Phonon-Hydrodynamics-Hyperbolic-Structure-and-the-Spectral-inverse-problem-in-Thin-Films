# 07 — Experiments (fc3-parametrised rate class, AlN, 300 K)

All scripts in `scripts/` (copies of the memory-gap helpers aln_common, aln_events, aln_geometry, memgap, fastops,
symmetry; nothing imported from the thesis tree). Results in `results/m553/` (5x5x3) and `results/m775/` (7x7x5).
Environments: phono3py env `D:\ResearchLab\envs\aln-phono3py-4.5.0` (phono3py/phonopy 4.5.0, symfc 1.7.3, Rust
kernels) for anything calling phono3py or symfc; research env for analysis and optimisation. Inputs hash-checked
at every phono3py run (aln_common.verify_inputs). Event lists: memory-gap exports (5x5x3:
`20261004-memory-gap/results/aln/events_m553_s0.1.npz`; 7x7x5: `.../review/cx_runs/events_m775_s0.1.npz`).
Units: tau in ps, kappa in W/(m K), lambda in THz (collision units).

## E1 — Reference from fc3 (task 1)  [p3_events.py, t1_reference.py]
phono3py's interaction strengths pp recomputed on the fixed event list (each event at phono3py's representative
triplet, 1.2 s per pass at 5x5x3, 21 s at 7x7x5); W = conv g0w /(4 s_p s_a s_b) (or /8 s_p s_a^2).
* W * pp reproduces the memory-gap event rates: max relative difference 1.3e-15 (5x5x3), 1.4e-15 (7x7x5).
* 5x5x3, orbit-averaged (memory-gap convention): kappa_xx 214.2999, kappa_zz 201.9670, tau_x 34.2530,
  tau_z 27.4465, tau_CS 19.1903 / 15.5168, lambda_min(C|H) 9.7643e-4 (F = 1 ceiling 81.499 ps). Identical to the
  memory-gap reference.
* 5x5x3, gauge-invariant block average (base point of task 3): kappa 214.2999 / 201.9039, tau 34.2530 / 27.4410.
* 7x7x5: kappa 243.027 / 272.334, tau 52.609 / 50.623, tau_CS 21.103 / 20.016, lambda_min 4.0353e-4 (ceiling
  197.2 ps). Identical to the memory-gap 7x7x5 values (52.6 / 50.6 ps).
  Files: results/m553/t1_reference_m553.json, results/m775/t1_reference_m775.json.

## E2 — Independent amplitude engine versus phono3py (task 5, amplitudes)  [amp_engine.py, t5_engine_check.py]
NumPy re-implementation of real-to-reciprocal (both conventions, all_shortest shortcut) and reciprocal-to-normal.
* All 41,131 allowed events, reference fc3: max relative difference 6.5e-12 (r0 average) and 2.5e-11 (single leg);
  max |difference| / max pp = 2.5e-15 over all 44,406 events.
* Witness fc3 (E9): 4.1e-10. Random symmetric fc3 (iso_0.10_s0, E8 check): model and phono3py agree (below).
* Finding: phono3py's pure-Python prototype uses a different pre-phase for some Umklapp triplets (09, item 1).

## E3 — Phase convention make_r0_average (task 2a)  [t1_reference.py]
| mesh | kappa_xx | kappa_zz | tau_x | tau_z | tau_CS x / z | lambda_min | lifetimes |
|---|---|---|---|---|---|---|---|
| 5x5x3 True | 214.30 | 201.97 | 34.253 | 27.446 | 19.19 / 15.52 | 9.764e-4 | — |
| 5x5x3 False | 209.95 (-2.0 %) | 197.94 (-2.0 %) | 33.080 (-3.4 %) | 26.211 (-4.5 %) | 18.80 / 15.21 | 1.0068e-3 | -13.5 % ... +8.8 % (rms 2.9 %) |
| 7x7x5 True | 243.03 | 272.33 | 52.609 | 50.623 | 21.10 / 20.02 | 4.035e-4 | — |
| 7x7x5 False | 238.44 (-1.9 %) | 264.60 (-2.8 %) | 50.503 (-4.0 %) | 45.514 (-10.1 %) | 20.70 / 19.45 | 4.108e-4 | -14.2 % ... +14.1 % (rms 2.0 %) |
Per-event rates change by factors 0.02 - 76 (0.1 % / 99.9 % quantiles, 5x5x3) although tau moves by 3 - 10 %.

## E4 — Real-space truncation (task 2b)  [variants_p3.py, t2_variants.py]
cut_zero_r: phono3py cutoff_fc3_by_zero (ASR broken). cut_sym_r: cut_zero_r projected onto the symfc cutoff-r
symmetric + ASR subspace (basis sizes 218 / 569 / 1,282 / 2,098 at 3.5 / 4.0 / 4.5 / 5.0 A). 3.2 and 3.5 A give
identical tensors (no cluster in between). eps = relative Frobenius distance to fc3_0.
5x5x3 (tau ratios to reference; dr = max |lifetime change|):
| variant | eps | kappa x / z | tau_x | tau_z | tau/tau_CS x / z | dr |
|---|---|---|---|---|---|---|
| cut_sym 2.0 | 0.097 | 1.094 / 1.468 | 47.35 (1.382) | 64.30 (2.343) | 2.26 / 2.82 | 0.75 |
| cut_sym 3.2 | 0.044 | 1.010 / 1.018 | 38.04 (1.111) | 29.23 (1.065) | 1.96 / 1.85 | 0.25 |
| cut_sym 4.0 | 0.036 | 0.969 / 0.974 | 32.77 (0.957) | 27.04 (0.985) | 1.76 / 1.79 | 0.45 |
| cut_sym 4.5 | 0.031 | 1.000 / 1.009 | 33.46 (0.977) | 28.09 (1.023) | 1.74 / 1.80 | 0.30 |
| cut_sym 5.0 | 0.022 | 0.993 / 0.999 | 33.33 (0.973) | 26.95 (0.982) | 1.75 / 1.74 | 0.19 |
| cut_zero 4.0 | 0.036 | 0.960 / 0.901 | 32.24 (0.941) | 21.70 (0.791) | 1.75 / 1.55 | 1.07 |
7x7x5:
| variant | eps | kappa x / z | tau_x | tau_z | tau/tau_CS x / z | dr |
|---|---|---|---|---|---|---|
| cut_sym 2.0 | 0.097 | 1.262 / 1.422 | 105.57 (2.007) | 173.03 (3.418) | 3.96 / 6.08 | 0.85 |
| cut_sym 3.2 | 0.044 | 1.028 / 1.074 | 58.87 (1.119) | 72.04 (1.423) | 2.71 / 3.35 | 0.33 |
| cut_sym 4.0 | 0.036 | 0.957 / 0.964 | 47.11 (0.895) | 48.51 (0.958) | 2.33 / 2.52 | 0.38 |
| cut_sym 4.5 | 0.031 | 0.980 / 0.975 | 49.98 (0.950) | 49.71 (0.982) | 2.42 / 2.55 | 0.20 |
| cut_sym 5.0 | 0.022 | 0.997 / 0.997 | 50.99 (0.969) | 51.16 (1.011) | 2.42 / 2.56 | 0.21 |
| cut_zero 3.2 | 0.044 | 0.973 / 1.005 | 42.63 (0.810) | 50.19 (0.991) | 2.08 / 2.50 | 2.69 |
Reference tau/tau_CS: 1.785 / 1.769 (5x5x3), 2.493 / 2.529 (7x7x5). Without re-imposing the sum rule, single
lifetimes change by factors up to 4 at 7x7x5.

## E5 — Random symmetry-preserving perturbations (task 2c)  [variants_p3.py, t2_variants.py]
iso: uniform direction in the 7,752-dim symmetric + ASR space; mul: P_sym(xi * fc3_0), xi iid N(0,1) per element;
both rescaled to relative Frobenius size eps. 5 seeds (5x5x3), 3 seeds (7x7x5). Check: (1.05) fc3_0 gives
tau and kappa x 1/1.1025 exactly (31.068 = 34.253/1.1025; 47.718 = 52.609/1.1025).
5x5x3, mean +- std [min, max] of tau / tau_ref; kappa ratio mean; max lifetime change range:
| family eps | tau_x ratio | tau_z ratio | kappa x / z | max |dr| |
|---|---|---|---|---|
| iso 1 % | 0.990 +- 0.004 [0.985, 0.994] | 1.019 +- 0.036 [0.987, 1.071] | 0.992 / 1.000 | 0.06 - 0.15 |
| iso 5 % | 0.760 +- 0.062 [0.682, 0.827] | 0.817 +- 0.070 [0.754, 0.936] | 0.872 / 0.855 | 0.71 - 1.74 |
| iso 10 % | 0.542 +- 0.043 | 0.544 +- 0.034 | 0.666 / 0.644 | 3.1 - 10.6 |
| iso 20 % | 0.295 +- 0.031 | 0.279 +- 0.034 | 0.382 / 0.357 | 7.8 - 28.7 |
| mul 1 % | 0.983 +- 0.040 [0.931, 1.027] | 0.978 +- 0.014 | 0.993 / 0.986 | 0.12 - 0.47 |
| mul 5 % | 0.827 +- 0.039 | 0.804 +- 0.049 | 0.896 / 0.877 | 0.79 - 1.53 |
| mul 10 % | 0.617 +- 0.063 | 0.516 +- 0.092 | 0.667 / 0.604 | 1.8 - 6.1 |
| mul 20 % | 0.323 +- 0.114 | 0.301 +- 0.173 | 0.342 / 0.269 | 6.2 - 57 |
7x7x5:
| family eps | tau_x ratio | tau_z ratio | kappa x / z | max |dr| |
|---|---|---|---|---|
| iso 1 % | 0.947 +- 0.019 | 0.983 +- 0.083 [0.934, 1.078] | 0.987 / 0.991 | 0.11 - 0.12 |
| iso 5 % | 0.653 +- 0.029 | 0.576 +- 0.050 | 0.858 / 0.797 | 1.1 - 2.5 |
| iso 10 % | 0.385 +- 0.008 | 0.323 +- 0.032 | 0.652 / 0.561 | 6.7 - 10.9 |
| iso 20 % | 0.204 +- 0.007 | 0.171 +- 0.007 | 0.380 / 0.317 | 17 - 30 |
| mul 1 % | 0.965 +- 0.036 | 0.971 +- 0.016 | 0.989 / 0.989 | 0.10 - 0.31 |
| mul 5 % | 0.655 +- 0.042 | 0.656 +- 0.075 | 0.862 / 0.834 | 0.76 - 1.85 |
| mul 10 % | 0.556 +- 0.073 | 0.477 +- 0.049 | 0.669 / 0.594 | 2.8 - 4.6 |
| mul 20 % | 0.366 +- 0.108 | 0.382 +- 0.109 | 0.421 / 0.360 | 5.2 - 40 |
Random perturbations add scattering on average (E|Phi0 + dPhi|^2 = |Phi0|^2 + E|dPhi|^2): tau and kappa fall;
tau/tau_CS falls from 1.78 towards 1.3 - 1.4. Per-mode lifetimes are the most sensitive observable (10 - 47 %
already at eps = 1 %). None of these variants is data-consistent beyond eps = 1 %.

## E6 — Exact linearisation, rank and fragility (5x5x3)  [linearize.py]
Rows: 2,292 allowed orbits (335 selection-rule orbits stay zero), 2,712 amplitude rows (block sizes 1/2/4/8:
1,948/318/20/6 orbits). Checks: L f0 = Phi0 to 7.6e-14; G c0 = Phi0 to 1.3e-11; linear model versus engine on an
explicit perturbed fc3 7.0e-15; finite-difference gradients (dlog tau, dK, dr) to 1e-7 relative (check_model.py).
Model versus phono3py (all events, orbit average) for iso_0.10_s0: tau_x 19.856483 vs 19.856483, tau_z 14.83220
vs 14.83183, kappa_z 132.617 vs 132.629 (difference = block versus orbit averaging).
* Rate Jacobian (2,292 x 7,752, relative): rank 2,113; the 179 null directions are the orbits tied by shared
  degenerate blocks (159 + 15 + 5). **Every independent block rate can be moved independently, at first order, by
  a symmetric fc3 perturbation.** Real rank of the amplitude map: 3,758 of 5,424.
* Fragility rho* = |Phi0| / (||G_m|| |c0|) (relative Frobenius change that silences amplitude m in its worst
  direction): min 6.6e-6, 1 % quantile 2.7e-5, median 2.6e-3, max 0.0215.

## E7 — Rigorous ceiling for Frobenius-bounded fc3 changes  [t3_ceiling.py]
tau <= 1/lambda_min(C(g_low(rho))|H), any symmetric fc3 with eps <= rho, no data needed:
rho = 0: 81.49 ps; 0.1 %: 480.6 ps (393 orbits individually silenceable); 0.2 %: 3,287 ps (850); >= 0.5 %:
vacuous (1,798 orbits silenceable at 0.5 %, all 2,292 at 5 %).

## E8 — Cluster structure of the symmetric basis  [inline check, results/m553/cluster_orbits.npz]
Each of the 8,098 symfc columns lies in one distance shell (spread of the largest pair distance 4e-14); the columns
form 456 symmetry-distinct cluster orbits (up to 27 columns each). Reference norm per orbit spans 0.058 - 1,166.
Shells (largest pair distance): on-site 2 orbits (31 % of ||fc3_0||^2), < 2 A 4 orbits (68 %), 2 - 3.2 A 24
orbits (0.7 %), 3.5 - 5 A 102 orbits (0.15 %), > 5 A 324 orbits with 5,780 columns (0.05 %).

## E9 — First-order data-preserving slope  [local_analysis.py]
Data Jacobian (117 lifetimes + kappa_xx, kappa_zz; relative) with respect to y = c/|c0|:
| space | n | rank | |grad log tau_x| | s_x | |grad log tau_z| | s_z | null-space fraction x / z |
|---|---|---|---|---|---|---|---|
| full | 7,752 | 112 | 156.4 | 21.4 | 307.5 | 27.4 | 0.137 / 0.089 |
| cutoff 5.0 A | 2,098 | 112 | 95.7 | 13.2 | 176.0 | 17.3 | 0.138 / 0.099 |
| cutoff 3.5 A | 218 | 112 | 33.7 | 2.52 | 53.6 | 3.33 | 0.075 / 0.062 |
s = |P_null grad log tau|: with the data held fixed to first order, tau ~ tau_0 exp(+- s eps). The static data
see 86 - 94 % of the tau sensitivity; the remainder is dominated by the long-range constants.

## E10 — Constrained range, full symmetric space (task 3)  [optimize_fc3.py, t3_summary.py]
Lifetimes of all 117 mode orbits within +-10 %, kappa_xx and kappa_zz within +-1 %, block-average base point
(tau 34.2530 / 27.4410 ps). Rigorous floor: tau >= K/|b|^2 >= (1 - dK) tau_CS = 19.00 ps (x), 15.36 ps (z).
| target | Frobenius bound | value found (ps) | ratio to base | eps used | min / max orbit-rate ratio | orbits suppressed 10x |
|---|---|---|---|---|---|---|
| max tau_x | rho <= 0.1 (t3_test) | 481.45 | 14.06 | 0.0967 | 1.5e-9 / 1.7e4 | 263 |
| max tau_x | none (t3_main) | 453.54 | 13.24 | 0.0962 | 4.8e-12 / 2.6e4 | 245 |
| min tau_x | none | 19.84 | 0.579 | 0.0826 | 3.0e-9 / 7.3e3 | 180 |
| max tau_z | none | 17,121.27 | 623.9 | 0.0956 | 5.2e-10 / 4.4e3 | 324 |
| min tau_z | none | 15.40 | 0.561 | 0.0852 | 1.8e-6 / 5.0e4 | 75 |
The unbounded max tau_x run stopped below the rho <= 0.1 run (a subset of its class): local optimiser (09, item 6).
The minima lie within 4.4 % (x) and 0.3 % (z) of the rigorous floor. Every witness reaches the lifetime tolerance
(max deviation 0.0998) and suppresses hundreds of event orbits by more than a factor 10.

## E11 — Constrained range, short-range subspace (cutoff 3.5 A, 218 parameters)  [t3_cut35]
Same constraints; the subspace contains the on-site, < 2 A and 2 - 3.2 A clusters (99.85 % of ||fc3_0||^2).
| rho | max tau_x | min tau_x | max tau_z | min tau_z |
|---|---|---|---|---|
| 0.05 | no feasible iterate | 26.86 | 43.78 | 20.86 |
| 0.1 | 74.02 | 25.22 | 57.27 | 19.41 |
| 0.2 | 91.16 | 24.06 | 81.98 | 18.34 |
| 0.5 | 97.35 (eps 0.29) | 24.02 | 89.73 (eps 0.29) | 18.30 |
| none | 97.37 (eps 0.29) | 24.01 | 91.27 (eps 0.31) | 18.36 |
Unbounded within the subspace: tau_x in [24.0, 97.4] ps (0.70 - 2.84 x base), tau_z in [18.3, 91.3] ps
(0.67 - 3.33 x); the unbounded optima use Frobenius changes of 29 - 31 %.

## E12 — Witness verification, full-space witnesses (task 5)  [witness_fc3.py, verify_witness.py]
fc3 tensors rebuilt from the optimiser coordinates; route B = phono3py 4.5.0 Rust interaction strengths for all
44,406 events, orbit-averaged rates, tau by dense eigendecomposition on H and by conjugate gradients; route C =
independent NumPy amplitude engine (E2), block average; route A = the linear optimisation model.
| witness | route A tau_x / tau_z | route B (phono3py) | route C (engine) | eig vs CG | engine vs Rust pp |
|---|---|---|---|---|---|
| max tau_x, rho <= 0.1 | 481.4494 / 35.2813 | 481.4494 / 35.2794 | 481.4494 / 35.2813 | 2e-13 | 4.1e-10 |
| max tau_z, unbounded | 42.5416 / 17,121.27 | 42.5416 / 17,120.62 | 42.5416 / 17,121.27 | 2e-12 | 1.3e-10 |
The z differences (4e-5) are the block versus orbit averaging of E1. lambda_min(C|H) falls to 5.27e-5 and
2.65e-6 THz (F = 1 ceilings 1,509 and 30,037 ps). Degeneracy-averaged lifetimes stay within 9.98 %; raw
orbit-average lifetimes of degenerate partners deviate by 58 % and 44 % (gauge artefact, 09 item 2).
Anatomy (09 item 5): on-site and < 2 A clusters change by 1.6 - 2.3 %, the 2 - 3.2 A shell by 34 %, the clusters
beyond 3.5 A by factors 1.1 - 5.0 shell-wise; the largest single cluster orbit of the tau_x witness changes by
42.9 times its reference norm.

## E13 — Calibration of fc3 plausibility with an independent workflow  [alamode_compare.py, alamode_transport.py]
Source: Phonon Olympics ALAMODE AlN cubic force constants, `Aluminum Nitride/ALAMODE/force_constants/3rd-order/
cubic_15nn_ndata200_ols.xml` at commit 0640f077 (ALM 1.4.2, 4x4x3 supercell, own DFT setup, a = 3.1135 A,
c = 4.9831 A versus 3.1130 / 4.9821 A for the phono3py workflow); SHA-256 2acff003...0f85d1; stored outside the
repository (D:\ResearchLab\datasets\phonon-olympics\0640f077...\AlN_alamode\). Mapping (script header): stored
(pair2, pair3) half completed by the swap; ALM image ordering (every sampled stored image is a shortest image);
polarity map z -> -z with factor (-1)^(#z), checked against the inversion map r -> -r + c/2 (agreement 1.2e-16,
see 09 item 12); 3x3x2 folding; Ry/bohr^3 -> eV/A^3. Folded tensor: sum rule 7.4e-13, symmetric-subspace residual
5.2e-14 (no ASR) and 7.9e-14 (with ASR).
Per-cluster-orbit relative difference delta_o = ||(eta_ALAMODE - eta_0)_o|| / ||eta_0,o||:
| shell (A) | orbits | share of ||fc3_0||^2 | shell rel. difference | orbit median | orbit max | norm ratio ALAMODE/ref |
|---|---|---|---|---|---|---|
| on-site | 2 | 30.7 % | 0.34 % | 0.31 % | 0.46 % | 0.997 |
| < 2 | 4 | 68.4 % | 0.36 % | 0.36 % | 0.56 % | 0.997 |
| 2 - 3.2 | 24 | 0.68 % | 7.6 % | 16 % | 100 % | 1.003 |
| 3.2 - 4 | 20 | 0.064 % | 27 % | 35 % | 103 % | 0.988 |
| 4 - 5 | 82 | 0.082 % | 42 % | 63 % | 447 % | 0.977 |
| 5 - 6 | 184 | 0.040 % | 80 % | 100 % | 1,479 % | 0.879 |
| > 6 | 140 | 0.006 % | 100 % (ALAMODE zero: 15th-neighbour cutoff) | 100 % | 100 % | 0 |
Frobenius distance 2.36 %; norm ratio 0.9972. The dominant clusters agree to better than 0.6 %; the small
long-range ones differ by order one.
Hybrid fc2 (phono3py) + fc3 (ALAMODE, folded, ASR-projected), same event list, orbit-average convention:
kappa_xx 214.53 (+0.11 %), kappa_zz 208.62 (+3.29 %); tau_x 33.996 ps (-0.75 %), tau_z 27.906 ps (+1.67 %);
tau_CS +0.11 % / +3.29 %; eig and CG agree to 2e-15; lambda_min 9.08e-4 THz. Degeneracy-averaged lifetimes:
-18.3 % ... +17.5 %, rms 4.9 %, 97.2 % of modes within 10 %. Event-rate ratios ALAMODE/reference: 0.011, 0.30,
0.98, 3.6, 43 at the 0.1, 1, 50, 99, 99.9 % quantiles. Selection-rule events stay zero (<= 2e-15 of the maximum).
A real second workflow thus moves single event rates by up to two orders of magnitude and lifetimes by up to 18 %,
but tau_mem by less than 2 %.

## E14 — First-order data-preserving slope at 7x7x5 (task 4)  [local_analysis.py, results/m775/local_analysis.json]
Data: 285 mode-orbit lifetimes + kappa_xx, kappa_zz (287 rows); linearisation of 5 October (23,313 amplitude rows,
G c0 = Phi0 to 5.1e-11), recovered from D:.
| space | n | rank | sv min / median / max | |grad log tau_x| | s_x | |grad log tau_z| | s_z | null fraction x / z |
|---|---|---|---|---|---|---|---|---|
| full | 7,752 | 276 | 3e-13 / 62.2 / 1,313 | 198.4 | 22.5 | 687.6 | 33.5 | 0.114 / 0.049 |
| cutoff 5.0 A | 2,098 | 276 | 2e-13 / 39.1 / 708 | 128.6 | 12.3 | 359.4 | 16.5 | 0.095 / 0.046 |
| cutoff 3.5 A | 218 | **218** | **0.118** / 11.1 / 223 | 45.6 | 7.5e-14 | 78.3 | 1.4e-13 | 1.7e-15 / 1.7e-15 |
Compared with 5x5x3 (E9): the full-space and 5 A slopes are mesh-stable (21.4 / 27.4 and 13.2 / 17.3 there).
In the 3.5 A subspace the data Jacobian has full column rank at 7x7x5 (rank 112 of 218 at 5x5x3, with 119 data):
the static data identify the short-range fc3 to first order, and tau has no first-order freedom at fixed data.
Conditioning: the smallest singular value 0.118 means a 1 % data change can move the short-range coordinates by up
to ~8.5 % (relative Frobenius). Local statement only (first order at the reference); no 7x7x5 optimisation was run.

## E15 — Per-cluster-orbit bounded class (task 3)  [orbit_scaling.py, t3_batch.py, t3_batch_min.py]
Every one of the 456 cluster orbits changes by at most eps relative to its reference norm; lifetimes +-10 %,
kappa +-1 %, no Frobenius bound; scaled coordinates (09 item 8). Runs of 5 October (C: mirror) and, for the four
minima at eps = 0.05 and 0.2, of 10 October. tau in ps (ratio to the block-average base 34.253 / 27.441):
| eps | min tau_x | max tau_x | min tau_z | max tau_z |
|---|---|---|---|---|
| 0.05 | 30.12 (0.88) | 39.92 (1.17) | 23.56 (0.86) | 33.73 (1.23) |
| 0.10 | 28.84 (0.84) | 44.80 (1.31)* | 22.35 (0.81) | 37.58 (1.37) |
| 0.20 | 27.05 (0.79) | 55.37 (1.62) | 20.85 (0.76) | 45.41 (1.65) |
(*) best of the multistart E16. Ranges are nested in eps as they must be. All witnesses saturate the lifetime
tolerance (0.0998) and the bound of at least one cluster orbit; they suppress 10 - 51 orbit rates by more than a
factor 10 (minimum rate ratio 1e-4 or below) although every cluster moves by at most eps: individual amplitudes
are fragile (E6), tau is not.
Verification (route B = phono3py Rust amplitudes, orbit average; route A = model): eps 0.1 max x 44.0793 / 44.0793,
max z 37.5909 / 37.5821; eps 0.2 max x 55.3689 / 55.3689, max z 45.4271 / 45.4135; eig and CG identical to the
printed digits; degeneracy-averaged lifetimes within 9.98 %. Orbit-average convention: kappa_z of the eps 0.2 z
witness is 1.01004 (block convention 1.0098): a convention effect of 4e-5, not a constraint violation.

## E16 — Multistart for the eps = 0.1 maxima  [t3_multistart.py]
Starts: zero; 'aim' = the Frobenius-class witness of E10 (481 ps / 17,121 ps) clipped orbit by orbit to 0.9 eps (its
largest orbit change is 42.9 times the orbit's reference norm); two random points in the box.
| target | zero | aim | rand1 | rand2 | spread |
|---|---|---|---|---|---|
| max tau_x | 44.08 | 44.79 | 43.68 | 44.80 | 2.6 % |
| max tau_z | 37.58 | 37.40 | 37.43 | 37.39 | 0.5 % |
The large-tau mechanism of the Frobenius class does not survive the per-cluster bound: the clipped witness
direction returns to the same maximum. Still inner bounds (no global certificate).

## E17 — Calibrated class, k = 2  [orbit_scaling.py --eps, t3_calibrated.py]
eps_o = max(2 delta_o, 0.005), delta_o from E13 (on-site and nearest-neighbour orbits: 0.6 - 1.1 %; 2 - 3.2 A median
32 %; beyond 6 A: 200 %; largest 2,958 %). The ALAMODE tensor itself lies in the k = 1 class. Weighted coordinates
(cond(M) = 2.4e4; first attempt stalled, 09 item 13).
| target | tau found (ps) | ratio | other axis | Frobenius size | route B (phono3py) |
|---|---|---|---|---|---|
| max tau_x | 127.46 | 3.72 | tau_z 31.13 | 4.4 % | 127.4634 / 31.1313 |
| max tau_z | 163.42 | 5.96 | tau_x 37.16 | 4.4 % | 37.1576 / 163.3993 |
| min tau_x | 21.80 | 0.64 | tau_z 32.86 | 4.1 % | 21.7990 / 32.8683 |
| min tau_z | 16.74 | 0.61 | tau_x 34.01 | 4.0 % | 34.0109 / 16.7412 |
All four saturate their largest per-orbit bound (ratio 0.999) and the lifetime tolerance. Anatomy (every witness
alike): on-site and < 2 A shells change by 0.6 - 0.7 %; 2 - 3.2 A by 14 %; 3.2 - 4 A by 51 - 53 %; 4 - 5 A by
76 - 81 %; 5 - 6 A by 132 - 146 %; beyond 6 A by 133 - 166 %. The memory time is carried by the small long-range
anharmonic constants (0.19 % of ||fc3_0||^2 beyond 3.2 A), precisely the ones on which two independent workflows
disagree. Minima lie within 15 % (x) and 9 % (z) of the rigorous floor (1 - dK) tau_CS.
Orbit-average convention: kappa_z ratio 1.01057 for the min tau_x witness and 0.98995 for the min tau_z witness
(block convention within 1 %): convention effects of 6e-4 and 5e-5.

## E18 — Calibrated class, k = 1 (bounds = the observed inter-workflow difference)  [t3_calibrated.py]
| target | tau found (ps) | ratio | other axis | Frobenius size | route B (phono3py) |
|---|---|---|---|---|---|
| max tau_x | 80.33 | 2.35 | tau_z 30.00 | 2.2 % | 80.3258 / 30.0010 |
| max tau_z | 59.04 | 2.15 | tau_x 34.05 | 2.2 % | 34.0546 / 59.0475 |
| min tau_x | 23.89 | 0.70 | tau_z 29.01 | 2.1 % | 23.8942 / 29.0183 |
| min tau_z | 18.77 | 0.68 | tau_x 35.27 | 2.1 % | 35.2741 / 18.7697 |
Anatomy: on-site and < 2 A 0.3 - 0.5 %; 2 - 3.2 A 7 %; 3.2 - 4 A 25 - 27 %; 4 - 5 A 38 - 40 %; 5 - 6 A 67 - 73 %;
beyond 6 A 66 - 83 %. 48 - 71 orbit rates suppressed by more than a factor 10. These witnesses are no farther from
the reference (Frobenius 2.1 - 2.2 %) than the ALAMODE tensor is (2.36 %), keep every lifetime within 10 % and both
kappa within 1 % (ALAMODE itself: lifetimes up to 18 %, kappa_zz +3.3 %), yet span tau in [0.68, 2.35] x the
reference. Orbit-average kappa ratios within 1.0099 / 0.99006 (convention effects <= 6e-5).

<!-- E19-E22 written on the C: mirror while D: was unavailable (10-11 October 2026), merged here. -->

Context: the counterexample review of 10 October (`review/counterexample-review.md`, D:) reported that the 5x5x3
witnesses of E15–E18 lose their memory effect and violate the static data on the 7x7x5 event set. E19 reproduces
that finding with the campaign's own code; E20 asks the question it leaves open.

## E19 — mesh transfer of the 5x5x3 witnesses, campaign route (10 October 2026)
* **Inputs.** 7x7x5 event list regenerated on C: from the hash-pinned Olympics inputs (`aln_events.py --mesh 7 7 5
  --sigma 0.1 --nsig 4 --T 300`, `ALN_DATASET` = C: copy, same SHA-256); 481,899 events. Witness tensors: the
  campaign's own rebuilt tensors `results/m553/witness/fc3_<name>.npy` (E12, E15, E17, E18).
* **Route.** `p3_events.py --mesh 7 7 5 --fc3 name=path` (phono3py 4.5.0 Rust amplitudes) -> `mesh_transfer_check.py`
  (fc3class.FC3Class, orbit-averaged rates, Cholesky). No code from the review is used.
* **Reference.** kappa 243.027 / 272.334 W/mK, tau 52.609 / 50.623 ps — identical to E1 (7x7x5).
* **Result** (`results/m775x/mesh_transfer.json`; ratios to the 7x7x5 reference; lifetimes degeneracy-averaged):

| witness (5x5x3 tau ratio) | kappa_xx | kappa_zz | tau_x | tau_z | lifetimes |
|---|---|---|---|---|---|
| cal1 max x (2.345) | 0.916 | 0.932 | **0.830** | 0.918 | -31 % ... +44 % |
| cal1 max z (2.151) | 1.019 | 0.965 | 0.898 | **1.081** | -30 % ... +54 % |
| cal1 min x (0.698) | 0.929 | 0.912 | 0.689 | 0.700 | -17 % ... +115 % |
| cal1 min z (0.684) | 0.931 | 0.838 | 0.719 | 0.539 | -66 % ... +261 % |
| cal2 max x (3.721) | 0.853 | 0.827 | **0.705** | 0.643 | -47 % ... +101 % |
| cal2 max z (5.953) | 0.994 | 0.703 | 0.786 | **0.527** | -37 % ... +338 % |
| orb 0.1 max x (1.287; E15 zero start, not the E16 best 1.308) | 1.018 | 0.998 | 1.087 | 1.138 | -29 % ... +21 % |
| orb 0.1 max z (1.370) | 1.023 | 1.031 | 1.037 | **1.347** | -23 % ... +22 % |

  Agreement with the review's independent operator (r8) and with phono3py's own LBTE at 7x7x5 to three to four
  digits. **Conclusion (corrected 11 October, see E24):** the calibrated maxima of E17–E18 are not properties of AlN fc3; they are produced by tuning
  amplitudes at the 75-point mesh's triplets. On 7x7x5 every witness violates kappa (by 1.8–30 %) and the lifetimes
  (by up to +338 %); the calibrated maxima reverse sign for x.
* **Linear model at 7x7x5.** `linearize_parts.py` (23,313 rows, exact; Phi0 = G c0 to 1.6e-10). Projected to the
  scaled coordinates of each class (G M^-1/2; formula checked on 5x5x3 against the stored files to 6.7e-16). At the
  cal1 witnesses the linear 7x7x5 model reproduces the phono3py-route ratios above to <= 5e-5
  (`results/joint_modelcheck.py`).

## E20 — joint-mesh optimisation (data on 5x5x3 AND 7x7x5)  [joint_mesh.py]
* **Question.** Does an fc3 exist that keeps the 300 K data on two samplings of the operator and still moves tau?
  One unknown u (scaled coordinates of the class) drives the exact linearised models of both meshes; all mode-orbit
  lifetimes (+-10 %) and kappa_xx, kappa_zz (+-1 %) are imposed on BOTH meshes; objective log tau on 7x7x5.
  Same augmented-Lagrangian scheme, margins and feasibility test as E15–E18; maxiter 300, 15 outer iterations,
  start u = 0. Classes: calibrated k = 1 with the 266 ALAMODE-cutoff orbits at the 0.5 % floor ("cal1nz", the
  review's construction rebuilt on C:: cond(M) 5.4e8, isometry 2.5e-10) and uniform per-cluster eps = 0.1.
* **Witnesses** (best feasible iterate; verified by phono3py amplitudes on both meshes, orbit-averaged convention,
  `results/joint/verify_m553.json`, `verify_m775.json`: kappa within 0.98–1.0099, degeneracy-averaged lifetimes
  within 9.98 %, tau equal to the optimiser's to <= 1e-5):

| run | tau ratio 7x7x5 (x, z) | tau ratio 5x5x3 (x, z) | Frobenius size | largest orbit change |
|---|---|---|---|---|
| cal1nz max x | **1.799**, 1.138 | 1.056, 1.047 | 1.74 % | 13.0 x orbit norm (bound active) |
| cal1nz max z | 0.975, **3.053** | 1.023, 1.070 | 1.75 % | 7.1 x |
| orb 0.1 max x | **1.326**, 1.089 | 1.057, 1.055 | 4.1 % | 0.0999 |
| orb 0.1 max z | 1.013, **1.690** | 1.009, 1.184 | 3.3 % | 0.0999 |

  Every inner L-BFGS-B solve stopped at maxiter (as in E15–E18, review section 3.1): inner bounds, not optima.
* **Reading.** Imposing the data on 5x5x3 does not prevent large tau changes on 7x7x5; the tau change concentrates
  on the objective mesh (5x5x3 tau moves by <= 18 %). Anatomy as in E18 (calibrated: on-site and < 2 A 0.5 %,
  2–3.2 A 7 %, 3.2–4 A 23 %, 4–5 A 37 %, 5–6 A 49–50 %).

## E21 — validation on a third mesh, 9x9x6 (not used in any optimisation)  [aln_events.py, p3_events.py, mesh_transfer_check.py]
* Event set 9x9x6, sigma 0.1 THz, 4 sigma, 300 K (`results/m996/`). Reference: kappa 265.76 / 256.90 W/mK, tau
  61.40 / 59.15 ps.
* **The reference memory time is not mesh-converged** at sigma = 0.1 THz: tau_x 34.25 / 52.61 / 61.40 ps and tau_z
  27.45 / 50.62 / 59.15 ps on 5x5x3 / 7x7x5 / 9x9x6 (kappa_xx 214.3 / 243.0 / 265.8, kappa_zz 202.0 / 272.3 / 256.9).
  The discretisation change of tau (+79 % / +116 %) exceeds every robust fc3 effect found below.
* Results (`results/joint/verify_m996.json`, `results/m996/transfer_controls.json`; ratios to the 9x9x6 reference):

| tensor | kappa_xx | kappa_zz | tau_x | tau_z | lifetimes (deg.-avg.) |
|---|---|---|---|---|---|
| ALAMODE (second workflow) | 0.9945 | 0.9937 | 0.978 | 0.932 | -28 % ... +19 % |
| E18 cal1 max x (5x5x3 only; 2.345 there) | 0.938 | 0.949 | 0.985 | 0.972 | -31 % ... +29 % |
| E20 cal1nz max x (1.799 on 7x7x5) | 0.952 | 0.965 | **1.183** | 0.958 | -26 % ... +56 % |
| E20 cal1nz max z (3.053 on 7x7x5) | 1.014 | 0.918 | 1.122 | **0.973** | -35 % ... +28 % |
| E20 orb 0.1 max x (1.326 on 7x7x5) | 0.994 | 0.997 | **1.182** | 1.074 | -18 % ... +21 % |
| E20 orb 0.1 max z (1.690 on 7x7x5) | 0.995 | 0.967 | 1.026 | **1.144** | -26 % ... +13 % |

* **Conclusion.** On a mesh that took no part in the optimisation every witness loses most or all of its effect
  (3.05 -> 0.97, 1.80 -> 1.18, 1.69 -> 1.14, 1.33 -> 1.18) and violates the data again (kappa by up to 8 %,
  lifetimes by up to +56 %). Two-mesh constraints do not produce sampling-robust witnesses: the optimiser exploits
  whichever triplets the objective mesh samples. No fc3 in either class has been found that changes tau robustly
  across samplings by more than about 20 % — and those ~18 % candidates already break the data on 9x9x6.
* Note: the +-10 % lifetime tolerance is narrower than the inter-workflow spread at every mesh (ALAMODE: 18 %,
  27 %, 28 % on 5x5x3, 7x7x5, 9x9x6), so the data class is tighter than first-principles precision.
* Missing: the 9x9x6 amplitudes of the E15 orb 0.1 max z witness (the computation was stopped for lack of memory and
  not repeated).

## E22 — continuation runs of the review on 5x5x3 (recorded from review/runs, D:)
Per-cluster eps = 0.1 (start: eps 0.2 witnesses halved; maxiter 600–1000): max tau_x 47.87 ps (1.398), max tau_z
39.40 ps (1.436), min tau_x 28.38 ps (0.829), min tau_z 21.97 ps (0.801) — the E15/E16 range [0.84, 1.31] x
[0.81, 1.37] is superseded by at least [0.83, 1.40] x [0.80, 1.44]. cal1nz on 5x5x3: max 1.858 (x) / 1.820 (z),
min 0.743 (x) / 0.728 (z). All 5x5x3-only statements; E19–E21 show they are sampling-specific.

## E23 — infrared structure of the memory moment (test T1 of the peer review; RTA level)  [ir_rta_moment.py, ir_cutoff.py, ir_branch.py, ir_branch_share.py]
* **Data.** phono3py 4.5.0 RTA per-mode output of the transport-baseline campaign (D:\...\20261002-aln-transport-
  baseline\runs; same hash-pinned inputs; 300 K): tetrahedron method on 9x9x5, 15x15x9, 31x31x17 (the benchmark's
  phono3py grid); Gaussian sigma = 0.05 / 0.1 / 0.2 THz on 15x15x9 and 19x19x11.
* **Quantity.** With C diagonal, x = b/r, so tau_mem(RTA) = N/K = sum kappa_mode tau / sum kappa_mode with
  tau = 1/(4 pi gamma): the kappa-weighted mean lifetime. (The full event operator differs; on 7x7x5 and 9x9x6 the
  peer review finds x_mu / (b_mu / r_mu) = 0.96–1.21 for the lowest modes, so the infrared part is RTA-like.)
* **Totals** (`results/ir/ir_rta.json`):

| mesh, integration | f_min (THz) | kappa_xx / kappa_zz | tau_mem,x | tau_mem,z |
|---|---|---|---|---|
| 9x9x5, tetrahedron | 2.30 | 231.2 / 212.4 | 60.2 ps | 52.2 ps |
| 15x15x9, tetrahedron | 1.30 | 244.2 / 227.0 | 69.7 ps | 64.2 ps |
| 31x31x17, tetrahedron | 0.69 | 270.5 / 247.1 | **135.4 ps** | **99.8 ps** |
| 15x15x9, sigma 0.05 / 0.1 | 1.30 | 256.4 / 254.1 (xx) | 85.9 / 82.5 ps | 84.2 / 79.2 ps |
| 19x19x11, sigma 0.05 / 0.1 / 0.2 | 1.07 | 269.3 / 266.8 / 260.7 (xx) | 124.8 / 119.9 / 106.8 ps | 82.5 / 80.4 / 76.9 ps |

  kappa changes by 11 % from 15x15x9 to 31x31x17; tau_mem doubles.
* **Cutoff test** (`ir_cutoff.json`): the part of N from modes above f_c, divided by the full K:
  f_c = 3 THz: x 31.10 / 31.62 / 31.70 ps, z 27.01 / 26.85 / 27.80 ps (15x15x9 tet / 19x19x11 sigma 0.05 / 31x31x17
  tet); f_c = 4 THz: 17.5 / 19.0 / 17.3 and 15.3 / 16.7 / 15.7 ps. **The mesh dependence of tau_mem lies entirely in
  acoustic modes below about 3 THz**, which on 31x31x17 supply 77 % (x) and 72 % (z) of tau_mem but 24 % / 23 % of K.
* **Concentration** (`ir_branch_share.json`; 31x31x17): below 2 THz the LA branch supplies 35.1 ps of tau_mem,x
  (absent on 15x15x9: no LA mode below 2 THz there), the two TA branches 14.6 and 18.5 ps. One LA star alone
  (q = (1/31, 0, 0), 6 points, f = 1.29 THz, tau = 1.06 ns) contributes 28.7 ps (21 %) of tau_mem,x. Low-q LA
  linewidths are irregular: gamma = 7.5e-5 THz at (1/31,0,0), 1.7e-4 at (0,0,1/17), 1.2e-4 at (1/31,1/31,0)
  (f = 2.22 THz), 7.9e-4 at (0,0,2/17) (f = 2.57 THz).
* **Low-frequency exponents** of gamma ~ f^n (31x31x17, `ir_branch_m313117.json`): TA branches n = 1.1–1.4 below
  2 THz (7–10 points), 1.2–2.1 over 1.5–4 THz; LA n = 2.3 (1.5–3 THz, 7 points) and 2.9 (2–4 THz). A power law
  gamma ~ f^n with constant group velocity and classical heat capacity gives N ~ integral f^(2-2n) df: divergent at
  f -> 0 iff n >= 3/2, while K converges for n < 3. On these data the TA branches lie below the threshold and the LA
  branch above it, but the windows are pre-asymptotic and few points carry them.
* **Status.** NUMERICALLY DEMONSTRATED (RTA, phono3py, up to 31x31x17): tau_mem is not converged on the benchmark's production
  grid 31x31x17 (kappa changes by 11 % from 15x15x9, tau_mem by 94 %), its mesh dependence is infrared (< 3 THz) and is carried by a
  handful of acoustic q-stars nearest Gamma, while the part above 3 THz is converged to about 2 %.
  UNRESOLVED: whether the bulk limit exists. PLAUSIBLE (review F1; LA exponent >= 2 in the accessible window): it does
  not, i.e. tau_mem of the bulk three-phonon BTE is infrared-divergent and must be regularised (film thickness,
  frequency cutoff, or a finite observation frequency).
* **Bearing on E15–E22.** The optimisation objective on every mesh is dominated by that mesh's few lowest acoustic
  modes (review M2a), which differ from mesh to mesh. This explains why witnesses do not transfer and makes the
  transfer test a weak discriminator between "artefact" and "physical ambiguity".

## E24 — numerical audit and peer review: corrections and new data (11 October 2026)
Sources: `review/numerical-audit.md` (independent operator, eigendecomposition with the null pair removed; scripts
`review/na_*.py`, outputs `review/na_out/`) and `review/peer-review.md` (`review/peer-review-checks/`). Every number of
E19–E21 is reproduced to <= 1e-13 by an independent operator.
* **Linear algebra is not the issue** (audit 2.1–2.3): cond(H) 148 / 210 / 351 (5x5x3 / 7x7x5 / 9x9x6), clean null
  vector, tau invariant to the regulariser to 1e-13 at gamma = tr C / n, 4-sigma cut converged to <= 0.06 %, events
  with tiny rates irrelevant (total |elasticity| <= 3e-6). On even c-meshes b . e = 4.7e-9 |b_z| (band-sticking plane
  q_z = 1/2); harmless at the campaign's gamma, should be projected out.
* **Mesh sequence of the event operator** (audit 2.6): tau_x / tau_z = 34.25 / 27.45 (5x5x3), 45.06 / 44.40 (6x6x4),
  52.61 / 50.62 (7x7x5), 61.40 / 59.15 (9x9x6), **73.57 / 70.05 ps (12x12x8)**; observed order of the family
  6x6x4 / 9x9x6 / 12x12x8 p = -0.07 ... -0.17 — no convergent regime; tau grows about like N^0.15. Fewer than 0.5 % of
  the modes (the Gamma-adjacent acoustic modes) carry 25–65 % of N. Within-mesh exponent of the lowest acoustic
  rates: 2.01 / 1.95 / 1.61 (7x7x5 / 9x9x6 / 12x12x8).
* **Smearing at fixed mesh** (audit 2.4, exact re-weighting to sigma' <= 0.1 THz): tau_x at sigma' = 0.05 / 0.10 THz:
  40.4 / 34.3 (5x5x3), 84.5 / 52.6 (7x7x5, -38 %), 64.7 / 61.4 (9x9x6), 74.8 / 73.6 (12x12x8). 5x5x3 and 7x7x5 are
  under-resolved in the smearing; on 7x7x5 witness ratios move with sigma' (cal1nz J max x 1.80 -> 2.17 at 0.08 THz).
* **Convention** (audit 2.5): on even c-meshes the orbit-average convention differs from the gauge-invariant block
  average by 1–5 % for z quantities (9x9x6: tau_z 59.15 vs 59.78 ps, kappa_zz 256.9 vs 251.8); witness ratios change
  by <= 1.8e-3.
* **Missing 9x9x6 transfers, computed by the audit** (audit 3.4; campaign route agrees to 1e-13):

| tensor (5x5x3; 7x7x5 ratio) | tau_x | tau_z | kappa_xx | kappa_zz | lifetimes |
|---|---|---|---|---|---|
| E15 orb 0.1 max z (1.370; 1.347) | 1.025 | **1.518** | 1.025 | 1.046 | -35 % ... +20 % |
| review c1_z_cont02 (1.436; 1.315) | 1.044 | **1.599** | 1.038 | 1.058 | -37 % ... +21 % |
| E15 orb 0.1 max x (1.287; 1.087) | **1.157** | 0.980 | 1.018 | 0.968 | -31 % ... +22 % |
| review c1_x_cont02 (1.387; 1.048) | **1.118** | 1.058 | 1.012 | 0.977 | -39 % ... +22 % |

  **The c-axis enhancement of the per-cluster eps = 0.1 class persists across all three meshes** (1.37 -> 1.35 ->
  1.52; 1.44 -> 1.32 -> 1.60): these tensors slow the long-wavelength c-axis TA modes (lambda_min(H) 0.63–0.66 of the
  reference on 9x9x6), which every mesh samples. What fails on the finer meshes is the data (kappa_zz +5–6 %,
  lifetimes to -37 %). The calibrated-class and E20 witnesses instead engineer one slow eigenpair among the
  objective mesh's Gamma-adjacent modes (one pair carries 58–88 % of N there; no such concentration on 9x9x6).
* **Two levers** (peer review M4, NUMERICAL OBSERVATION): on the objective mesh the gain is mostly off-diagonal
  (in-scattering "drag" on the lowest modes, log-share 0.37–0.64 of 0.43–0.71); the part that transfers is diagonal
  (lower low-frequency acoustic rates; E15 orb 0.1 max z on 7x7x5: 0.226 diagonal of 0.265).
* **Optimiser** (audit 5): the E20 best iterates are near first-order KKT points (residual 0.7–2.0 % of the
  objective gradient; 113–174 active orbit bounds); calibration on 5x5x3 maps a 4 % residual to a 7–9 % shortfall.
  Inner bounds a few per cent below the local maxima of their basins; single start; global optimality unknown.
* **Agreement bounds corrected:** linear model vs phono3py route <= 6e-5 (E19) and <= 1.8e-4 (E20), entirely the
  block- versus orbit-average convention (the amplitudes agree to 1.5e-11); E20 "tau equal to the optimiser's to
  <= 1e-5" holds for x (<= 2e-7), not for z (<= 3.3e-4).
* **Rigorous floor per mesh** (peer review): (1 - dK) tau_CS allows ratios down to 0.40 / 0.39 (7x7x5) and
  0.36 / 0.31 (9x9x6); no joint minimisation was run.
