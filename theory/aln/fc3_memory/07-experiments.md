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

## E25 — long-wavelength acoustic linewidths at fixed q, converged in the integration mesh  [ir_linewidths.py, run_e25.sh]
* phono3py 4.5.0 RTA, 300 K, linear tetrahedron, NAC as in the transport baseline; the same physical q evaluated on
  integration meshes 48x48x24, 96x96x48, 144x144x72 (`results/ir_lw/lw_m*.json`). Smoke test: the 31x31x17 values of
  the transport baseline are reproduced exactly (basal LA at 1.29 THz: gamma 7.54e-5 THz).
* **Convergence at fixed q** (48 / 96 / 144): basal LA at 0.834 THz 4.246 / 4.087 / 4.067e-5; basal LA at 1.664 THz
  1.262 / 1.250 / 1.248e-4; c-axis TA at 0.490 THz 4.804 / 4.841 / 4.838e-5; c-axis LA at 0.913 THz 6.21 / 6.58 /
  6.90e-5 (not converged); lowest basal TA at 0.454 THz 6.99 / 8.02 / 8.77e-5 (not converged; q = 1, 2, 3 grid steps),
  at 0.911 THz 1.601 / 1.543 / 1.533e-4 (converged at 4–6 steps). **Linewidths of modes within a few grid steps of
  Gamma are integration-limited** (the lowest basal TA branch is worst); the lowest stars of any uniform mesh — which
  carry tau_mem (E23, E24) — have linewidth errors of tens of per cent, in either direction (1/144 step: 1.62e-4
  where the converged trend gives below 1e-4).
* **Fixed sum rule:** the Olympics fc3 and fc2 satisfy the acoustic sum rule to 7e-14 / 4e-14 (max |fc3| 68.3);
  symmetrising them in phono3py changes no linewidth (all four digits). The low-q anomaly is not a sum-rule effect.
* **Low-frequency exponents from converged or converging values** (gamma ~ f^n): basal LA n = 1.2 (0.56–0.83 THz),
  1.6 (0.83–1.66), 1.5 (1.66–2.49); basal TA (upper) n = 1.26 (0.47–0.94); basal TA (lower) n <= 0.8 (0.45–0.91,
  using the 3-step value at 0.45 THz, which is still rising); c-axis TA n = 0.8 (0.16–0.49), 1.4 (0.49–0.98), 1.3
  (0.98–1.46); c-axis LA n about 2 (0.61–0.91, lower point at 2 steps only) — undetermined.
* **Bearing on F1.** Below about 1 THz every branch with a converged value has n = 0.8–1.3 < 3/2, for which the RTA
  memory moment converges (integrand f^(2-2n)); the n of 1.5–2.9 seen in 1.5–4 THz windows (E23) is pre-asymptotic.
  The c-axis LA branch is the exception still open. Status: infrared divergence of tau_mem NOT SUPPORTED by the
  converged linewidths; a finite bulk limit PLAUSIBLE but reached only on meshes far finer than 31x31x17, because the
  0.3–3 THz modes (n near 1) carry contributions decaying only like f_min^(3-2n) ~ f_min^0.4–1.

## E26 — boundary-regularised memory time (film regulator, RTA)  [ir_film.py]
* RTA per-mode data of E23; Matthiessen 1/tau_L = 1/tau_3ph + |v_z|/L (diffuse grey boundary for a c-oriented film of
  thickness L; HEURISTIC regulator, not a film solution of the BTE) or |v|/L. tau_mem(L) = sum c v^2 tau_L^2 / sum c
  v^2 tau_L (`results/ir/ir_film.json`).
* **Cross-plane (z), |v_z|/L**, meshes 9x9x5 / 15x15x9 / 19x19x11 (sigma 0.05) / 31x31x17: L = 30 nm 5.09 / 5.19 /
  5.23 / 5.27 ps; 100 nm 10.32 / 10.45 / 10.59 / 10.68; 300 nm 17.9 / 18.5 / 19.1 / 19.3; 1 um 29.5 / 32.4 / 34.7 /
  35.6; 3 um 40.2 / 46.8 / 53.0 / 56.0; bulk 52.2 / 64.2 / 82.5 / 99.8 ps.
* The film-regularised cross-plane memory time is mesh-converged to about 1 % (30–100 nm), 1–2 % (300 nm) and
  2.5 % (1 um) between the two finest meshes; it is not converged beyond a few micrometres. In-plane (x) with the
  |v_z|/L regulator stays infrared-controlled (in-plane modes do not see the boundaries).
* Status: NUMERICALLY DEMONSTRATED for this RTA + Matthiessen model: for films up to about 1 um the memory time is a
  well-defined, mesh-converged first-principles quantity (5.3 ps at 30 nm to 35.6 ps at 1 um, cross-plane, 300 K).

## E28 — third workflow: ShengBTE / thirdorder.py fc3 (Quantum ESPRESSO, PBEsol)  [shengbte_compare.py]
* Inputs downloaded from the benchmark repository at commit 0640f077 (`D:\ResearchLab\datasets\phonon-olympics\
  0640f077...\AlN_shengbte\`, SHA256SUMS; no licence in the source repository, not redistributed): FORCE_CONSTANTS_3RD
  (1,444 triplet blocks, eV/A^3), CONTROL (structure), QE inputs. **The QE inputs use PBEsol PAW pseudopotentials**
  (Al.pbesol-n-kjpaw_psl.1.0.0, N.pbesol-n-kjpaw_psl.1.0.0), confirming for this workflow the common functional
  stated by the peer review (M6). Lattice a = 3.116 A, c/a = 1.6007 (phono3py: 3.113 A, 1.6004); opposite polarity
  (axial Al-N bond along +c).
* **Mapping** (position-based, independent of the ALAMODE code): mirror z -> -z required — the unmirrored map fails
  the basis match (error 0.237 in fractional units), the mirrored map matches the basis to 1.4e-3 and every cluster
  atom to 6.7e-5. Folded tensor: permutation-symmetric exactly, sum rule 7e-7 relative, symmetric-space residual
  7e-12 (no ASR) / 2e-7 (with ASR). Largest pair distance 3.65 A (thirdorder cutoff): no constants beyond.
* **Differences by shell** (relative norm of the difference, `results_c/m553/shengbte_mirror/shengbte_compare.json`):

| shell | share of norm^2 | ALAMODE vs phono3py | ShengBTE vs phono3py | ShengBTE vs ALAMODE |
|---|---|---|---|---|
| on-site | 30.7 % | 0.34 % | **1.30 %** | 1.34 % |
| < 2 A | 68.4 % | 0.36 % | **1.61 %** | 1.49 % |
| 2–3.2 A | 0.68 % | 7.6 % | 12.4 % | 11.6 % |
| 3.2–4 A | 0.064 % | 27 % | 33 % | 27 % |
| 4–6 A | 0.12 % | 42–80 % | 100 % (zero) | 88–98 % |

  Frobenius distance to phono3py 4.1 % (ALAMODE 2.4 %). The dominant constants of the third workflow differ by about
  1.5 %, four times the phono3py–ALAMODE difference on which the calibrated class (E17–E18) was built; the
  single-pair calibration understated realistic precision on the dominant orbits.
* **Transport** (fc2 phono3py, fc3 replaced; `transfer_workflows.json` in m553/sbt_pp, m775x, m996):

| fc3 | mesh | kappa_xx | kappa_zz | tau_x | tau_z | lifetimes (deg.-avg.) |
|---|---|---|---|---|---|---|
| ShengBTE | 5x5x3 | 0.965 | 0.991 | 0.968 | 1.090 | -18 % ... +74 % |
| ShengBTE | 7x7x5 | 0.974 | 0.970 | 1.009 | 1.007 | -23 % ... +33 % |
| ShengBTE | 9x9x6 | 0.977 | 0.981 | 1.024 | 1.027 | -22 % ... +34 % |
| ALAMODE | 5x5x3 / 7x7x5 / 9x9x6 | 1.001 / 0.997 / 0.995 | 1.033 / 1.014 / 0.994 | 0.993 / 0.977 / 0.978 | 1.017 / 1.006 / 0.932 | up to -28 % / +19 % |

* **Reading** (NUMERICAL OBSERVATION, three workflows, one functional): independent first-principles fc3 tensors move
  tau_mem by -7 % ... +9 % on every mesh, while single-mode lifetimes move by up to 74 % and kappa by up to 3.5 %. The
  realistic inter-workflow spread of tau_mem (at fixed mesh) is an order of magnitude below the worst-case ranges of
  the per-orbit classes. Second functional: not yet available (no DFT code installed; planned via Quantum ESPRESSO).

## E25 status amended — model level of the infrared exponents (novelty protocol, `review/novelty-protocol.md` §4)
* The exponents of E25 are Fermi-golden-rule three-phonon linewidths (thermal phonons infinitely long-lived). For
  omega tau_th << 1 the damping of long-wavelength sound enters the Akhiezer regime (attenuation ∝ omega^2; Woodruff &
  Ehrenreich 1961, Simons 1964), and mode-coupling theory gives sound damping ∝ q^2 and a t^(-3/2) heat-current
  tail in 3D (Lepri, Livi & Politi, Phys. Rep. 377 (2003) 1, Eqs. 134, 141; Ernst, Hauge & van Leeuwen 1971). Both
  give n = 2 >= 3/2 asymptotically. The criterion "finite iff n < 3/2" is a one-line corollary of Vermeersch et al.,
  PRB 91, 085202 (2015), Eqs. 19–20, and is presented as such.
* **Amended status.** For the three-phonon golden-rule operator: finite bulk limit UNRESOLVED (converged exponents
  0.8–1.3 below 1 THz are consistent with a finite limit; c-axis LA undetermined). For the physical crystal: bulk
  heat-flux memory time expected infinite (LITERATURE SUPPORTED, n = 2 in the Akhiezer / hydrodynamic regime; the
  crossover near 1/(2 pi tau_th), OUR ESTIMATE GHz to tens of GHz at 300 K, lies far below every mesh used here and
  is not computed). The E25 line "divergence NOT SUPPORTED" is withdrawn; S7 of the report ("UNRESOLVED; divergence
  PLAUSIBLE") stands once the model level is stated.
* Consequence: identifiability must be posed for a regularised, mesh-converged quantity — the film-regularised
  cross-plane memory time tau_mem(L) of E26 (L <= 1 um removes the Akhiezer modes, whose mean free paths exceed L) or
  a frequency-cutoff memory time — not for the bulk tau_mem of the discretised operators (E10–E21).

## E30 — film-regularised memory time of the conserving event operator  [film_event.py, queue/film_e30.sh]
* C_L = C(g) + diag(|v_z| / (4 pi L)) (operator units THz; physical boundary rate |v_z|/L, Matthiessen as E26), now
  on the full conserving operator instead of RTA; C_L is positive definite, no energy-mode regulariser. L = inf
  reproduces E1 (5x5x3: 34.253 / 27.446 ps). Rates orbit-tied (FC3Class); `symmetry.orbits(allow_missing=True)` skips
  the 88 unexported selection-rule images of the 12x12x8 list (0 on the other meshes). Results `results_c/film/`.
* **Reference, cross-plane tau_mem,z(L) [ps]** (5x5x3 / 7x7x5 / 9x9x6): L = 30 nm 5.30 / 5.52 / 5.34; 100 nm
  11.08 / 11.50 / 10.89; 300 nm 17.7 / 19.8 / 18.8; 1 um 23.2 / 31.4 / 31.4; bulk 27.4 / 50.6 / 59.2. The full-operator
  values at 30–100 nm agree with the mesh-converged RTA values of E26 (5.27, 10.68 ps on 31x31x17) to 1–4 %; the
  spread over these coarse meshes is 4 % (30 nm, 100 nm), 11 % (300 nm), against a factor 2.2 for the bulk.
* **Inter-workflow spread** (ALAMODE, ShengBTE fc3; 5x5x3 / 7x7x5): tau_z(30 nm) -2.1 % ... +0.7 %, tau_z(100 nm)
  -1.4 % ... +1.0 %.
* **Bulk witnesses evaluated in the film** (ratios to the reference on the same mesh; witnesses optimised for the
  BULK memory time, E15–E20 and review continuations):

| witness | mesh | bulk tau_z | tau_z(1 um) | tau_z(300 nm) | tau_z(100 nm) | tau_z(30 nm) |
|---|---|---|---|---|---|---|
| orb 0.1 max z | 5x5x3 / 7x7x5 / 9x9x6 | 1.370 / 1.347 / 1.518 | 1.248 / 1.111 / 1.103 | 1.124 / 1.036 / 1.035 | 1.027 / 1.006 / 1.007 | 0.987 / 0.997 / 0.998 |
| c1_z_cont02 | 5x5x3 / 9x9x6 | 1.436 / 1.599 | 1.280 / 1.127 | 1.129 / 1.044 | 1.018 / 1.010 | 0.980 / 0.999 |
| orb 0.2 max z | 5x5x3 | 1.655 | 1.426 | 1.199 | 1.029 | 0.970 |
| cal1 max z | 5x5x3 / 7x7x5 | 2.151 / 1.081 | 1.711 / 1.054 | 1.299 / 1.018 | 1.016 / 0.995 | 0.935 / 0.989 |
| cal2 max z | 5x5x3 / 7x7x5 | 5.953 / 0.527 | 2.540 / 0.673 | 1.257 / 0.782 | 0.871 / 0.866 | 0.863 / 0.931 |
| cal1 min z | 5x5x3 / 7x7x5 | 0.684 / 0.539 | 0.737 / 0.712 | 0.814 / 0.838 | 0.918 / 0.924 | 1.005 / 0.971 |
| cal1nz joint max z | 5x5x3 / 7x7x5 | 1.070 / 3.053 | 1.042 / 1.153 | 1.014 / 0.894 | 0.991 / 0.896 | 0.986 / 0.945 |

  In-plane: the x witnesses keep most of their effect on tau_x(100 nm) (e.g. orb 0.1 max x 1.334 at 5x5x3, 1.218 at
  9x9x6), as expected — the |v_z|/L regulator does not remove in-plane long-wavelength modes.
* **Reading** (NUMERICAL OBSERVATION): the bulk-tau freedom of E10–E21 lives in long cross-plane mean free paths
  that a film of L <= 100 nm removes; the cross-plane film memory time responds to it by at most -14 % ... +3 %
  (30–100 nm) where the bulk value moves by factors up to 6. This is a property of these witnesses only: the
  attainable range of tau_z(L) itself requires optimising tau_z(L) directly (E31).

## E31 — attainable range of the film memory time (direct optimisation)  [robust_mesh.py --film]
* robust_mesh.py extended with `--film L`: objective tau_z(L) of C_L, constraints the bulk static data (all lifetimes
  within 10 %, kappa_xx, kappa_zz within 1 %), per-orbit class bounds. Gradient: dK_L/dg, dN_L/dg take the bulk form
  with C_L (D independent of fc3); finite-difference check at a random point 3.4e-10 (x) / 1.5e-11 (z).
* **5x5x3, per-cluster class eps = 0.1, L = 100 nm** (single start from the reference, 15 outer iterations, both
  feasible, violation <= 9e-8; `results_c/film_opt/m553_orb0.1_z_L100_{max,min}.json`): tau_z(100 nm) ratio
  **[0.909, 1.098]** (bulk tau_z at these optima 0.987 / 1.121; kappa_xx, kappa_zz at -0.98 % / +0.98 %, lifetimes at
  the 10 % bound, every active orbit at its bound). The same class gives bulk tau_z in [0.80, 1.37] (E15/E16, E22).
  Inner bounds (feasible iterates; global optimality unknown). Other L, classes and two-mesh robust runs:
  `queue/film_opt.sh` (in progress).

## E29 — persistent c-axis witnesses on 12x12x8 (fourth mesh, not used in any optimisation)  [mesh_transfer_check.py]
* Event list and amplitudes of the numerical audit (`review/na_scratch/events_m12128_s0.1.npz`, 10,733,499 events,
  13,821 modes; phono3py 4.5.0 amplitudes `review/na_scratch/m12128w/pp_r01_*.npz`), evaluated with the campaign's
  own operator (FC3Class, orbit-averaged; the 88 unexported selection-rule images skipped). Reference tau 73.57 /
  70.05 ps, kappa 280.44 / 267.35 (event operator), identical to the audit (S8). `results_c/m12128/mesh_transfer.json`.

| witness | tau_z ratio 5x5x3 / 7x7x5 / 9x9x6 / **12x12x8** | 12x12x8 tau_x | kappa_xx | kappa_zz | lifetimes (deg.-avg.) |
|---|---|---|---|---|---|
| orb 0.1 max z (E15) | 1.370 / 1.347 / 1.518 / **1.186** | 1.037 | **+3.1 %** | +1.4 % | **-25 % ... +30 %** |
| c1_z_cont02 (review continuation) | 1.436 / 1.32 / 1.599 / **1.205** | 1.054 | **+4.3 %** | +2.4 % | **-25 % ... +35 %** |

* Both persistent c-axis witnesses keep a +19–20 % bulk tau_z on the fourth mesh, but violate the static data there
  (kappa by 3–4 times the 1 % tolerance, lifetimes by 2.5–3.5 times the 10 % tolerance). The statement of the report
  "every witness violates the data on a mesh not used in its construction" now holds on 7x7x5, 9x9x6 and 12x12x8.

### E30 addendum — 12x12x8 (12 October 2026)
* `results_c/film/film_m12128.json`; the 88 skipped symmetry images carry rate <= 1e-18 of the largest (true
  selection-rule zeros), so skipping them is exact to that level.
* **Reference tau_mem,z(L) [ps], 5x5x3 / 7x7x5 / 9x9x6 / 12x12x8:** 30 nm 5.30 / 5.52 / 5.34 / **5.49**; 100 nm
  11.08 / 11.50 / 10.89 / **11.23**; 300 nm 17.7 / 19.8 / 18.8 / 19.7; 1 um 23.2 / 31.4 / 31.4 / 33.9; bulk 27.4 / 50.6 /
  59.2 / **70.1**. Over four meshes the 30 nm and 100 nm values stay within +-2 % and +-3 % of their mean; the bulk
  value grows by a factor 2.6. (NUMERICAL OBSERVATION; Gaussian sigma = 0.1 THz, orbit tying; the smearing and
  convention checks of E24 have not yet been repeated for the film quantity.)
* Persistent c-axis witnesses on 12x12x8: orb 0.1 max z tau_z ratio 0.997 (30 nm), 1.005 (100 nm), 1.025 (300 nm),
  1.057 (1 um), 1.186 (bulk); c1_z_cont02 0.998 / 1.009 / 1.035 / 1.077 / 1.205.

### E31 update — 5x5x3, per-cluster class eps = 0.1 (feasible optima)
* L = 30 nm: max **1.060** (bulk tau_z at that optimum 1.008); L = 300 nm: max **1.164** (bulk 1.256). Minima at 30 and
  300 nm and the 1 um runs were interrupted by the drive disconnection of 11 October 22:43 and restarted
  12 October 21:10 (`queue/film_opt.sh`, resumable).

### E31 update — 5x5x3 film ranges (12 October, 23:15; feasible optima, violation <= 1e-7 unless stated)
| class | L | min tau_z(L) (bulk tau_z there) | max tau_z(L) (bulk tau_z there) | bulk range of the class |
|---|---|---|---|---|
| per-cluster eps 0.1 | 30 nm | 0.937 (1.064) | 1.060 (1.008) | [0.80, 1.37] (E15/E16, E22) |
| per-cluster eps 0.1 | 100 nm | 0.909 (0.987) | 1.098 (1.121) | |
| per-cluster eps 0.1 | 300 nm | 0.882 (0.842) | 1.164 (1.256) | |
| per-cluster eps 0.1 | 1 um | 0.844 (0.814) | 1.261 (1.361) | |
| calibrated k = 1 (cal1) | 100 nm | **0.802** (1.310) | **1.247** (1.457) | [0.684, 2.151] (E18) |
| three-workflow (cal3w) | 100 nm | (running; 0.811 infeasible) | (running; 1.264 at violation 2e-4) | not run for bulk |

* **Reading.** The admissible range of tau_z(L) contracts relative to the bulk in every class, but its width depends on
  the precision class: about +-6 % (30 nm) to +-10 % (100 nm) for a uniform 10 % per-cluster bound, about -20 % ...
  +25 % at 100 nm when each orbit may move by the observed phono3py–ALAMODE difference. The film minimum of the
  calibrated class has a bulk tau_z of 1.31: film and bulk extremes are reached by different tensors.
* Caveat inherited from E19–E21: the calibrated-class extremes of E17–E18 were 5x5x3-specific. These film extremes
  are single-mesh too; the two-mesh runs (`m553m775_*`) and the 9x9x6 validation decide whether they are properties of
  AlN fc3 or of the 75-point mesh. First two-mesh iterate (orb 0.1, 100 nm, infeasible 1e-3): 1.068 / 1.044.
* 9x9x6 linearisation complete (`results_c/m996/merge.done`; float32 projections for orb, cal1nz, cal3w, cal3wnz;
  Phi0 = G c0 to 3.0e-6 relative, float32 storage).

## E32 — mesh transfer of the film-optimised tensors (E31, built on 5x5x3)  [witness_verify.py build, p3_events.py, film_event.py, mesh_transfer_check.py]
* Tensors rebuilt from the optimiser coordinates (`results_c/m553/witness/fc3_<class>_fL<L>_<sense>.npy`), amplitudes
  by phono3py on each mesh (`results_c/<mesh>/filmw_pp`), film operator as E30. **Check on 5x5x3:** the phono3py route
  reproduces the optimiser's film ratios (1.098 / 0.909 / 1.247 / 0.802 at 100 nm).
* **7x7x5** (`results_c/film/filmw_m775.json`, `filmw_m775_transfer.json`; ratios to the 7x7x5 reference):

| tensor (5x5x3 film ratio) | tau_z(L) on 7x7x5 | kappa_xx | kappa_zz | lifetimes |
|---|---|---|---|---|
| orb 0.1, 30 nm max (1.060) | 1.019 | 0.976 | 0.970 | -15 % ... +21 % |
| orb 0.1, 30 nm min (0.937) | 0.984 | 1.020 | 1.002 | -15 % ... +9 % |
| orb 0.1, 100 nm max (1.098) | 1.015 | 1.011 | 0.964 | -20 % ... +21 % |
| orb 0.1, 100 nm min (0.909) | 0.974 | 0.996 | 0.992 | -16 % ... +21 % |
| orb 0.1, 300 nm max (1.164) | 1.018 | 1.024 | 0.991 | -23 % ... +26 % |
| orb 0.1, 300 nm min (0.882) | 0.939 | 0.979 | 0.958 | -24 % ... +54 % |
| orb 0.1, 1 um max (1.261) | 1.096 | 1.025 | 1.021 | -22 % ... +28 % |
| orb 0.1, 1 um min (0.844) | 0.850 | 0.973 | 0.941 | -27 % ... +82 % |
| cal1, 100 nm max (1.247) | 0.992 | 0.978 | 0.920 | -34 % ... +61 % |
| cal1, 100 nm min (0.802) | 0.884 | 0.925 | 0.790 | -61 % ... +739 % |

* **Reading.** As for the bulk witnesses (E19), every single-mesh film extremum violates the static data on the next
  mesh (kappa by 2–21 %, lifetimes by 9–739 % beyond tolerance) and keeps only a small part of its effect: at 30–100 nm
  the transferred tensors move tau_z(L) by -3 % ... +2 % (per-cluster class) and -12 % ... -1 % (calibrated). The
  single-mesh film ranges of E31 are therefore upper estimates of the freedom; whether any of it survives is decided by
  the multi-mesh robust runs (E31 two-mesh, three-mesh to follow). 9x9x6 evaluation in progress.

### E31 update — remaining 5x5x3 runs and first two-mesh result (13 October, 01:30)
* Single-mesh 5x5x3, 100 nm, feasible (violation <= 1e-7): three-workflow cal3w **0.802 / 1.269** (bulk tau_z at the
  optimum 1.248 / 1.524); cal3wnz 0.830 / 1.209 (1.173 / 1.364); per-cluster eps = 0.2: 0.870 / 1.154 (1.025 / 1.237).
  Lane A (`m553_*`) complete: 16 runs, all feasible.
* **Two-mesh** orb 0.1, 100 nm, max (5x5x3 + 7x7x5 data and objective; `m553m775_orb0.1_z_L100_max`): first feasible
  iterate at outer 4 (violation 8.7e-5 against tolerance 1e-4): film ratio **1.071 (5x5x3) / 1.045 (7x7x5)**; bulk
  tau_z 1.063 / 0.935. Still iterating. Requiring the data on a second mesh halves the single-mesh maximum (1.098 to
  1.045 on the finer mesh) but does not remove it.

### E32 addendum — 9x9x6 (13 October)  [results_c/film/filmw_m996.json, filmw_m996_transfer.json]
No missing symmetry images on 9x9x6. Ratios to the 9x9x6 reference (tau_z(L) at the optimised L; data columns bulk):

| tensor | tau_z(30 nm) | tau_z(100 nm) | tau_z(L_opt) | kappa_xx | kappa_zz | lifetimes |
|---|---|---|---|---|---|---|
| orb 0.1, 30 nm max | 1.008 | 0.994 | 1.008 | 0.965 | 0.933 | -18 % ... +32 % |
| orb 0.1, 30 nm min | 0.991 | 0.996 | 0.991 | 1.021 | 1.046 | -14 % ... +17 % |
| orb 0.1, 100 nm max | 1.009 | 1.005 | 1.005 | 1.003 | 0.937 | -19 % ... +34 % |
| orb 0.1, 100 nm min | 0.995 | 0.987 | 0.987 | 0.995 | 1.016 | -15 % ... +20 % |
| orb 0.1, 300 nm max | 1.007 | 1.008 | 1.005 | 1.020 | 0.974 | -21 % ... +19 % |
| orb 0.1, 300 nm min | 1.000 | 0.986 | 0.955 | 0.974 | 0.987 | -16 % ... +29 % |
| orb 0.1, 1 um max | 1.001 | 1.010 | 1.083 | 1.028 | 1.029 | -27 % ... +18 % |
| orb 0.1, 1 um min | 1.005 | 0.991 | 0.907 | 0.970 | 0.988 | -16 % ... +32 % |
| cal1, 100 nm max | 1.008 | 0.991 | 0.991 | 0.961 | 0.885 | -34 % ... +80 % |
| cal1, 100 nm min | 0.952 | 0.894 | 0.894 | 0.906 | 0.783 | -16 % ... +317 % |

* **Reading.** Same pattern as on 7x7x5: every tensor violates the 9x9x6 data, and at 30–100 nm the per-cluster
  tensors keep -1.3 % ... +1 % of their effect. The one sizable residue (cal1 min, -11 % at 100 nm) comes with
  kappa_zz -22 % and lifetimes up to +317 %. It is not admissible and gives no evidence of freedom. The longer films
  (300 nm, 1 um) keep more (-9 % ... +8 %), as the bulk does.

## E33 — three-mesh robust film runs  [queue/film_opt3.sh, robust_mesh.py --film --meshes m553,m775,m996]
Launched 13 October, 01:30, one run at a time (9x9x6 model about 10 GB): orb 0.1 at 100 and 30 nm (max/min), then
cal1nz, cal3wnz and cal3w at 100 nm, then orb 0.1 at 300 nm. Outputs `results_c/film_opt/m3_*`. 12x12x8 validation of
the optima to follow.

## E34 — temperature dependence of kappa and tau_mem for the bulk witnesses on 9x9x6  [review/r10_temperature.py; results_c/m996/temperature_witnesses.json]
REF (9x9x6, bulk-regularised operator): kappa_xx = 3109 / 546 / 265.8 / 133.3 / 78.4 W/mK and tau_x = 1599 / 157.5 /
61.4 / 25.6 / 13.8 ps at 100 / 200 / 300 / 500 / 800 K. Witness ratios kappa_xx/kappa_zz and tau_x/tau_z:

| witness | 100 K | 300 K | 800 K |
|---|---|---|---|
| orb joint, x (orb_J_x) | k 1.013/1.005, tau 1.063/0.998 | k 0.994/0.997, tau 1.182/1.074 | k 0.988/0.991, tau 1.211/1.064 |
| cal1nz joint, x | k 0.920/0.940, tau 1.045/0.941 | k 0.952/0.965, tau 1.183/0.958 | k 0.957/0.972, tau 1.246/0.963 |
| orb E15, z | k 0.991/1.013, tau 0.931/1.095 | k 1.025/1.046, tau 1.025/1.518 | k 1.031/1.034, tau 1.062/1.481 |
| ALAMODE | k 0.974/0.981, tau 0.943/0.928 | k 0.995/0.994, tau 0.978/0.932 | k 1.001/0.999, tau 0.995/0.941 |
| ShengBTE | k 0.991/1.021, tau 1.023/1.024 | k 0.977/0.981, tau 1.024/1.027 | k 0.975/0.976, tau 1.024/1.014 |

* **Reading (NUMERICAL OBSERVATION).** The witness that fits the 9x9x6 data at 300 K (orb_J_x: kappa within 0.6 %)
  stays within 1.3 % of kappa from 100 to 800 K and keeps tau_x +6 % ... +21 %. A kappa(T) curve at the 1 % level would
  constrain it only marginally (1.3 % at 100 K). Temperature data therefore does not visibly restore identifiability of
  the bulk tau_mem. A decisive test needs an optimisation with kappa(T) constraints (not yet implemented).
* Cost (13 October, 03:50): one outer iteration takes about 1 h for the two-mesh runs and more than 2.5 h for the
  three-mesh runs (9x9x6 model 15 GB resident). The old stopping rule (violation < 1e-7) is not met by the two-mesh
  runs, which therefore go the full 15 outers. `robust_mesh.py` (backup `robust_mesh.py.bak_prestop`) now also stops
  after outer >= 3 when two consecutive iterates are feasible for the nominal tolerances and the robust log ratio moved
  by less than 1e-4 (`--rtol`). It takes an optional warm start (`--start file:U.npy`, or `PREFIX_start.npy` if present),
  used to seed three-mesh runs from the corresponding two-mesh optima. Applies to runs started after this change; the
  run in progress (`m3_orb0.1_z_L100_max`, zero start) keeps the old rule.

## E27 — RTA on finer meshes (48x48x26 complete; 56x56x30 running)  [rta_chunks.py; results_c/rta/m484826/assembled.json]
phono3py 4.5.0, 300 K, linear tetrahedron, NAC, same inputs as E23. 48x48x26 computed in 48 chunks of irreducible
points; the chunked assembly was checked on 9x9x5 against a single run (kappa identical to 4 digits).

| quantity | 31x31x17 (E23/E26) | 48x48x26 | change |
|---|---|---|---|
| f_min (THz) | 0.69 | 0.452 | |
| kappa_xx / kappa_zz (W/mK) | 270.5 / 247.1 | 272.8 / 250.6 | +0.8 % / +1.4 % |
| tau_mem,x / tau_mem,z, bulk (ps) | 135.4 / 99.8 | 152.8 / 124.1 | +13 % / +24 % |
| part above 3 THz, x / z (ps) | 31.70 / 27.80 | 31.87 / 28.37 | +0.5 % / +2.0 % |
| tau_z(L), 30 / 100 / 300 nm (ps) | 5.27 / 10.68 / 19.3 | 5.295 / 10.728 / 19.45 | +0.5 / +0.4 / +0.8 % |
| tau_z(L), 1 / 3 / 10 um (ps) | 35.6 / 56.0 / — | 36.0 / 57.6 / 83.8 | +1.1 / +2.9 % |

* **Reading (NUMERICALLY DEMONSTRATED, RTA + Matthiessen).** With 2.9 times more q-points than 31x31x17, kappa moves by
  about 1 %, the bulk memory time by a further 13–24 %. It keeps growing with the infrared cutoff, consistent with E23
  and with the expected Akhiezer divergence. The film-regularised cross-plane memory time is converged to 0.5 % at
  30–100 nm and 1 % at 1 um. 56x56x30 (4,640 irreducible points) started 13 October, 05:10.

### E31 update — two-mesh orb 0.1, 100 nm, max: final (13 October, 06:50)
Best feasible iterate at outer 14 (violation 2.6e-7): film tau_z(100 nm) **1.0712 (5x5x3) / 1.0455 (7x7x5)**; robust
ratio stable to 1e-6 over the last three outers. Bulk tau_z at the optimum 1.063 / 0.935 (it lowers the bulk value on
7x7x5 while raising the film value). kappa_xx, kappa_zz at the tolerance edge on 7x7x5 (+0.98 %, +0.98 %). Optimum
`m553m775_orb0.1_z_L100_max_u.npy`; 9x9x6 and 12x12x8 validation via the phono3py route (E32) is pending until memory
allows (the three-mesh run occupies 15 GB). The two-mesh minimum then started, under the new stopping rule.

### E31 update — two-mesh orb 0.1, 100 nm, min: final (13 October, 08:36)
Stopped by the new rule at outer 6 (7 outers, about 2 h). Best iterate feasible for the nominal tolerances at outer 5
(inner-margin violation 2e-4; kappa_zz on 7x7x5 0.99000, at the edge): film tau_z(100 nm) **0.9276 (5x5x3) / 0.9496
(7x7x5)**; bulk tau_z 1.010 / 1.170. The outer-6 iterate gives 0.9496 at violation 2.8e-5.
* **Two-mesh admissible range of tau_z(100 nm), per-cluster class eps = 0.1: [0.950, 1.046]** (on the finer mesh),
  against [0.909, 1.098] single-mesh (E31) and [0.974, 1.015] for the transferred single-mesh extremes (E32). The
  same tensors move the bulk tau_z by -6.5 % and +17 %: the film quantity is several times better constrained.

## E35 — smearing dependence of the film memory time (T2 for the film quantity)  [film_smearing.py; results_c/film_smear/]
Exact re-weighting of the stored sigma = 0.1 THz events to sigma' = 0.05 / 0.06 / 0.08 / 0.10 THz (as audit 2.4), film
operator of E30; ratios to REF at the same sigma'. sigma' = 0.1 reproduces E30 and the E32 transfer ratios exactly.

| | sigma' (THz) | 0.05 | 0.06 | 0.08 | 0.10 |
|---|---|---|---|---|---|
| 5x5x3 REF | tau_z 30 / 100 nm (ps) | 5.563 / 12.232 | 5.444 / 11.779 | 5.338 / 11.317 | 5.298 / 11.078 |
| 7x7x5 REF | tau_z 30 / 100 nm (ps) | 5.595 / 11.767 | 5.576 / 11.697 | 5.544 / 11.586 | 5.520 / 11.502 |
| 7x7x5 | ALAMODE / ShengBTE, 100 nm | 1.0043 / 0.9876 | 1.0040 / 0.9868 | 1.0039 / 0.9861 | 1.0042 / 0.9862 |
| 7x7x5 | orb 0.1 fL100 max / min (E32) | 1.0167 / 0.9731 | 1.0156 / 0.9731 | 1.0149 / 0.9736 | 1.0147 / 0.9743 |
| 7x7x5 | cal1 fL100 max / min (E32) | 0.9946 / 0.8948 | 0.9923 / 0.8911 | 0.9911 / 0.8870 | 0.9918 / 0.8844 |
| 5x5x3 | orb 0.1 fL100 max / min | 1.0955 / 0.9137 | 1.1039 / 0.9064 | 1.1059 / 0.9025 | 1.0982 / 0.9087 |
| 5x5x3 | cal1 fL100 max / min | 1.2634 / 0.8843 | 1.2792 / 0.8526 | 1.2730 / 0.8114 | 1.2468 / 0.8020 |

* **Reading (NUMERICAL OBSERVATION).** On 7x7x5 the reference film memory time changes by 1.4 % (30 nm) and 2.3 %
  (100 nm) over sigma' = 0.05–0.1 THz, against 5 % / 10 % on 5x5x3 (the coarse mesh is under-resolved in the
  smearing, as for the bulk in audit 2.4). Film ratios of the tensors and of the workflows move by <= 1 % on 7x7x5.
  The film-quantity conclusions are therefore not smearing artefacts on 7x7x5. On 5x5x3 the extremes move by up to
  8 % (cal1 min 0.80 to 0.88), another reason not to read the single-mesh 5x5x3 ranges as properties of AlN.
* The three workflows agree on tau_z(100 nm) to within -1.4 % ... +0.4 % on 7x7x5 at every sigma'.
* The 9x9x6 run was stopped by memory pressure (the three-mesh optimisation holds 15 GB) and is pending.

### E33 update — three-mesh orb 0.1, 100 nm, max: first feasible iterate (13 October, 12:30)
Outer 3 (about 12 h; violation 7.7e-5, feasible for the nominal tolerances on all three meshes): film tau_z(100 nm)
**1.058 (5x5x3) / 1.041 (7x7x5) / 1.037 (9x9x6)**; bulk tau_z 1.042 / 0.937 / 0.919; kappa within +1 %. Run continues
(old stopping rule). Sequence of the admissible maximum of tau_z(100 nm) on the finest mesh used: 1.098 (one mesh),
1.046 (two), 1.037 (three, provisional). Two-mesh 30 nm max: 1.0362 at outer 7, feasible.
* Two-mesh orb 0.1, 30 nm, max: final at outer 7 (new stopping rule; violation 1.3e-5): film tau_z(30 nm)
  **1.049 (5x5x3) / 1.036 (7x7x5)**; bulk tau_z 0.992 / 0.895 (single mesh: 1.060). 30 nm min running.
* Two-mesh orb 0.1, 30 nm, min: final at outer 6 (violation 2.5e-5): film tau_z(30 nm) **0.950 (5x5x3) / 0.965
  (7x7x5)**; bulk tau_z 1.064 / 1.180. **Two-mesh admissible range of tau_z(30 nm): [0.965, 1.036]** (single mesh
  [0.937, 1.060]); of tau_z(100 nm): [0.950, 1.046]. Three-mesh 100 nm max at outer 4: 1.0378 (feasible).
* Two-mesh orb 0.1, 300 nm, max: final at outer 6 (violation 2.9e-5): film tau_z(300 nm) **1.105 (5x5x3) / 1.058
  (7x7x5)**; bulk tau_z 1.151 / 1.028 (single mesh 1.164). 300 nm min running.
* Two-mesh orb 0.1, 300 nm, min: final (best feasible at outer 6, violation 1.2e-4 inner margin): film tau_z(300 nm)
  **0.910 (5x5x3) / 0.931 (7x7x5)**; bulk tau_z 0.902 / 1.057.

**Two-mesh admissible ranges of tau_z(L) on 7x7x5, per-cluster class eps = 0.1 (single-mesh 5x5x3 in brackets):**

| L | min | max |
|---|---|---|
| 30 nm | 0.965 (0.937) | 1.036 (1.060) |
| 100 nm | 0.950 (0.909) | 1.046 (1.098) |
| 300 nm | 0.931 (0.882) | 1.058 (1.164) |

The half-width grows with L (3.5 %, 4.8 %, 6.4 %) and is roughly half the single-mesh value at every L. The
calibrated two-mesh runs (cal1nz, 100 nm) follow.
* Two-mesh calibrated cal1nz, 100 nm, max: final at outer 7 (violation 1.3e-5): film tau_z(100 nm) **1.124 (5x5x3) /
  1.079 (7x7x5)**; bulk tau_z 1.154 / 0.914; kappa at the +1 % edge on both meshes. (Single-mesh comparison: cal1, not
  cal1nz, 1.247 on 5x5x3; cal3wnz 1.209.) cal1nz min running (0.90 at outer 1, infeasible).
* Two-mesh calibrated cal1nz, 100 nm, min: final at outer 7 (violation 5.1e-6): film tau_z(100 nm) **0.882 (5x5x3) /
  0.897 (7x7x5)**; bulk tau_z 1.131 / **1.730**. **Two-mesh calibrated range of tau_z(100 nm): [0.897, 1.079]** on
  7x7x5. The same data leave the bulk tau_z between 0.914 and 1.730 at these two optima alone. Lane B and lane A
  complete (`film_opt.done`).

### E27 update — RTA 56x56x30 (14 October; results_c/rta/m565630/assembled.json)

| quantity | 31x31x17 | 48x48x26 | 56x56x30 |
|---|---|---|---|
| f_min (THz) | 0.69 | 0.452 | 0.389 |
| kappa_xx / kappa_zz (W/mK) | 270.5 / 247.1 | 272.8 / 250.6 | 273.1 / 251.8 |
| tau_mem,x / tau_mem,z, bulk (ps) | 135.4 / 99.8 | 152.8 / 124.1 | 158.0 / 132.3 |
| part above 3 THz, x / z (ps) | 31.70 / 27.80 | 31.87 / 28.37 | 31.67 / 27.99 |
| tau_z(L), 30 / 100 / 300 nm (ps) | 5.27 / 10.68 / 19.3 | 5.295 / 10.728 / 19.45 | 5.299 / 10.737 / 19.48 |
| tau_z(L), 1 / 3 / 10 um (ps) | 35.6 / 56.0 / — | 36.0 / 57.6 / 83.8 | 36.1 / 58.1 / 85.2 |

* From 48x48x26 to 56x56x30 kappa moves by 0.1 / 0.5 %, the bulk memory time by +3.4 % / +6.6 % (still growing), the
  film value by 0.08 % (30–100 nm), 0.15 % (300 nm), 0.4 % (1 um), 1.7 % (10 um). The film-regularised cross-plane
  memory time is converged to better than 0.1 % at 30–100 nm (RTA + Matthiessen). 64x64x35 started 14 October, 02:58.
* 9x9x6 smearing check (E35) relaunched 03:45 after memory became available.

### E35 addendum — 9x9x6 (results_c/film_smear/m996.json)
| sigma' (THz) | 0.05 | 0.06 | 0.08 | 0.10 |
|---|---|---|---|---|
| REF tau_z 30 / 100 nm (ps) | 5.364 / 10.957 | 5.351 / 10.928 | 5.340 / 10.902 | 5.339 / 10.894 |
| ALAMODE / ShengBTE, 100 nm | 1.0021 / 0.9935 | 1.0021 / 0.9922 | 1.0022 / 0.9902 | 1.0023 / 0.9889 |
| orb 0.1 fL100 max / min | 1.0040 / 0.9866 | 1.0041 / 0.9863 | 1.0045 / 0.9863 | 1.0050 / 0.9867 |
| cal1 fL100 max / min | 0.9934 / 0.9015 | 0.9927 / 0.8976 | 0.9917 / 0.8940 | 0.9910 / 0.8939 |

* On 9x9x6 the reference film memory time changes by 0.5 % (30 nm) and 0.6 % (100 nm) over sigma' = 0.05–0.1 THz
  (7x7x5: 1.4 / 2.3 %; 5x5x3: 5 / 10 %); tensor and workflow ratios move by <= 0.8 %. **The film quantity is
  converged in the smearing on 9x9x6 at the 1 % level** (NUMERICAL OBSERVATION). The three workflows agree to
  -1.1 % ... +0.2 % at 100 nm.

## E36 — transfer of the two-mesh film optima (E31, 5x5x3 + 7x7x5) to 9x9x6 and 12x12x8  [queue/film_transfer2.sh]
Tensors `results_c/m553/witness/fc3_<class>_m2L<L><sense>.npy` rebuilt from the optimiser coordinates; amplitudes by
phono3py on each mesh (`<mesh>/film2_pp`); film operator as E30; static data by mesh_transfer_check.py.
* **7x7x5 (reproduction check, a data mesh):** film ratios equal the optimiser's (100 nm: 1.0456 / 0.9495 orb,
  1.0790 / 0.8965 cal1nz); all lifetimes within +-10 %, kappa within +-1 %. Pipeline confirmed.
* **9x9x6 (not used by the optimiser):**

| tensor | tau_z 30 / 100 / 300 nm on 9x9x6 | on 7x7x5 at L_opt | bulk tau_z | kappa_xx / kappa_zz | lifetimes |
|---|---|---|---|---|---|
| orb 30 nm max | **1.022** / 1.012 / 0.974 | 1.036 | 0.762 | 0.992 / 0.957 | -15 % ... +30 % |
| orb 30 nm min | **0.977** / 0.980 / 1.005 | 0.965 | 1.288 | 0.996 / 1.028 | -20 % ... +14 % |
| orb 100 nm max | 1.018 / **1.020** / 1.001 | 1.046 | 0.774 | 1.000 / 0.954 | -14 % ... +34 % |
| orb 100 nm min | 0.982 / **0.974** / 0.979 | 0.950 | 1.154 | 0.994 / 1.011 | -14 % ... +15 % |
| orb 300 nm max | 1.005 / 1.016 / **1.028** | 1.058 | 0.933 | 1.007 / 0.977 | -16 % ... +13 % |
| orb 300 nm min | 1.000 / 0.988 / **0.970** | 0.931 | 1.002 | 0.995 / 1.009 | -11 % ... +17 % |
| cal1nz 100 nm max | 1.029 / **1.033** / 1.005 | 1.079 | 0.696 | 1.009 / 0.926 | -18 % ... +89 % |
| cal1nz 100 nm min | 0.967 / **0.947** / 0.938 | 0.897 | 0.895 | 0.990 / 0.943 | -31 % ... +21 % |

* **Reading.** None of the two-mesh optima satisfies the 9x9x6 data: lifetimes leave the +-10 % band in every case
  (by 3–79 points) and kappa_zz by up to 7 %. The data on two meshes therefore still do not imply the data on a
  third. On 9x9x6 the tensors keep about half of their film effect (orb 100 nm: +2.0 % / -2.6 %; cal1nz: +3.3 % /
  -5.3 %), against one quarter or less for the single-mesh optima (E32). The bulk tau_z of the same tensors moves by
  -30 % ... +29 %. Each added data mesh roughly halves the film freedom (one mesh about +-10 %, two meshes +-5 %,
  three meshes: +3.8 % maximum at 100 nm, E33). Whether this sequence tends to zero or to a finite width is the open
  question, to be addressed by the 12x12x8 evaluation of the three-mesh optima.
* 12x12x8 evaluation of the two-mesh optima in progress.

### E33 update — 100 nm max converged in value; warm starts for the remaining runs (14 October, 06:40)

- m3_orb0.1_z_L100_max, outer 10 (t = 102 688 s): robust ratio 1.037929, violation 8.4e-7, feasible. Outers 7 to 10 change the
  ratio by less than 2e-5 in total. Per-mesh film ratios at 100 nm (tau_x, tau_z): 5x5x3 1.027 / 1.0583,
  7x7x5 0.967 / 1.0411, 9x9x6 0.956 / 1.0379, so the binding mesh for the tau_z objective is 9x9x6. Bulk
  tau_z ratios at this iterate: 1.042 / 0.936 / 0.917. The run was started before the --rtol rule and continues to outer 15
  (about 2.8 h per outer). Its tensor is written only at the end, so it is left to finish.
- The remaining runs of the queue start from staged initial points (PREFIX_start.npy, recorded as start_used):
  the two-mesh optima (E31) for orb 0.1 at 30, 100 (min), 300 nm and for cal1nz at 100 nm; the 5x5x3 optima for cal3w and
  cal3wnz. Together with --rtol 1e-4 this cuts each run from 15 outers to an expected 4 to 6.
  Caveat: the problem is nonconvex, and a warm start can select a different local optimum than a zero start. The max run
  above started from zero and reached 1.0379. Its min partner starts from the two-mesh min (bulk-feasible on 5x5x3 and
  7x7x5). If time allows, one pair will be repeated from zero to check this.

### E36 note — 12x12x8 film metrics rerun with a low-memory film_event.py (14 October, 07:50)

The first 12x12x8 film_event pass failed at the REF tensor (MemoryError while assembling C; the three-mesh optimisation holds
14.6 GB of 27.9 GB). film_event.film_metrics now assembles C on the (n+1)^2 bincount buffer with the dummy slot decoupled
(row and column 0, diagonal 1), rebuilds it for each L and factorises it in place; the bulk rank-one term is added by row
blocks. Peak memory drops from about four n^2 arrays to one. Check: on 7x7x5 (REF, orb0.1_m2L100max, cal1nz_m2L100min, all L,
kappa, tau, tauRTA) the new code reproduces results_c/film/film2_m775.json to a maximum relative difference of 4.4e-16.
Previous version kept as scripts_c/film_event.py.bak_lean. The rerun (queue/film_e36_m12128.sh) starts after the 12x12x8
mesh_transfer_check of film_transfer2.sh finishes.

Literature status of the film result: review/novelty-film.md (14 October). C1 (identifiability map, narrowing with the number
of data meshes), C2 (multi-mesh fc3 optimisation) and C3 (bulk non-convergence against converged film value) are each
PARTIALLY KNOWN, with the AlN-specific results POTENTIALLY NOVEL WITHIN THE SEARCH SCOPE. Two caveats from that report go
into the analysis: (i) adding data meshes can only shrink the admissible range for a fixed objective, so the informative
content is the rate and the limit of the narrowing, not the narrowing itself; (ii) the data are generated by the same
operator that is inverted (inverse-crime setting), so the ranges are inner bounds of the class at fixed discretisation.

### E36 addendum — 12x12x8 (14 October, 09:55)  [results_c/film/film2_m12128.json, film2_m12128_transfer.json]

Same eight tensors, 12x12x8 events (review/na_scratch/events_m12128_s0.1.npz, sigma 0.1 THz), REF amplitudes
m12128w/pp_r01_REF.npz; 88 events lack a symmetry image, with largest rate 1.0e-18 of the maximum, which is negligible.
The film metrics use the low-memory film_event.py (E36 note). Film tau_z ratio to REF at the optimised L, and bulk tau_z
ratio, on the fitted mesh and the two held-out meshes:

| tensor | L | tau_z(L) 7x7x5 | 9x9x6 | 12x12x8 | bulk 7x7x5 | 9x9x6 | 12x12x8 |
|---|---|---|---|---|---|---|---|
| orb 0.1 max | 30 | 1.0362 | 1.0215 | 1.0228 | 0.895 | 0.762 | 0.813 |
| orb 0.1 min | 30 | 0.9648 | 0.9771 | 0.9757 | 1.180 | 1.288 | 1.174 |
| orb 0.1 max | 100 | 1.0456 | 1.0196 | 1.0209 | 0.935 | 0.774 | 0.856 |
| orb 0.1 min | 100 | 0.9495 | 0.9741 | 0.9734 | 1.170 | 1.154 | 1.122 |
| orb 0.1 max | 300 | 1.0580 | 1.0284 | 1.0259 | 1.028 | 0.933 | 0.956 |
| orb 0.1 min | 300 | 0.9315 | 0.9699 | 0.9665 | 1.057 | 1.002 | 1.040 |
| cal1nz max | 100 | 1.0790 | 1.0332 | 1.0321 | 0.914 | 0.696 | 0.765 |
| cal1nz min | 100 | 0.8965 | 0.9469 | 0.9498 | 1.730 | 0.895 | 0.960 |

Static data on 12x12x8 (tolerances: lifetimes ±10 %, kappa ±1 %). All eight tensors violate them. Degeneracy-averaged
lifetime deviations: orb 30 max -13 %/+26 %, 30 min -16 %/+13 %, 100 max -15 %/+23 %, 100 min -18 %/+15 %,
300 max -13 %/+16 %, 300 min -17 %/+13 %, cal1nz max -17 %/+54 %, cal1nz min -20 %/+17 %.
kappa (xx, zz) ratios: orb 30 max 0.985/0.971, 30 min 1.003/1.009, 100 max 0.996/0.975, 100 min 1.000/1.004,
300 max 1.006/0.983, 300 min 0.995/1.010, cal1nz max 1.001/0.952, cal1nz min 1.004/0.957.

Findings (NUMERICAL OBSERVATION):
1. Each tensor's film value is mesh-converged on the held-out meshes. Between 9x9x6 and 12x12x8 the film ratio at the
   optimised L changes by at most 0.34 points, and at 30 and 100 nm by at most 0.14 points for the per-cluster tensors (cal1nz min 0.29). Over the same
   step the bulk ratio of the same tensors moves by up to 11.4 points (orb 30 min 1.288 -> 1.174; orb 100 max 0.774 -> 0.856), and cal1nz min moves from
   1.730 on 7x7x5 to 0.895 and 0.960. This is the E27 and E30 contrast again, now for perturbed tensors and the full
   operator: the film quantity is a property of the tensor, while the bulk quantity is a property of the discretisation.
2. The film excursion on the fitted 7x7x5 mesh overstates the converged excursion of the same tensor by a factor of 1.5
   to 2.5 (orb 100 max: +4.6 % on 7x7x5, +2.1 % converged; cal1nz min: -10.4 % and -5.0 %). Part of what the optimiser
   exploits is discretisation-specific, consistent with E32 and E36 on 9x9x6.
3. No tensor fitted on 5x5x3 + 7x7x5 satisfies the 12x12x8 data, so these converged excursions (about ±2-3 % per-cluster,
   +3.2/-5.0 % calibrated) are not admissible ranges. They are what tensors admissible on two coarse meshes do to the
   converged film quantity. The admissible range on fine meshes remains the open question, to be answered by E33 (three
   meshes, current maximum 1.0379 on 9x9x6) and its 12x12x8 validation.

### E33 update — three-mesh orb 0.1, 100 nm, max: final (14 October, 08:16)  [results_c/film_opt/m3_orb0.1_z_L100_max.json, _u.npy]

Stopped at outer 13 by the original rule (violation 6.6e-9 < 1e-7, ratio change < 1e-6), t = 110 159 s (30.6 h).
Robust ratio (smooth minimum of the film tau_z(100 nm) ratio over the three meshes) **1.03793**, binding mesh 9x9x6.
Per-mesh film ratios at 100 nm (tau_x, tau_z): 5x5x3 1.0266 / 1.0583, 7x7x5 0.9674 / 1.0411, 9x9x6 0.9562 / 1.0379.
Bulk tau ratios (x, z): 1.0305 / 1.0420, 0.9971 / 0.9357, 0.9813 / 0.9169. Data active: kappa_xx and kappa_zz at +0.98 % on
7x7x5 and 9x9x6, kappa_zz +0.98 % on 5x5x3, max lifetime deviation 9.98 % on all three meshes, orbit bound 0.999 eps.
Comparison for the same objective at 100 nm (max): one data mesh 1.098 (5x5x3), two 1.046 (7x7x5 binding), three 1.0379.
The min partner (warm start from the two-mesh min) started at 08:16. Held-out 12x12x8 validation of this tensor:
E37 (queue/m3_validate2.sh).

## E37 — held-out validation of the three-mesh film optima (E33) on 12x12x8  [queue/m3_validate2.sh; results_c/film3/]

Per tensor: witness_verify.py build (m553 basis) -> p3_events.py amplitudes on the 7x7x5, 9x9x6 and 12x12x8 event sets
-> film_event.py (low-memory) -> mesh_transfer_check.py. REF amplitudes as in E36. (A first copy of the queue,
queue/m3_validate.sh, had a positional-argument error in the output names. It was made inert with skip markers before it
wrote any final result. Its stray files results_c/film3/orb0.1_m3L100max_.json and _transfer.json are the 7x7x5
evaluation and are identical to the _m775 files.)

### orb 0.1, 100 nm, max (three-mesh optimum 1.03793)

| mesh | role | film tau_z(100 nm) ratio | bulk tau_z ratio | kappa (xx, zz) ratio | lifetime dev. min/max |
|---|---|---|---|---|---|
| 7x7x5 | fitted | 1.0412 (optimiser 1.0411) | 0.9357 | 1.0098 / 1.0098 | -9.98 % / +9.98 % |
| 9x9x6 | fitted | 1.0372 (optimiser 1.0379) | 0.9212 | 1.0097 / 1.0069 | -10.07 % / +10.15 % |
| 12x12x8 | held out | **1.0327** | 0.9497 | 1.0058 / 1.0108 | -11.6 % / +11.9 % |

On the fitted meshes, the phono3py route reproduces the optimiser's film value to 0.07 points. The 9x9x6 lifetime bound is
exceeded by at most 0.15 points, which is the difference between the optimiser's 9x9x6 model and the phono3py route. On the
held-out 12x12x8 mesh, the tensor keeps 86 % of its film excursion (+3.27 % against +3.79 %). The two-mesh optimum kept 46 %
(+2.09 % against +4.56 %, E36). The data are violated only marginally: lifetimes by at most 1.9 points beyond ±10 %, and
kappa_zz by 0.08 points beyond 1 %. For the two-mesh optima the excess was up to 16 points and 1.5 points (E36).
NUMERICAL OBSERVATION: with three data meshes, the maximiser is nearly admissible on the next finer mesh, and its film
excursion is nearly mesh-stable. This suggests that the admissible film range at 100 nm in this class converges to a finite
width of about +3 %, not to zero. The min partner, the other lengths and the 12x12x8 evaluation of the other classes are
needed before this is stated. mesh_transfer_check.py now also records the fraction of modes outside ±10 % and the 50/90/99 %
quantiles of the lifetime deviation (added after this tensor's check).

### E33 update — three-mesh orb 0.1, 100 nm, min: final (14 October, 18:58)  [results_c/film_opt/m3_orb0.1_z_L100_min.json, _u.npy]

Warm start from the two-mesh min (start_used file:...m3_orb0.1_z_L100_min_start.npy). It stopped at outer 4 by the --rtol rule
(feasible at outers 3 and 4, log-ratio change 1e-5). Violation 2.3e-5, t = 37 903 s (10.5 h; the max run took 30.6 h).
Robust ratio (smooth maximum of the film tau_z(100 nm) ratio over the three meshes) **0.96317**, binding mesh 9x9x6.
Per-mesh film ratios at 100 nm (tau_x, tau_z): 5x5x3 0.9719 / 0.9321, 7x7x5 0.9645 / 0.9483, 9x9x6 0.9661 / 0.9632.
Bulk tau ratios (x, z): 0.9633 / 1.0260, 0.9605 / 1.1924, 0.9781 / 1.1051. Active data: kappa_zz -0.98 % on all three
meshes, kappa_xx -0.98 % on 9x9x6, lifetimes 9.98 % on all three meshes, orbit bound 0.999 eps.

Admissible film tau_z(100 nm) range in the per-cluster eps = 0.1 class, as a function of the number of data meshes
(inner bounds; the binding mesh is the finest data mesh in each case):

| data meshes | min | max | width |
|---|---|---|---|
| 5x5x3 | 0.909 | 1.098 | 0.189 |
| 5x5x3 + 7x7x5 | 0.950 | 1.046 | 0.096 |
| 5x5x3 + 7x7x5 + 9x9x6 | 0.963 | 1.038 | 0.075 |

The width ratios of successive rows are 0.51 and 0.78, so the narrowing slows. Held-out 12x12x8 check of the min tensor:
E37.

### E37 — orb 0.1, 100 nm, min (three-mesh optimum 0.96317)

| mesh | role | film tau_z(100 nm) ratio | bulk tau_z ratio | kappa (xx, zz) ratio | lifetime dev. min/max | modes outside ±10 % | |dev| q50/q90/q99 |
|---|---|---|---|---|---|---|---|
| 7x7x5 | fitted | 0.9483 (optimiser 0.9483) | 1.1924 | 0.9911 / 0.9902 | -9.98 % / +9.98 % | 0 | 2.9 / 6.9 / 10.0 % |
| 9x9x6 | fitted | 0.9644 (optimiser 0.9632) | 1.1056 | 0.9900 / 0.9913 | -10.4 % / +10.0 % | 0.07 % | 2.7 / 7.4 / 10.0 % |
| 12x12x8 | held out | **0.9658** | 1.1019 | 0.9949 / 0.9884 | -15.6 % / +12.3 % | 0.34 % | 2.6 / 6.0 / 9.0 % |

On 12x12x8 the tensor keeps 93 % of its film excursion (-3.42 % against -3.68 %). It violates the data marginally:
kappa_zz by 0.16 points beyond 1 %, and lifetimes for 0.34 % of the modes (99 % of the modes stay within 9.0 %).

**Summary of E33 and E37 at 100 nm (per-cluster eps = 0.1).** Film tau_z(100 nm) ratio of the extremal tensors:

| data meshes | optimum (binding mesh) | held-out 12x12x8 | data on 12x12x8 |
|---|---|---|---|
| 2 (5x5x3, 7x7x5) | 0.950 / 1.046 | 0.973 / 1.021 | violated (lifetimes to -18 % / +23 %, kappa_zz -2.5 %) |
| 3 (+ 9x9x6) | 0.963 / 1.038 | 0.966 / 1.033 | nearly kept (lifetimes outside ±10 % for <= 0.34 % of modes, kappa within 1.2 %) |

NUMERICAL OBSERVATION (one length, one class): with three data meshes the extremal tensors are nearly admissible on the
finer held-out mesh, and their film excursions are stable to 0.3–0.5 points. The admissible film range at 100 nm
appears to converge to roughly [0.96, 1.04]: a finite width of about ±3.5 %, not zero. Statement status: SUPPORTED BY
EVIDENCE for this class at 100 nm. It still needs: the other lengths (30 and 300 nm, queued), the calibrated classes (queued),
an independent multistart (local optima), and the numerical audit.

## E38 — Block C: second functional (PBE) for fc3, launched 15 October  [scripts_blockc/; runs in WSL Ubuntu, /home/koussay/blockc]

Quantum ESPRESSO 6.7 (Ubuntu package, Open MPI, 4 ranks; WSL capped at 4 CPUs and 7 GB so the film queues keep their
memory). pslibrary 1.0.0 PAW pseudopotentials, the same family as the reference ShengBTE workflow (E28): Al/N
`pbe-n-kjpaw_psl.1.0.0` and `pbesol-n-kjpaw_psl.1.0.0` (SHA256 in WSL `pseudo/`). Numerics as the reference: ecutwfc 120 Ry,
ecutrho 600 Ry, conv_thr 1e-12 (primitive) / 1e-11 (supercells).
Design: PBE against PBEsol **in the same 3x3x2 supercell (72 atoms), Gamma point, phono3py displacements (0.03 A) with pair
cutoff 3.66 A**. The reference tensor's largest pair distance is 3.65 A, and a shell lies between 3.65 and 3.70 A, so the counts are
318 / 370 / 466 / 466 / 1254 supercells for cutoffs 3.2 / 3.65 / 3.7 / 4.0 / none. The functional difference is read at fixed
numerics. The PBEsol branch (reference structure) also validates the pipeline against the reference fc3. The 4x4x3 supercell (192
atoms, same count) was rejected on cost (about 15–20 times more per SCF).
Pipeline (resumable, skips finished outputs): PBE vc-relax (12x12x12 k) -> 370 + 1 PBE supercell SCFs -> 370 + 1 PBEsol
supercell SCFs -> PBEsol vc-relax (consistency check of the reference structure). Progress log: WSL `blockc/progress.log`.
Cost per supercell SCF: to be measured from the first runs. The full 742 SCFs are not expected to finish by 17 October.

### E38 correction (15 October, after 21 PBE supercells)
1. **Pair cutoff.** The cutoff of 3.66 A split the functionals. In PBEsol (a = 3.116 A) the third Al-N shell lies at 3.646/3.652 A and
   is kept (466 supercells). In PBE (a = 3.133 A) it lies at 3.666/3.672 A and was dropped (318 supercells). The reference tensor's
   largest pair, 3.65 A, is this shell. Corrected to **3.9 A** (next shell 4.38 / 4.40 A), giving 466 supercells for both
   functionals. The earlier count "370 at 3.65 A" cut through the shell's two sub-distances. phono3py numbering does not depend on
   the cutoff: the 22 finished PBE cells (disp-00000..00021) are identical in the 3.9 A set (checked: 0 position changes).
2. **Cost and cutoff test.** A 72-atom SCF at 120/600 Ry on 4 WSL ranks takes 2160–2570 s while the Windows queues load the CPU,
   so 932 SCFs would take about 24 days. Pipeline v2 runs a test whose rule was fixed before the runs (`ecut_check.py`): it
   recomputes dF = F(disp-00001) - F(perfect) at 80/480 and 60/480 Ry and accepts the lowest cutoff with
   max|dF - dF(120/600)| <= 0.5 % of max|dF(120/600)|. If a lower cutoff is accepted, the 120/600 outputs are kept in
   `pbe_e120/` and both branches run at the chosen cutoff. Otherwise 120/600 continues and the 22 cells are reused.
   Result: `ecut_test/result.json`.

### E33 update: three-mesh orb 0.1, 30 nm, max (final 15 October, 20:45)  [results_c/film_opt/m3_orb0.1_z_L30_max.json, _u.npy]

Restarted from scratch after the reboot of 15 October (no checkpoint). Stopped at outer 5, violation 1.0e-5, t = 64 344 s (17.9 h).
Robust ratio (smooth minimum of the film tau_z(30 nm) ratio over the three meshes) **1.0281**, binding mesh 9x9x6.
Per-mesh film ratios at 30 nm (tau_x, tau_z): 5x5x3 1.0152 / 1.0459, 7x7x5 0.9207 / 1.0359, 9x9x6 0.9028 / 1.0281.
Bulk tau ratios (x, z): 1.0094 / 0.9866, 0.9376 / 0.8956, 0.9137 / 0.8766. Active data: kappa at +0.98 % (5x5x3, 7x7x5),
max lifetime deviation 9.98 % on all three meshes, orbit bound 0.999 eps.

### E37: orb 0.1, 30 nm, max (three-mesh optimum 1.0281)

| mesh | role | film tau_z(30 nm) | film tau_z(100 nm) | bulk tau_z | kappa (xx, zz) | lifetime dev. min/max | modes outside ±10 % | abs. dev. q50/q90/q99 |
|---|---|---|---|---|---|---|---|---|
| 7x7x5 | fitted | 1.0359 (optimiser 1.0359) | 1.0366 | 0.8956 | 1.0098 / 1.0098 | -9.97 % / +9.98 % | 0 | 4.8 / 8.5 / 9.9 % |
| 9x9x6 | fitted | 1.0276 (optimiser 1.0281) | 1.0254 | 0.8776 | 1.0047 / 0.9963 | -10.0 % / +10.0 % | 0.17 % | 4.9 / 8.7 / 10.0 % |
| 12x12x8 | held out | **1.0276** | 1.0251 | 0.8981 | 0.9992 / 1.0014 | -11.4 % / +15.4 % | 0.89 % | 4.6 / 8.2 / 10.0 % |

On 12x12x8 the tensor keeps 98 % of its film excursion at 30 nm (+2.76 % against +2.81 %). Kappa stays within 0.15 %,
and 99 % of the modes stay within the ±10 % lifetime band. The bulk tau_z ratio moves by 2 points between 9x9x6 and 12x12x8,
and the film ratio by less than 0.01 points (the E27/E30 contrast again). NUMERICAL OBSERVATION: at 30 nm the three-mesh maximum
(+2.8 %) is smaller than at 100 nm (+3.8 %; held out +3.3 %). The 30 nm min run started on 16 October, 00:16.

### E38 second correction (16 October): the 120/600 production runs are invalid; production restarted at 60/480 Ry

1. **Failed cutoff test.** The rule of `ecut_check.py` rejected 80/480 and 60/480: the error in dF against 120/600 was 216 % and
   56 %. The 80/480 supercell, however, had an energy 3.48 Ry *below* the 120/600 one. A lower cutoff cannot give a lower
   variational energy, and the Gamma eigenvalues of the *perfect* supercell differed between the three runs by 0.1–1.3 eV, with
   different degeneracy patterns.
2. **Primitive-cell cutoff ladder** (`ecut_ladder/`; PBE, 4 atoms, one Al displaced by 0.0311 A along a, 6x6x4 k, nosym). Force on the
   displaced atom along x, in Ry/bohr: 60/480 -0.0183863, 80/320 -0.0183779, 80/480 -0.0183778, 100/500 -0.0183768. Energies -135.85422,
   -135.85535, -135.85535, -135.85542 Ry, monotone. **Forces are converged to 1e-5 Ry/bohr (5e-4 relative) at 60/480.**
3. **Folding test** (`fold_check/`). The perfect 3x3x2 supercell at Gamma must equal 18 x the primitive cell on a Gamma-centred 3x3x2 k-mesh.

| cutoff (Ry) | 18 x E_prim (Ry) | E_supercell, Gamma (Ry) | difference (Ry) | HOMO prim / supercell (eV) |
|---|---|---|---|---|
| 60/480 | -2445.14298906 | -2445.14298871 | **3.5e-7** | 6.9433 / 6.9433 |
| 80/480 | -2445.16340826 | -2448.99133383 | -3.83 | 6.9427 / 7.2739 |
| 120/600 | -2445.16602366 | -2445.51333121 | -0.35 | 6.9427 / 7.7402 |

   The Gamma-point 72-atom calculation of QE 6.7 (Ubuntu package, 4 MPI ranks) fails at 80/480 and 120/600: it converges to
   energies below the variational minimum. The failure is reproducible and silent: no warning, and SCF accuracy < 1e-11 Ry.
   The cause was not identified. All three cutoffs share the grids of the passing case only in part (dense grid of the
   60/480 run = smooth grid of the 120/600 run). **The 34 PBE supercells computed at 120/600 Ry (disp-00000..00033) are invalid.**
   They are moved to `pbe_e120_invalid/`. The ecut_check.py rule had used them as its reference, so its outcome is void.
4. **Restart at 60/480 Ry** (`ecut.choice.reason`). The two 60/480 supercells of the cutoff test (byte-identical inputs) are reused.
   Before production, a check fixed in advance (`kcheck.py`) compares disp-00001 between the Gamma-trick path and the complex path
   (K_POINTS automatic 1 1 1). The run continues only if |dE| < 1e-6 Ry and max|dF| < 1e-5 Ry/bohr. Cost about 15 min per
   supercell under the present load, so about 5 days for the PBE branch and about 10 days for both.
5. **Still to do.** E28's reference fc3 was produced at 120/600 with another QE version, which is unaffected by this failure. Before
   the functional comparison, the PBEsol branch at 60/480 is validated against the reference fc3 (supercell, k-sampling and
   cutoff differ, so expect a small but non-zero difference). A spot check of a few displaced cells with an independent QE build
   (7.x) is advisable.
   **kcheck result (16 October):** disp-00001 at 60/480, Gamma-trick vs complex path: |dE| = 1.2e-7 Ry, max|dF| = 2.1e-7 Ry/bohr,
   both below the limits fixed in advance (1e-6 Ry, 1e-5 Ry/bohr). PASS. PBE production at 60/480 continues from disp-00002.

### E27 update: RTA 64x64x35 (assembled 16 October, 13:52; results_c/rta/m646435/assembled.json)

| quantity | 48x48x26 | 56x56x30 | 64x64x35 |
|---|---|---|---|
| lowest frequency on the mesh (THz) | 0.452 | 0.389 | 0.336 |
| kappa x / z (W/mK) | 272.8 / 250.6 | 273.1 / 251.8 | 273.6 / 251.8 |
| bulk tau_mem x / z (ps) | 152.8 / 124.1 | 158.0 / 132.3 | 162.4 / 136.4 |
| tau_mem from modes above 3 THz, x / z (ps) | 31.87 / 28.37 | 31.67 / 27.99 | 31.58 / 27.78 |
| film tau_mem z at 30 / 100 / 300 nm (ps) | 5.30 / 10.73 / 19.45 | 5.30 / 10.74 / 19.48 | 5.30 / 10.75 / 19.51 |

The bulk memory time is still growing: +2.8 % (x) and +3.1 % (z) from 56x56x30 to 64x64x35. Kappa changes by at most 0.2 %, the
part above 3 THz by at most 0.8 %, and the film memory time by at most 0.15 % at 30–300 nm. In x the growth per step is linear in
the lowest mesh frequency within the three points (83 ps/THz on both steps). The z increments (130 and 77 ps/THz) are not
consistent with that. HEURISTIC, not a result: a linear extrapolation of x to f_min -> 0 would give about 190 ps. Whether the
bulk limit is finite remains UNRESOLVED (E25: it depends on the infrared exponent of the acoustic linewidths). The
statement of paper 1 is unaffected: only regularised (film or frequency-cutoff) memory times are quoted, and these are converged.

### E33/E37 update (17 October): orb 0.1, 30 nm, min, and cal1nz, 100 nm, max

**orb 0.1, 30 nm, min** (outer 5, 17.5 h, violation 2.7e-5). Robust ratio **0.9727**, binding mesh 9x9x6. Per-mesh film tau_z(30 nm):
5x5x3 0.9499, 7x7x5 0.9636, 9x9x6 0.9727. Held-out check:

| mesh | role | film tau_z(30 nm) | bulk tau_z | kappa (xx, zz) | lifetime dev. min/max | modes outside ±10 % | abs. dev. q50/q90/q99 |
|---|---|---|---|---|---|---|---|
| 7x7x5 | fitted | 0.9636 | 1.2103 | 0.9902 / 0.9902 | -9.9 % / +10.0 % | 0 | 4.1 / 7.9 / 10.0 % |
| 9x9x6 | fitted | 0.9732 (optimiser 0.9727) | 1.1570 | 0.9924 / 0.9997 | -10.3 % / +10.0 % | 0.07 % | 4.4 / 8.1 / 9.9 % |
| 12x12x8 | held out | **0.9730** | 1.1223 | 0.9990 / 0.9924 | -12.6 % / +12.3 % | 0.41 % | 3.9 / 7.3 / 9.5 % |

Admissible film tau_z ratio at 30 nm, three data meshes, per-cluster eps = 0.1: **[0.973, 1.028]**, and the held-out values are
**[0.973, 1.028]** (99–100 % of the excursion kept). At 100 nm the corresponding range is [0.963, 1.038], held out [0.966, 1.033].

**cal1nz, 100 nm, max** (calibrated class, outer 7, 17.3 h). Robust ratio **1.0630**, binding mesh 9x9x6. Per-mesh film tau_z(100 nm):
5x5x3 1.1104, 7x7x5 1.0745, 9x9x6 1.0630. Held-out 12x12x8: film **1.0472** (75 % of the excursion kept), bulk tau_z 0.9481,
kappa 1.0056 / 1.0091, lifetimes -14.5 % / +27.1 %, 1.9 % of the modes outside ±10 %, q99 11.1 %. For the calibrated class the
optimum is less mesh-stable, and its data are violated more on the held-out mesh than for orb 0.1 (two-mesh cal1nz max: +3.2 %
converged, E36). The min partner is running (started 17 October, 07:54).

NUMERICAL OBSERVATION: in the per-cluster class the three-mesh extremal tensors are mesh-stable at both 30 and 100 nm. The
admissible half-width is about 2.8 % at 30 nm and about 3.5 % at 100 nm. Remaining: the cal1nz min, cal3wnz, cal3w and
300 nm runs (queued), an independent multistart, and the numerical audit.
