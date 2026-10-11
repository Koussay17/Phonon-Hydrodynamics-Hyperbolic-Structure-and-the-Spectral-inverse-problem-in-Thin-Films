# Numerical audit — fc3-parametrised memory campaign (20261005-fc3-memory), claims S1, S5, S8 and the optimiser

**Date:** 11 October 2026. **Object:** `13-final-report.md` (version of 11 October) and `07-experiments.md` E1, E14,
E19–E22 (D:), with the scripts and results of the C: mirror (`scripts/`, `results/m553`, `results/m775x`,
`results/m996`, `results/joint`). **Mandate:** decide whether the numbers behind S1, S5, S8 and the optimiser
statements can be trusted (consistency, conditioning, regularisation, discretisation, convention, iteration error).
No campaign file outside `review/` was modified. Scripts: `review/na_*.py`; raw outputs: `review/na_out/`; large
scratch files (two new event lists, five 9x9x6 pp files) in `C:\Users\Koussay\AppData\Local\Temp\na_fc3\`
(regeneration commands in section 8). Environments: research Python (NumPy 2.5.3, SciPy 1.18.1, OpenBLAS, 4 threads)
for all analysis; phono3py 4.5.0 replica environment only for the event export at 6x6x4 / 12x12x8 and for the
9x9x6 amplitudes of four witness tensors (inputs hash-checked by `aln_common.verify_inputs`).

---

## 0. Verdicts

| claim | verdict | decisive numbers |
|---|---|---|
| **S1** Phi0 = G c0 to 1.6e-10 (7x7x5) | **RELIABLE** | all three stored scaled matrices (orb, cal1, cal1nz) reproduce Phi0 with row-wise max 1.6–1.7e-10 (rows of tiny \|Phi0\|), 1.0–1.1e-13 relative to max \|Phi0\|, 99.9 % quantile 7e-12 |
| **S1** linear model vs phono3py route <= 5e-5 | **RELIABLE WITH CAVEATS** | amplitude level: uc\|Phi\|^2 vs Rust pp on 20,141 non-degenerate rows <= 1.5e-11 (pp > 1e-6 max), <= 5.1e-11 (all); tau-ratio differences are up to **5.6e-5** (cal1, abs.; 7.2e-5 rel.) and **1.7e-4** (E20 cal1nz max z), entirely the block- vs orbit-average convention (reproduced to 1e-9 by block-averaging the phono3py pp). E20's "tau equal to the optimiser's to <= 1e-5" holds for x (<= 1.7e-7) but not for z (up to 7.4e-5 at 7x7x5, 3.3e-4 at 5x5x3) |
| **S5** tau ratios of E19 / E21 at 7x7x5 and 9x9x6 | **RELIABLE** | independent operator reproduces all 22 stored ratios to <= 1e-13; cond(C on e^perp) 137–457; null eigenvalue <= 1e-17, next eigenvalue >= 1.9e-4; regularisation-independent at gamma = tr C/n; cut 3.5 vs 4 sigma <= 0.1 %; convention <= 1.8e-3 |
| **S5** wording "every single-mesh witness tested on a finer mesh loses or reverses its memory change" | **NOT RELIABLE** | E19's own table: orb 0.1 max z keeps 1.347 of 1.370 on 7x7x5. New (this audit): on 9x9x6 the same tensor gives **tau_z x 1.518**, the review continuation c1_z_cont02 **1.599** (both routes agree to 1e-13; stable under cut and smearing changes); the data are violated there (kappa_zz +4.6 / +5.8 %, lifetimes down to -35 / -37 %) |
| **S6** "largest tau change of any witness on an unused mesh is 1.18 x" | **NOT RELIABLE** | contradicted on 7x7x5 by E19 (orb 0.1 max z: 1.347) and review 4.2 / 3.3 (1.330, 1.315), and on 9x9x6 by the values above (1.518, 1.599) |
| **S8** reference tau / kappa on 5x5x3 / 7x7x5 / 9x9x6 | **RELIABLE WITH CAVEATS** | reproduced to <= 3e-14; linear algebra converged (cond <= 351, gamma-scan spread <= 1.3e-13 on 5x5x3 / 7x7x5); 4-sigma cut converged (<= 0.06 %); **caveat:** at 9x9x6 the orbit-average convention differs from the gauge-invariant block average by +1.06 % (tau_z) and -1.97 % (kappa_zz) |
| **S8** "tau_mem is not mesh-converged; genuine discretisation effect" | **RELIABLE (understated)** | not a regularisation, conditioning, cut or round-off artefact. Two further meshes: 6x6x4 45.06 / 44.40 ps, 12x12x8 **73.57 / 70.05 ps** (+20 % / +18 % over 9x9x6). The family 6x6x4–9x9x6–12x12x8 is in no convergent power-law regime (observed order -0.07 … -0.17). N is carried by the Gamma-adjacent acoustic modes (<= 0.5 % of the modes carry 25–65 % of N). Classification: **UNDER-RESOLVED** in the mesh; linear algebra **CONVERGED** |
| **Optimiser** (all inner L-BFGS-B solves at maxiter) | **does not change any stated conclusion; RELIABLE WITH CAVEATS** | E20 best iterates are near first-order KKT points: KKT residual 0.7–2.0 % of the objective gradient (154–230 active constraints). Calibration on the 5x5x3 runs: residual 3.5–4.2 % went with a 7–9 % shortfall, 0.6 % with 0.8 %. E20 maxima are inner bounds a few per cent below the local maxima of their basin; global optimality unknown. **TOLERANCE-LIMITED (iteration cap), not unstable** |

---

## 1. Independent implementation used for the checks

`na_common.py` re-implements the operator from the written definitions only (02-assumptions A/B, aln_events.py
docstring), sharing no code with fastops / memgap / aln_geometry / fc3class: own mode indexing, own conserving
stoichiometry, own rates g = conv pp g0w / (4 s_p s_a s_b) (or / 8 s_p s_a^2), own symmetry tying (group average
(1/24) sum_R g(R alpha), which equals the orbit average for a group action; no union-find), own dense assembly, and
**C^+ b by full eigendecomposition with the null pair removed explicitly** (no regularisation). Checks of the inputs:
* rates from the pp files vs the exported event rates: 1.31e-15 (5x5x3), 1.29e-15 (7x7x5), 1.44e-15 (9x9x6);
* g0w is the Gaussian of the detuning to 3e-12 (5x5x3) and 4.4e-6 (7x7x5, 9x9x6; phono3py's own evaluation);
* my orbit partition equals the campaign's union-find partition exactly (7x7x5: 23,674 orbits; `na_linear.py`);
* reference values: tau and kappa agree with the campaign to 2.2e-15 (5x5x3), 2.2e-15 (7x7x5), 3.1e-14 (9x9x6).

---

## 2. S8 — the reference memory time on three (five) meshes

### 2.1 Linear algebra on H = e^perp (`na_spectral.py`; `na_out/spectral_m*.json`)

| mesh | n | lambda_null | lambda_min(H) | lambda_max | cond(H) | eig vs Cholesky (campaign) | gamma-scan spread of tau (gamma/(trC/n) = 1e-8 … 1e4) | rel. residual |
|---|---|---|---|---|---|---|---|---|
| 5x5x3 | 897 | 1.8e-17 | 9.764e-4 | 0.145 | 148 | 2.2e-15 | 1.3e-13 | <= 1.6e-15 (1.4e-12 at gamma = 1e4) |
| 7x7x5 | 2,937 | -2.4e-18 | 4.035e-4 | 0.085 | 210 | 2.2e-15 | 2.8e-14 | <= 1.6e-15 |
| 9x9x6 | 5,829 | -9.9e-18 | 2.877e-4 | 0.101 | 351 | 3.1e-14 | 7e-9 (x), **4.4e-3 (z)**; <= 4e-11 for gamma >= 1e-4 trC/n | <= 1.9e-15 |

* The energy null vector is clean: |C e_hat| / ||C|| <= 1.9e-16; |e_hat . q_null| = 1 - 2e-16; no negative or second
  near-zero eigenvalue (gap lambda_min(H) / |lambda_null| > 1e13). cond(H) grows mildly with the mesh (148 -> 351;
  6x6x4: 270; 12x12x8: lambda_min 2.19e-4, cond about 4.6e2 if lambda_max is about 0.1 as on the other meshes;
  lambda_max was not computed there). In double precision the expected relative error of N is
  of order cond^2 eps ~ 1e-11: **not floating-point limited**.
* **Regularisation.** With b orthogonal to e and ker C = span(e), any gamma > 0 gives C^+ b. On 9x9x6 b is *not*
  exactly orthogonal to e: b_z . e_hat = 4.7e-9 |b_z| (x: 6.9e-12), from the plane q_z = 1/2 of the even c-mesh, where
  the bands stick together in pairs with group velocities +-v_z (6_3 screw axis plus time reversal) and the pair
  cancellation in sum (d f)^2 v holds only to about 1e-9 (`na_bdote.py`: the whole residual 1.59e-2 of
  sum (d f)^2 v_z comes from that plane; |v_z| there up to 91 THz A). The spurious term (b . e_hat)^2 / gamma^2 in N
  is 1e-18 relative at the campaign's gamma = tr C / n, but it biases tau_z by **0.44 %** at gamma = 1e-8 tr C / n.
  The campaign's choice is safe; projecting b onto e^perp before the solve (as done here) removes the issue for any
  gamma. 5x5x3 and 7x7x5 (odd c-meshes) have b . e_hat <= 5e-17.

### 2.2 Where N comes from (lowest eigenpairs, mode frequencies)

N = sum_k (q_k . b)^2 / lambda_k^2 weights the slow end of the spectrum quadratically. Shares of N (eigenpairs
ordered by eigenvalue; dense eigendecomposition up to 9x9x6, shift-invert Lanczos with the Cholesky factor at
12x12x8):

| mesh | N_x: lowest 1 / 5 / 20 | largest single pair (x) | pairs for 50 % (x) | N_z: lowest 1 / 5 / 20 | largest single pair (z) | pairs for 50 % (z) |
|---|---|---|---|---|---|---|
| 5x5x3 | 0.00 / 0.52 / 0.58 | 0.32 | 5 | 0.00 / 0.00 / 0.25 | 0.33 | 24 |
| 6x6x4 | 0.00 / 0.00 / 0.47 | — | — | 0.10 / 0.27 / 0.48 | — | — |
| 7x7x5 | 0.00 / 0.39 / 0.64 | 0.29 | 16 | 0.00 / 0.00 / 0.51 | 0.31 | 18 |
| 9x9x6 | 0.00 / 0.00 / 0.52 | 0.20 | 20 | 0.00 / 0.30 / 0.30 | **0.30** | 32 |
| 12x12x8 | 0.00 / 0.00 / 0.27 | 0.12 (within the lowest 24) | — | 0.00 / 0.22 / 0.35 | 0.22 | — |

("lowest 1" is zero whenever the slowest mode carries no current in that direction, e.g. the c-axis TA modes for x.)
The slow eigenvectors are localised on the **Gamma-adjacent acoustic modes** (participation ratio 4–14, weight on
acoustic bands >= 0.97): first in-plane shell for x, the four TA modes at q = (0, 0, +-1/N_c) for z (6x6x4 is the
exception: its four slowest eigenpairs sit on the band-sticking pairs at the A point q = (0, 0, 1/2), 5.16 THz,
participation ratio 2.0–2.2, and the lowest five carry 27 % of N_z — the zone-face structure of 2.5). Split of
N = |x|^2 over modes: the four lowest c-axis TA modes carry 24.7 % (5x5x3, 3.70 THz), 36.4 % (7x7x5, 2.30 THz),
34.6 % (9x9x6, 1.93 THz) of N_z; the first in-plane acoustic shell carries 64 % (7x7x5) and 41 % (9x9x6) of N_x;
at 12x12x8 the 16 modes below 2 THz carry 25 % (x) and 24 % (z). Fewer than 0.5 % of the modes carry 25–65 % of N.
These modes, their frequencies (3.70 / 2.85 / 2.30 / 1.93 / 1.46 THz) and their rates move with the mesh:
lambda_min(H) = 9.76 / 4.71 / 4.04 / 2.88 / 2.19 x 1e-4 (fit lambda_min ~ f_min^1.55).

### 2.3 Stability of tau with respect to the rates (elasticities)
e_alpha = d log tau / d log g_alpha (exact formula with u = C^+ x; sum e_alpha = -1 checked to 3e-14):
sum |e_alpha| = 1.13 / 1.22 (5x5x3, x / z), 1.31 / 1.35 (7x7x5), 1.39 / 1.53 (9x9x6); largest single-orbit
|elasticity| 0.03–0.06. A relative error delta in every rate changes tau by at most about 1.5 delta. Events with rates
below 1e-3 x median (8.3 / 6.0 / 3.9 % of the events) carry a total |elasticity| of 8e-7 / 3e-6 / 3e-6: **events
with tiny rates are irrelevant**. The dominance of a few slow eigenpairs is a property of the operator, not a
numerical fragility.

### 2.4 Event cut and smearing at fixed mesh
* Cut (rates set to zero beyond n sigma): tau(3.5 sigma) / tau(4 sigma) - 1 = +0.03 / +0.06 % (5x5x3), +0.04 / +0.05 %
  (7x7x5), +0.05 / +0.04 % (9x9x6); 3 sigma: +0.2 … +0.4 %. The 4-sigma truncation is **CONVERGED** (residual
  < 0.01 % by extrapolation). Not the source of the mesh dependence.
* Smearing (exact Gaussian re-weighting of the stored events to sigma' <= 0.1 THz with a 4 sigma' cut):

| mesh | tau_x at sigma' = 0.05 / 0.07 / 0.10 THz | tau_z at 0.05 / 0.07 / 0.10 THz |
|---|---|---|
| 5x5x3 | 40.43 / 35.97 / 34.25 (-15 %) | 30.51 / 28.59 / 27.45 (-10 %) |
| 7x7x5 | 84.45 / 64.43 / 52.61 (**-38 %**) | 49.67 / 49.94 / 50.62 |
| 9x9x6 | 64.69 / 62.87 / 61.40 (-5 %) | 58.67 / 58.17 / 59.15 |
| 12x12x8 | 74.80 / 74.73 / 73.57 (-1.6 %) | 71.20 / 70.77 / 70.05 (-1.6 %) |

  5x5x3 and 7x7x5 are strongly under-resolved in the smearing (the delta functions are sampled by few near-resonant
  triplets); at 12x12x8 the smearing dependence has fallen to 1.6 % while the mesh dependence has not (2.6).

### 2.5 Convention: orbit average versus gauge-invariant block average (`na_blockavg.py`, `na_bigmesh.py ... block`)
The campaign's reported numbers use the orbit average of per-event rates, which is gauge dependent for events that
touch degenerate bands (09 item 2). The block average (sum of pp over ordered band triples of the degenerate blocks
divided by the block size; missing partners count as zero) was rebuilt from the pp files; it reproduces the
campaign's own block-average reference at 7x7x5 to 3.6e-9 (tau) / 4.6e-9 (kappa). Block minus orbit, relative:

| mesh | tau_x | tau_z | kappa_xx | kappa_zz | degenerate modes |
|---|---|---|---|---|---|
| 5x5x3 | 1e-9 | -2.0e-4 | 4e-10 | -3.1e-4 | 22 |
| 6x6x4 | -1.0e-3 | **+1.5e-2** | -2.3e-3 | **-5.1e-2** | 502 (432 on q_z = 1/2) |
| 7x7x5 | 1e-9 | +6.6e-5 | 1e-10 | -1.2e-4 | 38 |
| 9x9x6 | -1.7e-3 | **+1.06e-2** | -1.9e-3 | **-1.97e-2** | 1,090 (972 on q_z = 1/2) |
| 12x12x8 | -2.2e-4 | +5.7e-3 | -1.0e-3 | -8.2e-3 | 1,894 (1,728 on q_z = 1/2) |

On even c-meshes the band-sticking plane makes the orbit convention uncertain at the 1–5 % level for z quantities.
The S8 values at 9x9x6 therefore carry a convention uncertainty of about 1 % (tau_z 59.15 vs 59.78 ps) and 2 %
(kappa_zz 256.9 vs 251.8). Witness *ratios* are insensitive (section 3.3). The permutation of the velocity
eigenvalues among degenerate partners changes tau by <= 2e-10 (`na_gauge.py`; the random-rotation variants in
`gauge_m*.json` are not gauge transformations and were used only as a sensitivity diagnostic).

### 2.6 Mesh refinement (`na_bigmesh.py`, `na_refine.py`; events exported with the campaign's `aln_events.py`)

| mesh | modes | events | f_min (THz) | tau_x / tau_z (ps), orbit | block | kappa_xx / kappa_zz, orbit | N/N_q (x) | K/N_q (x) |
|---|---|---|---|---|---|---|---|---|
| 5x5x3 | 897 | 44,406 | 3.70 | 34.25 / 27.45 | 34.25 / 27.44 | 214.3 / 202.0 | 1.37e10 | 3.19e7 |
| 6x6x4 | 1,725 | 159,544 | 2.85 | 45.06 / 44.40 | 45.02 / 45.08 | 241.6 / 234.3 | 2.03e10 | 3.59e7 |
| 7x7x5 | 2,937 | 481,899 | 2.30 | 52.61 / 50.62 | 52.61 / 50.63 | 243.0 / 272.3 | 2.39e10 | 3.61e7 |
| 9x9x6 | 5,829 | 1,905,484 | 1.93 | 61.40 / 59.15 | 61.29 / 59.78 | 265.8 / 256.9 | 3.05e10 | 3.95e7 |
| 12x12x8 | 13,821 | 10,733,499 | 1.46 | **73.57 / 70.05** | 73.56 / 70.46 | 280.4 / 267.4 | 3.86e10 | 4.17e7 |

The lean solver reproduces the campaign at 9x9x6 to 1e-14 and the dense route at 6x6x4 to 1e-14.
* **Observed order** for the similar-aspect family 6x6x4 / 9x9x6 / 12x12x8 (h = 1/N, Q = Q_inf + C h^p): the ratio of
  successive differences is 1.342 (tau_x) / 1.353 (tau_z) in the orbit convention and 1.327 / 1.377 in the block
  convention, below the value ln(3/2)/ln(4/3) = 1.409 that any convergent power law (p > 0) must exceed; p = -0.14 /
  -0.12 (orbit), -0.17 / -0.07 (block). tau grows roughly like N^0.15 with no sign of saturation; **Richardson
  extrapolation is not applicable** and no converged value can be quoted. kappa behaves differently: p = 0.4 (xx),
  1.2–1.3 (zz); its extrapolations (about 390 / 293 W/mK) are themselves unreliable at such low orders.
* **Infrared structure** (`na_ir.py`). N per q-point grows with the mesh (x 1.37 -> 3.86 e10) while K per q-point grows
  slowly (3.19 -> 4.17 e7). Within-mesh fits of the Gamma-adjacent acoustic rates, r = diag C ~ f^p for f < 4 THz:
  p = 2.01 (7x7x5), 1.95 (9x9x6), 1.61 (12x12x8). HEURISTIC: with b ~ kT v near Gamma (equipartition) and r ~ q^p,
  the continuum RTA analogue N ~ int q^2 dq q^-2p diverges at small q when p >= 3/2, as q_min^(3 - 2p), i.e. as
  (mesh division)^(2p - 3) ~ N^0.2 for p = 1.6 — the same order as the observed growth of tau. The RTA analogue
  itself grows 30.1 / 40.4 / 49.2 / 59.0 / 72.7 ps (x). Whether tau_mem has a finite limit at sigma = 0.1 THz is
  therefore **not decidable from these meshes**: if the low-frequency exponent stays >= 3/2 the memory time is
  infrared-divergent (sigma > 0 cuts this off only when f_min ~ sigma, i.e. at about 100 divisions along c); if it
  drops below 3/2 the limit exists but is approached very slowly.

### 2.7 Verdict S8
The tabulated values are the correct values of the stated discrete operators (RELIABLE; caveat: 1–2 % convention
uncertainty in the z values at 9x9x6). The mesh dependence is a **genuine discretisation effect** — it is not due to
the regularisation (invariant to 1e-13 at the campaign's gamma), to ill-conditioning (cond <= 351), to round-off, or to
the 4-sigma cut (<= 0.06 %). It is carried by the Gamma-adjacent slow acoustic modes, whose frequencies and rates are
set by the mesh. The report understates it: tau_mem rises a further 20 % / 18 % from 9x9x6 to 12x12x8 (+115 % /
+155 % from 5x5x3), and the refinement sequence shows no convergent regime. "+79 % / +116 % from 5x5x3 to 9x9x6" is a
lower bound on the discretisation error of the 5x5x3 value, not an estimate of it.

---

## 3. S5 — witness tensors on other meshes

### 3.1 Reproduction (`na_transfer.py`; `na_out/transfer_m*.json`)
Every stored tau ratio, kappa ratio and degeneracy-averaged lifetime extremum of `mesh_transfer.json` (8 tensors,
7x7x5), `verify_m775.json` (4), `verify_m553.json` (4), `verify_m996.json` (4) and `transfer_controls.json` (2) is
reproduced by the independent operator: max |difference| of the tau ratios 1.0e-13 (7x7x5), 6.8e-14 (9x9x6),
7.5e-15 (5x5x3); lifetime extrema identical to all printed digits.

### 3.2 Conditioning and near-null structure of the witness operators
* cond(H) = 137–457 at 7x7x5, 318–455 at 9x9x6; lambda_min(H) >= 1.87e-4 for every tensor; null eigenvalue <= 1e-17.
  No witness creates a spurious near-null mode; the ratios are not floating-point or regularisation limited.
* But several ratios are carried by **one** slow eigenpair. On their objective mesh (7x7x5) the E20 witnesses lower
  lambda_min(H) to 0.46 (cal1nz max z), 0.63 (cal1nz max x), 0.72 / 0.80 (orb 0.1 z / x) of the reference, and a single
  eigenpair carries 88 % of N_z (cal1nz z; reference 31 %), 58 % of N_x (cal1nz x; reference 29 %), 71 % (orb z),
  43 % (orb x). On 9x9x6 the same tensors have lambda_min ratios 0.85–1.04 and no such concentration: the optimiser
  engineered a slow mode among the few Gamma-adjacent acoustic modes of the objective mesh. This is the mechanism of
  the sampling artefact, and it is resolved accurately (it is a physical slow mode of that discrete operator, well
  separated from zero).

### 3.3 Robustness of the ratios
* Cut 3.5 sigma vs 4 sigma: ratio changes <= 0.1 % (largest: cal1nz max z on 7x7x5, 3.0534 -> 3.0559); 3 sigma <= 0.2 %.
* Convention (block minus orbit, ratios): <= 1.7e-4 (7x7x5), <= 1.8e-3 (9x9x6), <= 1.4e-4 (5x5x3); every kappa
  violation of E19/E21 is unchanged in the gauge-invariant convention (e.g. cal1nz max z on 9x9x6: kappa_zz 0.918 in
  both conventions).
* Smearing at fixed mesh: on 7x7x5 some ratios move strongly (cal1nz max x 1.80 -> 1.92 -> 2.17 at sigma' = 0.10 /
  0.09 / 0.08 THz; cal1 max x 0.83 -> 0.79 -> 0.74; orb 0.1 J max x 1.33 -> 1.43), z ratios little (3.05 -> 3.04;
  1.69 -> 1.69). On 9x9x6 all ratios move by <= 3 %. The conclusion "the effect does not transfer" holds for every
  sigma' in [0.08, 0.10] THz and every cut in [3, 4] sigma.

### 3.4 Witnesses never evaluated on 9x9x6 (`na_untested.py`, `na_untested_variants.py`)
E21 records that the 9x9x6 amplitudes of the E15 per-cluster max z witness were not computed. They were computed here
(phono3py 4.5.0 Rust kernels, `p3_events.py`, campaign tensors `results/m553/witness/fc3_orb0.1_c*_max.npy` and the
review tensors `review/work/fc3_c1_*_cont02.npy`; the recomputed 9x9x6 reference pp is bit-identical to the
campaign's). Ratios to the 9x9x6 reference, independent operator; the campaign's own `mesh_transfer_check.py` gives
the same numbers to 1e-13:

| tensor (5x5x3 ratio; 7x7x5 ratio) | tau_x | tau_z | kappa_xx | kappa_zz | lifetimes (deg.-avg.) | tau_z at sigma' 0.08 / cut 3 sigma / block conv. |
|---|---|---|---|---|---|---|
| E15 orb 0.1 max z (1.370; 1.347) | 1.025 | **1.518** | 1.025 | 1.046 | -35 % … +20 % | 1.497 / 1.518 / 1.520 |
| review c1_z_cont02 (1.436; 1.315) | 1.044 | **1.599** | 1.038 | 1.058 | -37 % … +21 % | 1.577 / 1.599 / 1.601 |
| E15 orb 0.1 max x (1.308; 1.087) | **1.157** | 0.980 | 1.018 | 0.968 | -31 % … +22 % | — |
| review c1_x_cont02 (1.387; 1.048) | **1.118** | 1.058 | 1.012 | 0.977 | -39 % … +22 % | — |

The c-axis enhancement of the per-cluster eps = 0.1 class is **not** a sampling artefact: it persists on 7x7x5 and
grows on 9x9x6 (1.37 -> 1.35 -> 1.52; 1.44 -> 1.32 -> 1.60). Its mechanism is visible in the spectrum: lambda_min(H)
falls to 0.66 / 0.63 of the reference on 9x9x6 and the lowest five eigenpairs carry 53–54 % of N_z (reference 30 %) —
these tensors slow down the long-wavelength c-axis TA modes, which every mesh samples. What fails on the finer meshes
is the data: kappa_zz +4.6 / +5.8 % and lifetimes down to -35 / -37 % (the ALAMODE tensor: kappa within 0.6 %,
lifetimes -28 … +19 % on 9x9x6).

### 3.5 Verdict S5
The numbers of E19 and E21 are RELIABLE. The statement S5 as worded ("every single-mesh witness tested on a finer
mesh loses or reverses its memory change") is NOT RELIABLE: it is contradicted by E19's own row for orb 0.1 max z and
by 3.4. S6 ("largest tau change of any witness on an unused mesh is 1.18 x") is contradicted (1.35 on 7x7x5 in E19;
1.52 and 1.60 on 9x9x6). Section 1 of the report ("All large changes found (factors 1.3–6) are artefacts of the
sampling on which they were optimised") holds for the calibrated-class and E20 witnesses but not for the per-cluster
c-axis witnesses. What survives without exception is: **every witness violates the static data on a mesh not used
in its construction.** The headline sentence of section 1 ("no fc3 change … moves tau_mem by more than about 20 % in a
way that survives a change of Brillouin-zone sampling *while keeping the static data*") remains literally true, but
only because of the data clause; the tau change itself does survive for the c-axis per-cluster witnesses.

---

## 4. S1 — exact linearisation at 7x7x5 (`na_linear.py`; `na_out/linear_m775.json`, `linear_m553.json`)

* **Phi0 = G c0.** The unscaled parts lin_G_*.npy are no longer on disk; the check was made on the stored scaled
  matrices G_s = G M^-1/2 that E20 actually used: G_s (M^-1/2)^-1 c0 vs Phi0 (solve residual <= 4e-15). Row-wise max
  relative error 1.65e-10 (cal1), 1.59e-10 (cal1nz), 1.69e-10 (orb) — the stated 1.6e-10 — attained on rows with
  small |Phi0|; relative to max |Phi0| 1.0–1.1e-13; 99.9 % quantile 6.7–7.2e-12. 5x5x3: 1.2–1.4e-11 row-wise.
* **Amplitude level at the witnesses** (sharper than tau ratios): uc |Phi_row(u)|^2 of the linear model vs phono3py's
  Rust pp of the rebuilt tensor at the same representative event, all 20,141 non-degenerate rows, 8 tensors:
  max relative difference 1.5e-11 for pp > 1e-6 max (5.1e-11 over all rows); absolute difference / max pp 1.8e-14.
  Representative triplets and bands match for 100 % of the rows. The optimiser history values are reproduced to 3e-15.
* **tau / kappa ratios** model (block convention) minus phono3py route (orbit convention):
  cal1 witnesses (E19): x <= 6.6e-7, z **-6.0e-6, -5.59e-5, -3.82e-5, -3.86e-5** (relative up to 7.2e-5) — the stated
  "<= 5e-5" is exceeded slightly (5.6e-5); E20 witnesses at 7x7x5: z up to **-1.74e-4** (cal1nz max z); at 5x5x3: z up
  to -1.43e-4. Block-averaging the phono3py pp reproduces every one of these differences to 1e-9, so they are the
  convention, not a model error.
* E20 "tau equal to the optimiser's to <= 1e-5": tau in ps, phono3py route vs optimiser: x <= 1.7e-7 everywhere;
  objective z of cal1nz max z 9.0e-6, but orb 0.1 max z 7.4e-5 (7x7x5) and all z values on 5x5x3 2.2–3.3e-4.

**Verdict S1: RELIABLE WITH CAVEATS** — the linearisation is exact to 1e-11 at the amplitude level; the quoted
agreement bounds should read "<= 6e-5 (cal1, E19), <= 1.8e-4 (E20), block- versus orbit-average convention".

---

## 5. Optimiser: inner solves at maxiter (`na_kkt.py`)

* Histories (`results/joint/*.json`): all 60 inner solves stopped at nit = 300; feasibility from outer iteration
  4–7 on; over the last six outer iterations tau varies by 1.7e-5 – 1.1e-4 (relative); final violation 6e-8 – 8e-7;
  penalty mu up to 3.3e6. The plateau alone does not prove stationarity (review 3.1).
* **First-order optimality** at the reported best iterates, with the campaign's linearised models (memory-mapped):
  KKT residual rho = min_{lambda >= 0} ||grad log tau - J_act^T lambda|| / ||grad log tau|| over the constraints within
  2e-3 of their bound.

| run | active (lifetimes 5x5x3 / 7x7x5, orbit bounds) | rho | multiplier mass: 7x7x5 lifetimes / kappa / orbit bounds |
|---|---|---|---|
| cal1nz max x (1.799) | 17 / 35, 158 | 0.020 | 0.44 / 0.15 / 0.37 |
| cal1nz max z (3.053) | 20 / 33, 174 | 0.007 | 0.48 / 0.24 / 0.21 |
| orb 0.1 max x (1.326) | 16 / 23, 113 | 0.017 | 0.53 / 0.17 / 0.17 |
| orb 0.1 max z (1.690) | 22 / 22, 123 | 0.017 | 0.59 / 0.16 / 0.11 |

* Calibration with the same optimiser on 5x5x3 (per-cluster eps = 0.1, max tau_x): E15 iterate rho = 0.042 (ratio
  1.287), E16 best start 0.035 (1.308), review continuation 0.006 (1.387), review polish 0.0017 (1.398). A residual of
  3.5–4 % went with a 7–9 % shortfall to the later local value; 0.6 % with 0.8 %.
* **Effect on the conclusions: none.** (i) Every E19–E21 statement is about specific tensors, evaluated after the
  fact (sections 3.1, 4). (ii) The E20 maxima are within a few per cent of the local maxima of their basins (rho
  0.7–2 %); they are inner bounds, as stated, and global optimality is unknown (one start per run). (iii) The only
  conclusion an incomplete search can weaken is the absence statement S6 / section 1, and that statement is limited
  more by what was searched (no run maximised a ratio on a mesh while imposing the data on every mesh used for
  validation) and by section 3.4 than by the iteration cap.

---

## 6. Further observations

1. **Symmetry closure of exported event sets.** At 12x12x8 the exported list misses 88 symmetry images (44 events
   lack an image; all have rates <= 2.7e-15 x median, total share 2.5e-22): selection-rule zeros that phono3py returns
   as tiny positive numbers for one member of an orbit and exactly zero for another. Harmless numerically, but
   `symmetry.orbits` (FC3Class) raises "event set not invariant under symmetry" on such a list; a convergence study
   with the campaign code needs a tolerance (or treats missing images as zero rate, as done here).
2. **Even c-meshes** (6x6x4, 9x9x6, 12x12x8) contain the band-sticking plane q_z = 1/2; there the orbit convention is
   uncertain at the 1–5 % level for z quantities (2.5) and b is orthogonal to e only to 1e-9 (2.1). A mesh study
   should use the block (gauge-invariant) convention and project b onto e^perp.
3. **E14** (read, not re-run): rank 218 / 218 in the 3.5 A subspace with sigma_min = 0.118 against sigma_max = 223
   (condition 1.9e3) is an unambiguous numerical rank; the full-space values (sigma_min 3e-13) are rank deficient as
   stated.

---

## 7. Corrections suggested for the report

* S5: "Every single-mesh witness tested on a finer mesh violates the data there; the calibrated-class witnesses also
  lose or reverse their memory change, whereas the per-cluster c-axis enhancement persists (1.35 on 7x7x5; 1.52 and
  1.60 on 9x9x6 for E15 orb 0.1 max z and its continuation) with kappa_zz +3 … +6 % and lifetimes down to -37 %."
* S6: replace "1.18 x" by "1.60 x (c1_z_cont02 on 9x9x6), with the data violated there"; restrict "1.18" to the E20
  two-mesh witnesses.
* Section 1 / 4: "All large changes found are artefacts of the sampling" -> "all large changes violate the static data
  on a mesh not used in their construction; the calibrated-class and two-mesh effects are in addition sampling
  artefacts".
* S8: add the 6x6x4 / 12x12x8 values, the observed orders (-0.07 … -0.17: no convergent regime) and the infrared
  anatomy; state the 1–2 % convention uncertainty at 9x9x6; "+79 % / +116 %" is a lower bound on the discretisation
  error of the 5x5x3 value. "Must be quoted with a convergence study" should mention that tau_mem may be
  infrared-divergent (2.6) — a convergence study may not converge.
* S1 / E19 / E20: "<= 5e-5" -> "<= 6e-5 (E19), <= 1.8e-4 (E20)"; "tau equal to the optimiser's to <= 1e-5" -> "x:
  <= 2e-7; z: <= 3.3e-4 (block versus orbit convention)".

---

## 8. Files and reproduction

Scripts (`review/`): `na_common.py` (independent operator), `na_spectral.py` (2.1–2.4), `na_bdote.py` (2.1),
`na_transfer.py` (3.1–3.3), `na_linear.py` (4), `na_gauge.py` (2.5), `na_blockavg.py` (2.5, 3.3), `na_bigmesh.py`
(2.6; options: cut, block convention, smearing), `na_check664.py` (dense cross-check), `na_ir.py` (2.6),
`na_refine.py` (2.6), `na_kkt.py` (5), `na_closure.py` (6.1), `na_untested.py`, `na_untested_variants.py` (3.4).
Outputs: `review/na_out/*.json` (one file per script and mesh; event-list metadata `events_m664_s0.1.meta.json`,
`events_m12128_s0.1.meta.json`; `untested_m996_*.json`).

```
<research env> -B na_spectral.py m553|m775|m996            <research env> -B na_transfer.py m553|m775|m996
<research env> -B na_linear.py m775|m553                    <research env> -B na_gauge.py m553|m775|m996
<research env> -B na_blockavg.py m553|m775|m996             <research env> -B na_refine.py
<phono3py env> -B scripts/aln_events.py --mesh 6 6 4 | 12 12 8 --sigma 0.1 --nsig 4 --T 300 --out <tmp>/events_m664_s0.1.npz | events_m12128_s0.1.npz
<research env> -B na_bigmesh.py <events> <out.json> [none|nsig] [orbit|block] [sigma']
<research env> -B na_ir.py <out.json> <events> ...
<research env> -B na_kkt.py cal1nz_x_max | cal1nz_z_max | orb01_x_max | orb01_z_max
<research env> -B na_kkt.py custom <u.npy> orb 0.1 0 m553 m553 <out.json>
<phono3py env> -B scripts/p3_events.py --mesh 9 9 6 --events results/m996/events_m996_s0.1.npz --outdir <tmp>/m996w --fc3 REF --fc3 orb0.1_c2_max=results/m553/witness/fc3_orb0.1_c2_max.npy --fc3 orb0.1_c0_max=... --fc3 c1_z_cont02=review/work/fc3_c1_z_cont02.npy --fc3 c1_x_cont02=...
<research env> -B na_untested.py ; na_untested_variants.py ; scripts/mesh_transfer_check.py (campaign route cross-check)
```
Peak memory of any single run: about 2.5 GB (12x12x8 dense operator, 1.53 GB, factorised in place).
