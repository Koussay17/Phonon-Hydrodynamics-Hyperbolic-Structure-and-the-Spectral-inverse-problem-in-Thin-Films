# 07 — Experiments (all scripts in `scripts/`, raw data in `results/`)

Environment: research Python 3.14.7 (numpy 2.5.3, scipy 1.18.1, python-flint 0.9.0, cvxpy 1.9.3
with Clarabel); phono3py 4.5.0 environment for the AlN export only. Machine: 6 cores, 28 GB.
Units in the Debye models: lattice constant 1, energies in units of c_T, rates in model units;
tau in the inverse rate unit. AlN: tau in ps (tau[ps] = (N/K)/(4 pi), operator in THz).
Revised after the proof audit (`review/proof-audit.md`); certified lower bounds are rounded down.

## E1. Exact Ben-Tal–Teboulle constant, 1D N = 9 (`bt_exact.py`, `results/bt_exact_d1_N9.json`)
n = 16 modes, m = 26 exact events, d = 15. All C(26,14) = 9,657,700 subsets enumerated (331 s);
3,566,333 bases. Condition numbers separate cleanly: nonsingular subsets log10 cond in [1, 8),
singular ones >= 16. **M = 16284.0067**; second method (rational kernel of the integer
stoichiometric matrix, |y_I| = |P_H D phi_I| / |b^T D phi_I| at 50 digits): 16284.00672356293
(float value agrees to 5e-11). K0 = 17541.29 -> M^2 K0 = 4.65e12, versus tau_ref = 169.1 and the
lifetime-constrained tau_max(found) = 627.5 (unconstrained rates; 509.8 with time reversal).
The condition-number cut used here to declare subsets singular is heuristic. The value was
certified independently in the proof audit (`review/certify_M_1d9.py`): an integer-rank bound on
S_I plus a time-reversal certificate classify all 9,657,700 subsets rigorously (3,566,333 bases,
the same count; 8,430 rank-14 subsets with b^T D phi = 0 are exactly singular). The exact
recomputation with energies exactly 2 pi k/N gives M = 16284.006723562 (3e-14 relative to the
float-energy value above). Next hyperplane: 6725.96.

## E2. Lower bounds on M by basis pivoting (`memgap.bt_local_search`, `results/bt_local_search.json`)
Each value is an explicit basic solution, hence a valid lower bound on M; the bases are genuine
(audit `review/check_pivot_bases.py`: cond(B_I) = 2e5–3e9, 50-digit re-solves agree to <= 1.2e-8).
The search is strongly seed-dependent: with this seed it recovers the exact 1D N = 9 value, another
seed stops at 6725.96, and a re-run for 3D N = 3 reaches 1.77e6 (above the 2.05e5 listed here).
These numbers are lower bounds, not estimates of M; the json stores neither the subsets nor their
conditioning. M >= 2.3e2 (1D N=15), 4.1e2 (N=21), 3.6e4 (N=27), 2.3e4 (2D N=3), 2.6e5 (2D N=5),
2.0e5 (3D N=3). sqrt(tau_ref/K0) = 0.06–0.10 in all cases. Resulting vacuity factor
M^2 K0/tau_ref (audit `check_bt_numbers.json`): >= 7.88e6 (1D N=15), >= 3.15e7 (1D N=21), and
>= 2.7e10 in the other geometries (1D N=9: 2.75e10).

## E3. Inner tau interval, Debye event sets (`scan.py`, `results/scan/`, `results/scan_summary.json`)
Reference rates: Klemens-type (02-assumptions B). Fixed: diag C = r, K(0) = K0 (x current).
Multi-start (per configuration 4 seeds x 13 starts: reference, lognormal sigma 0.7/1.5/3, LP
vertices of the lifetime polytope mixed toward the K-minimiser) x both signs, multiplicative
projected gradient with exact Newton restoration (feasibility <= 1e-12). tie: none / tr (time
reversal) / full (cubic point group incl. inversion).

| config | n | m | tau_RTA | tau_ref | tau_min (found) | tau_max >= | max/min >= |
|---|---|---|---|---|---|---|---|
| 1D N=9 tr | 16 | 26 | 113.3 | 169.1 | 104.9 | 509 | 4.86 |
| 1D N=21 tr | 40 | 124 | 673.4 | 725.7 | 210.6 | 1.05e6 | 4.99e3 |
| 1D N=33 tr | 64 | 294 | 1711 | 1497 | 351.9 | 6.05e6 | 1.72e4 |
| 1D N=63 tr | 124 | 1032 | 6409 | 4720 | 659.8 | 1.20e8 | 1.82e5 |
| 2D N=3 full | 16 | 32 | 25.2 | 24.19 | 24.19 | 24.19 | 1 (slice is a point) |
| 2D N=5 full | 48 | 196 | 65.8 | 89.8 | 59.04 | 4.36e3 | 74.0 |
| 2D N=7 full | 96 | 544 | 88.1 | 128.9 | 63.81 | 3.60e4 | 565 |
| 2D N=9 full | 160 | 1284 | 104.9 | 135.8 | 59.58 | 1.89e5 | 3.18e3 |
| 3D N=3 full | 52 | 290 | 33.6 | 35.3 | 28.72 | 67.0 | 2.33 |
| 3D N=5 full | 248 | 4050 | 72.6 | 88.7 | 45.77 | 1.61e5 | 3.52e3 |
| 3D N=5 none | 248 | 4050 | 72.6 | 88.7 | 38.20 | 2.12e6 | 5.55e4 |
Lower bounds (tau_max, ratio) are rounded down; found minima are rounded up. Full table (32
configurations): `results/scan_summary.json`. The best maximum is typically found by 1 of ~50
starts (column "hits"): tau_max values are lower bounds; tau_min values are reproduced by most starts.

## E4. Independent verification of witnesses (`verify_witness.py`, `results/verify_debye_witnesses.json`)
12 witnesses (min and max of 1D N=21/27 tr, 2D N=5 full/tr, 3D N=3 full/tr): Cholesky,
eigendecomposition on H, matrix-free conjugate gradients, and arb ball arithmetic (128 bits).
Agreement >= 9 digits; arb radii < 1e-20; diagonal residual <= 1e-12; K residual <= 1e-12;
condition numbers on H up to 4.3e6. Every witness is admissible (the certified solve of
C + gamma e_hat e_hat^T proves ker C = span(e)). The untied 1D N = 9 maximum (627.5) is not among
these 12 witnesses.

## E5. Lower end of the interval vs the Cauchy–Schwarz bound tau_CS = K0/|b|^2
tau_min(found)/tau_CS (from `results/scan_summary.json`; tie in brackets): 1.648 (1D N=9 tr),
1.652 (1D N=15 tr), 1.261 (N=21 tr), 1.188 (N=33 tr), 1.103 (N=63 tr); 1.684 (2D N=5 full),
1.093 (2D N=5 tr), 1.318 (2D N=7 full), 1.200 (2D N=9 full); 1.206 (3D N=3 full), 1.012 (3D N=3
none), 1.242 (3D N=5 full), 1.036 (3D N=5 none); 1.00034 (AlN x), 1.00001 (AlN z). Maximum over
all 32 stored configurations: 1.6839 (2D N=5, full symmetry). The proved lower bound is nearly
attained: the minimiser makes b nearly an eigenvector of C. The rigorous single-mode bound
(P4(b) in 13) is < 1 in model units and is not useful at these sizes.

## E6. Convex side: minimum DC response over the lifetime polytope (`results/K_lower_sdp.json`)
K_lo(r) = min over P(r) of the extended-value K (13, P5) = infimum over admissible g: SDP with the
Schur-complement LMI [[C(g), b],[b^T, t]] >= 0 (Clarabel) versus mirror-descent minimiser (convex
problem): 13309.64 / 13309.64 (1D N=9), 2524.886 / 2524.887 (2D N=3), 12164.53 / 12164.55 (2D N=5),
4606.77 / 4606.78 (3D N=3); relative differences 1.0e-7, 2.3e-7, 2.0e-6, 3.3e-6, i.e. agreement to
<= 3.3e-6 where the SDP converged. 1D N=15: SDP "optimal_inaccurate", 45602 vs 45790 (4.1e-3).
In all cases K_RTA/3 < K_lo < K0 (L3).

## E7. Infrared scaling without optimisation (`ir_scaling.py`, `results/ir_scaling_d{1,2,3}.json`)
Prescribed r = A eps^alpha, reference = KL projection of Klemens rates. Local exponents of tau
between the two finest meshes:
| d | alpha | meshes | tau_RTA exponent | tau_ref exponent | RTA prediction |
|---|---|---|---|---|---|
| 2 | 0.5 | 5..25 | 0.07 | 1.01 | converges |
| 2 | 1.0 | 5..25 | 0.33 | 0.69 | log N |
| 2 | 1.5 | 5..25 | 0.95 | 0.65 | N^1 |
| 3 | 1.0 | 3..11 | 0.33 | 0.32 | converges |
| 3 | 1.5 | 3..11 | 0.69 | 0.61 | log N |
| 3 | 2.0 | 3..11 | 1.29 | 1.16 | N^1 |
tau_RTA follows the moment criterion qualitatively (converging / log / power); the event operator
tracks tau_RTA in 3D for alpha >= 1.5 (tau_ref/tau_RTA = 1.04–1.15 at N >= 5) but exceeds it at
small alpha. K0/N^d still drifts (exponents 0.6–1.1): the meshes are pre-asymptotic, so these
numbers SUPPORT but do not confirm the exponent law for the RTA sums, and they say nothing
two-sided about the event operator (its mesh exponents are unjustified; 13, T7(ii)). The exact
RTA lattice sums of the audit (`review/check_ir.py`, d = 3 up to N = 161) give the asymptotic
laws: convergence (alpha = 1), log N (alpha = 3/2), N^1 (alpha = 2), N^d/log N (alpha = d = 3).
In 1D, K0/N and tau_ref grow linearly even at alpha = 0.25 (vanishing umklapp measure; 09, item 7).

## E8. AlN event geometry, 5x5x3, sigma = 0.1 THz, 300 K (`aln_events.py`, `aln_geometry.py`)
897 modes, 44,406 events (56 % umklapp, 178 repeated daughters), 2627 event orbits and 117 mode
orbits under 6mm x time reversal (24 operations). Validation vs phono3py's symmetrised Omega'
(local copy of the baseline check): kappa_xx 214.30 vs 214.28, kappa_zz 201.97 vs 202.82 W/(m K);
tau_mem,x 34.25 vs 34.21 ps, tau_mem,z 27.45 vs 27.99 ps; energy residual 4.5e-16. The z values of
this comparison used the symmetrised Omega' without degeneracy averaging; the counterexample
review's degeneracy-averaged rebuild gives kappa_zz 202.18 W/(m K) and tau_z 27.54 ps (x unchanged:
214.28, 34.21 ps), closer to the event operator; C^+ b differs from Omega' by 0.15 % (x) and
0.79 % (z). Reference at 7x7x5 (review, same export script): kappa 243.03 / 272.33 W/(m K),
tau_ref 52.61 / 50.62 ps, tau_CS 21.10 / 20.02 ps, lambda_min 4.04e-4 THz. All AlN numbers below
are 5x5x3 numbers.
Selection-rule zeros (review `cx_forbidden_anatomy.json`, `cx_selection_rules.json`): 3,275 events
(7.38 %, 335 orbits) have reference rates <= 4.35e-15 of the maximum (allowed events >= 3.07e-9);
95.6 % of them share a nontrivial unitary little group (17.2 % of allowed events): evidenced, not
proved event by event, to be space-group selection rules. 2,292 allowed orbits remain.
tau_RTA (same diagonal): 30.12 ps (x), 22.52 ps (z). tau_CS = K0/|b|^2: 19.1903 ps (x), 15.5168 ps (z).
Full row rank of the constraint maps at the reference (audit `check_dimension_count.py`): W 897
(sigma_min/sigma_max 1.96e-2), [W; grad K_ab] 903 (9.0e-6), W_gamma 117 (2.76e-2),
[W_gamma; grad K_xx; grad K_zz] 119 (5.9e-6): local dimensions 43,503 and 2,508 (13, section 2).

## E9. AlN inner interval (`aln_scan.py`, `results/aln/scan_m553_s0.1/`, `scan_m553_summary.json`)
Fixed: all 897 mode lifetimes, K_xx and K_zz (symmetric rates; K_yy = K_xx and off-diagonals
vanish by symmetry) or all six K_ab (unconstrained rates). 10 starts per (direction, sign) for
symmetric rates, 5 for unconstrained; up to 3000 steps.

| rates | direction | tau_min (ps) | tau_max >= (ps) | reference (ps) |
|---|---|---|---|---|
| symmetric | x (basal) | 19.1969 (10/10 starts) | 6906.7 (arb, well-conditioned); 7163.5 (float); 3.575e7 (near-singular, see below) | 34.25 |
| symmetric | z (c axis) | 15.5169 (10/10) | 959,689.9 | 27.45 |
| unconstrained | x | 19.1905 | 1.140e5 | 34.25 |
| unconstrained | z | 15.5171 | 1.201e5 | 27.45 |
Inner gap ratios (rounded down): x >= 359.78 (6906.7/19.197), z >= 61,848.
**Status after the counterexample review.** All maxima in this table (and the minima) activate
selection-rule-forbidden events (E8): the 6906.7 ps witness draws up to 99.998 % of some mode
lifetimes from them, the 959,689.9 ps witness up to 32.9 %. They are SUPERSEDED as AlN statements by
the review's witnesses with those events fixed at zero (`review/cx_runs/nofor_*`, verified
`cx_verify_newwitness.json`, certified with 256-bit arb and exact energy conservation
`cx_certify_new.json`): **tau_x = 7,038.17 ps** (arb radius ~1e-66, K_xx residual 3.4e-14, lifetimes
5.8e-13, lambda_min 2.8e-6, condition 4.2e4) and **tau_z = 4,322,183.47 ps** (radius < 1e-60, K_zz
residual 2.9e-11, lambda_min 1.1e-8, condition 1.1e7); factors over the reference >= 205 and
>= 1.57e5. The unconstrained rows break basal isotropy (tau_x = 114,019 ps with tau_y = 514 ps) and
are not AlN statements.
Transfer to the exact data (audit `review/check_witness_ift.py`): the x-max witness (6906.75 ps)
has residual 1.6e-14, a full-row-rank constraint Jacobian (119 rows, sigma_min = 0.0102), Newton
correction <= 1.6e-12 in log-rates and first-order |dtau/tau| <= 6e-13; the x-min witness: correction
<= 2.6e-11, |dtau/tau| <= 4e-12.

Verification (`results/aln/verify_m553.json`): z-max, z-min and x-min agree across Cholesky,
eigendecomposition and CG to <= 1e-10 (condition numbers <= 3e6). The x-max witness at 3.58e7 ps
has lambda_min on H = 2.5e-15 (condition 6.4e13); the three float evaluations agree only to 1e-3
(35.75e6, 35.73e6, 35.76e6 ps) and its K residual is 1.1e-10: it is reported only with the
ball-arithmetic result (`certify_aln.py`, 256-bit arb, `results/aln/certify_arb_m553.json`):
tau_x = 35,759,344.509 ps (radius < 1e-52 ps), K_xx/K0 - 1 = 1.1e-11, K_zz exact to 1e-16.
Because the witness sits next to a 2.5e-15 eigenvalue, tau may be sensitive at O(1) to the
1e-11 data mismatch, and the arb balls are built from the float stoichiometric coefficients
l_p, l_d (energy conservation to ~1e-16 only). The statement certified is therefore "tau_x = 3.58e7 ps
is attained for data within 1.1e-11 (relative) of the AlN data, for an event geometry whose energy
conservation holds to float precision". The robust certified bound uses the well-conditioned
witness full_c0_max_08 (arb): tau_x = 6906.7472 ps, K_xx/K0 - 1 = 1.1e-14 (the 7163.6 ps witness
_09 is float-verified only).

Mechanisms (`full_c0_*` witnesses): x-max suppresses 73 % of the events (rates below 1e-6 of the
reference; umklapp share among suppressed events equals the overall share) and creates an almost
exact hidden invariant in the 5–10 THz band (eigenvalue 2.5e-15 versus 9.8e-4 for the reference;
current overlap (b.u)^2/|b|^2 = 6.6e-19; 99.99 % of N). The well-conditioned 6906.7 ps witness
family makes low-frequency acoustic modes nearly conserved (83 % of the eigenvector weight below
5 THz, eigenvalue 2.5e-6, for the 7163.5 ps member). The minimiser cancels the soft-mode response:
the share of |x|^2 below 5 THz drops from 0.53 (reference) to 0.06; the ratio x_mu r_mu/b_mu on the
lowest x-coupled modes is 1.13–1.21 at the reference and 0.294 at the x-minimiser (audit
`check_ir.py`), so the RTA-like infrared hypothesis H4 is operator-specific.

## E10. Results of the independent counterexample review (`review/counterexample-review.md`)
Recorded here as experiments by an independent code base (`review/cx_*.py`); not re-run.
* **Factor-F prior (13, P7).** Rigorous ceiling tau <= F/lambda_min(C_ref|H) = 81.50 F ps (5x5x3),
  197.2 F ps (7x7x5). Box-constrained inner values (g_ref/F <= g <= F g_ref, symmetric rates,
  lifetimes and K_xx, K_zz fixed; `cx_verify_box.json`): F = 2: x [27.96, 45.63], z [20.71, 44.16] ps;
  F = 10: x [20.73, 153.48], z [15.72, 162.60] ps; F = 100: x [19.44, 480.12], z [15.517, 1421.2] ps.
* **Smooth amplitude prior** g = g_ref exp P(eps_p, eps_a, eps_b): degree <= 10 pins tau
  (34.2530 / 27.4465 ps); degrees 12–20, |P| up to 30: maxima found <= 1.80x the reference
  (43.19 ps x, 49.43 ps z).
* **Additional data fixed** (300 K lifetimes and K_xx, K_zz always fixed): kappa(T) at
  100–800 K -> tau_z >= 131,918 ps (gap persists); accumulation below 5/10/15 THz ->
  tau_z >= 61,254 ps (persists); K(z) at 1/ns, 1/(100 ps), 1/(10 ps) -> maxima found 62.9 ps (x),
  33.9 ps (z) (inner evidence, not a bound; slow modes of 0.1 % weight remain); lifetimes at five
  temperatures -> inconclusive (stiff constraints; first-order method crept to 35.9 / 28.4 ps).
* **L3 and the linewidth convention.** No counterexample with r = diag C (minimum of
  K s_max/K_RTA = 1.0000 over random, adversarial and property-based searches). With phono3py-type
  self-energy linewidths in place of diag C the bound fails for repeated-daughter networks: chain
  1, 2, 4, 8, 16 gives K/K_RTA(Gamma) ~ 0.3333 < 1/2 exactly; AlN: 1.2500 (x) vs 1.2572 with diag C.
* **Numerical artifacts:** none found (two solution routes; arb with exact energy conservation;
  1e-10 data perturbations move the robust witness by < 1e-9 relative).
