# 13 — Final report: identifiability gap of the memory moment in fixed event geometries

**Date:** 4 October 2026; revised 4–5 October 2026 after the proof audit (`review/proof-audit.md`)
and the independent counterexample review (`review/counterexample-review.md`).
Definitions and declared modelling choices: `02-assumptions.md`. Experiments: `07-experiments.md`.
Negative results and bugs: `09-failed-approaches.md`.

## 0. Declared class (every statement below is about this class)
* **Unknowns:** arbitrary nonnegative event rates g on the FIXED, energy-allowed three-phonon event
  set (momentum and energy selection; for AlN, rates symmetric under 6mm x time reversal). **No
  amplitude model**: no smoothness of |Phi|^2, no fc3 parametrisation, no per-event prior.
* **Data:** the 300 K mode lifetimes (r = diag C), the energy invariant, and the DC tensor
  (K_xx, K_zz for AlN).
* "Realistic event geometry" refers to the event SET only, not to the rates.
The gap reported here is a property of this class. It is not a statement about AlN's actual
memory; with a microscopic amplitude model the rates, and hence tau_mem, are determined.

## 1. Answer
**Within the declared class the gap is large and one-sided.**
* **Lower end essentially identified:** tau_mem >= tau_CS = K(0)/|b|^2 (proved); the minimum found
  lies within a factor <= 1.69 of tau_CS in every configuration tested (maximum 1.684, 2D N = 5,
  full cubic symmetry; 1.00034 for AlN x and 1.00001 for AlN z at 5x5x3).
* **Upper end not identified:** certified feasible witnesses exceed the reference memory by
  factors 10^2–10^5, and the inner ratio tau_max/tau_min grows with mesh size even with the full
  crystal symmetry imposed (Debye 2D: 1, >= 74.0, >= 565, >= 3.18e3 for N = 3, 5, 7, 9; 3D:
  >= 2.33, >= 3.52e3 for N = 3, 5).
* **AlN (5x5x3 mesh, sigma = 0.1 THz, 300 K; symmetric rates; selection-rule-forbidden events kept
  at zero; all 897 lifetimes and kappa_xx, kappa_zz fixed)** — certified witnesses from the
  counterexample review (256-bit ball arithmetic, exact energy conservation):
  basal: tau_x >= 7,038.17 ps (reference 34.25 ps, tau_CS = 19.190 ps; factor >= 205 over the
  reference); c axis: tau_z >= 4,322,183.47 ps (reference 27.45 ps, tau_CS = 15.517 ps; factor
  >= 1.57e5). The minima found are 19.197 ps (x) and 15.517 ps (z) in the class without the
  selection-rule zeros, and <= 19.44 ps (x), <= 15.517 ps (z) for box-constrained rates that keep
  the zeros (F = 100, section 2, P7); inner ratios >= 362 (x) and >= 2.78e5 (z).
  These numbers are mesh properties: at 7x7x5 the reference memory is 52.6 ps (x) and 50.6 ps (z)
  (kappa 243.03 / 272.33 W/(m K)); no large-tau optimisation was run there.
* **Superseded witnesses.** The previously quoted AlN maxima (6,906.7 ps x; 959,689.9 ps z;
  3.58e7 ps x near-singular) activate events whose three-phonon vertex is numerically zero
  (3,275 events, 7.4 %, 335 orbits; rates <= 4.4e-15 of the maximum versus >= 3.1e-9 for allowed
  events). 95.6 % of these events have a nontrivial common little group (mirror, glide, or C6v on
  Gamma–A) against 17.2 % of allowed events: the selection-rule interpretation is EVIDENCED
  statistically, not proved event by event. The 6,906.7 ps witness draws up to 99.998 % of some
  mode lifetimes from such events. Its factors over the reference (201.6 x, 3.49e4 z) are close to
  the replacements, which do not use those events. The 3.58e7 ps witness is a mode of weight 1.1e-6
  at about 32 s and is not used.
* **Unsymmetric rates are not AlN statements.** Witnesses with unconstrained rates break basal
  isotropy (tau_x = 114,019 ps but tau_y = 514 ps); they are not quoted as AlN properties.
* **Cone bound.** tau <= M^2 K(0) (Ben-Tal–Teboulle) is sharp over the cone when the lifetimes are
  free (T2); with lifetimes fixed it can still be attained in degenerate cases (current parallel to
  an event vector, section 2); for physical full-support currents the lifetime-constrained question
  is open. Quantitatively vacuous: M^2 K(0)/tau_ref >= 7.8e6 in every geometry computed.

## 2. Statements and status
Setting as in 02-A: C(g) = sum g a a^T on H = e^perp, x = C^+ b, K = b^T x, N = |x|^2, tau = N/K;
P(r) = {g >= 0 : diag C(g) = r}; F(r, K0) = admissible g in P(r) with K = K0. Standing hypotheses of
T1–T2: g admissible, b in H \ {0}, the full event set spans H (kernel dimension 1, checked).

**T1 (cone bound; Ben-Tal & Teboulle, Linear Algebra Appl. 139, 165 (1990); Cauchy–Binet proof
in note 18) — PROVED.** x/K lies in the convex hull of the basic solutions y_I of the ACTIVE bases
(a_alpha^T y = 0 for alpha in I, |I| = d-1, b^T y = 1); |x| <= M K, tau <= M^2 K, M = max_I |y_I|,
rate independent. Extends to non-admissible g with b in range C(g) (audit section 2.8).

**T2 (sharpness without lifetimes) — PROVED (standard weight-concentration argument).**
sup{tau(g) : g admissible, K(g) = K0} = M^2 K0. Proof: let I* attain M; g^eps = 1 on I*, eps
elsewhere; in the Cauchy–Binet weights c_J = eps^{|J \ I*|} with |J \ I*| >= 1 for J != I*, so
x/K -> y_{I*}; rescale to K = K0 (degree -1 homogeneity). The supremum is approached, not attained
in general: the rescaled rates on I* grow without bound and the lifetimes of the modes touched by
I* tend to zero, so T2 says nothing about F(r, K0). Related: Stewart, Linear Algebra Appl. 112, 189
(1989); O'Leary, Linear Algebra Appl. 132, 115 (1990); Hanke & Neumann, Linear Algebra Appl. 190,
137 (1993). No novelty is claimed.

**Lifetime-constrained sharpness — counterexample and open question.** Exact example (audit,
`review/check_cone_bounds.py`, instance B; 7 modes, 10 events): b = a_beta (current parallel to an
event vector); for every g supported on the basis J* + {beta}, x/K = y_{J*} exactly. With
r = (915/16, 7, 9/2, 29/9, 17/16, 2891/640, 153/80) and K0 = 4 the maximum over F(r, K0) is
attained and equals M^2 K0 = 44367253/1074614. A physical current b = D(v eps) has full support;
for such currents it is OPEN whether sup_F tau < M^2 K0. For symmetric rates T2 does not transfer
unless I* is a union of orbits.

**P3 (sub-network invariant form of M) — PROVED.** ker S_I^T = span{eps, phi_I},
|y_I| = |P_H D phi_I| / |b^T D phi_I| (phi_I rational for integer stoichiometry). For 1D N = 9 the
maximum is certified (audit: integer-rank and time-reversal certificates over all 9,657,700
subsets): M = 16284.006723562.

**L3 (DC response bounded below by RTA) — PROVED; sharp; novelty not checked.** Let s_max be the
largest number of DISTINCT modes touched by an active event (3 for p -> a + b, 2 for p -> 2a, 4 for
four-phonon events) and **r = diag C** (single-mode relaxation rates including the self-coupling
of repeated slots). Then K(g) >= K_RTA(r)/s_max, K_RTA(r) = sum_mu b_mu^2/r_mu. Proof: b_mu^2 <=
r_mu D_mu by Cauchy–Schwarz, D_mu = sum_{alpha touches mu} g_alpha (a_alpha^T x)^2, and
sum_mu D_mu <= s_max K. Sharp: exact networks give K/K_RTA = 1/3 (11 modes) and 1/2 (chain of
repeated-daughter events). **The bound is about diag C, not about self-energy linewidths:**
phono3py's linewidth counts a repeated-daughter self-coupling with half the weight that diag C
does, and for the chain with energies 1, 2, 4, 8, 16 (events 2->1+1, ..., 16->8+8, s_max = 2) a
rational current gives K/K_RTA(Gamma) = 7383516975/22150549363 ~ 0.3333 < 1/2, while
K/K_RTA(diag C) ~ 0.5406 >= 1/2 (counterexample review, `cx_attack5_chain_exact.json`). For AlN the
difference is negligible: K/K_RTA(Gamma) = 1.2500 (x), 1.4203 (z) versus 1.2572, 1.4319 with diag C.

**P4 (lower bounds on the memory) — PROVED.** (a) tau >= K/|b|^2, equality iff b is an eigenvector
of C on H. (b) N >= 2 b^T v - |C v|^2 for every v; with v = t e_mu:
tau_min >= max_mu b_mu^2/(rho_mu r_mu)^2 / K0 (valid, numerically weak here).

**P5 (convex structure; SDP/Schur) — PROVED with the extended-value K** (02-A):
{g in P(r) : K(g) <= K0} = {g in P(r) : [[C(g), b], [b^T, K0]] >= 0}; this set contains
non-admissible g with b in range C(g), so the identity is false for K defined only on admissible g.
K_lo(r) = min over P(r) of the extended K is an SDP and equals the infimum over admissible g when
P(r) contains an admissible point. SDP and mirror descent agree to <= 3.3e-6 where the SDP
converged. No convexity of tau is claimed.

**P6 (geometric bound on the gap) — PROVED.** tau_max/tau_min <= M^2 |b|^2 (vacuous numerically).

**P7 (ceiling under a factor-F amplitude prior) — PROVED** (counterexample review, two lines;
checked here). Let g_ref be admissible and g_alpha >= g_ref,alpha / F for every event. Then
C(g) - C(g_ref)/F is PSD, so lambda_min(C(g)|H) >= lambda_min(C(g_ref)|H)/F, and since x is in H,
tau(g) = |x|^2/(x^T C x) <= 1/lambda_min(C(g)|H). Hence **tau <= F / lambda_min(C_ref|H)**, with or
without lifetime and K data, symmetric or not. AlN: lambda_min = 9.7643e-4 THz, so
tau <= 81.50 F ps at 5x5x3 (197.2 F ps at 7x7x5), and tau_max/tau_min <= 4.247 F (x), 5.253 F (z).
Conversely, a witness of memory tau must suppress some event by a factor >= tau/81.50 ps (86.4 for
the 7,038.17 ps witness, 53,034 for the 4,322,183 ps witness). Inner values under the box
g_ref/F <= g <= F g_ref (symmetric rates, all lifetimes, K_xx and K_zz fixed; local optima):

| F | tau_x found [min, max] (ps) | tau_z found [min, max] (ps) | ceiling (ps) |
|---|---|---|---|
| 1 | 34.25 | 27.45 | — |
| 2 | [27.96, 45.63] | [20.71, 44.16] | 163.0 |
| 10 | [20.73, 153.48] | [15.72, 162.60] | 815.0 |
| 100 | [19.44, 480.12] | [15.517, 1421.2] | 8150 |

**Smooth amplitude prior — NUMERICAL (inner evidence, not a bound).** Rates g = g_ref exp P with P
a polynomial in the three phonon energies (symmetric in the daughters): up to degree 10 tau is
pinned (degree 6: the reference is isolated; degree 10: tau_x = 34.2530, tau_z = 27.4465 ps);
degrees 12–20 with |P| up to 30 gave no witness above 1.80x the reference (max 43.19 ps x,
49.43 ps z).

**T7 (infrared criterion).** (i) RTA moments — PROVED: M_m = int b^2 r^{-m} finite iff
m alpha < d (log at equality), under r = a(Omega) q^alpha with a bounded above and below and
int b0(Omega)^2 a^{-m} dOmega finite and positive; RTA mesh laws: converges (2 alpha < d), log N
(2 alpha = d), N^{2 alpha - d} (d/2 < alpha < d), N^d/log N (alpha = d). (ii) Event operator —
DERIVED UNDER H4 (|x_mu| >= c |b_mu|/r_mu at small q, c uniform in the mesh; a property of a given
operator, not of the class): N -> infinity for 2 alpha >= d, and tau >= |x|/|b| = sqrt(N)/|b|
gives tau -> infinity with no condition on K; mesh lower bound tau >~ N^{alpha - d/2}. Two-sided
mesh exponents for the event operator are UNJUSTIFIED. On the AlN 5x5x3 mesh x_mu r_mu/b_mu is
1.13–1.21 at the reference and 0.294 at the x-minimiser; no statement is made about the continuum
memory of AlN. (iii) M^2 K0 >= tau_ref diverges whenever tau_ref does (under H4).

**Dimension count — VERIFIED AT THE REFERENCE.** Full row rank of [W; grad K] at the AlN 5x5x3
reference: 903 (unconstrained) and 119 (orbits) rows, relative smallest singular values 9.0e-6 and
5.9e-6; local dimensions 43,503 and 2,508. With the 3,275 forbidden events fixed at zero the
symmetric class has 2,292 allowed orbits.

**Gap numbers — NUMERICALLY DEMONSTRATED (inner bounds).** Every quoted extremum is attained by an
explicit feasible witness; the AlN replacement witnesses are certified in 256-bit ball arithmetic
(tau_x: radius ~1e-66, K_xx residual 3.4e-14, lifetimes 5.8e-13, condition 4.2e4; tau_z: radius
< 1e-60, K_zz residual 2.9e-11, condition 1.1e7). Maxima are lower bounds on the supremum.

## 3. Additional data (counterexample review, inner evidence)
* kappa(T) at 100, 200, 300, 500, 800 K (with T-independent |Phi|^2 delta) fixed in addition:
  witnesses persist (tau_z >= 131,918 ps). The mode accumulation function below 5, 10, 15 THz
  fixed in addition: tau_z >= 61,254 ps. Neither closes the gap.
* K(z) at three Laplace rates (1/ns, 1/(100 ps), 1/(10 ps)) fixed in addition: the maxima found
  drop to 62.9 ps (x) and 33.9 ps (z). This is not a bound: a slow mode of small weight w at time T
  is invisible in finitely many K(z) yet adds w T to tau_mem.
* Lifetimes at five temperatures: inconclusive (constraint set numerically stiff; the rates are
  not identified, the memory was neither shown nor refuted to be).
* The large-tau witnesses differ strongly from the reference in K(z) at 1/(10 ns)–1/(100 ps)
  (22–53 % lower), in the accumulation function and in kappa(T) at other temperatures: they are
  observable in principle.

## 4. Fragile points and scope
1. The gap exists only without an amplitude prior (P7, smooth-prior evidence): a per-event factor
   F caps tau at 81.5 F ps (5x5x3); the 2–5 orders of magnitude need per-event suppressions of
   10^2–10^4 and more.
2. Global optimality of tau_max is not established; all maxima are lower bounds.
3. The selection-rule reading of the zero-vertex events is statistical, not proved event by event.
4. AlN numbers are 5x5x3 numbers. Reference comparison with phono3py's Omega': kappa_xx 214.30 vs
   214.28, tau_x 34.25 vs 34.21 ps. For z the earlier comparison (202.82 W/(m K), 27.99 ps) used the
   symmetrised Omega' without degeneracy averaging; the review's degeneracy-averaged rebuild gives
   202.18 W/(m K) and 27.54 ps, closer to the event operator (201.97, 27.45); the residual
   difference was not traced further. Solutions C^+ b differ from Omega' by 0.15 % (x), 0.79 % (z).
5. Debye meshes are pre-asymptotic for the exponent law.

## 5. Conclusion
Static transport data at one temperature (lifetimes and the DC tensor) determine the shortest
compatible Cattaneo time, tau_CS = K(0)/|b|^2, to within a factor 1.69 in every geometry tested,
but they leave the memory itself undetermined by orders of magnitude in the event-cone class.
What pins tau_mem is the microscopic amplitude structure: a factor-F prior per event confines it
to [tau_CS, 81.5 F ps] at 5x5x3, and smooth polynomial amplitudes kept it within 1.8x the
reference. The physically decisive open computation is therefore the fc3-parametrised rate class
(rates generated by symmetric anharmonic force constants, with the selection rules built in),
together with the lifetime-constrained sharpness question for full-support currents.

## 6. Reproduction
```
python -B bt_exact.py 1 9                         # E1
python -B scan.py tasks_A.json 6 ; python -B summarize_scan.py     # E3
python -B ir_scaling.py 3 1.0,1.5,2.0,phys 3,5,7,9,11 ../results/ir_scaling_d3.json   # E7
<phono3py env> -B aln_events.py --mesh 5 5 3 --sigma 0.1 --nsig 4 --out ../results/aln/events_m553_s0.1.npz
python -B aln_geometry.py ../results/aln/events_m553_s0.1.npz     # E8 validation
python -B aln_scan.py ../results/aln/events_m553_s0.1.npz m553_s0.1 4 3000   # E9
python -B certify_aln.py <witness.npz> ...        # arb certification
```
Audit and review checks: `review/` (not modified here); replacement witnesses
`review/cx_runs/nofor_*`, certification `review/cx_certify_new.json`, all runs
`review/cx_summary.json`. No file in the thesis repository was edited (one bytecode cache
regeneration, see 09-10).
