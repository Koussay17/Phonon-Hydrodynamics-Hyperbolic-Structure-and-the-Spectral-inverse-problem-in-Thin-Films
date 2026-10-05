# 02 — Definitions, assumptions and declared modelling choices

## A. Mathematical setting (all results are stated in this class)
* **Mode space.** R^n with the Euclidean inner product in entropy coordinates y = D^{-1} dn,
  D = diag(d_mu), d_mu = sqrt(N_mu (1 + N_mu)), N_mu the Bose occupation at temperature T > 0.
  Mode energies eps_mu > 0 (zero-frequency Gamma modes excluded).
* **Invariant.** e = D eps (energy). H = e^perp, d = n - 1. No other invariant is assumed: the
  event vectors of the full event set span H (checked numerically for every geometry used: kernel
  dimension 1).
* **Events.** Stoichiometric vectors s_alpha (parent -1, daughters +1; repeated daughter +2),
  event vectors a_alpha = D^{-1} s_alpha in H. s(alpha) = |supp a_alpha| is the number of DISTINCT
  modes touched by the event: 3 for p -> a + b (a != b), 2 for a repeated daughter p -> 2a; s_max is
  its largest value over active events (the constant of L3 in 13). Energy conservation
  s_alpha^T eps = 0 holds exactly: in integer arithmetic for the exact 1D lattice model; by the
  declared conserving stoichiometry (section C) for tolerance events.
* **Operators.** C(g) = sum_alpha g_alpha a_alpha a_alpha^T, g >= 0. g is *admissible* if the
  active event vectors span H, equivalently ker C(g) = span(e).
* **Current.** b = D (v * eps) (one column per Cartesian direction); b is orthogonal to e on a
  time-reversal-symmetric full grid (checked numerically, |b.e|/(|b||e|) < 1e-16).
* **Moments.** For admissible g: x = C^+ b (the unique solution in H of C x = b), K = b^T x,
  N = |x|^2, tau_mem = N/K = -K'(0)/K(0), where K(z) = b^T (z + C)^{-1} b on H.
  tau_mem is the relaxation time of the single (Cattaneo) pole K(0)/(1 + z tau) that matches the
  response and its first derivative at z = 0. Matches note 18 (tau_mem = |C^+ b|^2 / b^T C^+ b).
* **Extended-value K (convex side only).** For every g >= 0: K(g) = b^T C(g)^+ b if b is in
  range C(g), and +infinity otherwise (equivalently sup_y {2 b^T y - y^T C(g) y}). It coincides with
  K on admissible g, is convex and lower semicontinuous on R^m_{>=0}, and is the function for which
  {K <= K0} equals the Schur-complement LMI set (13, P5). Statements about tau and F use admissible g.
* **Data.** r = diag C(g) (mode relaxation rates = inverse lifetimes in this convention, including
  the self-coupling of repeated daughter slots; NOT the self-energy linewidth of RTA codes, which
  weights a repeated-daughter self-coupling by one half — L3 in 13 holds for diag C and can fail for
  self-energy linewidths), the invariant (energy only), and K0 = K(g) (and, for AlN, the DC tensor
  components K_xx, K_zz with symmetric rates, or all six components with unconstrained rates; the
  latter break basal isotropy and are not AlN statements).
* **Feasible set.** F(r, K0) = {g >= 0 admissible : diag C(g) = r, K(g) = K0};
  tau_min = inf_F tau_mem, tau_max = sup_F tau_mem. Symmetric version F_G: g constant on orbits of a
  group G of mode permutations preserving the geometry (time reversal; cubic point group for the
  Debye meshes; 6mm x time reversal, 24 operations, for AlN).
* **Inner/outer.** Values found by optimisation are attained by explicit feasible witnesses and
  give an INNER interval [tau_min_found, tau_max_found] contained in [tau_min, tau_max]; proved
  bounds give an OUTER interval. The gap ratio tau_max/tau_min is therefore bounded BELOW by the
  inner ratio.

## B. Debye-type event sets (declared model)
* Simple-cubic reciprocal lattice, N^d mesh, N odd (unique zone representatives, no
  self-conjugate zone-boundary points).
* Two branches with linear isotropic dispersion eps = c_b |kappa| (Euclidean norm of the zone
  representative), c = (1, 2) ("T", "L"); v = c_b kappa/|kappa|. Units: lattice constant 1,
  energies in units of c_T; kT = 0.5 * max(eps) unless stated.
* d = 1: exact integer resonance (energies c_b |kappa| in units of 2 pi/N). Umklapp
  T + T -> L requires k1 + k2 = 2N/3 mod N: exists for N divisible by 3.
* d >= 2: resonance tolerance |eps_p - eps_a - eps_b| <= 0.5 * c_T * (2 pi/N) with conserving
  stoichiometry (section C).
* Reference ("physical-model") rates: Klemens-type amplitude Gamma = eps_p eps_a eps_b / eps_max^3
  (1/2 for repeated daughters), resonance weight 1 (d = 1) or 1/(2 tol) (box), Bose factor
  sqrt(N_p(1+N_a)(1+N_b)(1+N_p)N_a N_b), normalisation 1/N^d. Emergent infrared exponent
  alpha ~ 2 (fitted, d = 1).
* Prescribed-exponent data: r = A eps^alpha (A fixed by the total rate of the reference); the
  reference rates are the KL projection of the Klemens rates onto {W g = r}.

## C. Conserving stoichiometry for tolerance events (declared)
For an accepted event with detuning Delta = eps_p - eps_a - eps_b != 0:
s = -l_p e_p + l_d (e_a + e_b), l_p = sqrt((eps_a + eps_b)/eps_p), l_d = 1/l_p,
so s^T eps = 0 exactly and l_p l_d = 1 (parent-daughter products unchanged). This is a
conserving, positive, rank-one approximation; it is NOT the exact Jacobian of the detuned Bose
bracket (note 19, section 5) and carries no error bound for inverse moments.

## D. AlN event geometry (declared)
* Inputs: Phonon Olympics AlN force constants (commit 0640f077; POSCAR, BORN, fc2 5x5x3,
  fc3 3x3x2), SHA-256 verified at every run; phono3py 4.5.0 production settings
  (make_r0_average = True, Rust kernels, NAC), Gaussian sigma = 0.1 THz, 300 K.
* Events: decay channel of phono3py's triplets q0 + q1 + q2 = G from the parent's row:
  p = (q0, j0) -> a = (-q1, j1) + b = (-q2, l), kept if |w0 - w1 - w2| <= 4 sigma and pp*g0 > 0;
  each unordered event once. The third phonon is computed as q2 = -q0 - q1 (mod G); phono3py's
  internal map returns the third phonon of the little-group representative triplet (see 09).
* Rates: g = conv*pp*g0 / (4 sinh(x_p/2) sinh(x_a/2) sinh(x_b/2)) (distinct daughters),
  conv*pp*g0 / (8 sinh(x_p/2) sinh(x_a/2)^2) with s = -e_p + 2 e_a (repeated daughter); this is
  the rank-one decomposition of phono3py's physical operator Omega' = D + C1 - (C0 + C2) J.
  Units THz; physical rates 4 pi g (1/ps); tau_mem[ps] = (N/K)/(4 pi).
* Reference rates: orbit average over the 24 operations (removes band-gauge noise in degenerate
  subspaces; changes kappa and tau_mem by < 1e-5 relative).
* Validation against phono3py's symmetrised Omega' at the same mesh and sigma (5x5x3):
  kappa_xx 214.30 vs 214.28 W/(m K), tau_mem,x 34.25 vs 34.21 ps; kappa_zz 201.97 vs 202.82,
  tau_mem,z 27.45 vs 27.99 ps; diagonal median relative deviation 8e-4; energy conserved to
  4.5e-16 (Omega' itself: 8.9e-4).
* Scope: a coarse-mesh, Gaussian-broadened, conserving event operator derived from AlN data.
  It is a realistic event GEOMETRY with realistic reference rates; it is not a converged AlN
  collision operator (mesh 5x5x3 lacks the infrared; kappa_LBTE at 31x31x17 is 302/286 W/(m K)).

* Selection rules: 3,275 of the 44,406 events (7.4 %, 335 orbits) have a numerically zero vertex
  (reference rates <= 4.35e-15 of the maximum versus >= 3.07e-9 for the others), with the
  statistical signature of space-group selection rules (95.6 % share a nontrivial mirror, glide or
  C6v little group; evidenced, not proved event by event). The 6mm x time-reversal tying does not
  impose these zeros; the AlN headline witnesses fix them at zero (2,292 allowed orbits).

## E. What the "unknown rates" class does and does not include
* **Declared class.** Arbitrary nonnegative rates on the fixed, energy-allowed event set (for AlN:
  symmetric under 6mm x time reversal, selection-rule zeros kept at zero), with NO amplitude model;
  data = the 300 K lifetimes (diag C) plus the DC tensor. Every gap statement is a statement about
  this class. "Realistic event geometry" refers to the event set only, not to the rates.
* Included: any nonnegative rates on the fixed event set (optionally symmetric). This is the
  model-agnostic class relevant to inferring a closure from static transport data without a
  microscopic amplitude model.
* Not included: smoothness of |Phi|^2 in the momenta, or parametrisation by a finite set of
  anharmonic force constants. With a first-principles fc3 the rates are determined (no
  identifiability problem); with fewer amplitude parameters than lifetimes the diagonal alone can
  pin the rates. The gap reported here is the gap of the event-cone class, not of fc3-parametrised
  models. Quantitatively (13, P7): a per-event prior g >= g_ref/F confines tau_mem to
  [tau_CS, F/lambda_min(C_ref|H)] = [tau_CS, 81.5 F ps] at 5x5x3, and smooth polynomial amplitudes
  (degree <= 20) gave no witness above 1.8x the reference.
