# Branch B-X — counterexample search, Workstream B (rank and injectivity of the 1D front-face response)

Date: 11 November 2026. Scripts and raw data: `derivations/B-counterexamples/` (all paths below relative to it).
Status: complete for this branch (numerical attack; no proofs).

## 0. Verdicts

All results refer to the physical impedance (transmission-line form of `src/forward_model.py`); see §1.1.
Units: ξ₁ = b_s = 1 at the point, θ written (ξ₁, b/b_s, τ_R/ξ₁², τ_ℓ/ξ₁²); see §1.2.

| Claim | Verdict | Exact counterexample | Practical failure (quantified) |
|---|---|---|---|
| B1 full-band rank 4 off Σ = {τ_R = τ_ℓ} ∪ {τ = 0} | SURVIVES search; Σ as stated is too large: the faces τ_ℓ = 0 and τ_R = 0 have rank 4; rank 3 on τ_R = τ_ℓ; rank 2 on τ_R = τ_ℓ, b = b_s | none | s₄ ∝ \|ln τ_R/τ_ℓ\| (exponent 1.0000); s₃, s₄ ∝ \|d\| along b = b_s; s₄ ∝ (τ/ξ₁²)^(−1/2) in the ballistic corner (§2.5) |
| B2 injectivity off τ_R = τ_ℓ | SURVIVES search | none (smallest distinct full-band misfit 1.2·10⁻²) | finite bands: GK ≡ Cattaneo with the same τ_R − τ_ℓ to 1.2·10⁻⁸ (thin band) and 2–7·10⁻³ (FDTR-like band); at b = b_s, ξ₁ × 0.21 at 2.3·10⁻³ (§3.3–3.4) |
| B3 N = 2 frequencies, rank 4 generically | SURVIVES as a generic local statement; BROKEN as an identification statement | θ = (1, 1, 0.5, 0.1), ωξ₁² = (0.2, 2): five other θ' reproduce both values of Z to 60 digits (§4.4); singular pairs lie on curves, e.g. (4.217, 6.188) for θ = (1, 2, 0.5, 0.1) | thick pairs: s₄/s₁ ∝ e^{−2 Re A} (2·10⁻⁶⁵ example, §4.3) |
| B4 interface resistance, generic rank 5 | SURVIVES search | none; R = 0 face rank 5; rank 4 on τ_R = τ_ℓ; rank 3 on τ_R = τ_ℓ, b = b_s, R = 0 | s₅ ∝ \|d\| near resonance; 10⁻⁵ in the ballistic / high-R corner (§5) |

## 1. Model, units, structure

### 1.1 The impedance formula of the brief

The brief writes Z = (1/B)(B tanh A + b_s√p)/(B + b_s√p tanh A). This is the reciprocal ratio of the
transmission-line form; it tends to b_s√p/B² as e → 0 instead of 1/(b_s√p) (`t00_verify_model.py`:
|Z b_s√p − 1| = 0.75 at e → 0 for b = 2b_s). It equals the physical formula with y → 1/y, i.e. with the
reflection coefficient ρ → −ρ (`s00_identities.py`, item 1, exact). The brief has since been corrected to
the physical form.

**All results of this report use the physical form**

Z = (Z_s + tanh A / B)/(1 + B Z_s tanh A) = (1/B)(B + b_s√p tanh A)/(b_s√p + B tanh A),  Z_s = R + 1/(b_s√p),

identical to `src/forward_model.response` (200 random stacks, 60 frequencies each, 10³–10¹¹ s⁻¹:
max |Z/Z_ref − 1| = 8.6·10⁻¹⁴) and to note 16. The swapped form was only used for one side scan
(`b1_scan.py brief`), kept as a side note in §7.

Stable form used everywhere: E = e^{−2A}, y = B Z_s, ρ = (y − 1)/(y + 1), Z = (1/B)(1 + ρE)/(1 − ρE).

### 1.2 Nondimensionalisation

b_s is fixed and known; it is used as the effusivity unit (β = b/b_s, Z b_s). This is a choice of units,
not a parameter reduction: β remains free. The time scaling
(ξ₁, τ_R, τ_ℓ, ω) → (λ^{1/2}ξ₁, λτ_R, λτ_ℓ, ω/λ), Z → λ^{1/2} Z maps the data set "all ω > 0" onto
itself, so the full-band rank and injectivity at θ are those at the rescaled point; we put ξ₁ = 1 at the
point examined (t_R = τ_R/ξ₁², t_ℓ = τ_ℓ/ξ₁²). The Jacobian is always taken in all four parameters,
ξ₁ included. For a finite frequency set the frequencies rescale too (ω ξ₁² is what matters).
For B4, R is measured in units of ξ₁/b_s.

### 1.3 Jacobian structure (exact; `s00_identities.py`)

With u = A ∂Z/∂A, v = B ∂Z/∂B (A fixed), w = u − v, g_X = p/(2(1 + τ_X p)):

∂Z/∂ln ξ₁ = u, ∂Z/∂ln b = v, ∂Z/∂τ_R = g_R w, ∂Z/∂τ_ℓ = −g_ℓ w;  AB = ξ₁ b p, A/B = (ξ₁/b)(1+τ_R p)/(1+τ_ℓ p).

Consequences:
* τ_R = τ_ℓ: the τ-columns are opposite, (0,0,1,1) is in the kernel: rank ≤ 3.
* τ_R = τ_ℓ and b = b_s: y ≡ 1, ρ ≡ 0, Z ≡ 1/B for every ξ₁, u ≡ 0 and v = −w: rank ≤ 2.
* Off resonance: rank < 4 ⟺ q := (u+v)/(2(u−v)) ∈ span_R{1, g_R, g_ℓ}, which implies that
  Q := −v/u = (1−2q)/(1+2q) is a real rational function of p (degree ≤ 2), where (exact identity)
  Q = sinh A [cosh A (1 + y²) + 2y sinh A] / (A (1 − y²)).
  Obstruction (sketch, not a proof): along p = iω → ∞, Q grows like e^{2 Re A} when τ_R, τ_ℓ > 0;
  for Cattaneo (τ_ℓ = 0) Re A → ξ₁/(2√τ_R) and Q ≈ sinh 2A/(2A) has infinitely many sign changes of
  Re Q; for τ_R = 0, A → ∞ at p = −1/τ_ℓ (essential singularity of sinh A). None is compatible with a
  rational Q (the same argument gives rank exactly 3 on the resonance when b ≠ b_s). The numerical attacks
  below found no exception.

Verification of the implementation: analytic Jacobian against 50-digit numerical differentiation
(both variants, with R): max relative error 1.2·10⁻³⁰; double against 40-digit: 2.5·10⁻¹².

Metrics. Data = relative response (Re, Im of dZ/Z = d ln|Z| + i d arg Z), weights d ln ω.
"raw" singular values (log parameters; ξ₁² ∂/∂τ for a vanishing time; R₀ = ξ₁/b for R = 0) and
"equilibrated" ones (columns scaled to unit norm, s₄ ∈ [0, 1], inverse collinearity index).
Equilibration must not be used to decide an exact rank when a column vanishes identically (it turns
the zero column of the invisible-film set into O(1) noise); exact ranks below use raw values.
Full band: characteristic frequencies {1/ξ₁², min(β², β⁻²)/ξ₁², 1/τ_R, 1/τ_ℓ, τ_ℓ/(τ_R ξ₁²)} ± 1 decade,
on a grid that resolves the thermal-wave resonances of Cattaneo-like films
(steps ≤ 0.1 in ln ω and ≤ min(0.2, (1−|ρE|)/3) in Im 2A; `zmodel.resolved_grid`).


## 2. B1 — full-band rank

### 2.1 Random scans (`b1_scan.py`)

Log-uniform θ (ξ₁ = b_s = 1), one tenth each on the strata τ_ℓ = 0, τ_R = 0, β = 1.
* Pass 1 (log grid, 24 points per decade, ± 4 decades; β ∈ [10⁻³, 10³], t ∈ [10⁻⁵, 10⁵], 20 000 points):
  smallest equilibrated s₄ off resonance (|ln τ_R/τ_ℓ| > 0.5): 4.1·10⁻⁷ at (β, t_R, t_ℓ) = (872, 7.1·10³, 9.8·10⁴).
  This grid aliases the thermal-wave resonances and over-weights the asymptotic high-frequency region
  (see 2.5 e); its Cattaneo-stratum values (median 7·10⁻⁷) are an artefact of that, not a property of the model.
* Pass 2 (resolved grid, ± 1 decade; β ∈ [10⁻³, 10³], t ∈ [10⁻⁴, 10⁴]; 19 848 points, 152 dropped because
  > 6·10⁴ nodes were needed): equilibrated s₄ medians 0.041 (interior), 0.136 (τ_ℓ = 0), 0.077 (τ_R = 0),
  0.025 (β = 1). Off resonance (|ln τ_R/τ_ℓ| > 0.1): minimum 3.5·10⁻⁶ at (720, 8.6·10⁻³, 9.8·10³)
  (ballistic, high-effusivity corner); restricted to β ∈ [10⁻², 10²], t ∈ [10⁻², 10²]: minimum 3.0·10⁻⁴,
  1 % quantile 9.4·10⁻⁴, median 0.065 (2 222 points).

### 2.2 Exact versus practical (`b1_mpcheck.py`)

The smallest values were recomputed in 50- and 80-digit arithmetic on the same grids:
(886, 9.7·10⁴, 0): s₄ = 1.368907·10⁻¹¹ at 50 and at 80 digits (double: 1.368907·10⁻¹¹);
(872, 7145, 98173): s₃ = 1.7946·10⁻⁶, s₄ = 7.39356·10⁻⁷ at both precisions; pass-2 minima (unit weights on
the resolved nodes, ≤ 2 500 nodes): (720, 8.6·10⁻³, 9.8·10³) 5.223098·10⁻⁶, (15.4, 8.8·10³, 0)
5.773991·10⁻⁶, (421, 0, 9.8·10³) 7.546568·10⁻⁶, identical in double and at 50 digits. An exact
deficiency would give s₄ ~ 10⁻⁵⁰ resp. 10⁻⁸⁰. Every small value found is a property of the parameter point
(practical), not a rank drop. (A sampled Jacobian of rank 4 certifies full-band rank 4 whatever the grid;
grid resolution only matters for the size of s₄.)

### 2.3 Adversarial minimisation (`b1_adversarial.py`)

Nelder–Mead on log₁₀ s₄ (both metrics) over log₁₀β ∈ [−2, 2], log₁₀ t ∈ [−3, 3] with the hard constraint
|ln τ_R/τ_ℓ| ≥ c (c = 0.1, 1), 60 starts each, and on the two boundary strata (30 starts each).
Result: every local minimum with s₄ < 10⁻² lies on a face of the box (asymptotic corners β → 10², t → 10³
or t_R → 10⁻³) or on the resonance constraint; the few interior local minima are well conditioned:
s₄ = 0.16 (τ_ℓ = 0, β = 1.10, t_R = 0.76), 0.01 (β = 0.023, t_R = 0.11, t_ℓ = 0.19), 0.01 (τ_R = 0,
β = 1.06, t_ℓ = 0.13); 50-digit values coincide with double. Lowest values reached: 10⁻⁴·⁵⁴ (c = 0.1,
on the constraint at the box corner), 10⁻⁴·¹⁹ (c = 1, box corner β = 100, t_R = 10⁻³, t_ℓ = 993).
No hidden interior degeneracy was found.
(An earlier version of this search clipped both times to the same box face, which lands exactly on
τ_R = τ_ℓ and produced spurious minima of 10⁻¹⁸; replaced by penalised, unclipped variables.)

### 2.4 Exact rank on the strata (`b1_strata_mp.py`, `b1_resonance_scaling.py` (c)); 50 digits, raw values

| stratum (30 random points each, β ∈ [10^-1.5, 10^1.5], t ∈ [10⁻², 10²]) | s₃/s₁ | s₄/s₁ | rank |
|---|---|---|---|
| Cattaneo τ_ℓ = 0 < τ_R (column ξ₁² ∂/∂τ_ℓ) | 1.2·10⁻³ – 3.9·10⁻² | 1.0·10⁻⁴ – 1.7·10⁻² | 4 |
| τ_R = 0 < τ_ℓ | 3.1·10⁻³ – 0.30 | 1.1·10⁻³ – 8.7·10⁻² | 4 |
| b = b_s, τ_R ≠ τ_ℓ | 8.7·10⁻³ – 0.23 | 1.3·10⁻³ – 0.13 | 4 |
| high-frequency matching β = √(τ_R/τ_ℓ) (ρ → 0 as ω → ∞) | 9.0·10⁻³ – 0.25 | 1.7·10⁻³ – 9.7·10⁻² | 4 |
| resonance τ_R = τ_ℓ, β ≠ 1 | 1.1·10⁻² – 0.29 | ≤ 4·10⁻⁵² | 3 |
| invisible film τ_R = τ_ℓ, b = b_s | ≤ 4.4·10⁻⁵² | ≤ 2.2·10⁻⁵² | 2 |

Correction to the claim: the boundary strata τ_ℓ = 0 and τ_R = 0 (with τ_R ≠ τ_ℓ) are **not** degenerate.
The map is analytic across τ = 0 on the imaginary axis (1 + τp ≠ 0 for p = iω), the one-sided Jacobian in
the linear time coordinate has rank 4. In log-time coordinates the boundary sits at −∞ and the column
τ ∂Z/∂τ vanishes linearly in τ; a rank drop "at τ = 0" seen in log coordinates is a coordinate artefact.
The exact degenerate set is Σ = {τ_R = τ_ℓ} (rank 3) ⊃ {τ_R = τ_ℓ, b = b_s} (rank 2, ∂Z/∂ξ₁ ≡ 0).
The matched case b = b_s is not degenerate off resonance.

### 2.5 Scaling laws near Σ and in the corners (`b1_resonance_scaling.py`, `b1_scan.py`, rays)

(a) Fourier resonance, d = ln(τ_R/τ_ℓ) → 0: raw s₄ = C|d|, fitted exponent 1.0000 in all 12 cases
(β ∈ {0.1, 0.5, 2, 10}, τ̄/ξ₁² ∈ {10⁻², 1, 10²}, d from 10⁻¹ to 10⁻⁹). C ≈ 0.21 at τ̄ = 10⁻²ξ₁² for every
β, 0.09–0.21 at τ̄ = ξ₁², and smaller in the ballistic high-effusivity corner (C = 1.9·10⁻² at β = 2,
τ̄ = 100ξ₁²; 4.1·10⁻³ at β = 10, τ̄ = 100ξ₁²).
The near-null direction is the common shift of the two times.
(b) Invisible-film set {d = 0, β = 1}: on the resonance, s₃ = C'|β − 1| (exponent 1.000; C' = 0.83, 0.13,
0.83 for τ̄/ξ₁² = 10⁻², 1, 10²); along b = b_s, s₃ and s₄ are both linear in |d| (exponents 1.000;
s₄/|d| = 6.3·10⁻³, 3.3·10⁻², 3.7·10⁻² for τ̄/ξ₁² = 10⁻², 1, 10², against C ≈ 0.21 off b = b_s at
τ̄ = 10⁻²ξ₁²).
Two singular values vanish at the codimension-2 set, both linearly.
(c) Ballistic corner τ/ξ₁² → ∞ (Knudsen number ~ √(τ/ξ₁²) → ∞): s₄ ∝ (τ/ξ₁²)^(−1/2) (rays β = 0.1, 1, 10 with
τ_R/τ_ℓ = 10 and β = 10 with τ_R/τ_ℓ = 1/10; the factor 10^{1/2} per decade is reached for t ≳ 10² at β ≥ 1
and approached at β = 0.1); for τ_ℓ > τ_R also s₃ ∝ (τ/ξ₁²)^(−1/2) (practical rank 2).
Mechanism: the non-Fourier information sits at ω ~ 1/τ ≪ 1/ξ₁², where the film is thermally thin and
enters only through O(ξ₁√ω) corrections.
(d) β → ∞ or β → 0 at fixed times: s₄ saturates (≈ 0.04 equilibrated at τ_R = 10ξ₁², τ_ℓ = ξ₁²); the limits
are the adiabatic-back (coth A/B) and isothermal-back (tanh A/B) films, which are themselves of rank 4.
(e) Window dependence of any "full-band conditioning" number: for Cattaneo films the reflections never
decay (Re A → ξ₁/(2√τ_R)), and above 1/τ_R the sensitivities to ln ξ₁ and ln τ_R are dominated by the
transit time ξ₁√τ_R = e/v, growing like ω; with log weights the equilibrated s₄ falls as the window is
extended (β = 2, τ_R = 0.5ξ₁²: 0.066, 0.006, 8·10⁻⁵ for 1, 2, 4 decades of padding). For GK films the
high-frequency region sees only b√(τ_ℓ/τ_R) (b, τ_R, τ_ℓ columns collinear there). The full-band rank is
unaffected; a full-band condition number is not a well-defined quantity without a declared band.
(f) In the full-band problem ξ₁ → 0 at fixed times is the ballistic corner (c) and ξ₁ → ∞ is t → 0: for
β ≠ 1 the raw s₄ stays ≈ 0.4 down to t = 10⁻⁴ (rays β = 0.1, 10, τ_R/τ_ℓ = 10); at b = b_s it falls like t
(5.9·10⁻⁴ at t_ℓ = 10⁻⁴), the approach to the corner {b = b_s, τ_R = τ_ℓ = 0} of the invisible-film set.
The band-limited versions of these limits are in §6.

### 2.6 Verdict B1

SURVIVES the search: no θ with τ_R ≠ τ_ℓ and full-band rank < 4 was found (20 000 + 19 848 scanned,
adversarial minimisation, 180 strata points at 50 digits), and the analytic criterion of §1.3 has an
obstruction in every regime. The claimed exceptional set must be corrected (boundary strata are rank 4;
rank 3 on the resonance, rank 2 on resonance ∩ {b = b_s}). PRACTICALLY ILL-CONDITIONED near Σ (linear
laws (a), (b)) and in the ballistic corner (law (c)).

## 3. B2 — injectivity

### 3.1 Search design (`b2_injectivity.py`)

Eight targets θ (ξ₁ = b_s = 1): T1 GK (β, t_R, t_ℓ) = (2, 0.5, 0.1); T2 GK with τ_ℓ > τ_R (0.5, 0.1, 0.5);
T3 Cattaneo (2, 0.5, 0); T4 matched effusivity (1, 0.5, 0.1); T5 near resonance (2, 0.5e^{0.05}, 0.5e^{−0.05});
T6 ballistic (10, 100, 10); T7 low effusivity (0.1, 10, 1); T8 τ_R = 0 (2, 0, 0.3).
Four bands: full (characteristic frequencies ± 1 decade, resonances resolved) and three windows of
4.3 decades (10 kHz–200 MHz) with upper edge ω_b ξ₁² = 10⁻² (thin), 2 (mid; this is the position of the
500 nm AlN/sapphire scenario, ξ₁² ≈ 3 ns), 10² (thick). For θ' in each stratum (GK, Cattaneo, τ_R = 0,
Fourier): 40 (GK) or 20 multistart trust-region least-squares fits of the RMS relative complex misfit,
starts log-uniform in ξ₁' ∈ [0.1, 10], β' ∈ [0.01, 100], t' ∈ [10⁻⁴, 10⁴]. Plus the explicit swap
τ_R ↔ τ_ℓ with (ξ₁, β) refitted ("swap2p"), then all four refitted from there ("swap4p").
"Far" = separation > 0.05 in some log-parameter (or a different stratum).

### 3.2 Full band (the function-level question)

Smallest misfit of a far minimum: T5 (near resonance, d = 0.1) 1.18·10⁻² at θ' = (1.22, 1.81, 7.04, 8.08);
T4 (b = b_s) 2.05·10⁻² at θ' = (0.112, 0.988, 3.21, 0.72); all others ≥ 8.7·10⁻² (T6) and typically
0.2–0.4. The swap τ_R ↔ τ_ℓ with (ξ₁, b) refitted leaves misfits 2.7·10⁻² (T5) to 0.72 (T7) (T1 0.52,
T2 0.33, T4 0.29, T6 0.20); refitting all four from the swapped point returns to θ itself (T1, T5) or stalls
at ≥ 8.7·10⁻² (T2, T4, T6, T7). No pair with a misfit approaching round-off was found. (Refinement of the
best pairs: §3.4.)

### 3.3 Finite bands: explicit practical non-injectivity

| target | band | θ' (ξ₁', β', t_R', t_ℓ') | RMS misfit | mechanism |
|---|---|---|---|---|
| T1 GK (1, 2, 0.5, 0.1) | thin | (1.000, 2.000, 0.3996, 0) Cattaneo | 1.2·10⁻⁸ | only τ_R − τ_ℓ at O(ωτ) |
| T2 (1, 0.5, 0.1, 0.5) | thin | (1.000, 0.500, 0, 0.3996) τ_R = 0 | 4.4·10⁻⁸ | idem, τ_ℓ − τ_R |
| T4 (1, 1, 0.5, 0.1) | thin | (1.014, 1.000, 0.395, 0) Cattaneo | 9.8·10⁻⁹ | idem |
| T5 near resonance | thin | (1.000, 2.000, 0.0498, 0) Cattaneo | 7.1·10⁻⁹ | idem |
| T1 GK | mid | (1.073, 1.937, 0.345, 0) Cattaneo | 6.7·10⁻³ | idem, compensated by ξ₁, b |
| T2 | mid | (1.128, 0.530, 0, 0.316) τ_R = 0 | 6.0·10⁻³ | idem |
| T4 b = b_s | mid | (0.209, 0.973, 1.771, 0.536) GK | 2.3·10⁻³ | at b = b_s only ξ₁(τ_R − τ_ℓ) at leading order |
| T5 near resonance | mid | (1.032, 1.972, 0.0242, 0) Cattaneo | 2.0·10⁻³ | near-Fourier |

For orientation: a relative complex misfit of 2·10⁻³ at every frequency corresponds, e.g., to 0.08° of phase
together with 0.14 % of amplitude.
Matched effusivity is the worst case for the thickness: constraining ξ₁' = e^D ξ₁ and refitting the other
three (`practical_pairs.py` (a); 12 starts, upper bounds of the profile) gives, in the mid band,
misfits 2.3·10⁻³, 2.6·10⁻³, 1.6·10⁻³ for D = −1.5, −1, −0.2 at b = b_s (a factor 4.5 in ξ₁ costs 0.23 %;
thicker films are harder to mimic: 1.1·10⁻², 1.8·10⁻² for D = +0.5, +1), against 3.3·10⁻², 3.0·10⁻²,
1.1·10⁻² at b = 2b_s; for a nearly Fourier matched film (Cattaneo, τ_R = 0.05ξ₁²) the profile stays below
8.8·10⁻⁴ for D ∈ [−1.5, 0.2] and reaches 2.2·10⁻³ at D = +1. On the full band the b = b_s profile rises to
2–8·10⁻².
Low-frequency gauge τ_X → τ_X + c (both times shifted, ξ₁ and b unchanged; β = 2, τ_R/τ_ℓ = 5, mid band):
misfit ≈ 1.3·10⁻² (c/τ_ℓ)(ω_b τ_R)² (ω_b τ_R = 1, 0.1, 0.01: 1.32·10⁻³, 1.26·10⁻⁵, 1.22·10⁻⁷ for
c = 0.1 τ_ℓ; linear in c).

### 3.4 Refinement of the best pairs (`b2_refine_pairs.py`)

Each pair was re-minimised (Levenberg–Marquardt, tolerances 10⁻¹⁵), its misfit recomputed in 40-digit
arithmetic, and evaluated on bands widened by 1 and 2 decades on each side:

| target, band, θ' | band misfit (refined / 40 digits) | +1 decade | +2 decades | full band |
|---|---|---|---|---|
| T1 thin, Cattaneo (1, 2, 0.3996, 0) | 1.194·10⁻⁸ / 1.194·10⁻⁸ | 6.7·10⁻⁶ | 1.9·10⁻³ | 0.65 |
| T1 mid, Cattaneo (1.073, 1.937, 0.345, 0) | 6.70·10⁻³ / 6.70·10⁻³ | 0.13 | 0.61 | 0.59 |
| T2 mid, τ_R = 0 (1.128, 0.530, 0, 0.316) | 6.04·10⁻³ / 6.04·10⁻³ | 0.11 | 0.26 | 0.32 |
| T4 mid, GK (0.209, 0.973, 1.771, 0.536) | 2.26·10⁻³ / 2.26·10⁻³ | 3.7·10⁻² | 6.4·10⁻² | 8.4·10⁻² |
| T5 mid, Cattaneo (1.032, 1.972, 0.0242, 0) | 2.00·10⁻³ / 2.00·10⁻³ | 4.1·10⁻² | 0.28 | 6.3·10⁻² |
| T8 thin, GK (0.973, 2.034, 297.6, 296.6) | 4.06·10⁻⁶ / 4.06·10⁻⁶ | 4.7·10⁻⁴ | 1.5·10⁻² | 0.75 |
| T4 full, GK (0.112, 0.988, 3.214, 0.720) | 2.05·10⁻² / 2.05·10⁻² | 2.5·10⁻² | 2.7·10⁻² | 2.05·10⁻² |
| T5 full, GK (1.221, 1.806, 7.037, 8.076) | 1.18·10⁻² / 1.18·10⁻² | 1.19·10⁻² | 1.19·10⁻² | 1.18·10⁻² |
| T1 full, GK (25.56, 1.697, 0.279, 0.074) | 0.204 / 0.204 | 0.227 | 0.229 | 0.203 |

Every misfit stalls at a positive value (no further decrease, identical at 40 digits) and grows when the
band is widened: these are practical, not exact, counterexamples. The full-band minima are genuine
distinct local minima of the misfit with positive value.

### 3.5 Verdict B2

SURVIVES the search as a function-level statement: no distinct θ, θ' with equal full-band response
(smallest full-band misfit of a distinct local minimum 1.2·10⁻², near the resonance). An asymptotic
argument supports it (sketch): Z on iω determines Z on (0, ∞) by analytic continuation; as p → +∞,
E = e^{−2A} → 0 (all strata with τ_R > 0), Z → 1/B with 1/B² = (1+τ_R p)/(b² p (1+τ_ℓ p)) rational, which
fixes (b, τ_R, τ_ℓ) off resonance, and the e^{−2A} correction fixes ξ₁ unless ρ ≡ 0 (the invisible film);
for τ_R = 0, Z ~ coth(ξ₁/√τ_ℓ)/(b√τ_ℓ p) separates that stratum. PRACTICALLY NON-INJECTIVE in finite
bands: in any band with ωτ ≪ 1 a GK film is reproduced by the Cattaneo (or τ_R = 0) film with the same
τ_R − τ_ℓ up to O((ωτ)²), and at b = b_s the thickness trades against τ_R − τ_ℓ.

## 4. B3 — two frequencies

M₂(θ; ω₁, ω₂) = [Re J(ω₁); Im J(ω₁); Re J(ω₂); Im J(ω₂)], J = ∂ln Z/∂(ln ξ₁, ln b, ln τ_R, ln τ_ℓ).

### 4.1 When is M₂ singular for every pair? (exact criterion)

For fixed θ, det M₂ ≡ 0 in (ω₁, ω₂) iff the 2-planes K(ω) = ker_R J(ω) ⊂ R⁴ pairwise intersect, i.e.
either (case 1) they share a line — the full-band rank is < 4 — or (case 2) they lie in a common
hyperplane: a fixed real covector n lies in span_R{Re J(ω), Im J(ω)} for every ω. Writing the columns as
[m, w, g_R w, g_ℓ w] (m = (u+v)/2) and solving Re(ζ) = n₂, Re(ζ g_R) = n₃, Re(ζ g_ℓ) = n₄ for ζ = γw,
consistency off resonance forces n₃ = n₄ = c, n₂ = 2c(τ_R + τ_ℓ), ζ = 2c(1+τ_R p)(1+τ_ℓ p)/p
(`s00_identities.py` item 5), and case 2 becomes

  Re Φ(iω) = const for all ω > 0,  Φ(p) = (1 + τ_R p)(1 + τ_ℓ p) q(p)/p,  q = (u+v)/(2(u−v)).

In the thermally thick regime q → −1/2 exponentially, so Re Φ → −(τ_R + τ_ℓ)/2: case 2 holds
asymptotically for pairs that are both thick (explains 4.3). At low frequency Re Φ is large and positive
(e.g. 63 at ω = 10⁻³ for (β, t_R, t_ℓ) = (2, 0.5, 0.1)), so case 2 fails there.

### 4.2 Tests (`b3_two_freq.py`)

* 1 500 random (θ, ω₁, ω₂) (β ∈ [10⁻², 10²], t ∈ [10⁻³, 10³], ω ξ₁² ∈ [10⁻³, 10³]) at 50 digits: Hadamard
  ratio |det M₂|/∏‖col‖ median 5.5·10⁻⁴, 10 % quantile 3.3·10⁻⁸; s₄/s₁ median 3.0·10⁻⁴.
* Ball arithmetic (Arb, 200 bits): det M₂ ≠ 0 certified at 200/200 random points
  (β ∈ [10^-1.5, 10^1.5], t ∈ [10⁻², 10²], ω ξ₁² ∈ [10⁻², 10²]); each certificate also certifies full-band
  rank 4 at that θ (B1) and, by continuity, rank 4 on a neighbourhood of the pair.
* Case-2 defect over 300 random θ (spread of Re Φ on 400 frequencies, normalised): minimum 2.3·10⁻⁴ at
  (β, t_R, t_ℓ) = (1.24, 473, 27.3) (ballistic corner), median 0.21; never zero.
* Adversarial: minimising over θ (|ln τ_R/τ_ℓ| ≥ 0.1 or 1) the best Hadamard ratio over all pairs of a
  25-point grid drives θ to the box corner (β ≈ 10², t ≈ 10³): 10⁻⁵·⁸; no interior zero.
* Non-generic pairs exist (`b3_singular_pairs.py`, 30 digits, sign map of det M₂ on a 33 × 33 grid,
  ω ξ₁² ∈ [10⁻², 10²]): off-diagonal sign changes for 4 of 5 θ tested, i.e. singular curves in the (ω₁, ω₂)
  plane. Explicit singular pair: θ = (1, 2, 0.5, 0.1), (ω₁, ω₂)ξ₁² = (4.2170, 6.1881): s₄/s₁ = 3.5·10⁻³²
  (zero at 30 digits). For the Cattaneo film (1, 2, 0.5, 0) no off-diagonal sign change was found on that
  grid. Besides these curves, the diagonal ω₁ = ω₂ is trivially singular. This is consistent with
  "generic" and shows the qualifier is necessary.

### 4.3 Practical failure

If both frequencies lie where the film is thermally thick, ∂Z/∂ξ₁ ∝ e^{−2A} and
s₄/s₁ ∝ e^{−2 Re A(min ω)}. Explicit example: θ = (1, 29.10, 37.24, 0.1792), ω₁ξ₁² = 60.48, ω₂ξ₁² = 757.2:
s₄/s₁ = 2.1432·10⁻⁶⁵ at 100 and 200 digits (2.8·10⁻⁶⁴ at 50 digits is round-off), with
e^{−2 Re A(ω₁)} = 3.2·10⁻⁶⁶. In the random sample at least 1 % of the pairs are of this type (Hadamard
ratio below 10⁻³⁰⁰). Band laws for many frequencies are in §6.

### 4.4 Global aspect: two frequencies do not identify θ (`b3_global_N2.py`) — EXACT COUNTEREXAMPLE

With N = 2 the data map R⁴ → R⁴ is square; rank 4 gives local injectivity only. Solving
Z(θ'; ω_k) = Z(θ; ω_k), k = 1, 2, from 300 log-uniform starts and refining every root with a 40-digit Newton
iteration (mp.findroot):

θ = (ξ₁, β, t_R, t_ℓ) = (1, 1, 0.5, 0.1) (b = b_s), (ω₁, ω₂) ξ₁² = (0.2, 2.0): six distinct exact solutions,
polished by Newton in 60-digit arithmetic (residuals 10⁻⁵⁵ to 10⁻⁶²; 40-digit values in
`data/b3_global_N2_roots40.json`):

| ξ₁'/ξ₁ | b'/b_s | τ_R'/ξ₁² | τ_ℓ'/ξ₁² | full-band RMS misfit to θ |
|---|---|---|---|---|
| 1 | 1 | 0.5 | 0.1 | 0 (the target) |
| 0.0707235524 | 0.944355050 | 4.66562510 | 0.748602149 | 0.11 |
| 0.712987928 | 11.3841859 | 487.225744 | 4.57596138·10⁻³ | 2.8 |
| 1.00596395 | 16.2128067 | 977.124091 | 2.30888665·10⁻³ | 2.8 |
| 2.35462321 | 38.2506010 | 5386.77066 | 4.23041312·10⁻⁴ | 2.8 |
| 2.74928911 | 44.6832602 | 7346.63697 | 3.10374191·10⁻⁴ | 2.9 |

All lie in the open parameter domain, off the resonance, and differ from θ elsewhere in frequency
(e.g. 1.3·10⁻² and 2.3·10⁻² relative at ωξ₁² = 0.6 and 1 for the first one), consistent with B2. For four other (θ, pair) cases only the target was
found from 300 starts (multiplicity depends on θ and on the pair). With N = 3 (two cases) only the target
was found.

### 4.5 Verdict B3

As stated (generic rank 4 with N = 2): SURVIVES (rank 4 for generic θ and generic pairs; 200/200
certified; the exceptional pairs form curves, explicit example in 4.2).
Read as an identifiability statement ("two frequencies suffice"): BROKEN — explicit exact counterexample
in 4.4 (five parameter sets distinct from θ reproducing Z(ω₁), Z(ω₂) to 60-digit working precision). Local rank does not give
global uniqueness for the square system; N = 3 showed no second solution in the cases tried.
Not uniform even locally: exponentially ill-conditioned for thick pairs (4.3), algebraically for pairs in
the thin or low-ωτ regimes (laws of §6).

## 5. B4 — interface resistance R (`b4_resistance.py`)

Z_s = R + 1/(b_s√p), θ₅ = (ξ₁, b, τ_R, τ_ℓ, R). New column ∂Z/∂R = 4E/((1 − ρE)²(y + 1)²), y = B Z_s;
exact relation v = −Z + Z_s ∂Z/∂R (not a linear dependence: Z_s depends on p). Window extended to the
resistance frequency 1/(R b_s)².

* Random scan, 3 000 points (β ∈ [10⁻², 10²], t ∈ [10⁻³, 10³], R ∈ [10⁻³, 10³] ξ₁/b_s, one sixth Cattaneo);
  143 Cattaneo points with unresolvable resonance combs (> 6·10⁴ nodes) dropped. Equilibrated s₅:
  median 0.021, minimum 1.4·10⁻⁵ at (β, t_R, t_ℓ, R) = (0.0167, 686, 823, 0.024) (near resonance, ballistic);
  the six smallest re-evaluated at 50 digits coincide with double to all printed digits (practical).
  (A first pass that did not drop unresolved grids gave 1.2·10⁻¹¹ at aliased Cattaneo points; discarded.)
* Exact strata (50 digits, raw s/s₁): generic (2, 0.5, 0.1, 0.3): s₅ = 1.2·10⁻² (rank 5);
  resonance: s₅ = 4·10⁻⁵⁷ (rank 4); resonance with b = b_s and R > 0: rank 4 — the resistance makes the
  matched Fourier film visible again; resonance with b = b_s and R = 0: s₄, s₅ ~ 10⁻⁵² (rank 3);
  R = 0 boundary (one-sided, column (ξ₁/b)∂/∂R): rank 5 (s₅ = 1.5·10⁻²); b = b_s, R = 0, off resonance:
  rank 5; Cattaneo + R: rank 5 (2.1·10⁻³); τ_R = 0 + R: rank 5 (1.4·10⁻³); Fourier film + R: rank 4
  (three Fourier parameters plus the resonance direction).
* Near resonance (β = 2, τ̄ = 0.5ξ₁², R = 0.3): s₅ = 7.4·10⁻²|d|, exponent 0.9999.
* R → 0: raw s₅ ∝ R in log-R coordinates (coordinate artefact; equilibrated value stays 0.08).
  R → ∞: raw s₅ stays at 0.31 because the full band reaches down to ω ~ 1/(R b_s)², where R is still
  visible against the substrate impedance; in a fixed band the substrate and R decouple as R → ∞
  (adiabatic back) and the ln R column vanishes.
* N = 3 frequencies: a non-zero 5×5 minor certified in ball arithmetic at 100/100 random points.
* Adversarial (24 starts, |ln τ_R/τ_ℓ| ≥ 1 hard): 13 minima on the box, 9 on the constraint, 2 interior with
  s₅ ≥ 10⁻²·⁷; lowest 10⁻⁵·¹⁸ at the box corner (t_ℓ → 10³, R ~ 10²–10³).

Thin-film hierarchy with R (`s02_thin_R.py`): first-order coefficient c₁ = √p (R + ξ₁/β − βξ₁), i.e.
R + e/λ − eρc/b_s² in physical units; the cumulative rank grows by one per order and reaches 5 at ε⁵
(band exponents 1/2, 1, 3/2, 2, 5/2).

Verdict B4: SURVIVES (generic rank 5; exact degenerate set {τ_R = τ_ℓ}, rank 4, and
{τ_R = τ_ℓ, b = b_s, R = 0}, rank 3; the face R = 0 is not degenerate). Practically ill-conditioned in the
same places as B1, and in thin films, where R enters at leading order only through R + e/λ − eρc/b_s².

## 6. Practical non-identifiability in a finite band: scaling laws (`practical_bands.py`, `s01_thin_hierarchy.py`)

Band [ω_b/2·10⁴, ω_b] (the 4.3 decades of 10 kHz–200 MHz); ŝ_k = s_k/√(ln 2·10⁴): a unit step in the
log-parameters along the k-th singular vector changes Z by ŝ_k (RMS relative over the band).

**Thin film** (ω_b ξ₁² → 0, times fixed relative to ξ₁²): ŝ_k ∝ (ω_b ξ₁²)^{κ_k} with fitted
κ = 0.495, 0.997, 1.504, 2.499 (β = 2; GK in either order, Cattaneo, τ_R = 0, and β = 0.5 GK give the same
exponents to ± 0.01). Symbolic expansion in ε (√p → ε√p, τp → ε²τp) gives the order at which the
cumulative rank of the coefficient map reaches 1, 2, 3, 4: ε¹, ε², ε³, ε⁵, i.e. κ = 1/2, 1, 3/2, 5/2
exactly. Leading coefficients (W = b_s√p Z):
c₁ = −√p ξ₁(β² − 1)/β, c₂ = p ξ₁²(β² − 1), c₃ = −p^{3/2} ξ₁[ξ₁²(3β⁴ − 4β² + 1) − 3(τ_R − τ_ℓ)]/(3β);
the times enter at ε³ only through τ_R − τ_ℓ, the remaining time combination at ε⁵. The weakest direction
is the common shift dτ_R = dτ_ℓ (singular vector (0, 0, −0.196, −0.981) in log coordinates for τ_R = 5τ_ℓ).
At b = b_s: c₁ = c₂ = 0, c₃ = p^{3/2} ξ₁(τ_R − τ_ℓ), c₄ = −p² ξ₁²(τ_R − τ_ℓ); cumulative-rank orders
1, 3, 4, 5, i.e. κ = 1/2, 3/2, 2, 5/2 (numerical: 0.496, 1.493, 1.994, 2.499). At leading order only the
product ξ₁(τ_R − τ_ℓ) is seen: the thickness trades against the time difference (the mechanism of the
b = b_s pairs of §3.3).

**Low ωτ at fixed band position** (ω_b ξ₁² = 2, β = 2, τ_R/τ_ℓ = 5, times scaled by ε): ŝ₁, ŝ₂ constant,
ŝ₃ ≈ 0.082 ω_b τ_R, ŝ₄ ≈ 6.8·10⁻³ (ω_b τ_R)² (exponents 1.000, 2.001). Near the resonance inside this band
(β = 2, τ̄ = 0.05ξ₁², ω_b τ̄ = 0.1): ŝ₄ = 3.0·10⁻⁴ |d|, exponent 1.000.

**Thick film** (band bottom ω_a ξ₁² → ∞): ŝ₄ ∝ e^{−2 Re A(ω_a)} (fitted exponent 0.98), singular vector
→ ξ₁ alone; e.g. ŝ₄ = 4.6·10⁻²⁵ at ω_a ξ₁² = 316.

**Reading for the 500 nm AlN/sapphire scenario** (β ≈ 2, ξ₁² ≈ 3 ns, band ω ξ₁² ≈ 2·10⁻⁴ – 4): this is the
mid band. For τ_R in the 10–100 ps range (ω_b τ_R ≈ 0.013–0.13), the laws above (τ_R/τ_ℓ = 5) give
ŝ₃ ≈ 10⁻³–10⁻² and ŝ₄ ≈ 10⁻⁶–10⁻⁴:
the separate values of τ_R and τ_ℓ are not resolvable at 0.1°/1 % noise; at best τ_R − τ_ℓ is, marginally.
These are local Fisher-type statements with the four film parameters free and no instrument or interface
nuisance; Workstream A must use its own nuisance sets and the actual τ_R(500 nm).

## 7. Attacks tried, attacks not performed, limitations

Tried: random log-uniform scans (≈ 40 000 points over β ∈ [10⁻³, 10³], t ∈ [10⁻⁵, 10⁵], including both
boundary strata and b = b_s); adversarial minimisation of the smallest singular value with hard distance
constraints; 50/80-digit recomputation of every small value; ball-arithmetic certificates of non-zero
minors (300 points); exact strata ranks at 50 digits (180 points); near-resonance and near-invisible-film
power laws; ballistic, high/low-effusivity, thin and thick limits; multistart global search for distinct
parameter sets with equal response in 4 bands × 4 strata × 8 targets, explicit τ_R ↔ τ_ℓ swap, refinement
and band widening; the square N = 2 system (300 starts × 5 cases, 60-digit polishing); N = 3; interface
resistance (scan, strata, certificates, adversarial); complex rank versus real rank (equal at all points
tested); the physical formula against the project forward model.

Not performed / limitations:
* No proof. The analytic obstruction of §1.3 and the asymptotic argument of §3.5 are sketches; the
  function-level claims rest on numerical evidence plus these sketches. Hidden degeneracies on sets of
  measure zero, or outside the scanned ranges (β beyond 10^±3, t beyond 10^±5), are not excluded.
* Ball-arithmetic certificates are pointwise (300 points), not over parameter regions.
* The multistart searches are heuristic: in the thin band some minima were reached only through the swap
  path, not from the 40 random starts, so the misfits of §3.2–3.3 are upper bounds on the best attainable.
* Full-band conditioning numbers depend on the declared window (§2.5 e); only ranks and band-limited
  numbers are meaningful.
* Cattaneo films with t_R ≳ 10³ have resonance combs needing > 6·10⁴ nodes; 152 (B1) and 143 (B4) such
  points were excluded from the resolved scans (their rank is covered by the Cattaneo-stratum analysis at
  smaller t_R and by the analytic argument).
* Only the 1D model of note 7 with the "local" interface temperature; the entropic or kinetic temperature
  conventions, finite beams, transducer layers and instrument nuisances were not attacked.
* N = 2 global counterexample: one (θ, pair) case with multiplicity 6; how often multiple roots occur over
  (θ, ω₁, ω₂) was not mapped (4 other cases showed a single root from 300 starts).

Side note (non-physical form). The swapped formula of the original brief gives qualitatively similar rank
statistics (`b1_scan.py brief`: interior median equilibrated s₄ = 0.16, minimum 2·10⁻⁵ near resonance;
τ_R = 0 stratum minimum 0.046). It is not the physical model and is not used elsewhere.

## 8. Files (`derivations/B-counterexamples/`)

| file | content |
|---|---|
| `zmodel.py` | impedance, analytic Jacobian (double and mpmath), resolved frequency grids |
| `t00_verify_model.py` | checks against `src/forward_model.response`, 50-digit derivatives, limits |
| `s00_identities.py`, `s01_thin_hierarchy.py`, `s02_thin_R.py` | symbolic identities, thin-film hierarchies |
| `b1_scan.py`, `b1_mpcheck.py`, `b1_adversarial.py`, `b1_strata_mp.py`, `b1_resonance_scaling.py` | B1 |
| `b2_injectivity.py`, `b2_refine_pairs.py`, `practical_pairs.py` | B2 |
| `b3_two_freq.py`, `b3_singular_pairs.py`, `b3_global_N2.py`, `b3_global_N2_roots40.py` | B3 |
| `b4_resistance.py` | B4 |
| `practical_bands.py`, `practical_bands_beta1.py` | band-limited laws |
| `run_chain2.sh` | sequential batch used for the second series |
| `data/*.json, *.npy, *.log` | raw results of every run |
