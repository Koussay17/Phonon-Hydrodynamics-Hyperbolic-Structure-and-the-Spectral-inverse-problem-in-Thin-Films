# Branch B-S — computer-algebra route: rank, degenerate locus, gauge and expansions of the 1D front-face impedance

Campaign `20261111-fdtr-gauge-memory`, Workstream B (00-question.md), assumptions B1–B3 (02-assumptions.md).
Scripts: `derivations/B-symbolic/` (Python with SymPy, mpmath, python-flint/Arb, Pint; Wolfram Language `.wl` files).
Raw outputs: `derivations/B-symbolic/results/`. All numbers quoted below are copied from those files.

Wolfram computations ran on a Wolfram Language 15.0.1 kernel through remote evaluation (the local Wolfram Engine
reports a licence/activation problem); the `.wl` files are self-contained and run unchanged with `wolframscript -file`.
Compute limits respected: single-threaded SymPy/mpmath/Arb jobs, NumPy capped at two threads, grid arrays sized below
1 GB.

---

## 0. Summary

| # | Question | Result | Status |
|---|---|---|---|
| 0 | Closure formula of the brief | The brief's expression has numerator and denominator swapped: it equals 1/(B²Z). Correct: Z = (1/B)(B + S tanh A)/(S + B tanh A), S = b_s√p | **DISCREPANCY FOUND** (corrected; the correction was acknowledged in the campaign during the run; every result below uses the correct form) |
| 1 | Generic rank, finite N | Real Jacobian [Re; Im] (2N × 4): rank 2 for N = 1, rank 4 for N ≥ 2. Complex Jacobian (N × 4): rank min(N, 4). **Minimal N for rank 4: N = 2 (real), N = 4 (complex)** | VERIFIED (50-digit SVD, two derivative methods, Wolfram; Arb-certified minors at 40 rational points; generic by analyticity) |
| 2 | Degenerate locus Σ (function level, all ω > 0) | **Σ = {τ_R = τ_ℓ}** (contains the Fourier point). Rank 4 off Σ; rank 3 on Σ ∩ {b ≠ b_s}, null vector (0,0,1,1); rank 2 on Σ ∩ {b = b_s}, null vectors (1,0,0,0), (0,0,1,1). The candidates b = b_s, τ_R = 0, τ_ℓ = 0 are **not** degenerate by themselves | DERIVED, two independent proofs (§3.3, §3.4), all symbolic steps machine-checked (SymPy + Wolfram), certified numerics on every stratum, numerical search found nothing else. Not yet audited |
| 2' | N = 2 minor | det[Re r₁; Im r₁; Re r₂; Im r₂] = (τ_R − τ_ℓ)/(4ξ₁b) · \|W(iω₁)\|² \|W(iω₂)\|² · D̃ ; the only parameter-only factor is (τ_R − τ_ℓ)/(ξ₁b) | VERIFIED (exact rational evaluation + Wolfram symbolic) |
| 3 | Gauge | e∂_e + 2a∂_a + λ∂_λ + ℓ∂_ℓ annihilates Z; kernel of the composite map (e, a, λ, ℓ, τ_R) ↦ Z is exactly this direction off Σ (dimension 2 on Σ, 3 on Σ ∩ {b = b_s}); the only unit rescaling acting on the parameters and fixing (ω, b_s, Z) is the length rescaling | VERIFIED |
| 4 | Expansions | Small p: Puiseux series in √p, coefficients z₋₁ … z₅ in closed form; relaxation times enter first at order p through τ_R − τ_ℓ only, and separate at order p². Large p: power series of 1/B (determines b, τ_R, τ_ℓ) plus exponentially small reflections (the only carriers of ξ₁, b_s, R) for τ_R, τ_ℓ > 0; stratum-specific forms for τ_ℓ = 0, τ_R = 0, τ_R = τ_ℓ | VERIFIED (SymPy and Wolfram: exact agreement of the small-p coefficients, identical large-p coefficients; 50-digit truncation tests, all ratios → 1) |
| 5 | Interface resistance R | **Generic rank 5**; minimal N = 3 (real), 5 (complex). Σ_R = {τ_R = τ_ℓ}: rank 4 if R > 0; rank 3 if R = 0, with the extra null vector (1, 0, 0, 0, −(1 − b²/b_s²)/b) (thickness ↔ resistance trade-off) | DERIVED + VERIFIED (same methods) |
| + | By-product (claim B2) | θ ↦ Z(iω, ·) is **injective on {τ_R ≠ τ_ℓ}**, by explicit inversion of the small-p coefficients; fibres on Σ: the half-lines {(ξ₁, b, τ, τ), τ ≥ 0}; on Σ ∩ {b = b_s}: {(ξ₁, b_s, τ, τ)}, ξ₁ > 0, τ ≥ 0 | DERIVED (constructive), symbolic steps verified, end-to-end numerical inversion; not yet audited |

---

## 1. Model, conventions and the closure (task 0)

**DEFINITIONS.** θ = (ξ₁, b, τ_R, τ_ℓ) ∈ (0,∞)² × [0,∞)², substrate effusivity b_s > 0 known, optional interface
resistance R ≥ 0. With principal square roots,

ρ(p) = √((1 + τ_R p)/(1 + τ_ℓ p)), A = ξ₁√p ρ (= σe), B = b√p/ρ (= λ_eff σ), S = b_s√p,
Z_s = R + 1/S, w = B Z_s, t = tanh A, Γ = (w − 1)/(w + 1), E = e^(−2A).

On p = iω, ω > 0: arg[p(1 + τ_R p)/(1 + τ_ℓ p)] = π/2 + arctan ωτ_R − arctan ωτ_ℓ ∈ (0, π), so the single square root
of the brief equals √p · ρ, and A, B are analytic on D = ℂ \ (−∞, 0] (the Möbius map p ↦ (1+τ_R p)/(1+τ_ℓ p) sends
non-real p to non-real values). Re A > 0 and Im A > 0 on iℝ₊.

**Closure.** With the quadrupole M = [[cosh A, sinh A/B], [B sinh A, cosh A]] and θ₀ = (M₀₀Z_s + M₀₁)/(M₁₀Z_s + M₁₁):

Z = (1/B)(w + t)/(1 + w t) = (1/B)(1 + ΓE)/(1 − ΓE); with R = 0: **Z = (1/B)(B + S t)/(S + B t)**.

Thin film (t → 0): Z → Z_s = 1/S (bare substrate). Thick film (t → 1): Z → 1/B.

**DISCREPANCY.** The brief's expression Z_b = (1/B)(B t + S)/(B + S t) satisfies identically **Z_b = 1/(B² Z)** and
Z_b − Z = −(B − S)(B + S)/[B(B + S t)(S + B t) cosh² A]. It coincides with Z only where B² ≡ S² (b = b_s and τ_R = τ_ℓ)
or asymptotically for a thick film; its thin-film limit is S/B² instead of 1/S, which is unphysical.

Evidence (results/00_closure.txt, 00_closure_wl.txt):
- SymPy: quad − closed = 0, Z_b Z B² − 1 = 0, the multiple-reflection form, both limits (symbolic, exact).
- Wolfram: the same three identities = 0; limits {1/s, s/bb²}.
- mpmath, 50 digits, 40 random rational points on iω (incl. τ = 0 boundaries): |Z_quad − Z|/|Z| ≤ 4.3e-51; with R
  (20 points) ≤ 2.3e-51; |Z_b − Z|/|Z| ranges from 6.6e-13 (a thick-film point, ω ξ₁² ≈ 58) to O(1). At the Wolfram
  spot point the relative difference is 3.018 in both systems.
- The project's own implementation (src/forward_model.py `response`, double precision, five AlN/sapphire-like stacks
  incl. Cattaneo, GK, R ≠ 0, 10 kHz–200 MHz): agreement with the closed form ≤ 5.0e-16; distance to Z_b ≥ 4.1e-3.

---

## 2. Structure of the Jacobian (key lemma)

**LEMMA 1 (columns).** Let U := A ∂_A Z, V := B ∂_B Z (at fixed Z_s), Y := ∂_{Z_s} Z, g_R = p/(2(1 + τ_R p)),
g_ℓ = p/(2(1 + τ_ℓ p)). Then

U = A(1 − t²)(1 − w²)/(B(1 + wt)²), V = −t(1 + 2wt + w²)/(B(1 + wt)²), Y = (1 − t²)/(1 + wt)²,

∂_{ξ₁}Z = U/ξ₁, ∂_b Z = V/b, ∂_{τ_R}Z = g_R (U − V), ∂_{τ_ℓ}Z = −g_ℓ (U − V), ∂_R Z = Y.

*Proof.* A is linear in ξ₁ and free of b; B is linear in b and free of ξ₁; ∂ ln A/∂τ_R = −∂ ln B/∂τ_R = g_R and
∂ ln A/∂τ_ℓ = −∂ ln B/∂τ_ℓ = −g_ℓ (from ln A = ln ξ₁ + ½ ln p + ½ ln X, X = ρ²); chain rule. ∎

Verification: SymPy (results/02_sigma_symbolic.txt, A1–A2: all differences cancel to 0; A3: the assembled columns
agree with direct SymPy derivatives of the full expression to 6.0e-51 at 30 random points, incl. τ_ℓ = 0 and R = 0);
Wolfram (results/02_sigma_wl.txt, i1–i2: `Simplify` → 0 on p > 0, residual 0 to ~196 digits at p = iω).

Immediate consequences: on τ_R = τ_ℓ, col₃ + col₄ = (g_R − g_ℓ)(U − V) = 0 (rank ≤ 3); if moreover b = b_s (and R = 0),
w ≡ 1 so U ≡ 0, col₁ = 0 (rank ≤ 2).

---

## 3. Function-level rank and the degenerate locus Σ (task 2)

### 3.1 Reduction
Z(p; ξ₁, b, τ_R, τ_ℓ; b_s, R) = (ξ₁/b_s) Ẑ(p ξ₁²; b/b_s, τ_R/ξ₁², τ_ℓ/ξ₁², R b_s/ξ₁), Ẑ := Z(·; 1, β, T_R, T_ℓ; 1, R̂).
The change of variables is a diffeomorphism for ξ₁ > 0 and p ↦ pξ₁² maps iℝ₊ onto itself, so the function-level rank
depends only on (β, T_R, T_ℓ[, R̂]); this is used for the numerical search.

### 3.2 Statement
**THEOREM 1 (DERIVED; not yet audited).** Fix b_s > 0 and R = 0. For θ ∈ (0,∞)² × [0,∞)² let J(θ) be the real-linear
map c ↦ Σ_k c_k ∂_kZ(iω; θ), ω > 0 (equivalently on any subset of iℝ₊ accumulating in iℝ₊ or at ω = 0). Then
(i) rank J = 4 iff τ_R ≠ τ_ℓ;
(ii) rank J = 3, ker J = span{(0,0,1,1)} on {τ_R = τ_ℓ, b ≠ b_s};
(iii) rank J = 2, ker J = span{(1,0,0,0), (0,0,1,1)} on {τ_R = τ_ℓ, b = b_s}.
The ranks over ℝ and ℂ coincide (the ∂_kZ are real on p > 0).

### 3.3 Proof A (transcendence)
Assume Σ_k c_k ∂_kZ = 0 on iℝ₊; by the identity theorem it holds on D. By Lemma 1, multiplying by B(1 + wt)²,

e₀ + e₁ t + e₂ t² = 0, e₀ = (c₁/ξ₁ + h) A(1 − w²), e₁ = −(c₂/b − h)(1 + w²), e₂ = −e₀ − 2(c₂/b − h) w,

with h := c₃g_R − c₄g_ℓ (SymPy A4: coefficients checked exactly). The e_j lie in the algebraic function field
K = ℂ(p, √p, ρ).
**LEMMA 2.** t = tanh A is transcendental over K. *Proof.* Otherwise e^{2A} = (1 + t)/(1 − t) is algebraic over ℂ(p).
A is a non-constant algebraic function with a pole reachable inside D along which Re A → +∞ (p → +∞ if τ_R > 0 or
τ_R = τ_ℓ; p → −1/τ_ℓ from Im p > 0 if τ_R = 0 < τ_ℓ). There |e^{2A}| grows faster than any power, whereas an algebraic
function grows at most like a power (Puiseux). ∎
Hence e₀ = e₁ = e₂ = 0. Since w² → b²/b_s² at p → 0, 1 + w² ≢ 0, so h ≡ c₂/b; h(0) = 0 gives c₂ = 0 and h ≡ 0, i.e.
c₃(1 + τ_ℓ p) ≡ c₄(1 + τ_R p): c₃ = c₄ = 0 unless τ_R = τ_ℓ, where c₃ = c₄ (determinant of this 2 × 2 system: τ_R − τ_ℓ,
SymPy C). Finally e₀ = (c₁/ξ₁)A(1 − w²) ≡ 0 forces c₁ = 0 unless w² ≡ 1 ⇔ b²(1 + τ_ℓ p) ≡ b_s²(1 + τ_R p) ⇔ b = b_s and
τ_R = τ_ℓ. This is exactly (i)–(iii). ∎

### 3.4 Proof B (exact, from the small-p Puiseux coefficients; independent of Lemma 2)
qZ is analytic in q = √p near 0 jointly with θ, so Z = Σ_{n≥−1} z_n(θ) qⁿ with z_n rational in θ (computed exactly,
§6.1), and any function-level relation Σc_k∂_kZ ≡ 0 implies Σc_k∂_k z_n = 0 for every n. Exact minors
(results/04_expansions.txt, 06_inversion_series.txt, 06b_sigma_minors.txt):
- det ∂(z₀, z₁, z₂, z₄)/∂θ = −2ξ₁⁴(b² − b_s²)(τ_R − τ_ℓ)/(b⁴b_s³) (b ≠ b_s, τ_R ≠ τ_ℓ ⇒ rank 4);
- det ∂(z₀, z₂, z₃, z₄)/∂θ at b = b_s = −2ξ₁⁴(τ_R − τ_ℓ)²/b_s⁵ (b = b_s, τ_R ≠ τ_ℓ ⇒ rank 4);
- on τ_R = τ_ℓ: det ∂(z₀, z₁, z₂)/∂(ξ₁, b, τ_R) = 2ξ₁³(b² − b_s²)/(b³b_s³) (rank ≥ 3 if b ≠ b_s);
- on τ_R = τ_ℓ, b = b_s: det ∂(z₀, z₂)/∂(b, τ_R) = −2ξ₁²/b_s³ (rank ≥ 2).
With the upper bounds of Lemma 1 this reproves (i)–(iii) by exact rational arithmetic. (det ∂(z₀, z₁, z₂, z₃)/∂θ ≡ 0:
z₋₁ … z₃ depend on the times only through τ_R − τ_ℓ, so four coefficients are not enough; z₄ is needed.)

### 3.5 Strata (all verified)

| Stratum | Rank | Kernel | Upper bound | Lower bound (certificate) |
|---|---|---|---|---|
| τ_R ≠ τ_ℓ (incl. τ_R = 0 or τ_ℓ = 0, and b = b_s) | 4 | 0 | — | exact minor (§3.4); Arb 4 × 4 minors at 40 points (§4) |
| τ_R = τ_ℓ = τ ≥ 0, b ≠ b_s (incl. Fourier τ = 0) | 3 | (0,0,1,1) | col₃ + col₄ = 0 | exact 3 × 3 minor; Arb 3 × 3 minors at 6 points |
| τ_R = τ_ℓ, b = b_s | 2 | (1,0,0,0), (0,0,1,1) | col₁ = 0, col₃ + col₄ = 0 | exact 2 × 2 minor; Arb 2 × 2 minors at 5 points |

Physical reading: on Σ the film is exactly a Fourier layer (Fourier resonance), the common time is invisible; with
matched effusivity as well, Γ ≡ 0 and Z = 1/(b√p) carries no thickness information.

Verification of the strata (results/02_sigma_ball.txt, 02_sigma_wl.txt, 02_sigma_symbolic.txt):
- Arb (320-bit, forward-mode differentiation in ball arithmetic, no symbolic derivative): certified minors of the
  expected size at every stratum point, next-size minors contain 0, |col₃ + col₄| ≤ 1e-93 and |col₁| ≤ 1e-94
  (b = b_s) on the strata.
- Wolfram: i3 = (∂_{τ_R} + ∂_{τ_ℓ})Z|_{τ_ℓ=τ_R} → 0, i4 = ∂_{ξ₁}Z|_{τ_ℓ=τ_R, b=b_s} → 0 (symbolic); normalised singular
  values: generic 4 nonzero; τ_R = τ_ℓ: σ₄/σ₁ = 8.8e-78; τ_R = τ_ℓ, b = b_s: σ₃/σ₁ = 5.6e-116; Fourier 3 resp. 2;
  b = b_s off Σ, Cattaneo, τ_R = 0: rank 4; null vectors confirmed to 40 digits.
- SymPy, 50 digits: |∂_{τ_R}Z + ∂_{τ_ℓ}Z| on τ_R = τ_ℓ = 0.0 at 20 random points.

### 3.6 Numerical search for other degeneracies (results/02_sigma_numeric.txt, 02_sigma_corner.txt)
Grid β ∈ 10^[−2,2] (17 values, incl. β = 1), T_R, T_ℓ ∈ {0} ∪ 10^[−3,3] (26 values each): 11 492 points, 65 frequencies
P = ωξ₁² ∈ [1e−4, 1e4], real 130 × 4 Jacobian.
- On the diagonal: max σ₄/σ₁ = 4.8e-16 (rank ≤ 3); β ≠ 1: min σ₃/σ₁ = 1.3e-4; β = 1 with natural column scaling:
  σ₃/σ₁ ≤ 2.9e-16 (rank 2); at 40 digits σ₃/σ₁ ~ 1e-42.
- Near the diagonal, T_ℓ = T_R(1 + ε): σ₄/σ₁ ∝ ε, fitted slopes 0.85–1.00 (nine (β, T) pairs), as required by the
  factor (τ_R − τ_ℓ).
- Off the diagonal (|log₁₀ T_R/T_ℓ| ≥ 0.25, and the boundaries T = 0): min σ₄/σ₁ = 1.19e-7 at β = 31.6, T_R = 10³,
  T_ℓ = 0; recomputed at 50 digits (1.1865e-7, identical): genuinely nonzero. τ_R = 0 boundary: min 8.5e-5.
  No interior zero was found.
- Artefact identified and removed: on {β = 1, T_R = T_ℓ} the ξ₁ column is identically zero and column equilibration
  rescales its round-off to unit norm (σ₃/σ₁ = 0.53 in the equilibrated double-precision run); natural scaling and
  40-digit arithmetic give rank 2 as proved.

### 3.7 Conditioning remarks (not degeneracies)
- Near Σ the fourth singular value is O(|τ_R − τ_ℓ|) (factor of the determinant, §4.2).
- Cattaneo corner, τ_R ≳ ξ₁² (T_R = 1 … 10⁵, τ_ℓ = 0): in the weakest direction (δ ln ξ₁, δ ln b, δ ln τ_R) the ratio
  δ ln ξ₁ : δ ln τ_R is ½ : −1 at every point computed (transit time ξ₁√τ_R = e/v fixed); the b component tends to −½
  as T_R → ∞ (already at all T_R for β = 31.6; between −0.19 and −1.1 for β ≤ 1 at moderate T_R), so asymptotically the
  wave impedance b/√τ_R and the heat capacity ξ₁b are fixed as well and only the damping per transit ξ₁/√τ_R → 0
  changes. Equilibrated σ₄/σ₁ ≈ 4e-4 (T_R = 1) to 1e-8 (T_R = 10⁵) for the band [1e−4, 1e4] (β = 31.6); the rank is 4
  (Theorem 1; certified Cattaneo points 24–27 of §4).

---

## 4. Finite frequency sets (task 1)

### 4.1 Ranks (results/01_rank_mp.txt, 01_rank_ball.txt, 01_rank_wl.txt, 01_compare_wl.txt)
40 exact rational points: 24 random (ξ₁ ∈ [0.1, 5], b, b_s ∈ [0.1, 10], T_R, T_ℓ ∈ [0.01, 20], ωξ₁² ∈ [0.01, 100]),
4 Cattaneo (τ_ℓ = 0), 4 with τ_R = 0, 4 with b = b_s, 4 AlN/sapphire SI points (e = 500 nm, λ = 321 or 60 W/mK,
C = 2.41e6 J/m³K, sapphire 35 W/mK, 3.03e6 J/m³K, τ_R = 0.1–10 ns, τ_ℓ = 0–3 ns, 8 frequencies 10 kHz–200 MHz).
- Derivatives: (a) SymPy symbolic, (b) mpmath numerical differentiation at 100 digits of the quadrupole product
  (no symbolic input). Max relative entry difference 2.3e-49; singular values agree to ~1e-50.
- 50-digit SVD of the column-equilibrated matrices. Every point:

| N | real rank (2N × 4) | complex rank (N × 4) | smallest σ_min/σ₁ (real / complex) over the 40 points |
|---|---|---|---|
| 1 | 2 | 1 | 0.35 / — |
| 2 | 4 | 2 | 1.0e-5 / 0.13 |
| 3 | 4 | 3 | 7.7e-5 / 3.6e-4 |
| 4 | 4 | 4 | 7.1e-4 / 3.2e-5 |
| 8 | 4 | 4 | 1.4e-3 / 1.9e-4 |

  An exact zero would appear below 1e-45 at this precision: the gap between the smallest nonzero singular value and the
  noise floor is ≥ 40 decades at every point.
- Wolfram (8 of the points, symbolic D of the quadrupole product, N[·, 60], `SingularValueList`): agreement with mpmath
  to 2.2e-19 (the printed precision) for all 64 matrices.
- Arb certificates (forward-mode AD in 320-bit balls) at all 40 points: real N = 1 2 × 2 minor, real N = 2 4 × 4
  determinant, complex N × N minors for N = 1…4, all with balls excluding 0 (radii ~1e-91); the identity
  det_real = −¼ det_complex(±iω) holds in ball arithmetic at every point.

**Minimal N.** Real Jacobian: N = 2 (N = 1 gives at most 2). Complex Jacobian: N = 4 (N ≤ 3 gives at most N).
Genericity: the 4 × 4 determinant is real-analytic in (θ, ω₁, ω₂) on a connected domain and certified nonzero at 40
points, so its zero set is a proper analytic subset (nowhere dense, measure zero).

### 4.2 Factorisation of the N = 2 minor (task 2, "common factors")
For four arbitrary points p₁…p₄ (SymPy B1: exact rational evaluation at 200 random points, 0 failures; Wolfram:
`Together` of the difference = 0, symbolic):

det[∂_kZ(p_m)] = −(τ_R − τ_ℓ)/(ξ₁ b) · Π_m W(p_m) · det[1, V/W, g_R, κ](p_m), W = U − V, κ = p²/(2(1 + τ_R p)(1 + τ_ℓ p)),

and det[Re r₁; Im r₁; Re r₂; Im r₂] = −¼ det[r₁; r̄₁; r₂; r̄₂] (SymPy B2, generic complex rows), with r̄(iω) = r(−iω). Hence

**det J_real(ω₁, ω₂) = (τ_R − τ_ℓ)/(4ξ₁b) · |W(iω₁)|² |W(iω₂)|² · D̃(iω₁, −iω₁, iω₂, −iω₂).**

The only factor depending on θ alone is (τ_R − τ_ℓ)/(ξ₁b): Σ is the only parameter-only degeneracy of the two-frequency
design. The remaining zeros are frequency-design zeros (W(iω_j) = 0 or D̃ = 0).
N = 2 on the grid of §3.6 (Hadamard ratio |det|/Π‖col‖, best of 8 random pairs P₁, P₂ ∈ 10^[−2,2]): min off Σ
8.1e-9 (corner β = 100, T ~ 10³, where the pairs lie far above 1/T), Cattaneo boundary 7.4e-4, τ_R = 0 boundary 2.4e-7;
on Σ ≤ 8.9e-17.

### 4.3 What is proved and what is supported
- Complex N = 4 and real N = 4: for **every** θ ∉ Σ, generic frequency sets give rank 4 (Theorem 1 + the elementary
  lemma: functions linearly independent on a set E admit points of E with a nonzero determinant; then analyticity).
- Real N = 2 (and 3): rank 4 for generic (θ, ω₁, ω₂) (certified + analyticity) and at every grid point off Σ
  (numerical). "For every θ ∉ Σ" with N = 2 is NUMERICALLY SUPPORTED, not proved. A proof route: the identity
  D̃ ≡ 0 on the conjugate-pair slice reduces to an algebraic relation between tanh A(p) and tanh A(−p), excluded by
  Ax's theorem (functional Schanuel) since A(p), A(−p) are ℚ-linearly independent modulo constants — a sketch only.

---

## 5. Gauge (task 3)  (results/03_gauge.txt, 03_gauge_wl.txt)

Physical parameters q = (e, a, λ, ℓ, τ_R), ρc = λ/a eliminated: ξ₁ = e/√a, b = λ/√a, τ_ℓ = κℓ²/a (κ = 3 historical,
κ = 1 + α in general).
- D = e∂_e + 2a∂_a + λ∂_λ + ℓ∂_ℓ annihilates ξ₁, b, τ_R, τ_ℓ exactly (SymPy), hence Z. On the full composite
  expression: Wolfram `Simplify` → 0; SymPy 50-digit check at 25 random points: |DZ|/|e∂_eZ| ≤ 1.4e-41. Without the
  ℓ term the residual is ≥ 0.026 relative (SymPy) and 3.9 at the Wolfram spot point — the ℓ term is required.
- Finite form: (e, a, λ, ℓ) → (μe, μ²a, μλ, μℓ) leaves (ξ₁, b, τ_ℓ) unchanged; Z unchanged (Wolfram symbolic 0; SymPy 0
  at 20 random points). Equivalently in (e, λ, ρc, ℓ, τ_R): kernel (e, λ, −ρc, ℓ, 0), i.e. ρc → ρc/μ.
- dφ/dq (4 × 5) has rank 4 with kernel spanned by (e, 2a, λ, ℓ, 0) (SymPy, Wolfram).
- Composite map q ↦ (Z(iω_j))_{j≤8}, 50 digits, log-scaled columns: off Σ exactly one zero singular value (~6e-52, the
  others ≥ 4.7e-3); on τ_R = τ_ℓ two (~1e-51); on τ_R = τ_ℓ, b = b_s three. So **off Σ the gauge is the only
  non-identifiability of the physical parameters**; on Σ add the motion along the diagonal (τ_R and ℓ together with
  τ_R = κℓ²/a); on Σ ∩ {b = b_s} add the thickness direction.
- Dimensional analysis (Pint for the dimensions, SymPy and, independently, hand-written exponents in Wolfram for the
  linear algebra). Exponents (m, kg, s, K): e (1,0,0,0); a (2,0,−1,0); λ (1,1,−3,−1); ρc (−1,1,−2,−1); ℓ (1,0,0,0);
  τ_R (0,0,1,0); b_s (0,1,−5/2,−1); ω (0,0,−1,0); Z (0,−1,3,1). Unit rescalings fixing ω, b_s and Z: the 2-parameter
  group {T = 1, Θ = M}; its action on ln(e, a, λ, ρc, ℓ, τ_R) is (1, 2, 1, −1, 1, 0)·ln L (the gauge) and 0 for M.
  **No other dimensional direction acts on the parameters.**

---

## 6. Expansions (task 4)  (results/04_expansions.txt, 04_expansions_wl.txt, 04_compare_wl.txt)

### 6.1 Small p (q = √p): Z = Σ_{n≥−1} z_n qⁿ (convergent), R = 0
With R_f := ξ₁/b = e/λ (film resistance) and C_f := ξ₁b = ρc e (film heat capacity per area), δ := τ_R − τ_ℓ:

| term | coefficient | determines |
|---|---|---|
| p^(−1/2) | z₋₁ = 1/b_s | nothing about the film (substrate) |
| p⁰ | z₀ = ξ₁(1/b − b/b_s²) = R_f − C_f/b_s² | R_f − C_f/b_s² |
| p^(1/2) | z₁ = −(ξ₁²/b_s)(1 − b²/b_s²) = −(C_f/b_s) z₀ | C_f (if b ≠ b_s); with z₀: (ξ₁, b) |
| p¹ | z₂ = ξ₁δ/b − ξ₁³(3b⁴ − 4b²b_s² + b_s⁴)/(3b b_s⁴) | τ_R − τ_ℓ (first appearance of the times) |
| p^(3/2) | z₃ = −ξ₁²δ/b_s + ξ₁⁴(b² − b_s²)(3b² − 2b_s²)/(3b_s⁵) | τ_R − τ_ℓ again (ξ₁ when b = b_s) |
| p² | z₄ = −ξ₁τ_ℓδ/b + G₄(ξ₁, b, δ) | τ_ℓ(τ_R − τ_ℓ): first coefficient separating τ_R and τ_ℓ |

(z₄ and z₅ in full in the results file; (∂_{τ_R} + ∂_{τ_ℓ})z_n ≡ 0 for n ≤ 3, = −ξ₁δ/b for n = 4.)
With interface resistance: z₀ = R + ξ₁(1/b − b/b_s²), z₁ = −(ξ₁/b_s³)(2Rb b_s² − ξ₁(b² − b_s²)), further terms in the file.
SymPy and Wolfram give the same polynomials; the exact rational values of z₋₁ … z₅ at a test point coincide in both
systems (14/14, with and without R).

### 6.2 Large p along iω, by stratum
- **(G) τ_R, τ_ℓ > 0, τ_R ≠ τ_ℓ.** Z = 1/B + (2/B)ΓE/(1 − ΓE): power series of 1/B plus exponentially small terms,
  |E| = e^(−2 Re A), Re A ~ ξ₁√(τ_R/τ_ℓ)√(ω/2). With K = √(τ_R/τ_ℓ)/b, a_R = 1/τ_R, a_ℓ = 1/τ_ℓ:
  1/B ~ K p^(−1/2) [1 + (a_R − a_ℓ)/(2p) − (a_R − a_ℓ)(a_R + 3a_ℓ)/(8p²) + (a_R − a_ℓ)(a_R² + 2a_Ra_ℓ + 5a_ℓ²)/(16p³) + …].
  The first three coefficients determine b√(τ_ℓ/τ_R) (high-frequency effusivity), 1/τ_R − 1/τ_ℓ, and a_R + 3a_ℓ, hence
  (b, τ_R, τ_ℓ) when τ_R ≠ τ_ℓ. ξ₁, b_s and R appear only in the exponentially small reflection, leading term
  (2Γ_∞K/√p) exp(−2ξ₁√(τ_R/τ_ℓ)√p)·(1 + O(p^(−1/2))), A = ξ₁√(τ_R/τ_ℓ)[√p + (a_R − a_ℓ)/(2√p) + …],
  Γ_∞ = (b√(τ_ℓ/τ_R) − b_s)/(b√(τ_ℓ/τ_R) + b_s) for R = 0 (Γ → 1 if R > 0). Rate ξ₁√(τ_R/τ_ℓ) = e/√a_∞,
  a_∞ = κℓ²/τ_R the GK high-frequency diffusivity.
- **(F) τ_R = τ_ℓ (incl. Fourier).** 1/B = 1/(b√p) exactly: one power term; transseries
  Z = (1/(b√p))[1 + 2Γe^(−2ξ₁√p) + 2Γ²e^(−4ξ₁√p) + 2Γ³e^(−6ξ₁√p) + …], Γ = (b − b_s)/(b + b_s): determines b, then Γ, then ξ₁.
- **(C) Cattaneo, τ_ℓ = 0 < τ_R.** No decaying power series: with t = tanh A exact (bounded on iω since
  Re A → ξ₁/(2√τ_R) > 0), Z = t√τ_R/b + (1 − t²)/(b_s√p) + t(2b²(t² − 1) + b_s²)/(2b b_s²√τ_R p) − b²t²(t² − 1)/(b_s³τ_R p^(3/2)) + …,
  A = ξ₁√τ_R p + ξ₁/(2√τ_R) − ξ₁/(8τ_R^(3/2) p) + …. The leading oscillatory term carries √τ_R/b (wave impedance⁻¹),
  ξ₁√τ_R = e/v (transit time) and ξ₁/(2√τ_R) (damping per transit); the next ones bring b_s.
- **(N) τ_R = 0 < τ_ℓ.** A → A_∞ = ξ₁/√τ_ℓ: power series from p⁻¹,
  Z = coth A_∞/(b√τ_ℓ p) − b_s csch²A_∞/(b²τ_ℓ p^(3/2)) + c₄ p⁻² + c₅ p^(−5/2) + … (c₄, c₅ symbolic in
  results/04_expansions_wl.txt; SymPy agrees at the test point to 1.6e-33); unlike case (G), ξ₁ enters the power series:
  the first two coefficients already involve ξ₁, b, τ_ℓ (through A_∞, b√τ_ℓ) and b_s.

Numerical truncation tests on p = iω (50 digits): see §6.3.

### 6.3 Truncation tests (results/04b_truncation.txt)
Test point ξ₁ = 7/10, b = 23/10, b_s = 11/10, τ_R = 3/10, τ_ℓ = 1/8 (R = 0 or 1/5); coefficients from exact-rational
series evaluated at 60 digits, Z at 50 digits. Ratio = |error of the truncation| / |first omitted term|, which must
tend to 1:

| expansion | first omitted term | ratio along iω |
|---|---|---|
| small p, R = 0 | z₅ q⁵ | 0.911, 0.971, 0.991, 0.997, 0.9991 (ω = 1e−2 … 1e−6) |
| small p, R = 1/5 | z₅ q⁵ | 0.938, 0.980, 0.994, 0.998, 0.9994 |
| large p, generic | y₄ p^(−9/2) | 0.99997, 0.9999997, 0.999999997, 1.0 (ω = 1e3 … 1e6); \|Z − 1/B\|/\|Z\| = 2.9e-22 at ω = 1e3 (e^(−2 Re A) = 9.7e-22), at the 50-digit floor beyond |
| large p, Cattaneo | C₄(tanh A) u⁴, u = p^(−1/2) (oscillatory, \|tanh A\| between 0.58 and 1.74) | 0.963, 0.949, 0.991, 0.9998, 0.9988 (ω = 1e3 … 1e7) |
| large p, τ_R = 0 | c₆ u⁶ | 1.0018, 1.0006, 1.0002, 1.00006 (ω = 1e3 … 1e6) |

(The numerical block at the end of 04_expansions.py converted the coefficients to double precision, which produced a
spurious ratio at ω = 1e5; it is superseded by 04b_truncation.py.)

---

## 7. Interface resistance (task 5)  (results/05_resistance.txt, 05_scaling.txt, 05_extra.txt)

θ = (ξ₁, b, τ_R, τ_ℓ, R), Z_s = R + 1/(b_s√p). Column ∂_RZ = Y (Lemma 1). Repeating Proof A with the extra column:
e₁ forces c₂ = 0, h ≡ 0; e₀ = (c₁/ξ₁)A(1 − w²) + c₅B ≡ 0 ⇔ c₁[ρ² − b²(R√p + 1/b_s)²] + c₅b ≡ 0. The odd part in √p is
−2c₁Rb²√p/b_s (SymPy C), so for R > 0: c₁ = c₅ = 0. For R = 0: c₁(ρ² − b²/b_s²) + c₅b ≡ 0 forces c₁ = 0 unless ρ ≡ 1
(τ_R = τ_ℓ), where c₅ = −c₁(1 − b²/b_s²)/b. Hence (DERIVED):

| Stratum | Rank | Kernel |
|---|---|---|
| τ_R ≠ τ_ℓ (any R ≥ 0, any b) | 5 | 0 |
| τ_R = τ_ℓ, R > 0 | 4 | (0,0,1,1,0) |
| τ_R = τ_ℓ, R = 0, b ≠ b_s | 3 | (0,0,1,1,0), (1, 0, 0, 0, −(1 − b²/b_s²)/b) |
| τ_R = τ_ℓ, R = 0, b = b_s | 3 | (0,0,1,1,0), (1,0,0,0,0) |

The second null vector on {τ_R = τ_ℓ, R = 0} (SymPy A5b: U/ξ₁ = ((1 − w²)/b)Y exactly): a thin Fourier slab δe added at
the interface changes Z_s by δe(1/λ − ρc/b_s²), a frequency-independent shift, because p Z_s² = 1/b_s² for a Fourier
substrate; this is a resistance δR = δξ₁(1/b − b/b_s²). It needs both the Fourier film and R = 0 at the base point.

Verification:
- 16 random rational points (4 Cattaneo), 50-digit SVD, symbolic vs numerical derivatives (quadrupole product with
  [[1, R], [0, 1]]) agree to 1.5e-48. Real ranks 2, 4, 5, 5, 5, 5 for N = 1, 2, 3, 4, 5, 8; complex ranks min(N, 5).
  **Minimal N: 3 (real), 5 (complex).** Smallest σ₅/σ₁: real 6.7e-5 (N = 3), complex 6.6e-7 (N = 5).
- Arb at all 16 points: real N = 2 4 × 4 minor, real N = 3 5 × 5 minor, complex N = 5 determinant certified ≠ 0.
- Strata (two random rational points each): singular values and certified minors as in the table; null-vector
  residuals ≤ 8e-51 relative.
- Small-p coefficient map with R: det ∂(z₀…z₄)/∂θ = −(2ξ₁⁵/(3b³b_s⁴))(τ_R − τ_ℓ)·Q, Q = 3b_s²(τ_R − τ_ℓ) + 2Rb³ξ₁ − R²b²b_s².
  On Q = 0 this truncated map is singular, but the function-level rank is still 5 (Arb-certified 5 × 5 minors at two
  points on Q = 0); Q = 0 is not a degeneracy.
- Grid (9 β × 14 T_R × 14 T_ℓ × 7 R̂, 65 frequencies): on Σ with R > 0 max σ₅/σ₁ = 3.6e-16; on Σ with R = 0, β ≠ 1
  max σ₄/σ₁ = 4.0e-16; off Σ min σ₅/σ₁ = 1.7e-11 at (β = 100, T_R = 0, T_ℓ = 1e−3, R̂ = 0), i.e. 1e−3 from the Fourier
  corner, where two singular values vanish at R = 0 (05_scaling: σ₅ ∝ ε near Σ at R = 0; σ₄ ∝ R̂ on Σ for β ≠ 1).
  Off Σ, R → 0 is regular (σ₅/σ₁ = 0.030, 0.055 for β = 0.1, 1; 1.1e-8 for β = 100, Arb-certified rank 5 there).

---

## 8. By-product: injectivity off Σ (claim B2)  (results/06_inversion_series.txt)

Z(iω), ω > 0, determines its Puiseux coefficients at ω → 0⁺, hence z₀ … z₄. Exactly (SymPy I1–I2):
- if z₀ ≠ 0 (⇔ b ≠ b_s): C_f = −b_s z₁/z₀, R_f = z₀ + C_f/b_s², ξ₁ = √(R_fC_f), b = √(C_f/R_f); δ = (z₂ − F₂(ξ₁, b)) b/ξ₁;
  if δ ≠ 0, τ_ℓ = −(z₄ − G₄(ξ₁, b, δ)) b/(ξ₁δ), τ_R = τ_ℓ + δ;
- if z₀ = 0 (⇔ b = b_s): z₂ = ξ₁δ/b_s, z₃ = −ξ₁²δ/b_s, z₄ = −ξ₁δ(3τ_ℓ − 2ξ₁²)/(3b_s); if δ ≠ 0, ξ₁ = −z₃/z₂, δ = b_s z₂/ξ₁,
  τ_ℓ from z₄;
- δ is determined by Z in both cases, so points of Σ and of its complement never share a response.
Hence **F is injective on {τ_R ≠ τ_ℓ}**; the fibre through θ ∈ Σ, b ≠ b_s, is {(ξ₁, b, τ, τ) : τ ≥ 0}; through
θ ∈ Σ, b = b_s, it is {(ξ₁', b_s, τ, τ) : ξ₁' > 0, τ ≥ 0} (Z = 1/(b_s√p)).
The inversion formulas reproduce θ exactly at random rational points (0 failures). End to end, the coefficients
z₋₁ … z₅ were extracted from values of Z alone (trapezoidal Cauchy integrals on |q| = 0.05/ξ₁, 64 nodes, branch
√p = q, 50 digits) and inverted at six random points, one of them with b = b_s: maximum reconstruction errors between
1.5e-32 and 3.3e-46.
Status: DERIVED (constructive); requires the same audit as Theorem 1. It answers B2 in the function-level setting
(B2 of 00-question.md) with b_s known and R = 0.

---

## 9. Independence of the checks; limitations

- Closed form vs model: the closed form was checked against the matrix product (SymPy, Wolfram, mpmath) and against the
  project's code; all rank computations use the closed form or the matrix product, never the brief's expression.
- Derivatives: three independent constructions — SymPy symbolic, Wolfram symbolic (of the matrix product), mpmath
  numerical differentiation (no symbolic input), plus Arb forward-mode differentiation (no symbolic input).
- Rank: floating SVD (50 digits) is evidence; the Arb minors are proofs of the lower bounds at the stated points; the
  upper bounds on Σ are exact identities (SymPy and Wolfram independently). Theorem 1 has two proofs sharing only the
  closed form and Lemma 1: Proof A (transcendence, analytic) and Proof B (exact rational minors of the Puiseux map).
- Correlated elements: SymPy and Wolfram series agree exactly, but both expand the same closed form; the closed form
  itself is independently confirmed (§1). Lemma 2 uses a standard growth argument written here, not machine-checked.
- Not covered: the N = 2 statement for every θ ∉ Σ (numerical only); injectivity with R > 0 (needs more coefficients or
  the large-p data); graded films; any lateral (finite-beam) model.

---

## 10. Files
`fdtr1d.py` (shared model: mpmath, Arb and SymPy versions); `00_closure.py/.wl`; `01_rank.py`, `01_rank_ball.py`,
`01_make_wl.py`, `01_rank.wl`, `01_compare_wl.py`; `02_sigma_symbolic.py`, `02_sigma.wl`, `02_sigma_ball.py`,
`02_sigma_numeric.py`, `02_sigma_corner.py`; `03_gauge.py/.wl`; `04_expansions.py/.wl`, `04_compare_wl.py`,
`04b_truncation.py`;
`05_resistance.py`, `05_scaling.py`, `05_extra.py`; `06_inversion_series.py`, `06b_sigma_minors.py`.
Results: `results/*.txt` (`*_wl.txt` are the Wolfram outputs), points in `results/01_points.json`, grids in
`results/02_grid*.npz`.
