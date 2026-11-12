# Proof audit — Supplemental Material, Sections S2 and S4 (draft of 12 November 2026)

**Date:** 12 November 2026. **Object:** `paper/manuscript/supplement.tex` (compiled `supplement.pdf`; equations
(S1)–(S7), Tables S1–S2, Lemma S1), checked against `sections/identifiability.tex` (Lemma 2, Theorem 3, hypotheses
(H1)–(H5); compiled `main.pdf`), `sections/kinetics.tex` (Eq. (15)), the earlier audit
`checks/proof_audit_sections_II_IV.md` (§1.13 lattice lemma, §4 Lemma 2), `checks/verify_manuscript_algebra.py`
(block "Lemma 2 sharpness") and `checks/audit/audit_lemma2.py`, `checks/audit/audit_ir_lattice.py`. The `.tex` files
were not modified.

**Scope.** Only material that had not been audited before: the proof of Lemma S1 (a)–(c), the two regularity
counterexamples, the mesh laws and the four laws for τ_mem^RTA with their constants, the event-operator bound (S7);
the sufficiency proof of the equality networks, Table S1, the identity min_j K₀/K_RTA = 1/λ_max(G) and the bound (S2);
and the agreement of the main-text statements of Lemma 2 and Theorem 3 (i)–(iii) with what S2 and S4 prove.

**Method.** Every step was re-derived by hand. The fragile steps were rechecked by computations that differ from
those cited in the supplement: other lattices (hexagonal, skew, bcc), angular functions with jumps on a closed null
set of directions, adversarial error terms, an explicit null-space proof of the stoichiometric rank, a bordered-system
solve for K₀, Bose weights at 50 digits, exact real-root isolation of det(λR − C), an exact positive-definiteness
certificate for the strict inequality, implicit differentiation at δ = 0, closed forms for A₁ and A₂, QUADPACK
quadrature with algebraic endpoint weights for the continuum moments, and a smoothed-cutoff (Poisson) evaluation of
the lattice constant S₂. Scripts are in `checks/audit_supp/` (§14). There are 37 checks and all pass. One check of the
mesh-law script first failed because its own criterion ("deviation strictly decreasing from N₁ = 401") met an exact
tie after a pre-asymptotic crossing of 1; the criterion was corrected to the last three meshes and rerun. It did not
involve the manuscript.

**Classification.** VALID · VALID UNDER EXTRA ASSUMPTION (stated) · INCOMPLETE (the conclusion is true but the
argument as written does not reach it) · UNJUSTIFIED · FALSE (counterexample given) · UNKNOWN.

---

## 0. Verdict at a glance

| Item | Verdict | Defect | Repair |
|---|---|---|---|
| Lemma S1, counting constant (6/ρ₀)^d | **VALID** | — | — |
| Lemma S1, Lebesgue criterion and (S6) | **VALID** | "closed" is never used, only "zero measure" | none needed |
| Lemma S1(a), uniform bound and δ-split, (S4) | **VALID** | — | — |
| Lemma S1(b), Cesàro argument and error terms | **VALID** | — | — |
| Lemma S1(c), dominated convergence, positivity | **VALID** | — | — |
| Regularity counterexamples | **VALID** for d ≥ 2 | the first needs d ≥ 2 (for d = 1 every direction is a lattice direction and A = 0) | "(d ≥ 2)" |
| Mesh laws, mα ≥ d | **VALID** | — | — |
| Mesh law, mα < d | **VALID**; one step **INCOMPLETE** as written | the part k < δ of the *optical* branches is not covered by Lemma S1(a) with exponent mα, and (H5) gives Riemann integrability only on k ≥ δ | one clause: Lemma S1(a) with m = 0 bounds it by C″δ^d |
| Four τ_mem^RTA laws with constants | **VALID** | — | optional: α > d gives τ ~ N₁^α S₂/S₁ |
| Event-operator bound (S7) and Theorem 3(ii) | **VALID** under (H4) and uniform (H1) | — | — |
| Theorem 3(iii), main text | **INCOMPLETE statement** | "under the hypotheses of (ii)" omits 2α ≥ d, which (ii) needs for τ_mem → ∞ | "under the hypotheses of (ii) and 2α ≥ d" |
| Table S2 (numbers) | **VALID** (reproduced independently) | the "extrapolated 0.99785" keeps the uncorrected M₂ correction; a three-term fit gives 0.9990 | describe it as such |
| S2, equality criterion and sufficiency (S1) | **VALID** | — | — |
| S2, use of Eq. (15) for 1 ↔ 3 events | **INCOMPLETE** (citation) | Eq. (15) is derived for three-phonon events | one sentence: the same linearization gives ζ = W⁻¹s for 1 ↔ 3 and 2 ↔ 2 events (checked, §6.7) |
| Table S1 (three exact networks) | **VALID** | — | — |
| min_j K₀/K_RTA = 1/λ_max(G) | **VALID**; the step K₀ = uᵀG⁺u is asserted without reason | G⁺ ≠ R^{1/2}C⁺R^{1/2} in general | one line (§7.2) |
| Bound (S2) and the four minima | **VALID** (certified enclosures) | — | optional: the bound holds for every positive weight vector |
| Self-energy counterexample (main text cites S2) | claim **TRUE** (verified here); reference **INCOMPLETE** | S2 marks it PENDING, the main text cites it as given | write the two-line counterexample (§8) or mark the sentence pending |
| Sentence before Lemma 2 ("bounded below by the RTA response") | **imprecise** | the factor 1/ς_max is missing; K₀ < K_RTA occurs (Table S1) | "bounded below by a fixed fraction of the RTA response" |
| Lemma 2 statement vs S2 | **faithful**; label PROVED justified | — | optional scope clause (§10.1) |
| Theorem 3(i) statement vs S4 | **faithful**; label PROVED justified | — | optional remark on S_m (§10.2) |
| Theorem 3(ii) statement vs S4 | **faithful**; label DERIVED UNDER ASSUMPTIONS justified | — | — |

**No mathematical statement of S2 or S4 is false.** The defects are a missing hypothesis in the main-text statement of
Theorem 3(iii), a dangling cross-reference to a counterexample that is still pending in S2, and short gaps in the
writing.

---

## 1. Lemma S1 (lattice sums)

Setting: B invertible, m ≥ 0, α > 0, s := mα; Φ ≥ 0 bounded and continuous outside a closed null set F ⊂ S^{d−1},
A = ∫Φ dΩ > 0; φ(n) = Φ(Bn/|Bn|)|Bn|^{−s}; |ε_N(n)| ≤ ϑ(|Bn|/N), ϑ bounded on (0, k₀), ϑ(t) → 0.

| # | Step | Class. | Reason |
|---|---|---|---|
| 1.1 | ρ₀ = min_{n≠0}\|Bn\| > 0 | VALID | the lattice 𝖡ℤ^d is discrete (B invertible) |
| 1.2 | #{0 < \|Bn\| < R} ≤ (3R/ρ₀)^d for R ≥ ρ₀ | VALID | The open balls of radius ρ₀/2 about the lattice points are disjoint and lie in the ball of radius R + ρ₀/2 ≤ 3R/2, so the count (origin included) is ≤ (2R/ρ₀ + 1)^d ≤ (3R/ρ₀)^d. Checked on square, hexagonal, skew 2D and bcc lattices: the largest ratio is 0.588 |
| 1.3 | Σ_l empty unless 2^{l+1} > ρ₀; Σ_{n∈Σ_l}\|Bn\|^{−s} ≤ C₁2^{l(d−s)}, C₁ = (6/ρ₀)^d (S5) | VALID | Σ_l ⊂ {0 < \|Bn\| < 2^{l+1}} and \|Bn\|^{−s} ≤ 2^{−ls} (s ≥ 0). With R = 2^{l+1} ≥ ρ₀, the count is ≤ (6·2^l/ρ₀)^d. Largest ratio on the four lattices, s ∈ {0, d/2, d, 3d/2}: 0.548 |
| 1.4 | φ̃ is bounded on D = {1 ≤ \|Bu\| < 2} and continuous outside a null set | VALID | u ↦ Bu/\|Bu\| is continuous on u ≠ 0, so φ̃ is continuous wherever Bu/\|Bu\| ∉ F. In polar coordinates v = Bu the cone over F has measure ∫₁²ρ^{d−1}dρ·σ(F) = 0, and B⁻¹ preserves null sets. φ̃ ≤ sup Φ on D |
| 1.5 | Lebesgue criterion; lattice sums are Riemann sums of the cube partition u + [0,h)^d | VALID | φ̃·1_D is bounded with compact support. Its discontinuities lie in the cone over F and in ∂D (two ellipsoids), a null set, so φ̃·1_D is Riemann integrable on a box. For any partition with mesh → 0 the Darboux sums squeeze every tagged sum, the corner-tagged grid sums included. Only "zero measure" is used; the closedness of F is never used (harmless) |
| 1.6 | c* = \|det B\|⁻¹A∫₁²ρ^{d−1−s}dρ > 0 | VALID | v = Bu, du = dv/\|det B\|, polar coordinates; A > 0 |
| 1.7 | φ(n) = 2^{−ls}φ̃(2^{−l}n), n ∈ Σ_l ⇔ 2^{−l}n ∈ D, hence (S6) | VALID | homogeneity of degree −s. Checked on the hexagonal lattice with an angular function that jumps on four directions, s = 1, 2, 3: deviation of a_l/2^{l(d−s)} from c* below 3·10⁻⁵ at l = 11 |
| 1.8 | (a), first claim, C = C₁/(1 − 2^{s−d}) | VALID | The points with \|Bn\| < R lie in shells with 2^l < R. Geometric series with ratio 2^{d−s} > 1, top term 2^{L(d−s)} < R^{d−s}. Then R = δN. The bound holds for every δ > 0 and N, the empty case included. Largest ratio on the hexagonal and skew lattices: 0.10 |
| 1.9 | (a), δ-split: Riemann sum of φ̃ over {δ ≤ \|Bu\| < k₀} with h = 1/N; small ball ≤ Cδ^{d−s} sup Φ; N → ∞, then δ → 0 | VALID | The limsup and liminf are within Cδ^{d−s}sup Φ + \|det B\|⁻¹Aδ^{d−s}/(d−s) of the target for every δ. Checked: N^{s−d}Σφ/target − 1 = −6.9·10⁻⁴ (hexagonal, jumps) and −1.6·10⁻⁴ (skew) |
| 1.10 | (a), error bound in (S4) | VALID | \|N^{s−d}Σε_N\|Bn\|^{−s}\| ≤ ϑ̄N^{s−d}Σ\|Bn\|^{−s}, and the latter converges by 1.9 with Φ ≡ 1, A = \|S^{d−1}\|. Checked with ε_N = +0.2(t/k₀)^{1/2}: the limit equals the integral prediction to 4·10⁻⁵ and lies below the bound (S4). (S4) is not used for the mesh laws, which need only 1.8 |
| 1.11 | (b), Σ_{l<L_N}a_l = c*L_N + o(L_N) | VALID | Cesàro: a_l → c* and the shells with l < l_min are empty. "Only finitely many shells with l < 0 are nonempty" is a roundabout way of saying the sum starts at l_min |
| 1.12 | (b), top shell ≤ C₁ sup Φ; L_N log 2 = log N + O(1) | VALID | 1.3 with s = d |
| 1.13 | (b), error terms o(log N) | VALID | Below δN: at most ε′C₁ per shell, with at most log₂N + log₂k₀ − l_min + 1 nonempty shells, uniformly in δ < k₀. Between δN and k₀N: at most log₂(k₀/δ) + 2 shells, each ≤ ϑ̄C₁. Checked adversarially with ϑ(t) = 0.5/(1 + \|log t\|): the error sum grows like log log N, and its ratio to log N falls (0.0216 → 0.0158 for N = 2⁶ … 2¹⁴). The slope of S_N in log N equals \|det B\|⁻¹A to 2·10⁻⁵ |
| 1.14 | (c), S_Φ converges absolutely | VALID | Σ_l C₁2^{l(d−s)} is geometric with ratio < 1, starting at l_min |
| 1.15 | (c), error terms → 0 by dominated convergence | VALID | Indicator form Σ_{n≠0}1_{\|Bn\|<k₀N}ε_N(n)\|Bn\|^{−s}. Dominant ϑ̄\|Bn\|^{−s} is summable. For each fixed n, ϑ(\|Bn\|/N) → 0. Checked with ϑ(t) = t^{0.3}: decay like N^{−0.3} |
| 1.16 | (c), S_Φ > 0 | VALID | σ({Φ > 0}) > 0 since A > 0. Removing the null F leaves a point ω where Φ is continuous with Φ(ω) > 0, so Φ > 0 on an open U ∋ ω. Q^d∖{0} gives the same directions as Z^d∖{0}, and they are dense. ω ↦ Bω/\|Bω\| is a homeomorphism of S^{d−1} with inverse ω ↦ B⁻¹ω/\|B⁻¹ω\|. Checked with Φ = 1 on the arc 0.10 < θ < 0.12 only: S_N converges to 1.2560·10⁻⁷ > 0, and the first contributing point is n = (9, 1) |
| 1.17 | d = 3 | VALID | bcc reciprocal lattice, jump on a plane of directions: log slope 1.0002 × \|det B\|⁻¹A at s = d, convergence at s = 6 |

**Lemma S1: VALID.** Recommended (cosmetic): drop "closed" or say it is not needed; in (b) say "the shells with
l < l_min are empty".

## 2. Regularity counterexamples

| # | Step | Class. | Reason |
|---|---|---|---|
| 2.1 | Φ = 0 on the lattice directions, 1 elsewhere: A = \|S^{d−1}\| > 0, φ ≡ 0 | VALID **for d ≥ 2** | The set of lattice directions is countable, hence null on S^{d−1} only for d ≥ 2. For d = 1, S⁰ = {±1} consists of lattice directions and A = 0. Φ is discontinuous everywhere (both level sets are dense), so it violates the hypothesis of Lemma S1 |
| 2.2 | J = J₀(k̂) with J₀² = Φ, r = k^α near Γ satisfy (H1)–(H3); the mesh laws fail | VALID | (H1): r₀ = 1. (H2): J − J₀ = 0 near Γ, J bounded, A_m = \|S^{d−1}\|. (H3) is arranged away from Γ. Every mesh point near Γ has a lattice direction, so J = 0 there: S_m = 0 contradicts S_m ∈ (0, ∞) for mα > d, the log law fails at mα = d, and the limit misses the ball part for mα < d |
| 2.3 | J = 0 at the mesh points with k ≥ k₀ of all meshes | VALID | The countable union of meshes is null, so (H1)–(H3) and every 𝓜_m are unchanged, while the part k ≥ k₀ of every mesh sum is removed. The mα < d law fails whenever that part of 𝓜_m is positive |
| 2.4 | "(H5) excludes both constructions" | VALID | 2.2 violates both parts of (H5): Φ_m is nowhere continuous, and J²r⁻ᵐ is nowhere continuous on δ ≤ k < k₀. 2.3 violates the second part: wherever the original integrand is continuous and positive, the modified one is discontinuous on a set of positive measure |
| 2.5 | (H5) "admits the discontinuities of frequency-sorted acoustic branches … closed set of directions of zero measure" | PHYSICAL INTERPRETATION; consistent | At conical acoustic axes the one-sided group velocities of a sorted branch depend on the direction of approach. Isolated directions and crossing curves are closed and null. Scope remark on S_m in §10.2 |

## 3. Mesh laws (Theorem 3(i)) and the four laws for τ_mem^RTA

| # | Step | Class. | Reason |
|---|---|---|---|
| 3.1 | Σ_ac J²r⁻ᵐ = [Φ_m + ε_ac]k^{−mα}, sup_k̂\|ε_ac\| → 0, additive | VALID | As in the earlier audit (§1.2): with the o(1) of (H1) at most 1/2, ε_ac is bounded by n_ac sup J²(2/r₋)^m + sup Φ_m. Then ϑ(t) := sup_{\|k\|=t}\|ε_ac\| is admissible in Lemma S1 |
| 3.2 | The acoustic part k < k₀ equals N₁^{mα−d}V_cell⁻¹S_{N₁} | VALID | V = N₁^dV_cell. The ball k < k₀ lies inside the zone with k₀ < ρ₀/2, so each mesh point in it is its own representative, with no double counting |
| 3.3 | Optical branches and k ≥ k₀ contribute O(1) | VALID | ≤ n_b sup J² c⁻ᵐ/V_cell by (H3) |
| 3.4 | mα > d: 𝓜_m(N₁) = N₁^{mα−d}S_m/V_cell[1 + o(1)] | VALID | Lemma S1(c) and S_m > 0, plus the O(1) of 3.3 |
| 3.5 | mα = d: (2π)^{−d}A_m log N₁ + o(log N₁) | VALID | Lemma S1(b) and \|det B\|V_cell = (2π)^d |
| 3.6 | mα < d: part k < δ ≤ C′δ^{d−mα} uniformly; remaining Riemann sums converge; limsup ≤ 2C′δ^{d−mα} | VALID; **INCOMPLETE as written** | Lemma S1(a) with exponent mα covers the acoustic branches only. On the *optical* branches near Γ the integrand is bounded, but (H5) gives Riemann integrability only on k ≥ δ, so their part k < δ must also be bounded. Lemma S1(a) with m = 0 (#{k_ν < δ}/V ≤ C(0)δ^d/V_cell) does this, as does the trivial bound for the integral. Repair: "for all branches", or one clause citing Lemma S1(a) with m = 0. The Riemann-sum step on k ≥ δ is standard, with boundary points of the zone O(N₁^{d−1}) in number |
| 3.7 | τ laws: 2α < d → 𝓜₂/𝓜₁; 2α = d: (2π)^{−d}A₂ log N₁/𝓜₁; d/2 < α < d: N₁^{2α−d}S₂/(V_cell𝓜₁); α = d: (2π)^d S₂N₁^d/(V_cell A₁ log N₁) | VALID | Ratios of 3.4–3.6 with m = 1, 2. 𝓜₁ > 0 by (H2). Dimensions check: each constant is a time. Not stated but immediate: α > d gives τ = N₁^α S₂/S₁[1 + o(1)] |

**Independent reproduction (paper's test model, Table S2).**
* A₁ = π/√1.6 and A₂ = 1.3π/1.6^{3/2} in closed form.
* Continuum moments by QUADPACK with the weight k^{1−mα}: 𝓜₂/𝓜₁ = 0.6291353667 (α = 0.6) and (2π)⁻²A₂/𝓜₁ = 0.3885388430 (α = 1).
* Smoothed-cutoff evaluation of Σ_{n≠0}Φ₂(n̂)\|n\|^{−3} = 3.4851299133, stable to 10⁻¹³ in R. It agrees with the sharp-cutoff value 3.48512991.
* Quadrant-folded mesh sums reproduce every entry of Table S2:
  * 2.606·10⁻³ and exponent 0.792;
  * 0.996570, increasing;
  * 2.505·10⁻³, exponent 0.5000, ratio 1.02571, Richardson value 0.99785;
  * 1.000867, 1.001582 and 1.0572.

**Presentation remark (Table S2, α = 1.5).** The value 0.99785 comes from a Richardson step that removes only the
N₁^{−1/2} correction of 𝓜₁. 𝓜₂(N₁)/(N₁S₂) − 1 is still +2.5·10⁻³ at N₁ = 3201 and decays more slowly than N₁^{−1}. A
three-term fit 1 + aN₁^{−1/2} + bN₁^{−1}log N₁ on N₁ = 801, 1601, 3201 gives 0.9990. The proved limit is 1. Suggest
"0.998 (Richardson in N₁^{−1/2}; the remainder is the 𝓜₂ correction, 2.5·10⁻³ at N₁ = 3201)". Also: "exact lattice
sums" in the caption means direct double-precision sums. "Direct lattice sums" is more accurate. The sentence "𝓜₁(N₁),
which converges like N₁^{−(d−α)} for α < d" is stated generally but observed at one α. It is the generic rate (its
coefficient is a lattice-zeta constant of Φ₁ that could vanish), so it should be scoped to the test, as in the main
text.

**Independent extension (hexagonal lattice, jumps on lattice directions).** J² and r₀ jump on θ ∈ {0, π/3, π, 4π/3}.
Two of these directions are lattice directions. The values there are fixed by exact integer tests. Results:
* 𝓜₁, 𝓜₂ converge for α = 0.6, with errors 1.6·10⁻⁴ and 2.1·10⁻³ at N₁ = 3201.
* The log slope of 𝓜₂ at α = 1 is 1.0012 × (2π)⁻²A₂.
* 𝓜₂V_cell/(N₁S₂) = 1.0015 at α = 1.5.
* Changing r₀ only on the two lattice directions θ = 0, π changes S₂ by −5.9 % (§10.2).

## 4. Event-operator bound (S7) and Theorem 3(ii)–(iii)

| # | Step | Class. | Reason |
|---|---|---|---|
| 4.1 | Uniform (H1) ⇒ r_ν ≤ r₀k_ν^α[1 + δ₁(k_ν)], δ₁ independent of N₁ | VALID | This is the definition in Theorem 3(ii). Take δ₁ ≥ 0 and the max over branches, and k₀ also below the radius of (H4) (min of the two radii) |
| 4.2 | (S7) | VALID | r_ν^{−2} ≥ r₀^{−2}k^{−2α}(1 + δ₁)^{−2}. Summing over the acoustic branches at each k gives V⁻¹k^{−2α}[Φ₂ + ε′], with ε′ = Σ_ac r₀^{−2}[J²(1 + δ₁)^{−2} − J₀²] → 0 uniformly (J → J₀ uniformly, J bounded, r₀ ≥ r₋). ε′ is bounded once k₀ is taken with δ₁ ≤ 1/2 |
| 4.3 | Right side ~ (2π)^{−d}A₂ log N₁ (2α = d), N₁^{2α−d}S₂/V_cell (2α > d) | VALID | Lemma S1(b), (c) with Φ = Φ₂ (bounded by (H1)–(H2), continuous outside a closed null set by hypothesis), ε_N(n) = ε′(Bn/N₁), A₂ > 0 |
| 4.4 | "Only the continuity of Φ₂ … enters" | VALID | No Riemann integrability away from Γ is used |
| 4.5 | (H4) and Prop. 2 ⇒ (ii) | VALID | N ≥ Σ_{ac,k_ν<k₀}x_ν² ≥ c²Σj_ν²r_ν^{−2}, and \|j\|² ≤ n_b sup J²/V_cell. Then τ ≥ N^{1/2}/\|j\| ≳ N₁^{α−d/2} for 2α > d and ≳ (log N₁)^{1/2} for 2α = d |
| 4.6 | Theorem 3(iii): "Under the hypotheses of (ii), M²_{N₁}K_{0,N₁} → ∞" | **INCOMPLETE statement** | Theorem 2 gives τ ≤ M²K₀. The divergence of τ is established in (ii) only *if 2α ≥ d*, which is not among "the hypotheses of (ii)" (the "Suppose that …" list). For 2α < d the statement is unproved. It also uses implicitly that 𝒞_{N₁} = 𝒞(g_{N₁}) ∈ 𝔈(Z_{N₁}), which "the event family on the mesh" implies. Repair: "Under the hypotheses of (ii) and 2α ≥ d, …" |

## 5. Status labels of Theorem 3 (main text) against S4

* (i) PROVED. The integral criterion (main text) uses (H1)–(H3). The mesh laws (S4) use (H5) in addition. The
  statement in the main text matches S4 exactly, including S_m ∈ (0, ∞) and the four τ laws. The label is justified
  once 3.6 is completed with one clause.
* (ii) DERIVED UNDER ASSUMPTIONS, "namely (H4) and the uniform form of (H1)". These are exactly the unproved properties
  used in §4. The hypotheses (H2) and the continuity of Φ₂ are proved inputs. The label is faithful. The defect of the
  earlier audit (§1.14 there) is repaired.
* (iii) same label. It is faithful apart from the missing 2α ≥ d (4.6).

---

## 6. S2: equality networks

| # | Step | Class. | Reason |
|---|---|---|---|
| 6.1 | Necessity: equality ⇒ j = 𝒞y with 𝒞y = ς diag(r)y, every active event of size ς, ζ_{γν}y_ν constant on supports | VALID | Lemma 2. Event term g_γ[ς Σζ²y² − (ζᵀy)²] ≥ 0 by Cauchy–Schwarz on the support. It vanishes iff ζ_{γν}y_ν is constant on the support and (s_γ − ς)Σ_{supp}ζ²y² = 0. The connectivity argument needs e to have no zero entry, which (K1) gives |
| 6.2 | Sufficiency (S1): ζ_γᵀy = ςc_γ, ζ²_{γν}y_ν = ζ_{γν}c_γ ⇒ 𝒞y = ς diag(r)y | VALID | Event by event: ζ_γ(ζ_γᵀy) = ς diag(ζ_γ∘ζ_γ)y. Off the support both sides vanish. Inactive events do not enter. Checked symbolically, event by event, for symbolic weights (rates factor out) |
| 6.3 | j ∈ ℋ∖{0}, K₀ = yᵀ𝒞y, K_RTA = ς²yᵀRy, ratio 1/ς | VALID | 𝒞𝒞⁺𝒞 = 𝒞. j = ςRy ≠ 0 because r_ν > 0. Admissibility makes j ∈ range 𝒞 = ℋ |
| 6.4 | ζ_γ = W⁻¹s_γ, y = Wψ ⇒ ζ_{γν}y_ν = s_{γν}ψ_ν; "all positive weights and rates" | VALID | Admissibility for every positive rate vector follows from rank n − 1 and W invertible |
| 6.5 | Table S1 ranks n − 1 | VALID | Proved here by hand through the null space of the map x ↦ (s_γᵀx)_γ, which is span{ε} in each case. Chain: x_{ν+1} = 2x_ν. ς = 3: the four events from mode 1 with mode 4 give x₅ = … = x₈ = x₁ − x₄, the event 1 → 5 + 6 gives x₁ = 2x₄, and so x ∝ (2,2,2,1,…,1). ς = 4: likewise x₄ = … = x₉, x₁ = x₂ = x₃ = 3x₄. Confirmed by exact null-space computation |
| 6.6 | ψ choices: chain ψ_{ν+1} = −2ψ_ν; parents +1, daughters −1 | VALID | Each event of the ς = 3, 4 networks has distinct daughters, so s_{γν}ψ_ν = −1 on every support. Checked |
| 6.7 | "[Eq. (15)]" for the 1 ↔ 3 events of the ς = 4 network | **INCOMPLETE** (citation) | Eq. (15) is derived in the main text for p → a + b. The same linearization of G[f_p∏(1+f_a) − (1+f_p)∏f_a] under detailed balance gives −P sᵀW^{−2}δf with s = −δ_p + δ_a + δ_b + δ_c (and likewise for 2 ↔ 2), hence ζ = W⁻¹s. This was rederived here. Repair: one sentence |
| 6.8 | Pure 2 ↔ 2 networks: phonon number W𝟏 is invariant, ∉ span{e} unless all energies are equal | VALID | ζ_γᵀW𝟏 = s_γᵀ𝟏 = 0, and W𝟏 ∥ Wε iff ε ∥ 𝟏. The "unless" is needed: for 4 modes of equal energy, the three 2 ↔ 2 events have rank 3 = n − 1 |
| 6.9 | Verification paragraph | faithful | The block "Lemma 2 sharpness" of `verify_manuscript_algebra.py` does what is stated: energy conservation, rank, symbolic identity, rational ratio with K₀ = jᵀ(𝒞 + eeᵀ/\|e\|²)⁻¹j |

**Independent checks** (`audit_s2_networks.py`):
* K₀/K_RTA = 1/ς exactly from the bordered system [[𝒞, e], [eᵀ, 0]] for four random rational weight/rate sets per
  network.
* With Bose weights at T = 0.7, 3 and 20 (energy units), the ratio is 1/ς to 10⁻⁵¹.
* det(ςR − 𝒞) = 0 and ςR − 𝒞 ⪰ 0 exactly, so λ_max(G) = ς.
* The equality criterion ("a y with constant ζ_{γν}y_ν on all supports exists" ⇔ det(ςR − 𝒞) = 0) holds on 99 random
  admissible uniform networks, 50 of them with equality.

## 7. S2: networks with events of different sizes

| # | Step | Class. | Reason |
|---|---|---|---|
| 7.1 | u = R^{−1/2}j ∈ range G ⇔ j ∈ ℋ | VALID | range G = R^{−1/2} range 𝒞 |
| 7.2 | K₀ = uᵀG⁺u, K_RTA = \|u\|² | VALID; **reason missing** | G⁺ ≠ R^{1/2}𝒞⁺R^{1/2} in general (a congruence does not commute with the Moore–Penrose inverse), so the step needs a line. If Gz = u, then 𝒞(R^{−1/2}z) = R^{1/2}u = j. Hence K₀ = jᵀR^{−1/2}z = uᵀG⁺u (any solution gives the same jᵀx for j ∈ range 𝒞). Equivalently, the variational principle with y = R^{−1/2}z |
| 7.3 | min_j K₀/K_RTA = 1/λ_max(G), attained | VALID | Rayleigh quotient of G⁺ on range G, attained at the top eigenvector. Checked: K₀/K_RTA at j* = 𝒞y* (y* the top generalized eigenvector) equals 1/λ_max to 10⁻²⁵. 300 random j never undercut it |
| 7.4 | λ_max(G) ≥ yᵀ𝒞y/yᵀRy for every y ≠ 0 | VALID | z = R^{1/2}y; y need not lie in ℋ |
| 7.5 | Mixed family: ς_max = 4, rank 5, admissible | VALID | Null space span{ε}. The six three-phonon events alone already have rank 5 |
| 7.6 | yᵀ𝒞y = 16 + 10δ, yᵀRy = 4 + 12δ | VALID | Rederived by hand and symbolically. **With y = Wψ both forms are independent of the weights**, because ζ_γᵀy = s_γᵀψ and r_νw_ν² = Σ_γg_γs²_{γν}. The bound (S2) therefore holds for every positive weight vector, Bose weights at any T included. The text states it for unit weights only |
| 7.7 | Strictness 1 < min 4K₀/K_RTA | VALID | Lemma 2. Exact certificate: every leading principal minor of 4R − 𝒞 is a polynomial in δ with nonnegative coefficients, for example det(4R − 𝒞) = 16δ³(25519δ³ + 36634δ² + 17285δ + 2678). So 4R − 𝒞 ≻ 0 for every δ > 0 |
| 7.8 | (S2) upper bound and limit 1 | VALID | From 7.4 and 7.6. Implicit differentiation of the reduced determinant det(λR − 𝒞)/δ² at the simple root (4, 0) gives 4/λ_max = 1 + (1339/568)δ + O(δ²) = 1 + 2.3574δ + O(δ²), against the slope 38/16 = 2.375 of the bound. The bound is sharp to 0.7 % at first order |
| 7.9 | The minima 1.2198, 1.0234, 1.00236, 1.000024 | VALID | Exact real-root isolation of det(λR − 𝒞) (width 10⁻³⁰) gives enclosures 1.2198290594650, 1.0234090535814, 1.0023557393965 and 1.0000235737781. They round to the stated values, lie in (1, bound] and decrease |

## 8. Self-energy linewidths (main text cites S2; S2 is pending)

The main text after Lemma 2 says: "the inequality can fail for it (an exact counterexample is given in the
Supplemental Material, Sec. S2)". The S2 paragraph is `\PENDING`. **The claim is true.** Take an equality network and
its equality current j = 𝒞y. If r̃ ≤ r with r̃_ν < r_ν on some mode with j_ν ≠ 0, then K̃_RTA = Σj²/r̃ > K_RTA, so
K₀ = K_RTA/ς < K̃_RTA/ς. For the chain of Table S1, counting the self-coupling 4g_γ/w²_ν of the repeated daughter with
half weight lowers r_ν on modes 1–4. Every mode carries current (j_ν = 2r_νw_νψ_ν ≠ 0). Exact arithmetic with three
random weight and rate sets gives K₀/(K̃_RTA/2) = 0.722, 0.686 and 0.688. **Status of the reference: INCOMPLETE.** The
paragraph can be written from the lines above.

## 9. Main-text statements against S2

* **Lemma 2.** "Equality … iff 𝒞(g)y = ς diag(r)y …; this requires every active event to touch exactly ς_max modes.
  The constant is attained for ς_max = 2, 3 and 4. When some active event touches fewer than ς_max modes … the
  inequality is strict, but the constant cannot be improved." Each sentence is proved: §6 (Table S1 with the
  sufficiency proof) and §7 (family (S2)). The label PROVED is justified.
* **Sentence before Lemma 2:** "the dc response is bounded below by the RTA response built from the same diagonal".
  This is imprecise: the bound is K_RTA/ς_max, and Table S1 has K₀ = K_RTA/ς_max < K_RTA. Repair: "by the RTA response
  built from the same diagonal, up to the factor 1/ς_max".
* The proof of Lemma 2 cites S2 for "Exact networks attain … and a family of mixed networks approaches it". S2 delivers
  both.

## 10. Scope remarks (not defects)

**10.1 Which currents.** The attainment in Lemma 2 is shown for currents j = 𝒞Wψ ∈ ℋ of abstract stoichiometric
networks. Physical currents are odd under time reversal. With the parent/daughter pattern ψ = ±1, the equality vector
is even under any time reversal that maps parents to parents. Whether time-reversal-odd currents of a
time-reversal-symmetric admissible network attain 1/ς_max is not addressed. It is UNKNOWN, analogous to the open point
of Remark 5. Optional clause: "for currents in ℋ".

**10.2 The constant S_m depends on values on a null set.** For mα > d, S_m = ΣΦ_m(Bn/\|Bn\|)\|Bn\|^{−mα} is a sum over
lattice directions, a null set. Its dominant terms come from the shortest lattice vectors, which in symmetric crystals
lie along the high-symmetry directions where acoustic branches are degenerate. These are the points of the
discontinuity set that (H5) admits. The mesh law remains correct with the values of Φ_m that the mesh data take there
(the uniform (H1)–(H2) fix them as limits along each lattice ray). The hexagonal test shows that changing r₀ only on
θ = 0, π changes S₂ by −5.9 %, and that the mesh sums follow the values actually used. Physically, the leading
constant for mα > d depends on how degenerate modes are treated on those rays (basis or averaging convention). The
growth exponent does not. Optional sentence after (H5).

**10.3 Bose weights.** §6.4 and §7.6 hold for every positive weight vector, so all S2 constructions are valid with
the physical weights w_ν = [f⁰(1 + f⁰)]^{1/2} at any temperature.

---

## 11. Earliest failures in reading order

1. Main text, sentence before Lemma 2: missing factor 1/ς_max (wording, §9).
2. Main text, text after Lemma 2: cross-reference to a pending S2 paragraph (claim true, §8).
3. Main text, Theorem 3(iii): missing hypothesis 2α ≥ d (§4.6) — **the earliest gap in a formal statement**.
4. S2: Eq. (15) cited for 1 ↔ 3 events (§6.7); K₀ = uᵀG⁺u unexplained (§7.2).
5. S4: optical part k < δ in the mα < d law (§3.6); d ≥ 2 in the first counterexample (§2.1); presentation of the
   α = 1.5 extrapolation in Table S2 (§3).

None of these propagates into another proof, and no conclusion changes.

## 12. Minimal repairs (consolidated)

* Theorem 3(iii): "Under the hypotheses of (ii) and 2α ≥ d, M²_{N₁}K_{0,N₁} → ∞ …".
* Before Lemma 2: "… is bounded below by the RTA response built from the same diagonal, up to the factor 1/ς_max".
* S2: write the self-energy paragraph (§8), or mark the main-text parenthesis as pending.
* S2, after "[Eq. (15)]": "(the same linearization gives ζ_γ = W⁻¹s_γ for 1 ↔ 3 and 2 ↔ 2 events)".
* S2, before "gives K₀ = uᵀG⁺u": "if Gz = u then 𝒞(diag(r)^{−1/2}z) = j, so that".
* S2, mixed family: optionally "for all positive weights (take y = Wψ)".
* S4, mα < d bullet: "the part k < δ of 𝓜_m(N₁), on all branches, is at most C′δ^{d−mα} (Lemma S1(a), with exponent
  mα on the acoustic and 0 on the optical branches)".
* S4, first counterexample: "(d ≥ 2)".
* Table S2, α = 1.5: describe the extrapolation (§3); caption "direct lattice sums"; scope "converges like
  N₁^{−(d−α)}" to the test.
* Optional: "closed" in Lemma S1 is unused; add α > d: τ_mem^RTA = N₁^αS₂/S₁[1 + o(1)].

## 13. Candidates for machine checking

* **Sufficiency identity (S1) together with the equality criterion of Lemma 2**, event by event: if \|supp ζ\| = ς and
  ζ_νy_ν = c on supp ζ, then ζ(ζᵀy) = ς diag(ζ∘ζ)y, and conversely yᵀ(ς diag(ζ∘ζ) − ζζᵀ)y = 0 forces this. It is short,
  central (it carries "attained" and "strict" in Lemma 2) and finite-dimensional. It is recommended for Lean 4 + Mathlib
  (`Finset.sum`, `inner_mul_le_norm_mul_norm`), together with the inequality 𝒞 ⪯ ς_max diag 𝒞 already recommended in
  the earlier audit. The check was not run here: the memory limit of this audit (≤ 2 GB beside a running 15 GB job) is
  below what a Mathlib import typically needs.
* Lemma S1 is now written out and audited, and is costly to formalize. A written proof suffices.

## 14. Check files (`paper/manuscript/checks/audit_supp/`, research Python, single thread, < 400 MB, about 1 min in total)

| file | content | result |
|---|---|---|
| `audit_s2_networks.py/.out.txt` | Table S1: ranks through null spaces, event-wise symbolic identity, bordered-system K₀, Bose weights at 50 digits, det(ςR − 𝒞) = 0 with ςR − 𝒞 ⪰ 0. Equality criterion on 99 random networks. Mixed family: weight independence, exact positive-definiteness certificate, certified minima, slope 1339/568. Minimizer of K₀/K_RTA. Self-energy counterexample | 18/18 |
| `audit_s4_lattice_lemma.py/.out.txt` | Lemma S1 on square, hexagonal, skew and bcc lattices with jump functions and adversarial errors: (S5) constants, (S6), (a) with (S4), (b) with a log log N error, (c) with an arc-supported Φ, d = 3 | 6/6 |
| `audit_s4_mesh_laws.py/.out.txt/.json` | Table S2 reproduced by closed forms, QUADPACK and smoothed-cutoff S₂. Hexagonal mesh laws with jumps on lattice directions; dependence of S₂ on null-set values | 13/13 |

Run each with `python -B <file> > <file>.out.txt` from `checks/audit_supp/`.

## 15. Overall verdict

* **Lemma S1 is a valid proof.** The counting constant (6/ρ₀)^d, the Lebesgue-criterion step, the δ-split of (a), the
  Cesàro argument and the error bookkeeping of (b), and the dominated convergence and positivity argument of (c) are
  all correct. The constants are confirmed on lattices and angular functions that the paper's own test does not
  cover, discontinuous ones included. The two regularity counterexamples are valid (the first for d ≥ 2) and show that
  (H5) is not removable.
* **The mesh laws and the four τ_mem^RTA laws are proved under (H1)–(H3) and (H5)**, with one clause to add for the
  optical branches near Γ in the mα < d case. The label PROVED of Theorem 3(i) is justified.
* **The event-operator bound (S7) is valid.** Theorem 3(ii) is a valid derivation under (H4) and the uniform (H1), and
  its label names exactly these. Theorem 3(iii) needs "2α ≥ d" added to its statement.
* **S2 is valid.**
  * The sufficiency proof of the equality networks is correct.
  * The stoichiometric reduction ζ_{γν}y_ν = s_{γν}ψ_ν makes Table S1 independent of weights and rates.
  * The three networks are admissible (ranks proved by hand) and attain 1/ς_max exactly for ς_max = 2, 3, 4.
  * min_j K₀/K_RTA = 1/λ_max(G) holds; the step K₀ = uᵀG⁺u deserves its one-line reason.
  * The bound 1 < min 4K₀/K_RTA ≤ (16 + 48δ)/(16 + 10δ) is certified exactly. It holds for all positive weights and is
    sharp to first order within 0.7 %.
* The main-text statements of Lemma 2 and Theorem 3(i)–(ii) match what S2 and S4 prove, with honest status labels.
  The open items are Theorem 3(iii) (missing 2α ≥ d), the wording before Lemma 2, and the self-energy counterexample.
  The main text cites that counterexample as given, but S2 still marks it pending; it is true and two lines long (§8).
