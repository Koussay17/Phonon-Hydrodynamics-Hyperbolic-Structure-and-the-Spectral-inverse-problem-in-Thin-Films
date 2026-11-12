import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

/-!
# Theorem 1, Steps 4–6: from the film invariants to the parameters

Campaign `20261111-fdtr-gauge-memory`, `proofs/B-analytic.md`, §4 (Theorem 1, Proof A).

The memory ratio of the film law is `g(p) = (1 + τ_R p)/(1 + τ_ℓ p)`.  This file checks the
finite-dimensional part of the global identifiability argument:

* `crossmul_iff_coeffs`, `step5`, `step5_iff`: **Step 5**, the polynomial identity
  `(1 + a p)(1 + d p) = (1 + c p)(1 + b p)` for all `p`, where `(a, b) = (τ_R, τ_ℓ)` and
  `(c, d) = (τ_R', τ_ℓ')`.  It holds iff `(c, d) = (a, b)` or `a = b ∧ c = d`.
* `step5_g`, `step5_two_points`: the same conclusion from `g = g'` on `p > 0` (or at two distinct
  positive points) for nonnegative times.
* `step4`: **Step 4** (algebraic part): equal `z = ξ₁² p g` and `B² = b² p / g` on `p > 0` give
  `ξ₁ = ξ₁'`, `b = b'`, `g = g'`.
* `mem_fibre_iff_of_not_D2`: Steps 4–5 combined, as a statement about the fibre `𝔉(θ)`.
* `thm1_of_steps`: the case assembly of Theorem 1, (a) ⇔ (c), from the analytic steps taken as
  hypotheses.
* `hasDerivAt_g_tR`, `hasDerivAt_g_tl`, `g_tangent_resonance`: `∂_{τ_R} g + ∂_{τ_ℓ} g = 0` on
  `τ_R = τ_ℓ` (the kernel direction `(0,0,1,1)` of Theorem 2).

The analytic steps (identity principle, two-load argument, local inversion of `cosh √z`) are not
formalised here; they enter `thm1_of_steps` only as hypotheses.
-/

namespace FdtrGauge

/-! ## Step 5: the rational identity -/

/-- The cross-multiplied identity `(1 + a p)(1 + d p) = (1 + c p)(1 + b p)` for all real `p`
is equivalent to the two coefficient identities `a - b = c - d` and `a d = c b`. -/
theorem crossmul_iff_coeffs (a b c d : ℝ) :
    (∀ p : ℝ, (1 + a * p) * (1 + d * p) = (1 + c * p) * (1 + b * p)) ↔
      (a - b = c - d ∧ a * d = c * b) := by
  constructor
  · intro h
    have h1 := h 1
    have h2 := h (-1)
    exact ⟨by linear_combination (h1 - h2) / 2, by linear_combination (h1 + h2) / 2⟩
  · rintro ⟨h1, h2⟩ p
    linear_combination p * h1 + p ^ 2 * h2

/-- **Theorem 1, Step 5.** If `(1 + a p)(1 + d p) = (1 + c p)(1 + b p)` for all `p`, then either
`a ≠ b` and `(c, d) = (a, b)`, or `a = b` and `c = d`. -/
theorem step5 {a b c d : ℝ}
    (h : ∀ p : ℝ, (1 + a * p) * (1 + d * p) = (1 + c * p) * (1 + b * p)) :
    (a ≠ b → c = a ∧ d = b) ∧ (a = b → c = d) := by
  obtain ⟨h1, h2⟩ := (crossmul_iff_coeffs a b c d).1 h
  refine ⟨fun hab => ?_, fun hab => by linarith⟩
  have hδ : (a - b) * (d - b) = 0 := by linear_combination h2 - b * h1
  have hd : d = b := by
    rcases mul_eq_zero.1 hδ with h | h
    · exact absurd (sub_eq_zero.1 h) hab
    · exact sub_eq_zero.1 h
  exact ⟨by linarith, hd⟩

/-- **Theorem 1, Step 5, fibre form.** The cross-multiplied identity holds for all `p` iff the
two pairs coincide or both lie on the resonance. -/
theorem step5_iff (a b c d : ℝ) :
    (∀ p : ℝ, (1 + a * p) * (1 + d * p) = (1 + c * p) * (1 + b * p)) ↔
      ((c = a ∧ d = b) ∨ (a = b ∧ c = d)) := by
  constructor
  · intro h
    by_cases hab : a = b
    · exact Or.inr ⟨hab, (step5 h).2 hab⟩
    · exact Or.inl ((step5 h).1 hab)
  · rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩) p <;> ring

/-- The memory ratio `g(p) = (1 + τ_R p)/(1 + τ_ℓ p)` of the film law. -/
noncomputable def g (tR tl p : ℝ) : ℝ := (1 + tR * p) / (1 + tl * p)

lemma one_add_mul_pos {t p : ℝ} (ht : 0 ≤ t) (hp : 0 ≤ p) : 0 < 1 + t * p := by positivity

lemma g_pos {tR tl p : ℝ} (hR : 0 ≤ tR) (hl : 0 ≤ tl) (hp : 0 ≤ p) : 0 < g tR tl p :=
  div_pos (one_add_mul_pos hR hp) (one_add_mul_pos hl hp)

/-- On the resonance `τ_R = τ_ℓ = t ≥ 0`, `g ≡ 1` on `p ≥ 0`. -/
lemma g_diag {t p : ℝ} (ht : 0 ≤ t) (hp : 0 ≤ p) : g t t p = 1 :=
  div_self (one_add_mul_pos ht hp).ne'

/-- Equality of `g` at a point `p ≥ 0` is the cross-multiplied identity at `p`. -/
lemma g_eq_iff_crossmul {tR tl tR' tl' p : ℝ} (hl : 0 ≤ tl) (hl' : 0 ≤ tl') (hp : 0 ≤ p) :
    g tR tl p = g tR' tl' p ↔ (1 + tR * p) * (1 + tl' * p) = (1 + tR' * p) * (1 + tl * p) := by
  unfold g
  rw [div_eq_div_iff (one_add_mul_pos hl hp).ne' (one_add_mul_pos hl' hp).ne']

/-- Two distinct positive points already force the cross-multiplied identity for all `p`
(the difference of the two sides is `p (α + β p)`). -/
theorem crossmul_of_two_points {tR tl tR' tl' p₁ p₂ : ℝ} (hl : 0 ≤ tl) (hl' : 0 ≤ tl')
    (hp₁ : 0 < p₁) (hp₂ : 0 < p₂) (hne : p₁ ≠ p₂)
    (h₁ : g tR tl p₁ = g tR' tl' p₁) (h₂ : g tR tl p₂ = g tR' tl' p₂) :
    ∀ p : ℝ, (1 + tR * p) * (1 + tl' * p) = (1 + tR' * p) * (1 + tl * p) := by
  have e₁ := (g_eq_iff_crossmul hl hl' hp₁.le).1 h₁
  have e₂ := (g_eq_iff_crossmul hl hl' hp₂.le).1 h₂
  -- `e_i : p_i (α + β p_i) = 0` with `α = tR + tl' - tR' - tl`, `β = tR tl' - tR' tl`
  have hα₁ : (tR + tl' - tR' - tl) + (tR * tl' - tR' * tl) * p₁ = 0 := by
    have : p₁ * ((tR + tl' - tR' - tl) + (tR * tl' - tR' * tl) * p₁) = 0 := by
      linear_combination e₁
    exact (mul_eq_zero.1 this).resolve_left hp₁.ne'
  have hα₂ : (tR + tl' - tR' - tl) + (tR * tl' - tR' * tl) * p₂ = 0 := by
    have : p₂ * ((tR + tl' - tR' - tl) + (tR * tl' - tR' * tl) * p₂) = 0 := by
      linear_combination e₂
    exact (mul_eq_zero.1 this).resolve_left hp₂.ne'
  have hβ : tR * tl' - tR' * tl = 0 := by
    have : (tR * tl' - tR' * tl) * (p₁ - p₂) = 0 := by linear_combination hα₁ - hα₂
    exact (mul_eq_zero.1 this).resolve_right (sub_ne_zero.2 hne)
  have hα : tR + tl' - tR' - tl = 0 := by linear_combination hα₁ - p₁ * hβ
  intro p
  linear_combination p * hα + p ^ 2 * hβ

/-- **Theorem 1, Step 5, two-point form.** For nonnegative times, `g = g'` at two distinct
positive points gives the conclusion of Step 5. -/
theorem step5_two_points {tR tl tR' tl' p₁ p₂ : ℝ} (hl : 0 ≤ tl) (hl' : 0 ≤ tl')
    (hp₁ : 0 < p₁) (hp₂ : 0 < p₂) (hne : p₁ ≠ p₂)
    (h₁ : g tR tl p₁ = g tR' tl' p₁) (h₂ : g tR tl p₂ = g tR' tl' p₂) :
    (tR ≠ tl → tR' = tR ∧ tl' = tl) ∧ (tR = tl → tR' = tl') :=
  step5 (crossmul_of_two_points hl hl' hp₁ hp₂ hne h₁ h₂)

/-- **Theorem 1, Step 5, rational-function form.** For nonnegative times, `g = g'` on `p > 0`
gives: off the resonance the times coincide; on the resonance the other point is on the
resonance too. -/
theorem step5_g {tR tl tR' tl' : ℝ} (hl : 0 ≤ tl) (hl' : 0 ≤ tl')
    (h : ∀ p : ℝ, 0 < p → g tR tl p = g tR' tl' p) :
    (tR ≠ tl → tR' = tR ∧ tl' = tl) ∧ (tR = tl → tR' = tl') :=
  step5_two_points hl hl' one_pos two_pos (by norm_num) (h 1 one_pos) (h 2 two_pos)

/-! ## Step 4 (algebraic part) -/

/-- A real quadratic polynomial vanishing at `p = 1, 2, 3` has zero constant term. -/
lemma const_coeff_of_three_points {A B C : ℝ}
    (h1 : A + B * 1 + C * 1 ^ 2 = 0) (h2 : A + B * 2 + C * 2 ^ 2 = 0)
    (h3 : A + B * 3 + C * 3 ^ 2 = 0) : A = 0 := by
  linear_combination 3 * h1 - 3 * h2 + h3

/-- **Theorem 1, Step 4 (algebraic part).** If the invariants `z = ξ₁² p g(p)` and
`B² = b² p / g(p)` of two points of `Θ` agree on `p > 0`, then `ξ₁ = ξ₁'`, `b = b'` and `g = g'`
on `p > 0`.  (The informal proof uses `g(0) = 1`; here the constant term of a quadratic identity
in `p` plays that role.) -/
theorem step4 {xi b tR tl xi' b' tR' tl' : ℝ} (hxi : 0 < xi) (hb : 0 < b) (hxi' : 0 < xi')
    (hb' : 0 < b') (hl : 0 ≤ tl) (hR' : 0 ≤ tR') (hl' : 0 ≤ tl')
    (hz : ∀ p : ℝ, 0 < p → xi ^ 2 * p * g tR tl p = xi' ^ 2 * p * g tR' tl' p)
    (hB : ∀ p : ℝ, 0 < p → b ^ 2 * p / g tR tl p = b' ^ 2 * p / g tR' tl' p) :
    xi = xi' ∧ b = b' ∧ ∀ p : ℝ, 0 < p → g tR tl p = g tR' tl' p := by
  -- (i) on p > 0: ξ² g = ξ'² g'
  have hz' : ∀ p : ℝ, 0 < p → xi ^ 2 * g tR tl p = xi' ^ 2 * g tR' tl' p := by
    intro p hp
    have := hz p hp
    have hp0 : p ≠ 0 := hp.ne'
    apply mul_left_cancel₀ hp0
    linear_combination this
  -- (ii) cross-multiplied form at p > 0
  have key : ∀ p : ℝ, 0 < p →
      xi ^ 2 * ((1 + tR * p) * (1 + tl' * p)) = xi' ^ 2 * ((1 + tR' * p) * (1 + tl * p)) := by
    intro p hp
    have h := hz' p hp
    have d1 := (one_add_mul_pos hl hp.le).ne'
    have d2 := (one_add_mul_pos hl' hp.le).ne'
    unfold g at h
    field_simp at h
    linear_combination h
  -- (iii) constant term: ξ² = ξ'²
  have hsq : xi ^ 2 = xi' ^ 2 := by
    have k1 := key 1 one_pos
    have k2 := key 2 two_pos
    have k3 := key 3 (by norm_num)
    have := const_coeff_of_three_points (A := xi ^ 2 - xi' ^ 2)
      (B := xi ^ 2 * (tR + tl') - xi' ^ 2 * (tR' + tl))
      (C := xi ^ 2 * (tR * tl') - xi' ^ 2 * (tR' * tl))
      (by linear_combination k1) (by linear_combination k2) (by linear_combination k3)
    linarith
  have hxi_eq : xi = xi' := by
    have h : (xi - xi') * (xi + xi') = 0 := by linear_combination hsq
    rcases mul_eq_zero.1 h with h | h
    · linarith
    · linarith
  -- (iv) g = g' on p > 0
  have hg : ∀ p : ℝ, 0 < p → g tR tl p = g tR' tl' p := by
    intro p hp
    have h := hz' p hp
    rw [hxi_eq] at h
    exact mul_left_cancel₀ (pow_ne_zero 2 hxi'.ne') h
  -- (v) b² = b'² from B² = B'² at p = 1
  have hb_eq : b = b' := by
    have h := hB 1 one_pos
    rw [hg 1 one_pos] at h
    have hgpos := g_pos hR' hl' zero_le_one
    have hb2 : b ^ 2 = b' ^ 2 := by
      field_simp at h
      linarith
    have h' : (b - b') * (b + b') = 0 := by linear_combination hb2
    rcases mul_eq_zero.1 h' with h' | h'
    · linarith
    · linarith
  exact ⟨hxi_eq, hb_eq, hg⟩

/-! ## Fibres of Theorem 1 -/

/-- A parameter point `θ = (ξ₁, b, τ_R, τ_ℓ)`. -/
@[ext] structure Param where
  xi : ℝ
  b : ℝ
  tR : ℝ
  tl : ℝ

/-- `Θ = {ξ₁ > 0, b > 0, τ_R ≥ 0, τ_ℓ ≥ 0}`. -/
def Theta : Set Param := {θ | 0 < θ.xi ∧ 0 < θ.b ∧ 0 ≤ θ.tR ∧ 0 ≤ θ.tl}

/-- The Fourier resonance `Δ = {τ_R = τ_ℓ}`. -/
def Delta : Set Param := {θ | θ ∈ Theta ∧ θ.tR = θ.tl}

/-- `D₂ = Δ ∩ {b = b_s}`. -/
def D2 (bs : ℝ) : Set Param := {θ | θ ∈ Delta ∧ θ.b = bs}

/-- The fibre `𝔉(θ)` of Theorem 1: `{θ}` off `Δ`; `{(ξ₁, b, t, t) : t ≥ 0}` on `Δ ∖ D₂`;
`{(ξ, b_s, t, t) : ξ > 0, t ≥ 0} = D₂` on `D₂`. -/
def Fib (bs : ℝ) (θ : Param) : Set Param :=
  {θ' | θ' ∈ Theta ∧
    ((θ ∉ Delta ∧ θ' = θ) ∨
      (θ ∈ Delta ∧ θ ∉ D2 bs ∧ θ'.xi = θ.xi ∧ θ'.b = θ.b ∧ θ'.tR = θ'.tl) ∨
      (θ ∈ D2 bs ∧ θ' ∈ D2 bs))}

/-- **Steps 4–5 combined.** For `θ ∈ Θ ∖ D₂` and `θ' ∈ Θ`: `ξ₁' = ξ₁`, `b' = b` and `g' = g` on
`p > 0` iff `θ' ∈ 𝔉(θ)`. -/
theorem mem_fibre_iff_of_not_D2 {bs : ℝ} {θ θ' : Param} (hθ : θ ∈ Theta) (hθ' : θ' ∈ Theta)
    (hD2 : θ ∉ D2 bs) :
    (θ'.xi = θ.xi ∧ θ'.b = θ.b ∧ ∀ p : ℝ, 0 < p → g θ'.tR θ'.tl p = g θ.tR θ.tl p) ↔
      θ' ∈ Fib bs θ := by
  obtain ⟨_, _, hR, hl⟩ := hθ
  obtain ⟨_, _, hR', hl'⟩ := hθ'
  constructor
  · rintro ⟨hxi, hb, hg⟩
    refine ⟨⟨by assumption, by assumption, hR', hl'⟩, ?_⟩
    have h5 := step5_g hl hl' (fun p hp => (hg p hp).symm)
    by_cases hres : θ.tR = θ.tl
    · have hmem : θ ∈ Delta := ⟨⟨by assumption, by assumption, hR, hl⟩, hres⟩
      exact Or.inr (Or.inl ⟨hmem, hD2, hxi, hb, h5.2 hres⟩)
    · have hnot : θ ∉ Delta := fun h => hres h.2
      obtain ⟨h1, h2⟩ := h5.1 hres
      exact Or.inl ⟨hnot, Param.ext hxi hb h1 h2⟩
  · rintro ⟨_, h | h | h⟩
    · obtain ⟨_, rfl⟩ := h
      exact ⟨rfl, rfl, fun _ _ => rfl⟩
    · obtain ⟨hmem, _, hxi, hb, hres'⟩ := h
      refine ⟨hxi, hb, fun p hp => ?_⟩
      rw [hres', hmem.2, g_diag hl' hp.le, g_diag hl hp.le]
    · exact absurd h.1 hD2

/-- **Theorem 1, (a) ⇔ (c): case assembly.**  `Same θ θ'` stands for "equal responses on a
frequency set with an accumulation point".  The hypotheses are the analytic steps of Proof A:
`hS4` (Steps 1–4, `θ ∉ D₂`), `hS6` (Step 6 with Step 2: `θ ∈ D₂` forces `θ' ∈ D₂`) and `hc`
((c) ⇒ (b): points of the fibre give equal responses). -/
theorem thm1_of_steps (bs : ℝ) (Same : Param → Param → Prop)
    (hS4 : ∀ θ θ', θ ∈ Theta → θ' ∈ Theta → θ ∉ D2 bs → Same θ θ' →
      θ'.xi = θ.xi ∧ θ'.b = θ.b ∧ ∀ p : ℝ, 0 < p → g θ'.tR θ'.tl p = g θ.tR θ.tl p)
    (hS6 : ∀ θ θ', θ ∈ Theta → θ' ∈ Theta → θ ∈ D2 bs → Same θ θ' → θ' ∈ D2 bs)
    (hc : ∀ θ θ', θ ∈ Theta → θ' ∈ Fib bs θ → Same θ θ') :
    ∀ θ θ', θ ∈ Theta → θ' ∈ Theta → (Same θ θ' ↔ θ' ∈ Fib bs θ) := by
  intro θ θ' hθ hθ'
  refine ⟨fun h => ?_, hc θ θ' hθ⟩
  by_cases hD2 : θ ∈ D2 bs
  · exact ⟨hθ', Or.inr (Or.inr ⟨hD2, hS6 θ θ' hθ hθ' hD2 h⟩)⟩
  · exact (mem_fibre_iff_of_not_D2 hθ hθ' hD2).1 (hS4 θ θ' hθ hθ' hD2 h)

/-! ## The kernel direction `(0,0,1,1)` on the resonance -/

lemma hasDerivAt_affine (p x : ℝ) : HasDerivAt (fun t : ℝ => 1 + t * p) p x := by
  have h := ((hasDerivAt_id' x).mul_const p).const_add (1 : ℝ)
  simpa using h

/-- `∂g/∂τ_R = p/(1 + τ_ℓ p)`. -/
theorem hasDerivAt_g_tR (tR tl p : ℝ) :
    HasDerivAt (fun t => g t tl p) (p / (1 + tl * p)) tR := by
  unfold g
  exact (hasDerivAt_affine p tR).div_const (1 + tl * p)

/-- `∂g/∂τ_ℓ = -(1 + τ_R p) p/(1 + τ_ℓ p)²` where `1 + τ_ℓ p ≠ 0`. -/
theorem hasDerivAt_g_tl {tR tl p : ℝ} (h : 1 + tl * p ≠ 0) :
    HasDerivAt (fun t => g tR t p) (-(1 + tR * p) * p / (1 + tl * p) ^ 2) tl := by
  unfold g
  have h2 := (hasDerivAt_const tl (1 + tR * p)).div (hasDerivAt_affine p tl) h
  convert h2 using 1
  ring

/-- On the resonance `τ_R = τ_ℓ = t` (with `1 + t p ≠ 0`), the two partial derivatives of `g`
add to zero: `(0,0,1,1)` is tangent to the level sets of `g`. -/
theorem g_tangent_resonance {t p : ℝ} (h : 1 + t * p ≠ 0) :
    ∃ dR dl : ℝ, HasDerivAt (fun s => g s t p) dR t ∧ HasDerivAt (fun s => g t s p) dl t ∧
      dR + dl = 0 := by
  refine ⟨_, _, hasDerivAt_g_tR t t p, hasDerivAt_g_tl h, ?_⟩
  field_simp
  ring

end FdtrGauge
