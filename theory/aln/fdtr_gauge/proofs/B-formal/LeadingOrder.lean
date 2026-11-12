import ResearchProofs.FdtrGauge.MomentMinors
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Data.Finset.Powerset

/-!
# Theorem N2, step (3): the leading term of the two-frequency determinant

Campaign `20261111-fdtr-gauge-memory`, `proofs/B-analytic.md`, §8 (Theorem N2).

With `s = e^{iπ/4} σ` and frequencies `σ₁ = σ`, `σ₂ = λ σ`, the Cauchy–Binet expansion reads
`det 𝒥₂ = Σ_K σ^{|K|} t_K(λ) G_K`, `|K| = Σ_{k ∈ K} k`, where `t_K(λ)` is the `4 × 4`
determinant with rows `cos(kπ/4)`, `sin(kπ/4)`, `λ^k cos(kπ/4)`, `λ^k sin(kπ/4)` (`k ∈ K`).
This file checks the finite facts used in step (3):

* `t_0123`, `t_0124`, `t_0234`: the three `t_K` printed in `07w_wolfram_crosscheck.out`
  (written with `√2/2 = 1/√2`), and `t_0124_ne_zero`, `t_0234_ne_zero` for `λ > 0`, `λ ≠ 1`.
* `index_sets_le_nine`, `index_sum_six`, `index_sum_seven`, `index_sum_ge_six`: the
  four-element index sets with `|K| ≤ 9` are exactly the seven listed in the proof; `|K| ≥ 6`,
  with equality only for `{0,1,2,3}`; `|K| = 7` only for `{0,1,2,4}`; all listed sets except
  `{0,2,3,4}` contain `{0,1}`.
* `leading_coeff_ne_zero`: off the resonance (`ξ₁, b, b_s > 0`), `t_{0124}(λ) G_{0124} ≠ 0` if
  `b ≠ b_s`; `matched_low_order`: if `b = b_s`, every minor with `|K| ≤ 9` other than
  `G_{0234}` vanishes and `t_{0234}(λ) G_{0234} ≠ 0`.

Not formalised: the Cauchy–Binet expansion of `det 𝒥₂` itself (an infinite sum), the remainder
estimate of step (4), and Lemma F.  `G_K` for a set `K` is `G` of its increasing enumeration.
-/

namespace FdtrGauge.LeadingOrder

open Real FdtrGauge.Moments

/-- The matrix of `t_K(λ)`. -/
noncomputable def tmat (K : Fin 4 → ℕ) (lam : ℝ) : Matrix (Fin 4) (Fin 4) ℝ :=
  Matrix.of ![fun j => Real.cos (K j * π / 4), fun j => Real.sin (K j * π / 4),
    fun j => lam ^ K j * Real.cos (K j * π / 4), fun j => lam ^ K j * Real.sin (K j * π / 4)]

/-- `t_K(λ) = det tmat`. -/
noncomputable def t (K : Fin 4 → ℕ) (lam : ℝ) : ℝ := (tmat K lam).det

lemma cos_0 : Real.cos ((0 : ℝ) * π / 4) = 1 := by simp
lemma sin_0 : Real.sin ((0 : ℝ) * π / 4) = 0 := by simp
lemma cos_1 : Real.cos ((1 : ℝ) * π / 4) = √2 / 2 := by
  rw [show (1 : ℝ) * π / 4 = π / 4 by ring]; exact Real.cos_pi_div_four
lemma sin_1 : Real.sin ((1 : ℝ) * π / 4) = √2 / 2 := by
  rw [show (1 : ℝ) * π / 4 = π / 4 by ring]; exact Real.sin_pi_div_four
lemma cos_2 : Real.cos ((2 : ℝ) * π / 4) = 0 := by
  rw [show (2 : ℝ) * π / 4 = π / 2 by ring]; exact Real.cos_pi_div_two
lemma sin_2 : Real.sin ((2 : ℝ) * π / 4) = 1 := by
  rw [show (2 : ℝ) * π / 4 = π / 2 by ring]; exact Real.sin_pi_div_two
lemma cos_3 : Real.cos ((3 : ℝ) * π / 4) = -(√2 / 2) := by
  rw [show (3 : ℝ) * π / 4 = π - π / 4 by ring, Real.cos_pi_sub, Real.cos_pi_div_four]
lemma sin_3 : Real.sin ((3 : ℝ) * π / 4) = √2 / 2 := by
  rw [show (3 : ℝ) * π / 4 = π - π / 4 by ring, Real.sin_pi_sub, Real.sin_pi_div_four]
lemma cos_4 : Real.cos ((4 : ℝ) * π / 4) = -1 := by
  rw [show (4 : ℝ) * π / 4 = π by ring]; exact Real.cos_pi
lemma sin_4 : Real.sin ((4 : ℝ) * π / 4) = 0 := by
  rw [show (4 : ℝ) * π / 4 = π by ring]; exact Real.sin_pi

/-- `t_{0123}(λ) = λ (λ − 1)² (λ² + 1)/2`. -/
theorem t_0123 (lam : ℝ) : t ![0, 1, 2, 3] lam = lam * (lam - 1) ^ 2 * (lam ^ 2 + 1) / 2 := by
  rw [t, det_fin_four']
  simp only [tmat, Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons,
    Nat.cast_zero, Nat.cast_one, Nat.cast_ofNat, cos_0, sin_0, cos_1, sin_1, cos_2, sin_2, cos_3,
    sin_3, pow_zero, pow_one]
  -- every Leibniz term carries one `√2` from each of the columns `k = 1` and `k = 3`
  have h2 : √2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  linear_combination (lam * (lam - 1) ^ 2 * (lam ^ 2 + 1) / 4) * h2

/-- `t_{0124}(λ) = λ (λ − 1)² (λ + 1)(λ² + 1)/√2` (here `1/√2 = √2/2`). -/
theorem t_0124 (lam : ℝ) :
    t ![0, 1, 2, 4] lam = √2 * lam * (lam - 1) ^ 2 * (lam + 1) * (lam ^ 2 + 1) / 2 := by
  rw [t, det_fin_four']
  simp only [tmat, Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons,
    Nat.cast_zero, Nat.cast_one, Nat.cast_ofNat, cos_0, sin_0, cos_1, sin_1, cos_2, sin_2, cos_4,
    sin_4, pow_zero, pow_one]
  ring

/-- `t_{0234}(λ) = λ² (λ − 1)² (λ + 1)(λ² + 1)/√2`. -/
theorem t_0234 (lam : ℝ) :
    t ![0, 2, 3, 4] lam = √2 * lam ^ 2 * (lam - 1) ^ 2 * (lam + 1) * (lam ^ 2 + 1) / 2 := by
  rw [t, det_fin_four']
  simp only [tmat, Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons,
    Nat.cast_zero, Nat.cast_ofNat, cos_0, sin_0, cos_2, sin_2, cos_3, sin_3, cos_4, sin_4,
    pow_zero]
  ring

theorem t_0124_ne_zero {lam : ℝ} (h0 : 0 < lam) (h1 : lam ≠ 1) : t ![0, 1, 2, 4] lam ≠ 0 := by
  rw [t_0124]
  have hs : 0 < √2 := by positivity
  have hl : (lam - 1) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.2 h1)
  have : √2 * lam * (lam - 1) ^ 2 * (lam + 1) * (lam ^ 2 + 1) ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hs.ne' h0.ne') hl) (by positivity))
    positivity
  exact div_ne_zero this two_ne_zero

theorem t_0234_ne_zero {lam : ℝ} (h0 : 0 < lam) (h1 : lam ≠ 1) : t ![0, 2, 3, 4] lam ≠ 0 := by
  rw [t_0234]
  have hs : 0 < √2 := by positivity
  have hl : (lam - 1) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.2 h1)
  have : √2 * lam ^ 2 * (lam - 1) ^ 2 * (lam + 1) * (lam ^ 2 + 1) ≠ 0 := by
    apply mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero hs.ne' (by positivity)) hl)
      (by positivity))
    positivity
  exact div_ne_zero this two_ne_zero

/-! ## Index sets with `|K| ≤ 9` -/

/-- The seven index sets of step (3). -/
def lowSets : Finset (Finset ℕ) :=
  {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 2, 5}, {0, 1, 3, 4}, {0, 1, 2, 6}, {0, 1, 3, 5}, {0, 2, 3, 4}}

lemma lowSets_eq :
    ((Finset.range 10).powersetCard 4).filter (fun K => K.sum id ≤ 9) = lowSets := by
  decide

/-- **The four-element index sets with `|K| ≤ 9` are exactly the seven listed.** -/
theorem index_sets_le_nine (K : Finset ℕ) (hK : K.card = 4) :
    K.sum id ≤ 9 ↔ K ∈ lowSets := by
  constructor
  · intro hs
    have hsub : K ⊆ Finset.range 10 := by
      intro k hk
      have h1 : id k ≤ K.sum id := Finset.single_le_sum (fun i _ => Nat.zero_le (id i)) hk
      rw [Finset.mem_range]
      exact Nat.lt_succ_of_le (le_trans h1 hs)
    rw [← lowSets_eq, Finset.mem_filter, Finset.mem_powersetCard]
    exact ⟨⟨hsub, hK⟩, hs⟩
  · intro h
    rw [← lowSets_eq, Finset.mem_filter] at h
    exact h.2

/-- `|K| ≥ 6` for every four-element `K ⊂ ℕ`. -/
theorem index_sum_ge_six (K : Finset ℕ) (hK : K.card = 4) : 6 ≤ K.sum id := by
  by_contra h
  have hmem := (index_sets_le_nine K hK).1 (by omega)
  have : ∀ L ∈ lowSets, 6 ≤ L.sum id := by decide
  exact h (this K hmem)

/-- `|K| = 6` only for `K = {0,1,2,3}`. -/
theorem index_sum_six (K : Finset ℕ) (hK : K.card = 4) (hs : K.sum id = 6) : K = {0, 1, 2, 3} := by
  have hmem := (index_sets_le_nine K hK).1 (by omega)
  have : ∀ L ∈ lowSets, L.sum id = 6 → L = {0, 1, 2, 3} := by decide
  exact this K hmem hs

/-- `|K| = 7` only for `K = {0,1,2,4}`. -/
theorem index_sum_seven (K : Finset ℕ) (hK : K.card = 4) (hs : K.sum id = 7) :
    K = {0, 1, 2, 4} := by
  have hmem := (index_sets_le_nine K hK).1 (by omega)
  have : ∀ L ∈ lowSets, L.sum id = 7 → L = {0, 1, 2, 4} := by decide
  exact this K hmem hs

/-- All listed sets except `{0,2,3,4}` contain `{0,1}`. -/
theorem lowSets_contain_01 : ∀ L ∈ lowSets, L ≠ {0, 2, 3, 4} → ({0, 1} : Finset ℕ) ⊆ L := by
  decide

/-! ## The leading coefficient -/

/-- **Generic case.** Off the resonance with `b ≠ b_s`, the coefficient of `σ⁷`,
`t_{0124}(λ) G_{0124}`, is nonzero for `λ > 0`, `λ ≠ 1`. -/
theorem leading_coeff_ne_zero {bs xi b tR tl lam : ℝ} (hbs : 0 < bs) (hxi : 0 < xi) (hb : 0 < b)
    (hm : b ≠ bs) (hres : tR ≠ tl) (h0 : 0 < lam) (h1 : lam ≠ 1) :
    t ![0, 1, 2, 4] lam * G bs ![0, 1, 2, 4] xi b tR tl ≠ 0 := by
  refine mul_ne_zero (t_0124_ne_zero h0 h1) ?_
  rw [G_0124 hb.ne' hbs.ne']
  have h2 : b ^ 2 - bs ^ 2 ≠ 0 := by
    intro h
    have h' : (b - bs) * (b + bs) = 0 := by linear_combination h
    rcases mul_eq_zero.1 h' with h' | h'
    · exact hm (sub_eq_zero.1 h')
    · linarith
  have h3 : tl - tR ≠ 0 := sub_ne_zero.2 (Ne.symm hres)
  have h4 : xi ^ 4 ≠ 0 := pow_ne_zero 4 hxi.ne'
  have h5 : b ^ 4 * bs ^ 3 ≠ 0 := by positivity
  exact div_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne_zero h4) h2) h3) h5

/-- **Matched case `b = b_s`.** Every minor with `|K| ≤ 9` other than `G_{0234}` vanishes, and the
coefficient of `σ⁹`, `t_{0234}(λ) G_{0234}`, is nonzero off the resonance. -/
theorem matched_low_order {bs xi tR tl lam : ℝ} (hbs : 0 < bs) (hxi : 0 < xi) (hres : tR ≠ tl)
    (h0 : 0 < lam) (h1 : lam ≠ 1) :
    G bs ![0, 1, 2, 3] xi bs tR tl = 0 ∧ G bs ![0, 1, 2, 4] xi bs tR tl = 0 ∧
      G bs ![0, 1, 2, 5] xi bs tR tl = 0 ∧ G bs ![0, 1, 3, 4] xi bs tR tl = 0 ∧
      G bs ![0, 1, 2, 6] xi bs tR tl = 0 ∧ G bs ![0, 1, 3, 5] xi bs tR tl = 0 ∧
      t ![0, 2, 3, 4] lam * G bs ![0, 2, 3, 4] xi bs tR tl ≠ 0 := by
  refine ⟨G_01_matched hbs.ne' _ _, G_01_matched hbs.ne' _ _, G_01_matched hbs.ne' _ _,
    G_01_matched hbs.ne' _ _, G_01_matched hbs.ne' _ _, G_01_matched hbs.ne' _ _, ?_⟩
  refine mul_ne_zero (t_0234_ne_zero h0 h1) ?_
  rw [G_0234_matched hbs.ne']
  have h2 : xi ^ 4 ≠ 0 := pow_ne_zero 4 hxi.ne'
  have h3 : (tR - tl) ^ 2 ≠ 0 := pow_ne_zero 2 (sub_ne_zero.2 hres)
  have h4 : bs ^ 5 ≠ 0 := pow_ne_zero 5 hbs.ne'
  exact div_ne_zero (mul_ne_zero (mul_ne_zero (by norm_num) h2) h3) h4

end FdtrGauge.LeadingOrder
