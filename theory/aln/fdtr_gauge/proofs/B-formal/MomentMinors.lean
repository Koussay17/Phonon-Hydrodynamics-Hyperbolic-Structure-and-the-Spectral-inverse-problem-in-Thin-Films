import ResearchProofs.FdtrGauge.MomentCoefficients
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.LinearAlgebra.Matrix.Nondegenerate
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# Minors of the moment Jacobian (Theorem 2, Proof C; Theorem N2, step (3))

Campaign `20261111-fdtr-gauge-memory`, `proofs/B-analytic.md`, §5 and §8.

`G_K = det [∇c_k]_{k ∈ K}` for `K = (k₀, k₁, k₂, k₃)`, with the coefficients and gradients of
`MomentCoefficients.lean` (taken as given rational functions; see the header there).

* `G_0123`: `G_{0123} ≡ 0`.
* `G_0124`: `G_{0124} = 2 ξ₁⁴ (b² − b_s²)(τ_ℓ − τ_R)/(b⁴ b_s³)`.
* `G_0234`: the general form of `G_{0234}`, and `G_0234_matched`:
  `G_{0234}|_{b = b_s} = −2 ξ₁⁴ (τ_R − τ_ℓ)²/b_s⁵`.
* `G_01_matched`: at `b = b_s`, `∇c₀ ∥ ∇c₁`, so every `G_{01k₂k₃}` vanishes.
* `G_0125`, `G_0134`, `G_0126`, `G_0135`: the remaining minors with `|K| ≤ 9` in
  `08_low_frequency_two_point.out` (a).
* `moment_minor_ne_zero`, `moment_kernel_trivial`: off the resonance `τ_R ≠ τ_ℓ`
  (`ξ₁ > 0`, `b > 0`, `b_s > 0`, any real times, so including `τ_R = 0` or `τ_ℓ = 0`), a vector
  annihilated by `∇c₀, …, ∇c₄` is zero.  This is the algebraic content of Proof C.
* `resonance_direction`: on `τ_R = τ_ℓ`, `(0, 0, 1, 1)` is annihilated by every `∇c_k`,
  `k = 0, …, 6` (consistency with the kernel of Theorem 2).

Not formalised: Lemma 3 (iv) (`D_c Z ≡ 0` implies `D_c c_k = 0`), i.e. the passage from the
moment Jacobian to the function-level Jacobian.
-/

set_option linter.style.longLine false
set_option linter.unusedSimpArgs false

namespace FdtrGauge.Moments

open Matrix

/-- Laplace expansion of a `4 × 4` determinant along the first row. -/
theorem det_fin_four' {R : Type*} [CommRing R] (A : Matrix (Fin 4) (Fin 4) R) :
    A.det =
      A 0 0 * (A 1 1 * (A 2 2 * A 3 3 - A 2 3 * A 3 2) - A 1 2 * (A 2 1 * A 3 3 - A 2 3 * A 3 1)
                + A 1 3 * (A 2 1 * A 3 2 - A 2 2 * A 3 1))
      - A 0 1 * (A 1 0 * (A 2 2 * A 3 3 - A 2 3 * A 3 2) - A 1 2 * (A 2 0 * A 3 3 - A 2 3 * A 3 0)
                + A 1 3 * (A 2 0 * A 3 2 - A 2 2 * A 3 0))
      + A 0 2 * (A 1 0 * (A 2 1 * A 3 3 - A 2 3 * A 3 1) - A 1 1 * (A 2 0 * A 3 3 - A 2 3 * A 3 0)
                + A 1 3 * (A 2 0 * A 3 1 - A 2 1 * A 3 0))
      - A 0 3 * (A 1 0 * (A 2 1 * A 3 2 - A 2 2 * A 3 1) - A 1 1 * (A 2 0 * A 3 2 - A 2 2 * A 3 0)
                + A 1 2 * (A 2 0 * A 3 1 - A 2 1 * A 3 0)) := by
  rw [Matrix.det_succ_row_zero]
  simp only [Fin.sum_univ_succ, Fin.sum_univ_zero, det_fin_three, submatrix_apply,
    Fin.succ_zero_eq_one, Fin.succ_one_eq_two]
  simp [Fin.succAbove, Fin.lt_def]
  ring_nf

lemma coeff_0 (bs : ℝ) : coeff bs 0 = c0 bs := rfl
lemma coeff_1 (bs : ℝ) : coeff bs 1 = c1 bs := rfl
lemma coeff_2 (bs : ℝ) : coeff bs 2 = c2 bs := rfl
lemma coeff_3 (bs : ℝ) : coeff bs 3 = c3 bs := rfl
lemma coeff_4 (bs : ℝ) : coeff bs 4 = c4 bs := rfl
lemma coeff_5 (bs : ℝ) : coeff bs 5 = c5 bs := rfl
lemma coeff_6 (bs : ℝ) : coeff bs 6 = c6 bs := rfl

/-- `G_K = det [∇c_{K 0}; ∇c_{K 1}; ∇c_{K 2}; ∇c_{K 3}]`. -/
noncomputable def G (bs : ℝ) (K : Fin 4 → Fin 7) (xi b tR tl : ℝ) : ℝ :=
  Matrix.det (Matrix.of fun i => grad (coeff bs (K i)) xi b tR tl)

lemma G_eq (bs : ℝ) (K : Fin 4 → Fin 7) (xi b tR tl : ℝ) :
    G bs K xi b tR tl =
      Matrix.det (Matrix.of ![grad (coeff bs (K 0)) xi b tR tl, grad (coeff bs (K 1)) xi b tR tl,
        grad (coeff bs (K 2)) xi b tR tl, grad (coeff bs (K 3)) xi b tR tl]) := by
  unfold G
  congr 1
  ext i j
  fin_cases i <;> rfl

/-- Evaluation of `G_K` from four explicit rows. -/
lemma det_rows (r0 r1 r2 r3 : Fin 4 → ℝ) :
    Matrix.det (Matrix.of ![r0, r1, r2, r3]) =
      r0 0 * (r1 1 * (r2 2 * r3 3 - r2 3 * r3 2) - r1 2 * (r2 1 * r3 3 - r2 3 * r3 1)
                + r1 3 * (r2 1 * r3 2 - r2 2 * r3 1))
      - r0 1 * (r1 0 * (r2 2 * r3 3 - r2 3 * r3 2) - r1 2 * (r2 0 * r3 3 - r2 3 * r3 0)
                + r1 3 * (r2 0 * r3 2 - r2 2 * r3 0))
      + r0 2 * (r1 0 * (r2 1 * r3 3 - r2 3 * r3 1) - r1 1 * (r2 0 * r3 3 - r2 3 * r3 0)
                + r1 3 * (r2 0 * r3 1 - r2 1 * r3 0))
      - r0 3 * (r1 0 * (r2 1 * r3 2 - r2 2 * r3 1) - r1 1 * (r2 0 * r3 2 - r2 2 * r3 0)
                + r1 2 * (r2 0 * r3 1 - r2 1 * r3 0)) := by
  rw [det_fin_four']
  simp

/-! ## The minors used in Proof C and Theorem N2 -/

/-- **`G_{0123} ≡ 0`.** -/
theorem G_0123 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 2, 3] xi b tR tl = 0 := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c2 hb hbs, grad_c3 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- **`G_{0124} = 2 ξ₁⁴ (b² − b_s²)(τ_ℓ − τ_R)/(b⁴ b_s³)`.** -/
theorem G_0124 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 2, 4] xi b tR tl = 2 * xi ^ 4 * (b ^ 2 - bs ^ 2) * (tl - tR) / (b ^ 4 * bs ^ 3) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c2 hb hbs, grad_c4 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- **`G_{0234}`, general form** (`07_moment_jacobian.out`). -/
theorem G_0234 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 2, 3, 4] xi b tR tl =
      2 * xi ^ 4 * (tR - tl) * (2 * b ^ 4 * xi ^ 2 - 5 * b ^ 2 * bs ^ 2 * xi ^ 2 - 3 * bs ^ 4 * tR
        + 3 * bs ^ 4 * tl + 3 * bs ^ 4 * xi ^ 2) / (3 * b ^ 4 * bs ^ 5) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c2 hb hbs, grad_c3 hb hbs, grad_c4 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- **`G_{0234}|_{b = b_s} = −2 ξ₁⁴ (τ_R − τ_ℓ)²/b_s⁵`.** -/
theorem G_0234_matched {bs xi tR tl : ℝ} (hbs : bs ≠ 0) :
    G bs ![0, 2, 3, 4] xi bs tR tl = -2 * xi ^ 4 * (tR - tl) ^ 2 / bs ^ 5 := by
  rw [G_0234 hbs hbs]
  field_simp
  ring

/-- A `4 × 4` determinant whose first two rows are supported in the second column vanishes. -/
lemma det_rows_col1 (a a' : ℝ) (r2 r3 : Fin 4 → ℝ) :
    Matrix.det (Matrix.of ![![0, a, 0, 0], ![0, a', 0, 0], r2, r3]) = 0 := by
  rw [det_rows]
  simp

/-- **Matched case.** At `b = b_s`, `∇c₀ = (0, −2ξ₁/b_s², 0, 0)` and `∇c₁ = (0, 2ξ₁²/b_s², 0, 0)`
are parallel, so `G_K = 0` whenever `{0, 1} ⊂ K`. -/
theorem G_01_matched {bs xi tR tl : ℝ} (hbs : bs ≠ 0) (k2 k3 : Fin 7) :
    G bs ![0, 1, k2, k3] xi bs tR tl = 0 := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, coeff_0, coeff_1]
  have e0 : -xi * (bs ^ 2 + bs ^ 2) / (bs ^ 2 * bs ^ 2) = -2 * xi / bs ^ 2 := by
    field_simp
    ring
  have e1 : 2 * bs * xi ^ 2 / bs ^ 3 = 2 * xi ^ 2 / bs ^ 2 := by
    field_simp
  have h0 : grad (c0 bs) xi bs tR tl = ![0, -2 * xi / bs ^ 2, 0, 0] := by
    rw [grad_c0 hbs hbs, e0, sub_self, neg_zero, zero_mul, zero_div]
  have h1 : grad (c1 bs) xi bs tR tl = ![0, 2 * xi ^ 2 / bs ^ 2, 0, 0] := by
    rw [grad_c1 hbs hbs, e1, sub_self, mul_zero, zero_mul, zero_div]
  rw [h0, h1]
  exact det_rows_col1 _ _ _ _

/-- `G_{0125}` (`08_low_frequency_two_point.out` (a)). -/
theorem G_0125 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 2, 5] xi b tR tl =
      2 * xi ^ 5 * (b - bs) * (b + bs) * (tR - tl) / (b ^ 3 * bs ^ 4) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c2 hb hbs, grad_c5 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- `G_{0134}` (`08_low_frequency_two_point.out` (a)). -/
theorem G_0134 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 3, 4] xi b tR tl =
      2 * xi ^ 5 * (b - bs) * (b + bs) * (tR - tl) / (b ^ 3 * bs ^ 4) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c3 hb hbs, grad_c4 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- `G_{0126}` (`08_low_frequency_two_point.out` (a)). -/
theorem G_0126 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 2, 6] xi b tR tl =
      -4 * xi ^ 4 * (b - bs) * (b + bs) * (tR - tl) * (2 * b ^ 2 * xi ^ 2 - 3 * bs ^ 2 * tl
        - bs ^ 2 * xi ^ 2) / (3 * b ^ 4 * bs ^ 5) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c2 hb hbs, grad_c6 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-- `G_{0135}` (`08_low_frequency_two_point.out` (a)). -/
theorem G_0135 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    G bs ![0, 1, 3, 5] xi b tR tl =
      -2 * xi ^ 6 * (b - bs) * (b + bs) * (tR - tl) / (b ^ 2 * bs ^ 5) := by
  rw [G_eq]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons, coeff_0, coeff_1, coeff_2,
    coeff_3, coeff_4, coeff_5, coeff_6]
  rw [grad_c0 hb hbs, grad_c1 hb hbs, grad_c3 hb hbs, grad_c5 hb hbs, det_rows]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.cons_val_three, Matrix.head_cons, Matrix.tail_cons]
  field_simp
  ring

/-! ## Rank four off the resonance (Proof C) -/

/-- **Off the resonance one of `G_{0124}`, `G_{0234}` is nonzero** (`ξ₁ > 0`, `b > 0`, `b_s > 0`;
any real `τ_R ≠ τ_ℓ`): `G_{0124} ≠ 0` if `b ≠ b_s`, and `G_{0234} ≠ 0` if `b = b_s`. -/
theorem moment_minor_ne_zero {bs xi b tR tl : ℝ} (hbs : 0 < bs) (hxi : 0 < xi) (hb : 0 < b)
    (hres : tR ≠ tl) :
    G bs ![0, 1, 2, 4] xi b tR tl ≠ 0 ∨ G bs ![0, 2, 3, 4] xi b tR tl ≠ 0 := by
  have hd : tR - tl ≠ 0 := sub_ne_zero.2 hres
  by_cases hm : b = bs
  · right
    subst hm
    rw [G_0234_matched hbs.ne']
    have h1 : xi ^ 4 ≠ 0 := pow_ne_zero 4 hxi.ne'
    have h2 : (tR - tl) ^ 2 ≠ 0 := pow_ne_zero 2 hd
    have h3 : b ^ 5 ≠ 0 := pow_ne_zero 5 hb.ne'
    exact div_ne_zero (mul_ne_zero (mul_ne_zero (by norm_num) h1) h2) h3
  · left
    rw [G_0124 hb.ne' hbs.ne']
    have h2 : b ^ 2 - bs ^ 2 ≠ 0 := by
      intro h
      have h' : (b - bs) * (b + bs) = 0 := by linear_combination h
      rcases mul_eq_zero.1 h' with h' | h'
      · exact hm (sub_eq_zero.1 h')
      · linarith
    have h3 : tl - tR ≠ 0 := sub_ne_zero.2 (Ne.symm hres)
    have h1 : xi ^ 4 ≠ 0 := pow_ne_zero 4 hxi.ne'
    have h4 : b ^ 4 * bs ^ 3 ≠ 0 := by positivity
    exact div_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne_zero h1) h2) h3) h4

/-- **Proof C (algebraic part).** Off the resonance, a direction `c ∈ ℝ⁴` with
`D_c c_k = ∇c_k · c = 0` for `k = 0, …, 4` is zero. -/
theorem moment_kernel_trivial {bs xi b tR tl : ℝ} (hbs : 0 < bs) (hxi : 0 < xi) (hb : 0 < b)
    (hres : tR ≠ tl) (c : Fin 4 → ℝ)
    (hc : ∀ k : Fin 7, (k : ℕ) ≤ 4 → grad (coeff bs k) xi b tR tl ⬝ᵥ c = 0) : c = 0 := by
  have key : ∀ K : Fin 4 → Fin 7, (∀ i, (K i : ℕ) ≤ 4) → G bs K xi b tR tl ≠ 0 → c = 0 := by
    intro K hK hG
    apply Matrix.eq_zero_of_mulVec_eq_zero hG
    ext i
    simp only [Matrix.mulVec, Matrix.of_apply, Pi.zero_apply]
    exact hc (K i) (hK i)
  rcases moment_minor_ne_zero hbs hxi hb hres with h | h
  · exact key _ (fun i => by fin_cases i <;> decide) h
  · exact key _ (fun i => by fin_cases i <;> decide) h

/-! ## Consistency on the resonance -/

/-- On `τ_R = τ_ℓ = t`, `(0, 0, 1, 1)` is annihilated by `∇c_k` for every `k = 0, …, 6`. -/
theorem resonance_direction {bs xi b t : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) (k : Fin 7) :
    grad (coeff bs k) xi b t t 2 + grad (coeff bs k) xi b t t 3 = 0 := by
  fin_cases k
  · change grad (c0 bs) xi b t t 2 + grad (c0 bs) xi b t t 3 = 0
    rw [grad_c0 hb hbs]; simp
  · change grad (c1 bs) xi b t t 2 + grad (c1 bs) xi b t t 3 = 0
    rw [grad_c1 hb hbs]; simp
  · change grad (c2 bs) xi b t t 2 + grad (c2 bs) xi b t t 3 = 0
    rw [grad_c2 hb hbs]; simp; ring
  · change grad (c3 bs) xi b t t 2 + grad (c3 bs) xi b t t 3 = 0
    rw [grad_c3 hb hbs]; simp; ring
  · change grad (c4 bs) xi b t t 2 + grad (c4 bs) xi b t t 3 = 0
    rw [grad_c4 hb hbs]; simp; field_simp; ring
  · change grad (c5 bs) xi b t t 2 + grad (c5 bs) xi b t t 3 = 0
    rw [grad_c5 hb hbs]; simp; field_simp; ring
  · change grad (c6 bs) xi b t t 2 + grad (c6 bs) xi b t t 3 = 0
    rw [grad_c6 hb hbs]; simp; field_simp; ring

end FdtrGauge.Moments
