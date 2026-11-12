import Mathlib.LinearAlgebra.Matrix.Determinant.Basic
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity

/-!
# The two-load kernel lemma (Theorem 1, Steps 1–3; Theorem 2, Proof A)

Campaign `20261111-fdtr-gauge-memory`, `proofs/B-analytic.md`, §4–§5.

On the `s`-plane (`p = s²`, `S = b_s s`), write `F = F(s)`, `F̃ = F(-s)` and `v = (C, X, Y)` for
the entries of the film quadrupole `M = [[C, X], [Y, C]]`, `C² − X Y = 1`.  The identities
`F (Y + C S) = C + X S` and `F̃ (Y − C S) = C − X S` read `𝖫 v = 0` with

    𝖫 = [[F S − 1, −S, F], [−F̃ S − 1, S, F̃]].

All statements are pointwise in `s` and hold over any field of characteristic zero (`ℂ` in the
application; the rank criterion is false in characteristic 2):

* `twoLoad_minor01/02/12`: the three `2 × 2` minors (CAS item 3 of `01_symbolic_identities`).
* `twoLoad_minors_eq_zero_iff`: for `S ≠ 0`, all minors vanish iff `F S = 1` and `F̃ S = −1`.
* `quadrupole_mem_ker`: the true `(C, X, Y)` is in the kernel.
* `ker_line`: if `𝖫` has rank 2, its kernel is a line.
* `step3_pointwise`: two kernel vectors with `C² − X Y = 1` satisfy `v' = λ v`, `λ² = 1`.
* `thm2_pointwise`: a kernel vector `u` with `2 C u_C − X u_Y − Y u_X = 0` vanishes.
* `F_sub_inv_S`, `F_eq_inv_S_iff`: `F − 1/S = X (S² − B²)/(S (Y + C S))`, `B² = Y/X` (Step 2).

The last section checks the algebraic endgame of Theorem 2, Proof A, as identities on `p > 0`
(`thm2_endgame`, `thm2_kernel_iff`, `thm2_D2_kernel_iff`).

Not formalised: the passage between pointwise and meromorphic statements (rank 2 off a discrete
set, constancy of `λ` on a connected set, `C → 1` as `s → 0`), and the derivation of the identities
in `p` of the last section from `D_c M ≡ 0` and from the factorised form of `D_c Z` on `D₂`.
-/

namespace FdtrGauge

open Matrix

section TwoLoad

variable {K : Type*} [Field K]

/-- The two-load matrix `𝖫` of Theorem 1, Step 1. -/
def twoLoad (F Ft S : K) : Matrix (Fin 2) (Fin 3) K :=
  !![F * S - 1, -S, F; -Ft * S - 1, S, Ft]

/-- The `2 × 2` minor of a `2 × 3` matrix on the columns `i, j`. -/
def minor2 (M : Matrix (Fin 2) (Fin 3) K) (i j : Fin 3) : K :=
  (M.submatrix id ![i, j]).det

lemma minor2_eq (M : Matrix (Fin 2) (Fin 3) K) (i j : Fin 3) :
    minor2 M i j = M 0 i * M 1 j - M 0 j * M 1 i := by
  simp [minor2, Matrix.det_fin_two]

theorem twoLoad_minor01 (F Ft S : K) :
    minor2 (twoLoad F Ft S) 0 1 = S * (S * (F - Ft) - 2) := by
  rw [minor2_eq]; simp [twoLoad]; ring

theorem twoLoad_minor02 (F Ft S : K) :
    minor2 (twoLoad F Ft S) 0 2 = 2 * F * Ft * S + F - Ft := by
  rw [minor2_eq]; simp [twoLoad]; ring

theorem twoLoad_minor12 (F Ft S : K) : minor2 (twoLoad F Ft S) 1 2 = -S * (F + Ft) := by
  rw [minor2_eq]; simp [twoLoad]; ring

/-- **Step 2 (rank criterion).** For `S ≠ 0` the three minors vanish simultaneously iff
`F S = 1` and `F̃ S = −1`, i.e. `F = 1/S`, `F̃ = −1/S`. -/
theorem twoLoad_minors_eq_zero_iff [CharZero K] {F Ft S : K} (hS : S ≠ 0) :
    (minor2 (twoLoad F Ft S) 0 1 = 0 ∧ minor2 (twoLoad F Ft S) 0 2 = 0 ∧
        minor2 (twoLoad F Ft S) 1 2 = 0) ↔ (F * S = 1 ∧ Ft * S = -1) := by
  rw [twoLoad_minor01, twoLoad_minor02, twoLoad_minor12]
  constructor
  · rintro ⟨h01, -, h12⟩
    have hsum : F + Ft = 0 := by
      have : S * (F + Ft) = 0 := by linear_combination -h12
      exact (mul_eq_zero.1 this).resolve_left hS
    have hdiff : S * (F - Ft) - 2 = 0 := (mul_eq_zero.1 h01).resolve_left hS
    exact ⟨by linear_combination (hdiff + S * hsum) / 2,
      by linear_combination (S * hsum - hdiff) / 2⟩
  · rintro ⟨h1, h2⟩
    have hsum : F + Ft = 0 := by
      have : (F + Ft) * S = 0 := by linear_combination h1 + h2
      exact (mul_eq_zero.1 this).resolve_right hS
    exact ⟨by linear_combination S * h1 - S * h2, by linear_combination 2 * F * h2 - hsum,
      by linear_combination -S * hsum⟩

lemma twoLoad_mulVec_apply_zero (F Ft S : K) (v : Fin 3 → K) :
    (twoLoad F Ft S *ᵥ v) 0 = (F * S - 1) * v 0 + -S * v 1 + F * v 2 := by
  simp [twoLoad, Matrix.mulVec, dotProduct, Fin.sum_univ_three]

lemma twoLoad_mulVec_apply_one (F Ft S : K) (v : Fin 3 → K) :
    (twoLoad F Ft S *ᵥ v) 1 = (-Ft * S - 1) * v 0 + S * v 1 + Ft * v 2 := by
  simp [twoLoad, Matrix.mulVec, dotProduct, Fin.sum_univ_three]

/-- The kernel condition `𝖫 (C, X, Y) = 0`, row by row. -/
theorem twoLoad_ker_iff (F Ft S C X Y : K) :
    twoLoad F Ft S *ᵥ ![C, X, Y] = 0 ↔
      ((F * S - 1) * C + -S * X + F * Y = 0 ∧ (-Ft * S - 1) * C + S * X + Ft * Y = 0) := by
  constructor
  · intro h
    have h0 := congrFun h 0
    have h1 := congrFun h 1
    rw [twoLoad_mulVec_apply_zero] at h0
    rw [twoLoad_mulVec_apply_one] at h1
    exact ⟨h0, h1⟩
  · rintro ⟨h0, h1⟩
    ext i
    fin_cases i
    · rw [Pi.zero_apply, ← h0]; exact twoLoad_mulVec_apply_zero F Ft S _
    · rw [Pi.zero_apply, ← h1]; exact twoLoad_mulVec_apply_one F Ft S _

/-- **Step 1.** The true `(C, X, Y)` solves `𝖫 v = 0` with `F = (C + X S)/(Y + C S)` and
`F̃ = (C − X S)/(Y − C S)` (the value at `−s`; `C, X, Y` are even in `s`). -/
theorem quadrupole_mem_ker {C X Y S : K} (hp : Y + C * S ≠ 0) (hm : Y - C * S ≠ 0) :
    twoLoad ((C + X * S) / (Y + C * S)) ((C - X * S) / (Y - C * S)) S *ᵥ ![C, X, Y] = 0 := by
  rw [twoLoad_ker_iff]
  constructor
  · field_simp
    ring
  · field_simp
    ring

/-- The components `(n₀, n₁, n₂) = (m₁₂, −m₀₂, m₀₁)` of the cross product of the rows. -/
def n0 (F Ft S : K) : K := -S * (F + Ft)
/-- See `n0`. -/
def n1 (F Ft S : K) : K := -(2 * F * Ft * S + F - Ft)
/-- See `n0`. -/
def n2 (F Ft S : K) : K := S * (S * (F - Ft) - 2)

/-- If `𝖫` has rank 2 (and `S ≠ 0`), the cross product of the rows is nonzero. -/
lemma cross_ne_zero [CharZero K] {F Ft S : K} (hS : S ≠ 0) (hrk : ¬(F * S = 1 ∧ Ft * S = -1)) :
    n0 F Ft S ≠ 0 ∨ n1 F Ft S ≠ 0 ∨ n2 F Ft S ≠ 0 := by
  by_contra hcon
  simp only [not_or, not_not] at hcon
  obtain ⟨h0, h1, h2⟩ := hcon
  apply hrk
  rw [← twoLoad_minors_eq_zero_iff hS, twoLoad_minor01, twoLoad_minor02, twoLoad_minor12]
  unfold n0 at h0
  unfold n1 at h1
  unfold n2 at h2
  exact ⟨h2, by linear_combination -h1, h0⟩

/-- A kernel vector is parallel to the cross product of the rows. -/
lemma ker_cross {F Ft S C X Y : K} (hv : twoLoad F Ft S *ᵥ ![C, X, Y] = 0) :
    n0 F Ft S * X = n1 F Ft S * C ∧ n0 F Ft S * Y = n2 F Ft S * C ∧
      n1 F Ft S * Y = n2 F Ft S * X := by
  obtain ⟨h1, h2⟩ := (twoLoad_ker_iff F Ft S C X Y).1 hv
  unfold n0 n1 n2
  exact ⟨by linear_combination Ft * h1 - F * h2, by linear_combination -(S * h1) - S * h2,
    by linear_combination (-Ft * S - 1) * h1 - (F * S - 1) * h2⟩

/-- **Kernel lemma.** If `S ≠ 0` and `𝖫` has rank 2, then any kernel vector `(C', X', Y')` is a
multiple of a nonzero kernel vector `(C, X, Y)`. -/
theorem ker_line [CharZero K] {F Ft S C X Y C' X' Y' : K} (hS : S ≠ 0)
    (hrk : ¬(F * S = 1 ∧ Ft * S = -1)) (hv : twoLoad F Ft S *ᵥ ![C, X, Y] = 0)
    (hw : twoLoad F Ft S *ᵥ ![C', X', Y'] = 0) (hv0 : ¬(C = 0 ∧ X = 0 ∧ Y = 0)) :
    ∃ μ : K, C' = μ * C ∧ X' = μ * X ∧ Y' = μ * Y := by
  obtain ⟨a01, a02, a12⟩ := ker_cross hv
  obtain ⟨b01, b02, b12⟩ := ker_cross hw
  set m0 := n0 F Ft S
  set m1 := n1 F Ft S
  set m2 := n2 F Ft S
  rcases cross_ne_zero hS hrk with h | h | h
  · -- `m0 ≠ 0`: `X = m1 C/m0`, `Y = m2 C/m0`, so `C ≠ 0`
    have hC : C ≠ 0 := by
      intro hC0
      apply hv0
      refine ⟨hC0, ?_, ?_⟩
      · have : m0 * X = 0 := by rw [a01, hC0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
      · have : m0 * Y = 0 := by rw [a02, hC0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
    refine ⟨C' / C, by field_simp, ?_, ?_⟩
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination C * b01 - C' * a01
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination C * b02 - C' * a02
  · -- `m1 ≠ 0`: `C = m0 X/m1`, `Y = m2 X/m1`, so `X ≠ 0`
    have hX : X ≠ 0 := by
      intro hX0
      apply hv0
      refine ⟨?_, hX0, ?_⟩
      · have : m1 * C = 0 := by rw [← a01, hX0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
      · have : m1 * Y = 0 := by rw [a12, hX0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
    refine ⟨X' / X, ?_, by field_simp, ?_⟩
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination X' * a01 - X * b01
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination X * b12 - X' * a12
  · -- `m2 ≠ 0`: `C = m0 Y/m2`, `X = m1 Y/m2`, so `Y ≠ 0`
    have hY : Y ≠ 0 := by
      intro hY0
      apply hv0
      refine ⟨?_, ?_, hY0⟩
      · have : m2 * C = 0 := by rw [← a02, hY0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
      · have : m2 * X = 0 := by rw [← a12, hY0, mul_zero]
        exact (mul_eq_zero.1 this).resolve_left h
    refine ⟨Y' / Y, ?_, ?_, by field_simp⟩
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination Y' * a02 - Y * b02
    · field_simp
      apply mul_left_cancel₀ h
      linear_combination Y' * a12 - Y * b12

/-- **Theorem 1, Step 3 (pointwise).** At a point where `𝖫` has rank 2, two kernel vectors with
`C² − X Y = 1` satisfy `(C', X', Y') = λ (C, X, Y)` with `λ² = 1`. -/
theorem step3_pointwise [CharZero K] {F Ft S C X Y C' X' Y' : K} (hS : S ≠ 0)
    (hrk : ¬(F * S = 1 ∧ Ft * S = -1)) (hv : twoLoad F Ft S *ᵥ ![C, X, Y] = 0)
    (hw : twoLoad F Ft S *ᵥ ![C', X', Y'] = 0) (hdet : C ^ 2 - X * Y = 1)
    (hdet' : C' ^ 2 - X' * Y' = 1) :
    ∃ μ : K, C' = μ * C ∧ X' = μ * X ∧ Y' = μ * Y ∧ μ ^ 2 = 1 := by
  have hv0 : ¬(C = 0 ∧ X = 0 ∧ Y = 0) := by
    rintro ⟨rfl, rfl, rfl⟩
    norm_num at hdet
  obtain ⟨μ, hC, hX, hY⟩ := ker_line hS hrk hv hw hv0
  refine ⟨μ, hC, hX, hY, ?_⟩
  rw [hC, hX, hY] at hdet'
  linear_combination hdet' - μ ^ 2 * hdet

/-- **Theorem 2, Proof A (pointwise).** If `𝖫 u = 0`, `𝖫 v = 0`, `C² − X Y = 1` and `u` satisfies
the linearised constraint `2 C u_C − X u_Y − Y u_X = 0`, then `u = 0` (here `u = D_c v`). -/
theorem thm2_pointwise [CharZero K] {F Ft S C X Y uC uX uY : K} (hS : S ≠ 0)
    (hrk : ¬(F * S = 1 ∧ Ft * S = -1)) (hv : twoLoad F Ft S *ᵥ ![C, X, Y] = 0)
    (hu : twoLoad F Ft S *ᵥ ![uC, uX, uY] = 0) (hdet : C ^ 2 - X * Y = 1)
    (hlin : 2 * C * uC - X * uY - Y * uX = 0) : uC = 0 ∧ uX = 0 ∧ uY = 0 := by
  have hv0 : ¬(C = 0 ∧ X = 0 ∧ Y = 0) := by
    rintro ⟨rfl, rfl, rfl⟩
    norm_num at hdet
  obtain ⟨μ, hC, hX, hY⟩ := ker_line hS hrk hv hu hv0
  have hμ : μ = 0 := by
    have h2 : 2 * μ = 0 := by
      rw [hC, hX, hY] at hlin
      linear_combination hlin - 2 * μ * hdet
    simpa using h2
  exact ⟨by rw [hC, hμ, zero_mul], by rw [hX, hμ, zero_mul], by rw [hY, hμ, zero_mul]⟩

/-- **Theorem 1, Step 2.** `F − 1/S = X (S² − B²)/(S (Y + C S))` with `B² = Y/X`. -/
theorem F_sub_inv_S {C X Y S : K} (hS : S ≠ 0) (hX : X ≠ 0) (hD : Y + C * S ≠ 0) :
    (C + X * S) / (Y + C * S) - 1 / S = X * (S ^ 2 - Y / X) / (S * (Y + C * S)) := by
  field_simp
  ring

/-- Consequently `F = 1/S` at a point iff `S² = B²` there. -/
theorem F_eq_inv_S_iff {C X Y S : K} (hS : S ≠ 0) (hX : X ≠ 0) (hD : Y + C * S ≠ 0) :
    (C + X * S) / (Y + C * S) = 1 / S ↔ S ^ 2 = Y / X := by
  rw [← sub_eq_zero, F_sub_inv_S hS hX hD, div_eq_zero_iff, mul_eq_zero, sub_eq_zero]
  constructor
  · rintro ((h | h) | h)
    · exact absurd h hX
    · exact h
    · exact absurd h (mul_ne_zero hS hD)
  · intro h
    exact Or.inl (Or.inr h)

end TwoLoad

/-! ## Theorem 2, Proof A: the algebraic endgame on `p > 0` -/

/-- A real quadratic vanishing at `p = 1, 2, 3` is zero. -/
lemma quad_eq_zero_of_three {A B C : ℝ} (h1 : A + B * 1 + C * 1 ^ 2 = 0)
    (h2 : A + B * 2 + C * 2 ^ 2 = 0) (h3 : A + B * 3 + C * 3 ^ 2 = 0) :
    A = 0 ∧ B = 0 ∧ C = 0 :=
  ⟨by linear_combination 3 * h1 - 3 * h2 + h3,
    by linear_combination (-5 * h1 + 8 * h2 - 3 * h3) / 2,
    by linear_combination (h1 - 2 * h2 + h3) / 2⟩

/-- `D_c log g = c₃ p/(1 + τ_R p) − c₄ p/(1 + τ_ℓ p)` for `c = (c₁, c₂, c₃, c₄)`. -/
noncomputable def dlogg (tR tl c3 c4 p : ℝ) : ℝ := c3 * p / (1 + tR * p) - c4 * p / (1 + tl * p)

/-- **Theorem 2, Proof A, case `θ ∉ D₂` (endgame).** The identities `D_c log z ≡ 0` and
`D_c log B² ≡ 0` on `p > 0`, i.e. `2 c₁/ξ₁ + D_c log g = 0` and `2 c₂/b − D_c log g = 0`, force
`c₁ = c₂ = 0`, `c₃ = c₄` and `c₃ (τ_ℓ − τ_R) = 0`. -/
theorem thm2_endgame {xi b tR tl c1 c2 c3 c4 : ℝ} (hxi : xi ≠ 0) (hb : b ≠ 0) (hR : 0 ≤ tR)
    (hl : 0 ≤ tl)
    (hz : ∀ p : ℝ, 0 < p → 2 * c1 / xi + dlogg tR tl c3 c4 p = 0)
    (hB : ∀ p : ℝ, 0 < p → 2 * c2 / b - dlogg tR tl c3 c4 p = 0) :
    c1 = 0 ∧ c2 = 0 ∧ c3 = c4 ∧ c3 * (tl - tR) = 0 := by
  -- cross-multiplied form of `hz` (multiplied by `ξ₁`)
  have key : ∀ p : ℝ, 0 < p →
      2 * c1 + (2 * c1 * (tR + tl) + xi * (c3 - c4)) * p +
        (2 * c1 * (tR * tl) + xi * (c3 * tl - c4 * tR)) * p ^ 2 = 0 := by
    intro p hp
    have h := hz p hp
    have d1 : 1 + tR * p ≠ 0 := by positivity
    have d2 : 1 + tl * p ≠ 0 := by positivity
    unfold dlogg at h
    field_simp at h
    linear_combination h
  obtain ⟨hA, hB', hC⟩ := quad_eq_zero_of_three (key 1 one_pos) (key 2 two_pos)
    (key 3 (by norm_num))
  have hc1 : c1 = 0 := by linarith
  rw [hc1] at hB' hC
  have hc34 : c3 = c4 := by
    have : xi * (c3 - c4) = 0 := by linear_combination hB'
    exact sub_eq_zero.1 ((mul_eq_zero.1 this).resolve_left hxi)
  have hc3 : c3 * (tl - tR) = 0 := by
    have : xi * (c3 * tl - c4 * tR) = 0 := by linear_combination hC
    have h' := (mul_eq_zero.1 this).resolve_left hxi
    rw [← hc34] at h'
    linear_combination h'
  have hdl : dlogg tR tl c3 c4 1 = 0 := by
    have d1 : 1 + tR * 1 ≠ 0 := by positivity
    have d2 : 1 + tl * 1 ≠ 0 := by positivity
    unfold dlogg
    rw [← hc34]
    field_simp
    linear_combination hc3
  have hc2 : c2 = 0 := by
    have h := hB 1 one_pos
    rw [hdl, sub_zero, div_eq_zero_iff] at h
    rcases h with h | h
    · linarith
    · exact absurd h hb
  exact ⟨hc1, hc2, hc34, hc3⟩

/-- **Theorem 2 on `Θ ∖ D₂`: the kernel.** The set of `c = (c₁, c₂, c₃, c₄)` satisfying both
identities on `p > 0` is `{0}` off the resonance and `ℝ (0, 0, 1, 1)` on it. -/
theorem thm2_kernel_iff {xi b tR tl : ℝ} (hxi : xi ≠ 0) (hb : b ≠ 0) (hR : 0 ≤ tR) (hl : 0 ≤ tl)
    (c1 c2 c3 c4 : ℝ) :
    ((∀ p : ℝ, 0 < p → 2 * c1 / xi + dlogg tR tl c3 c4 p = 0) ∧
        (∀ p : ℝ, 0 < p → 2 * c2 / b - dlogg tR tl c3 c4 p = 0)) ↔
      (c1 = 0 ∧ c2 = 0 ∧ c3 = c4 ∧ (tR ≠ tl → c3 = 0)) := by
  constructor
  · rintro ⟨hz, hB⟩
    obtain ⟨h1, h2, h34, h3⟩ := thm2_endgame hxi hb hR hl hz hB
    refine ⟨h1, h2, h34, fun hne => ?_⟩
    exact (mul_eq_zero.1 h3).resolve_right (sub_ne_zero.2 (Ne.symm hne))
  · rintro ⟨rfl, rfl, rfl, h⟩
    have hdl : ∀ p : ℝ, dlogg tR tl c3 c3 p = 0 := by
      intro p
      by_cases hres : tR = tl
      · subst hres; simp [dlogg]
      · simp [dlogg, h hres]
    exact ⟨fun p _ => by simp [hdl p], fun p _ => by simp [hdl p]⟩

/-- **Theorem 2 on `D₂`: the kernel.** With
`D_c Z = −(1 − e^{−2A})/(2B) · (2 c₂/b − (c₃ − c₄) p/(1 + t p))` and the prefactor `≢ 0`, the
kernel is `{c₂ = 0, c₃ = c₄}` (`c₁` free). -/
theorem thm2_D2_kernel_iff {b t : ℝ} (hb : b ≠ 0) (ht : 0 ≤ t) (c2 c3 c4 : ℝ) :
    (∀ p : ℝ, 0 < p → 2 * c2 / b - (c3 - c4) * (p / (1 + t * p)) = 0) ↔ (c2 = 0 ∧ c3 = c4) := by
  constructor
  · intro h
    have e : ∀ p : ℝ, 0 < p → 2 * c2 + (2 * c2 * t - b * (c3 - c4)) * p = 0 := by
      intro p hp
      have h' := h p hp
      have d : 1 + t * p ≠ 0 := by positivity
      field_simp at h'
      linear_combination h'
    have e1 := e 1 one_pos
    have e2 := e 2 two_pos
    have hβ : 2 * c2 * t - b * (c3 - c4) = 0 := by linear_combination e2 - e1
    have hc2 : c2 = 0 := by linarith
    rw [hc2] at hβ
    have : b * (c3 - c4) = 0 := by linarith
    exact ⟨hc2, sub_eq_zero.1 ((mul_eq_zero.1 this).resolve_left hb)⟩
  · rintro ⟨rfl, rfl⟩ p _
    simp

end FdtrGauge
