import Mathlib.Analysis.Calculus.Deriv.Add
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Inv
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.Ring

/-!
# Low-frequency Laurent coefficients of the front-face response and their gradients

Campaign `20261111-fdtr-gauge-memory`, `proofs/B-analytic.md`, §5 (Theorem 2, Proof C) and §8
(Theorem N2).  Perfect contact (`R = 0`).

With `p = s²` and `F(s) = Z(s²)`, `F(s) = 1/(b_s s) + Σ_{k ≥ 0} c_k s^k` near `s = 0`.
**Here the coefficients `c₀, …, c₆` are taken as given rational functions of
`θ = (ξ₁, b, τ_R, τ_ℓ)`** (the closed forms of `derivations/B-analytic/07_moment_jacobian.out`,
which agree with `laurent.py` and with the forms displayed in Proof C).  They are not derived from
`Z` in this file; `RiccatiHierarchy.lean` checks that they solve the Riccati hierarchy in `ξ₁`
that the Laurent coefficients of `Z` obey.

`grad f θ` is the vector of the four partial derivatives of `f` at `θ`, each defined with `deriv`
of the restriction of `f` to one variable.  The lemmas `grad_c0`, …, `grad_c6` evaluate them for
`b ≠ 0`, `b_s ≠ 0`.  The coefficients are polynomial in `ξ₁, τ_R, τ_ℓ`, so at boundary points
(`τ = 0`) these are the two-sided derivatives of the polynomial extension, which is how Proof C
uses them.
-/

set_option linter.unusedSimpArgs false
set_option linter.style.longLine false
set_option linter.unusedVariables false

namespace FdtrGauge.Moments

/-- Derivative evaluation for the explicit rational expressions of this file: differentiate with
the lambda-form `deriv` lemmas (side conditions by `fun_prop` or from the hypotheses), then clear
denominators. -/
macro "deriv_tac" : tactic => `(tactic| (
  simp (disch := first | fun_prop | simp [*]) only [deriv_div_const, deriv_fun_div, deriv_fun_mul,
    deriv_neg, deriv_neg'', deriv_fun_add, deriv_fun_sub, deriv_id'', deriv_const', deriv_fun_pow,
    deriv_const_mul_id, deriv_const_mul_field, deriv_mul_const_field] <;>
  (try field_simp) <;>
  ring))

/-! ## The coefficients (07_moment_jacobian.out, verbatim) -/

/-- The Laurent coefficient `c_0` (perfect contact). -/
noncomputable def c0 (bs xi b tR tl : ℝ) : ℝ :=
  -xi*(b - bs)*(b + bs)/(b*bs^2)

/-- The Laurent coefficient `c_1` (perfect contact). -/
noncomputable def c1 (bs xi b tR tl : ℝ) : ℝ :=
  xi^2*(b - bs)*(b + bs)/bs^3

/-- The Laurent coefficient `c_2` (perfect contact). -/
noncomputable def c2 (bs xi b tR tl : ℝ) : ℝ :=
  -xi*(3*b^4*xi^2 - 4*b^2*bs^2*xi^2 - 3*bs^4*tR + 3*bs^4*tl + bs^4*xi^2)/(3*b*bs^4)

/-- The Laurent coefficient `c_3` (perfect contact). -/
noncomputable def c3 (bs xi b tR tl : ℝ) : ℝ :=
  xi^2*(3*b^4*xi^2 - 5*b^2*bs^2*xi^2 - 3*bs^4*tR + 3*bs^4*tl + 2*bs^4*xi^2)/(3*bs^5)

/-- The Laurent coefficient `c_4` (perfect contact). -/
noncomputable def c4 (bs xi b tR tl : ℝ) : ℝ :=
  -xi*(15*b^6*xi^4 - 30*b^4*bs^2*xi^4 - 20*b^2*bs^4*tR*xi^2 + 20*b^2*bs^4*tl*xi^2 +
    17*b^2*bs^4*xi^4 + 15*bs^6*tR*tl + 10*bs^6*tR*xi^2 - 15*bs^6*tl^2 - 10*bs^6*tl*xi^2 -
    2*bs^6*xi^4)/(15*b*bs^6)

/-- The Laurent coefficient `c_5` (perfect contact). -/
noncomputable def c5 (bs xi b tR tl : ℝ) : ℝ :=
  xi^2*(45*b^6*xi^4 - 105*b^4*bs^2*xi^4 - 75*b^2*bs^4*tR*xi^2 + 75*b^2*bs^4*tl*xi^2 +
    77*b^2*bs^4*xi^4 + 45*bs^6*tR*tl + 60*bs^6*tR*xi^2 - 45*bs^6*tl^2 - 60*bs^6*tl*xi^2 -
    17*bs^6*xi^4)/(45*bs^7)

/-- The Laurent coefficient `c_6` (perfect contact). -/
noncomputable def c6 (bs xi b tR tl : ℝ) : ℝ :=
  -xi*(315*b^8*xi^6 - 840*b^6*bs^2*xi^6 - 630*b^4*bs^4*tR*xi^4 + 630*b^4*bs^4*tl*xi^4 +
    756*b^4*bs^4*xi^6 + 420*b^2*bs^6*tR*tl*xi^2 + 714*b^2*bs^6*tR*xi^4 -
    420*b^2*bs^6*tl^2*xi^2 - 714*b^2*bs^6*tl*xi^4 - 248*b^2*bs^6*xi^6 + 105*bs^8*tR^2*xi^2 -
    315*bs^8*tR*tl^2 - 420*bs^8*tR*tl*xi^2 - 126*bs^8*tR*xi^4 + 315*bs^8*tl^3 +
    315*bs^8*tl^2*xi^2 + 126*bs^8*tl*xi^4 + 17*bs^8*xi^6)/(315*b*bs^8)

/-- The family `k ↦ c_k`, `k = 0, …, 6`. -/
noncomputable def coeff (bs : ℝ) : Fin 7 → ℝ → ℝ → ℝ → ℝ → ℝ :=
  ![c0 bs, c1 bs, c2 bs, c3 bs, c4 bs, c5 bs, c6 bs]

/-- `∇f(θ) = (∂_{ξ₁} f, ∂_b f, ∂_{τ_R} f, ∂_{τ_ℓ} f)(θ)`. -/
noncomputable def grad (f : ℝ → ℝ → ℝ → ℝ → ℝ) (xi b tR tl : ℝ) : Fin 4 → ℝ :=
  ![deriv (fun y => f y b tR tl) xi, deriv (fun y => f xi y tR tl) b,
    deriv (fun y => f xi b y tl) tR, deriv (fun y => f xi b tR y) tl]

/-! ## The forms displayed in Proof C (`proofs/B-analytic.md`, §5) -/

theorem c0_eq_proofC {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    c0 bs xi b tR tl =
      xi * (bs ^ 2 - b ^ 2) / (b * bs ^ 2) := by
  unfold c0
  field_simp
  ring

theorem c1_eq_proofC {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    c1 bs xi b tR tl =
      xi ^ 2 * (b ^ 2 - bs ^ 2) / bs ^ 3 := by
  unfold c1
  field_simp
  ring

theorem c2_eq_proofC {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    c2 bs xi b tR tl =
      xi * (tR - tl) / b - xi ^ 3 * (3 * b ^ 4 - 4 * b ^ 2 * bs ^ 2 + bs ^ 4) / (3 * b * bs ^ 4) := by
  unfold c2
  field_simp
  ring

theorem c3_eq_proofC {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    c3 bs xi b tR tl =
      xi ^ 2 * (3 * bs ^ 4 * (tl - tR) + (3 * b ^ 4 - 5 * b ^ 2 * bs ^ 2 + 2 * bs ^ 4) * xi ^ 2) / (3 * bs ^ 5) := by
  unfold c3
  field_simp
  ring

/-! ## Partial derivatives -/

theorem deriv_c0_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c0 bs y b tR tl) xi =
      -(b - bs)*(b + bs)/(b*bs^2) := by
  simp only [c0]
  deriv_tac

theorem deriv_c0_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c0 bs xi y tR tl) b =
      -xi*(b^2 + bs^2)/(b^2*bs^2) := by
  simp only [c0]
  deriv_tac

theorem deriv_c0_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c0 bs xi b y tl) tR =
      0 := by
  simp only [c0]
  deriv_tac

theorem deriv_c0_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c0 bs xi b tR y) tl =
      0 := by
  simp only [c0]
  deriv_tac

theorem deriv_c1_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c1 bs y b tR tl) xi =
      2*xi*(b - bs)*(b + bs)/bs^3 := by
  simp only [c1]
  deriv_tac

theorem deriv_c1_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c1 bs xi y tR tl) b =
      2*b*xi^2/bs^3 := by
  simp only [c1]
  deriv_tac

theorem deriv_c1_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c1 bs xi b y tl) tR =
      0 := by
  simp only [c1]
  deriv_tac

theorem deriv_c1_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c1 bs xi b tR y) tl =
      0 := by
  simp only [c1]
  deriv_tac

theorem deriv_c2_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c2 bs y b tR tl) xi =
      -(3*b^4*xi^2 - 4*b^2*bs^2*xi^2 - bs^4*tR + bs^4*tl + bs^4*xi^2)/(b*bs^4) := by
  simp only [c2]
  deriv_tac

theorem deriv_c2_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c2 bs xi y tR tl) b =
      -xi*(9*b^4*xi^2 - 4*b^2*bs^2*xi^2 + 3*bs^4*tR - 3*bs^4*tl - bs^4*xi^2)/(3*b^2*bs^4) := by
  simp only [c2]
  deriv_tac

theorem deriv_c2_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c2 bs xi b y tl) tR =
      xi/b := by
  simp only [c2]
  deriv_tac

theorem deriv_c2_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c2 bs xi b tR y) tl =
      -xi/b := by
  simp only [c2]
  deriv_tac

theorem deriv_c3_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c3 bs y b tR tl) xi =
      2*xi*(6*b^4*xi^2 - 10*b^2*bs^2*xi^2 - 3*bs^4*tR + 3*bs^4*tl + 4*bs^4*xi^2)/(3*bs^5) := by
  simp only [c3]
  deriv_tac

theorem deriv_c3_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c3 bs xi y tR tl) b =
      2*b*xi^4*(6*b^2 - 5*bs^2)/(3*bs^5) := by
  simp only [c3]
  deriv_tac

theorem deriv_c3_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c3 bs xi b y tl) tR =
      -xi^2/bs := by
  simp only [c3]
  deriv_tac

theorem deriv_c3_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c3 bs xi b tR y) tl =
      xi^2/bs := by
  simp only [c3]
  deriv_tac

theorem deriv_c4_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c4 bs y b tR tl) xi =
      -(15*b^6*xi^4 - 30*b^4*bs^2*xi^4 - 12*b^2*bs^4*tR*xi^2 + 12*b^2*bs^4*tl*xi^2 +
        17*b^2*bs^4*xi^4 + 3*bs^6*tR*tl + 6*bs^6*tR*xi^2 - 3*bs^6*tl^2 - 6*bs^6*tl*xi^2 -
        2*bs^6*xi^4)/(3*b*bs^6) := by
  simp only [c4]
  deriv_tac

theorem deriv_c4_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c4 bs xi y tR tl) b =
      -xi*(75*b^6*xi^4 - 90*b^4*bs^2*xi^4 - 20*b^2*bs^4*tR*xi^2 + 20*b^2*bs^4*tl*xi^2 +
        17*b^2*bs^4*xi^4 - 15*bs^6*tR*tl - 10*bs^6*tR*xi^2 + 15*bs^6*tl^2 + 10*bs^6*tl*xi^2 +
        2*bs^6*xi^4)/(15*b^2*bs^6) := by
  simp only [c4]
  deriv_tac

theorem deriv_c4_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c4 bs xi b y tl) tR =
      xi*(4*b^2*xi^2 - 3*bs^2*tl - 2*bs^2*xi^2)/(3*b*bs^2) := by
  simp only [c4]
  deriv_tac

theorem deriv_c4_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c4 bs xi b tR y) tl =
      -xi*(4*b^2*xi^2 + 3*bs^2*tR - 6*bs^2*tl - 2*bs^2*xi^2)/(3*b*bs^2) := by
  simp only [c4]
  deriv_tac

theorem deriv_c5_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c5 bs y b tR tl) xi =
      2*xi*(45*b^6*xi^4 - 105*b^4*bs^2*xi^4 - 50*b^2*bs^4*tR*xi^2 + 50*b^2*bs^4*tl*xi^2 +
        77*b^2*bs^4*xi^4 + 15*bs^6*tR*tl + 40*bs^6*tR*xi^2 - 15*bs^6*tl^2 - 40*bs^6*tl*xi^2 -
        17*bs^6*xi^4)/(15*bs^7) := by
  simp only [c5]
  deriv_tac

theorem deriv_c5_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c5 bs xi y tR tl) b =
      2*b*xi^4*(135*b^4*xi^2 - 210*b^2*bs^2*xi^2 - 75*bs^4*tR + 75*bs^4*tl +
        77*bs^4*xi^2)/(45*bs^7) := by
  simp only [c5]
  deriv_tac

theorem deriv_c5_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c5 bs xi b y tl) tR =
      -xi^2*(5*b^2*xi^2 - 3*bs^2*tl - 4*bs^2*xi^2)/(3*bs^3) := by
  simp only [c5]
  deriv_tac

theorem deriv_c5_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c5 bs xi b tR y) tl =
      xi^2*(5*b^2*xi^2 + 3*bs^2*tR - 6*bs^2*tl - 4*bs^2*xi^2)/(3*bs^3) := by
  simp only [c5]
  deriv_tac

theorem deriv_c6_xi {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c6 bs y b tR tl) xi =
      -(315*b^8*xi^6 - 840*b^6*bs^2*xi^6 - 450*b^4*bs^4*tR*xi^4 + 450*b^4*bs^4*tl*xi^4 +
        756*b^4*bs^4*xi^6 + 180*b^2*bs^6*tR*tl*xi^2 + 510*b^2*bs^6*tR*xi^4 -
        180*b^2*bs^6*tl^2*xi^2 - 510*b^2*bs^6*tl*xi^4 - 248*b^2*bs^6*xi^6 + 45*bs^8*tR^2*xi^2 -
        45*bs^8*tR*tl^2 - 180*bs^8*tR*tl*xi^2 - 90*bs^8*tR*xi^4 + 45*bs^8*tl^3 +
        135*bs^8*tl^2*xi^2 + 90*bs^8*tl*xi^4 + 17*bs^8*xi^6)/(45*b*bs^8) := by
  simp only [c6]
  deriv_tac

theorem deriv_c6_b {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c6 bs xi y tR tl) b =
      -xi*(2205*b^8*xi^6 - 4200*b^6*bs^2*xi^6 - 1890*b^4*bs^4*tR*xi^4 + 1890*b^4*bs^4*tl*xi^4 +
        2268*b^4*bs^4*xi^6 + 420*b^2*bs^6*tR*tl*xi^2 + 714*b^2*bs^6*tR*xi^4 -
        420*b^2*bs^6*tl^2*xi^2 - 714*b^2*bs^6*tl*xi^4 - 248*b^2*bs^6*xi^6 -
        105*bs^8*tR^2*xi^2 + 315*bs^8*tR*tl^2 + 420*bs^8*tR*tl*xi^2 + 126*bs^8*tR*xi^4 -
        315*bs^8*tl^3 - 315*bs^8*tl^2*xi^2 - 126*bs^8*tl*xi^4 - 17*bs^8*xi^6)/(315*b^2*bs^8) := by
  simp only [c6]
  deriv_tac

theorem deriv_c6_tR {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c6 bs xi b y tl) tR =
      xi*(30*b^4*xi^4 - 20*b^2*bs^2*tl*xi^2 - 34*b^2*bs^2*xi^4 - 10*bs^4*tR*xi^2 +
        15*bs^4*tl^2 + 20*bs^4*tl*xi^2 + 6*bs^4*xi^4)/(15*b*bs^4) := by
  simp only [c6]
  deriv_tac

theorem deriv_c6_tl {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    deriv (fun y => c6 bs xi b tR y) tl =
      -xi*(30*b^4*xi^4 + 20*b^2*bs^2*tR*xi^2 - 40*b^2*bs^2*tl*xi^2 - 34*b^2*bs^2*xi^4 -
        30*bs^4*tR*tl - 20*bs^4*tR*xi^2 + 45*bs^4*tl^2 + 30*bs^4*tl*xi^2 +
        6*bs^4*xi^4)/(15*b*bs^4) := by
  simp only [c6]
  deriv_tac

/-- `∇c_0`. -/
theorem grad_c0 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c0 bs) xi b tR tl =
      ![-(b - bs)*(b + bs)/(b*bs^2),
        -xi*(b^2 + bs^2)/(b^2*bs^2),
        0,
        0] := by
  rw [grad, deriv_c0_xi hb hbs, deriv_c0_b hb hbs, deriv_c0_tR hb hbs,
    deriv_c0_tl hb hbs]

/-- `∇c_1`. -/
theorem grad_c1 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c1 bs) xi b tR tl =
      ![2*xi*(b - bs)*(b + bs)/bs^3,
        2*b*xi^2/bs^3,
        0,
        0] := by
  rw [grad, deriv_c1_xi hb hbs, deriv_c1_b hb hbs, deriv_c1_tR hb hbs,
    deriv_c1_tl hb hbs]

/-- `∇c_2`. -/
theorem grad_c2 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c2 bs) xi b tR tl =
      ![-(3*b^4*xi^2 - 4*b^2*bs^2*xi^2 - bs^4*tR + bs^4*tl + bs^4*xi^2)/(b*bs^4),
        -xi*(9*b^4*xi^2 - 4*b^2*bs^2*xi^2 + 3*bs^4*tR - 3*bs^4*tl - bs^4*xi^2)/(3*b^2*bs^4),
        xi/b,
        -xi/b] := by
  rw [grad, deriv_c2_xi hb hbs, deriv_c2_b hb hbs, deriv_c2_tR hb hbs,
    deriv_c2_tl hb hbs]

/-- `∇c_3`. -/
theorem grad_c3 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c3 bs) xi b tR tl =
      ![2*xi*(6*b^4*xi^2 - 10*b^2*bs^2*xi^2 - 3*bs^4*tR + 3*bs^4*tl + 4*bs^4*xi^2)/(3*bs^5),
        2*b*xi^4*(6*b^2 - 5*bs^2)/(3*bs^5),
        -xi^2/bs,
        xi^2/bs] := by
  rw [grad, deriv_c3_xi hb hbs, deriv_c3_b hb hbs, deriv_c3_tR hb hbs,
    deriv_c3_tl hb hbs]

/-- `∇c_4`. -/
theorem grad_c4 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c4 bs) xi b tR tl =
      ![-(15*b^6*xi^4 - 30*b^4*bs^2*xi^4 - 12*b^2*bs^4*tR*xi^2 + 12*b^2*bs^4*tl*xi^2 + 17*b^2*bs^4*xi^4 + 3*bs^6*tR*tl + 6*bs^6*tR*xi^2 - 3*bs^6*tl^2 - 6*bs^6*tl*xi^2 - 2*bs^6*xi^4)/(3*b*bs^6),
        -xi*(75*b^6*xi^4 - 90*b^4*bs^2*xi^4 - 20*b^2*bs^4*tR*xi^2 + 20*b^2*bs^4*tl*xi^2 + 17*b^2*bs^4*xi^4 - 15*bs^6*tR*tl - 10*bs^6*tR*xi^2 + 15*bs^6*tl^2 + 10*bs^6*tl*xi^2 + 2*bs^6*xi^4)/(15*b^2*bs^6),
        xi*(4*b^2*xi^2 - 3*bs^2*tl - 2*bs^2*xi^2)/(3*b*bs^2),
        -xi*(4*b^2*xi^2 + 3*bs^2*tR - 6*bs^2*tl - 2*bs^2*xi^2)/(3*b*bs^2)] := by
  rw [grad, deriv_c4_xi hb hbs, deriv_c4_b hb hbs, deriv_c4_tR hb hbs,
    deriv_c4_tl hb hbs]

/-- `∇c_5`. -/
theorem grad_c5 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c5 bs) xi b tR tl =
      ![2*xi*(45*b^6*xi^4 - 105*b^4*bs^2*xi^4 - 50*b^2*bs^4*tR*xi^2 + 50*b^2*bs^4*tl*xi^2 + 77*b^2*bs^4*xi^4 + 15*bs^6*tR*tl + 40*bs^6*tR*xi^2 - 15*bs^6*tl^2 - 40*bs^6*tl*xi^2 - 17*bs^6*xi^4)/(15*bs^7),
        2*b*xi^4*(135*b^4*xi^2 - 210*b^2*bs^2*xi^2 - 75*bs^4*tR + 75*bs^4*tl + 77*bs^4*xi^2)/(45*bs^7),
        -xi^2*(5*b^2*xi^2 - 3*bs^2*tl - 4*bs^2*xi^2)/(3*bs^3),
        xi^2*(5*b^2*xi^2 + 3*bs^2*tR - 6*bs^2*tl - 4*bs^2*xi^2)/(3*bs^3)] := by
  rw [grad, deriv_c5_xi hb hbs, deriv_c5_b hb hbs, deriv_c5_tR hb hbs,
    deriv_c5_tl hb hbs]

/-- `∇c_6`. -/
theorem grad_c6 {bs xi b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) :
    grad (c6 bs) xi b tR tl =
      ![-(315*b^8*xi^6 - 840*b^6*bs^2*xi^6 - 450*b^4*bs^4*tR*xi^4 + 450*b^4*bs^4*tl*xi^4 + 756*b^4*bs^4*xi^6 + 180*b^2*bs^6*tR*tl*xi^2 + 510*b^2*bs^6*tR*xi^4 - 180*b^2*bs^6*tl^2*xi^2 - 510*b^2*bs^6*tl*xi^4 - 248*b^2*bs^6*xi^6 + 45*bs^8*tR^2*xi^2 - 45*bs^8*tR*tl^2 - 180*bs^8*tR*tl*xi^2 - 90*bs^8*tR*xi^4 + 45*bs^8*tl^3 + 135*bs^8*tl^2*xi^2 + 90*bs^8*tl*xi^4 + 17*bs^8*xi^6)/(45*b*bs^8),
        -xi*(2205*b^8*xi^6 - 4200*b^6*bs^2*xi^6 - 1890*b^4*bs^4*tR*xi^4 + 1890*b^4*bs^4*tl*xi^4 + 2268*b^4*bs^4*xi^6 + 420*b^2*bs^6*tR*tl*xi^2 + 714*b^2*bs^6*tR*xi^4 - 420*b^2*bs^6*tl^2*xi^2 - 714*b^2*bs^6*tl*xi^4 - 248*b^2*bs^6*xi^6 - 105*bs^8*tR^2*xi^2 + 315*bs^8*tR*tl^2 + 420*bs^8*tR*tl*xi^2 + 126*bs^8*tR*xi^4 - 315*bs^8*tl^3 - 315*bs^8*tl^2*xi^2 - 126*bs^8*tl*xi^4 - 17*bs^8*xi^6)/(315*b^2*bs^8),
        xi*(30*b^4*xi^4 - 20*b^2*bs^2*tl*xi^2 - 34*b^2*bs^2*xi^4 - 10*bs^4*tR*xi^2 + 15*bs^4*tl^2 + 20*bs^4*tl*xi^2 + 6*bs^4*xi^4)/(15*b*bs^4),
        -xi*(30*b^4*xi^4 + 20*b^2*bs^2*tR*xi^2 - 40*b^2*bs^2*tl*xi^2 - 34*b^2*bs^2*xi^4 - 30*bs^4*tR*tl - 20*bs^4*tR*xi^2 + 45*bs^4*tl^2 + 30*bs^4*tl*xi^2 + 6*bs^4*xi^4)/(15*b*bs^4)] := by
  rw [grad, deriv_c6_xi hb hbs, deriv_c6_b hb hbs, deriv_c6_tR hb hbs,
    deriv_c6_tl hb hbs]

end FdtrGauge.Moments
