import ResearchProofs.FdtrGauge.MomentCoefficients
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.LinearCombination

/-!
# The Riccati hierarchy: an independent characterisation of the Laurent coefficients

Campaign `20261111-fdtr-gauge-memory`; supports `proofs/B-analytic.md`, §5 (Proof C) and §8.
This route is not used in `B-analytic.md`; it checks the coefficient formulas there by a method
independent of the series division of `laurent.py` and of the CAS series of
`07_moment_jacobian.py` / `07w_wolfram_crosscheck.wl`.

1. `hasDerivAt_Zfilm`, `hasDerivAt_Zfilm_model`: the impedance seen through a uniform film,
   `Z(x) = (L cosh κx + sinh κx / B)/(B L sinh κx + cosh κx)` (Lemma 1 with `A = κ x`), satisfies
   `Z' = (κ/B)(1 − B² Z²)`, `Z(0) = L`.  With Lemma 0 (`κ B = b p`, `κ/B = g/b`) this is
   `∂_{ξ₁} Z = g/b − b p Z²`.
2. With `p = s²` and `s F(s) = Σ_k q_k s^k` (`q₀ = 1/b_s`, `q_{k+1} = c_k`), comparing powers
   of `s` gives the triangular hierarchy
       `∂_{ξ₁} q_{n+1} = [n even] γ_{n/2}/b − b Σ_{j=0}^{n} q_j q_{n−j}`,
   `γ_m` the Taylor coefficients of `g(s²)` in `s²` (`gam_series`), with initial data
   `q_{n+1}(0) = 0` (perfect contact, `F|_{ξ₁=0} = 1/(b_s s)`).
   `hierarchy`, `hierarchy_init`: the coefficients `c₀, …, c₆` of `MomentCoefficients.lean`
   satisfy it; `hierarchy_unique`: it has no other solution.

Not formalised (step 2, first half): that the Laurent coefficients of `F` are differentiable in
`ξ₁` and obey the hierarchy, i.e. termwise differentiation of the Laurent series in `ξ₁` and the
identification of coefficients.  Under that standard analytic input the coefficients used in the
minor identities are the Laurent coefficients of the model.
-/

set_option linter.unusedSimpArgs false

namespace FdtrGauge.Riccati

open FdtrGauge.Moments

/-! ## 1. The Riccati equation of the front-face impedance -/

/-- Impedance through a uniform film with `A = κ x`, film admittance factor `B` and load `L`. -/
noncomputable def Zfilm (κ B L x : ℂ) : ℂ :=
  (L * Complex.cosh (κ * x) + Complex.sinh (κ * x) / B) /
    (B * L * Complex.sinh (κ * x) + Complex.cosh (κ * x))

theorem Zfilm_zero (κ B L : ℂ) : Zfilm κ B L 0 = L := by
  simp [Zfilm]

/-- **Riccati equation in the thickness variable:** `Z' = (κ/B)(1 − B² Z²)`. -/
theorem hasDerivAt_Zfilm {κ B L x : ℂ} (hB : B ≠ 0)
    (hD : B * L * Complex.sinh (κ * x) + Complex.cosh (κ * x) ≠ 0) :
    HasDerivAt (Zfilm κ B L) (κ / B * (1 - B ^ 2 * Zfilm κ B L x ^ 2)) x := by
  have hA : HasDerivAt (fun y : ℂ => κ * y) κ x := by
    simpa using (hasDerivAt_id x).const_mul κ
  have hc : HasDerivAt (fun y => Complex.cosh (κ * y)) (Complex.sinh (κ * x) * κ) x :=
    (Complex.hasDerivAt_cosh (κ * x)).comp x hA
  have hs : HasDerivAt (fun y => Complex.sinh (κ * y)) (Complex.cosh (κ * x) * κ) x :=
    (Complex.hasDerivAt_sinh (κ * x)).comp x hA
  have hN := (hc.const_mul L).add (hs.div_const B)
  have hDd := (hs.const_mul (B * L)).add hc
  have h := hN.div hDd hD
  unfold Zfilm
  convert h using 1
  simp only [Pi.add_apply]
  field_simp

/-- **The same with the model invariants** (Lemma 0: `A B = ξ₁ b p`, `A/B = (ξ₁/b) g` with
`A = κ ξ₁`): `∂_{ξ₁} Z = g/b − b p Z²`. -/
theorem hasDerivAt_Zfilm_model {κ B L x b p g : ℂ} (hB : B ≠ 0) (hb : b ≠ 0)
    (hD : B * L * Complex.sinh (κ * x) + Complex.cosh (κ * x) ≠ 0)
    (hκB : κ * B = b * p) (hκg : κ = g / b * B) :
    HasDerivAt (Zfilm κ B L) (g / b - b * p * Zfilm κ B L x ^ 2) x := by
  convert hasDerivAt_Zfilm hB hD using 1
  have e : κ / B = g / b := by rw [hκg]; field_simp
  have e2 : κ / B * B ^ 2 = b * p := by rw [← hκB]; field_simp
  linear_combination (-1 : ℂ) * e + Zfilm κ B L x ^ 2 * e2

/-! ## 2. The hierarchy for the Laurent coefficients -/

/-- Taylor coefficients of `g(s²) = (1 + τ_R s²)/(1 + τ_ℓ s²)` in powers of `s²`. -/
def gam (tR tl : ℝ) : ℕ → ℝ
  | 0 => 1
  | m + 1 => (tR - tl) * (-tl) ^ m

/-- `(1 + τ_ℓ p) Σ_{m ≤ M} γ_m p^m = 1 + τ_R p − (τ_R − τ_ℓ)(−τ_ℓ)^M p^{M+1}`: the `γ_m` are the
Taylor coefficients of `g` at `p = 0`. -/
theorem gam_series (tR tl p : ℝ) (M : ℕ) :
    (1 + tl * p) * ∑ m ∈ Finset.range (M + 1), gam tR tl m * p ^ m =
      1 + tR * p - (tR - tl) * (-tl) ^ M * p ^ (M + 1) := by
  induction M with
  | zero => simp [gam]; ring
  | succ M ih =>
    rw [Finset.sum_range_succ, mul_add, ih]
    simp only [gam]
    ring

/-- `q₀ = 1/b_s`, `q_{k+1} = c_k` (`k = 0, …, 6`): the coefficients of `s F(s)`, as functions of
`ξ₁` (the value `0` for `k ≥ 8` is never used). -/
noncomputable def q (bs b tR tl : ℝ) : ℕ → ℝ → ℝ
  | 0 => fun _ => 1 / bs
  | 1 => fun xi => c0 bs xi b tR tl
  | 2 => fun xi => c1 bs xi b tR tl
  | 3 => fun xi => c2 bs xi b tR tl
  | 4 => fun xi => c3 bs xi b tR tl
  | 5 => fun xi => c4 bs xi b tR tl
  | 6 => fun xi => c5 bs xi b tR tl
  | 7 => fun xi => c6 bs xi b tR tl
  | _ => fun _ => 0

/-- Coefficient of `s^n` in `g(s²)/b − b s² F(s)²` when `s F = Σ_k Q_k s^k`. -/
noncomputable def rhs (b tR tl : ℝ) (Q : ℕ → ℝ → ℝ) (n : ℕ) (xi : ℝ) : ℝ :=
  (if n % 2 = 0 then gam tR tl (n / 2) / b else 0) -
    b * ∑ j ∈ Finset.range (n + 1), Q j xi * Q (n - j) xi

lemma hasDerivAt_of_deriv {f : ℝ → ℝ} {x d : ℝ} (hf : DifferentiableAt ℝ f x)
    (h : deriv f x = d) : HasDerivAt f d x :=
  h ▸ hf.hasDerivAt

set_option linter.flexible false in
set_option maxHeartbeats 1000000 in
-- The order-6 identity expands products of c₀, …, c₅; the default heartbeat budget is too small.
/-- **The coefficients `c₀, …, c₆` solve the hierarchy** (`n = 0, …, 6`). -/
theorem hierarchy {bs b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) (n : ℕ) (hn : n ≤ 6) (xi : ℝ) :
    HasDerivAt (q bs b tR tl (n + 1)) (rhs b tR tl (q bs b tR tl) n xi) xi := by
  interval_cases n
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c0; fun_prop) ?_
    simp only [q, deriv_c0_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c1; fun_prop) ?_
    simp only [q, deriv_c1_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c2; fun_prop) ?_
    simp only [q, deriv_c2_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0 c1
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c3; fun_prop) ?_
    simp only [q, deriv_c3_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0 c1 c2
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c4; fun_prop) ?_
    simp only [q, deriv_c4_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0 c1 c2 c3
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c5; fun_prop) ?_
    simp only [q, deriv_c5_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0 c1 c2 c3 c4
    field_simp
    ring
  · refine hasDerivAt_of_deriv (by simp only [q]; unfold c6; fun_prop) ?_
    simp only [q, deriv_c6_xi hb hbs]
    simp [rhs, q, gam, Finset.sum_range_succ]
    unfold c0 c1 c2 c3 c4 c5
    field_simp
    ring

/-- **Initial data** of perfect contact: `q_{n+1}(0) = 0`. -/
theorem hierarchy_init (bs b tR tl : ℝ) (n : ℕ) (hn : n ≤ 6) : q bs b tR tl (n + 1) 0 = 0 := by
  interval_cases n <;> simp [q, c0, c1, c2, c3, c4, c5, c6]

/-- **Uniqueness.** Any family of functions of `ξ₁` with `Q₀ = 1/b_s` that satisfies the hierarchy
on `ℝ` for `n = 0, …, 6` with the initial data `Q_{n+1}(0) = 0` coincides with `q`, i.e.
`Q_{k+1} = c_k` for `k = 0, …, 6`. -/
theorem hierarchy_unique {bs b tR tl : ℝ} (hb : b ≠ 0) (hbs : bs ≠ 0) (Q : ℕ → ℝ → ℝ)
    (h0 : ∀ x, Q 0 x = 1 / bs)
    (hQ : ∀ n ≤ 6, ∀ x, HasDerivAt (Q (n + 1)) (rhs b tR tl Q n x) x)
    (hinit : ∀ n ≤ 6, Q (n + 1) 0 = 0) :
    ∀ k ≤ 7, Q k = q bs b tR tl k := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro hk
    rcases k with _ | n
    · funext x
      simp [q, h0]
    · have hn : n ≤ 6 := by omega
      have hlow : ∀ j ≤ n, Q j = q bs b tR tl j := fun j hj => ih j (by omega) (by omega)
      have hrhs : ∀ x, rhs b tR tl Q n x = rhs b tR tl (q bs b tR tl) n x := by
        intro x
        have hs : ∑ j ∈ Finset.range (n + 1), Q j x * Q (n - j) x =
            ∑ j ∈ Finset.range (n + 1), q bs b tR tl j x * q bs b tR tl (n - j) x := by
          apply Finset.sum_congr rfl
          intro j hj
          rw [Finset.mem_range] at hj
          rw [hlow j (by omega), hlow (n - j) (by omega)]
        simp only [rhs, hs]
      have hD : ∀ x, HasDerivAt (fun y => Q (n + 1) y - q bs b tR tl (n + 1) y) 0 x := by
        intro x
        have h := (hQ n hn x).sub (hierarchy (tR := tR) (tl := tl) hb hbs n hn x)
        rw [hrhs x, sub_self] at h
        exact h
      have hconst := is_const_of_deriv_eq_zero (fun x => (hD x).differentiableAt)
        (fun x => (hD x).deriv)
      funext x
      have h := hconst x 0
      simp only [hinit n hn, hierarchy_init bs b tR tl n hn, sub_zero] at h
      linarith

end FdtrGauge.Riccati
