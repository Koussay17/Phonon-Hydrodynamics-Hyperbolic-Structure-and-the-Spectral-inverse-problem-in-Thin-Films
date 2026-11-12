"""Workstream B, formal check: emit ResearchProofs/FdtrGauge/MomentCoefficients.lean.

The Laurent coefficients c_0..c_6 (perfect contact) are taken verbatim from
derivations/B-analytic/07_moment_jacobian.out (re-derived in 01_coefficients.py by series division
and, independently, by integrating the Riccati hierarchy).  The stated partial derivatives are
SymPy's; Lean re-proves every one of them from the definitions (`deriv` of the restriction to one
variable), so this script cannot introduce an unchecked statement.
Run:  python 03_generate_lean_coefficients.py <output .lean path>
"""
import sys

import sympy as sp

xi, b, bs, tR, tl = sp.symbols('xi b bs tR tl')
loc = {'xi': xi, 'b': b, 'bs': bs, 'tR': tR, 'tl': tl}
C_TXT = [
    "-xi*(b - bs)*(b + bs)/(b*bs**2)",
    "xi**2*(b - bs)*(b + bs)/bs**3",
    "-xi*(3*b**4*xi**2 - 4*b**2*bs**2*xi**2 - 3*bs**4*tR + 3*bs**4*tl + bs**4*xi**2)/(3*b*bs**4)",
    "xi**2*(3*b**4*xi**2 - 5*b**2*bs**2*xi**2 - 3*bs**4*tR + 3*bs**4*tl + 2*bs**4*xi**2)/(3*bs**5)",
    "-xi*(15*b**6*xi**4 - 30*b**4*bs**2*xi**4 - 20*b**2*bs**4*tR*xi**2 + 20*b**2*bs**4*tl*xi**2"
    " + 17*b**2*bs**4*xi**4 + 15*bs**6*tR*tl + 10*bs**6*tR*xi**2 - 15*bs**6*tl**2"
    " - 10*bs**6*tl*xi**2 - 2*bs**6*xi**4)/(15*b*bs**6)",
    "xi**2*(45*b**6*xi**4 - 105*b**4*bs**2*xi**4 - 75*b**2*bs**4*tR*xi**2 + 75*b**2*bs**4*tl*xi**2"
    " + 77*b**2*bs**4*xi**4 + 45*bs**6*tR*tl + 60*bs**6*tR*xi**2 - 45*bs**6*tl**2"
    " - 60*bs**6*tl*xi**2 - 17*bs**6*xi**4)/(45*bs**7)",
    "-xi*(315*b**8*xi**6 - 840*b**6*bs**2*xi**6 - 630*b**4*bs**4*tR*xi**4 + 630*b**4*bs**4*tl*xi**4"
    " + 756*b**4*bs**4*xi**6 + 420*b**2*bs**6*tR*tl*xi**2 + 714*b**2*bs**6*tR*xi**4"
    " - 420*b**2*bs**6*tl**2*xi**2 - 714*b**2*bs**6*tl*xi**4 - 248*b**2*bs**6*xi**6"
    " + 105*bs**8*tR**2*xi**2 - 315*bs**8*tR*tl**2 - 420*bs**8*tR*tl*xi**2"
    " - 126*bs**8*tR*xi**4 + 315*bs**8*tl**3 + 315*bs**8*tl**2*xi**2 + 126*bs**8*tl*xi**4"
    " + 17*bs**8*xi**6)/(315*b*bs**8)",
]
C = [sp.sympify(t, locals=loc) for t in C_TXT]
VARS = [('xi', xi), ('b', b), ('tR', tR), ('tl', tl)]


def lean(e):
    return sp.sstr(e).replace('**', '^')


def wrap(s, indent=4, width=96):
    """Break a long Lean expression after ' + ' / ' - ' (the operator stays at the line end)."""
    import re
    toks = re.split(r'(?<= [+-] )', s)
    out, line = [], ''
    for tok in toks:
        if len(line) + len(tok) > width and line:
            out.append(line.rstrip())
            line = ' '*indent + tok
        else:
            line += tok
    out.append(line)
    return '\n'.join(out)


HEADER = r'''import Mathlib.Analysis.Calculus.Deriv.Add
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
'''


def main(path):
    lines = [HEADER]
    lines.append('/-! ## The coefficients (07_moment_jacobian.out, verbatim) -/\n')
    for k, ck in enumerate(C):
        lines.append(f'/-- The Laurent coefficient `c_{k}` (perfect contact). -/')
        lines.append(f'noncomputable def c{k} (bs xi b tR tl : ℝ) : ℝ :=')
        lines.append(wrap('  ' + lean(ck)) + '\n')
    lines.append('/-- The family `k ↦ c_k`, `k = 0, …, 6`. -/')
    lines.append('noncomputable def coeff (bs : ℝ) : Fin 7 → ℝ → ℝ → ℝ → ℝ → ℝ :=')
    lines.append('  ![c0 bs, c1 bs, c2 bs, c3 bs, c4 bs, c5 bs, c6 bs]\n')
    lines.append('/-- `∇f(θ) = (∂_{ξ₁} f, ∂_b f, ∂_{τ_R} f, ∂_{τ_ℓ} f)(θ)`. -/')
    lines.append('noncomputable def grad (f : ℝ → ℝ → ℝ → ℝ → ℝ) (xi b tR tl : ℝ) : Fin 4 → ℝ :=')
    lines.append('  ![deriv (fun y => f y b tR tl) xi, deriv (fun y => f xi y tR tl) b,')
    lines.append('    deriv (fun y => f xi b y tl) tR, deriv (fun y => f xi b tR y) tl]\n')
    lines.append('/-! ## The forms displayed in Proof C (`proofs/B-analytic.md`, §5) -/\n')
    proof_forms = [
        ("c0", "xi * (bs ^ 2 - b ^ 2) / (b * bs ^ 2)"),
        ("c1", "xi ^ 2 * (b ^ 2 - bs ^ 2) / bs ^ 3"),
        ("c2", "xi * (tR - tl) / b - xi ^ 3 * (3 * b ^ 4 - 4 * b ^ 2 * bs ^ 2 + bs ^ 4) / (3 * b * bs ^ 4)"),
        ("c3", "xi ^ 2 * (3 * bs ^ 4 * (tl - tR) + (3 * b ^ 4 - 5 * b ^ 2 * bs ^ 2 + 2 * bs ^ 4) * xi ^ 2) / (3 * bs ^ 5)"),
    ]
    for name, form in proof_forms:
        lines.append(f'theorem {name}_eq_proofC {{bs xi b tR tl : ℝ}} (hb : b ≠ 0) (hbs : bs ≠ 0) :')
        lines.append(f'    {name} bs xi b tR tl =')
        lines.append(f'      {form} := by')
        lines.append(f'  unfold {name}')
        lines.append('  field_simp')
        lines.append('  ring\n')
    lines.append('/-! ## Partial derivatives -/\n')
    for k, ck in enumerate(C):
        for name, v in VARS:
            d = sp.factor(sp.diff(ck, v))
            pt = {'xi': 'xi', 'b': 'b', 'tR': 'tR', 'tl': 'tl'}[name]
            args = {'xi': 'y b tR tl', 'b': 'xi y tR tl', 'tR': 'xi b y tl', 'tl': 'xi b tR y'}[name]
            lines.append(f'theorem deriv_c{k}_{name} {{bs xi b tR tl : ℝ}} (hb : b ≠ 0) (hbs : bs ≠ 0) :')
            lines.append(f'    deriv (fun y => c{k} bs {args}) {pt} =')
            lines.append(wrap('      ' + lean(d), indent=8) + ' := by')
            lines.append(f'  simp only [c{k}]')
            lines.append('  deriv_tac\n')
    for k, ck in enumerate(C):
        ds = [lean(sp.factor(sp.diff(ck, v))) for _, v in VARS]
        lines.append(f'/-- `∇c_{k}`. -/')
        lines.append(f'theorem grad_c{k} {{bs xi b tR tl : ℝ}} (hb : b ≠ 0) (hbs : bs ≠ 0) :')
        lines.append(f'    grad (c{k} bs) xi b tR tl =')
        lines.append('      ![' + ',\n        '.join(ds) + '] := by')
        lines.append(f'  rw [grad, deriv_c{k}_xi hb hbs, deriv_c{k}_b hb hbs, deriv_c{k}_tR hb hbs,')
        lines.append(f'    deriv_c{k}_tl hb hbs]\n')
    lines.append('end FdtrGauge.Moments\n')
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(lines))


if __name__ == '__main__':
    main(sys.argv[1])
