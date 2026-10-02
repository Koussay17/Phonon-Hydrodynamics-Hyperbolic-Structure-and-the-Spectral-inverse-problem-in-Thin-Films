# Final hostile numerical review (4th red review) — 2 October 2026

Scope: note `C:\Users\Koussay\these\notes\23_Integration_resonante.tex`, archive
`theory/aln/resonance_measure/*` (A, B, C, C2, D) and runner `scripts/reproduce_resonance_measure.py`,
with focus on the 2 Oct fixes (log-space `reaction_data`; `-O` refusal). Repository read-only
(hashes identical before/after all runs). The proof, counterexample and peer red reports and the
`red_counter_*`/`red_proof_*` artefacts were NOT read. The 28 Sep numerical audit script
(`R_numerics_audit.py`) was skimmed only after all my own tests were complete, to check coverage.

Environment: research Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.3.0, python-flint 0.9.0.
Implementation under test: current `C2_weak_form.py` sha256 `2c4fa15a…`; frozen pre-fix copy
`933e7ae5…` (manifest hashes verified). References: (R1) python-flint Arb balls, definitional
`A-B`, `(A-B)/(log A-log B)`, adaptive precision, relative radius < 1e-30; (R2) mpmath 60+ digits,
identity route `B·expm1(-d)` with exact `d`. R1 and R2 agree to ≤ 1.4e-49 on all 30,021 inputs.
Errors are measured at the exact floating-point inputs; u = 2^-53 = 1.11e-16;
κ_y = Σ_j |ξ_j ∂y/∂ξ_j| / |y| (analytic, mpmath).

## Verdict

The 2 Oct log-space fix is numerically sound on the requested domain ξ ∈ [1e-12, 1e3]³: zero
non-finite outputs in 30,021 stressed events (frozen code: 231 NaN/inf plus 175 finite but 100 %-wrong
values, silent underflow to 0, with κ < 1e4), Λ relative error ≤ 5.73e-14, F backward stable
(error ≤ 2.2·u·(κ + |log base| + 3)), exact resonance gives F = 0 exactly, and the harness
reproduces 13/13 checks with `-O`/`-OO`/`PYTHONOPTIMIZE` refused. I found **no defect affecting a
reported number**. I found one **confirmed test-design defect** (the revised C2 thermal-identity
assertion is non-discriminating) and several **fragilities**: silent overflow/NaN in the selected branch
outside [1e-12,1e3]; complex-step accuracy limits; weak flux checks in the new regression;
underflow of Λ and of the capacity weights at the large-ξ end; a provenance gap.

## Artefacts (all in `experiments/`, nothing overwritten)

| File | Content |
|---|---|
| `red_numerics_final_common.py` | read-only module loader, Arb and mpmath references |
| `red_numerics_final_stress.py/.json/_raw.npz/.log` | Task 1 + Task 4 (30,021 events, seed 20261002) |
| `red_numerics_final_complex_step_v2.py/.json/.log` | Task 2 (Arb derivative reference). v1 files superseded, see "Own errors" |
| `red_numerics_final_exprel.py/.json` | Task 3 |
| `red_numerics_final_harness.py/.json/.log`, `red_numerics_final_harness_run/`, `red_numerics_final_optimize_*` | Task 5 |
| `red_numerics_final_checklist.py/.json`, `…_checklist_v2.py/.json`, `…_rootcond.py/.json` | Task 6 |

## 1. Randomized stress of F and Λ (Task 1) — OK

Sets: S1 12,000 independent log-uniform legs; S2 6,000 near resonance (|d| relative 1e-17…1e-1 and
absolute 1e-16…1e-1, both signs); S3 3,000 exact resonance (d = 0 verified with Fractions);
S4 3,000 with |d| ∈ [1e-15, 1e-13]; S5 3,000 with all legs in [1e-12, 1e-6]; S6 3,000 with
subnormal outputs (both branches and near resonance); S7 21 regression/corner cases including
(1,400,400) and the eight C2 cases.

| Quantity | Result |
|---|---|
| Non-finite F or Λ (new) | 0 / 30,021; no RuntimeWarning emitted |
| Non-finite (frozen) | 231 (225 in S1, 6 in S7), e.g. NaN at (1,400,400) |
| Λ, normal range (27,108 samples) | max rel. error 5.73e-14, median 1.4e-15 (κ_L ≤ 1e3 everywhere) |
| F, normal range, by κ_F | κ ≤ 10: 1.62e-14; ≤ 1e3: 5.70e-14; ≤ 1e4: 3.36e-13; ≤ 1e8: 3.3e-9 |
| Backward-stability ratio | max rel/(u(κ+\|log base\|+3)) = 2.22 (F), 2.15 (Λ); rel/(uκ) ≤ 73 |
| Exact resonance (3,000 + 67 others) | F returned exactly 0 in every case |
| Subnormal outputs (5,681) | \|err\| ≤ 0.90 × (u(κ+\|log base\|+3)\|ref\| + 1 subnormal ulp) |
| Huge logs (\|log base\| > 300) | Λ rel. error ≤ 5.7e-14 = 1.03·u·\|log base\| at worst |
| Affinity `b@α` vs exact d | ≤ 2.3e-13 absolute (2 ulp of 1e3) |
| Branch formulas compared near d≈0 | relative gap ≤ 5.7e-14 (continuous at the switch) |
| Same, large \|d\| (S1) | gap up to 2.5e-6: the NON-selected formula loses accuracy (subnormal base); selecting max(A,B) is necessary |
| Structure | F·d_code ≤ 0 in 30,021/30,021: per-event entropy production −F·d ≥ 0 always |

Near resonance the relative error of F is unbounded (up to 100 %, 324 cases with F = 0 while
|F_exact| > 1e-300, 35 cases off by a factor 2) because d = ξ_p − ξ_a − ξ_b is formed in floating point;
these have κ_F ≈ 3.7e16 and rel/(uκ) ≈ 0.24, i.e. inherent conditioning, not instability (already stated in C2
§5). The log-space evaluation costs accuracy relative to the frozen product where the latter was finite:
for κ < 10 the new code reaches 1.6e-14 (ξ→0+, \|log B\| ≈ 80) versus ≤ 8e-16. Harmless against every
1e-12 tolerance (≥ 17× margin). Classification: **OK — FLOATING-POINT LIMITED at the conditioning level.**

## 2. Complex-step derivatives, h = 1e-28 (Task 2)

Reference: rigorous Arb, definitional ∂F = A∂logA − B∂logB, ∂Λ = (F∂d − d∂F)/d², limit formula at d = 0.
5,305 inputs (random, a d-sweep through 0 and across the series cutoff for six base scales, exact
resonance, an underflow sweep, ξ→0+).

* **F through `reaction_data`: OK.** Normwise relative error ≤ 2.5e-14 whenever |F| ≳ 1e-280;
  ≤ 6.9e-15 across d = 0 and the branch switch; ≤ 2.1e-14 at exact resonance.
* **Λ through `reaction_data`: FRAGILITY.** Error ≈ (2–4)·u/|d| on the expm1 side of the cutoff:
  4.6e-11 for 1e-5 ≤ |d| < 3e-5, 6.2e-12 (3e-5…1e-4), 2.8e-12 (1e-4…1e-3), 1.8e-13 (1e-3…1e-2),
  ≤ 2.5e-14 elsewhere (including the series side and d = 0). Cause: NumPy complex division
  Im(expm1(z)/z) = ε(e^x − expm1(x)/x)/x cancels. A complex-step Jacobian of K/Λ checked with the harness's
  own 1e-11 tolerance would fail spuriously near |d| ≈ 1e-5. The harness never complex-steps
  `reaction_data` (all its Jacobian checks use `direct_products=True`), so the fix's complex-step
  validity was not tested by the harness.
* **Underflow of the step: FRAGILITY.** For |F| ∈ [1e-290, 1e-280] the error is 1.1e-7 (F) and 2.6e-5 (Λ); below
  1e-290 the derivative is lost entirely (error 1.0). This occurs inside [1e-12, 1e3]: for example ξ_p ≳ 645 with d < 0.
* **Harness's own path (direct Bose products): FRAGILITY outside the tested states.** Error up to
  7.2e-5 at exact resonance with ξ ≈ 1e-12, 3.6e-5 for ξ→0+, 4.7e-10 at ξ ≈ 1e-7 (imaginary-part cancellation
  ∝ 1/n). At the harness states (ξ ≥ 0.69) it is 1.6e-16 and 3.0e-16, so the reported checks are valid.

## 3. `exprel_stable` (Task 3) — OK (derivative fragility above)

Values compared with mpmath at 60 digits: series side |z| < 1e-5 ≤ 1.00 u; real expm1 side ≤ 2.48 u;
window 0.95e-5…1.05e-5 ≤ 2.34 u; complex arguments at 7 angles ≤ 1.00 u (series) and ≤ 4.05 u (expm1).
The truncation term z⁶/7! is 2e-34 at the cutoff. The complex-step derivative is ≤ 2.9e-16 on the series side
but 8.3e-11 at x = 1.0038e-5, 5.1e-11 for x in [1.1e-5, 1e-4) and 4.4e-12 for x in [1e-4, 1e-3);
observed/(2u/|x|) ≤ 3.74. `expm1(z)/z` overflows prematurely for z ∈ (709.78, 716.36) (true value
representable); only the non-selected branch ever receives z > 0, so this is harmless.

## 4. Can `np.errstate` hide a genuine overflow/NaN in the SELECTED branch? (Task 4) — YES; FRAGILITY

Inside [1e-12, 1e3]³ it does not: 0 / 30,021 selected-branch evaluations trapped under
`errstate(over/invalid/divide='raise')`, and no non-finite output occurred. The bound is max(A,B) ≤ ~1e36.
Constructed failures inside the note's declared open domain ξ > 0, or just outside it. All are silent
with the errstate in place. "Without errstate" means the same arithmetic run with warnings enabled.

| ξ | Code output | True value | Without errstate |
|---|---|---|---|
| (2e-160, 1e-160, 1e-160), d = 0 | F = NaN, Λ = inf | F = 0, Λ = 5e479 | overflow + invalid warnings |
| (2e-160 + 4 ulp, 1e-160, 1e-160) | F = −inf | F = −6.32e304 (representable) | overflow warning |
| (3e-110, 1e-110, 1e-110) | F = −inf | F = −3.33e219 (representable) | overflow warning |
| (−0.5, 1, 1), domain violation | NaN, NaN | — (frozen: finite garbage F = −5.84) | invalid warning |
| (1, 0, 1) | NaN, inf | — | divide-by-zero still warned (not suppressed) |

The threshold for ξ = (2x, x, x) is x ≈ 1.41e-103. This is physically irrelevant, but the note explicitly does not
prove domain invariance, so a trajectory leaving ξ > 0 now yields silent NaN. A minimal fix is to evaluate each
branch only on its mask (boolean indexing, so no errstate is needed), or to assert finiteness and min ξ > 0 on
the selected outputs.

## 5. Fresh harness run and optimize refusal (Task 5) — OK

`python -B scripts/reproduce_resonance_measure.py --output-dir experiments/red_numerics_final_harness_run`:
exit 0, 6.9 s, **13/13 checks true** (A_mass, A_energy_kernel, A_entropy_identity,
B_exact_scalar_calculus, B_interpolation_negative_control, B_energy_kernel, C_assertions,
D_independent_poisson_sums, D_positive_heating, D_regular_heating_order, D_small_drift_wrong_measure,
C2_finite_weak_form_assertions, C2_separated_variables_finite_and_accurate). Every log has empty stderr.
The copied scripts are hash-identical to the repository, and the JSONs record the matching `script_sha256`.
Refusals: the runner under `-O`, `-OO`, `PYTHONOPTIMIZE=1` and `PYTHONOPTIMIZE=2` exits 1 with
"do not run under python -O". A, B, C and C2 (copies) exit 1 under `-O` and under `PYTHONOPTIMIZE=1`.
An AST census shows each guard precedes the first assert. Minor issues:
(i) D-limits-check has no guard; it has 0 asserts and its checks are evaluated in the runner, so this is
harmless, but "every script refuses" is literally false. (ii) The runner refuses only after creating the output
directory and copying 5 scripts, which leaves a half-populated directory that `exist_ok=False` then blocks.

## 6. Original red-numerics checklist (Task 6)

**(a) C2 tolerance change.** The frozen module reproduces the recorded failure bit for bit:
ρ = −1e-4, E = 8.1378800375857028e-11, PS = 8.1378800376373305e-11, |E−PS| = 5.163e-22 > 8.14e-23.
The current module gives identical values. The original predicate fails for every |ρ| ≤ 1e-4.
The diagnosis is **confirmed**: |E − PS − T| ≤ 1.57·u|E| for all ten ρ, where T = Σ ωΛΔ(a − βΔ) is summed exactly.
Also |a − βΔ| ≤ 7.1e-16, against an a-priori rounding bound of 1.9e-14 (27× margin).
**CONFIRMED DEFECT (test design; no reported number affected):** the revised predicate
|E−PS| ≤ Σ ωΛ|Δ||a−βΔ| + 32ε(|E|+|PS|) is a triangle inequality on its own residual. It passes on
deliberately wrong states whose thermal identity is violated by 30 % (α = 1.3βe), 1,580 % and 1,136 % (random α);
the original predicate rejects all three. Its only remaining sensitivity is to F ≠ −Λa: it detects a 1e-12
relative corruption of F but not 1e-14. The Cauchy–Schwarz assertion passes on all wrong states (vacuous by
construction). A non-vacuous replacement that passes on the actual data is to assert |a−βΔ| ≤ γ_n Σ|b||α| and
|E − PS − T| ≤ c·u|E|.

**(b) Physical vs surrogate energy — OK.** The two-node physical drifts agree with 50-digit mpmath to
1.1e-15 and 1.5e-15 (0.0098043592758994659, 0.0087339718071175044). The positive-square identity gap is 1e-52,
and the surrogate drift is exactly 0.0.

**(c) Root coverage and regular-root conditioning — OK (minor).** A closed-form edge census predicts C's root
counts on all 7 meshes (0,0,0,0,2,2,2 at d = 0.9999). The mass error at N ≥ 256 is 6.5e-15–1.5e-13,
within |cot q|·xtol = 3.5e-13. Approaching criticality with C's own `scan_roots`, the relative mass error is
≈ (0.1–0.4)·u/(1−d): 1.7e-13 (1−d = 1e-4), 7.4e-10 (1e-8), 1.3e-5 (1e-12), 7 % (1−d = u). With offset 0.37,
both roots are missed for 1−d ≤ 1e-6 (N = 1024) and ≤ 1e-8 (N = 16384), consistent with the note.
C's float "exact" near-critical mass carries a 1.26e-13 error from the cancellation in 1−d². This is cosmetic
(sqrt((1−d)(1+d)) gives 1e-17) and no assertion depends on it.

**(d) Quadrature — CONVERGED.** SciPy `quad` against mpmath tanh-sinh at 30 digits, 15 integrals
(regular σ = 0.2 and 0.00625, critical σ = 0.1 and 0.0015625, empty σ = 0.2, 0.025 and 0.0125):
relative error ≤ 3.3e-15 for the 14 non-tiny integrals, and every QUADPACK estimate covers the true error. The empty-case leakage at
σ = 0.0125 (7.6e-15) lies below `epsabs` = 2e-13 yet is accurate to 3.4e-12; its assertion passes with
mass/bound = 0.019, which is a weak check. D: the i0e continuum agrees with `besseli` to ≤ 3.6e-16; the heating
mesh sums (N = 8192, σ/h ≥ 16) agree with the continuum to ≤ 3.2e-16; the observed order is 2.0017 and
drift/(βσ²A0) = 1.0004.

**(e) Note claim not covered by the harness — verified.** The fold constant
2^{1/4}Γ(5/4)/√(2π) = 0.43001999366226 matches direct quadrature. E_σ/(σ^{3/2}·const) → 0.9999928 at σ = 1e-5,
with observed order 1.49997.

**(f) Weak or vacuous checks — FRAGILITY.** Every assertion has margin ≥ 20× (the smallest is the C2
separated-variable mobility check, 4.8e-14 against 1e-12). The following checks are non-discriminating and should
not count as independent evidence:
* A_mass: the weights are built from the same analytic coarea formula.
* A_energy_kernel: b·e = O(u) holds by construction at analytic roots.
* A_entropy_identity: z^T C C^{-1} g = z^T g holds for any invertible C.
* D_positive_heating: every summand is ≥ 0 pointwise.
* C2 rhs_vs_K_alpha: F = −Λd holds by construction.
* The C2 revised thermal identity and the C2 Cauchy–Schwarz assertion (see (a)).

The new separated-variable regression has two weak spots. The `or flux_ref == 0.0` escape is taken in 2 of
8 cases ((5,2,3) and (3,1,2), both exact resonance; the flux is unchecked there, although it is exactly 0).
The `max(|ref|, 1e-300)` floor understates the (705,700,4) flux error by about 1e6: it reports 5.6e-20 against a true
4.8e-14, so the effective tolerance is ~9e-7 relative. None of its cases lies at 0 < |d| ≪ 1.
The Decimal reference uses `repr` strings; the inputs differ from the floats only at 1e-3, which is negligible.

**(g) Large-ξ underflow — FRAGILITY.** Λ underflows to 0 for ξ = (2x, x, x) with x > 372.6,
inside [1e-12, 1e3]. This contradicts the stated Λ_r > 0, although no NaN results. The capacity weight n(1+n) is
subnormal for ξ > 708.4 and zero for ξ > 709.79. `bose()` = 1/expm1 returns 0 for ξ > 709.78, although n is
representable down to ξ ≈ 745. The volume side was not converted to log space, so numerically M can lose positive
definiteness when every volume ξ ≳ 708.

**(h) Provenance — FRAGILITY.** The repository was unchanged by all runs, and the manifest hashes are verified.
The campaign-folder copies (`experiments/C2_weak_form.py` 933e7a…, A c114d0…, B 46bd24…, C b3edc1…) are the frozen
pre-fix scripts. README, RESUME, 05 and 07 record only pre-fix replays (11 and 12 checks). The only post-fix full
replay I found is mine (13 checks, `experiments/red_numerics_final_harness_run`). The note is dated
28 Sep but describes the 2 Oct fix.

## Classification summary

| # | Finding | Class |
|---|---|---|
| 1 | Log-space F/Λ on [1e-12,1e3]³: 0 non-finite, Λ ≤ 5.7e-14, F backward stable, exact zeros at d = 0 | OK |
| 2 | Near-resonance relative accuracy of F limited by forming d (κ up to 3.8e16) | OK (inherent) |
| 3 | Log-space accuracy cost ≤ ~\|log base\|·u (1.6e-14 for κ < 10, versus 8e-16 frozen) | OK |
| 4 | errstate silences genuine selected-branch overflow (ξ ≲ 1.4e-103) and domain-violation NaN | FRAGILITY |
| 5 | Complex-step ∂Λ loses up to 5 digits for 1e-5 ≤ \|d\| ≲ 1e-2; not exercised by the harness | FRAGILITY |
| 6 | Complex step (h = 1e-28) underflows for \|F\| ≲ 1e-280 (derivative lost below 1e-290) | FRAGILITY |
| 7 | Harness's direct-product Jacobian check is inaccurate for ξ ≪ 1 (up to 7e-5) | FRAGILITY |
| 8 | Revised C2 thermal-identity assertion is non-discriminating (passes at 30–1580 % violation) | CONFIRMED DEFECT (test design) |
| 9 | Diagnosis of the C2 failure (a − βΔ rounding) is quantitatively correct | OK |
| 10 | Weak flux checks in the separated-variable regression (escape and 1e-300 floor) | FRAGILITY |
| 11 | Λ and capacity weights underflow at large ξ; `bose` premature underflow | FRAGILITY |
| 12 | exprel values ≤ 4.05 u (real and complex) | OK |
| 13 | Harness 13/13; -O/-OO/PYTHONOPTIMIZE refused; D unguarded but assert-free | OK |
| 14 | Quadrature, mesh, Poisson, fold constant, physical/surrogate | OK / CONVERGED |
| 15 | No post-fix replay recorded in the campaign; campaign copies are pre-fix | FRAGILITY (provenance) |

## Own errors corrected during this review

1. `red_numerics_final_complex_step.py` (v1) computed norms in float, so ~1e-218 values underflowed when squared;
   its identity-route reference also cancels for tiny components. These produced a spurious "route gap 1.0" and the
   D4 zeros. Superseded by v2 with Arb and mpmath norms. The normwise F and Λ numbers agree between v1 and v2.
2. `red_numerics_final_checklist.py` (v1) placed breakpoints outside [−π, π] for σ = 0.2, giving a spurious 1e-4
   mesh error. Corrected in `red_numerics_final_checklist_v2.py`, where the error is 3.9e-17.

VERDICT: FIX VERIFIED — log-space F/Λ backward stable and finite on [1e-12,1e3]³; harness 13/13 with
optimize refused; no reported number found wrong; 1 confirmed test-design defect (vacuous revised C2 identity
assertion) and 8 fragilities (silent selected-branch overflow/NaN outside the stress domain, complex-step limits,
weak regression checks, large-ξ underflow, provenance). Numerically reliable within the stated bounded synthetic scope.
