# 10-proof-status


## Workstream B (11 November)
Proof audit `branches/B-proof-audit.md`. PROVED: Theorem 1 (global identifiability and fibres), Theorem 2 (rank table; Σ = Δ),
Corollary 2.1, Theorems 1′/2′ (interface resistance), Theorem G (gauge), Proposition TW (two-frequency twin, re-certified).
PROVED WITH REPAIRS: Theorem N2, Proposition NR (D1, D2; repairs written in proofs/B-analytic.md addendum).
GAP: Proposition L structural part (D7) → CONJECTURE. Formal verification of the polynomial identities: in progress.
Status conditional on Assumption M (1D, local interface temperature, no extra GK boundary condition), exact noise-free data.

## Workstream A (12 November)
Numerical audit `branches/A-numerical-audit.md`: RELIABLE WITH CAVEATS. Claim as admitted: in the note-16 Au/AlN/sapphire
twin, the static-data range of the film memory time is not resolvable by FDTR when the instrument and interface nuisances
are free (N ⊇ N2), by ≥ 4.2× (0.01°) and ≥ 34× (0.1°) against the half-width; marginal for instrument-only nuisances at
500 nm and 0.01°; margins halve for the full width and shrink with smaller interface resistances. Hostile review pending.

## After hostile review (12 November) — status changes
- Workstream B theorems are now stated as **conditional on the local interface-temperature convention (c_n = 0)**. Under the
  entropic convention implemented in the twin, the referee finds rank 4 on Δ and D₂ (resonant GK film differs from Fourier by
  0.089°). Degenerate-set statements are therefore convention-dependent; classification of c_n ≠ 0 pending (revision R-B1).
- Novelty of B downgraded: partially anticipated (DPL literature; Ordóñez-Miranda & Alvarado-Gil 2010; classical
  identifiability methods); full text of the 2010 paper still unread.
- Workstream A headline restricted pending revision R-A1: computed widths only (30–300 nm), 0.1° verdict as the robust claim,
  sensitivity box over interface resistances and regulator normalisation; δ = τ_R − τ_ℓ range open (R-A2).

## Formal verification (12 November) — `proofs/B-formal.md`, Lean sources `proofs/B-formal/`
FORMALLY VERIFIED (Lean 4.34.1, Mathlib v4.34.1, 68 statements, axioms propext/Classical.choice/Quot.sound only, no sorry):
Theorem 1 step 5 (rational identity, fibre form) and the algebraic part of step 4; the 2×3 two-load kernel lemma (pointwise);
Theorem 2 Proof A endgame (kernel {0} off Δ, line on Δ∖D₂, D₂ case); the moment minors G₀₁₂₃ ≡ 0, G₀₁₂₄, G₀₂₃₄ (general b and
b = b_s); N2 step 3 (t_K values, index-set enumeration, nonzero leading coefficient); Riccati hierarchy uniquely determining
c₀…c₆. Not formalised: the analytic steps (Lemmas 0–3, termwise differentiation of Laurent coefficients, meromorphic parts
of the kernel argument, Lemmas T/LW, Theorems 1′/2′, NR, TW, G). Status of Theorems 1, 2, N2: PROVED (audited), finite
core FORMALLY VERIFIED. Remark for the write-up: step 4 uses g(0) = 1 with p = 0 outside D; harmless for rational functions.

## Revision R-B (12 November) — `proofs/B-revision.md` (not yet audited)
DERIVED (audit pending): interface-convention classification (local κ = 0: Σ = Δ; any tied ε ≠ 0 incl. entropic and
Beardo-kinetic: no degeneracy on Δ, only Π ∩ {R = 0}; free c_n: τ_n = c_nC identifiable, resonance degeneracy returns);
twin observable is T(0), Cattaneo layers convention-free (Proposition O, from the code). B-4: known transducer with unknown
R₁, R₂ — rank 6 / 5 / 4 (local) as the referee found; R₁ never in a kernel; Theorems T1–T3, T4(i)–(iv) and TW3
(three-frequency certified twin). CONJECTURED: N = 3 on the hypersurface; entropic T4″.

## Revision R-A (12 November) — `07-experiments.md` R-A1…R-A6
- Width at 500 nm computed (inner bounds): ≥ 6.13 % on the data meshes, 5.0–5.4 % held out. Regulator normalisation c
  enters only through L/c (exact).
- Cattaneo projection (τ_ℓ = 0): not resolvable at 0.1° for every declared nuisance set over the whole sensitivity box
  (σ/w ≥ 7, ≤ 0.005 bit) — NUMERICALLY DEMONSTRATED, model-conditional. At 0.01°: ≥ 3.3 on the note-16 stack; flips only at
  500 nm for R₁ ≤ 5e-9 with c ≤ 3/4.
- Kinetic nonlocal time of the same regularised operator: τ_ℓ = |e|²(u·Ĉ⁻¹u)/K² (DERIVED, audit pending; reduces to 4τ_c/5 for
  the grey Callaway case); τ_ℓ/τ_mem = 1.8 / 3.0 / 6.3 / 9.2 at 30/100/300/500 nm (NUMERICALLY SUPPORTED, RTA converged
  < 0.5 %). Static-data range of δ = τ_R − τ_ℓ is 12–27× wider than that of τ_mem at 300–500 nm (inner bounds, 9×9×6), and
  would be resolvable at 0.01° (1D, local convention). "The Cattaneo comparison is conservative": CONTRADICTED under the
  kinetic closure (supported only under the grey closure). Audit and counterexample attack pending.
