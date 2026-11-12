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
