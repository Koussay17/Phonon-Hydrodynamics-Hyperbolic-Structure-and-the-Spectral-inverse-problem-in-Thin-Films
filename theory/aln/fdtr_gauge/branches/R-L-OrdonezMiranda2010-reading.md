# Reading note — Ordóñez-Miranda & Alvarado-Gil, J. Heat Transfer 132, 061302 (2010), doi:10.1115/1.4000748

Full text read (9 pages, 061302-1 to -9), 12 November 2026. The PDF is kept in the local reference library, not in the
repository.

## What the paper does
- **Model.** First-order (Taylor-expanded) dual-phase-lag law q + τ_q ∂_t q = −k(∇T + τ_T ∂_t ∇T) (their Eq. 2), i.e. the 1D
  GK law of note 7 with τ_q = τ_R, τ_T = τ_ℓ. Finite layer (thickness l) on a semi-infinite layer, modulated surface flux at
  x = 0 (Dirac source), front-face temperature. Interface: temperature jump R·(flux) and flux continuity (Eqs. 10a,b) — the
  local convention with a resistance; the substrate also obeys DPL (lags τ_q0, τ_T0).
- **Equal lags.** States explicitly (p. 3) that for τ_q = τ_T the DPL predictions coincide with Fourier at every frequency.
  This is our Δ ⊂ Σ, in exactly our setting (with credit to the classical DPL literature).
- **Method.** Only in the two limits ε₀/ε ≪ 1 or ≫ 1 (substrate effusivity much lower or much higher than the layer's) and
  R ε₀ ≪ 1, where the substrate drops out (Eqs. 14, 15). A normalised signal Γ = Θ₂/Θ₁ (Eq. 16) — the ratio of the two limiting
  configurations, i.e. the same layer on two different substrates — whose phase zeros, extrema and amplitude-ones give
  polynomial equations (19)–(23) for (τ_q, τ_T).
- **Known thickness.** Explicit: "for a layer of known thickness and thermal diffusivity (known cutoff frequency) … determine
  its phase lags" (Sec. 3.1.1); later "for a layer of known thickness, … determine its phase lags and thermal diffusivity
  simultaneously" from three characteristic frequencies (Sec. 3.3, Eqs. 35–37), then the effusivity from a non-normalised
  amplitude (Eq. 38) and k = ε√α.
- **Uniqueness.** One sentence (Sec. 3.1.1): positivity and validity of Eq. 19 for all zeros "ensure that the time phase lags
  are unique", with a sketched argument (perturb τ_q, τ_T and show the right-hand side changes). This is a local uniqueness
  claim for their characteristic-frequency equations, not a theorem on the response map. Approximate formulas (23), (35–37)
  are tested numerically (errors up to ~50–95 % on τ_q for the first zeros).
- **Not in the paper.** No rank or Jacobian analysis; no general effusivity ratio (only the two extreme limits); no matched
  effusivity case; no statement for unknown thickness or the depth-scaling gauge; no function-level injectivity; no
  degenerate-set classification beyond τ_q = τ_T; no resistance identifiability (R is assumed negligible); no finite-frequency
  non-uniqueness result; single-configuration identifiability not addressed (their normalisation uses two substrates).

## Verdict for workstream B
PARTIAL PRIOR ART, now confirmed from the full text:
- Anticipated: the equal-lag Fourier equivalence in this exact configuration; a constructive determination of (τ_q, τ_T, α, ε)
  for a layer of known thickness in two extreme-effusivity regimes with a two-substrate normalisation; a sketched local
  uniqueness claim for those formulas.
- Not anticipated (potentially novel within the literature searched): Σ = Δ exactly, for every effusivity ratio including
  b = b_s; global injectivity of θ = (ξ₁, b, τ_R, τ_ℓ) ↦ Z from a single configuration with unknown thickness (identifying the
  gauge-invariant combinations), with explicit fibres on Δ; the interface-resistance extension (rank 5, Theorem 2′);
  two-frequency local-but-not-global identifiability (certified twin).
Consistency: their need for a known thickness is exactly the gauge of Theorem G — with l known the gauge is broken, and our
theorem then gives (a, λ, ρc, τ_q, τ_T) as well.
Required in the paper: cite this work for the equal-lag reduction in the layered setting and for the constructive
known-thickness determination; state the differences above.
