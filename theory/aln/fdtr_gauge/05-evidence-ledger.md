# 05 — Evidence ledger

| ID | Claim | Branch | Method | Evidence | Assumptions | Status |
|---|---|---|---|---|---|---|
| E-B1 | Jacobian of θ ↦ Z has rank 4 iff τ_R ≠ τ_ℓ; rank 3 on τ_R = τ_ℓ, rank 2 on τ_R = τ_ℓ and b = b_s; b = b_s, τ_R = 0, τ_ℓ = 0 not degenerate | AB-A | Low-frequency coefficient determinants: det ∂(L₀,L₁,L₂,L₄) = −2ξ₁⁴(b²−b_s²)(τ_R−τ_ℓ)/(b⁴b_s³); at b = b_s det ∂(L₀,L₂,L₃,L₄) = 2ξ₁⁴(τ_R−τ_ℓ)²/b_s⁵; SVD (σ_min ∝ \|δ\|) | B1–B3 | DERIVED, audit pending; awaits B-M, B-S, B-X |
| E-B2 | Injective off Σ: θ recovered explicitly from the low-frequency coefficients (incl. b = b_s); fibres on Σ: line (ξ₁, b, τ, τ), quarter-plane (ξ₁′, b_s, τ, τ) if b = b_s | AB-A | Explicit inversion, 50-digit check | as above | DERIVED, audit pending |
| E-B3 | At first order in ωτ the times enter only through δ = τ_R − τ_ℓ; FDTR band ωτ ≤ 0.035, so in practice (ξ₁, b, δ) only; τ_R and τ_ℓ columns whitened cosine −0.9994 | AB-A | Expansion + Fisher | 1D | NUMERICALLY SUPPORTED |
| E-A1 | Phase change from τ_R → τ_R(1+w): δφ = w ωτ_R W(ω) + O(w², (ωτ)²); peak 0.022°/w at 63 MHz for 500 nm (range 0.014–0.034°/w), 0.012°/w at 300 nm, 0.0032°/w at 100 nm; Au/AlN contact resistance suppresses it ~30× (bare AlN 1.13°/w) | AB-A | Exact 1D quadrupole; matches twin kernel at k = 0 to 1e-15; Gaussian twin within 1–6 % at 40–200 MHz | A1–A6, 1D | NUMERICALLY SUPPORTED (1D) |
| E-A2 | 1D detectability w* (1σ, 0.01°): N0 0.14, N1 0.18, N2 0.29, N3 1.8, N4 1.5e3; at 0.1°: ×~9. Admissible half-widths 0.028–0.08: not resolvable in any declared set; closest N0 most favourable case 0.090 vs 0.08 | AB-A | Fisher, whitened | 1D, inner-bound widths | NUMERICALLY SUPPORTED (1D); awaits A-C full twin + numerical audit |

## Brief error found
The impedance formula in the B briefs had numerator and denominator swapped (thin-film limit b_s√p/B² instead of
1/(b_s√p)). Found by branch AB-A; correct form Z = (1/B)(B + b_s√p tanh A)/(b_s√p + B tanh A) (src/forward_model.py,
src/twin_fdtr.py). Corrections sent to B-M, B-S, B-X on 11 November; their results are accepted only on the correct form.
