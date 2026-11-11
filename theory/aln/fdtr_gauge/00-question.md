# 00 — Question

Two workstreams of paper 1 (paper/PRIOR_ART.md, "Results still needed before submission", items 3 and 4).

## Workstream A — FDTR consequence of the film memory-time gap

**Question.** The static data (mode lifetimes within ±10 %, κ within ±1 %, three data meshes) leave the
film-regularised memory time τ_mem(L) of the conserving closure of AlN at 300 K undetermined within the admissible
ranges of fc3-memory E33/E37 (ratio to the reference, three data meshes; held-out 12×12×8 in brackets):

| class | 30 nm | 100 nm | 300 nm |
|---|---|---|---|
| orb 0.1 | [0.973, 1.028] ([0.973, 1.028]) | [0.963, 1.038] ([0.966, 1.033]) | [0.949, 1.054] ([0.954, 1.041]) |
| cal1nz / cal3wnz / cal3w | — | up to [0.921, 1.080] ([0.934, 1.059]) | — |

Is that gap visible in a frequency-domain thermoreflectance (FDTR) measurement of the Au/AlN/sapphire scenario
(500 nm AlN, 10 kHz – 200 MHz, phase noise 0.01–0.1°, amplitude noise 1 %), and under which nuisance-parameter sets?

**Formulation.** τ_mem = −K′(0)/K(0) is the Cattaneo relaxation time of the conserving closure (paper/PRIOR_ART.md,
central contribution 1). The film quantity is mapped to τ_R of the AlN layer in `src/twin_fdtr.py`. Observable: the
probe-averaged complex response H(f; θ) of `gaussian_response`, with instrument gain, phase offset and delay
(`instrument`). Quantities to produce:
1. Δφ(f) = φ(f; τ_R(1+w)) − φ(f; τ_R) for the admissible relative changes w, with all other parameters fixed.
2. The detectability threshold w*(noise, nuisance set, band): the smallest |w| for which the τ_R direction is
   resolved at 1σ (and 2σ) in the local Fisher analysis, with nuisance parameters freed (sets declared in advance:
   N0 none; N1 instrument gain/phase/delay; N2 + interface resistances R_Au/AlN, R_AlN/sapphire; N3 + λ_z and C of AlN;
   N4 + τ_ℓ (GK) and Au properties; N5 + beam radii).
3. Comparison of w* with the admissible half-widths w(L); a statement of the form "the static-data gap of τ_mem(L) is
   (is not) resolvable by FDTR at noise σ with nuisance set N".

**What counts as a solution.** Converged numbers (Hankel quadrature, frequency grid, finite-difference step of the
Jacobian, all checked) for w* at each (σ, N) and band, plus Δφ(f) curves, with the 500 nm base value τ_R(500 nm)
derived from the campaign's operator, not assumed.

**Falsification.** The main hypothesis H_A ("the gap is far below FDTR resolution once instrument and interface
nuisances are free") is false if w* ≤ the admissible half-width for any realistic N ⊇ N2 at σ_φ = 0.1°.

## Workstream B — rank and gauge of the 1D non-Fourier front-face response

**Setting** (note 7). One film of thickness e on a semi-infinite substrate; 1D transfer, constant coefficients, zero
initial state, no interior source, unit flux absorbed at the front face, response = front-face temperature. Film law
Cattaneo or Guyer–Krumhansl with effective conductivity λ_eff(p) = (λ + 3ℓ²pρc)/(1 + τ_R p) (historical α = 2; the
conserving α = 1/3 only changes the numerical factor in τ_ℓ = (1+α)ℓ²/a). The film enters only through
σe = ξ₁ √(p(1+τ_R p)/(1+τ_ℓ p)) and λ_eff σ = b √(p(1+τ_ℓ p)/(1+τ_R p)), ξ₁ = e/√a, b = √(λρc).

**Claim to prove or refute (B1).** With the substrate effusivity b_s > 0 known (Fourier substrate), the map
F: θ = (ξ₁, b, τ_R, τ_ℓ) ∈ (0,∞)² × [0,∞)² ↦ Z(·; θ), Z(p) the front-face impedance on p = iω, ω > 0, has
Jacobian of rank four at every θ outside an explicit proper analytic subset Σ, and the rank on Σ is computed.
Candidates for Σ: the Fourier resonance τ_R = τ_ℓ; the matched-effusivity case b = b_s; boundaries τ = 0.
**Claim B2 (global).** F is injective on the complement of Σ (or the exact fibres are described).
**Claim B3 (gauge).** The response depends on (e, a, λ, ρc, ℓ, τ_R) only through θ; the residual one-parameter gauge
e → μe, a → μ²a, λ → μλ, ρc → μ⁻¹ρc, ℓ → μℓ is the only non-identifiability off Σ; it is broken by any of: known e,
known λ or ρc, a lateral wavenumber (finite beams with known radii and in-plane properties), an independent velocity
with a Debye relation, an interior source. Each breaking mechanism to be proved or demoted.
**Extensions.** Interface resistance R film/substrate (one more parameter; rank five?); a known Fourier transducer on
top; finitely many frequencies (how many suffice generically).

**What counts as a solution.** A written proof (THEOREM with hypotheses), the exact degenerate set, an independent
CAS verification of the rank and of Σ, and an audit. **Falsification.** A θ off the proposed Σ with rank < 4, or two
distinct θ off Σ with the same Z.
