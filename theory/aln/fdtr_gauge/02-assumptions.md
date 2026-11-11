# 02 — Assumptions

## Workstream A
- A1. The film memory time of E30/E33 uses a **heuristic Matthiessen regulator** d_α = |v_z,α|/(4πL) in the THz-linewidth units of the operator (physical rate |v_z|/L) added to the
  conserving event operator. It is not a film solution of the BTE. Every A result is conditional on this regulator.
- A2. τ_mem(L) (z, cross-plane) is identified with the Cattaneo τ_R of the AlN layer of the twin. τ_ℓ (nonlocal time)
  is not fixed by τ_mem and is treated as a separate parameter (nuisance, or zero for the Cattaneo variant).
- A3. The admissible ranges are **inner bounds** (local optimiser; multistart E39 running). A detectability verdict
  "below resolution" is only as strong as the width it is compared with; if E39 widens the ranges, A is re-run.
- A4. The ranges exist at 30/100/300 nm, not at 500 nm. For the 500 nm scenario: (i) the base τ_R(500 nm) is computed
  from the same operator (REF, film_event route on 12×12×8, and RTA on the converged mesh); (ii) the width at 500 nm is
  bracketed by evaluating the existing witnesses at 500 nm (inner) and by the thickness trend; the decisive output is
  w*(L) compared with w(L) at the lengths where w is known. The FDTR scenario is also evaluated at 100 and 300 nm films.
- A5. Stack values are illustrative (note 16): Au 150 W/mK, 2.49e6 J/m³K, 80 nm; sapphire 35 W/mK, 3.03e6; R 1e-8
  and 2e-8 m²K/W; AlN C = 2.41e6 J/m³K; AlN λ_z from the film κ_z(L) of the campaign operator. Beam radii as in note 16.
- A6. Noise Gaussian, independent across frequencies (1 % on ln|H|, 0.01–0.1° phase), 80 log-spaced frequencies
  10 kHz – 200 MHz unless stated. Local (Fisher) analysis; Bayesian/profile checks where the Fisher result is marginal.
- A7. Synthetic only. No statement about a real sample.

## Workstream B
- B1. 1D, linear, constant coefficients, zero initial state, no interior source, front-face flux and detection,
  semi-infinite Fourier substrate with known b_s; film Cattaneo/GK with the quadrupole of note 7 (reproduces Camacho
  de la Rosa et al. to machine precision); "local" temperature convention at the interface (note 7's model).
- B2. Response known for all ω > 0 (function-level identifiability) unless a finite frequency set is stated.
- B3. Parameter domain: ξ₁, b > 0; τ_R, τ_ℓ ≥ 0 (boundary strata Fourier/Cattaneo treated separately).
