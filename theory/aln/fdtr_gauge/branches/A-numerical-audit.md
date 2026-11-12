# Branch A-N — independent numerical audit of branch A-C (workstream A)

12 November 2026. Object: 07-experiments.md entries A1–A7, scripts `experiments/A/*.py`, raw outputs
`experiments/A/out/*.json`, forward model `C:\Users\Koussay\these\src\twin_fdtr.py` (imported read-only, not modified).
Audit scripts and raw outputs: `experiments/A_audit/` and `experiments/A_audit/out/` (list and commands in section 10).
Every number is SYNTHETIC and conditional on the twin physics and noise model of 02-assumptions.md (A1–A7).

**Claim audited (E-A3).** In the Au(80 nm)/AlN/sapphire twin, w* = σ(ln τ_R) exceeds the admissible relative width of
the film memory time (2.8–8 %) at 30/100/300/500 nm for all nuisance sets N0–N5 and phase noise 0.01–0.1°
(e.g. 500 nm, N2, 0.01°: 0.340; 300 nm, N0, 0.01°: 0.246); largest admissible phase change 1.6e-3°.

## 0. Verdict

**RELIABLE WITH CAVEATS.** Numerically the A4 table is **CONVERGED**: it is not limited by quadrature, finite-difference
step, roundoff or conditioning. An implementation that shares no code with the twin or with `indep.py` (other language,
other layer formulation, other quadrature, other differentiation, other Fisher algebra) reproduces all 270 finite-beam
entries to ≤ 7.2e-6 relative (N0–N3: ≤ 4.3e-7), and a 50-digit reference confirms the Jacobian. The base values A1 are
reproduced exactly from the source files. No unit, weighting or noise-scaling error exists; nothing numerical can
lower w* by an order of magnitude.

The caveats are about reading and design, not numerics:
1. "By large factors for all nuisance sets" is too strong at 500 nm for N0 and N1. On the declared 80-point
   log-spaced design the 1σ ratios σ/w are 1.9–2.8 (N0) and 2.6–3.7 (N1) at 0.01°. With the c-optimal allocation of the
   same 80 measurements in the same band they fall to 0.70–1.01 (N0) and 1.3–1.9 (N1). "Large factors" is correct for
   N ⊇ N2: ≥ 4.2 at 0.01° and ≥ 34 at 0.1° on the declared design, ≥ 2.4 and ≥ 22 with the optimal design.
2. σ is compared with the half-width w(L). Telling the two ends of the admissible range apart needs σ ≤ 2w, which halves
   every margin. At 500 nm, N0, 0.01° the ratio becomes 0.96–1.4.
3. The size of the margin is a property of the illustrative stack, above all of the interface resistances (note 16:
   R1 = 1e-8, R2 = 2e-8 m²K/W). With R1/R2 = 3e-9/5e-9, σ drops ×5 (N0) and ×3.5 (N2). With all favourable changes
   computed (better interfaces, 20 nm Au, 320 points, 0.1 % amplitude noise), N2 at 500 nm reaches σ/w = 1.2–1.7 (0.01°)
   and 2.0–2.9 (0.1°): still not resolvable, but not by "more than an order of magnitude" (section 7, F8).
4. Model-level inputs are untouched: the base τ_R are 12×12×8 event values (probably 2–4 % low, section 6); the 500 nm
   width is a trend bracket; the finite-beam GK rows have no independent reference (section 8, F7).

None of these reverses the falsification test of 00-question.md (N ⊇ N2, σ_φ = 0.1°). On the declared stack and design,
σ/w = 34–50 at 500 nm, 100 at 300 nm, 590–1300 at 100 nm and 4.5e4 at 30 nm; with the optimal design, ≥ 22 at 500 nm;
with every favourable change, ≥ 2.0. H_A survives.

## 1. Independent implementations

| layer | twin (branch production) | `indep.py` (branch check) | this audit |
|---|---|---|---|
| layer physics | global mode-matching system per (k, p), LU with equilibration | impedance recursion upward | ABCD (quadrupole) product top-down, each layer scaled by e^{-γd}, Julia (`aud_model.jl`); unscaled cosh/sinh at 50 digits, mpmath (`aud3_mp_check.py`) |
| Hankel average | trapezoid in x = ln(√η k), step 0.04, x ∈ [−22, 2] | composite Gauss–Legendre in u = √η k on [0, 7] | adaptive Gauss–Kronrod (7, 15), own code, in s = η k² on [0, 120], componentwise rtol 1e-13 on H and on every derivative; tanh-sinh on [0, ∞) at 50 digits |
| derivatives | central differences of ln H, Richardson h = 2e-3 | forward AD (JAX) | dual numbers (ForwardDiff), all 13 thermal/beam parameters in one pass; 1e-20 central difference under the integral at 50 digits |
| Fisher algebra | SVD of the column-normalised whitened Jacobian (`common.fisher`) | the same routine | orthogonal projection on the nuisance span (Householder QR, re-orthogonalised); 256-bit BigFloat inverse; Cholesky of the Jacobi-scaled normal matrix and a 40-digit inverse (Python, `aud4_compare.py`) |
| base values | `a1_base_values.py` | — | own re-assembly of the phono3py chunks and REF blocks (`aud1_base_values.py`) |

Physics re-derived for the audit model (not copied). Cattaneo layer in Hankel space: Q_z' = −(Cp + λ_r k²/g)T,
T' = −g Q_z/λ_z, g = 1 + τ_R p, hence γ² = (Cpg + λ_r k²)/λ_z and Y = λ_z γ/g. Probe average
H = (1/2π)∫k e^{−ηk²} Z dk with η = (w_p² + w_r²)/8. Its steady limit on a half-space,
1/(√(2π) λ √(w_p² + w_r²)), is the standard Gaussian-beam heating (A/(2√π Λ w) for equal radii), so the normalisation is
right. Instrument columns: ∂/∂g = 1, ∂/∂φ0 = i (rad), ∂/∂t_d = −2πif·1e-9 (t_d in ns). Internal checks: scaled vs
unscaled ABCD ≤ 5.6e-16 relative; dual numbers vs central differences (h = 1e-6) at a node agree to the FD error; GK15
integrates x^22 to 4e-17.

## 2. Reproduction of the A4 table (finite beams 10/8 µm, Cattaneo, λ_r = RTA κ_x(L), 80 points, 10 kHz – 200 MHz)

| L | set | σ_φ (°) | branch σ(ln τ_R) | audit σ(ln τ_R) | relative difference |
|---|---|---|---|---|---|
| 500 nm | N0 | 0.01 | 0.154002 | 0.154002 | +5e-12 |
| 500 nm | N1 | 0.01 | 0.205122 | 0.205122 | +6e-11 |
| 500 nm | N2 | 0.01 | 0.340379 | 0.340379 | −2e-10 |
| 500 nm | N3 | 0.01 | 2.25388 | 2.25388 | −1e-9 |
| 500 nm | N2 | 0.1 | 2.73316 | 2.73316 | −1e-10 |
| 500 nm | N5 | 0.1 | 56.6743 | 56.6743 | −3e-9 |
| 500 nm, RTA base | N2 | 0.01 | 0.296775 | 0.296775 | +7e-11 |
| 300 nm | N0 | 0.01 | 0.246161 | 0.246161 | +3e-10 |
| 300 nm | N2 | 0.01 | 0.662085 | 0.662085 | +7e-10 |
| 300 nm | N2 | 0.1 | 5.36486 | 5.36486 | +6e-10 |
| 100 nm | N0 | 0.01 | 0.938393 | 0.938393 | +3e-10 |
| 100 nm | N2 | 0.1 | 47.4065 | 47.4065 | +1e-9 |
| 30 nm | N0 | 0.01 | 5.25262 | 5.25262 | −4e-9 |
| 30 nm | N2 | 0.1 | 1259.692 | 1259.692 | −2e-8 |
| 500 nm, 1–200 MHz | N2 | 0.01 | 0.532389 | 0.532389 | −3e-10 |
| 500 nm, 10 kHz–50 MHz | N1 | 0.01 | 3.55678 | 3.55678 | +2e-9 |

* All 270 finite-beam entries (5 base points × 6 sets × 3 bands × 3 noise levels): max |relative difference| 7.2e-6, at
  30 nm, N5, 10 kHz–50 MHz (σ = 5.6e4, normalised condition number 5e4). N0–N3: ≤ 4.3e-7. Full band, all sets: ≤ 1.6e-6.
* 1D limit (N0–N4, four L): ≤ 5.7e-6. 1D GK rows (τ_ℓ = 0.5 and 2 τ_R, `aud_gk1d.jl`): N0–N3 ≤ 3.5e-7; N4
  (σ ≈ 2e3–2e7, nearly degenerate) ≤ 5.3e-3, i.e. the branch's FD noise amplified by the degeneracy, as A2 states.
* A3 phase changes (exact, nonlinear) are reproduced to the four printed digits: 500 nm, w = ±8 %: 1.603e-3°;
  30 nm, ±2.8 %: 1.554e-5°; max |dφ/d ln τ_R| = 5.55e-4 / 3.403e-3 / 1.225e-2 / 2.004e-2° at 200 / 156 / 83 / 65 MHz
  (30/100/300/500 nm); asymmetry |Δφ(+8 %)| − |Δφ(−8 %)| ≤ 8e-5 relative.
* Order-of-magnitude check by hand. At 500 nm the τ_R term is a reactive correction iωτ_R L/λ_z to the film resistance.
  The transducer heat capacity shunts it (pC_Au d_Au Z_below ≈ 3i at 65 MHz), so a lumped estimate gives |δ ln H| ≈ 2.3e-4
  rad per unit w, i.e. 0.013°, against 0.020° computed. With 18.6 points per decade and sensitivity ∝ f below the peak,
  Σ h² ≈ 11 effective points, so σ_N0 ≈ 1.745e-4/(3.5e-4·√11) ≈ 0.15. The magnitude of the table is physical, not an
  artefact.

## 3. Error budget of the branch numerics

**Jacobian (finite differences).** Twin FD through an own wrapper (`aud2_twin_checks.py`, not `common.py`) against the
dual-number Jacobian:

| | 500 nm | 30 nm |
|---|---|---|
| central h = 0.1 / 0.03 / 0.01 / 1e-3, τ_R column (max-normalised) | 1.72e-3 / 1.55e-4 / 1.72e-5 / 1.72e-7 | 1.66e-3 / 1.50e-4 / 1.66e-5 / 1.67e-7 |
| observed order (h = 0.1 → 0.01) | 2.00 | 2.00 |
| central h = 1e-4 / 1e-5 (roundoff branch) | 3.9e-9 / 3.7e-8 | 1.7e-7 / 1.7e-6 |
| production Richardson h = 2e-3, τ_R column | 6.7e-10 | 2.8e-8 |
| production, other columns (R1, λ_z, C_Au, w-scale) | ≤ 8e-12 | ≤ 7e-11 |
| effect on σ, N0–N5, 0.01 and 0.1° | ≤ 4.3e-9 | ≤ 1.6e-6 |

The A2 rows are confirmed. The 30 nm floor is FD roundoff: |J_τ| ≤ 1.6e-5 there, against an ln H roundoff of ~1e-15
divided by 2h.

**50-digit reference** (8 grid frequencies, 30 and 500 nm). The audit model agrees with it to ≤ 9e-16 in ln H and
≤ 5.5e-15 relative in J_τ at every frequency; the self-check of the 50-digit quadrature is ≤ 1.5e-51. The twin agrees to
2.0–2.2e-14 in ln H (see below) and 1.2e-10 – 5.2e-5 relative in J_τ per frequency. The 5.2e-5 occurs at 30 nm,
10 kHz, where |J_τ| is 3.3e-4 of its maximum: absolute 3e-13, 1.7e-8 of the column maximum, weight 1e-7 in the Fisher
sum.

**Hankel quadrature of the twin** against the audit reference, max |Δ ln H| over 80 frequencies (same at all four L):

| rule | 30 nm | 500 nm |
|---|---|---|
| step 0.32 / 0.16 / 0.08 / 0.04 / 0.02 / 0.01 (x_low −22) | 2.3e-6 / 6.5e-13 / 2.2e-14 / 2.2e-14 / 2.2e-14 / 1.6e-13 | 2.8e-6 / 7.4e-13 / 2.2e-14 / 2.2e-14 / 2.2e-14 / 1.6e-13 |
| x_low −14 / −30 (step 0.04) | 1.5e-12 / 2.2e-14 | 2.0e-12 / 2.2e-14 |
| x_high 1.5 / 3 instead of 2 (rebuilt from the twin's own nodes) | 2.3e-9 / 2.2e-14 | 2.2e-9 / 2.2e-14 |
| twin default per-frequency step | 3.6e-14 | 1.1e-14 |

Convergence in the step is exponential, as claimed (step 0.32 → 0.16: six orders). **Finding F3.** The floors 2.2e-14 at
steps 0.02–0.08 and 1.6e-13 at step 0.01 are not quadrature error. They are a constant real offset of ln|H| (the same at
every frequency; the imaginary parts agree to 3e-16). `trapezoid_nodes` takes the step as x[1] − x[0] of a `linspace`
that starts at −22, and that difference carries a rounding error: h_used/h_true − 1 = −2.13e-14 at steps 0.02–0.08 and
+1.56e-13 at step 0.01. The weight sum is off by −2.1e-14 and +1.6e-13. A constant factor cancels in ln(H₊/H₋) and is
absorbed by the gain, so no σ is affected. The 2.2e-14 "agreement" of twin and `indep.py` in A2 is this offset; the two
independent models agree with the 50-digit value to ≤ 1e-15. Adaptive rule of the audit: rtol 1e-10 against 1e-13
changes σ by ≤ 1.2e-12 (all sets, L, noise levels).

**Fisher algebra and conditioning.** QR projection (audit) against the branch SVD: as in section 2. QR against the
256-bit BigFloat inverse: ≤ 5e-13. Cholesky of the normal matrix against a 40-digit inverse: ≤ 2e-8 at 30 nm (the
normal equations square the condition number) and ≤ 4e-13 at 500 nm. Normalised condition numbers (column-normalised
whitened Jacobian):

| L, σ_φ | N1 | N2 | N3 | N4 | N5 |
|---|---|---|---|---|---|
| 100 nm, 0.1° | 6.96 | 13.0 | 207 | 3.59e3 | 3.62e3 |
| 500 nm, 0.1° | 2.22 | 7.66 | 147 | 210 | 325 |
| 30 nm, 0.01° | 4.64 | 109 | 1.69e3 | 4.85e4 | 4.94e4 |

The branch statement (2–13, 150–210, 210–3.6e3 at 100 and 500 nm) is confirmed. At 30 nm the condition number reaches
4.9e4. That is harmless in double precision but not mentioned in A4.

**Classification: CONVERGED.** The dominant numerical error is the FD roundoff of the τ_R column at 30 nm: ≤ 1.6e-6 in
σ, five or more orders below every modelling sensitivity.

## 4. Parameterisation and the reading w* = σ(ln τ_R)

* All positive parameters are log-parameters; g is the log gain; φ0 is in rad; t_d is in ns with an exact column. The
  audit derives and codes these columns separately and reproduces N1–N5 to ≤ 1.6e-6, so the instrument enters correctly.
* Phase units. The whitening divides Im J (rad) by σ_φ in rad. A degree/radian slip would be visible:
  * reading 0.01 as rad gives σ_N0(500 nm) = 4.56 (×30; the 1 % amplitude then dominates);
  * a Jacobian in degrees against noise in rad gives 0.00269 (÷57);
  * only the correct conversion reproduces the branch value 0.154002 (300 nm: 8.35 / 0.0043 against 0.24616).

  Phase-only / amplitude-only at 500 nm are 0.15407 (0.01°) / 5.320 (1 %), as in A6 (1.54 at 0.1°, 5.32). The phase
  carries ≥ 99.9 % of the N0 information at 0.01°, and 94–97 % of it comes from above 10 MHz.
* **Linear response in τ_R** (noiseless profile, nuisances refitted by Gauss–Newton with dual-number Jacobians,
  `aud_model2.jl`). Δχ²(w) against (w/σ_F)²:

| L, set, σ_φ | σ_F | Δχ²(−2.8 %) / (+2.8 %) | (w/σ_F)², 2.8 % | Δχ²(−8 %) / (+8 %) | (w/σ_F)², 8 % | τ_R = 0 (w = −1): Δχ² / pred | w = +3: Δχ² / pred |
|---|---|---|---|---|---|---|---|
| 500, N0, 0.01 | 0.1540 | 0.0330576 / 0.0330565 | 0.0330570 | 0.269867 / 0.269840 | 0.269853 | 42.19 / 42.16 | 378.9 / 379.5 |
| 500, N1, 0.01 | 0.2051 | 0.0186293 / 0.0186373 | 0.0186333 | 0.152016 / 0.152201 | 0.152109 | 23.59 / 23.77 | 219.0 / 213.9 |
| 500, N2, 0.01 | 0.3404 | 0.00676666 / 0.00676720 | 0.00676693 | 0.0552340 / 0.0552465 | 0.0552402 | 8.620 / 8.631 | 78.04 / 77.68 |
| 500, N2, 0.1 | 2.733 | 1.04942e-4 / 1.04959e-4 | 1.04951e-4 | 8.5654e-4 / 8.5694e-4 | 8.5674e-4 | 0.1335 / 0.1339 | 1.216 / 1.205 |
| 500, N3, 0.01 | 2.254 | 1.54096e-4 / 1.54567e-4 | 1.54332e-4 | 1.25437e-3 / 1.26535e-3 | 1.25985e-3 | 0.1864 / 0.1969 | 2.082 / 1.772 |
| 300, N0, 0.01 | 0.2462 | 0.0129385 / 0.0129380 | 0.0129383 | 0.105624 / 0.105613 | 0.105619 | 16.51 / 16.50 | 148.2 / 148.5 |
| 300, N2, 0.01 | 0.6621 | 1.78852e-3 / 1.78848e-3 | 1.78850e-3 | 0.0146004 / 0.0145996 | 0.0146000 | 2.282 / 2.281 | 20.51 / 20.53 |
| 100, N2, 0.1 | 47.41 | 3.48861e-7 / 3.48842e-7 | 3.48852e-7 | 2.84800e-6 / 2.84754e-6 | 2.84777e-6 | 4.454e-4 / 4.450e-4 | 3.993e-3 / 4.005e-3 |

  At the admissible changes (|w| ≤ 8 %) the profile equals (w/σ_F)² to ≤ 6.1e-4 relative for N0–N2 and ≤ 4.4e-3 for N3.
  It is symmetric in w, not in ln(1 + w). At |w| = 8 % the two signs differ by ≤ 1.2e-3 (N0–N2) and 8.7e-3 (N3); a
  response linear in ln τ_R would give 5.8 % at |w| = 2.8 % and 17 % at |w| = 8 %. So σ(ln τ_R) is the resolvable
  relative change w at 1σ, as the branch reads it.

  Over w ∈ [−1, 3] the deviation is ≤ 2.4 % for N0–N2 (worst N1 at w = 3). N3 departs at large w:
  * −5.3 % at w = −1, where the profile lies below the Fisher value, so the Fisher number is slightly optimistic for
    the Fourier alternative;
  * +5.6 % at w = 1 and +17.5 % at w = 3, where it lies above the Fisher value (conservative).

  The w = −1 values reproduce A7 (42.2, 23.6, 8.62, 0.133, 0.186; 16.5, 2.28; 4.4e-4), and the ±2.8 / ±8 % values
  reproduce A5.
* Linearity of H itself: max |ΔH − w τ_R ∂H/∂τ_R| / max |w τ_R ∂H/∂τ_R| at w = −1 / −0.5 / +1 / +3 is
  2.4e-4 / 1.2e-4 / 2.4e-4 / 7.1e-4 (30 nm), 6.4e-4 / 3.2e-4 / 6.4e-4 / 1.9e-3 (100 nm), 2.5e-3 / 1.3e-3 / 2.5e-3 / 7.7e-3
  (300 nm), 5.6e-3 / 2.8e-3 / 5.7e-3 / 1.75e-2 (500 nm). This confirms "H linear in τ_R over 0 ≤ τ_R ≤ 4τ_R,0" (to
  1.8 %). The quadratic term is the (ωτ_R)² correction, ωτ_R ≤ 0.031 at 200 MHz.

## 5. Frequency grid

σ·√N, audit model, on two kinds of grid over 10 kHz – 200 MHz: N log-spaced points with both endpoints (the branch's
grid), and N points at the centres of N equal cells in ln f:

| L, set, σ_φ | grid | N = 20 | 40 | 80 | 160 | 320 | 640 | 1280 |
|---|---|---|---|---|---|---|---|---|
| 100 nm, N0, 0.01° | endpoints | 8.031 | 8.260 | 8.393 | 8.465 | 8.502 | 8.521 | 8.531 |
| | centres | 8.512 | 8.534 | 8.539 | 8.540 | 8.541 | 8.541 | 8.541 |
| 100 nm, N2, 0.1° | endpoints | 324.1 | 377.3 | 424.0 | 454.8 | 472.5 | 482.1 | 487.0 |
| | centres | 648.5 | 523.5 | 499.5 | 493.9 | 492.5 | 492.2 | 492.1 |
| 500 nm, N0, 0.01° | endpoints | 1.3914 | 1.3832 | 1.3774 | 1.3742 | 1.3725 | 1.3716 | 1.3712 |
| | centres | 1.3748 | 1.3717 | 1.3710 | 1.3708 | 1.3708 | 1.3708 | 1.3708 |
| 500 nm, N2, 0.1° | endpoints | 24.08 | 24.29 | 24.45 | 24.54 | 24.60 | 24.63 | 24.64 |
| | centres | 24.72 | 24.67 | 24.66 | 24.66 | 24.66 | 24.66 | 24.66 |

The two grids bracket a common limit, the continuum information density. The endpoint grid approaches it as O(1/N): for
100 nm, N2, the distances 68 / 37 / 19.5 / 10 / 5 for N = 80 … 1280 halve per doubling. This is the Euler–Maclaurin
end-point term, i.e. the full weight of the 200 MHz point, where the N2 information sits at 100 nm. At N = 80, σ√N of
the declared design differs from the continuum limit by:
* +0.5 % at 500 nm, N0 (the sensitivity peaks at 65 MHz and is lower at the end point);
* −0.9 % at 500 nm, N2, 0.1°;
* −1.7 % at 100 nm, N0;
* −14 to −15 % at 100 nm, N2.

A negative sign means the end point makes the design more informative per point than the continuum density.

The branch's numbers (8.03 … 8.50, 324 … 473, 1.391 … 1.373, 24.08 … 24.60) are reproduced, and its reading is correct.
The drift is a property of the design, not a discretisation error, and σ for the declared 80-point design needs no grid
extrapolation. The independent code confirms the plain sum over independent points (no frequency quadrature weight), so
at fixed per-point noise σ ∝ N^{-1/2}.

## 6. Base values (A1)

* **RTA chunks** (`results_c/rta/m484826, m565630, m646435`, own assembly code). The archived `assembled.json` film
  memory times (x and z, 30 nm – 10 µm) and κ are reproduced to ≤ 2.2e-16 on all three meshes. τ_z(500 nm) =
  25.4507 / 25.5095 / 25.5667 ps and κ_z(500 nm) = 169.455 / 169.445 / 169.344 W/mK. The q_z = 0 plane has weight 1/N_z
  and max |v_z| = 2.0e-29 m/s on every mesh.
* **Convergence wording (F4).** The per-step changes of τ_z(500 nm), +0.231 % and +0.224 %, are confirmed. They do not
  decrease (equal increments for N_x = 48 → 56 → 64), so the three meshes are not in an asymptotic regime and no power law
  fits them. Extrapolation with p = 2 / p = 1 gives 25.75 / 25.97 ps (+0.7 / +1.6 %). "+0.22 %" is a step size, not an
  error bar; a 1–2 % mesh uncertainty is the honest statement.
* **Event operator, REF blocks.** The REF block is identical in all twelve `film3/*_m12128.json` files (spread 0). At
  30/100/300 nm, τ_z = 5.48689 / 11.2264 / 19.6993 ps and κ_z = 62.551 / 115.621 / 167.248 W/mK, as used.
* **500 nm.** The ratio method is reproduced: τ_R = 24.99710 ps, λ_z = 188.5691 W/mK, λ_r = 215.8905 W/mK (RTA κ_x). The
  event(12×12×8)/RTA(64×64×35) τ_z ratio is 1.00965 at 300 nm, 0.93601 at 1000 nm and 0.97772 interpolated at 500 nm.
  The sentence "0.973 at 300 nm and 0.936 at 1000 nm" in A1 misquotes the 300 nm value (the arithmetic is right).
  Alternatives: log-log 24.81, PCHIP 25.14, RTA direct 25.57 ps. A same-mesh estimate (event/RTA on the 12×12×8
  diagonal, 1.0393 at 300 nm and 0.9963 at 1000 nm, interpolated and applied to RTA 64×64×35) gives 26.10 ps (+4.4 %).
* **Mesh bias of the base values (caveat).** RTA on the 12×12×8 diagonal lies −2.0 / −1.9 / −2.9 / −6.0 % below RTA
  64×64×35 at 30 / 100 / 300 / 1000 nm. If the event operator carries a similar bias, the converged event τ_z are about
  2–4 % larger than the base values at 30–500 nm. σ ∝ 1/τ_R for N0–N2 (A6; audit scenario τ_R × 1.044 in section 7), so
  σ would be 2–4 % lower. This is immaterial for the verdict, but it should be stated with E-A4.
* **Admissible widths.** From the result files: three-mesh [0.9727, 1.0281] / [0.9632, 1.0379] / [0.9491, 1.0537]; held
  out on 12×12×8 [0.9730, 1.0276] / [0.9658, 1.0327] / [0.9536, 1.0409] (30/100/300 nm), as in 00-question.md. The 300 nm
  half-width used, 0.052, is the mean of 5.1 % and 5.4 %; with 0.054 the 300 nm ratios fall by 4 %. At 1000 nm the
  held-out 300 nm witnesses give 1.0534 / 0.9354, so the 500 nm bracket 0.055–0.08 is plausible as an inner bound.

## 7. Search for effects that could lower w* by an order of magnitude

σ(ln τ_R) at 500 nm unless stated, on the declared design and noise unless stated (`aud_model2.jl`,
`aud5_design_bound.py`). Factors in brackets are relative to the declared case.

| candidate | effect | conclusion |
|---|---|---|
| degree/radian slip in the phase weighting | ×30 or ÷57 (section 4) | excluded: only the correct conversion reproduces the branch |
| per-point noise against a frequency weight | plain sum over independent points; σ ∝ N^{-1/2} | excluded as an error; 80 → 320 points halves σ (N2, 0.01°: 0.340 → 0.171) |
| missing information | phase and amplitude both used; phase carries ≥ 99.9 % of the N0 information at 0.01° | nothing omitted within the data model |
| allocation of the same 80 measurements inside the band (c-optimal design on the 80 grid frequencies as candidates; equivalence gap ≤ 1e-6, so φ* ≥ φ(1 − gap) certifies σ to 5e-7; a continuous optimum between grid points gains ≤ 0.12 % for N0) | gain 1.3–2.9: N0 0.0558, N1 0.1045, N2 0.195 at 0.01°; N2 1.76 at 0.1°. 300 nm: N0 0.0912, N2 0.379. 100 nm: N2 2.22 | ≤ 3×; N0/N1 at 500 nm reach w (F1) |
| amplitude noise 0.1 % instead of 1 % | N2: 0.273 (0.01°), 0.645 (0.1°, ×4.2) | large at 0.1°, small at 0.01° |
| interface resistances R1/R2 = 3e-9/5e-9 m²K/W (note 16: 1e-8/2e-8) | N0 0.030 (×5.1), N2 0.098 (×3.5) at 0.01°; N2 0.90 (×3.0) at 0.1° | the largest single lever |
| R1 = R2 = 1e-9 m²K/W | N0 0.0134 (×11.5); N2 0.098 (0.01°), 0.84 (0.1°) | N0 gains an order of magnitude; N2 saturates near ×3.5 because the R's are free |
| 20 nm Au transducer | N0 0.056 (×2.8), N2 0.215 (×1.6) at 0.01°; N2 1.95 at 0.1° | ≤ 3× |
| beams 3/2.4 µm or 30/24 µm | N2 0.70 or 0.30 (0.01°) | smaller beams are worse |
| f_max = 1 GHz at the same point density (model-internal) | N2 0.178 (0.01°), 1.69 (0.1°) | ≤ 2× |
| calibration priors φ0 ± 0.01°, t_d ± 1 ps (and R's ± 1 %) | N2 at 0.1°: 2.02 (1.98) | ≤ 1.4× |
| base τ_R +4.4 % (mesh-corrected, section 6) | σ × 0.958 | negligible |
| all favourable at once (R 3e-9/5e-9, 20 nm Au, 320 points, σ_A 0.1 %) | N0 0.0085, N1 0.083, N2 0.093 (0.01°); N2 0.159 (0.1°). 300 nm: N2 0.346 / 0.509 | N2/w = 1.2–1.7 (0.01°) and 2.0–2.9 (0.1°) at 500 nm; 6.7 / 9.8 at 300 nm |
| averaging n independent sweeps (σ ∝ n^{-1/2}), declared stack | 2σ ≤ w needs, at 500 nm, n = 15–31 (N0, 0.01°), 72–153 (N2, 0.01°), 4.7e3–9.9e3 (N2, 0.1°); at 300 nm, 600–650 (N2, 0.01°) | as in the A-C synthesis |

**Reading.** Within the declared scenario (stack of note 16, 80 log-spaced points, σ_A = 1 %), nothing numerical or
statistical lowers w* by a factor of ten. That factor is available only through the experiment or the stack:
* ≥ 100 repetitions;
* σ_φ ≤ 0.001°;
* interfaces several times better than the illustrative values (N0 only);
* combinations of these.

Even with every favourable change computed here, N ⊇ N2 stays above w(L) at 1σ: at 500 nm, 1.2–1.7 at 0.01° and
2.0–2.9 at 0.1°.

The A-C synthesis says "H_A survives with a margin of more than an order of magnitude at every thickness". That is
specific to the declared stack. At 500 nm, N2, 0.1° the margin shrinks from 34–50 to 11–16 with R1/R2 = 3e-9/5e-9, and
to 2.0–2.9 with all favourable changes. The A-C list of "what would change the verdict" covers phase noise, the size of
τ_R and the GHz band, but not the interface resistances; they should be added (F8).

## 8. Findings to record

* **F1 (wording of E-A3).** "Exceeds by large factors for all nuisance sets N0–N5" holds for N ⊇ N2, not for N0/N1 at
  500 nm (σ/w = 1.9–2.8 for N0, 2.6–3.7 for N1 at 0.01°; 0.70–1.01 and 1.3–1.9 with a c-optimal allocation of the same
  80 points). Suggested wording: "w*₁σ exceeds w(L) for every declared set; by ≥ 4 (0.01°) and ≥ 34 (0.1°) for
  N ⊇ N2. In the idealised N0/N1 cases at 500 nm it is within a factor 2–4 on the declared design and comparable to w
  with an optimised design."
* **F2 (half-width).** σ is compared with the half-width. Discriminating the ends of the admissible range needs σ ≤ 2w,
  which halves all ratios. At 500 nm, N0, 0.01°: 1.4 / 0.96 against 2w. State which comparison is meant.
* **F3 (quadrature floor).** The constant offset −2.13e-14 (+1.56e-13 at step 0.01) comes from the step x[1] − x[0] in
  `trapezoid_nodes`. It is not quadrature error and has no effect on σ. A2 should not present 2.2e-14 / 1.6e-13 as
  quadrature errors. (Using h = (x_high − x_low)/(n − 1) would remove it; not needed for any result.)
* **F4 (A1 text).** The "0.973 at 300 nm" ratio is a misquote (1.0097). The RTA mesh change "+0.22 %" is a step change,
  not an error estimate (increments not decreasing; extrapolation +0.7 … +1.6 %). The 12×12×8 base values are probably
  2–4 % low (σ 2–4 % high).
* **F5 (independence of indep.py).** `indep.py` imports neither the twin nor any of its functions. It uses only the
  constants of `common.py`, through a function-level import. However, the A2/A4 cross-check passed its Jacobian through
  the branch's own `weighted`/`fisher` routines, so it did not test the noise weighting or the marginal-variance algebra.
  This audit closes that gap: separate QR-projection and Cholesky/40-digit routines agree with the branch to ≤ 4.3e-7
  (N0–N3).
* **F6 (conditioning).** The A4 condition numbers are confirmed. At 30 nm the normalised condition number reaches 4.9e4
  (N5); harmless (QR against 256-bit: 5e-13).
* **F7 (not covered).** The finite-beam GK rows (isotropic GK film, slip walls, N4 with τ_ℓ free) have only the branch's
  Richardson self-consistency; no independent implementation exists. The verdict does not depend on them; they should be
  labelled as such. The in-plane-τ variant of A2 and the λ_r = λ_z variant were not recomputed.
* **F8 (missing sensitivity: interface resistances).** The τ_R signal is the reactive part of the film resistance, so
  its relative size scales with the share of the film in the resistance below the transducer. R1 = 1e-8 and
  R2 = 2e-8 m²K/W (note 16) make that share small. Lower resistances raise the sensitivity strongly:
  * N0 at 500 nm: ×5.1 with 3e-9/5e-9, ×11.5 with 1e-9/1e-9;
  * N2: about ×3.5 in both cases;
  * N2 at 0.1°: ×3.0.

  The A-C list of what would change the verdict should include the interface resistances. A statement such as "margin
  of more than an order of magnitude" should be tied to the note 16 values.

## 9. Status

E-A3 can be kept as NUMERICALLY DEMONSTRATED (model-conditional) after the numerical audit, with the F1/F2 wording. E-A4
(base values) is NUMERICALLY SUPPORTED as stated, with the mesh caveat of section 6. E-A5 (τ_R > 0 against Fourier only
at 300–500 nm, N0–N2, 0.01°) is confirmed by the profiles of section 4 (w = −1 column).

Numerical classification: CONVERGED.
* No element is UNDER-RESOLVED, TOLERANCE-LIMITED or BOUNDARY-CONTAMINATED. Domain truncation: x_low −14 → 2e-12,
  x_high 1.5 → 2e-9; the production values −22 and 2 sit below the 2e-14 floor.
* The only floating-point-limited quantities are the twin's constant normalisation (F3, 2e-14) and the FD τ_R column at
  30 nm (2.8e-8). Neither appears in σ above 1.6e-6.

## 10. Reproduction

Research Python (`ultimate-research`), Julia 1.13.1 with `--project=C:/Users/Koussay/ResearchLab/envs/julia-research`,
one thread each, at most two processes at a time.

| step | command (in `experiments/A_audit/`) | time | output |
|---|---|---|---|
| base values | `python aud1_base_values.py` | 2 min | out/aud1_base_values.json |
| inputs | `python aud0_inputs.py` | 1 s | out/aud_inputs.json, out/aud_inputs.jl |
| independent model, all tables | `julia --project=... aud_model.jl` | 77 s | out/aud_julia_core.json |
| twin checks (quadrature, FD Jacobian) | `python aud2_twin_checks.py` | 105 s | out/aud2_twin_checks.json |
| 50-digit reference | `python aud3_mp_check.py` | 2.5 min | out/aud3_mp_check.json |
| comparisons, own Fisher algebra | `python aud4_compare.py` | 10 s | out/aud4_compare.json |
| profiles, linearity, grid, units, scenarios | `julia --project=... aud_model2.jl` | 150 s | out/aud_julia_part2.json |
| c-optimal design bound | `python aud5_design_bound.py` | 1 min | out/aud5_design_bound.json |
| 1D GK rows | `julia --project=... aud_gk1d.jl` | 30 s | out/aud_gk1d.json |
