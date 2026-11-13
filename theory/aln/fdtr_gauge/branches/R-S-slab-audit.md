# Branch R-S-A — independent numerical audit of the full-operator kinetic slab (R-S, R-S1…R-S6)

13 November 2026. Object: `07-experiments.md` R-S, R-S1…R-S6 and the R-S synthesis; scripts `experiments/R_slab/`
(rs_core.py, rs1_validate.py, rs2_slab.py, rs2b_fits.py, rs2c_interp_check.py, rs3_margins.py, rs5_trend.py,
rs6_proxy.py, rs7_band.py) and their outputs `experiments/R_slab/out/`; the earlier RTA-slab audit
`branches/R-A-kinetic-audit.md` (au5/au6/au8). Audit scripts `experiments/R_slab_audit/sa*.py`, raw outputs
`experiments/R_slab_audit/out/` (list in §7). Research Python, two BLAS threads, one process at a time, peak working
set ≤ 1.4 GB (7×7×5 runs ≤ 1.0 GB; 9×9×6 run 1.0 GB under a free-memory guard: abort before any heavy step if free
RAM < 1.5 GB or own private memory > 1.9 GB); the optimisation job and the Quantum ESPRESSO run on the same machine
were not touched. Every number is SYNTHETIC and conditional on the event operator, the black diffuse walls and the
twin's 1D FDTR description (02-assumptions A1–A7).

**Claim audited.** A cross-plane slab with the full orbit-tied event operator (7×7×5) and diffuse black walls, solved
modally (A_e φ = λ²Bφ after the parity reduction q_z → −q_z and Schur elimination of the q_z = 0 plane), gives the
FDTR-visible δ = τ_R − τ_ℓ of −31.70 ps (300 nm) and −61.00 ps (500 nm) for REF; a static-data half-range of 5.2 ps
(300 nm) and 11.1–11.3 ps (500 nm) over the 12 witness tensors; margins σ/half-range at 500 nm, 0.01° of N1 0.59,
N2 1.01, N2lam 2.34, N3 5.97 ("not resolvable once λ_z is free"). Stated validation: modal vs adding-doubling 6.6e-12,
energy balance ≤ 8e-10, AAA interpolation effect ≤ 1e-4 ps, RTA reduction vs au6 5.9e-11.

## 0. Verdicts

| # | claim | verdict | decisive numbers |
|---|---|---|---|
| 1 | modal solution correct and complete (parity map, Schur elimination, diffusive mode, choice of decaying modes, conditioning) | **RELIABLE — CONVERGED** | an independent solver (no parity reduction, no Schur complement, no eigendecomposition, no invariant imbedding; z-discretised exact characteristics, cubic source, GMRES) reproduces the stored two-ports of REF on 7×7×5 at 6 frequencies (10 kHz–200 MHz) × 2 thicknesses to **≤ 1.7e-9** after Richardson extrapolation, with observed order **3.97–3.98**; the frequency-dependent part Im Y11 to ≤ 2.7e-9 of itself; the range-setting witness cal3w min to ≤ 1.6e-9; the RTA reduction to ≤ 2.5e-9; on 5×5×3 to ≤ 5.6e-11. cond(Φ) = 9.1–17, cond(P±) = 3.7–14 (column-normalised); the diffusive-mode polish changes Y by ≤ 1.2e-10 |
| 2 | REF δ = −31.70 / −61.00 ps (N3; N2lam −31.84 / −62.70) | **CONVERGED for the 7×7×5 model; UNDER-RESOLVED in the q-mesh** | independent fits −31.6953 / −61.0035 (N3), −31.8436 / −62.6951 ps (N2lam); single minimum of the misfit profile over δ ∈ [−400, +150] ps; numerical error of δ < 1e-4 ps (δ moves 0.021 ps per 1e-6 relative two-port error). Mesh, 500 nm (N3): −37.63 (5×5×3) → −61.00 (7×7×5) → **−67.36 ps (9×9×6, new)**; N2lam −40.60 → −62.70 → −69.35; steps shrink ×0.27–0.30 (heuristic limit ≈ −70 / −72 ps). 300 nm: −28.4 → −31.7 (9×9×6 not run) |
| 3 | half-range 5.2 ps (300 nm), 11.1–11.3 ps (500 nm) | **RELIABLE WITH CAVEATS** (arithmetic and interpolation verified; inner bounds; not mesh-converged) | 5.255 / 5.218 and 11.324 / 11.084 ps (N2lam / N3) recomputed; AAA on the range-setting witness against exact 80-frequency data: ≤ 1.6e-4 ps; leave-out cross-validation on all 13 tensors ≤ 0.008 ps. 500 nm, other meshes: 40.6 ps on 5×5×3 (film witnesses admissible there), **8.73 / 8.58 ps on 9×9×6** (the two end-setting witnesses; proxy 9.42 / 9.07) |
| 4 | margins 0.59 / 1.01 / 2.34 / 5.97 (500 nm, 0.01°); "not resolvable once λ_z is free" | **numbers RELIABLE**; conclusion **RELIABLE WITH CAVEATS** (7×7×5 and 9×9×6, these witnesses); **N2 = 1.01 INCONCLUSIVE** on 7×7×5 (at the threshold; 1.17 on 9×9×6) | σ = 6.73 / 11.41 / 26.51 / 66.22 ps, margins 0.594 / 1.008 / 2.341 / 5.974 (Cauchy-integral Fisher, own stack and fits); direct test 1.942 (N2lam), 5.149 (N3). Mesh: same extraction gives 0.55 / 1.47 on 5×5×3 (resolvable for N2lam), 2.34 / 5.97 on 7×7×5, **2.81 / 7.50 on 9×9×6**: the verdict flips between 5×5×3 and 7×7×5 and strengthens from 7×7×5 to 9×9×6. A flip would need a visible half-range ×2.34 / ×5.97 (7×7×5) or ×2.81 / ×7.50 (9×9×6): δ-optimised witnesses are the open risk |
| 5 | stated validation figures | **CONFIRMED by an independent route** (the stated self-checks share the Schur complement, and the doubling check also the parity blocks, with the production solver) | see §1–§3 |

**Overall: RELIABLE WITH CAVEATS.** The modal slab solution is correct and complete; nothing numerical in the chain
solver → two-port → AAA → GK fit → Fisher moves δ, the half-ranges or the margins by more than 2e-4 ps or 1e-3
relative. On 7×7×5 the conclusion "the static-data range of the visible δ is not resolvable at 0.01° once λ_z is free"
is reproduced (N2lam 2.34, N3 5.97; at 0.1° every set ≥ 5.6). A new full-operator slab on 9×9×6 (independent solver;
the paired reduction does not apply there) gives REF δ = −69.35 / −67.36 ps (N2lam / N3, 500 nm), a range contracted
by ×0.77 (half-range 8.73 / 8.58 ps from the two end-setting witnesses) and margins 2.81 / 7.50: the conclusion
strengthens with refinement. Limitations, all outside the numerics of the model: (i) the mesh — the same extraction
on 5×5×3 resolves the range (0.55), so the verdict rests on a three-mesh trend (0.55 → 2.34 → 2.81) whose last step
is favourable, not on a converged value; the REF δ itself is not converged (+10 % from 7×7×5 to 9×9×6, steps
shrinking); (ii) the ranges are inner bounds of witnesses not optimised for δ (a flip needs ×2.3–2.8 for N2lam);
(iii) N2 = 1.01 on 7×7×5 is not distinguishable from the threshold (1.17 on 9×9×6); (iv) the 9×9×6 data carry a small
non-reciprocity from the q_z = 1/2 plane (≤ 8.5e-4 in Y12 + Y21 at 200 MHz; 0.025 ps in δ). Status proposal: "not
resolvable with λ_z free (λ_z free, C known or free; 500 nm, 0.01°; black walls; 12 non-δ-optimised witnesses)":
NUMERICALLY SUPPORTED on 7×7×5 and 9×9×6, mesh trend favourable from the second refinement on, not converged;
"N2 not resolvable": UNRESOLVED on 7×7×5 (marginal), supported on 9×9×6 (1.17).

## 1. The modal solver: premises and internals (task 1) — `sa3_modal_diag.py`, `out/sa3_modal_diag.json`

* **Parity map** (own operator assembly, §2): grid map 15 is the reflection q_z → −q_z; frequencies symmetric to
  7.1e-15, v_z antisymmetric to 2.1e-14 (relative to max |v_z|), the assembled operator invariant under the mode
  permutation to 2.0e-15; the map fixes exactly 585 modes, all with v_z = 0 (q_z = 0 plane), none moving. The
  reduction is exact on 7×7×5 (and on 5×5×3).
* **Where it fails.** On 9×9×6 (N_z = 6) the zone-boundary plane q_z = 1/2 is mapped onto itself. Its 972 modes are all
  moving (|v_z| up to 9064 m/s, median 1994 m/s), in 482 degenerate pairs with v_z = ±v (|v₁ + v₂| ≤ 1.4e-3 m/s) and
  8 modes not resolved into pairs at a 1e-6 THz tolerance. The paired reduction needs v_z(Pi) = −v_z(i) on self-mapped modes, which is false
  there (the ± assignment within a degenerate pair is a gauge choice that the stored map does not follow). This is
  why no 9×9×6 slab exists in R-S; it is not a sign of an error on 7×7×5. The solver of §2 needs no pairing.
* **Schur elimination of the v_z = 0 modes** f₀ = −(p + C₀₀)⁻¹C₀±f± is exact algebra; p + C₀₀ is non-singular because
  ker C = span{e} has components on moving modes. It is common to all three of the author's routes (paired, unpaired,
  doubling); the solver of §2 keeps the local modes explicitly (pointwise algebraic relation inside the iteration),
  so the agreement of §2 tests it.
* **Diffusive mode.** At 100 Hz the unpolished geev eigenvalue is μ_d = λ_d² = 883.7 + 5.584e6 i m⁻² (real part wrong
  by a factor 220, relative error 1.6e-4 of |μ_d|, inside the backward-stability bound ε_mach μ_max/|μ_d| ≈ 2.6e-3); the polished one, 4.02 + 5.584e6 i m⁻²,
  agrees with the closed-form expectation pC/K(1 + δp) to 3.6e-9 (2.3e-8 at 10 kHz). Its effect on the two-port is small because λ_d d ≪ 1:
  polished vs unpolished Y differ by ≤ 8.8e-11 (100 Hz), ≤ 1.2e-10 (10 kHz–200 MHz); on Im Y11 by 3.2e-8 (100 Hz),
  5e-10 (10 kHz), ≤ 3e-11 above 1 MHz. The polish is correct but not needed in the FDTR band.
* **Choice of decaying modes.** Each μ gives Re λ > 0 and both decaying exponentials e^{−λz}, e^{−λ(d−z)} are used, so
  no growing exponential appears; the solution space s'' = Ms is 2n-dimensional and complete when M is
  diagonalisable (μ ≠ 0 for p ≠ 0). The even/odd combination of the wall conditions is algebraically equivalent to
  the four-block system. Verified by §2 (no modal basis there).
* **Conditioning** (REF 7×7×5, 100 Hz, 10 kHz, 1.5 MHz, 200 MHz): column-normalised eigenvector matrix Φ:
  κ₂ = 9.06 / 9.06 / 9.07 / 17.1; boundary matrices P+ / P− (column-normalised): 6.5 / 3.7 (300 nm) and 6.0 / 3.8
  (500 nm) up to 1.5 MHz, 14.2 / 13.8 and 14.1 / 13.9 at 200 MHz; spread of λ 2.4e3 … 8.1e9 m⁻¹. The stored
  "pivot ratio" of P− (up to 1.8e4) reflects column scaling (the diffusive column is O(λ_d d)), not ill-conditioning.

## 2. Independent solver and convergence (task 1) — `sa_core.py`, `sa0*.py`, `sa1_verify775.py`

* **Method** (shares only the raw event list, the amplitudes and the symmetry maps with R-S). Own operator assembly
  (event vectors, orbit tying by the minimum image index, dense C = 4π·10¹² Σ g a aᵀ): bulk K = 272.3333751891106
  W/mK and δ_bulk = −1151.5698615 ps, identical to rs1 B5 (relative 0 and 6e-16); 23 674 orbits, 0 missing images;
  C e = 0 to 6.5e-16. Slab: V f′ + (p + C)f = 0 on all 2937 modes (no parity, local modes explicit), black walls;
  splitting C = C_R − K′ with C_R the energy-conserving RTA of diag C; inner RTA solve by exact characteristic
  integration of every mode with the total source interpolated by cubic Lagrange polynomials on a 4-node stencil
  (exact exponential moments ψ_j), the scalar local temperature from a dense (N+1)² system; outer GMRES on the nodal
  values of f (rtol 2e-12; true residual checked). Graded z mesh h(s) = min(h_min + γs, h_max) (0.02 nm, 0.08,
  4 nm; N = 202 / 254 at 300 / 500 nm), bisection levels 0 and 1; Richardson R = (16 Y₁ − Y₀)/15. A second-order
  (linear-source) variant serves as an independent discretisation.
* **Scheme checks.** Grey two-stream slab (closed form): fourth order 1.1e-8 → 6.8e-10 → 4.2e-11 → 2.7e-12 (×16 per
  bisection), second order ×4 per bisection. On 5×5×3 (REF): second order observed p = 1.83–1.95 (sa0b); interior
  h_max and grading γ control the error, the wall cell size does not.
* **5×5×3, fourth order, levels 0/1/2** (error relative to max|Y| of the stored modal two-port):

| case | level 0 | level 1 | level 2 | observed order | Richardson (1, 2) |
|---|---|---|---|---|---|
| 10 kHz, 300 nm | 6.1e-7 | 3.9e-8 | 2.5e-9 | 3.97 | 2.5e-11 |
| 10 kHz, 500 nm | 7.5e-7 | 4.8e-8 | 3.0e-9 | 3.97 | 3.3e-11 |
| 200 MHz, 300 nm | 9.0e-7 | 5.8e-8 | 3.6e-9 | 3.97 | 4.7e-11 |
| 200 MHz, 500 nm | 1.0e-6 | 6.5e-8 | 4.1e-9 | 3.96 | 5.6e-11 |

* **7×7×5, REF, fourth order, levels 0/1** (12 cases, 10 kHz, 123 kHz, 1.51 MHz, 18.5 MHz, 64.7 MHz, 200 MHz ×
  300 / 500 nm): level-0 error 6.4e-7 … 1.10e-6, level-1 error 4.1e-8 … 7.0e-8, error ratio 15.63–15.73 (observed
  order 3.97–3.98), **Richardson error 7.5e-10 … 1.72e-9** (largest at 200 MHz, 500 nm); energy temperatures
  Th ≤ 2.3e-10; Im Y11 and Im Y12 relative to themselves ≤ 2.7e-9 and ≤ 2.9e-9 (the δ-carrying part); GMRES 28–29
  iterations (level 0), 17–18 (level 1, warm start); reciprocity and mirror symmetry of the independent solution
  ≤ 1.6e-14. Second-order variant at level 1: 4.5e-5 … 3.0e-4 (consistent with its order).
* **Witness cal3w_m3L100min** (sets the low end of the range), 500 nm, 0.91 MHz and 200 MHz: 1.0e-9 and 1.6e-9.
* **RTA reduction** (inner solver alone with D = diag C and K′ = 0, i.e. exactly R − (Re)(Re)ᵀ/(eᵀRe)) against
  rs2_m775_film2_rta_all_REF.npz: ≤ 2.5e-9 (12 cases). A first version of this test used the diagonal of the
  rank-one-corrected matrix as the relaxation rates; that different operator gave a constant 1.6e-4 offset, which
  identified the slip (corrected, `rta_vs_modal_v2`).
* **Reading.** Two discretisations of a formulation that shares nothing with the modal route converge at their
  design orders to the stored two-ports; the remaining differences (≤ 1.7e-9 on 7×7×5, ≤ 6e-11 on 5×5×3) are at the
  level of the extrapolated discretisation error. The stated self-checks (unpaired 1.4e-11, doubling 6.6e-12, energy
  balance ≤ 8e-10) are therefore confirmed by a route that also tests the Schur elimination and the parity pairing.

## 3. Extraction and GK fit (task 2) — `sa_fdtr.py`, `sa2_fits.py`, `sa5_aaa_witness.py`

Own stack (quadrupole products), own wall-port coupling (cross-check against au6_slab.Z_stack_slab 4.1e-14 in ln H),
fits by variable projection (instrument terms projected out exactly) and Levenberg–Marquardt from 48–60 starts.

* **REF δ_fit (ps)**, wall port, note-16 stack (author's value in brackets):

| | N2lam 0.01° | N3 0.01° | N2lam 0.1° | N3 0.1° | Δχ² (0.01°) |
|---|---|---|---|---|---|
| 300 nm | −31.8436 (−31.84) | −31.6953 (−31.70) | −31.7709 (−31.77) | −31.6208 (−31.62) | 1.3e-5 / 7.2e-6 |
| 500 nm | −62.6951 (−62.70) | −61.0035 (−61.00) | −62.2222 (−62.22) | −60.6070 (−60.61) | 2.7e-3 / 1.9e-3 |

  All starts reach the minimum (spread ≤ 3e-4 ps) except 300 nm N3 0.1°, where 20 of 60 starts stop in a nuisance-space
  local minimum (Δχ² 1.06); the global minimum is the published one.
* **Global minimum.** Misfit profile χ²(δ) (nuisances re-minimised at each δ, 0.01°) on δ ∈ [−400, +150] ps: one local
  minimum (500 nm: −62.5 N2lam, −60.0 N3 on the 2.5-ps grid; 300 nm: −32.5), Δχ² at the ends 138 / 72 (N2lam) and
  10.3 / 16.7 (N3) at 500 nm. Curvature of the profile: σ = 25.9 ps (N2lam) and 66.0 ps (N3), against the Fisher
  values 26.5 and 66.2 (the profile varies τ_ℓ at fixed τ_R, Fisher τ_R at fixed τ_ℓ).
* **Frequency set** (500 nm, 0.01°, N2lam / N3): 40 even points −62.73 / −61.26, 40 odd −62.68 / −60.85, every fourth
  point −62.89 / −61.92. **Band**: 10 kHz–50 MHz −64.19 / −63.66, 10 kHz–10 MHz −65.49 / −65.65, 1–200 MHz
  −69.59 / −67.03, 10–200 MHz −77.01 / −72.76 (rs7 reproduced). The visible δ is a band-effective first-order
  parameter (the slab is not a first-order GK film); in-band resampling changes it by ≤ 0.9 ps, the band choice by up
  to 19 %.
* **τ_R of the base** × 0.5 / × 2: −62.99 / −62.13 (N2lam), −61.06 / −60.93 (N3): the data fix δ only.
* **Sensitivity to two-port errors.** A smooth relative perturbation ε(1 + 0.5i)g(f) of every entry moves δ_fit
  (N2lam, 500 nm) by 0.016 / 0.021 / 0.005 ps for ε = 1e-6 (g = cos, ramp, quartic in ln f), linearly in ε. With the
  verified two-port accuracy (§2) the numerical error of every δ_fit is < 1e-4 ps.
* **Interpolation.** AAA from the 24 stored frequencies: REF against exact 80-frequency data 2.7e-9 in Y, δ changes
  ≤ 2.3e-5 ps; **the range-setting witness cal3w min against exact 80-frequency data** (modal solver, verified in §2):
  ≤ 8.1e-9 in Y, ≤ 3.2e-9 in ln H, δ changes −6.7e-5 / −1.6e-4 ps (300 nm), +2.6e-6 / −1.5e-4 ps (500 nm). Leave-out
  cross-validation (every second in-band support point removed, all 13 tensors): ≤ 7.5e-6 at the removed points, δ
  changes ≤ 0.008 ps. A cruder interpolant (cubic spline in ln f of Y_W/Y_REF) fails its own validation (1e-3 in Y on
  the wall variants) and moves witness δ by 0.3–0.9 ps: the accuracy of AAA is needed and present.

## 4. Margins recomputed (task 3) — `sa2_fits.py` stage `wit`, `out/sa2_fits.json`

Fisher σ at the REF slab-fit base point by Cauchy-integral differentiation (16 points on a circle in the complex
parameter plane; 24 points at half the radius agree to 2.5e-10), marginal σ along τ_R at fixed τ_ℓ by orthogonal
projection; half-ranges from my own fits of the 13 AAA-interpolated two-ports.

| 500 nm | N1 | N2 | N2lam | N3 |
|---|---|---|---|---|
| σ(δ) 0.01° (ps) | 6.73 (6.732) | 11.41 (11.413) | 26.51 (26.512) | 66.22 (66.217) |
| margin 0.01° | **0.594** | **1.008** | **2.341** | **5.974** |
| margin 0.03° / 0.1° | 1.78 / 5.65 | 2.95 / 8.03 | 6.53 / 14.24 | 16.7 / 41.9 |
| direct test 1/√max Δχ² (0.01° / 0.1°) | 0.024 / 0.23 | 0.058 / 0.55 | 1.942 / 12.5 | 5.149 / 37.9 |

300 nm, 0.01°: 2.539 / 3.757 / 7.225 / 14.054; 0.1°: 23.8 / 30.1 / 46.8 / 93.6; direct 6.448 (N2lam), 12.99 (N3).
Half-ranges: 11.324 (N2lam) and 11.084 ps (N3) at 500 nm, ends set by orb 0.1 30 nm max (+10.35 ps) and cal3w min
(−12.30 ps); 5.255 / 5.218 ps at 300 nm. Every published margin and direct-test value is reproduced to the printed
digits. N1/N2 margins combine σ at known λ_z with the λ_z-free range (campaign convention); the direct test shows that
under N1/N2 the witnesses are distinguished through the interior conductivity (0.02–0.06), not through δ.

## 5. Mesh dependence and δ-optimised witnesses (task 4)

### 5.1 5×5×3 → 7×7×5 with the same extraction (`sa2_fits.py` stage `mesh`)

| | REF δ N2lam / N3 (ps) | half-range N2lam / N3 (ps) | σ N2lam / N3 (ps) | margin N2lam / N3 |
|---|---|---|---|---|
| 5×5×3, REF + 10 film witnesses (admissible on 5×5×3), 500 nm | −40.60 / −37.63 | 40.63 / 39.90 | 22.40 / 58.80 | **0.55 / 1.47** |
| 7×7×5, REF + the 2 film witnesses run there, 500 nm | −62.70 / −61.00 | 10.67 / 10.92 | 26.51 / 66.22 | 2.49 / 6.07 |
| 7×7×5, REF + 12 campaign witnesses, 500 nm | −62.70 / −61.00 | 11.32 / 11.08 | 26.51 / 66.22 | 2.34 / 5.97 |
| 5×5×3, 300 nm | −28.67 / −28.42 | 23.31 / 23.22 | 32.77 / 65.09 | 1.41 / 2.80 |
| 7×7×5, 300 nm, 2 film witnesses / 12 campaign witnesses | −31.84 / −31.70 | 3.92 / 3.95; 5.25 / 5.22 | 37.96 / 73.34 | 9.69 / 18.6; 7.22 / 14.05 |

The witness families differ between the rows (the film witnesses violate the 7×7×5 static data), so the rows are a
probe of the trend, not a convergence sequence. They show that the verdict is decided by the mesh: on 5×5×3 the
same extraction would call the 500 nm range resolvable with λ_z free (0.55); on 7×7×5 it is not (2.3–2.5). σ changes
by +18 % (N2lam); for the two tensors run on both meshes the visible half-range falls 34.6 → 10.7 ps (N2lam) and
34.1 → 10.9 ps (N3), ×0.31 (rs5, reproduced here for 7×7×5).

### 5.2 RTA reading, 7×7×5 → 9×9×6 (au8, binned RTA diagonal)

REF visible δ (N3) −84.41 → −80.88 ps (N2lam −87.09 → −83.26) while the regularised δ_M grows −211.5 → −232.5 ps; range
transfer 0.315 → 0.264; the end-setting witness orb 0.1 30 nm max keeps its visible excursion (−10.0 % → −9.9 %).
In the RTA reading the visible quantities are nearly mesh-stable between the two finer event meshes.

### 5.3 Full-operator slab on 9×9×6 (`sa4_m996.py`, `out/sa4_m996.json`) — new

* **Set-up.** Independent solver of §2 (no pairing needed), fourth order, level-0 mesh (N = 254 at 500 nm), 16 support
  frequencies (1 kHz, every sixth FDTR frequency, 200 MHz) and AAA to the 80 FDTR frequencies (every-sixth support
  costs ≤ 2e-4 ps in δ on 7×7×5, rs2c); 500 nm only. A second, coarser mesh (γ = 0.12, h_max = 6 nm, N = 178) changes
  the two-port by 2.7e-6 (10 kHz) and 3.9e-6 (200 MHz): level-0 error ≲ 3e-6, i.e. ≲ 0.06 ps in δ (§3 sensitivity).
  2 min per frequency, peak working set 1.0 GB, free RAM ≥ 2.2 GB at every guard check.
* **Operator** (own assembly): 95 384 orbits, 0 missing images, C e = 0 to 1.1e-15; regularised C + QDQ at 500 nm:
  τ_R = 23.697, τ_ℓ = 217.897, δ_reg = −194.20067 ps (au8: −194.20067); bulk K = 256.90 W/mK, δ_bulk = −2352 ps.
* **Data issue of the q_z = 1/2 plane.** The mode set carries a net equilibrium current Σ v_z c = 27.6 W m⁻²K⁻¹
  (|ê·b|/|b| = 4.73e-9), all of it from that plane, and the discrete model has no exact v_z-reversal symmetry there:
  |Y11 + Y22|/|Y11| = 5.1–6.5e-7 at all frequencies, |Y12 + Y21|/|Y11| = 6.8e-7 (1 kHz) rising to 8.5e-4 (200 MHz).
  A time-reversal-invariant crystal cannot be non-reciprocal; this is the gauge of the degenerate pairs noted in R-A1,
  not a solver effect (the same solver gives ≤ 1.6e-14 on 7×7×5). Effect on the visible δ: replacing the two-port by
  its reciprocal, mirror-symmetric part moves δ by 0.025 (N2lam) and 0.018 ps (N3).
* **REF** (500 nm, wall port, note-16, 0.01°): **δ_fit = −69.35 ps (N2lam), −67.36 ps (N3)**; 0.1°: −68.82 / −66.92;
  Δχ² = 4.2e-3 / 3.1e-3; interior λ_z,fit = 236.2 W/mK (0.92 K_bulk); λ_eff = G₀d = 162.1 W/mK; σ(δ) at the fitted
  base point: N1 6.08, N2 10.25, N2lam 24.57 ps (N2lam base), N3 64.35 ps (N3 base).

| REF, 500 nm | 5×5×3 | 7×7×5 | 9×9×6 |
|---|---|---|---|
| visible δ, N2lam / N3 (ps) | −40.60 / −37.63 | −62.70 / −61.00 | −69.35 / −67.36 |
| regularised δ (C + QDQ) (ps) | −62.46 | −153.54 | −194.20 |
| bulk δ (ps) | −139.8 | −1151.6 | −2351.7 |
| σ(δ) N2lam / N3, 0.01° (ps) | 22.40 / 58.80 | 26.51 / 66.22 | 24.57 / 64.35 |
| proxy c (N2lam / N3) | — | 2.12 / 2.17 | 2.17 / 2.22 |
| proxy half-range of the 12 witnesses, N2lam / N3 (ps) | — | 12.57 / 12.09 | 9.42 / 9.07 |
| visible half-range, N2lam / N3 (ps) | (40.6 / 39.9, film witnesses) | 11.32 / 11.08 | 8.73 / 8.58 (two end-setters, below) |

  The visible δ approaches a limit (steps +22.1 / +6.7 ps for N2lam, +23.4 / +6.4 ps for N3) while the regularised
  and bulk values keep growing: the long-mean-free-path carriers that the finer meshes add are wall-limited in the
  slab. A geometric continuation of the last two steps (ratio 0.27–0.30) gives ≈ −72 ps (N2lam) and ≈ −70 ps (N3), a
  heuristic estimate on three unequal meshes, not a converged value.
* **Proxy on 9×9×6** (δ_reg at L/c with c recalibrated on the 9×9×6 REF slab; one bordered solve per tensor): the
  same two tensors set the ends as on 7×7×5 (cal3w min −10.35 ps, orb 0.1 30 nm max +8.48 ps, N2lam), and the
  half-range contracts by ×0.75 (12.57 → 9.42 ps). With the 7×7×5 slope (0.90 / 0.91) the predicted visible
  half-range is ≈ 8.5 / 8.2 ps and the predicted margins at 0.01° are ≈ 2.9 (N2lam) and ≈ 7.8 (N3).

* **Witness slabs on 9×9×6** (same set-up; the two tensors that set the ends on 7×7×5 and, by the proxy, on 9×9×6):

| 500 nm, 0.01° | cal3w min: δ_fit; excursion from REF (N2lam / N3) | orb 0.1 30 nm max: δ_fit; excursion (N2lam / N3) | two-witness half-range (N2lam / N3) | margin N2lam / N3 |
|---|---|---|---|---|
| 7×7×5 (R-S, reproduced) | −74.99 / −72.95; −12.30 / −11.95 | −52.35 / −50.78; +10.35 / +10.22 | 11.32 / 11.08 | 2.34 / 5.97 |
| 9×9×6 (this audit) | −78.17 / −75.97; **−8.82 / −8.61** | −60.71 / −58.80; **+8.64 / +8.56** | **8.73 / 8.58** | **2.81 / 7.50** |

  Excursions shrink ×0.72 (cal3w min) and ×0.84 (orb 0.1 30 nm max) from 7×7×5 to 9×9×6; fit/proxy ratios 0.85 / 1.02
  (N2lam; 7×7×5: 0.80 / 1.06), so the proxy ordering of the ends holds and no other witness is predicted to exceed
  them (next proxy excursions −8.25 for cal1nz min and +7.71 for cal3w max): the two-witness half-range is an inner
  bound that equals the 12-witness value up to proxy error. Interior conductivity follows δ as on 7×7×5 (λ_z,fit 223.5
  / 236.2 / 241.1 W/mK for cal3w min / REF / orb 30 max). Coarse-mesh checks 2.5–3.9e-6 for both tensors; symmetrising
  the two-ports moves each δ by 0.023 ps; Δχ² ≤ 5.5e-3.
* **Margins on 9×9×6** (500 nm, σ at the 9×9×6 REF fit, two-witness half-range): **0.01°: N1 0.70, N2 1.17,
  N2lam 2.81, N3 7.50**; 0.03°: 2.08 / 3.44 / 7.85 / 21.0; 0.1°: 6.63 / 9.40 / 17.2 / 53.0. From 7×7×5 to 9×9×6 every
  margin grows (N2lam ×1.20, N3 ×1.26), N2 leaves the threshold, N1 stays below 1. Mesh sequence of the N2lam margin
  with the same extraction: 0.55 (5×5×3, film witnesses) → 2.34 (7×7×5) → 2.81 (9×9×6).

### 5.4 δ-optimised witnesses

No witness was optimised for δ. On 7×7×5 a reversal at 0.01° needs a visible half-range of 26.5 ps (N2lam, ×2.34)
or 66 ps (N3, ×5.97); at 0.1° ×14 / ×42. On 9×9×6 the factors are ×2.81 (24.6 ps) and ×7.50 (64.4 ps) at 0.01°.
Through the proxy of rs6 (δ_reg at L/c, c ≈ 2.1, slope 0.90, correlation 0.99) the 7×7×5 requirement is a
regularised half-range at L_eff ≈ 235 nm of ≈ 29 ps against 12.6 ps for the present witnesses (9×9×6: ≈ 27 ps
against 9.4 ps, L_eff ≈ 230 nm).
The present witnesses, optimised for τ_z at 30–300 nm, already move δ_reg by up to +41 % of REF at 500 nm, mostly
through the off-diagonal (normal-process) structure that the lifetime data do not constrain; nothing in the audit
excludes a factor 2.3–2.8 for a δ objective. The question is open and is the main residual risk of the N2lam verdict
(N3 needs 6–7.5, which is a different order of demand).

## 6. Points the record should state or change

1. R-S1 validation: add the independent route (§2: no parity, no Schur complement, no eigendecomposition; ≤ 1.7e-9
   on 7×7×5 with observed order 3.97, ≤ 6e-11 on 5×5×3) and the conditioning numbers (κ(Φ) ≤ 17, κ(P±) ≤ 14).
2. R-S2/R-S3: the δ values and half-ranges are exact for the 7×7×5 model (numerical error < 1e-4 ps; AAA verified on
   the range-setting witness, ≤ 1.6e-4 ps); state "UNDER-RESOLVED in the q-mesh" next to them, with the sequence at
   500 nm: REF δ (N3) −37.6 / −61.0 / −67.4 ps and half-range 11.08 → 8.58 ps on 7×7×5 → 9×9×6 (§5.3).
3. R-S4/synthesis: keep N2lam 2.34 and N3 5.97 (reproduced) and add the 9×9×6 values 2.81 / 7.50 (two end-setting
   witnesses); write N2 = 1.01 as "at the threshold" on 7×7×5 (1.17 on 9×9×6) rather than "not resolvable"; add that
   the same extraction on 5×5×3 gives 0.55 / 1.47, so the conclusion is carried by the trend 5×5×3 → 7×7×5 → 9×9×6
   (0.55 → 2.34 → 2.81, N2lam) and by inner-bound witnesses.
4. R-S6: "the full-operator 7×7×5 values are not mesh-converged; the trend in both operator families is a falling
   transfer" is confirmed by the 9×9×6 slab: the visible range contracts ×0.77; the proxy recalibrated on 9×9×6
   (c = 2.17 / 2.22) gives a 12-witness half-range of 9.42 / 9.07 ps and predicts the two end-setters' slab excursions
   to a fit/proxy ratio of 0.85–1.02, so it remains usable as the objective for δ-optimised witnesses.
5. The visible δ is band-effective (−61.0 full band, −72.8 ps in 10–200 MHz, N3); the half-range and the margins refer
   to the declared 10 kHz–200 MHz band.
6. 9×9×6: the paired reduction is not applicable (972 self-mapped moving modes on q_z = 1/2); the slab needs an
   unpaired solver; the data of that plane carry a net equilibrium current (Σ v_z c = 27.6 W m⁻²K⁻¹, |ê·b|/|b| =
   4.7e-9, R-A1) and make the slab non-reciprocal at up to 8.5e-4 (200 MHz), a gauge artefact of the degenerate pairs
   with a 0.025-ps effect on δ. A gauge-consistent (pair-averaged) operator would remove it.
7. Not re-derived independently here: the partial-specular and gray-transmission wall variants (their wall operators
   were checked against the doubling route on synthetic operators only, rs1 A3), the energy-port numbers and 300 nm on
   9×9×6. The black-wall, wall-port numbers that carry the verdict are covered.

## 7. Files

Scripts (`experiments/R_slab_audit/`), outputs (`experiments/R_slab_audit/out/`):
* `sa_core.py` — own operator assembly, bulk/regularised kinetic times, independent slab solver (orders 2 and 4).
* `sa_fdtr.py` — own 1D stack, wall/energy ports, variable-projection fits, misfit profiles, Cauchy-integral Fisher.
* `sa0_proto553.py`, `sa0b_meshstudy553.py`, `sa0c_order4_553.py` → `sa0_proto553.json`, `sa0b_meshstudy553.json`,
  `sa0c_order4_553.json`, `sa0c.log`: 5×5×3 prototype, mesh study, fourth-order convergence.
* `sa1_verify775.py` → `sa1_verify775.json`, `sa1.log`, `sa1_rta.log`: 7×7×5 operator, full-operator and RTA two-ports
  against the stored modal ones, witness spot check.
* `sa2_fits.py` → `sa2_fits.json`, `sa2.log`: REF fits, profiles, subsets, bands, sensitivity, interpolation, witness
  δ, half-ranges, Fisher margins, direct test, 5×5×3/7×7×5 mesh trend.
* `sa3_modal_diag.py` → `sa3_modal_diag.json`, `sa3.log`: parity map, 9×9×6 zone-boundary plane, conditioning,
  diffusive-mode polish.
* `sa5_aaa_witness.py` → `sa5_aaa_witness.json`, `sa5_exact80_cal3w_m3L100min.npz`, `sa5.log`: exact 80-frequency
  witness two-ports, AAA check, leave-out cross-validation.
* `sa4_m996.py` → `sa4_m996.json`, `sa4_m996_{REF,cal3w_m3L100min,orb0.1_m3L30max}_500.npz` (two-ports at the 16
  support frequencies), `sa4_op.log`, `sa4_slab500.log` + `sa4_slab500b.log` (run resumed after adding the memory
  guard; completed frequencies reused), `sa4_fit*.log`, `sa4_proxy.log`, `sa4_wit_*.log`: 9×9×6 operator, slabs,
  fits (raw and symmetrised two-ports), proxy (§5.3).

Run order: sa0, sa0b, sa0c, sa1 (op rta full wit), sa2 (ref interp wit mesh), sa3, sa5, sa4 (op; slab; fit; proxy;
wit cal3w_m3L100min 500; wit orb0.1_m3L30max 500; fit). Each step runs in ≤ 40 min on two threads (sa4 slab ≈ 35 min
per tensor at 500 nm; peak working set 1.0 GB; free-memory guard 1.5 GB).
