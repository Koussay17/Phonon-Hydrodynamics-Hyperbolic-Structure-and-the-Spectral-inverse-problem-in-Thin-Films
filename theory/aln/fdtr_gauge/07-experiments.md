# 07 — Experiments

Workstream A, branch A-C (twin numerics). Scripts in `experiments/A/`, raw outputs (JSON/CSV/logs) in
`experiments/A/out/`, figures in `figures/`. Python: research environment (NumPy 2.5.3, SciPy 1.18.1, JAX 0.11.2,
float64), at most three threads; the forward model `C:\Users\Koussay\these\src\twin_fdtr.py` is imported read-only
(repository commit e4287dd). Every number below is SYNTHETIC: no statement about a real sample (A7).

**Common set-up (A2–A7).** `au_aln_sapphire` of the twin: Au 150 W/mK, 2.49e6 J/m³K, 80 nm; R_Au/AlN = 1e-8,
R_AlN/sapphire = 2e-8 m²K/W; AlN of thickness L, C = 2.41e6 J/m³K, λ_z and τ_R from A1; sapphire 35 W/mK, 3.03e6.
AlN law: Cattaneo (τ_R, ℓ = 0; λ_r = RTA κ_x(L)) or GK (τ_ℓ = 0.5 τ_R or 2 τ_R, α = 1/3, ℓ = √(τ_ℓ λ/((1+α)C)),
GKBoundary default: slip C_s = 1, tangential 'slip', temperature 'local'; isotropic, λ_r = λ_z). Beams w_p = 10 µm,
w_r = 8 µm (note 16), or the 1D large-beam limit 2 Z(k=0)/(π(w_p² + w_r²)) of the same observable. 80 log-spaced
frequencies 10 kHz – 200 MHz (sub-bands 1–200 MHz and 10 kHz–50 MHz are subsets of these points). Noise σ(ln|H|) = 1 %,
σ_φ ∈ {0.01, 0.03, 0.1} deg, Gaussian, independent. Trapezoidal rule in x = ln(√η k), fixed step 0.04,
x_low = −22 (validated in A2). Instrument: H_obs = exp(g + iφ0 − 2πi f t_d) H.

## A1 — base values of the AlN layer at 30, 100, 300, 500 nm  [a1_base_values.py; out/a1_base_values.json]

**Regulator.** Same heuristic Matthiessen regulator in both routes: physical boundary rate |v_z|/L (diffuse grey
boundary of a c-oriented film). In the event operator it appears as d = |v_z|/(4πL) because the operator is in THz
linewidth units (1/τ = 4πΓ): `film_event.py`, `CL[di] += vz / (FOURPI * L * 10.0)`. The form d = |v_z|/(4L) written in
02-assumptions.md A1 lacks the factor π; the computations use the correct physical rate |v_z|/L.

**Routes.** (a) RTA (phono3py 4.5.0, 300 K) on the three finest meshes of fc3-memory E27, recomputed from the stored
chunks (`D:\...\20261005-fc3-memory\results_c\rta\m{484826,565630,646435}`) with the assembly formula of
`rta_chunks.py`, 500 nm added; the archived values at 30–10 000 nm are reproduced exactly (relative deviation 0).
|v_z| is invariant on the star of each irreducible point (6mm and time reversal), so the irreducible sum is exact.
(b) Conserving event operator, 12×12×8, REF tensor: the REF block of every `results_c/film3/*_m12128.json` (identical
in all twelve files, spread 0). 500 nm is not available on that route.

| L (nm) | τ_z event 12×12×8 (ps) | τ_z RTA 64×64×35 (ps) | κ_z event (W/mK) | κ_z RTA (W/mK) | κ_x RTA 64×64×35 | κ_x RTA, 1/N_z → 0 | κ_x event 12×12×8 |
|---|---|---|---|---|---|---|---|
| 30 | 5.487 | 5.303 | 62.55 | 59.56 | 112.9 | 98.9 | 148.4 |
| 100 | 11.226 | 10.747 | 115.62 | 105.96 | 160.0 | 149.3 | 192.1 |
| 300 | 19.699 | 19.511 | 167.25 | 150.31 | 200.0 | 192.3 | 228.3 |
| 500 | **25.00** (interp., see below) | **25.567** (direct) | **188.6** (interp.) | 169.34 | 215.9 | 209.6 | 242.2 (interp.) |
| 1000 | 33.936 | 36.257 | 214.36 | 192.31 | 234.0 | 229.7 | 256.5 |

* **RTA mesh convergence** (48×48×26 → 56×56×30 → 64×64×35): τ_z(L) changes by +0.07/+0.07 % (30 nm),
  +0.09/+0.09 % (100 nm), +0.16/+0.17 % (300 nm), **+0.23/+0.22 % (500 nm)**; κ_z(L) by ≤ 0.08 %.
  τ_z,RTA(500 nm) = 25.45 / 25.51 / 25.57 ps.
* **500 nm, event operator.** Three estimates: log-log interpolation 300→1000 nm 24.81 ps; event/RTA ratio
  interpolated log-log and multiplied by the direct RTA value at 500 nm 25.00 ps; PCHIP in (ln L, ln τ) through all six
  thicknesses 25.14 ps (spread 1.3 %). Interpolation error measured on the RTA curve, where 500 nm is known:
  log-log 300→1000 nm −0.74 %, PCHIP +0.05 %. Same for κ_z(500): 185.8 / 188.6 / 188.8 W/mK (log-log error on RTA
  −1.5 %). **Adopted: τ_R(500 nm) = 25.00 ps, λ_z = 188.6 W/mK** (ratio method; it interpolates only the smooth
  event/RTA ratio, 0.973 at 300 nm and 0.936 at 1000 nm). Sensitivity case: RTA direct (25.57 ps, 169.3 W/mK).
* **In-plane value λ_r: mesh artefact found.** On a uniform mesh the plane q_z = 0 has v_z = 0 exactly
  (C2z with time reversal maps (q_x, q_y, 0) to itself and v_z to −v_z; max |v_z| on the plane 2e-29 m/s), so the
  regulator never suppresses it. It carries the weight 1/N_z and a constant 944/N_z W/mK of κ_x(L) at every L
  (35.0 / 30.9 / 27.0 W/mK on the three meshes; 24 % of κ_x(30 nm) on 64×64×35). The zone-boundary plane q_z = 1/2
  does not have v_z = 0 (sticking bands of the 6_3 screw, |v_z| up to 9.1 km/s). Hence κ_x(L) converges like 1/N_z:
  linear extrapolation in 1/N_z over the three meshes fits to ≤ 0.02 W/mK and gives the column "1/N_z → 0". On
  12×12×8 the plane alone would carry ≈ 118 W/mK: the event-operator in-plane film values are strongly biased and
  are not used. τ_z and κ_z are unaffected (their weights contain v_z^2).
* **Adopted λ_r** (Cattaneo layer, anisotropic): RTA 64×64×35 κ_x(L) (112.9 / 160.0 / 200.0 / 215.9 W/mK);
  sensitivity cases λ_r = λ_z and the 1/N_z-extrapolated value. GK layers of the twin must be isotropic: λ_r = λ_z.
* **Base values used in A2–A7** (τ_R = event τ_z, λ_z = event κ_z, C = 2.41e6 J/m³K): 30 nm 5.487 ps / 62.55 W/mK;
  100 nm 11.226 ps / 115.62; 300 nm 19.699 ps / 167.25; 500 nm 25.00 ps / 188.6.
* Command: `python a1_base_values.py` (about 10 s).

## A2 — numerical verification  [a2_numerics.py, indep.py; out/a2_numerics.json, out/a2_numerics.log]

**Independent forward model (indep.py)**, written for this branch, sharing no code with the twin: thermal quadrupoles
in Hankel space by impedance recursion from the substrate upwards (Cattaneo layers with separate in-plane and
cross-plane relaxation times; 1D GK by λ_eff = λ(1 + τ_ℓ p)/(1 + τ_R p) with the local temperature convention), Gaussian
average by composite Gauss–Legendre in u = √η k on [0, 7] (56 panels × 16 nodes), derivatives by exact forward-mode
automatic differentiation (JAX, float64). The twin solves a global mode-matching system, integrates by the trapezoidal
rule in ln k and is differentiated by finite differences, so the two routes share only the physics.

| check | result |
|---|---|
| Hankel quadrature, twin fixed step 0.16 / 0.08 / 0.04 / 0.02 / 0.01 vs adaptive (epsrel 1e-12), max abs Δ ln H, 80 f | 7.4e-13 / 2.2e-14 / 2.2e-14 / 2.2e-14 / 1.6e-13 (Cattaneo and GK, 30 and 500 nm; x_low −30 instead of −22: 2.2e-14) |
| same, effect on σ(ln τ_R), all sets N0–N5, against step 0.01 | ≤ 4e-9 (500 nm), ≤ 8.4e-6 (30 nm) relative |
| FD Jacobian on the adaptive rule (h = 1e-2, nodes move between evaluations) | σ within 1.5e-4 relative of the fixed-rule value |
| FD step, τ_R column vs exact autodiff (Cattaneo 500 nm): central h = 0.1 / 0.03 / 1e-2 / 1e-3 / 1e-4 / 1e-5 | 1.7e-3 / 1.5e-4 / 1.7e-5 / 1.7e-7 / 5.5e-9 / 4.6e-8 (h² until roundoff) |
| production scheme, Richardson h = 2e-3: τ_R column error; max σ error over N0–N5 | 6.7e-10; 8.4e-9 (500 nm) — 2.8e-8; 1.6e-6 (30 nm) |
| GK (no exact reference): Richardson h = 2e-3 vs h = 1e-3 | σ within 1.8e-9 (500 nm, τ_ℓ = 2τ_R), 1.1e-6 (100 nm, τ_ℓ = 0.5τ_R), including N4 |
| 1D: twin kernel(k = 0) vs independent closed-form chain Au/R/AlN(Cattaneo or GK λ_eff)/R/sapphire, 10 kHz – 10 GHz | max abs Δ ln Z ≤ 1.9e-15 (all 12 cases) |
| 1D: σ(ln τ_R) twin FD vs independent autodiff, N0–N4 | ≤ 5.7e-6 (Cattaneo); GK ≤ 3.5e-7 for N0–N3, N4 (σ ≈ 1e5–1e7, nearly degenerate) ≤ 5.3e-3 |
| large-beam limit, twin Gaussian (w_p = 30 µm … 3 mm, w_r = 0.8 w_p) vs 1D, 500 nm | max abs Δ ln H 9.5e-2, 9.3e-3, 1.0e-3, 9.4e-5, 1.0e-5; observed order 1.93, 1.995, 1.9995, 1.9999 (Cattaneo); GK 1.95 … 2.0000 |
| finite beams, twin vs independent model, all L, λ_r = RTA or λ_z | max abs Δ ln H 2.2e-14; τ_R column ≤ 2.8e-8 relative; σ(ln τ_R) ≤ 1.6e-6 (30 nm), ≤ 8e-9 (500 nm) |
| every Cattaneo row of A4 (all sets, 3 bands, 3 noise levels) recomputed by autodiff | max relative difference 3.0e-5 (finite beams), 8.0e-4 (1D, 30 nm, N4 σ ≈ 8e3) |
| frequency grid, σ·√N for N = 20, 40, 80, 160, 320 log-spaced points | N0: 8.03, 8.26, 8.39, 8.47, 8.50 (100 nm, 0.01 deg); 1.391 … 1.373 (500 nm). N2, 0.1 deg: 324, 377, 424, 455, 473 (100 nm), 24.08 … 24.60 (500 nm). Twin = autodiff at N = 160 to 1.5e-9. |

* **Reading.** The numerical error of every σ(ln τ_R) reported in A4 is ≤ 3e-5 relative (≤ 1e-8 at 300–500 nm),
  far below every modelling sensitivity (percent level, A4 variants) and irrelevant to the verdict. The frequency-grid
  drift is not a discretisation error of an integral: the sum over points weights the 200 MHz end point fully, and at
  100 nm (N2) the τ_R information sits there, so σ√N changes as O(1/N) (+7.3 %, +3.9 % per doubling); it reports how
  the design (80 points) compares with denser designs, not a convergence failure.
* **In-plane relaxation time** (the twin uses one τ_R for q_z and q_r; physically the in-plane memory time of the film
  is τ_x(L) = 157 / 121 / 111 / 111 ps, RTA 64×64×35): with τ_r,in-plane = τ_x(L) in the independent model, ln H moves
  by ≤ 1.5e-4 (500 nm), and the cross-plane derivative alone gives σ(ln τ_R) (N0, 0.01 deg) 5.13 / 0.867 / 0.232 /
  0.146 instead of 5.25 / 0.938 / 0.246 / 0.154 (30/100/300/500 nm); N2 −11 % (500 nm) to −39 % (30 nm); N3 ≤ 0.3 %.
* Command: `python a2_numerics.py` (31 min, mostly the adaptive-rule finite differences).

## A3 — phase change Δφ(f) at the admissible relative changes of τ_R  [a3_dphi.py; out/a3_dphi.json, out/a3_dphi_curves.csv; figures/A3_dphi.png]

Set-up as stated at the top of this file. Δφ(f) = arg H(τ_R(1+w)) − arg H(τ_R), everything else fixed, exact (no
linearisation).

**Maximum |Δφ| (deg) over 10 kHz – 200 MHz, Cattaneo, beams 10/8 µm** (1D limit 4–12 % larger; GK variants within
+8 % of Cattaneo; max |Δφ(+w)| − max |Δφ(−w)| ≤ 6e-8 deg (≤ 1e-4 relative), i.e. Δφ is odd in w and linear in τ_R;
max |Δ ln|H|| ≤ 5.5e-5 at w = 8 %):

| L | τ_R (ps) | w = 2.8 % | 3.7 % | 5.2 % | 8 % | f at max | max abs dφ/d ln τ_R (deg) |
|---|---|---|---|---|---|---|---|
| 30 nm | 5.49 | **1.55e-5** | 2.05e-5 | 2.89e-5 | 4.44e-5 | 200 MHz | 5.6e-4 |
| 100 nm | 11.23 | 9.5e-5 | **1.26e-4** | 1.77e-4 | **2.72e-4** | 156 MHz | 3.4e-3 |
| 300 nm | 19.70 | 3.43e-4 | 4.53e-4 | **6.37e-4** | 9.80e-4 | 83 MHz | 1.23e-2 |
| 500 nm | 25.00 | 5.61e-4 | 7.41e-4 | 1.04e-3 | 1.60e-3 | 65 MHz | 2.00e-2 |
| 500 nm, RTA base | 25.57 | 6.16e-4 | 8.14e-4 | 1.14e-3 | 1.76e-3 | 65 MHz | — |

Bold: the half-width admissible at that thickness (00-question.md). The largest phase change at an admissible w is
**1.6e-3 deg** (500 nm, w = 8 %), six times below the smallest per-point phase noise considered (0.01 deg) and 60
times below 0.1 deg; at 30 nm it is 1.6e-5 deg. Sub-bands: the maximum lies inside 1–200 MHz for every L; inside
10 kHz–50 MHz it is 0.65 (30 nm) to 0.89 (500 nm) of the full-band value. Command: `python a3_dphi.py` (77 s).

## A4 — local Fisher analysis, nuisance sets N0–N5  [a4_fisher.py; out/a4_fisher.json, out/a4_wstar_table.csv; figures/A4_wstar.png]

**Method.** Whitened Jacobian of [ln|H_obs|/σ_A ; arg H_obs/σ_φ], σ_A = 1 %, σ_φ ∈ {0.01, 0.03, 0.1} deg, independent
Gaussian. Positive parameters in logarithms; instrument columns exact (g → 1, φ0 → i, t_d → −2πif). Thermal columns:
Richardson-extrapolated central differences (h = 2e-3) with the quadrature rule frozen. σ(ln τ_R) from the SVD of
the column-normalised Jacobian (cond > 1e10 flagged, not inverted). Sets as declared: N0 {τ_R}; N1 + g, φ0, t_d;
N2 + R_Au/AlN, R_AlN/sap; N3 + λ_z, C (AlN); N4 + τ_ℓ (GK only) + λ_Au, C_Au; N5 + a common factor on both radii.
H depends on (w_p, w_r) only through η = (w_p² + w_r²)/8, so freeing both radii is exactly degenerate (flagged,
cond 3–4e12); N5 frees their common factor. In the 1D limit the beam factor is a pure gain and N5 is flagged.
w*_1σ = σ(ln τ_R), w*_2σ = 2σ(ln τ_R). A5 shows that H is linear in τ_R itself over 0 ≤ τ_R ≤ 4τ_R,0, so these are
also the resolvable relative changes w = δτ_R/τ_R, in both directions and not only locally.

**w*_1σ = σ(ln τ_R), Cattaneo, λ_r = RTA κ_x(L), beams 10/8 µm, 80 frequencies 10 kHz – 200 MHz:**

| L | σ_φ (deg) | N0 | N1 | N2 | N3 | N4 | N5 | admissible w(L) |
|---|---|---|---|---|---|---|---|---|
| 30 nm | 0.01 | 5.25 | 12.4 | 188 | 526 | 6.7e3 | 6.8e3 | 0.028 |
| | 0.03 | 15.7 | 36.9 | 532 | 1.1e3 | 1.3e4 | 1.3e4 | |
| | 0.1 | 51.4 | 114 | 1.3e3 | 1.7e3 | 3.0e4 | 3.6e4 | |
| 100 nm | 0.01 | 0.938 | 3.64 | 5.68 | 11.6 | 734 | 786 | 0.037 (0.08 calibrated classes) |
| | 0.03 | 2.81 | 10.8 | 16.7 | 33.2 | 2.0e3 | 2.0e3 | |
| | 0.1 | 9.24 | 31.4 | 47.4 | 82.4 | 4.4e3 | 4.4e3 | |
| 300 nm | 0.01 | 0.246 | 0.470 | 0.662 | 2.63 | 84.4 | 87.1 | 0.052 |
| | 0.03 | 0.737 | 1.40 | 1.94 | 7.53 | 194 | 195 | |
| | 0.1 | 2.40 | 4.35 | 5.36 | 19.3 | 286 | 289 | |
| 500 nm | 0.01 | 0.154 | 0.205 | 0.340 | 2.25 | 8.98 | 9.00 | not computed; trend 0.055–0.08 |
| | 0.03 | 0.460 | 0.612 | 0.997 | 6.46 | 24.3 | 24.7 | |
| | 0.1 | 1.48 | 1.94 | 2.73 | 17.1 | 49.6 | 56.7 | |

* **Ratio w*_1σ / w(L) at the most favourable entry (N0, 0.01 deg): 188 (30 nm), 25 (100 nm), 4.7 (300 nm),
  1.9–2.8 (500 nm against the trend bracket).** At the falsification test of 00-question.md (N ⊇ N2,
  σ_φ = 0.1 deg): w*_1σ = 1.3e3 / 47 / 5.4 / 2.7 at 30 / 100 / 300 / 500 nm, i.e. 4.5e4 / 1.3e3 (590 against 8 %) /
  100 / 34–50 times the admissible half-width. No entry of any table is below w(L), at 1σ or 2σ.
* **Phase noise needed** for 2σ(ln τ_R) = ln(1 + w(L)), σ_A = 1 %, full band: N0 2.6e-5 / 1.9e-4 / 1.0e-3 /
  1.65e-3 deg (30/100/300/500 nm, w = 2.8/3.7/5.2/5.2 %); N2 7.3e-7 / 3.2e-5 / 3.8e-4 / 7.4e-4 deg; N3 2.5e-7 /
  1.6e-5 / 9.6e-5 / 1.1e-4 deg; with w = 8 % at 100 nm, N2 6.8e-5 deg. The required noise is 6 (500 nm, N0) to
  1e4 (30 nm, N2) times below 0.01 deg.
* **Bands** (subsets of the 80 points): the τ_R information sits above 1 MHz. 1–200 MHz reproduces N0 to 0.1 % and
  costs a factor 1.1–1.6 in N1/N2 at 300–500 nm; 10 kHz–50 MHz costs 1.6–2.3 in N0 and a factor 10–17 in N1 at
  300–500 nm (below 50 MHz the τ_R direction is nearly i ω × const, collinear with the delay). 500 nm, 0.01 deg,
  1–200 MHz: N0 0.154, N1 0.233, N2 0.532, N3 8.09; 10 kHz–50 MHz: 0.251, 3.56, 3.61, 11.8.
* **Degeneracy partners** (correlation of ln τ_R with the most correlated nuisance, full band, 0.1 deg): N1 the
  delay (+0.93 at 100 nm, +0.41 at 500 nm; dropping t_d from N1 gives 1.16 instead of 3.64 at 100 nm, 0.186 instead
  of 0.205 at 500 nm, 0.01 deg); N2 the delay at 100 nm (+0.87), φ0 at 500 nm (−0.73); N3 λ_z (−0.98 at 500 nm);
  N4 C_Au (+0.9997 at 100 nm, +0.90 at 500 nm). Variance inflation over N0 at 500 nm, 0.1 deg: 1.3 / 1.8 / 11.5 /
  33.5 / 38.3 (N1 … N5). Condition numbers of the normalised Jacobians (100 and 500 nm, 0.1 deg): 2–13 (N1, N2),
  150–210 (N3), 210–3.6e3 (N4, N5); none of the declared finite-beam sets is numerically unidentifiable; the τ_R
  direction is simply weak.
* **Variants** (relative change of σ(ln τ_R) against the primary row, N0–N3, both noise levels, all four L;
  complete table in out/summary_tables.md): λ_r = λ_z: N0 ≤ 2 %, N1 −2…+6 %, N2 −19 % (30 nm), −8 % (100 nm),
  ≤ 3 % (300–500 nm), N3 ≤ 2.5 %; λ_r = 1/N_z-extrapolated κ_x: ≤ 0.5 % except N2 at 30 nm (−5 %). 1D limit: N0 −2…−8 %,
  N1 −13…+13 %, N2 −42 % (30 nm), −27 % (100 nm), −10…−11 % (300–500 nm), N3 +5…−18 %. GK with τ_ℓ = 0.5 or 2 τ_R held
  fixed, compared with the isotropic Cattaneo row: N0 ±1.4 %, N1 ±7 %, N2 −2…−19 %, N3 −22…−53 % (the isotropic GK
  layer ties λ_r to λ_z, so lateral spreading informs λ_z); in 1D GK and Cattaneo agree to ±3 % (N0–N2), −3…+12.5 %
  (N3). 500 nm with the RTA base: −4…−23 % (N0 0.140 / 1.35, N2 0.297 / 2.41). Extra sets at 500 nm, 0.01 deg:
  N3 + λ_r 5.10; N5 + both thicknesses 10.1.
* **GK, τ_ℓ free (N4).** In 1D, τ_R and τ_ℓ enter only through λ_eff = λ(1 + τ_ℓ p)/(1 + τ_R p), degenerate at first
  order in ωτ: σ(ln τ_R) = 2.1e3 (500 nm) to 1.9e7 (30 nm) in 1D. With finite beams the twin gives 5.2 (500 nm) to
  8.5e3 (30 nm), smaller, because at k > 0 the
  slip boundary layers of the isotropic GK film make its in-plane conductance depend on ℓ (Poiseuille-type
  suppression): moving τ_R and τ_ℓ together at fixed τ_R − τ_ℓ changes ln Z(k) by up to 2e-3 per unit relative
  change at k = 1e5 m⁻¹, 100 kHz (2e4 times the pure-τ_R change), while at k = 0 the ratio is ωτ_ℓ exactly
  (7.05e-3 at 200 MHz, 100 nm). That information is a property of the isotropic GK twin with slip walls, not of the
  τ_R signature, and every N4 value is ≥ 5 anyway.
* **Independent check of every Cattaneo row**: exact forward-mode autodiff of the independent quadrupole model
  (indep.py) against the twin's finite differences; see A2 for the maximum differences.
* Command: `python a4_fisher.py` (740 s, single thread); the CSV is rewritten from the JSON by `python a4_csv.py`
  with w*_1σ = σ and w*_2σ = 2σ as relative changes of τ_R (the first run wrote exp(σ) − 1, a log-space reading that
  A5 shows to be inappropriate); `python summarize_a.py > out/summary_tables.md` prints every table.

## A5 — profile likelihood: is the Fisher number misleading?  [a5_profile.py; out/a5_profile.json; figures/A5_profile.png]

500 nm, Cattaneo, beams 10/8 µm, 80 frequencies. (a) Noiseless profile with the twin: data = twin at τ_R,0 = 25.0 ps;
τ_R fixed at τ_R,0(1 + w), nuisances of the set refitted (instrument by exact linear projection; R's, λ_z, C by
Levenberg–Marquardt, tolerances 1e-15). (b) Noisy data from the twin's `synthetic_measurement` (σ_A = 1 %,
σ_φ = 0.01 or 0.1 deg, seeds 20261111 + i, i < 200), joint maximum-likelihood fit of ln τ_R and the nuisances with the
independent model (exact Jacobian; ln τ_R bounded to [−12, 4]); three seeds refitted with the twin.

* **Noiseless profile = Fisher prediction in the relative change w**, for N0–N3 and both noise levels:
  Δχ²(±2.8 %) = 0.03306 / 0.03306 (N0, 0.01 deg; (0.028/0.1540)² = 0.03306), 0.006767 (N2, 0.01 deg; 0.006766),
  1.54e-4 / 1.55e-4 (N3, 0.01 deg; 1.54e-4), 1.05e-4 (N2, 0.1 deg; 1.05e-4); Δχ²(±8 %) = 0.05523 / 0.05525 (N2, 0.01 deg;
  0.05524). Symmetric in w, not in ln(1 + w): the response is linear in τ_R. Across w ∈ [−1, 3] the profile follows
  (w/σ_F)² to ≤ 2.4 % for N0–N2 (worst N1 at w = 3; Fourier film 42.19 against 42.17 for N0, 8.62 against 8.63 for N2,
  0.01 deg); N3 rises above it at large w (+6 % at w = 1, +18 % at w = 3), i.e. the Fisher number is conservative
  there. At the admissible changes (|w| ≤ 8 %) the deviation is ≤ 0.06 % (N0–N2) and ≤ 0.44 % (N3).
* **Noisy data, 200 seeds** (estimates of the relative change τ̂_R/τ_R,0 − 1; robust width = half the 16–84 % range):

| set, σ_φ | σ_F | robust width | fraction at τ_R → 0 (Gaussian-in-τ_R prediction) | ln τ̂_R coverage of ±σ_F | χ²/dof |
|---|---|---|---|---|---|
| N0, 0.01° | 0.154 | 0.148 | 0 (0.000) | 0.695 | 0.988 |
| N1, 0.01° | 0.205 | 0.174 | 0 (0.000) | 0.730 | 0.987 |
| N2, 0.01° | 0.340 | 0.331 | 0.005 (0.002) | 0.700 | 0.987 |
| N0, 0.1° | 1.48 | 1.26 (clipped) | 0.23 (0.25) | 0.720 | 0.989 |
| N1, 0.1° | 1.94 | 1.30 (clipped) | 0.27 (0.30) | 0.705 | 0.988 |
| N2, 0.1° | 2.73 | 1.93 (clipped) | 0.32 (0.36) | 0.670 | 0.989 |

  The estimator behaves as a Gaussian in τ_R with width σ_F, truncated at τ_R = 0 (at 0.1 deg a quarter to a third
  of the fits run to the Fourier boundary, as predicted). In ln τ_R the scatter has heavy tails (std 0.92 for N2 at
  0.01 deg; 5–6 at 0.1 deg) because of that truncation, so σ(ln τ_R) must be read as a relative change of τ_R. The
  twin refits of seeds 20261111–3 (N2, 0.01 deg) give ln τ̂_R = +0.105544, −0.186128, +0.477517, equal to the
  independent-model fits to 1e-7.
* **Reading.** The Fisher numbers of A4 are not misleading: they are exact for the admissible changes (linear
  regime) and describe the scatter of the maximum-likelihood estimator. Near or above σ_F ≈ 1 the data do not bound
  τ_R from below at all (the Fourier film is inside 1σ).
* Command: `python a5_profile.py` (36 min).

## A6 — large-beam Fisher limit, phase/amplitude split, size of τ_R, and a model-internal design sweep  [a6_design.py; out/a6_design.json]

* **Large-beam limit of the Fisher numbers** (Cattaneo, 500 nm, w_r = 0.8 w_p; σ/σ_1D − 1 at 0.1 deg): w_p = 10 / 30 /
  100 / 300 / 1000 / 3000 µm: N0 +5.5e-2, +6.8e-3, +6.2e-4, +6.9e-5, +6e-6, +1e-6 (≈ 1/w²); N1–N2 likewise; N3 +2.2e-1 …
  +1.1e-5; N4 −0.32, −0.38, −0.39, −0.094, +0.012, +0.0015. N4 converges only beyond 1 mm: in 1D the 80 nm Au film is
  thermally thin and λ_Au is nearly unidentifiable, so the small lateral-spreading information on λ_Au dominates it
  until the beams are very large. N0–N3 reach the 1D values as expected.
* **Phase / amplitude split** (σ_φ = 0.1 deg, σ_A = 1 %; gain removed from phase-only sets, φ0 and t_d from
  amplitude-only sets): 500 nm N0 both 1.48, phase only 1.54, amplitude only 5.32; N2 2.73 / 3.41 / 6.85; N3 17.1 /
  22.7 / 77.1. 300 nm N2 5.37 / 6.64 / 13.4. The phase carries most of the information; the 1 % amplitude adds 5–25 %.
* **Size of τ_R** (the regulator is heuristic; 500 nm, σ(ln τ_R) with τ_R multiplied by 0.5 / 1 / 2 / 5 / 10 / 35):
  N0 (0.01 deg) 0.308 / 0.154 / 0.077 / 0.0308 / 0.0154 / 0.0036; N2 (0.1 deg) 5.48 / 2.73 / 1.36 / 0.535 / 0.260 /
  0.059. For N0–N2, σ ∝ 1/τ_R within 5 % up to ×10 (the response is linear in τ_R while ωτ_R ≪ 1); N3 falls
  faster (17.1 → 0.97 at ×10). For the falsification test (N2, 0.1 deg) to be met at 500 nm against w = 5.5–8 %,
  τ_R(500 nm) would have to be about 30–35 times larger
  (≈ 0.9 ns, beyond the bulk value of the same operator, 70 ps on 12×12×8); at 300 nm even ×35 gives N2 (0.1 deg) 0.14,
  still 2.7 times w(300 nm).
* **Design sweep** (independent model, exact Jacobians; twin checked at three points to ≤ 1e-9; constant density
  18.6 points per decade from 10 kHz; MODEL-INTERNAL: surface absorption, Fourier Au without electron–phonon coupling
  or optical penetration — not valid in the GHz range, so this only indicates scales). σ(ln τ_R), N2, 0.01 / 0.1 deg:

| L | Au | f_max 0.2 GHz | 1 GHz | 5 GHz | 20 GHz |
|---|---|---|---|---|---|
| 300 nm | 80 nm | 0.658 / 5.33 | 0.467 / 4.02 | 0.174 / 1.70 | 0.157 / 1.53 |
| 300 nm | 20 nm | 0.741 / 6.08 | 0.121 / 1.01 | 0.067 / 0.644 | 0.054 / 0.521 |
| 500 nm | 80 nm | 0.338 / 2.72 | 0.178 / 1.68 | 0.112 / 1.09 | 0.105 / 1.02 |
| 500 nm | 20 nm | 0.214 / 1.94 | 0.070 / 0.635 | 0.049 / 0.477 | 0.042 / 0.407 |

  The information saturates once the Au film is thermally thick (above ≈ 5 GHz for 80 nm). Even a 20 GHz band with a
  20 nm transducer at 0.01 deg reaches only w*_1σ ≈ 4–5 % in N2 (comparable to w(L), not resolvable at 2σ); at
  0.1 deg nothing comes within a factor 5. N0 with 20 nm Au at 0.2 GHz: 0.086 (300 nm), 0.055 (500 nm), 0.01 deg.
* Command: `python a6_design.py` (6 min).

## A7 — is the memory time itself visible? Fourier film and large changes  [a7_fourier_limit.py; out/a7_fourier_limit.json]

Noiseless profile: data = twin at the Cattaneo base point; alternative τ_R = 0 (Fourier film), 2τ_R, 10τ_R; nuisances
of the set refitted (instrument by exact linear projection, thermal ones by Levenberg–Marquardt). Δχ²_min (in units
where 1 = 1σ, 4 = 2σ), with the linear prediction (w/σ_F)² in brackets for the Fourier film (w = −1):

| L | N0 0.01° | N0 0.1° | N1 0.01° | N2 0.01° | N2 0.1° | N3 0.01° |
|---|---|---|---|---|---|---|
| 30 nm | 0.036 (0.036) | 3.8e-4 | 6.5e-3 | 3e-5 | 0 | 0 |
| 100 nm | 1.14 (1.14) | 0.012 | 0.076 | 0.031 | 4.4e-4 | 7.0e-3 (7.4e-3) |
| 300 nm | 16.5 (16.5) | 0.174 | 4.53 | 2.28 (2.28) | 0.035 | 0.133 (0.145) |
| 500 nm | 42.2 (42.2) | 0.457 | 23.6 | 8.62 (8.63) | 0.133 | 0.186 (0.197) |

* The response is linear in τ_R itself: Δχ²(τ_R = 0) = Δχ²(2τ_R) = (1/σ_F)² to ≤ 0.8 % (N0–N2, all L and noise
  levels), and Δχ²(10τ_R) = (9/σ_F)² to ≤ 8 % (N0–N2; worst N1 at 500 nm); N3 departs from linearity at large changes
  (10τ_R: 27.6 against 11.7 at 300 nm), not at the admissible ones (A5).
* **Reading.** Even the existence of a film memory time (Cattaneo against Fourier) is resolved at 2σ only for
  L ≥ 300 nm, N ≤ N1 (300 nm) or N ≤ N2 (500 nm), and σ_φ = 0.01 deg. At σ_φ = 0.1 deg it is never resolved
  (Δχ² ≤ 0.46); at 30–100 nm no set resolves it at 2σ (Δχ² ≤ 1.14, the only value above 1 being N0 at 0.01 deg,
  100 nm). A 3–8 % change of τ_R is a further factor 12–36 smaller in w (Δχ² smaller by w²).
* Command: `python a7_fourier_limit.py` (17 min).

## Branch A-C — synthesis (for 05-evidence-ledger and 13-final-report)

* **Result (NUMERICALLY DEMONSTRATED within the twin model and the stated noise model).** The static-data gap of the
  film memory time is not resolvable by FDTR in the Au(80 nm)/AlN/sapphire scenario: for every thickness (30, 100,
  300, 500 nm), every nuisance set N0–N5, every phase noise 0.01–0.1 deg and every band, w*_1σ = σ(ln τ_R) exceeds the
  admissible half-width w(L), and w*_2σ = 2σ(ln τ_R) exceeds it by twice as much. Closest approach: N0 (no nuisance),
  0.01 deg, full band: w*_1σ / w(L) = 188 (30 nm), 25 (100 nm), 4.7 (300 nm), 1.9–2.8 (500 nm against the trend
  bracket 5.5–8 %). Falsification test of H_A (N ⊇ N2, σ_φ = 0.1 deg): w*_1σ / w(L) = 4.5e4, 1.3e3, 100, 34–50.
  H_A survives with a margin of more than an order of magnitude at every thickness.
* **Mechanism.** In the film τ_R appears as a reactive correction ≈ iωτ_R L/λ_z to the film resistance (an
  inductance τ_R L/λ_z): ωτ_R ≤ 0.03 at 200 MHz, the film resistance is 1.6–8 % of the resistance below the
  transducer, and the transducer heat capacity shunts the highest frequencies. Maximum |dφ/d ln τ_R| = 5.6e-4,
  3.4e-3, 1.2e-2, 2.0e-2 deg (30/100/300/500 nm); the largest admissible Δφ is 1.6e-3 deg. The response is linear
  in τ_R, and the τ_R direction is collinear, at O(1), with the other reactive terms the nuisance sets free:
  instrument delay (N1), film heat capacity and λ_z (N3), transducer heat capacity (N4); in GK, τ_ℓ.
* **Even the existence of a memory time** (Fourier film against Cattaneo) is resolved at 2σ only for L ≥ 300 nm
  with N ≤ N1/N2 at 0.01 deg, and never at 0.1 deg.
* **What would change the verdict** (σ scales as σ_φ, as 1/√(repetitions) for white noise, and as 1/τ_R):
  phase noise 3.8e-4 deg (300 nm, N2) or 1.0e-3 deg (300 nm, N0) instead of 0.01–0.1 deg; or a film memory time
  ≈ 30–35 times larger (500 nm, N2, 0.1 deg). A GHz band with a thin transducer (model-internal, A6) still leaves N2
  at ≈ 4–5 % (0.01 deg) and ≥ 0.4 (0.1 deg). **Averaging caveat:** the verdict treats σ_φ as the final per-frequency
  noise of the measurement (as the question states). If σ_φ were the noise of one 80-point sweep and n independent
  sweeps with white noise were averaged, 2σ resolution of w(L) would need n ≈ 150 / 70 (500 nm, N2, 0.01 deg, against
  5.5 / 8 %), 650 (300 nm, N2, 0.01 deg), 9.4e4 (100 nm, N2, 0.01 deg), 4.7e3–9.9e3 (500 nm, N2, 0.1 deg) and
  4.3e4 (300 nm, N2, 0.1 deg); N0 at 500 nm, 0.01 deg: 15–31. Only the 500 nm, 0.01 deg cases are within reach of
  averaging, and there the resolution then rests on phase stability and model fidelity at the 1e-3 deg level.
* **Numerical status.** Numerical errors of σ(ln τ_R) ≤ 3e-5 relative (two independent implementations, exact
  autodiff against Richardson differences, quadrature to 1e-13, 1D closed form to 2e-15, large-beam order 2);
  profile likelihood and 200-seed maximum-likelihood fits confirm the Fisher numbers (A5).
* **Conditional on**: the twin physics (surface absorption, Fourier Au without electron–phonon coupling, Cattaneo or
  GK film with the stated boundary conditions, scalar τ_R for both flux components), white independent noise, the
  identification τ_R = τ_mem,z(L) with the heuristic |v_z|/L regulator, and the admissible widths as inner bounds
  (02-assumptions A1–A6). Freeing further parameters can only increase the marginal σ(ln τ_R) (the marginal Fisher
  variance is monotone in the nuisance set); an omitted effect that changes the base response could move σ either
  way, but it would have to raise the τ_R sensitivity by the margins above. The two inputs that could flip the
  verdict are a much larger true τ_R (A6: ≈ 30–35 times at 500 nm for N2, 0.1 deg) or much wider admissible ranges
  (E39). Because the response is linear in τ_R, the N2, 0.1 deg verdict flips only if the admissible range reached
  τ_R ≥ 3.7 τ_R,0 (500 nm) or 6.4 τ_R,0 (300 nm), against +5 to +8 % today.

## R-A — revision branch for objections A-1, A-2, A-m1, A-m3 (12 November)

Scripts in `experiments/R_A/`, raw outputs (JSON, NPZ, logs) in `experiments/R_A/out/`. Research Python, two threads,
one heavy process at a time, peak working set ≤ 1.3 GB (the E39 multistart, 15 GB, ran on the same machine); nothing
was run on the 12×12×8 operator. Every number is SYNTHETIC and conditional on the twin physics, the noise model and
the regularised operator (02-assumptions A1–A7). Run order: `run_chain1.sh`, `run_chain2.sh` (after
`ra0_kinetic_tests.py`, `ra1_film_event_ext.py m996`, `ra1b_dense_check.py`).

## R-A1 — film metrics of the twelve three-mesh optima at 500 nm and at the regulator lengths L/c  [R_A/ra1_film_event_ext.py; out/ra1_film_m775.json, out/ra1_film_m996.json]

* **Method.** Copy of fc3-memory `scripts_c/film_event.py` (E30, low-memory assembly of E36; FC3Class imported
  read-only from the campaign scripts), applied to REF and the twelve E33 optima (orb 0.1: 30/100/300 nm max/min;
  cal1nz, cal3wnz, cal3w: 100 nm max/min) with their stored E37 amplitudes (`film2_pp`) on 7×7×5 and 9×9×6. Changes:
  (i) 26 lengths from 15 to 1000 nm: 500 nm, and L/c for c ∈ {1/2, 2/3, 3/4, 4/3, 3/2, 2} at L = 30, 100, 300,
  500 nm; (ii) z column only; (iii) the regulator enters as Q D Q + γ ê êᵀ (Q = 1 − ê êᵀ, γ = tr C/n), its
  energy-conserving form, which equals D on the odd sector; (iv) the A-2 quantities of R-A6.
* **Reproduction of E37** (τ_z, κ_z, τRTA_z at 30/100/300/1000 nm, all 13 tensors): 6.7e-16 on 7×7×5. On 9×9×6:
  2.5e-5 (τ_z), 3.1e-6 (κ_z), 0 (τRTA_z). The 9×9×6 difference is not numerical noise of either route: on 9×9×6 the
  current vector has an energy component |ê·b|/|b| = 4.7e-9 (4e-18 on 7×7×5), carried entirely by the 972 modes of
  the plane q_z = 1/2 (N_z = 6), where the 6₃ screw makes bands stick in degenerate pairs whose velocities are
  assigned arbitrarily; C + D (E37) and the conserving form treat this component differently. It is irrelevant for
  every ratio below (≤ 2.5e-5 against excursions ≥ 1e-3).
* **Identity (A-m3).** The rate c|v_z|/L depends on L and c only through L/c. Hence, for every operator,
  τ_mem(L; c) = τ_mem(L/c; 1), κ(L; c) = κ(L/c; 1), and the static-data half-width w(L; c) = w(L/c; 1): a different
  normalisation of the heuristic regulator is a different effective thickness. Checked on the RTA route to 5.6e-9
  (rounding of L/c).
* Cost: 7×7×5 setup 108 s, 414 s in total; 9×9×6 setup 446 s, 2036 s in total, peak 1.28 GB.

## R-A2 — RTA film route for the regulator normalisation c  [R_A/ra2_rta_film_c.py; out/ra2_rta_film_c.json]

RTA 64×64×35 (E27 chunks), 1/τ_L = 1/τ_3ph + c|v_z|/L, A1 formulas (A1 RTA values reproduced exactly; an
independent capacity-weight evaluation of τ_mem agrees with the mode-kappa route to 1.0e-7). Base values used in
R-A4: τ_R(L; c) = τ_R,event(L) · τ_RTA(L/c)/τ_RTA(L), λ_z(L; c) = λ_z,event(L) · κ_z,RTA(L/c)/κ_z,RTA(L),
λ_r(L; c) = κ_x,RTA(L/c).

| τ_R(L; c) (ps) | c = 1/2 | 2/3 | 3/4 | 1 | 4/3 | 3/2 | 2 |
|---|---|---|---|---|---|---|---|
| 30 nm | 8.33 | 7.04 | 6.55 | 5.49 | 4.56 | 4.22 | 3.47 |
| 100 nm | 16.39 | 14.02 | 13.15 | 11.23 | 9.56 | 8.94 | 7.56 |
| 300 nm | 28.37 | 24.43 | 22.96 | 19.70 | 16.88 | 15.84 | 13.55 |
| 500 nm | 35.45 | 30.77 | 28.99 | 25.00 | 21.49 | 20.18 | 17.30 |

(RTA direct at 500 nm: 36.26 / 25.57 / 17.69 ps for c = 1/2, 1, 2, as R6 of the review.) Check of the c-scaling
τ_z(L/c)/τ_z(L) against the event operator (REF), 500 nm: c = 1/2: 1.418 (RTA 64×64×35), 1.276 (7×7×5), 1.323 (9×9×6);
c = 2: 0.692, 0.743, 0.729; at 30–100 nm the routes agree to ≤ 3.2 %. The coarse event meshes underestimate the
c-dependence at long effective lengths (their 1000 nm values are not converged, A1), so the box uses the RTA
64×64×35 scaling.

## R-A3 — static-data half-width at 500 nm and at the lengths L/c (INNER BOUNDS, no new optimisation)  [R_A/ra3_widths.py; out/ra3_widths.json]

Ratios ρ = τ_z,W(L)/τ_z,REF(L) of the twelve E33 optima from R-A1. Robust excursion per witness: the one of smaller
magnitude on 7×7×5 and 9×9×6 if both have the same sign, otherwise 0 (5×5×3 film values are not available at these
lengths; 5×5×3 was never binding at the optimised lengths). Ends of the range: largest robust excursion of each sign
over the witnesses of a class; half-width w = mean of the two ends (campaign convention, 300 nm [0.9491, 1.0537] →
5.2 %). 12×12×8 (held out; stored at 30/100/300/1000 nm only): kept-fraction interpolation (κ = (ρ₁₂ − 1)/(ρ₉ − 1) at
the stored lengths, interpolated in ln L, applied to the direct 9×9×6 value), with the shape-of-9×9×6 and the plain
log-linear interpolations as checks. Validation where the value is known: predicting 9×9×6 from 7×7×5 with the
kept-fraction rule errs by −0.0002 and −0.0040 for the two binding witnesses at 500 nm (orb 0.1 300 nm max/min),
up to 0.013 for witnesses whose excursion is small or not monotone in L; plain log-linear interpolation of the
9×9×6 values themselves errs by up to 0.0067 at 500 nm, the shape rule by up to 0.026 (it is ill-conditioned when the
excursion changes slope between 300 and 1000 nm). A negative result worth keeping: neither ρ nor its 9×9×6 shape can
be interpolated between 300 and 1000 nm to better than about half a point.

**500 nm.** τ_z(500 nm) ratio to REF:

| witness | 7×7×5 | 9×9×6 | robust | 12×12×8 kept / shape / log-lin |
|---|---|---|---|---|
| orb 0.1, 30 nm max | 0.9941 | 0.9833 | −0.59 % | 0.9875 / 0.9885 / 0.9875 |
| orb 0.1, 30 nm min | 0.9960 | 1.0110 | 0 (signs differ) | 1.0120 / 1.0018 / 1.0037 |
| orb 0.1, 100 nm max | 1.0269 | 1.0267 | +2.67 % | 1.0350 / 1.0254 / 1.0233 |
| orb 0.1, 100 nm min | 0.9650 | 0.9730 | −2.70 % | 0.9703 / 0.9682 / 0.9703 |
| **orb 0.1, 300 nm max** | 1.0670 | 1.0628 | **+6.28 %** | 1.0490 / 1.0488 / 1.0462 |
| **orb 0.1, 300 nm min** | 0.9265 | 0.9402 | **−5.98 %** | 0.9422 / 0.9416 / 0.9459 |
| cal1nz max / min | 1.0453 / 0.9562 | 1.0512 / 0.9782 | +4.53 / −2.18 % | 1.0414 / 0.9601 (kept) |
| cal3wnz max / min | 1.0419 / 0.9589 | 1.0522 / 0.9711 | +4.19 / −2.89 % | 1.0332 / 0.9587 (kept) |
| cal3w max / min | 1.0643 / 0.9647 | 1.0604 / 0.9695 | +6.04 / −3.05 % | 1.0502 / 0.9470 (kept) |

* **w(500 nm) ≥ 6.1 %** on the data meshes (range [−5.98 %, +6.28 %], set by the orb 0.1 300 nm optima; the same for
  the union of all classes; the widest calibrated witness reaches +6.04 % / −3.05 %). Held out on 12×12×8:
  **5.0–5.4 %** (kept 5.34 %, shape 5.36 %, log-linear 5.02 %; upper end +4.6 … +4.9 %, lower −5.4 … −5.8 %;
  interpolation uncertainty about ±0.4 point). These are inner bounds: no witness was optimised at 500 nm.
* The power law through the optimised orb 0.1 half-widths (2.77 / 3.74 / 5.23 % at 30/100/300 nm; slope 0.28–0.31)
  gives 5.95–6.12 % at 500 nm, which the non-optimised witnesses already reach: the trend is not an upper estimate,
  and an optimisation at 500 nm would widen the range. Against the earlier trend bracket 5.5–8 % (A1/A4), the computed
  inner bound sits in its lower half.
* Reproduction at the optimised lengths: the R-A1 ratios equal the E33/E37 phono3py-route values (e.g. 300 nm max
  1.0520 on 9×9×6 against optimiser 1.0537, E37 1.0520).

**Half-width as a function of the effective length L/c** (inner bounds; robust 7×7×5/9×9×6; held-out 12×12×8 kept):

| L/c (nm) | 15 | 30 | 60 | 100 | 150 | 200 | 300 | 500 | 750 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|
| orb 0.1 | 2.28 | 2.72 | 3.19 | 3.64 | 3.83 | 4.01 | 5.04 | 6.13 | 6.63 | 6.51 % |
| all classes | 3.38 | 5.00 | 6.65 | 7.56 | 7.82 | 7.63 | 6.75 | 6.13 | 6.63 | 6.51 % |
| all classes, 12×12×8 | 2.50 | 3.99 | 5.44 | 6.27 | 6.84 | 6.90 | 6.33 | 5.40 | 5.80 | 5.90 % |

The calibrated classes (optimised at 100 nm) dominate the union between 30 and 300 nm, the orb 0.1 300 nm optima above.
With c ∈ [1/2, 2] the half-width relevant at L = 500 nm is w(250 … 1000 nm) = 4.6–6.6 % (orb 0.1) and 6.1–7.2 % (all).

## R-A4 — sensitivity box of the A verdict (objection A-1 (ii)–(iii), A-m3)  [R_A/ra4_box.py; out/ra4_box.json, out/ra4_box.npz, out/ra4_au_subbox.npy]

* **Box (declared before the run).** R1 = R_Au/AlN and R2 = R_AlN/sapphire on a 7 × 7 grid {2, 3, 5, 7, 10, 14, 20} ×
  10⁻⁹ m²K/W; c ∈ {1/2, 2/3, 3/4, 1, 4/3, 3/2, 2} with τ_R, λ_z, λ_r at (L, c) from R-A2 and w(L/c) from R-A3;
  λ_z variants 'event' (A1), 'rta' (τ_R and λ_z both from RTA 64×64×35, the A4 sensitivity case), λ_z × 0.8,
  λ_z × 1.2; L = 30, 100, 300, 500 nm; σ_φ = 0.01, 0.03, 0.1°; N1, N2, N3 as declared. 5 488 Jacobians; A4 set-up
  otherwise (finite beams 10/8 µm, Cattaneo AlN, 80 log-spaced frequencies 10 kHz–200 MHz, 1 % amplitude noise).
  Margin m = σ(ln τ_R)/w (1σ against the half-width, the A4 convention); m ≤ 1 is called a flip.
* **Implementations.** Production: the independent quadrupole model of branch A (`indep.py`, Gauss–Legendre Hankel
  rule) with the base values as traced arguments, exact forward-mode Jacobian (JAX), marginal σ by Householder-QR
  projection on the nuisance span (not the branch's SVD routine). Check: the twin (`src/twin_fdtr.py` through
  `common.py`, Richardson central differences) at the 32 corners (R1, R2 ∈ {2e-9, 2e-8}, c ∈ {1/2, 2}, four L),
  all sets and noise levels: max relative difference **9.9e-8**. A4 reproduced at the note-16 point to ≤ 2.2e-6.
  Peak memory 0.46 GB, 474 s.

**Margins** (widths robust 7×7×5/9×9×6; base value for the orb 0.1 class first and for the union of classes in
brackets where they differ; box minima for the union of classes). "base" = note-16 stack (R1/R2 = 1e-8/2e-8), c = 1, event values:

| L | set | 0.01°: base / box minimum (bits at the minimum) | 0.03°: base / min | 0.1°: base / min | flips at 0.01° |
|---|---|---|---|---|---|
| 30 nm | N1 | 455 (248) / 27.5 (3e-4) | 738 / 82 | 4.2e3 (2.3e3) / 253 | none |
| | N2 | 6.9e3 (3.8e3) / 92 | 1.1e4 / 268 | 4.6e4 (2.5e4) / 709 | none |
| | N3 | 1.9e4 (1.1e4) / 251 | 2.3e4 / 645 | 6.1e4 (3.3e4) / 1.7e3 | none |
| 100 nm | N1 | 100 (48) / 6.2 (6e-3) | 142 / 18 | 861 (415) / 57 | none |
| | N2 | 156 (75) / 17 (8e-4) | 221 / 51 | 1.3e3 (627) / 142 | none |
| | N3 | 318 (153) / 97 | 439 / 271 | 2.3e3 (1.1e3) / 400 | none |
| 300 nm | N1 | 9.3 (7.0) / 1.93 (0.062) | 21 / 5.7 | 86 (65) / 17 | none |
| | N2 | 13.1 (9.8) / 2.81 (0.030) | 29 / 8.3 | 106 (80) / 24 | none |
| | N3 | 52 (39) / 12.0 | 112 / 33 | 382 (286) / 79 | none |
| 500 nm | N1 | 3.35 / **0.78** (0.31) | 10.0 / 2.32 | 31.7 / 7.0 | 4.2 % of the box |
| | N2 | 5.55 / **0.92** (0.24) | 16.3 / 2.73 | 44.6 / 8.4 | 0.7 % of the box |
| | N3 | 36.8 / 9.2 | 105 / 27 | 278 / 75 | none |

* **Robust region** (m > 1 over the whole box): every (L, set, σ_φ) except 500 nm at 0.01° with N1 or N2; at 0.03° the
  minimum is 2.3 (500 nm, N1), at 0.1° it is 7.0 (500 nm, N1), 8.4 (N2), 75 (N3), and ≥ 17 at 300 nm.
* **Flip region** (500 nm, 0.01° only). N2: R1 ≤ 3e-9 m²K/W together with c ≤ 2/3 and λ_z − 20 % (minimum 0.92 at
  R1 = R2 = 2e-9, c = 1/2). N1: R1 ≤ 5e-9 with c ≤ 3/4 (minimum 0.78). Below m ≈ 1.5 at 500 nm and 0.01°, N2 needs
  R1 ≤ 5e-9 (G₁ ≥ 200 MW m⁻²K⁻¹). At c = 1 nothing flips (N2 minimum 1.25, N1 1.05). The minimum always sits at the
  smallest c (longest memory time, largest width w(1000 nm) = 6.5 %) and lowest λ_z. With the held-out 12×12×8 widths
  the N2 minimum is 1.01 (no flip), N1 0.86 (1.3 % of the box).
* **Au thickness** (sub-box c = 1, 'event', full (R1, R2) grid, d_Au = 20, 50, 80, 100 nm): no flip; minima at 500 nm
  1.20 (N1, 0.01°), 1.42 (N2, 0.01°), 14 (N3, 0.01°), 13 (N2, 0.1°), all at R1 = R2 = 2e-9; at 300 nm 1.78 (N1) and
  3.9 (N2) at 0.01°.
* **Reading.** The 0.1° statement holds over the whole declared box, by ≥ 7 at 500 nm and ≥ 17 at 300 nm. The 0.01°
  statement holds on the note-16 stack (≥ 3.35 for N1, ≥ 5.55 for N2 at 500 nm) and fails only in a corner of the box
  that combines a top resistance ≤ 3–5e-9 m²K/W, a regulator normalisation c ≤ 3/4 and a lower λ_z; it is
  stack-dependent and marginal there (m ≈ 0.8–1.5, ≤ 0.31 bit). The c-optimal allocation of the 80 points (audit §7,
  gain 1.3–2.9) and the 2σ reading (σ ≤ 2w, audit F2) are not folded in. Either one moves the 0.01° corner further
  towards resolvability; at 0.1° either one alone leaves the box minimum ≥ 2.4, both together ≈ 1.2.

## R-A5 — information measures (A-m1)  [R_A/ra_info.py; used in ra4_box.py]

Prior x ~ U[−w, w] on the relative change of τ_R, Gaussian measurement of standard deviation σ (the response is linear
in τ_R at these changes, A5): mutual information I (bits) and the fraction of the prior variance removed depend on
r = σ/w only. Quadrature in y (composite Gauss–Legendre, cancellation-free truncated-normal variance) against Monte
Carlo (4e5 samples, fixed seed): agreement within the Monte Carlo standard error for r = 0.1 … 10; large-r limit
I → 1/(6 ln 2 r²), removed → 1/(3r²) used for r > 100 (error < 4e-5); review R5 reproduced (r = 1: 0.2073 bit,
25.2 %; r = 4.25: 0.0132 bit, 1.8 %).

| r = σ/w | 0.5 | 1 | 2 | 3.35 | 5.55 | 13.1 | 44.6 |
|---|---|---|---|---|---|---|---|
| bits | 0.604 | 0.207 | 0.058 | 0.021 | 0.0078 | 0.0014 | 1.2e-4 |
| prior variance removed | 58.7 % | 25.2 % | 7.7 % | 2.9 % | 1.07 % | 0.19 % | 0.017 % |

On the note-16 stack at 500 nm: N1 0.01° removes 2.9 % of the prior variance (0.021 bit), N2 0.01° 1.07 % (0.0078 bit),
N2 0.1° 0.017 % (1.2e-4 bit); at 300 nm N2 0.01° 0.19 %. At the box minima (500 nm, 0.01°): 0.31 bit (N1), 0.24 bit
(N2); at 0.1° ≤ 0.005 bit everywhere in the box.

## R-A6 — objection A-2: the FDTR-visible time δ = τ_R − τ_ℓ, its normal/resistive content and its static-data range  [R_A/ra_kinetic.py, ra0_kinetic_tests.py, ra1_film_event_ext.py, ra1b_dense_check.py, ra1c_ref_rta.py, ra1d_freq_kernel.py, ra2b_rta_taul_meshes.py, ra5_delta.py, ra5b_observable.py; out/ra0_kinetic_tests.json, ra1b_dense_check.json, ra1c_ref_rta.json, ra1d_freq_kernel_m775.json, ra2b_rta_taul_meshes.json, ra5_delta.json, ra5b_observable.json]

**1. Is the normal/resistive split in the witness results?** No stored film result contains it (`results_c/film*/*.json`
hold τ, κ and τRTA only; no script of fc3-memory uses it). The event lists carry a per-event label ev_umk
(0 normal, 1 umklapp; BZ-representative test of `aln_events.py`, "reporting only") and the witnesses' amplitudes are
stored, so the split is computable without new optimisation and was computed here (R-A1). The label is not unique:
event orbits mixing normal and umklapp labels (BZ-boundary representatives) are 4 053 of 23 674 orbits (17.7 % of the
REF rate) on 7×7×5 and 30 415 of 95 384 (35.2 %) on 9×9×6; normal events carry 46 % (7×7×5) and 41 % (9×9×6) of the
REF rate. On 7×7×5 at 500 nm the label convention alone moves a gray normal time by −14 … +19 % (G1 31.0 ps, range
26.7–36.8 ps between "normal if any event of the orbit is normal" and "only if all are"; G2 11.9 ps, 9.9–14.8 ps).

**2. Two closures for τ_ℓ.**
* *Gray closure of note 14* (conventions §7.1): τ_ℓ = (4/5)τ_c, τ_c = (1/τ_N + 1/τ_R)⁻¹, τ_R = τ_mem, with a gray
  normal time from the normal events: G1 κ_z-weighted normal rate, G2 capacity-weighted normal rate, RQ Rayleigh
  quotient of the normal part of C on the even moment u = Q V x.
* *Kinetic closure* (DERIVATION, `ra_kinetic.py`): for the symmetrised linear BTE (p + ikV + Ĉ)f = 0 with Ĉe = 0 and
  b = Ve, eliminating the temperature with the energy balance gives the exact current kernel
  κ(k, p) = b·[p + Ĉ + ik QVQ]⁻¹b on e⊥ (Q = 1 − ê êᵀ). Its expansion κ = K[1 − τ_R p − L²k² + …] gives τ_R = |x|²/K
  (x = Ĉ⁻¹b; the memory time) and K L² = u·Ĉ⁻¹u, u = QVx; with a = K/|e|² from the same reduction the
  Guyer–Krumhansl time is τ_ℓ = L²/a = |e|²(u·Ĉ⁻¹u)/K², and the dispersion relation p|e|² + k²κ = 0 gives at first
  order in p the effective conductivity λ(1 − δp), δ = τ_R − τ_ℓ, i.e. exactly the 1D GK film of the twin. For the
  gray isotropic Callaway operator (one speed; normal events conserve energy and momentum, resistive ones energy)
  this reduces to τ_ℓ = (4/5)τ_c, the conserving closure of note 14; without the projection Q it gives (9/5)τ_c, the
  historical coefficient. In general τ_ℓ = γ_g τ_c,eff with the kinematic factor γ_g = ⟨v_z⁴⟩/⟨v_z²⟩² − 1
  (capacity-weighted; 4/5 for one isotropic speed). The regulator enters as the energy-conserving Q D Q (R-A1).
* *Checks.* Gray isotropic operator, τ_N/τ_R = 0.01–100: τ_ℓ/τ_c = 0.8 and the unprojected variant 1.8 to 4e-14
  (R-A0 T1). Random non-gray parity-symmetric operators with the same conservation laws: the k² coefficient of the
  exact kernel (complex solves at small k, Richardson) equals K L² to 1e-8 (T2); independence of the regularisation
  constant γ to 2e-15 (T3). AlN 7×7×5: a dense Householder reduction to e⊥ with LU reproduces τ_R and τ_ℓ of the
  production route to 4e-15 (REF and one optimum, five lengths, R-A1b); the direct small-k kernel of the AlN operator
  (R-A1d) checks K and the k² coefficient (values below).

**3. Reference values.** τ in ps; 9×9×6 unless stated; "RTA" = energy-conserving RTA kinetic time on the same
diagonal (event meshes) or on the converged phono3py mesh 64×64×35:

| L | τ_mem 7×7×5 / 9×9×6 | τ_ℓ kinetic 7×7×5 / 9×9×6 | τ_ℓ RTA 7×7×5 / 9×9×6 / 64×64×35 | τ_ℓ/τ_mem | τ_c,eff / τ_c,RQ | τ_N G1 / G2 / RQ | δ kinetic | δ gray (G1 / G2 / RQ) |
|---|---|---|---|---|---|---|---|---|
| 30 nm | 5.52 / 5.34 | 10.1 / 9.7 | 11.4 / 10.7 / 9.1 | 1.81 | 3.5 / 3.1 | 24.0 / 13.1 / 16.3 | −4.3 (−0.81 τ_mem) | +0.35 / 0.43 / 0.40 τ_mem |
| 100 nm | 11.50 / 10.89 | 31.8 / 32.4 | 35.9 / 34.3 / 32.3 | 2.97 | 11.6 / 6.0 | 28.0 / 13.1 / 23.2 | −21.5 (−1.97) | +0.42 / 0.56 / 0.46 |
| 300 nm | 19.84 / 18.82 | 104.8 / 118.6 | 131.9 / 135.2 / 134.5 | 6.30 | 42.4 / 11.4 | 32.0 / 13.1 / 36.1 | −99.7 (−5.3) | +0.50 / 0.67 / 0.47 |
| 500 nm | 24.58 / 23.70 | 178.1 / 217.9 | 234.9 / 256.1 / 266.1 | 9.20 | 78.0 / 15.4 | 33.8 / 13.1 / 44.8 | −194.2 (−8.2) | +0.53 / 0.72 / 0.48 |
| 1000 nm | 31.38 / 31.36 | 334.1 / 460.4 | 464.6 / 558.4 / 648.8 | 14.7 | 165 / 22.5 | 36.1 / 13.1 / 58.3 | −429 (−13.7) | +0.57 / 0.76 / 0.48 |

γ_g = 3.03 (7×7×5), 2.79 (9×9×6), 2.73 (64×64×35).
* **Convergence.** The RTA kinetic time on the three converged phono3py meshes 48×48×26 / 56×56×30 / 64×64×35 is
  9.16 / 9.12 / 9.08 ps (30 nm), 32.28 / 32.27 / 32.28 (100 nm), 133.9 / 134.2 / 134.5 (300 nm), 264.4 / 265.8 /
  266.1 (500 nm): τ_ℓ(L) is a converged film quantity like τ_mem(L) (R-A2b). On the event meshes the RTA value
  converges towards it (500 nm: 234.9 → 256.1 → 266.1), and the full-operator value moves with it (178 → 218 ps);
  full/RTA on the same mesh is 0.76 (7×7×5) and 0.85 (9×9×6). Mesh-corrected full value at 500 nm ≈ 226 ps.
* **What carries τ_ℓ** (RTA 64×64×35, 500 nm): modes with |v_z|/v_max = 0.2–0.5 (76 %) and bulk lifetimes 100–1000 ps
  (64 %), i.e. boundary-limited cross-plane carriers; the plane q_z = 0 carries 0.2 % (6.5 % at 30 nm). With an
  isotropic regulator |v|/L instead of |v_z|/L: τ_ℓ(500 nm) = 198.7 ps (−25 %), τ_mem 20.5 ps.
* **Why τ_ℓ ≫ τ_R, whereas the gray closure bounds τ_ℓ < (4/5)τ_R.** (i) γ_g = 2.7–3.0 instead of 0.8 (broad
  velocity distribution); (ii) the even moment u relaxes on τ_c,eff = τ_ℓ/γ_g, an inverse-rate (slow-mode) average,
  five times the rate average τ_c,RQ that a gray τ_c represents (500 nm: 78 against 15 ps), and longer than τ_mem
  itself at L ≥ 100 nm. The gray closure with any of the three normal times gives δ = +(0.35–0.76) τ_mem; the
  kinetic expansion of the same operator gives δ = −(0.8–14) τ_mem (30–1000 nm). In the kinetic closure the first-order non-Fourier
  signature of the film is nonlocal and has the opposite sign of the Cattaneo one.

**4. Static-data range of δ** (inner bounds; witness excursions robust over 7×7×5/9×9×6 as in R-A3; absolute
half-widths in units of τ_mem,0 with the 9×9×6 REF ratio |δ₀|/τ_mem,0):

| L | class | w(τ_mem) | w(τ_ℓ) | w(δ), kinetic | half-range of δ (kinetic) / τ_mem,0 | kinetic / w(τ_mem) | w(τ_N) G1 / G2 / RQ | half-range of δ (gray G1 / G2) / τ_mem,0 | gray / w(τ_mem) |
|---|---|---|---|---|---|---|---|---|---|
| 30 nm | orb 0.1 | 2.72 % | 2.0 % | 7.7 % | 0.063 | 2.3 | 3.0 / 4.3 / 4.2 % | 0.010 / 0.014 | 0.37 / 0.50 |
| | all | 5.00 % | 2.0 % | 9.4 % | 0.077 | 1.5 | 3.0 / 4.3 / 4.2 % | 0.023 / 0.029 | 0.46 / 0.57 |
| 100 nm | orb 0.1 | 3.64 % | 5.6 % | 9.8 % | 0.193 | 5.3 | 2.4 / 4.3 / 2.6 % | 0.019 / 0.027 | 0.53 / 0.73 |
| | all | 7.56 % | 5.6 % | 10.3 % | 0.202 | 2.7 | 2.4 / 4.3 / 3.7 % | 0.043 / 0.056 | 0.57 / 0.74 |
| 300 nm | orb 0.1 | 5.04 % | 10.9 % | 13.1 % | 0.692 | 13.7 | 1.6 / 4.3 / 5.4 % | 0.037 / 0.051 | 0.73 / 1.00 |
| | all | 6.75 % | 12.3 % | 15.2 % | 0.806 | 12.0 | 1.6 / 4.3 / 7.0 % | 0.045 / 0.057 | 0.67 / 0.85 |
| 500 nm | orb 0.1 | 6.13 % | 14.0 % | 15.5 % | **1.27** | **21** | 1.2 / 4.3 / 7.1 % | 0.046 / 0.062 | 0.75 / 1.01 |
| | all | 6.13 % | 17.2 % | 19.8 % | **1.63** | **27** | 1.2 / 4.3 / 8.8 % | 0.046 / 0.062 | 0.75 / 1.01 |

* The τ_ℓ and δ excursions of the witnesses shrink from 7×7×5 to 9×9×6 by 20–45 % (cal3w min at 500 nm: +41.5 →
  +23.1 %), as their bulk excursions do (E36), whereas their τ_mem excursions at the optimised lengths are stable
  between 9×9×6 and 12×12×8 (E37); no 12×12×8 value of τ_ℓ exists (no 12×12×8 operator computation in this revision).
  τ_ℓ(500 nm) of a witness follows its bulk memory time (9×9×6: 30 nm max τ_ℓ −13.6 %, bulk τ_z −12.2 %; 30 nm min
  +14.3 %, +15.7 %), so the kinetic δ inherits part of the freedom that leaves the bulk memory time undetermined.
* Gray closure: the witnesses move τ_N by 1.2–3.0 % (G1), 4.3 % (G2) and 2.6–8.8 % (RQ). With x = τ_N/τ_R,
  ∂δ/∂τ_R = 1 − 0.8x²/(1+x)² ∈ (1/5, 1] and (τ_N/τ_R) ∂δ/∂τ_N = −0.8x/(1+x)² ∈ [−1/5, 0), so
  |Δδ|/τ_R ≤ w_R + w_N/5; the computed gray half-ranges are 0.37–1.01 times the τ_mem half-range. Within the gray
  closure the A comparison is conservative (equal within 1 % at worst, G2 at 300–500 nm).

**5. Frequency dependence in the band** (R-A1d, 7×7×5, all 13 tensors, 300 and 500 nm). Up to 200 MHz the in-phase
part of δ_eff(ω) = −(λ_eff/λ₀ − 1)/(iω) stays within 2.3 % of its static value for every tensor at 300 and 500 nm (REF,
500 nm: −153.5 ps static, −150.7 ps at 200 MHz); the out-of-phase part, the second-order term, reaches 11–17 % of |δ|
at 200 MHz (REF, 500 nm: 23.9 ps). Witness/REF ratios of δ_eff are preserved in the band to ≤ 0.5 % (30 nm max at
500 nm: 0.814 static, 0.815 at 200 MHz). Above the band the first-order kernel departs (1 GHz: in-phase part −36 %).
The direct small-k kernel of the AlN operator (7×7×5, REF, 500 nm) reproduces K to 8.9e-16 and the k² coefficient
to 3.7e-5 (Richardson in k). The static pair (τ_mem, τ_ℓ) is an adequate first-order description of the regularised
kernel in the FDTR band.

**6. FDTR consequences (1D large-beam limit).** The twin stack in 1D (A1 base λ_z and C; the witnesses change only
the dynamic part), local interface temperature unless stated. Closed form checked against the
twin's k = 0 kernel to 1.8e-15 (local) and 1.6e-15 (entropic convention).
* σ(δ) at the kinetic base point (τ_mem, τ_ℓ = mesh-corrected kinetic value) equals σ(τ_R) of the Cattaneo film to
  within 6 % for N1/N2 and is 1.4 times larger for N3 (500 nm, 0.01°: 7.16 against 7.61 ps for N2, 69.5 against
  49.0 ps for N3); along −τ_ℓ at fixed τ_R it is within 8 % of the τ_R direction (first-order dependence on δ only).
  The response is linear in δ to 18 % (500 nm) and 14 % (300 nm) between δ = τ_R and the kinetic δ₀.
* **Margins** σ(δ)/half-range of δ, all classes, note-16 stack (9×9×6-robust ranges); in brackets the direct test of
  the observable with the frequency-dependent 7×7×5 kernel of each witness (R-A5b; m = 1/√max Δχ², largest single
  excursion; the 7×7×5 excursions are 20–45 % larger than on 9×9×6):

| L | set | kinetic 0.01° | 0.03° | 0.1° | Cattaneo (same 1D set-up) 0.01° / 0.1° |
|---|---|---|---|---|---|
| 300 nm | N2 | 0.68 (0.40) | 2.0 | 5.6 (3.3) | 8.7 / 71 |
| | N3 | 3.7 (2.3) | 10.6 | 27 (17) | 34 / 234 |
| 500 nm | N2 | **0.18** (0.10) | 0.52 | **1.46** (0.87) | 5.0 / 40 |
| | N3 | 1.7 (0.97) | 4.8 | 12.5 (8.0) | 32 / 228 |

  With R1 = R2 = 2e-9 m²K/W the 500 nm N2 margins fall to 0.05 (0.01°) and 0.49 (0.1°). With the entropic interface
  temperature (θ = (1 + τ_ℓp)T at both AlN interfaces; twin GKBoundary 'entropic') the direct-test margins for N1/N2
  are 1.7–4.4 times smaller still (500 nm, N2: 0.053 at 0.01°, 0.52 at 0.1°; 300 nm, N2: 0.091, 0.90); for N3 they
  change by ≤ 20 %. The frequency-dependent and the static kernels give the same margins to ≤ 6 % (500 nm, N2: 0.101
  against 0.106, 0.867 against 0.890).
* **Fourier test (A7 analogue, 1D, note-16, nuisances refitted by Gauss–Newton).** Δχ² between the kinetic GK film and
  a Fourier film: 500 nm N2 798 (0.01°) / 11.7 (0.1°), N3 11.7 / 0.21; 300 nm N2 82 / 1.2, N3 4.0 / 0.08; 100 nm N2
  0.18 / 0.002. The Cattaneo film in the same set-up: 500 nm N2 10.8 / 0.16 (A7 with finite beams: 8.62 / 0.133).
* **Reading.** Within the gray closure the comparison of A is conservative. Within the kinetic closure of the same
  regularised operator it is not: the FDTR-visible first-order time is the nonlocal one, the witnesses' range of δ is
  12–27 times the τ_mem range at 300–500 nm, FDTR would resolve it at 0.01° (N1/N2, 300–500 nm), and at 0.1° the N2
  margin falls to 0.9–1.5 at 500 nm (3.3–5.6 at 300 nm); with λ_z and C free (N3) it stays ≥ 8 at 0.1°. τ_ℓ = 0 is
  not justified for these films (τ_ℓ = 1.8 τ_mem at 30 nm, 9 τ_mem at 500 nm), and the GK variant τ_ℓ = 2τ_R of A4 is
  not outside the physical range: the bound τ_ℓ < (4/5)τ_R is a property of the gray one-speed closure. Under the
  kinetic closure the A7 statement also changes: the first-order non-Fourier response of a 300–500 nm film is visible
  at 0.1° when λ_z and C are known (N2), not when they are free (N3).
* **Caveats.** First order in k² (L_GK = (τ_ℓ a)^½ ≈ 130 nm at 500 nm, k L_GK ≈ 0.5 at 200 MHz); bulk-like kernel of
  the Matthiessen-regularised operator (no Knudsen layers; |v|/L instead of |v_z|/L changes τ_ℓ by −25 %); the
  interface-temperature convention enters at O(τ_ℓ p); FDTR side in 1D only; the τ_ℓ and δ excursions are inner
  bounds of a range never optimised for δ, and they shrink by 20–45 % from 7×7×5 to 9×9×6.
* Status: τ_ℓ formula DERIVED (exact first-order expansion; manufactured and direct checks). τ_ℓ(L) values NUMERICALLY
  SUPPORTED (RTA converged to < 0.5 % on three phono3py meshes; full operator on two event meshes, mesh-corrected).
  "A conservative under the gray closure": NUMERICALLY SUPPORTED. "A conservative under the kinetic closure":
  CONTRADICTED within the regularised model (1D FDTR side, inner-bound ranges).

## R-A synthesis (answers to A-1, A-2, A-m1, A-m3; supersedes the corresponding statements of the A-C synthesis)

* **A-1, width at 500 nm.** From the existing optima, without new optimisation: w(500 nm) ≥ 6.1 % on the data meshes
  ([−5.98 %, +6.28 %], orb 0.1 300 nm optima; the union of classes gives the same), 5.0–5.4 % held out on 12×12×8
  (interpolated). INNER BOUNDS. Against these the note-16 margins at 500 nm are 5.55 (N2, 0.01°) and 44.6 (N2, 0.1°)
  instead of 4.25 and 34 against the assumed 8 %; at 300 nm 13.1 / 106 (orb 0.1) and 9.8 / 80 (all classes). A 500 nm
  optimisation can only widen the range and lower these numbers; the trend of the optimised widths (5.95–6.12 %)
  is already reached and is not an upper estimate.
* **A-1, sensitivity box** (R-A4). σ(ln τ_R) > w holds over the whole declared box except at 500 nm and 0.01° with N1
  or N2 (flips in 4.2 % and 0.7 % of the box: R1 ≤ 3–5e-9 m²K/W with c ≤ 3/4 and low λ_z; minima 0.78 and 0.92,
  ≤ 0.31 bit). Minimum over the box at 0.1°: 7.0 (500 nm, N1), 8.4 (N2), 75 (N3); 17 / 24 / 79 at 300 nm; ≥ 56 at
  30–100 nm. The 0.1° statement is robust; the 0.01° statement is stack-dependent and marginal at 500 nm.
  The A-C synthesis wording "every nuisance set N0–N5" and "a margin of more than an order of magnitude at every
  thickness" is replaced by these numbers.
* **A-m3.** τ_mem(L; c) = τ_mem(L/c; 1) and w(L; c) = w(L/c; 1) for any operator (R-A1); τ_R(500 nm) = 17.3–35.4 ps for
  c = 2 … 1/2 (R-A2); widths at L/c in R-A3, margins over c in R-A4 (the minima sit at c = 1/2).
* **A-m1.** R-A5: on the note-16 stack FDTR removes ≤ 2.9 % of the prior variance (≤ 0.021 bit) at 500 nm and 0.01°,
  ≤ 0.034 % (≤ 2.4e-4 bit) at 0.1°; ≤ 0.31 bit anywhere in the box at 0.01°, ≤ 0.005 bit at 0.1°.
* **A-2** (R-A6). The normal/resistive split is not in the stored results but is computable from the stored inputs
  (per-event labels, ambiguous for 18–35 % of the rate). (i) Gray conserving closure (τ_ℓ = 4τ_c/5) with three gray
  normal times: δ = +(0.35–0.76) τ_mem and the visible range of δ is 0.37–1.01 times the τ_mem range: the comparison
  of A is conservative within that closure. (ii) Kinetic closure (exact first-order expansion of the same regularised
  operator; it reduces to (i) for a gray isotropic operator): τ_ℓ = 1.8 τ_mem (30 nm) to 9.2 τ_mem (500 nm), converged
  in RTA to < 0.5 %; δ < 0; its static-data range (inner bounds) is 1.5–2.3 (30 nm), 2.7–5.3 (100 nm), 12–14 (300 nm)
  and 21–27 (500 nm) times the τ_mem range in absolute terms. In 1D with the local interface temperature this range
  would be resolved at 0.01° (N1/N2, 300–500 nm; margin 0.18 at 500 nm); at 0.1° the N2 margin is 0.9–1.5 (500 nm)
  and 3.3–5.6 (300 nm), with N3 ≥ 8. The comparison of A is therefore not conservative beyond the Cattaneo
  projection, and τ_ℓ = 0 cannot be justified for these films.
* **Headline that the record supports.** "In the Cattaneo projection of the conserving closure (τ_ℓ = 0), the
  static-data range of the film memory time — per-orbit and calibrated fc3 classes, three data meshes, inner bounds,
  half-width 2.7–7.6 % at 30–500 nm — is not resolvable by FDTR in the declared Au/AlN/sapphire scenario at 0.1° phase
  noise for any declared nuisance set over interface resistances 2–20 × 10⁻⁹ m²K/W, regulator normalisations 1/2–2,
  λ_z ± 20 % and Au 20–100 nm (σ(ln τ_R)/w ≥ 7; ≤ 0.005 bit); at 0.01° it is not resolvable on the note-16 stack
  (≥ 3.3), and at 500 nm the margin falls to 0.8–1.5 for top resistances ≤ 5 × 10⁻⁹ m²K/W, below 1 only together with
  a regulator normalisation c ≤ 3/4." Required qualifier: the first-order
  FDTR-visible time of the same regularised operator is dominated by a nonlocal time τ_ℓ ≈ 6–9 τ_mem at 300–500 nm,
  whose static-data range is 12–27 times wider and which would be resolvable at 0.01°; the statement is about the
  Cattaneo projection, not about the film's full first-order response.
* **Open.** δ-optimised ranges (objective τ_mem − τ_ℓ in robust_mesh.py); τ_ℓ on 12×12×8 (the RTA diagonal is light);
  a film solution with diffuse walls in place of the |v_z|/L regulator; the interface-temperature convention implied
  by the closure (it enters at O(τ_ℓ p)); finite beams for the kinetic GK film.
