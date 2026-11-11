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
