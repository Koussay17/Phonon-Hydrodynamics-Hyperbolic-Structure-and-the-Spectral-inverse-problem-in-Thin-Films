# Branch R-A-K — independent audit of the kinetic nonlocal time and of the δ-range result (R-A6, objection A-2)

12 November 2026. Object: `07-experiments.md` R-A, R-A1…R-A6 and the R-A synthesis; `10-proof-status.md` (Revision R-A);
`12-peer-review.md` objection A-2; scripts `experiments/R_A/` (ra_kinetic.py, ra0_kinetic_tests.py, ra1_film_event_ext.py,
ra1b/c/d, ra2b_rta_taul_meshes.py, ra5_delta.py, ra5b_observable.py) and their outputs; upstream operator data of
fc3-memory (events, amplitudes, RTA chunks; read-only); `src/twin_fdtr.py` (read-only); note 14 and conventions §7.
Audit scripts `experiments/R_A_audit/au1…au7*.py`, raw outputs `experiments/R_A_audit/out/` (list in §9). Research
Python, two threads, one process at a time, peak working set below 1 GB (dense 7×7×5 operators ≈ 0.4 GB); the 12×12×8 operator was not touched and the 9×9×6
operator was not rebuilt (its stored per-witness values are used for the ranges only). Every number is SYNTHETIC and
conditional on the twin physics, the noise model and the declared operators.

**Claim audited.** For the film-regularised conserving event operator (regulator |v_z|/L), the first-order (k², p)
expansion gives τ_ℓ = |e|²(u·Ĉ⁻¹u)/K², u = QVx; it reduces to 4τ_c/5 for the grey isotropic Callaway operator; for AlN
τ_ℓ/τ_mem ≈ 1.8 / 3.0 / 6.3 / 9.2 at 30/100/300/500 nm (RTA version mesh-converged < 0.5 %); the static-data range of
δ = τ_R − τ_ℓ is 12–27× wider than that of τ_mem at 300–500 nm and would be resolvable by FDTR at 0.01° (1D, local
convention).

## 0. Verdicts

| # | sub-claim | verdict | decisive numbers |
|---|---|---|---|
| 1 | τ_ℓ = \|e\|²(u·Ĉ⁻¹u)/K², u = QVx, is the first-order coefficient and equals the twin's τ_ℓ = (1+α)ℓ²C/λ | **RELIABLE** (derivation correct; caveat on the even sector) | three independent routes: exact rational grey case; bordered-system kernel of random non-grey operators (K 1e-14, τ_R 4e-12, k² coefficient ≤ 3e-7, δ from the dispersion relation ≤ 2e-6); twin `Material.tau_l` and k = 0 kernel vs an independent quadrupole 7e-16. The projection Q is required. τ_ℓ depends on how the regulator acts on the even sector, which no static (odd-sector) datum fixes: +38 / +12 / +4.0 / +2.5 % at 30/100/300/500 nm between the campaign's Q D Q and relaxation to the wall temperature |
| 2 | grey limit 4τ_c/5 (9τ_c/5 without Q) | **RELIABLE** | exact in SymPy for all τ_N/τ_R; also δ = τ_R − 4τ_c/5 from the exact dispersion relation |
| 3 | τ_ℓ/τ_mem ≈ 1.8/3.0/6.3/9.2; RTA mesh-converged < 0.5 % | **RELIABLE WITH CAVEATS** (RTA: CONVERGED; full operator: UNDER-RESOLVED) | RTA reproduced exactly on 48×48×26 and 64×64×35 by three routes (moment formula; operator algebra; Cauchy-integral kernel coefficients, ≤ 2e-11); refinement steps −0.9 / 0.0 / +0.4 / +0.6 %. Full operator reproduced exactly on 7×7×5 (13 tensors, 5 lengths, ≤ 4e-15) by an independent assembly. The quoted ratios are 9×9×6 full-operator values; τ_ℓ moves +13 % (300 nm) and +22 % (500 nm) from 7×7×5 to 9×9×6; the mesh-corrected value carries ≈ ±10 % |
| 4 | static-data range of δ 12–27× that of τ_mem (300–500 nm) | **RELIABLE WITH CAVEATS** (arithmetic exact; UNDER-RESOLVED in mesh) | reproduced exactly (11.95 / 26.52, all classes) from the raw per-witness values; the robust rule returns the 9×9×6 excursions; τ_ℓ excursions shrink by a factor 0.70–0.72 per refinement at 300–500 nm; a geometric continuation of that trend gives 9.7× (300 nm) and 19.9× (500 nm). The ratio of two inner bounds is not a bound. N/U labels do not enter the kinetic closure. The δ-extreme witnesses also move the film κ_z by −4.4 … +3.2 % (anticorrelated with τ_ℓ) |
| 5 | resolvable at 0.01° (1D, local) | **NOT RELIABLE** (model-dependent; holds only in part of the readings) | own Fisher reproduces σ and margins to 3 digits (500 nm N2: 0.176 at 0.01°, 1.46 at 0.1°; N3: 1.71, 12.5). 500 nm, 0.01°, note-16: first-order regularised model 0.18 (N2) / 0.51 (λ_z free, C known) / 1.71 (N3); with the mesh continuation 0.27 / 0.77 / 2.60; with the k² truncation (N2) 0.39. Kinetic slab (visible range × 0.26 from the witnesses' own lifetimes): 1.82 (λ_z free) / 5.4 (N3), and resolvable only if λ_z is known to ≤ 0.5 %. In the slab the N1/N2 premise with λ_M has no meaning (GK description needs λ_z ≈ 1.24 λ_M). 300 nm: ≥ 5 in the slab reading. At 0.1° no reading resolves on the note-16 stack |
| 6 | first order in k² is adequate in the FDTR band ("k L_GK ≈ 0.5") | **NOT RELIABLE** above ≈ 20 MHz at 300–500 nm | exact diffusive root of the same regularised RTA kernel (64×64×35): the non-Fourier part of λ_eff departs from the first-order GK form by 6 / 28 / 47 / 67 % at 10/50/100/200 MHz (500 nm); c₄/(K L⁴) = 4.8, c₆/(K L⁶) = 32 (GK resummation: 1, 1); \|k\|Λ_max ≈ 2 at 200 MHz, outside the convergence radius. FDTR: the first-order film misfits the exact-kernel data by Δχ² = 70 (N2, 0.01°, note-16; 1748 favourable). R-A1d tested only the p-dependence at first order in k² |
| 7 | the Matthiessen-kernel film model for a quantity carried by boundary-limited carriers | **NOT RELIABLE** for τ_ℓ and δ (factor ≈ 3) | RTA slab with diffuse black walls, same mode set, exact modal solution (verified). The FDTR-visible δ is −76 ps at 500 nm and −36 ps at 300 nm against δ_M = −240 / −113.5 ps (ratio 0.32, converged in binning to ±0.002, stable in port convention and stack). A uniform lifetime change transfers at 0.34–0.38; the witnesses' own lifetime changes transfer their δ range at 0.26–0.31 (both event meshes). The GK film with the regularised τ_ℓ is rejected by Δχ² = 48 (λ_z free, C known, 0.01°), and only by 4.4 with C also free |

**Overall.** Items 1–2 hold. Item 3 holds for the RTA route and is an under-resolved 9×9×6 number for the full operator.
Item 4 is correct arithmetic on 12 non-optimised witnesses whose δ excursions are not mesh-converged. Item 5 does not
follow once the model is tested; inside the first-order regularised GK film the margins are reproduced. At 500 nm and
0.01° the regularised first-order model resolves the range for N1/N2 and also with λ_z free when C is known (0.51), but
not with C free (1.71). A kinetic treatment of the walls shrinks the visible range by a factor ≈ 3.8 (witnesses' own
lifetimes). After that, resolution at 0.01° on the note-16 stack requires λ_z known to ≤ 0.5 %, or interface
resistances near 2e-9 m²K/W. At 300 nm, and at 0.1° on the note-16 stack, no reading resolves it. The first-order k²
truncation is itself inadequate at 0.01°. The exact regularised kernel moves the margins in both directions depending
on the nuisance set (×1.3–1.8 larger for N2; for N3 at 500 nm ×0.67–0.94 on the note-16 stack and ×0.35–1.33 on the
favourable one, the latter through the amplitude channel).
The robust content of R-A6 is qualitative: in every model tested, the first-order FDTR-visible time of a 300–500 nm AlN
film is negative and several times larger in magnitude than τ_mem. The Cattaneo projection is therefore not a
conservative stand-in for the film's first-order response, and τ_ℓ = 0 is not justified. The size of that time (−36 /
−76 ps in the slab, −115 / −241 ps in the regularised kernel) and of its static-data range are model-dependent at the
factor-3 level.

## 1. The formula (sub-claims 1–2) — `au1_derivation.py`, `out/au1_derivation.json`

**Setting (ASSUMPTIONS).** Symmetrised linear BTE in one cross-plane coordinate, (p + ikV + Ĉ)f = s, Ĉ = Ĉᵀ ≥ 0,
ker Ĉ = span{e} (resistive processes present; the regulator guarantees it), parity P with PVP = −V, PĈP = Ĉ, Pe = e;
b = Ve (odd, b·e = 0), |e|² = C. Temperature θ = e·f/|e|² (the energy temperature, the T of the twin's energy equation).

**Derivation (independent of `ra_kinetic.py`).** With f = θe + g, g ⊥ e, the e-projection gives the energy balance
p|e|²θ + ik b·g = 0 (e·Ĉ = 0, e·Ve = 0). The Q-projection gives (p + Ĉ + ikQVQ)g = −ikθ b on e⊥ (QVe = b). Hence
q = b·g = −ikθ κ(k, p), κ = b·(p + Ĉ + ikQVQ)⁻¹b on e⊥, exactly. The Neumann series around Ĉ on e⊥ (invertible there)
gives κ = K − pKτ_R − k²(u·Ĉ⁻¹u) + O(p², pk², k⁴). Here K = b·x, x = Ĉ⁻¹b, τ_R = x·x/K, u = QVQx = QVx. The odd term
x·Vx vanishes by parity. The twin's 1D GK film, τ_R q_t + q = −λT_z + L²q_zz with L² = (1+α)ℓ² (the 1D reduction of
ℓ²[∇²q + α∇∇·q]), has the kernel λ/(1 + τ_Rp + L²k²) = λ[1 − τ_Rp − L²k² + …]. Matching gives λ = K, L² = u·Ĉ⁻¹u/K and
τ_ℓ = (1+α)ℓ²C/λ = L²/a = |e|²(u·Ĉ⁻¹u)/K². On the dispersion relation p|e|² + k²κ = 0, k² = −p/a + O(p²), so
λ_eff = −p|e|²/k² = K[1 − (τ_R − τ_ℓ)p + O(p²)]. That is the twin's λ(1+τ_ℓp)/(1+τ_Rp) to first order. DERIVED.

**The projection.** Ĉ⁻¹ exists only on e⊥, and Vx has the energy component ê·Vx = K/|e| ≠ 0. That component is the
temperature mode, already carried by the energy balance. Relaxing it at a collision rate (the "unprojected" variant)
treats the conserved energy mode as a relaxing one. In the grey case this produces the historical coefficient 9/5. u
must be projected, and only onto e⊥ (no momentum projection: ker Ĉ = span{e}).

**Checks.**
* A (exact rational arithmetic, grey isotropic Callaway, orthonormal Legendre basis, exact through O(k⁶)):
  K = τ_R/3 (v = 1), τ_R,kin = τ_R, τ_ℓ = 4τ_Nτ_R/[5(τ_N+τ_R)] = (4/5)τ_c, unprojected (9/5)τ_c. The k² coefficient of
  the exact kernel equals u·Ĉ⁻¹u. The exact dispersion relation gives δ = τ_R − (4/5)τ_c. The k⁴ coefficient is
  c₄ = K L⁴ (55τ_N + 28τ_R)/[28(τ_N + τ_R)], i.e. 1 to 1.96 times the GK value K L⁴.
* B (random non-grey parity-symmetric operators, 60–240 modes, built in the parity-adapted basis, with and without a
  normal-type part that conserves an odd momentum not parallel to the current). The exact kernel comes from the bordered
  system [[p + ikV + Ĉ, −e], [eᵀ, 0]] (source along e only, θ = 1; no projection, no regularisation constant), and the
  closed formula from a KKT solve. Agreement: K 1.2e-14, τ_R 4.4e-12, k² coefficient (Richardson) ≤ 3.0e-7, δ from the
  Newton-solved diffusive root of the dispersion relation ≤ 1.9e-6. τ_ℓ/τ_R ranges 0.21–0.55 in these cases.
* C (twin convention). `twin_fdtr.Material.tau_l` = (1+α)ℓ²C/λ to 1e-16. The twin's k = 0 kernel of
  Au/R1/GK AlN/R2/sapphire equals an independent quadrupole product with λ_eff = λ(1+τ_ℓp)/(1+τ_Rp) to 6.8e-16 (local)
  and 7.3e-16 (entropic, θ = (1+τ_ℓp)T at both AlN faces).

**Even-sector completion (a caveat the record does not state).** Only the odd sector of the regulator is fixed by
the static data (on b and Ĉ⁻¹b every energy-conserving completion equals D). τ_ℓ also needs Ĉ⁻¹ on the even vector u,
where the completions differ. For two of them, the quadratic forms differ by a rank-one term ≥ 0: for g ∈ e⊥,
g·(QDQ)g − g·[D − (De)(De)ᵀ/(eᵀDe)]g = (eᵀDg)²/(eᵀDe). Hence τ_ℓ(QDQ) ≤ τ_ℓ(relaxation to the wall temperature).
The latter is the completion that corresponds to diffuse re-emission of boundary-scattered phonons at the flux-weighted
local temperature (gain ∝ D e). Computed:
* 7×7×5 full operator, REF (`au3_event775.py`): τ_ℓ(phys)/τ_ℓ(QDQ) = 1.379 / 1.121 / 1.040 / 1.025 / 1.012 at
  30/100/300/500/1000 nm. At 300–500 nm the witness excursion ranges change by ≤ 0.5 point for τ_ℓ and ≤ 0.9 point for δ; at 30 nm the δ range
  narrows (−10.4…+9.1 % → −6.9…+6.0 %).
* RTA 64×64×35 (`au2_rta.py`, τ_ℓ in ps at 30/100/300/500/1000 nm): V1 one local temperature (campaign) 9.08 / 32.28 /
  134.47 / 266.13 / 648.79; V2 separate local temperatures for the bulk and the wall processes 11.72 / 36.82 / 141.72 /
  274.93 / 659.56 (+29 / +14 / +5.4 / +3.3 / +1.7 %); V3 bulk local temperature + QDQ wall 9.30 / 33.98 / 138.08 / 270.76
  / 654.66; V4 QDQ for both 8.97 / 31.59 / 125.49 / 245.19 / 593.71.
At 300–500 nm the completion moves τ_ℓ by −8 … +5 % (V4 … V2; full operator +2.5 … +4 %), at 30–100 nm by up to
+29 … +38 %. State the completion with the number.

**Convergence radius.** The k-series of κ is controlled by the longest cross-plane mean free path of the regularised
operator, Λ_max = max |v_z|/rate (RTA: poles of the static kernel at k = ±i/Λ_i). For 64×64×35, Λ_max = 297 nm (300 nm)
and 491 nm (500 nm), against L_GK = (u·Ĉ⁻¹u/K)^½ = 91 nm and 136 nm. The series is therefore a series in (kΛ_max)², not
in (kL_GK)² (§5.1).

## 2. τ_ℓ values (sub-claim 3) — `au2_rta.py`, `au3_event775.py`

**RTA route, two converged phono3py meshes, three routes.** (i) The moment formula ⟨τ(v_z²τ − s)²⟩/s², s = ⟨v_z²τ⟩.
(ii) The operator algebra with Woodbury solves on e⊥ (γ-independence 7e-16, |ê·x|/|x| ≤ 4e-15). (iii) Taylor coefficients of the exact
driven kernel |e|²(b·A⁻¹e)/(e·A⁻¹e)/(−ik), extracted by a Cauchy integral on |k|²Λ_max² = 1/4 (64 points).
Agreement: K ≤ 7e-16, k² coefficient ≤ 1.5e-13, τ_R ≤ 2.1e-11.

| L (nm) | τ_mem 48×48×26 / 64×64×35 | τ_ℓ 48×48×26 | τ_ℓ 64×64×35 | step | campaign (R-A2b) |
|---|---|---|---|---|---|
| 30 | 5.295 / 5.303 | 9.16 | 9.08 | −0.9 % | 9.16 / 9.08 |
| 100 | 10.728 / 10.747 | 32.28 | 32.28 | 0.0 % | 32.28 / 32.28 |
| 300 | 19.446 / 19.511 | 133.88 | 134.47 | +0.4 % | 133.88 / 134.47 |
| 500 | 25.451 / 25.567 | 264.42 | 266.13 | +0.6 % | 264.42 / 266.13 |
| 1000 | 36.005 / 36.257 | 641.5 | 648.79 | +1.1 % | 641.5 / 648.79 |

RTA ratios τ_ℓ/τ_mem (64×64×35): 1.71 / 3.00 / 6.89 / 10.41. "RTA mesh-converged to < 0.5 %" is correct at
100–500 nm (0.9 % at 30 nm over two refinements). CONVERGED.

**Full operator, 7×7×5, independent assembly.** Modes, event vectors, rates and orbit tying were re-derived from the raw
event list. The event-image graph under the 24 stored maps has 23 674 connected components, as in the campaign, with 0
missing images. C = Σ g a aᵀ is assembled as a COO sum, and the solves use the KKT bordered system with LU (no padded
Cholesky, no Householder). For all 13 tensors × 5 lengths: τ_z ≤ 2.0e-15, τ_ℓ ≤ 3.8e-15, κ_z ≤ 1.6e-15 against
`ra1_film_m775.json`; energy residual ‖Ce‖/‖diag(C)e‖ ≤ 8.7e-16.

**Mesh status of the quoted ratios.** 1.81 / 2.97 / 6.30 / 9.20 are τ_ℓ/τ_mem of the full operator on 9×9×6 (e.g.
217.9/23.70 at 500 nm). They are not mesh-converged values. τ_ℓ moves from 7×7×5 to 9×9×6 by −4 / +2 / +13 / +22 /
+38 % at 30/100/300/500/1000 nm, and the full/RTA ratio on the same mesh moves 0.76 → 0.85 at 500 nm. The mesh-corrected
226 ps assumes the 9×9×6 full/RTA ratio is final. Its uncertainty is not smaller than the last observed change of that
ratio (≈ ±10 %). Defensible reading: τ_ℓ/τ_mem ≈ 9 ± 1 at 500 nm and ≈ 6–7 at 300 nm within the regularised model
(QDQ completion), plus the completion spread of §1.

## 3. The δ range (sub-claim 4) — `au4_delta_range.py`, `out/au4_delta_range.json`

* **Arithmetic.** Recomputed from the per-witness values of `ra1_film_m775/m996.json` with the campaign rule (smaller
  excursion if same sign, else 0; ends = largest robust excursion of each sign; w = mean of the ends). Ratios to
  w(τ_mem), all classes: 1.53 / 2.67 / 11.95 / 26.52 / 58.38 at 30/100/300/500/1000 nm, identical to `ra5_delta.json`
  to machine precision. δ_kin = τ_z − τ_ℓ is internally consistent (max deviation 6e-14 ps).
* **Same witnesses, same meshes.** The 12 witnesses and the two meshes are the same for τ_mem and δ. Because the δ
  excursions are larger on 7×7×5, the robust rule selects the 9×9×6 value for every δ end (9×9×6-only reading = robust
  reading: 11.95 and 26.52). The 7×7×5-only reading gives 14.53 and 36.32.
* **Not mesh-converged.** Median excursion ratio 9×9×6/7×7×5 for τ_ℓ: 0.705 (300 nm) and 0.721 (500 nm); for τ_mem
  0.724 and 0.936. No τ_ℓ exists on 12×12×8. If the next refinement shrinks the δ excursions by the same factor as the
  last one (a pessimistic continuation, applied per witness), the ratios become 9.7× (300 nm) and 19.9× (500 nm). The
  absolute half-range of δ falls from 40.6 to 26.7 ps (500 nm) and from 15.9 to 10.8 ps (300 nm); the τ_mem
  half-width under the same continuation is 5.4 % and 5.7 %.
* **Inner bounds.** No witness was optimised for δ or at 500 nm. Both ranges are inner bounds, so their ratio is not a
  bound in either direction; "12–27× wider" is a property of these 12 tensors. For the resolvability direction the
  inner-bound nature is favourable (a wider true range is easier to resolve); the mesh trend is not.
* **Normalisation.** Absolute half-range = w(δ) |δ₀/τ₀|(9×9×6) τ_R,0(A1) = 40.6 ps at 500 nm. With the δ₀ of the FDTR base
  point (25.0 − 226.4 ps) it is 40.0 ps (1.6 %). At 30 nm the two normalisations differ by 60 % (0.42 vs 0.26 ps); this
  is irrelevant at the margins found there.
* **Normal/umklapp labels.** The kinetic τ_ℓ and δ use the full operator only; the label ambiguity (18–35 % of the rate
  in mixed orbits) affects only the grey-closure comparison (τ_N G1/G2/RQ).
* **Film conductivity co-varies.** The witnesses that set the δ ends move the film κ_z(L) (robust): cal3w min δ +23.1 %
  with κ_z −4.4 %; cal3w max δ −16.6 % with κ_z +3.2 %; orb 30 nm max −15.1 % / +2.3 %; class half-width of κ_z
  3.8–3.9 % at 300–500 nm. R-A6 holds λ_z at the A1 base for every witness ("the witnesses change only the dynamic
  part"). Under N1/N2 (λ_z known) this treats tensors whose film conductivities differ by up to 4 % as sharing one known
  λ_z. Under N3 the issue disappears.

## 4. FDTR consequence (sub-claim 5) — `au5_fisher1d.py`, `out/au5_fisher1d.json`

Own implementation: NumPy impedance recursion (no JAX, no campaign module). Thermal Jacobian by Richardson central
differences, step error 3.1e-9 (h vs h/4). Instrument columns are exact; nuisances are profiled by QR projection; the
same stack, frequencies, noise and nuisance sets as R-A6; base point A1 τ_R and λ_z with the mesh-corrected τ_ℓ.

| L | set | σ(δ) via τ_R / via −τ_ℓ (ps), 0.01° | Cattaneo σ(τ_R) | margin 0.01° / 0.1° (campaign) | margin, own | direct nonlinear test (static GK witnesses, Levenberg–Marquardt refit) 0.01° / 0.1° |
|---|---|---|---|---|---|---|
| 300 | N2 | 10.84 / 10.25 | 11.56 | 0.68 / 5.6 | 0.683 / 5.62 | 0.59 / 4.92 |
| 300 | N3 | 59.5 / 57.9 | 45.8 | 3.7 / 27 | 3.748 / 27.1 | 3.48 / 25.8 |
| 500 | N1 | 4.67 / 4.47 | 4.72 | — | 0.115 / 1.10 | 0.096 / 0.92 |
| 500 | N2 | 7.16 / 6.67 | 7.61 | 0.18 / 1.46 | 0.176 / 1.463 | 0.141 / 1.20 |
| 500 | N3 | 69.5 / 71.4 | 49.0 | 1.7 / 12.5 | 1.711 / 12.52 | 1.59 / 12.6 |

Kinetic/Cattaneo σ ratio: 0.94–0.99 (N1/N2) and 1.30–1.46 (N3), as stated in R-A6. The τ_R and −τ_ℓ directions agree
within 3–10 %. The extreme witness in the direct test is cal3w_m3L100min. The numbers of R-A6 are reproduced; the
issues are of scope:
* "Resolvable at 0.01°" holds for N1/N2. At 500 nm it also holds with λ_z free if C is known (σ = 20.5 ps, m = 0.51;
  300 nm: 1.49). With λ_z and C both free (N3), which review objection A-3 calls the realistic minimum, the campaign's
  own margins are 1.7 (500 nm) and 3.7 (300 nm), i.e. not resolvable. The δ signature is strongly correlated with λ_z
  and C (σ grows ×2.9 from N2 to λ_z free and ×9.7 to N3 at 500 nm).
* The N1/N2 premise is that every witness has the same known λ_z, while the δ-extreme witnesses differ in film κ_z by
  up to 4.4 % (§3).
* Range continuation (§3): the 500 nm N2 margin at 0.01° becomes 0.27, the 300 nm one 1.00.

## 5. Validity of the first-order regularised description

### 5.1 k² truncation (sub-claim 6) — `au2_rta.py`, `au5b_ktrunc.py`

The exact diffusive root of p|e|² + k²κ(k, p) = 0 is computed with κ to all orders in k (Newton with continuation from
low frequency) on the same regularised RTA operator, at the FDTR frequencies. λ_eff = −p|e|²/k² is compared with the
first-order GK form and with the first-order-in-k² form that is exact in p (the R-A1d quantity). The two meshes agree to
≤ 0.1 % at every frequency, so no isolated pole of the discrete mode set is being tracked.

| 500 nm, 64×64×35 | 10 MHz | 50 MHz | 100 MHz | 200 MHz |
|---|---|---|---|---|
| \|k\| L_GK | 0.13 | 0.29 | 0.41 | 0.58 |
| λ_eff/λ₀ exact | 0.9991+0.0150i | 0.9807+0.0680i | 0.9437+0.1108i | 0.8774+0.1475i |
| first-order GK | 1.0000+0.0151i | 1.0006+0.0756i | 1.0024+0.1511i | 1.0097+0.3020i |
| deviation of the non-Fourier part (λ_eff − 1) | 6 % | 28 % | 47 % | 67 % |
| exact vs first order in k², exact in p | 0.1 % | 2.4 % | 7.9 % | 21 % |
| δ_eff = −(λ_eff/λ₀ − 1)/(iω) (ps) | −239−14i | −217−61i | −176−90i | −117−98i |

At 300 nm the deviation of the non-Fourier part is 2.7 / 13 / 25 / 44 % at 10/50/100/200 MHz. Static kernel
coefficients: c₄/(K L⁴) = 4.26 (300 nm) and 4.81 (500 nm); c₆/(K L⁶) = 24.9 and 32.4. The GK resummation has 1 and 1;
the grey Callaway operator has 1–1.96. The growth reflects |k|Λ_max ≈ 2 at 200 MHz. R-A1d's "in-phase part of δ_eff within
2.3 % up to 200 MHz" kept the kernel at first order in k² and therefore did not test this.

FDTR consequence within the regularised bulk model (REF data from the exact kernel, 1D stack, nuisances refitted).

| L | stack | first-order GK film vs exact-kernel data, Δχ² N1 / N2 / N3 at 0.01° (0.1°) | δ returned by a GK fit (τ_ℓ free) vs static δ (ps), N2 0.01° | transfer of a witness-like change (τ_3ph ×(1±0.1)) N2 / N3, 0.01° | Δχ²(+10 %) exact / first order, N2 0.01° |
|---|---|---|---|---|---|
| 300 | note-16 | 7.3 / 4.3 / 0.07 (0.09 / 0.06 / 0.003) | −94.0 vs −115.0 | 0.69 / 0.88 | 0.12 / 0.26 |
| 300 | favourable | 132 / 54 / 0.26 | −82.3 vs −115.0 | 0.67 / 0.27 | 0.56 / 1.38 |
| 500 | note-16 | 102 / 70 / 6.6 (1.4 / 1.1 / 0.12) | −187.8 vs −240.6 | 0.68 / 0.36 | 0.93 / 1.98 |
| 500 | favourable | 2776 / 1748 / 7.3 | −166.7 vs −240.6 | 0.50 / ill-determined | 7.2 / 27.7 |

Reading. At 0.01° the first-order GK film is not an adequate model of the regularised kinetic medium at 500 nm. A
first-order analysis would report a δ 18–31 % smaller than the static one. A change of the operator that moves the static
δ by Δ moves the reported δ by 0.5–0.7 Δ (N1/N2) or 0.3–0.9 Δ (N3). The distinguishability Δχ² halves.

### 5.2 Matthiessen kernel vs a kinetic slab (sub-claim 7) — `au6_slab.py`, `au6_test_twostream.py`

**Model.** RTA with energy conservation (one local temperature), bulk three-phonon lifetimes of 64×64×35 and no
regulator. Walls are diffuse and black: T⁺(0) = T_A, T⁻(d) = T_B. Modes are binned in (log τ, log |v_z|), with v_z = 0
modes as local classes, preserving Σcτ and Σcv²τ per bin. The solution is exact for the binned set: modal expansion
S = Φ[a e^{−λz} + b e^{−λ(d−z)}] from the eigenpairs of diag(g²/v²) − diag(2g/(τv²))1βᵀ, with the diffusive eigenvalue
polished on the secular equation in cancellation-free form; there is no z-discretisation. The Matthiessen film
parameters (λ_M, τ_R,M, τ_ℓ,M, V1 completion) are computed on the same classes, so slab and GK film describe the same
material.

**Verification.** Single class: equal to the closed-form two-stream slab to 3e-14 at 1 kHz–10 GHz, static limit
G = 2cvΛ/(d+2Λ) to 8e-13. Two classes plus one local class: equal to an independent Chebyshev collocation of the
first-order system to ≤ 4e-5, inside the collocation's own N = 120/160 difference of 2–5e-4. Energy balance
q(0) − q(d) = pU to ≤ 1e-8 for the AlN classes. Large-d limit G d/λ_bulk = 0.70 / 0.93 / 0.99 at 1 / 10 / 100 µm.

**Static.** At 500 nm the black-wall slab has R = 3.37e-9 m²K/W against d/λ_M = 2.95e-9 (+14 %; λ_eff,slab = 148
against λ_M = 170 W/m K). At 300 nm: 2.36e-9 against 1.99e-9 (λ 127 against 151).

**FDTR test** (slab inserted in the 1D stack, R1, R2 on the wall-emission temperatures; fitted with the twin's GK film,
local convention; nb = number of bins per axis):

| | nb = 12 | nb = 20 | nb = 30 | nb = 45 |
|---|---|---|---|---|
| classes | 95 | 234 | 468 | 920 |
| 500 nm δ_M (ps) | −228.2 | −234.7 | −238.0 | −239.5 (unbinned −240.6) |
| 500 nm δ fitted, λ_z free (N2lam) / N3, 0.01°, note-16 | — / −70.5 | −75.7 / −73.6 | −77.3 / −75.1 | −78.1 / −75.9 |
| ratio δ_fit/δ_M (N3) | 0.309 | 0.314 | 0.316 | 0.317 |
| transfer (τ_3ph ×(1±0.1)), N3 | 0.38 / 0.37 | 0.37 / 0.35 | 0.36 / 0.34 | (not run) |
| 300 nm δ_M; δ fitted N3 (ratio) | −108.4; −32.8 (0.30) | −111.8; −35.1 (0.31) | −113.5; −36.2 (0.32) | (not run) |

The N3 δ_fit sequence −70.5 → −73.6 → −75.1 → −75.9 ps has differences 3.1, 1.5, 0.8 (geometric limit ≈ −76.7 ps,
ratio 0.319 against the unbinned δ_M). Binning error ±0.002 in the ratio: CONVERGED.

Further results at nb = 30 (500 nm, note-16, 0.01°):
* With λ_z free and C fixed (N2lam), the GK film with the regularised τ_ℓ,M is rejected (Δχ² = 48.2; N3 4.4). With τ_ℓ
  free the fit is perfect (Δχ² = 0.005) at λ_z = 210.5 W/m K, i.e. 1.24 λ_M, with the wall resistance absorbed into R1,
  R2. A Fourier film gives Δχ² = 12.5 (N2lam) / 1.67 (N3). A Cattaneo film with τ_R free fits as well (Δχ² 0.03) with
  τ_R = −80 ps. The slab shows one first-order signature, δ ≈ −75 ps, in line with the first-order δ-only dependence of
  Result B.
* With λ_z fixed at λ_M (N2) the slab cannot be fitted (Δχ² = 5171 with τ_ℓ,M; 367 with τ_ℓ free). Fixing λ_z at the
  slab's own static value is worse (13104). In a film that behaves like the slab, the N1/N2 premise "λ_z known" has no
  value of λ_z that makes the GK film consistent.
* Port convention: with R1, R2 on the energy temperatures at the walls, δ_fit = −74.0 (N2lam) / −72.8 (N3) ps. Stack:
  favourable resistances give −75.7 / −61.6 ps. The ratio 0.26–0.33 holds in every case at 0.01°.
* 300 nm: δ_fit = −36.8 (N2lam) / −36.2 (N3) ps against −113.5 ps; GK with τ_ℓ,M rejected by Δχ² = 6.6 (N2lam).

Reading. The |v_z|/L bulk rate gives the boundary-limited carriers, which carry 76 % of τ_ℓ at 500 nm (R-A6), a bulk
nonlocal response with Λ_z ~ L. In a slab these carriers fly wall to wall: their contribution appears as wall
resistance and as a larger interior conductivity, not as a gradient correction. Within RTA with diffuse black walls, the
regularised kernel overstates the FDTR-visible first-order time by a factor ≈ 3.1 (300–500 nm). It overstates its
sensitivity to changes of the lifetimes by a factor ≈ 2.8 for a uniform change and 3.2–3.8 for the witnesses' own
lifetime changes (§5.4). Its static conductivity is off by 14–24 % depending on what
it is compared with. The slab is itself a model (RTA, black walls, no normal-process structure, no interface
transmission), but it treats the walls explicitly, which the regulator does not. The two descriptions of the same
material and the same film disagree by a factor ≈ 3 in exactly the quantity the claim is about, so the regularised value
of τ_ℓ is not a reliable estimate of the film's visible nonlocal time.

Also at nb = 30, the fitted wall resistances rise by 4.5–5.0e-10 m²K/W each, symmetric at both walls and independent
of the stack. Their sum equals R_slab − d/λ_fit to four digits (500 nm: 9.955e-10 against 9.960e-10). The slab is two
Knudsen-layer resistances plus an interior film of conductivity 206–210 W/m K (bulk κ_z 252, λ_M 170) whose first-order
time is δ ≈ −36 / −76 ps.

### 5.3 Full operator with the exact kernel, 7×7×5 witnesses — `au3b_exact_root775.py`, `au3c_spline_check.py`

Exact diffusive root of the independently assembled 7×7×5 operator (QDQ completion, dense complex-symmetric solves,
secant in k², 12 frequencies 1–200 MHz, spline in ln f to the 80 FDTR frequencies) for REF and the four witnesses that
set the δ ends (cal3w max/min, orb 0.1 30 nm max, cal1nz min). The static times reproduce the campaign's 7×7×5 values
(REF 300 nm 19.838 / 104.75 ps).
* Deviation of the non-Fourier part of λ_eff at 200 MHz from the first-order GK form: 0.25–0.39 (300 nm), 0.40–0.62
  (500 nm; REF 0.48). It is largest for the long-τ_ℓ witnesses and smaller than on the converged RTA mesh (0.44 / 0.67),
  because 7×7×5 lacks long-MFP carriers.
* First-order film vs exact-kernel REF data, Δχ² N1 / N2 / N3 at 0.01°: 2.1 / 1.2 / 0.03 (300 nm, note-16), 22.8 / 14.1 /
  0.86 (500 nm, note-16), 862 / 500 / 0.56 (500 nm, favourable).
* Witness distinguishability, margin ratio m_exact/m_first-order (REF vs witness, both through the same kernel type):

| | N2 0.01° / 0.1° | N3 0.01° / 0.1° |
|---|---|---|
| 300 nm, note-16 | 1.26–1.48 / 1.24–1.50 | 1.15–1.21 / 1.14–1.28 |
| 300 nm, favourable | 1.65–1.82 / 1.70–2.12 | 0.70–1.47 / 0.16–2.03 |
| 500 nm, note-16 | 1.30–1.45 / 1.28–1.48 | 0.67–0.94 / 0.57–0.98 |
| 500 nm, favourable | 1.56–1.70 / 1.58–1.79 | 0.35–1.33 / 0.06–0.23 |

  Under N2 the exact kernel makes the witnesses less distinguishable, consistent with the RTA transfer 0.5–0.7 of §5.1.
  Under N3 at 500 nm it can make them more distinguishable. The exact witness differences change the real part of λ_eff
  above ≈ 50 MHz (amplitude channel), which is less degenerate with λ_z and C than the first-order phase signature.
  Check (AU3c): the spline has leave-one-out errors ≤ 0.2 % at interior nodes, and cubic and PCHIP interpolation give the
  same Δχ² to 0.2 %. On the favourable stack the N3 Δχ² of cal1nz min (1.71) is 98 % amplitude (max |Δ ln|H|| =
  0.6 %); on the note-16 stack it is phase-dominated. Within the regularised model this effect is real, but it lives at
  |k|Λ_max ≈ 2, where §5.2 shows that the bulk kernel does not describe the film.
* Reading: the first-order margins are not conservative or optimistic in a fixed direction. They are wrong by factors
  0.35–1.8 at 0.01° (0.06–2.1 including 0.1°) depending on the nuisance set and stack, already within the regularised
  model.

### 5.4 Witness-specific transfer to the slab — `au8_witness_slab.py`

The slab test of §5.2 used a uniform lifetime scaling as a proxy for a witness. Here each witness's own single-mode
lifetimes enter both films. These are the RTA diagonal of its orbit-tied operator; orbits are labelled by the minimum
event index over the 24 images, 0 missing images; the REF diagonal reproduces `ra1c_ref_rta.json` exactly on both
meshes. Slab and Matthiessen film are evaluated on mode classes with membership fixed by REF (no binning noise in
differences), and each slab is fitted by the REF GK film with τ_ℓ free.

| mesh | L | δ-range transfer slab/Matthiessen, N2lam / N3 | witness T_W (N3), all | range-setting witnesses (cal3w min, cal3w max, orb 30 max, cal1nz min) | RTA-diagonal vs full-operator δ excursions: corr / slope |
|---|---|---|---|---|---|
| 7×7×5 | 300 | 0.301 / 0.300 | 0.12–0.60 | 0.12, 0.24, 0.30, 0.16 | 0.93 / 0.36 |
| 7×7×5 | 500 | 0.314 / 0.315 | 0.16–1.8 ¹ | 0.18, 0.24, 0.32, 0.16 | 0.93 / 0.33 |
| 9×9×6 | 300 | 0.281 / 0.281 | 0.14–0.34 | 0.14, 0.19, 0.29, 0.15 | 0.98 / 0.64 |
| 9×9×6 | 500 | 0.263 / 0.264 | 0.15–0.49 | 0.17, 0.21, 0.28, 0.17 | 0.98 / 0.62 |

¹ orb 0.1 300 nm optima, whose δ_M excursion is near zero (ill-conditioned ratio; irrelevant for the range ends).

Reading. The witnesses that set the δ range transfer least: 0.14–0.21 for the calibrated 100 nm optima. Their δ range
shrinks to 0.26–0.31 of its regularised value. On 9×9×6 the RTA diagonal carries 62–64 % of the full-operator δ
excursions, with correlation 0.98. The off-diagonal remainder is not covered by this test.

## 6. Combined margins — `au7_summary.py`, `au7b_table.py`, `au9_lambda_prior.py`

m = σ(δ)/(visible half-range of δ), 1D, local convention; m ≤ 1 = resolvable at 1σ. Readings: A = campaign
(first-order regularised GK film, campaign half-range); B = A with the geometric mesh continuation of the range
(×1/0.657 at 500 nm, ×1/0.681 at 300 nm); C = B with the k²-truncation transfer (§5.1; N1/N2); D = kinetic slab (σ at
the slab-fitted base point; visible half-range = campaign half-range × witness-specific range transfer of §5.4,
9×9×6); E = D with the mesh continuation. In the slab only N2lam (λ_z free, C known) and N3 are meaningful; with λ_z
fixed at λ_M the slab cannot be fitted.

| L | stack | σ_φ | A: N1 / N2 / N2lam / N3 | B: N2 / N2lam / N3 | C: N1 / N2 | D: N2lam / N3 | E: N2lam / N3 |
|---|---|---|---|---|---|---|---|
| 500 | note-16 | 0.01° | 0.11 / 0.18 / 0.51 / 1.71 | 0.27 / 0.77 / 2.60 | 0.24 / 0.39 | 1.82 / 5.38 | 2.77 / 8.18 |
| 500 | note-16 | 0.1° | 1.10 / 1.46 / 3.24 / 12.5 | 2.23 / 4.93 / 19.0 | 2.44 / 3.62 | 11.4 / 38.9 | 17.4 / 59.2 |
| 500 | favourable | 0.01° | 0.04 / 0.05 / 0.25 / 2.59 | 0.08 / 0.38 / 3.93 | 0.12 / 0.16 | 0.87 / 8.42 | 1.32 / 12.8 |
| 500 | favourable | 0.1° | 0.40 / 0.49 / 2.21 / 15.1 | 0.74 / 3.35 / 22.9 | 1.32 / 1.80 | 7.46 / 48.9 | 11.3 / 74.3 |
| 300 | note-16 | 0.01° | 0.49 / 0.68 / 1.49 / 3.75 | 1.00 / 2.18 / 5.50 | 0.99 / 1.45 | 5.20 / 11.7 | 7.64 / 17.2 |
| 300 | note-16 | 0.1° | 4.59 / 5.62 / 9.93 / 27.1 | 8.26 / 14.6 / 39.9 | 9.79 / 12.5 | 34.4 / 82.1 | 50.5 / 121 |
| 300 | favourable | 0.01° | 0.25 / 0.33 / 1.59 / 8.33 | 0.48 / 2.34 / 12.2 | 0.54 / 0.72 | 5.55 / 25.2 | 8.14 / 37.1 |
| 300 | favourable | 0.1° | 2.09 / 2.69 / 10.5 / 40.2 | 3.95 / 15.4 / 59.0 | 6.59 / 9.61 | 34.6 / 114 | 50.8 / 168 |

**Required knowledge of λ_z** (AU9; N2 plus a Gaussian prior of relative width s on λ_z, C known). At 500 nm, 0.01°,
note-16, the regularised range is resolved for every s ≤ 30 % (m = 0.36 at s = 1 %). The slab-visible range (10.7 ps)
needs s ≤ 0.48 % (m = 1.32 at 1 %, 1.78 at 5 %). On the favourable stack the slab-visible range is resolved even with
λ_z free (m = 0.87). At 300 nm the slab-visible range is not resolved even with λ_z exact (m = 2.5 note-16, 1.24
favourable). At 0.1° the slab-visible range is not resolved in any case, even with λ_z exact (m ≥ 1.9).

Reading. "Resolvable at 0.01°" survives in readings A–C at 500 nm. In the kinetic-slab readings D–E it survives only
for low interface resistances (favourable stack, λ_z free: 0.87, 1.32 with the mesh continuation) or for λ_z known to
better than 0.5 %. At 300 nm it holds only with λ_z and C known (N1/N2) in readings A–C, and on the note-16 stack it is
marginal there (B 1.00, C 1.45). At 0.1° it fails in every reading on the note-16 stack.

## 7. Answers to the requested points

1. **Re-derivation and convention.** The expansion is correct. It is the exact e⊥ reduction followed by a Neumann
   series, with the odd term removed by parity. τ_ℓ = |e|²(u·Ĉ⁻¹u)/K² is the twin's τ_ℓ = (1+α)ℓ²C/λ = L²/a because the
   twin's 1D GK film has the kernel λ/(1 + τ_Rp + (1+α)ℓ²k²) and T is the energy temperature. The energy projection is
   required and sufficient: u = QVx, no momentum projection. The unprojected form relaxes the conserved energy mode and
   gives 9τ_c/5 in the grey case; the projected form gives 4τ_c/5. Both are reproduced exactly with independent code.
   One point is missing from the record: the value depends on the even-sector completion of the regulator.
2. **Independent recomputation.** RTA τ_ℓ on 48×48×26 and 64×64×35 is reproduced by three routes to machine
   precision (§2). The 7×7×5 event operator is reproduced by an independent assembly to ≤ 4e-15 for all 13 tensors and
   5 lengths.
3. **δ range.** It uses consistent witnesses and meshes (robust rule = 9×9×6 for δ). The normalisation is consistent
   within 2 % at 300–500 nm. N/U labels play no role in the kinetic closure. It is an inner-bound ratio, not a bound,
   and it is not mesh-converged (excursions shrink ×0.70–0.72 per refinement; 9.7×/19.9× under a geometric
   continuation). The witnesses' film κ_z co-varies by up to 4.4 %.
4. **FDTR consequence.** It is reproduced by an independent 1D Fisher implementation to three digits (σ(δ), margins,
   direct nonlinear test). It holds within the first-order regularised GK film at 500 nm with C known (also with λ_z
   free: 0.51), but not with C free (N3: 1.71). In the kinetic-slab reading, which uses the witnesses' own lifetimes, the
   visible range is 0.26× and the 500 nm note-16 margin becomes 1.82 with λ_z free. Resolution would need λ_z known to
   ≤ 0.5 % or interface resistances near 2e-9 m²K/W. At 300 nm and at 0.1° no slab reading resolves (§6).
5. **Validity.** The k² truncation fails in the upper half of the band at 300–500 nm (non-Fourier part off by 28–67 % at
   50–200 MHz; first-order model misfit Δχ² = 70 at 0.01°). The Matthiessen-kernel film model overstates the visible
   nonlocal time by a factor ≈ 3 compared with a kinetic slab of the same material.

## 8. Suggested record changes

* R-A6 §2 and `10-proof-status.md`: keep "τ_ℓ formula DERIVED", add "(independently re-derived and checked; the value
  depends on the energy-conserving completion of the regulator on the even sector: QDQ vs relaxation to the wall
  temperature +2.5 / +4 / +12 / +38 % at 500 / 300 / 100 / 30 nm)".
* "τ_ℓ(L) values NUMERICALLY SUPPORTED (RTA converged < 0.5 %)": RTA route CONVERGED (reproduced); the quoted
  full-operator ratios 1.81/2.97/6.30/9.20 are 9×9×6 values, UNDER-RESOLVED (+13 % / +22 % from 7×7×5 to 9×9×6 at
  300 / 500 nm), mesh-corrected values ±10 %.
* R-A6 §5, "the static pair (τ_mem, τ_ℓ) is an adequate first-order description of the regularised kernel in the FDTR
  band": CONTRADICTED above ≈ 20 MHz at 300–500 nm (§5.1). R-A1d tested the p-dependence only.
* R-A6 §4 / synthesis, "static-data range of δ 12–27× wider": restate as "12–27× on the twelve non-optimised witnesses,
  9×9×6 values, not mesh-converged (≈ 10–20× under a geometric continuation of the mesh trend), not a bound;
  within the regularised first-order kinetic closure".
* "Would be resolvable at 0.01° (1D, local)": not supported as stated. Replace by: "within the first-order regularised
  GK film the margin at 0.01° is 0.18 / 0.51 / 1.71 (500 nm; N2 / λ_z free with C known / N3) and 0.68 / 1.49 / 3.75
  (300 nm). The first-order film misdescribes the regularised kernel itself at 0.01° (Δχ² = 70, 500 nm, N2). In a
  kinetic slab with diffuse black walls (RTA) the visible δ is ≈ 0.32 of the regularised value and the witnesses' δ range
  is ≈ 0.26–0.31 of it. The 500 nm margin with λ_z free then becomes 1.8 (note-16) and 0.87 (R1 = R2 = 2e-9), and
  resolution on the note-16 stack would need λ_z known to ≤ 0.5 %. The N1/N2 premise with the regularised λ_z has no
  consistent meaning in the slab." Status: UNRESOLVED (model-dependent at the factor-3 level; resolvable only for
  favourable interfaces or sub-percent knowledge of λ_z, at 500 nm and 0.01°).
* "The Cattaneo comparison is conservative: CONTRADICTED under the kinetic closure": the qualitative point survives
  (in every model tested the first-order visible time of a 300–500 nm film is not τ_mem and is larger in magnitude). The
  quantitative support (12–27×, resolvable at 0.01°) does not. Suggested status: SUPPORTED qualitatively, magnitude
  UNRESOLVED.
* Required qualifier of the headline: replace "τ_ℓ ≈ 6–9 τ_mem … would be resolvable at 0.01°" by "the first-order
  FDTR-visible time of the same films is not τ_mem; its value and static-data range depend on the treatment of the
  boundary-limited carriers at the factor-3 level (regularised kernel δ ≈ −240 ps, kinetic slab δ ≈ −75 ps at 500 nm,
  RTA), so the Cattaneo-projection verdict does not extend to the film's full first-order response, and no
  resolvability statement for that response is established".
* Open items that would settle it: a slab (or Monte Carlo) solution with the full operator and realistic interface
  transmission; δ-optimised witnesses; τ_ℓ on 12×12×8 (RTA diagonal light, full operator not).

## 9. Files

Scripts (`experiments/R_A_audit/`), outputs (`experiments/R_A_audit/out/`):
* `au1_derivation.py` → `au1_derivation.json`, `au1.log`: exact grey case (SymPy), random non-grey operators (bordered
  kernel, KKT closed form, dispersion root), twin convention.
* `au2_rta.py` → `au2_rta.json`, `au2.log`: RTA τ_ℓ on 48×48×26 and 64×64×35, completions V1–V4, Cauchy-integral
  kernel coefficients (c₄, c₆), exact diffusive root vs first-order forms.
* `au3_event775.py` → `au3_event775.json`, `au3.log`: independent 7×7×5 assembly, orbit tying, KKT solves, QDQ vs
  wall-temperature completion, witness excursions.
* `au3b_exact_root775.py` → `au3b_exact_root775.json`, `au3b.log`: exact diffusive root of the full 7×7×5 operator for
  REF and four range-setting witnesses; FDTR distinguishability exact vs first order.
* `au3c_spline_check.py` → `au3c_spline_check.json`: interpolation check and amplitude/phase split of the AU3b N3 result.
* `au4_delta_range.py` → `au4_delta_range.json`, `au4.log`: δ-range reconstruction, mesh trend, alternative readings,
  κ_z co-variation.
* `au5_fisher1d.py` → `au5_fisher1d.json`, `au5.log`: independent 1D Fisher, margins, direct nonlinear test.
* `au5b_ktrunc.py` → `au5b_ktrunc.json`, `au5b.log`: FDTR consequence of the k² truncation, transfer test.
* `au6_slab.py`, `au6_test_twostream.py` → `au6_slab_{12,20,30}.json`, `au6_slab_45_500_eta0.json`,
  `au6_test_twostream.json`, logs: kinetic slab, verification, FDTR fits.
* `au8_witness_slab.py` → `au8_witness_slab.json`, `au8.log`: witness-specific slab transfer (RTA diagonal, both
  event meshes).
* `au7_summary.py`, `au7b_table.py` → `au7_summary.json`, `au7b_table.json`: combined margins (table of §6).
* `au9_lambda_prior.py` → `au9_lambda_prior.json`, `au9.log`: required prior precision on λ_z.

Run order: au1, au2, au3, au4, au5, au5b, au6 (12, 20, 30; then `45 500 eta0`), au6_test_twostream, au3b, au3c, au8,
au7, au7b, au9. Each runs in ≤ 40 min on two threads; peak working set below 1 GB (dense 7×7×5 operators ≈ 0.4 GB).
