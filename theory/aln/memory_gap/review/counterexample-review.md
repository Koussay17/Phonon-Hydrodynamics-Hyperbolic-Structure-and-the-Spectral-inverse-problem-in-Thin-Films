# Counterexample review — identifiability gap of the memory moment

Campaign `20261004-memory-gap`, 4 October 2026. Subject: `13-final-report.md` (with `02-assumptions.md`,
`07-experiments.md`, `scripts/`, `results/`). Headline under attack: for linear energy-conserving phonon kinetics in
entropy coordinates, C(g) = Σ g_α a_α a_αᵀ (g ≥ 0, fixed event set), fixing the lifetimes diag C = r, the energy
invariant and the DC response K = bᵀC⁺b leaves τ_mem = bᵀ(C⁺)²b / bᵀC⁺b pinned from below (τ ≥ τ_CS = K/|b|²,
nearly attained) but unidentified from above. AlN 5×5×3, 300 K, σ = 0.1 THz, 44,406 events, 6mm×TR-symmetric rates:
basal 19.197 ps to ≥ 6906.7 ps (reference 34.25 ps), c axis 15.517 ps to ≥ 959,689.9 ps (reference 27.45 ps).
Also L3: K ≥ K_RTA/s_max.

Independence. Every attack below was designed and run without consulting `review/proof-audit.md` (not read).
All operators, orbits, moments, optimisers and certifications are new code (`review/cx_*.py`) written from the
definitions in `02-assumptions.md`. The campaign code was used only (i) to decode the orbit numbering of the stored
witnesses (`cx_common.their_ev_orb`; every decoded witness is re-checked for orbit constancy against independently
computed orbits, max deviation 7e-16), (ii) unchanged (`scripts/aln_events.py`) to export a 7×7×5 event geometry
into `review/cx_runs/`, (iii) as a copy with one flag changed (`cx_aln_events_symfc.py`). Python: research
environment (numpy 2.5.3, scipy 1.18.1, python-flint 0.9.0, hypothesis 6.168). Raw outputs: `review/cx_runs/`,
`review/cx_*.json`; table of all optimisation runs: `review/cx_summary.json` (`cx_summarize.py`).

---------------------------------------------------------------------------------------------------------------

## 0. Verdict in brief

* **Inside its declared class the claim survives every attack.** The class is arbitrary non-negative rates on
  the fixed event set, no amplitude model, and data = 300 K lifetimes plus the DC tensor. The witnesses are
  genuine, not numerical artifacts. They reproduce by two independent solution routes and in 256-bit ball
  arithmetic with exact energy conservation. A 1e-10 perturbation of the harmonic data does not move them. The
  lower end is pinned by the proved τ_CS.
* **The headline AlN witnesses are not physically admissible.** All of them switch on events whose three-phonon
  vertex is numerically zero, with the statistical signature of crystal-symmetry selection rules (§3.2). The 6906.7 ps basal witness takes up to 99.998 %
  of some mode lifetimes from such events. This breaks the specific witnesses, not the conclusion:
  selection-rule-respecting witnesses, certified here, reach τ_x ≥ 7038.17 ps and τ_z ≥ 4,322,183 ps.
* **The gap exists only without an amplitude prior** (scope qualification, quantitative). Rigorous ceiling: if
  no event rate is suppressed below 1/F of a reference operator's rate, then τ ≤ F/λ_min(C_ref|H) =
  **81.50·F ps** (5×5×3). With a factor-2 prior the whole admissible interval lies in [τ_CS, 163.0] ps. The
  reported witnesses require some events suppressed by factors ≥ 85 (basal) and ≥ 1.2e4 (c axis). A smooth
  amplitude prior (polynomial in the three energies, degree ≤ 20) gave no witness above 1.8× the reference.
* **Additional data.**
  * Independently computed observables discriminate the witnesses strongly: GHz-range K(z), mode accumulation,
    and κ(T) with temperature-independent |Φ|².
  * Fixing the κ(T) curve does not identify the memory either: verified c-axis witnesses at ≥ 1.3e5 ps match
    κ_xx and κ_zz at five temperatures.
  * Nor does the accumulation function (τ_z ≥ 61,254 ps).
  * Fixing K(z) at three Laplace rates (1/ns to 1/(10 ps)) cut the best maxima found to 62.9 ps (x) and
    33.9 ps (z). This is inner evidence only; a slow mode of small weight remains in principle possible.
  * Multi-temperature lifetimes could not be tested conclusively (numerical stiffness).
* **L3 survives as stated** (r = diag C). It is false if r is read as the self-energy linewidth of RTA codes
  when repeated-daughter events are present (exact counterexample, §6).

---------------------------------------------------------------------------------------------------------------

## 1. Summary table

| # | Attack | Verdict | Decisive numbers | Scripts / data |
|---|---|---|---|---|
| 1a | Positivity, detailed balance, time reversal, momentum selection, orbit symmetry of witnesses | SURVIVES | g ≥ 0; C symmetric in entropy coordinates; orbit-constant to 7e-16 | `cx_verify_reference.py/.json` |
| 1b | Crystal-symmetry selection rules (vertex zeros) | **BROKEN** for the headline witnesses; conclusion SURVIVES with replacements | 3,275 events (7.4 %) with rates ≤ 4.4e-15 of max; witness 6906.7 ps takes up to 99.998 % of lifetimes from them; certified selection-rule witnesses 7038.17 ps (x), 4,322,183 ps (z) | `cx_forbidden_anatomy`, `cx_selection_rules`, `cx_verify_newwitness`, `cx_certify_new` |
| 1c | Factor-F prior on each orbit rate | **SCOPE QUALIFICATION** (rigorous) | τ ≤ 81.50 F ps; found max 45.6 / 153.5 / 480.1 ps (x) and 44.2 / 162.6 / 1421 ps (z) for F = 2 / 10 / 100 | `cx_attack1_box.py`, `cx_jobs.py`, `cx_verify_box.json` |
| 1d | Smooth amplitude prior (polynomial in energies) | **SCOPE QUALIFICATION** (inner evidence) | degree ≤ 10: reference isolated; degree 12–20: max found ≤ 1.80× reference | `cx_attack1_smooth*.py` |
| 1e | Nature of the hidden slow modes / sum rules | SURVIVES (no sum rule violated); physical reading weakened | slow modes are not quasi-momenta (overlap ≤ 0.30); 3.58e7 ps witness = weight 1.1e-6 at 32 s | `cx_attack1_slowmode` |
| 1f | Non-symmetric ("without symmetry") witnesses | **BROKEN** as physical statements | τ_x = 114,019 ps but τ_y = 514 ps in a hexagonal crystal | `cx_forbidden_unconstrained.json` |
| 2 | Reference operator vs phono3py Ω′; mesh | SURVIVES; SCOPE (unconverged mesh) | C⁺b differs by 0.15 % (x), 0.79 % (z); 7×7×5: τ_ref 52.6 / 50.6 ps | `cx_attack2_reference`, `cx_attack2_mesh775` |
| 3 | Numerical artifacts | SURVIVES | arb, exact conservation: 35,759,344.509 ps, unchanged; 1e-10 data perturbation: ±1.1 % | `cx_attack3_arb`, `cx_certify_new` |
| 4 | Additional measurable data | κ(T) curve and accumulation: gap SURVIVES. K(z) at 3 Laplace rates: found maxima collapse (inner evidence only). Lifetimes at 5 T: INCONCLUSIVE (optimiser stiffness) | κ(T)-matched τ_z ≥ 131,918 ps; accumulation-matched τ_z ≥ 61,254 ps; K(z)-matched max found 62.9 / 33.9 ps; 5-T lifetimes: ≤ 35.9 / 28.4 ps found | `cx_attack4_observables`, `cx_verify_multiT`, `cx_verify_extras`, `cx_runs/mK5_*`, `Kz_*`, `acc_*`, `mT5*` |
| 5 | L3 counterexamples (4-phonon, repeated indices, zero modes) | SURVIVES; SCOPE (r must be diag C) | min K·s_max/K_RTA = 1.0000 over all searches; exact chain with self-energy linewidth: 0.3333 < 1/2 | `cx_attack5_*` |

---------------------------------------------------------------------------------------------------------------

## 2. Independent reconstruction (basis for all attacks)

`cx_common.py` rebuilds the event operator from `results/aln/events_m553_s0.1.npz`:
n = 897, m = 44,406, 2,627 event orbits, 117 mode orbits (my orbit computation, 24 operations), umklapp 56.1 %.
Energy conservation of the event vectors max |a·e|/(|a||e|) = 2.1e-17; b·e/(|b||e|) ≤ 6e-17.

Reference (orbit-averaged rates), reproduced:

| | κ (W/m K) | τ_mem (ps) | τ_CS (ps) | τ_RTA (ps) | K/K_RTA |
|---|---|---|---|---|---|
| x | 214.2999 | 34.2530 | 19.1903 | 30.1186 | 1.2572 |
| z | 201.9670 | 27.4465 | 15.5168 | 22.5213 | 1.4319 |

Smallest eigenvalue of C_ref on H = e⊥: λ_min = 9.7643e-4 THz (1/(4πλ_min) = **81.4986 ps**), condition number 148.5.

---------------------------------------------------------------------------------------------------------------

## 3. Attack 1 — physicality of the witnesses (most important)

### 3.1 Elementary physical constraints — SURVIVES
Every witness has g ≥ 0. Detailed balance and time reversal are built into the symmetric entropy-coordinate
form. The event set is fixed by momentum and energy selection, so normal events conserve crystal momentum by
construction. The symmetric witnesses are constant on the 6mm × TR orbits computed here (max relative deviation
7e-16). None of these is violated.

### 3.2 Selection-rule zeros — headline witnesses BROKEN, conclusion survives
**Finding.** The reference rates fall into three populations: 40,527 events between 1e-9 and 1 of the maximum,
and two clusters near 1e-20 and 1e-31. The raw phono3py rates show a clean gap. The largest "forbidden" rate is
4.35e-15 of the maximum; the smallest allowed one is 3.07e-9. Neither the Gaussian weight (≥ e⁻⁸ inside 4σ) nor
the Bose factors can produce a 1e-15 suppression, so these 3,275 events (7.38 %, 335 orbits) have a numerically
zero vertex |Φ|².

**Evidence that these are space-group selection rules** (`cx_selection_rules.py`):
* I reconstructed the 24 q-maps as integer matrices and kept the 12 unitary C6v operations (q_z preserved).
* 95.6 % of the forbidden events have a nontrivial common unitary stabilizer of all three wave vectors (a
  mirror or glide, or the full C6v on the Γ–A axis). Only 17.2 % of allowed events do.
* Among the 10,206 events that share such a little group, 30.7 % vanish. This is the parity signature of
  selection rules.
* Not proved event by event. The decisive check would be irreps of the eigenvectors, or a random symmetric
  fc3. The cheap test of re-exporting with `symmetrize_fc=True` was uninformative: it gives bit-identical
  rates, because the flag does not act on read-in force constants.

The 6mm×TR orbit tying of the campaign's "symmetric" class does not impose these zeros, so that class still
admits forbidden processes.

**Reliance of the witnesses on forbidden events** (`cx_forbidden_anatomy.json`, `cx_forbidden_unconstrained.json`):

| witness | τ (ps) | active forbidden events | max share of a mode's lifetime from forbidden events | modes > 10 % | effect of zeroing them |
|---|---|---|---|---|---|
| full_c0_max_08 (headline basal) | 6906.7 | 178 | 0.99998 | 188 | lifetimes off by up to 100 %, K_zz +833 % |
| full_c0_max_09 | 7163.6 | 145 | 1.000 | 211 | C becomes singular |
| full_c0_max_04 (near-singular) | 3.576e7 | 6 | 0.155 | 6 | lifetimes off 15.5 % |
| full_c2_max_01 (headline c axis) | 959,689.9 | 24 | 0.329 | 6 | lifetimes off 33 %, K_zz +21 % |
| none_c0_max (no symmetry) | 114,019 | — | 0.816 | 12 | — |
| none_c2_max (no symmetry) | 120,195 | — | 0.536 | 3 | — |
| x-min / z-min | 19.197 / 15.517 | 24 / 86 | 0.12 / 0.55 | 12 / 84 | — |

**Repair test.** I reran the maximisation with forbidden rates fixed at exactly 0: same data, symmetric rates,
log-rates in [−60, 60], my gradient-projection optimiser `cx_opt.py`. Results (`cx_runs/nofor_*`), verified
independently (`cx_verify_newwitness.json`) and certified in 256-bit arb with exact energy conservation
(`cx_certify_new.json`):

* basal: **τ_x = 7038.173225 ps** (arb radius ≈ 1e-66; K_xx residual 3.4e-14; lifetimes 5.8e-13; λ_min 2.8e-6,
  cond 4.2e4), from the campaign witness pattern. Other starts give 5076.7, 4436.8 and, from the reference
  itself, 2526.3 ps.
* c axis: **τ_z = 4,322,183.474 ps** (arb radius < 1e-60; K_zz residual 2.9e-11; λ_min 1.1e-8, cond 1.1e7),
  reached from the plain reference start. A start from the campaign witness pattern gives 23,456.2 ps.

So the qualitative conclusion survives the selection rules. In the restricted class the c-axis lower bound is
even larger than claimed (1.57e5 × reference, 2.79e5 × τ_CS). The specific witnesses quoted in 13-final-report
must be replaced.

### 3.3 Factor-F prior — SCOPE QUALIFICATION (rigorous ceiling plus inner values)
**Proposition (proved here; two lines).** Let g_ref be admissible. If g_α ≥ g_ref,α/F for every event
(one-sided; nothing else is assumed), then:
* C(g) − C(g_ref)/F = Σ (g_α − g_ref,α/F) a_α a_αᵀ ⪰ 0, so λ_min(C(g)|H) ≥ λ_min(C(g_ref)|H)/F;
* since x = C⁺b ∈ H, τ(g) = |x|²/(xᵀC x) ≤ 1/λ_min(C(g)|H).

Hence

  **τ(g) ≤ F / λ_min(C_ref|H) = 81.50·F ps** (AlN 5×5×3; 197.2·F ps at 7×7×5, §4),

and, with τ ≥ τ_CS, τ_max/τ_min ≤ 4.247·F (x) and 5.252·F (z). This holds for symmetric and non-symmetric
rates, with or without the lifetime and K data. Conversely, a witness with value τ must suppress at least one
event by a factor ≥ τ/81.50 ps. The values are 84.7 (6906.7 ps), 86.4 (7038.2 ps), 11,776 (959,689.9 ps),
53,034 (4,322,183 ps) and 4.4e5 (3.58e7 ps). The witnesses actually suppress 49–87 % of all events below
1e-6 of the reference (`cx_verify_reference.json`).

**Inner values under the box g_ref/F ≤ g ≤ F·g_ref** (symmetric rates, all lifetimes, K_xx and K_zz fixed):
* method: gradient projection with active set and exact Newton restoration, from the reference, from random
  points (these failed to restore except at F = 2), and by homotopy from the selection-rule witnesses;
* every listed witness was re-verified: inside the box, data ≤ 2e-12, two solution routes agree
  (`cx_verify_box.json`).

| F | τ_x found [min, max] (ps) | τ_z found [min, max] (ps) | rigorous ceiling (ps) | found max/τ_CS (x, z) | max/τ_CS allowed by ceiling (x, z) |
|---|---|---|---|---|---|
| 1 | 34.25 | 27.45 | — | 1.78, 1.77 | — |
| 2 | [27.96, 45.63] | [20.71, 44.16] | 163.0 | 2.38, 2.85 | 8.5, 10.5 |
| 10 | [20.73, 153.48] | [15.72, 162.60] | 815.0 | 8.0, 10.5 | 42.5, 52.5 |
| 100 | [19.44, 480.12] | [15.517, 1421.2] | 8150 | 25.0, 91.6 | 425, 525 |
| ∞ (claim) | [19.197, ≥ 6906.7] | [15.517, ≥ 959,689.9] | M²K0 | 360, 61,848 | — |

Inner maxima are local optima and may be beaten. The ceiling is rigorous but not shown to be tight. The lower
end stays pinned near τ_CS for every F ≥ 10. Reading: a factor-2 amplitude uncertainty, typical of DFT |Φ|²,
leaves the memory within a factor ≤ 4.8 (x) and ≤ 5.9 (z) of the reference (rigorous). The "2–6 orders of magnitude" require
amplitude priors looser than a factor 10²–10⁴ per event.

### 3.4 Smooth amplitude prior — SCOPE QUALIFICATION (inner evidence, not proof)
Class: g = g_ref·exp P(ω_p, ω_a, ω_b), with P a polynomial symmetric under a↔b; the data are as in the claim.
This is the class that 13-final-report §4 names as missing. Energy conservation makes ω_p nearly a function of
ω_a + ω_b, so the numerical rank saturates (121 at degree 10, 155 at degree 12, 197 at 16, 201 at 20, out of
2,627 orbits).

* Degree 6 (49 coefficients < 119 constraints): the reference is isolated, so τ is pinned exactly.
* Degree 10: τ_x 34.2530 and τ_z 27.4465 (|P| ≤ 0.004); no movement.
* Degrees 12–20, active-set projected gradient (`cx_attack1_smooth_pg.py`):

| degree | bound on \|P\| | τ_x max found (ps) | τ_z max found (ps) |
|---|---|---|---|
| 12 | ln 10 / ln 100 / 30 | 34.48 / 34.85 / 40.23 | 27.65 / 27.80 / 30.08 |
| 16 | ln 10 / ln 100 / 30 | 35.22 / 36.64 / 47.65 | 28.08 / 29.37 / 35.35 |
| 20 | 30 | 43.19 | 49.43 |

Even when amplitudes may vary by e^±30, smoothness in the energies kept τ below 1.8× the reference.

### 3.5 What the large-τ witnesses are physically — no sum rule violated, reading weakened
Spectral anatomy (`cx_verify_reference.json`, `cx_attack1_slowmode.json`):
* 6906.7 ps: a doubly degenerate mode at λ = 3.84e-6 THz (relaxation 20.7 ns) carries 33.2 % of K_xx and
  99.6 % of N. Weight: 63 % below 5 THz, 19 % at 5–10 THz, 18 % above 10 THz. Its best squared overlap with
  crystal momentum restricted to f < f_c is 0.30, so it is not a phonon-hydrodynamic quasi-momentum.
* 959,689.9 ps: λ = 4.4e-8 THz (1.81 μs) carries 53 % of K_zz; momentum overlap ≤ 0.10.
* 3.58e7 ps: λ = 2.46e-15 THz (relaxation 3.2e13 ps ≈ 32 s) carries 1.1e-6 of K and 99.993 % of N. This
  witness is a first-moment artifact: negligible spectral weight at an absurd time.

Energy remains the only exact invariant, and the conductivity f-sum rule holds for any C, so no conservation
law or sum rule is violated. The near-invariants are engineered by suppressing 49–87 % of the events by more
than 10⁶. This is admissible in the class and implausible for a P6₃mc anharmonic vertex, but plausibility is
not a theorem.

### 3.6 "Without symmetry" witnesses — BROKEN as physical statements
The unconstrained basal witness has τ_x = 114,019 ps but τ_y = 514.0 ps, and τ_z = 35.3 ps. Its rates break
the 6mm symmetry by a relative norm of 0.50. A hexagonal crystal requires basal isotropy of the whole response
function, so τ_xx = τ_yy for any physical operator. The interval "x [19.19, ≥1.14e5] ps" is therefore a
statement about a symmetry-violating class only. It also uses forbidden events (share up to 0.82).

---------------------------------------------------------------------------------------------------------------

## 4. Attack 2 — reference consistency and mesh

**Ω′ comparison — SURVIVES** (`cx_attack2_reference.py/.json`). I rebuilt Ω′ = D + C1 − (C0+C2)J from the
stored channels: degeneracy-averaged, symmetrised, 897 modes. Against the event operator C(g_ref):
* diagonal: median relative deviation 7.2e-4, maximum 0.51;
* off-diagonal: Frobenius relative difference 7.0 %;
* lowest eigenvalues on H agree to 0.1 % (9.7643e-4 vs 9.7707e-4, …);
* solutions: |x_ev − x_Ω′|/|x_Ω′| = 0.15 % (x) and 0.79 % (z).

The Ω′ moments do not depend on how its non-conserved energy direction is handled. Projection onto e⊥, the raw
inverse and spectral removal of the eigenvalue 1.06e-6 agree to 1e-12, because b's overlap with that mode is
≤ 1e-20. The agreement is structural, not coincidental.

Minor discrepancy: my Ω′ gives κ_xx 214.277 and τ_x 34.214 ps (matching 07-experiments), but κ_zz 202.185 and
τ_z 27.535 ps, not the quoted 202.82 and 27.99. The event operator (201.967, 27.446) is closer to this rebuild
than the report states.

**Mesh — SCOPE** (`cx_attack2_mesh775.py/.json`). I exported 7×7×5 with the campaign's export script:
n = 2,937, m = 481,899, 23,674 event orbits, 285 mode orbits, 4.61 % forbidden.

| | κ_xx / κ_zz (W/m K) | τ_ref x / z (ps) | τ_CS x / z (ps) | λ_min (THz) | factor-F ceiling |
|---|---|---|---|---|---|
| 5×5×3 | 214.30 / 201.97 | 34.25 / 27.45 | 19.19 / 15.52 | 9.76e-4 | 81.5·F ps |
| 7×7×5 | 243.03 / 272.33 | 52.61 / 50.62 | 21.10 / 20.02 | 4.04e-4 | 197.2·F ps |

The AlN memory numbers are mesh properties, not AlN properties: τ_ref changes by 54–84 %. The forbidden-event
structure persists, and the room left by a factor-F prior grows under refinement (ceiling/τ_CS: 4.2F → 9.3F for
x). A large-τ optimisation at 7×7×5 was not run.

---------------------------------------------------------------------------------------------------------------

## 5. Attacks 3 and 4

### 5.1 Attack 3 — numerical artifacts: SURVIVES
* **Two routes.** Cholesky of C + γêêᵀ and eigendecomposition on an explicit Householder basis of H reproduce
  every stored witness: 6906.747185 / 6906.747185, 959,689.9306 / 959,689.9304, 19.196925, 15.516870 ps. The
  near-singular witness gives 3.5746e7 / 3.5854e7 (0.3 % float spread at cond 6.4e13), as the campaign
  reported.
* **Exact energy conservation.** Ball arithmetic (256 bits) with l_p = √((f_a+f_b)/f_p) evaluated inside arb,
  so that sᵀε = 0 exactly (`cx_attack3_arb.py`). Near-singular witness: 35,759,344.50888929 ps (radius 2.5e-53).
  This is identical to all printed digits with the float-coefficient certification, so float energy
  conservation does not inflate τ. Robust witness: 6906.747184851676 ps.
* **Data perturbation.** Frequencies perturbed by 1e-10 (relative, random) and propagated consistently to D,
  l_p, b and e: near-singular witness 3.5675e7, 3.6140e7, 3.5497e7 ps (±1.1 %); robust witness 6906.747182 ps.
  The hidden invariant comes from a pattern of exact zeros, not from rounding.
* **Feasibility.** Lifetimes hold to ≤ 5e-13 and K to ≤ 6e-14 for the well-conditioned witnesses. The 3.58e7
  witness has the stated 1.1e-11 K_xx mismatch.

### 5.2 Attack 4a — observables the campaign does not fix (`cx_attack4_observables.py/.json`)
Here z is the Laplace variable in operator units, z = z_phys/(4π).

| quantity | reference x | witness 6906.7 (x) | reference z | witness 959,690 (z) |
|---|---|---|---|---|
| K(z)/K at z_phys = 1/(10 ns), 1/ns, 1/(100 ps) | 0.9966, 0.9674, 0.7662 | 0.7731, 0.6555, 0.5016 | 0.9973, 0.9736, 0.7988 | 0.4712, 0.4583, 0.3931 |
| accumulation f < 5 THz, f < 10 THz | 0.239, 0.777 | 3.071, 2.402 | 0.117, 0.898 | 0.086, **8.500** |
| second moment bᵀCb / reference | 1 | 1.806 | 1 | 1.860 |
| κ(T)/κ_ref(T), T-independent \|Φ\|²: 100, 200, 500, 800 K | 1 | 69.2, 4.70, 0.90, 1.54 (κ_x) | 1 | 1476, 97.4, 29.3, 61.1 (κ_z) |
| implied lifetimes at 100 K / reference | 1 | 0.19 – 4.70 | 1 | 0.14 – 3.06 |

The witnesses are therefore not invisible. Their K(z) lies 22–53 % below the reference at Laplace rates
1/(10 ns)–1/(100 ps). They also predict mode-resolved accumulations with ±800 % cancellations and κ(T) values off
by up to 10³ at other temperatures.

### 5.3 Attack 4b — information in multi-temperature lifetimes (rank)
With |Φ_α|²δ_α temperature independent, lifetimes at temperatures T_k give linear constraints on the 2,292
allowed orbit rates (`cx_attack4_observables.json`):

| temperature set | rows | numerical rank at relative tolerance 1e-3 / 1e-6 / 1e-13 |
|---|---|---|
| 300 K | 117 | 117 / 117 / 117 |
| 100, 200, 300, 500, 800 K | 585 | 289 / 490 / 585 |
| 10 values in 50–1500 K | 1,170 | 344 / 635 / 1,074 |
| 40 values in 20–5000 K | 4,680 | 406 / 777 / 1,635 |

Temperature dependence cannot pin |Φ|²: even 40 exact temperatures leave at least 657 undetermined rate
combinations, and at 0.1 % data accuracy at least 1,886.

### 5.4 Attack 4c — re-optimisation with additional data fixed
**κ(T) curve** (`cx_runs/mK5_*`, verified in `cx_verify_multiT.json`). Data: lifetimes at 300 K, K_xx and
K_zz at 100, 200, 300, 500 and 800 K, T-independent |Φ|²δ, selection rules enforced, 3000 s budget. The gap
does **not** collapse:

| start | τ_x (ps) | τ_z (ps) |
|---|---|---|
| reference | 309.4 | 12,703 |
| witness pattern | 675.1 | **131,918** |

All ten K values match to ≤ 2e-12. The c-axis witness has 46 % of K_zz in modes slower than 1 ns; its
implied lifetimes at 100 K are 0.11–3.9× the reference. Runs were stopped by the time limit, so these are lower
bounds.

**Other data, re-optimised** (`cx_runs/Kz_*`, `acc_*`, `mT5*`, `mT5tol6_*`; verified in `cx_verify_extras.json`).
All runs also fix the 300 K lifetimes and K_xx, K_zz; values are inner bounds within a 1500 s / 1200 s budget.

| additional data fixed | class | start | τ_x found (ps) | τ_z found (ps) | remark |
|---|---|---|---|---|---|
| K(z) of x and z at z_phys = 1/ns, 1/(100 ps), 1/(10 ps) (6 values) | claim's | witness pattern | 62.85 | 33.87 | slow modes persist (29 ns, 13 ns) with only 0.10 % / 0.05 % of K |
| accumulation below 5, 10, 15 THz, x and z (6 values) | claim's | witness pattern | 1117.6 | 61,254 | 19 % / 48 % of K in modes slower than 1 ns |
| κ(T) curve: K_xx, K_zz at 100–800 K | selection rules | witness / reference | 675.1 / 309.4 | 131,918 / 12,703 | table above |
| lifetimes and K at 100, 200, 300, 500, 800 K (595 constraints) | selection rules | reference | 35.93 | 28.44 | data matched to 1e-6; still creeping at about 0.002 ps per iteration |
| the same, witness or homotopy starts | selection rules | — | — | — | restoration failed: the witness patterns imply 100 K lifetimes off by 0.1–5× |

* **K(z) at three Laplace rates** was the most effective discriminator found within the budget. It drops the
  found maxima from ≥ 6906.7 / ≥ 959,690 ps to 62.9 / 33.9 ps, i.e. 1.8× and 1.2× the reference. This is
  inner evidence, not a bound. Spectrally, K(z) = Σ_k κ_k/(1 + z T_k) and τ = Σ_k κ_k T_k / Σ_k κ_k, with
  T_k = 1/λ_k and κ_k = (b·u_k)²/λ_k. A mode carrying a fraction w of K at time T is invisible in K(z) when
  zT ≫ 1, apart from the deficit w·K the other modes must cover, yet it adds w·T to τ. So finitely many K(z)
  cannot bound τ_mem in a general spectral class. The K(z)-constrained witnesses found here retain exactly such
  slow modes (29 ns), held at 0.1 % weight. Whether the event class admits much larger weight × time products
  under K(z) constraints was not settled.
* **Accumulation data** do not collapse the gap (τ_z ≥ 61,254 ps).
* **Lifetimes at five temperatures: inconclusive.**
  * The stacked constraint Jacobian has 17 singular values below 1e-8 of its maximum (smallest 2.7e-10
    relative).
  * With an exact tangent projection and truncated restoration, my first-order method crept from 34.25 to
    35.93 ps (x) and from 27.45 to 28.44 ps (z) in 20 minutes.
  * The exact feasible set has dimension ≥ 1,697 (§5.3), so the rates are certainly not identified. Whether
    τ_mem is identified under multi-temperature lifetime data is neither shown nor refuted here.

---------------------------------------------------------------------------------------------------------------

## 6. Attack 5 — L3 (K ≥ K_RTA/s_max): SURVIVES as stated; SCOPE on the meaning of r

**Random and adversarial search** (`cx_attack5_L3.py/.json`):
* networks: exact-energy integer networks of 4–11 modes with three-phonon events, repeated daughters
  (p → 2a), four-phonon events mixed with three-phonon ones (p+q → a+b, 2p → a+b), and conserving
  non-integer coefficients;
* rates: log-normal up to σ = 8 (near-singular C);
* currents: random, and adversarial (the generalised eigenvector minimising K/K_RTA on H);
* plus a Hypothesis property test with 400 examples.

Result:
* the minimum of K·s_max/K_RTA is exactly 1.0000 in every family (adversarial currents; sharpness confirmed);
* 0 Hypothesis failures;
* pure 2↔2 networks are never admissible, because phonon number becomes an extra invariant;
* zero modes are outside the entropy-coordinate setting (D_μ = ∞).

The Cauchy–Schwarz double-counting proof holds line by line: b_μ² ≤ r_μ D_μ and Σ_μ D_μ ≤ s_max·xᵀCx = s_max K.

**The linewidth convention is load-bearing** (`cx_attack5_rse_aln.py/.json`, `cx_attack5_chain_exact.py/.json`).
* On the AlN data, (diag C − Γ_phono3py)/Γ_phono3py = 0.503 × (repeated-daughter share of diag C), with
  correlation 0.992 over 66 modes and up to 3.8 %. phono3py's self-energy linewidth Γ counts a repeated-daughter
  self-coupling with half the weight that diag C does.
* L3 is a statement about r = diag C, not about the linewidth an RTA code reports.
* Exact counterexample to the Γ reading: chain with energies 1, 2, 4, 8, 16 and events 2→1+1, 4→2+2, 8→4+4,
  16→8+8 (s_max = 2, D = I, g = 1), with the rational current b ⟂ e given in the json:
  * K/K_RTA(Γ-type) = 7383516975/22150549363 ≈ 0.33333 < 1/2;
  * K/K_RTA(diag C) = 10336923765/19122114037 ≈ 0.54057 ≥ 1/2;
  * the infimum over b of the Γ-type ratio is 1/3.
* For AlN the difference is negligible: K/K_RTA(Γ) is 1.2500 (x) and 1.4203 (z), against 1.2572 and 1.4319.

---------------------------------------------------------------------------------------------------------------

## 7. Wording the claim must carry to survive

1. **Class.** "Within the event-cone class — arbitrary non-negative rates on the energy-allowed three-phonon event
   set, with no amplitude model …" must stay in every summary sentence. "Realistic event geometries" may refer
   to the event set only, not to the rates.
2. **Selection rules.** Either restrict the event set to symmetry-allowed vertices and replace the AlN witnesses,
   or state that the class contains selection-rule-forbidden processes (7.4 % of events). Selection-rule-
   respecting replacements certified here: τ_x ≥ 7038.17 ps (data to 3.4e-14) and τ_z ≥ 4,322,183 ps (data to
   2.9e-11).
3. **Prior dependence (rigorous).** "Any prior that bounds the suppression of each event rate to a factor F
   relative to a reference operator confines τ_mem to [τ_CS, F/λ_min(C_ref)] — here [τ_CS, 81.5·F] ps. The 2–5
   orders of magnitude require per-event suppressions of 10²–10⁴ and more. Under a smooth-in-energy amplitude
   prior no witness above 1.8× the reference was found."
4. **Magnitude.** "10²–10⁶" should read 10²–10⁵. The 10⁶ rests on the near-singular 3.58e7 ps witness, which is
   a 1e-6-weight, 32-second mode and which 13-final-report itself excludes. Robust well-conditioned factors over
   the reference are: 202 (x, campaign), 205 (x, this review), 3.5e4 (z, campaign) and 1.57e5 (z, this review;
   cond 1.1e7).
5. **Symmetry.** The "without symmetry" intervals violate basal isotropy (τ_x/τ_y = 222) and should not be quoted
   as AlN properties.
6. **Observability.** The large-τ witnesses differ strongly from the reference in K(z) at Laplace rates
   1/(10 ns)–1/(100 ps), in the
   accumulation function, and in κ(T) at other temperatures. Fixing κ(T) at five temperatures still leaves
   τ_z ≥ 1.3e5 ps.
7. **L3.** "r = diag C (single-mode relaxation rates including self-coupling of repeated slots); the bound may
   fail for self-energy linewidths when repeated-daughter events are present."
8. **Mesh.** The AlN interval endpoints are 5×5×3 values; τ_ref changes by 54–84 % at 7×7×5.

---------------------------------------------------------------------------------------------------------------

## 8. Attacks not performed or incomplete

* **fc3-parametrised class** (rates generated by symmetric anharmonic force constants): not computed. The rank of
  the fc3 → |Φ|² map relative to the 2,292 allowed orbits is unknown. This is the physically decisive class.
* **Event-by-event proof of the selection rules** (eigenvector irreps or a random symmetric fc3): not done.
  The statistical evidence is in §3.2.
* **Global optima.** No maximum, mine or the campaign's, is shown to be global. Random starts mostly failed to
  restore inside tight boxes. The factor-F ceiling is rigorous but its tightness is untested (no SDP or other
  relaxation attempted).
* **7×7×5 optimisation** not run; only the reference and the ceiling were computed there.
* **Complex-frequency response** (FDTR phase) not computed; only real Laplace frequencies were used.
* **Multi-temperature lifetimes**: inconclusive. The constraint set is numerically stiff (17 near-dependent
  combinations) and my first-order method crept by < 5 %. A second-order or interior-point method, or
  inequality data at realistic 1–10 % accuracy, are needed.
* **K(z)-invisible slow modes**: whether the event class admits slow modes of tiny weight and huge time while
  K(z) is fixed was not constructed (the spectral argument says it is not excluded in general).
* **Second-moment constraint** not imposed in a re-optimisation. It is a single linear constraint per direction;
  the witnesses change it by ×1.8.
* **proof-audit.md** was not read; no reconciliation with that audit is attempted here.

---------------------------------------------------------------------------------------------------------------

## 9. Verdict

The mathematical result survives: within the declared event-cone class, static data (300 K lifetimes and the
DC tensor) pin τ_mem from below (τ_CS, nearly attained) and not from above. Certified witnesses respecting the
crystal selection rules exceed the reference by 205× (basal) and 1.57×10⁵ (c axis). Adding the κ(T) curve at
five temperatures does not restore identifiability. No numerical artifact was found.

The physical reading does not survive as written.
* The quoted AlN witnesses use symmetry-forbidden processes.
* The non-symmetric witnesses break basal isotropy.
* The gap is entirely a consequence of having no amplitude prior. A per-event factor-F prior caps τ at
  81.5·F ps (rigorous), and a smooth amplitude prior kept τ within 1.8× the reference.
* The witnesses are strongly discriminated by GHz conductivity dispersion, accumulation, and multi-temperature
  lifetimes. Re-optimising with K(z) fixed at three Laplace rates found nothing above 1.8× the reference (no
  bound). With the κ(T) curve or the accumulation function fixed, witnesses of 1.3×10⁵ and 6.1×10⁴ ps persist.

L3 is correct for r = diag C and false for self-energy linewidths in networks with repeated daughters.

Status: claim **SURVIVES WITH SCOPE QUALIFICATIONS**. The headline witnesses are **BROKEN** (to be replaced). The
"10⁶" is **CONTRADICTED** as a robust figure; use 10⁵. L3 **SURVIVES**: its argument checks line by line
and no counterexample was found in scope; the r-convention must be made explicit.
