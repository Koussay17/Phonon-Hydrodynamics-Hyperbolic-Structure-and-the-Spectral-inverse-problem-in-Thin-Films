# Counterexample review — fc3-parametrised memory campaign (20261005-fc3-memory)

**Date:** 10 October 2026. **Object:** claims C1–C4 of the campaign as stated in `00-question.md`, `02-assumptions.md`,
`07-experiments.md` (E1–E18) and `09-failed-approaches.md` (D: copies of 10 October). **Mandate:** find errors,
artefacts and counterexamples; report what survives. All scripts and raw outputs of this review are in this folder
(`review/`); no campaign file outside `review/` was modified. Large intermediates are in `review/work/`.

Environments: research Python (NumPy 2.5.3, SciPy 1.18.1) for analysis and optimisation; phono3py 4.5.0 replica
environment for everything that calls phono3py (inputs hash-checked by `aln_common.verify_inputs`).

---

## 0. Verdicts

| claim | verdict | decisive numbers |
|---|---|---|
| **C3** mapping (polarity, images, folding, units) | **CONFIRMED** | independent re-derivation reproduces the folded tensor to 2.1e-16; FC2 mapped with the same conventions agrees with phono3py's fc2 to 0.72 % |
| **C3** "delta_o are physical differences, not artefacts" | **REFUTED in part** | on 266 of 456 orbits (58 %) the folded ALAMODE tensor is identically zero: delta_o = 1 is ALAMODE's cubic cutoff (largest cluster 5.678 A), not a measured disagreement; these orbits carry 45 % of the k = 1 box (sum of squared radii) |
| **C4** route B reproduces route A to <= 0.03 % | **CONFIRMED** (reproduced with an independent operator code) | REF and witnesses reproduced to every printed digit |
| **C4** "verified independently" (data consistency of the witnesses) | **REFUTED as worded** | under phono3py's own LBTE operator the witnesses violate the +-1 % kappa tolerance: kappa_xx +2.21 % (k = 1, max tau_x), +4.16 % (k = 2), kappa_zz +1.54 % / +3.29 %; linewidths up to +13 %. Cause identified: finite-sigma (detuned-event) weighting convention |
| **C1** range [0.84, 1.31] (x), [0.81, 1.37] (z) at eps = 0.1 | **REFUTED** (upper ends; lower ends see section 3) | feasible, independently verified witnesses at 47.514 ps (1.387, x) and 39.396 ps (1.436, z) versus the claimed 44.80 / 37.58 ps |
| **C2** numbers [0.70, 2.35] / [0.68, 2.15] (k = 1), [0.64, 3.72] / [0.61, 5.96] (k = 2) | **CONFIRMED as 5x5x3 event-operator inner bounds** | reproduced; same witnesses |
| **C2** evidence: the witnesses as physical (mesh-robust) counterexamples | **REFUTED** | the same witness tensors on the 7x7x5 event set: tau_x 0.83 x reference (k = 1 max x), 0.71 (k = 2 max x), tau_z 0.53 (k = 2 max z); kappa off by -8 % to -30 %; lifetimes up to +338 %. Confirmed with phono3py's own LBTE at 7x7x5. Also kappa(100 K) off by +23 % to +217 % |
| **C2** interpretation "static data do not pin tau_mem better than a factor ~2 at realistic fc3 precision; memory carried by small long-range constants" | **UNRESOLVED (unsupported)** — holds only as a statement about the 5x5x3 event operator | its only evidence fails the mesh, convention and temperature tests; whether mesh- and convention-robust tensors with a factor ~2 exist was not tested; 45 % of the k = 1 box is a cutoff artefact (removing it: max ratios 1.86 / 1.82 instead of 2.35 / 2.15) |

Details, numbers and the attacks that failed follow.

---

## 1. C3 — the ALAMODE mapping

### 1.1 Mapping mechanics: CONFIRMED
* **Independent re-derivation** (`r3_mapping.py`; ElementTree parsing, Cartesian displacement vectors, mirror applied
  to vectors, nearest-atom matching in the 3x3x2 supercell, per-axis sign product; no code shared with
  `alamode_compare.py`): the folded raw tensor agrees with `alamode/fc3_alamode_raw.npy` to
  **2.1e-16** (max abs, relative); sum rule 7.4e-13; matching error of atom positions 8.5e-6 (fractional).
* **Polarity**: the axial Al-N bond points along **+1.9035 A** in the ALAMODE cell and **-1.9032 A** in the phono3py
  cell: the z -> -z map is required. Without it the N sublattice does not match (residual 0.238 in fractional units).
* **Image ordering**: the ALM ordering ((0,0,0) first, then ia, ja, ka in -1..1, ia slowest) is consistent with the
  data (e.g. FC2 `pair2="3 1 14"` is the c-axis neighbour below, value equal to its partner above). Every FC3 entry
  carries exactly **one** image (no multi-image clusters: all ALAMODE clusters are shorter than 5.68 A, inside the
  4x4x3 Wigner-Seitz cell), so the fold has no multiplicity ambiguity.
* **Units**: on-site spot values -67.72 (ALAMODE) vs -68.24 eV/A^3 (reference), 32.31 vs 32.38.
* **Second tensor order (FC2 test)** (`r12_fc2_mapping.py`): ALAMODE `harmonic443.xml` mapped with the same
  conventions onto the phono3py 5x5x3 fc2: relative difference **0.72 %** over 1,404 pairs < 4.5 A; the basal
  nearest-neighbour block reproduces the signs of the odd-z components (xz -1.734 vs -1.739, yz 1.001 vs 1.004).

### 1.2 Are the delta_o measured differences? REFUTED for 58 % of the orbits
`r1_orbit_anatomy.py`, `r3_mapping.py`:
* The folded ALAMODE tensor is **identically zero on 266 of the 456 cluster orbits** (|eta_A,o| <= 1e-9 |eta_0,o|),
  so delta_o = 1 there by construction. Per label shell: 2-3.2 A: 2 of 24; 3.2-4 A: 4 of 20; 4-5 A: 16 of 82;
  5-6 A: 104 of 184; > 6 A: 140 of 140. E13 attributes the cutoff zeros to "> 6 A" only; this is inaccurate
  (126 cutoff-zero orbits lie below 6 A by the campaign's labels).
* Cause: ALAMODE's cubic cutoff. The largest ALAMODE cluster diameter is **5.678 A**. The 266 zero orbits are
  **exactly** the orbits whose consistent (single-embedding) diameter is >= 5.718 A; every nonzero orbit has a
  consistent diameter <= 5.677 A. delta_o = 1 on these orbits measures the reference norm, not a disagreement.
* These orbits hold **45 %** of the k = 1 box (sum over orbits of squared absolute radii max(delta_o, 0.005) |eta_0,o|).
  Shell shares of that box: 4-5 A 25 %, 5-6 A 45 %, > 6 A 11 %.
* Where both workflows have constants, they are correlated, not noise: median per-orbit cosine between eta_A and
  eta_0 is 0.99996 (on-site), 0.99998 (< 2 A), 0.992 (2-3.2 A), 0.955 (3.2-4 A), 0.88 (4-5 A), 0.66 (5-6 A;
  10 % quantile 0.19). The long-range differences on nonzero orbits are therefore partly physical, increasingly noisy
  beyond 5 A.
* **Shell labels are not cluster diameters.** The label is the largest of three independently chosen
  shortest-image distances. For **60** orbits no single embedding of the triplet realises the three distances; e.g.
  orbits 182 and 208 (label 3.113 A, "2-3.2 A shell") are collinear triplets i, i + a1, i + 2a1 = i - a1 whose
  physical diameter is 2a = **6.23 A** (both ALAMODE-zero, delta_o = 1). Re-binning the witness anatomies by the
  consistent diameter changes the shell numbers little (k = 1 max tau_x: "2-3.2 A" 7.1 % -> 7.0 %), so this is a
  labelling inaccuracy with small quantitative effect.

### 1.3 Folding and aliasing
Folding is the correct emulation of what a 3x3x2 displacement calculation would see **if** ALAMODE's constants were
the infinite-crystal constants; it does not inflate delta_o. What inflates the long-range delta_o is (i) the cutoff
(section 1.2) and (ii) the fact that the 3x3x2 reference folds in all physical clusters beyond 5.68 A while the folded
ALAMODE tensor cannot. Whether the reference's own long-range constants are noise could not be decided from one
pair of workflows; the cosines above say they are mostly signal up to 5 A.

---

## 2. C4 — independence of the witness verification

### 2.1 The stated agreement: CONFIRMED with an independent code
`r5_event_operator.py` builds the event operator from the written definitions (own orbit tying from `rot_maps`,
own dense assembly, own eigensolver; no campaign module imported). From the phono3py interaction strengths of the
campaign's pp files it reproduces: REF kappa_xx 214.29993549566626 (campaign 214.29993549566584), tau_x
34.252990795285 ps, kappa_zz 201.96695921348, tau_z 27.446471103010 ps; witness cal1_c0_max tau_x 80.32581 ps
(route B 80.3258), kappa ratio 1.0098001; cal2_c2_max tau_z 163.3993 ps (route B 163.3993). Route B is therefore a
correct evaluation of the campaign's operator.

### 2.2 Data consistency is not convention-independent: REFUTED as "verified independently"
Route B shares the event list, the conserving stoichiometry and the per-event Bose weights with the optimisation
model; it verifies the amplitude model, not the operator. A genuinely independent route is phono3py's own LBTE
(`r4_p3_lbte.py`: reducible collision matrix on the full 5x5x3 mesh, all triplets, phono3py's assembly and
symmetrisation; tau computed from the captured Omega and X). It reproduces the reference (kappa 214.277 / 202.191
vs 214.300 / 201.967; tau 34.217 / 27.538 ps vs 34.253 / 27.446 ps), and gives for the witnesses
(ratios to the phono3py reference; event-operator value in brackets):

| witness | kappa_xx | kappa_zz | tau_x | tau_z | linewidth change (deg.-avg.) |
|---|---|---|---|---|---|
| ALAMODE | 1.0014 (1.0011) | 1.0330 (1.0329) | 0.9935 | 1.0171 | -18.3 % ... +17.6 % |
| cal1 max x | **1.0221** (1.0098) | 0.9897 | 2.429 (2.345) | 1.092 | -10.5 % ... +11.1 % |
| cal1 max z | 1.0065 | **1.0154** (1.0099) | 0.992 | 2.181 (2.151) | -10.2 % ... +10.6 % |
| cal1 min x | 0.9909 | 1.0094 | 0.699 | 1.056 | -10.0 % ... +10.1 % |
| cal1 min z | **1.0110** (1.0096) | 0.9892 | 1.035 | 0.683 | -10.1 % ... +10.3 % |
| cal2 max x | **1.0416** (1.0098) | 0.9958 | 4.066 (3.721) | 1.132 | -10.6 % ... +11.5 % |
| cal2 max z | 1.0095 | **1.0329** (1.0096) | 1.085 | 6.328 (5.953) | -10.1 % ... +11.6 % |
| cal2 min x | 0.9912 | 1.0099 | 0.638 | 1.197 | -10.1 % ... +13.1 % |
| orb 0.1 max x | **1.0122** (1.0098) | 0.9901 | 1.297 | 1.062 | -10.4 % ... +10.2 % |
| orb 0.2 max x | **1.0154** | 0.9901 | 1.647 | 1.085 | -10.6 % ... +10.8 % |

The memory enhancement survives (it is even larger), but the kappa tolerance (+-1 %) is violated by 0.2 to 3.2
percentage points beyond the boundary the witnesses saturate. The reference itself differs between the two operators
by 1e-4 (xx) and 1.1e-3 (zz): the witnesses amplify the convention dependence by one to two orders of magnitude.

**Mechanism (identified).** Excluded: the 4-sigma event cut (phono3py with `sigma_cutoff = 4`: identical ratios,
1.0221 / 4.066); phono3py's internal options (`is_kappa_star` False, irreducible solver: identical, 219.006 W/mK);
degeneracy averaging (phono3py-style row/column averaging applied to the event operator leaves every kappa ratio
unchanged to 2e-4, `r6b`); permutation asymmetry of |Phi|^2 (the three perspectives of each event agree to
<= 6e-5 on degenerate-block sums, `r7`). Identified: the weighting of **detuned events** (finite sigma). An event
operator that differs from the campaign's only by using phono3py's Bose-factor forms for detuned events
(`r6d_detuning_convention.py`; identical for exactly conserving events) reproduces phono3py's numbers:
kappa_xx ratio 1.02212 (phono3py 1.02207) for cal1 max x, 1.04173 (1.04163) for cal2 max x, kappa_zz 1.03272
(1.0329) for cal2 max z, tau_x ratio 2.4297 (2.4293). The campaign's conserving rescaling of the stoichiometry
and the symmetric sinh weights are one admissible discretisation of the smeared LBTE; phono3py's is another; the
witnesses are data-consistent in the first and not in the second.

---

## 3. C1 — uniform per-cluster bound eps = 0.1

### 3.1 Optimiser budget
In every campaign run (t3_orb, t3_ms, t3_cal) **every** inner L-BFGS-B solve stopped at `maxiter = 200`
(`nit = 200` in all 15 outer iterations), while the penalty parameter rose to 8e5 - 3e6. The tau plateau across
outer iterations reflects the stiffening penalty, not a verified stationary point. The 2.6 % agreement of four
starts is therefore weak evidence about the maximum.

### 3.2 Counterexamples to the stated upper ends
Start: the eps = 0.2 witnesses scaled by 1/2 in the scaled coordinates (every orbit then sits exactly at 0.0999),
same optimiser and settings except `maxiter = 600` (`r2_opt.py`, runs `runs/c1_*_cont02`):

| target | claimed (E15/E16) | found | ratio | verification (phono3py pp + independent operator) |
|---|---|---|---|---|
| max tau_x | 44.80 ps (1.308) | **47.514 ps** | **1.387** | tau_x 47.5144, kappa 1.00980 / 0.99021, deg.-avg. lifetimes within 9.980 %, per-orbit change on the rebuilt tensor 0.09990, sum rule 1.5e-15 |
| max tau_z | 37.58 ps (1.370) | **39.396 ps** (orbit average 39.408) | **1.436** | kappa 1.00980 / 1.00994 (orbit-average convention), lifetimes within 9.980 %, per-orbit 0.09990 |

PENDING: polish (maxiter 1000) and minima continuation.

The stated range "[0.84, 1.31] (x), [0.81, 1.37] (z)" is therefore not the attainable range; the upper ends are at
least 1.39 and 1.44, and no evidence bounds them further.

### 3.3 Transfer of the C1 witnesses
7x7x5 event set (section 4.2 method): orb 0.1 max x: tau_x 1.087, kappa_xx +1.8 %; orb 0.1 max z: tau_z 1.347,
kappa_zz +3.1 %, kappa_xx +2.3 %; new max x: tau_x 1.048, kappa +1.0 % / +1.4 %, lifetimes -36 % ... +23 %; new
max z: tau_z 1.315, kappa +3.3 % / +2.4 %, lifetimes -19 % ... +26 %. The c-axis enhancement of the per-cluster
class partly survives the mesh change; the data constraints do not.

---

## 4. C2 — calibrated class

### 4.1 Numbers: CONFIRMED as 5x5x3 inner bounds
All four k = 1 and k = 2 maxima, and the minima, are reproduced by the independent operator from the campaign's
pp files (section 2.1). The rebuilt tensors are symmetric, satisfy the sum rule and the per-orbit bounds.

### 4.2 Mesh transfer: the witnesses are artefacts of the 5x5x3 sampling
An fc3 tensor is mesh independent; a witness for a statement about AlN must keep the data and the memory effect on a
finer sampling of the same smeared operator. Each witness tensor was evaluated on the 7x7x5 event set of the
memory-gap campaign (481,899 events, sigma 0.1 THz, 4 sigma, 300 K) with phono3py's Rust amplitudes
(`p3_events.py`, 37 s per tensor) and the independent operator (`r8_mesh_transfer.py`), ratios to the 7x7x5
reference (243.03 / 272.33 W/mK, 52.61 / 50.62 ps, identical to E1):

| witness (5x5x3 tau ratio) | 7x7x5 kappa_xx | 7x7x5 kappa_zz | 7x7x5 tau_x | 7x7x5 tau_z | 7x7x5 lifetimes (deg.-avg.) |
|---|---|---|---|---|---|
| ALAMODE (0.99 / 1.02) | 0.9965 | 1.0144 | 0.977 | 1.006 | -27 % ... +17 % (rms 3.8 %) |
| k = 1 max x (**2.345**) | 0.916 | 0.932 | **0.830** | 0.918 | -31 % ... +44 % |
| k = 1 max z (**2.151**) | 1.019 | 0.965 | 0.898 | **1.081** | -30 % ... +54 % |
| k = 1 min x (0.698) | 0.929 | 0.912 | 0.689 | 0.700 | -17 % ... +115 % |
| k = 1 min z (0.684) | 0.931 | 0.838 | 0.719 | 0.539 | -66 % ... +261 % |
| k = 2 max x (**3.721**) | 0.853 | 0.827 | **0.705** | 0.643 | -47 % ... +101 % |
| k = 2 max z (**5.953**) | 0.994 | 0.703 | 0.786 | **0.527** | -37 % ... +338 % |
| k = 1, cutoff orbits at floor, max x (1.858) | 0.933 | 1.019 | **0.872** | 1.330 | -29 % ... +31 % |
| k = 1, cutoff orbits at floor, max z (1.820) | 1.039 | 1.002 | 0.997 | **1.184** | -29 % ... +48 % |

**Independent confirmation** with phono3py's own LBTE at 7x7x5 (`r4_p3_lbte_mesh.py`): k = 1 max x kappa_xx 0.9157,
tau_x 0.8306; k = 1 max z kappa_zz 0.9659, tau_z 1.0841; k = 2 max x kappa 0.854 / 0.827, tau_x 0.707; k = 2 max z
kappa_zz 0.703, tau_z 0.527 (agreement with the event operator to about three digits).

The maxima found on 5x5x3 (2.35, 2.15, 3.72, 5.95) **reverse** on 7x7x5 for x (0.83, 0.71) and nearly vanish or reverse
for z (1.08, 0.53), while the 300 K data are violated by 3.5 % - 30 % in kappa and by up to +338 % in single-mode
lifetimes. A physical alternative tensor (ALAMODE) changes by similar amounts on both meshes. The witnesses achieve
their 5x5x3 effect by fine-tuning amplitudes at the 75-point mesh's triplets (48 - 71 orbit rates suppressed by more
than 10x; smallest single-event |Phi|^2 ratio to the reference 1e-5 - 1e-3 for k = 1 and 3e-9 - 3e-6 for k = 2,
against 1.3e-3 for ALAMODE), which a finer mesh does not reproduce.

### 4.3 Temperature transfer
Constraints were imposed at 300 K only. `r10_temperature.py` (5x5x3, same events, |Phi|^2 and Gaussian weights,
Bose factors at T): kappa ratios at 100 / 200 / 500 / 800 K:
ALAMODE 0.980 / 0.996 / 1.005 / 1.006 (xx), 1.031 / 1.036 / 1.032 / 1.031 (zz);
k = 1 max x: **1.234** / 1.057 / 0.980 / 0.969 (xx); k = 1 max z: **1.826** / 1.091 / 0.980 / 0.972 (zz);
k = 2 max x: **1.662** (xx at 100 K); k = 2 max z: **3.169** (zz at 100 K, tau_z x 22);
new C1 max z: 1.302 (zz at 100 K). A kappa(T) curve excludes every witness; ALAMODE stays within 3.6 %.

### 4.4 Cutoff-induced part of the class
Same optimiser, k = 1 bounds except on the 266 ALAMODE-zero orbits, which are held at the 0.5 % floor
(`r2_scaling.py` variant `cal1nz`, cond(M) = 5.4e8, isometry 1.9e-10; runs `c2_nz_*`):
max tau_x **63.66 ps (1.858)** instead of 80.33 (2.345); max tau_z **49.93 ps (1.820)** instead of 59.04 (2.152);
PENDING minima. Verified by phono3py pp + independent operator (63.6571 / 49.9377 ps, lifetimes within 9.980 %,
kappa within 1 %). About 40 % of the excess over the reference comes from orbits on which no inter-workflow
difference exists; a factor ~1.8 remains on the orbits with a measured difference (5x5x3 statement).

### 4.5 Does the ALAMODE data violation undermine or strengthen C2?
ALAMODE vs reference: lifetimes -18 % ... +18 % (5x5x3; -27 % ... +17 % at 7x7x5), kappa_zz +3.3 % (+1.4 %). The data
tolerance of the class (+-10 %, +-1 %) is therefore tighter than the spread of the data between the two
workflows. Formally this strengthens the inner-bound reading (a looser, realistic tolerance can only widen the
range). It weakens the reading of the k = 1 box as "realistic precision": the one realistic alternative tensor
moves tau by < 2 % on both meshes and at all temperatures, while the box's adversarial corners move it by 2.3x on
5x5x3 and in the opposite direction on 7x7x5. The box (any direction inside every orbit, any combination across
orbits, centred on the reference, 45 % of its radius from a cutoff) is a worst-case construction calibrated by a
single sample; it is not evidence about plausible fc3 uncertainty.

### 4.6 Status of the interpretation
"Static data do not pin tau_mem to better than a factor ~2 at realistic fc3 precision": **not supported as a
physical statement** — the only evidence offered are witnesses that are artefacts of the 5x5x3 sampling (4.2, the
witnesses as physical counterexamples are REFUTED), are data-consistent only in one finite-sigma convention (2.2),
and are excluded by kappa(T) (4.3). This review does not show the opposite (that the data pin tau_mem); the
physical question is UNRESOLVED. The statement **survives as a property of the 5x5x3 event operator** (inner bound). "Memory carried by small long-range constants": correct
as anatomy of these 5x5x3 witnesses (about 43 % of their squared fc3 change lies on orbits beyond ALAMODE's cutoff;
restricting those orbits removes about 40 % of the excess), but not established as physics. Whether tensors exist
that keep the data on both meshes (and in both conventions) and still move tau by a factor ~2 is **UNRESOLVED**.

---

## 5. Block versus orbit average and other convention checks
* Block (optimiser) vs orbit (phono3py route) averaging: kappa differences <= 1.4e-4 for all witnesses checked
  (e.g. new max z witness kappa_zz 1.00980 block vs 1.00994 orbit); tau differences <= 0.03 %. No error.
* phono3py-style degeneracy averaging of the event operator: kappa ratios unchanged to <= 2e-4 (`r6b`). Raw
  lifetimes of degenerate partners in the orbit average differ by up to -17 % / +28 % (gauge effect, item 2 of 09);
  degeneracy-averaged lifetimes stay within 9.98 %.
* The convention that matters is the detuned-event weighting (2.2): 1.2 - 3.2 percentage points in kappa.

## 6. Attacks that failed (claims survived)
* Polarity, image ordering, image multiplicity, units, symmetry and sum rule of the folded tensor (1.1).
* Selection-rule events in the witnesses: max |Phi|^2 on the 3,275 zero-vertex events <= 5.4e-15 of the maximum for
  every witness (they stay zero).
* Permutation symmetry of phono3py's |Phi|^2 on the event list: <= 6e-5 (REF and witnesses).
* Rebuilt witness tensors: re-projection residual 6e-16 (symmetric), sum rule 1e-15, per-orbit bounds satisfied when
  measured on the tensor itself (independent of W and M^-1/2).
* Stoichiometric rescaling alone ('raw' vs 'cons'): -0.3 points in the kappa ratio, opposite sign to the phono3py
  shift; not the mechanism.

## 7. Not performed / open
* Re-optimisation with the data imposed on both meshes, or in phono3py's convention: not run (the 7x7x5
  linearisation exists on D: but a joint model is about 3 GB and ~1.5 s per evaluation). Existence of mesh- and
  convention-robust witnesses with large tau: UNRESOLVED.
* Global optimality: no certificate for any maximum (campaign's or this review's).
* Second calibration sample (ShengBTE) not mapped; delta_o remains a single-sample spread.
* Behaviour at smaller sigma or converged meshes: not tested; the convention effect of 2.2 is O(Delta/kT) and is
  expected to shrink with sigma.

## 8. Files
Scripts (this folder): `r1_orbit_anatomy.py`, `r2_opt.py`, `r2_scaling.py`, `r3_mapping.py`, `r4_p3_lbte.py`,
`r4_p3_lbte_mesh.py`, `r4b_p3_variants.py`, `r5_event_operator.py`, `r6_save_omega.py`, `r6_compare_operators.py`,
`r6b_degeneracy_convention.py`, `r6c_compare_avg.py`, `r6d_detuning_convention.py`, `r7_pp_permutation.py`,
`r8_mesh_transfer.py`, `r9_build.py`, `r10_temperature.py`, `r11_anatomy.py`, `r12_fc2_mapping.py`.
Outputs: `r1_orbit_anatomy.json/.npz`, `r1_witness_anatomy.json`, `work/r3_mapping.json`, `work/r3_orbit_diameters.npz`,
`work/r12_fc2_mapping.json`, `r4/` (phono3py LBTE 5x5x3), `r4m775/` (phono3py LBTE 7x7x5), `r5/` (independent
operator, convention tests), `r7/` (permutation test), `m553v/` (new witnesses, 5x5x3 phono3py pp), `m775/`
(7x7x5 phono3py pp and transfer tables), `r10/` (temperature transfer), `runs/` (optimiser logs, json, `_y.npy`),
`work/fc3_*.npy` (rebuilt tensors), `work/scaled_cal1nz_*` (modified scaling).
