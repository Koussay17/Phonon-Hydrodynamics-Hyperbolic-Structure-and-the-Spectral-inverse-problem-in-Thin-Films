# Peer review — final report of 11 October 2026 (`13-final-report.md`)

**Date:** 11 October 2026. **Object:** `13-final-report.md` (version of 11 October 2026), read together with
`00-question.md`, `02-assumptions.md`, `07-experiments.md` (E1–E22), `09-failed-approaches.md` (items 14–17),
`review/counterexample-review.md` and `superseded_docs/13-final-report_2026-10-10_withdrawn.md`.
**Mandate:** hostile journal-style review: test every conclusion in sections 1, 2 and 4 against its evidence, and
judge the fairness of the mesh-transfer methodology, the meaning of "no robust change above ~20 %", the conclusion
drawn for paper 1, and the missing controls and literature.

**What was recomputed.** The event operator was re-implemented from the written definitions (02-assumptions A–B,
`memgap.py` docstring) without importing any campaign module (`review/peer-review-checks/referee_check.py`), using
raw per-event rates (no orbit averaging). It was evaluated on the campaign's stored phono3py interaction strengths
(C: mirror `results/m553`, `results/m775x`, `results/m996`). Additional diagnostics are infrared structure, spectral
location of each witness effect, a diagonal/off-diagonal decomposition and the rigorous floor per mesh. Scripts and
raw outputs are listed in section 11. No new phono3py run was made, and peak memory stayed below 2 GB. No campaign
file outside `review/` was modified.

---

## 1. Summary of the claims under review
1. **Headline (section 1).** No fc3 change "within first-principles precision" has been found that moves tau_mem
   by more than about 20 % in a way that survives a change of Brillouin-zone sampling while keeping the static
   data. All large changes found (factors 1.3–6) are "artefacts of the sampling on which they were optimised".
2. **Identifiability (S7).** Whether the data pin tau_mem is unresolved.
3. **Discretisation (S8).** tau_mem is not mesh-converged. Its "discretisation error" (+79 % / +116 %, 5x5x3 to
   9x9x6) "exceeds every robust fc3 effect found".
4. **Supporting statements.** S1–S4 cover reproduction, controllability, rigorous bounds and single-mesh ranges.
   S5 says single-mesh witnesses lose or reverse their effect on a finer mesh. S6 says the largest change on an
   unused mesh is 1.18x. S9–S12 cover calibration, local identifiability, sensitivity and convention dependence.
5. **Paper 1 (section 4).** The free-rate gap is "an artefact of ignoring the amplitude structure". Single-mesh fc3
   ambiguities are sampling artefacts. "tau_mem from first principles is plausible but must be quoted with a
   convergence study in mesh and smearing", followed by identifiability "on the converged operator".

## 2. Verification performed (numbers reproduce)
| quantity | report | this review (independent operator, raw rates) |
|---|---|---|
| REF 5x5x3 kappa_xx / kappa_zz; tau_x / tau_z | 214.300 / 201.967; 34.2530 / 27.4465 | 214.300 / 201.968; 34.2530 / 27.4470 |
| REF 7x7x5 | 243.027 / 272.334; 52.609 / 50.623 | 243.027 / 272.345; 52.609 / 50.622 |
| REF 9x9x6 | 265.76 / 256.90; 61.40 / 59.15 | 265.763 / 256.965; 61.397 / 59.146 |
| E19 (8 tensors on 7x7x5) | `results/m775x/mesh_transfer.json` | all ratios to <= 1e-5 |
| E20 (4 joint witnesses, both meshes) | `verify_m553.json`, `verify_m775.json` | all ratios to <= 1e-4 |
| E21 (9x9x6, 6 tensors) | `verify_m996.json`, `transfer_controls.json` | all ratios to <= 1e-4 (e.g. 1.1829, 0.9729, 1.1819, 1.1438; ALAMODE 0.9783 / 0.9324) |

The reference differences (<= 2.5e-4 relative) are the documented orbit-average versus raw-rate convention
(09 item 2). The 7x7x5 ALAMODE ratios were confirmed from `review/m775/r8_mesh_transfer_cons.json`.

One labelling error was found. The E19 row "orb 0.1 max x (1.308)" is the E15 zero-start witness: its 5x5x3 ratio
is 1.287 (44.08 ps; `results/m553/witness/verify_orb0.1_c0_max.json`), not the E16 multistart best (44.80 ps).

**The numerical content is reproducible.** The objections below concern interpretation, test design and physics,
not arithmetic.

## 3. Strongest parts
* The 10 October conclusion was withdrawn after an adversarial review, and the withdrawal is documented (09 items
  14–17). This is the correct scientific response.
* The engineering is careful:
  * hash-pinned inputs;
  * an amplitude engine independent of phono3py's Rust kernels (agreement 6.5e-12);
  * exact linearisations that reproduce phono3py-route ratios to <= 5e-5;
  * agreement of Cholesky, eigendecomposition and conjugate gradients.

  This review reproduced every number it checked.
* The ALAMODE mapping was validated by two maps related by a crystal symmetry. This caught a sign error that the
  symmetry residuals and the sum rule could not detect (09 item 12), which is a model practice.
* Every optimisation result is consistently labelled as an inner bound, and S7 is kept UNRESOLVED.
* Validating on a mesh unused in the optimisation (E21) is the right instinct. The second-workflow control at three
  meshes is the most informative physical datum in the campaign: tau moves by <= 7 % while single lifetimes move by
  up to 28 %.

---

## 4. Potentially fatal issue

### F1. A finite continuum limit of tau_mem is not established, and the campaign's own operators are consistent with an infrared divergence
**Claims attacked.**
* S8 and section 1: "its discretisation error (+79 % / +116 % from 5x5x3 to 9x9x6)".
* Section 4: "must be quoted with a convergence study in mesh and smearing".
* Section 4: "The natural next steps are the convergence study of tau_mem itself, and only then the identifiability
  question on the converged operator".

**Evidence** (reference fc3; `peer-review-checks/ref_m*.json`, `referee_ir.json`, `referee_floor.json`).
1. **tau_mem is carried by the lowest-frequency acoustic modes.**
   * On 7x7x5 and 9x9x6, the 2 % lowest-frequency modes carry 66–77 % of N = |C^+ b|^2 but only 26–36 % of K.
   * Four modes alone carry 36 % (7x7x5) and 35 % (9x9x6) of N_z. These are the doubly degenerate transverse
     acoustic modes at the smallest c-axis wavevector, +-q_c (f_min = 2.30 and 1.93 THz), and they carry no N_x.
2. **The diagonal rates follow Herring's low-frequency law.** r = diag C scales as f^1.95 on 7x7x5 (2.3–4.3 THz,
   46 modes) and f^2.06 on 9x9x6 (1.9–4.9 THz, 136 modes).
3. **Asymptotic argument.**
   * For these modes the operator is RTA-like: x = C^+ b satisfies x_mu / (b_mu / r_mu) = 0.96–1.21.
   * In entropy coordinates b_mu = d_mu eps_mu v_mu tends to k_B T v_mu as eps_mu goes to 0. The per-mode
     contribution to N therefore scales as omega^-4, while the 3D density of modes scales as omega^2.
   * Hence N ~ ∫ omega^-2 d omega ~ 1/omega_min diverges, while K ~ ∫ d omega converges.
   * This is the phonon form of a long-time tail. In 3D the heat-current autocorrelation decays as t^-3/2, so its
     integral (K) is finite and its first moment (N) is not. Equivalently, kappa(Omega) - kappa(0) ~ (i Omega)^(1/2),
     so tau_mem = -K'(0)/K(0) does not exist.
4. **The three meshes follow this law, not a convergent one.**
   * tau_x = 4.6 + 110/f_min[THz] ps reproduces all three values to 0.5 % (pairwise slopes 112 and 105 ps THz).
   * An exact three-point fit tau = tau_inf - c f_min^p gives p = -0.80 for x (divergent). For z it gives
     p = +0.03 with tau_inf = 1.6 ns, which is no usable limit.
   * tau/tau_CS grows 1.78 -> 2.49 -> 2.74 (x) and 1.77 -> 2.53 -> 3.23 (z), while tau_CS = K/|b|^2 changes by only
     17–18 %. The growth is in N, not K.

**Status.** PLAUSIBLE, not established (three non-nested meshes at fixed sigma, plus standard asymptotics). The
mechanism is textbook: Herring, Phys. Rev. 95, 954 (1954) for low-energy phonons; Pomeau & Résibois, Phys. Rep. 19,
63 (1975), doi:10.1016/0370-1573(75)90019-8 for long-time tails and divergent higher transport coefficients. The
omega^(1/2) and Lévy phenomenology appears in Koh & Cahill, PRB 76, 075207 (2007) and in Vermeersch, Carrete, Mingo
& Shakouri, PRB 91, 085202 (2015). The argument uses only r ~ omega^2 and the 3D mode density. It does not obviously
require sigma -> 0, so refinement at fixed sigma may diverge as well.

**Consequences if confirmed.**
* The +79 % / +116 % is the cutoff dependence of a divergent quantity, not a discretisation error. No convergence
  study converges bulk tau_mem.
* "Identifiability on the converged operator" is ill-posed for the absolute tau_mem. Well-posed alternatives are:
  * the tau ratio, which would tend to the ratio of infrared coefficients;
  * the infrared coefficient itself;
  * tau_mem with a physical regulator, the film thickness L of the thesis setting. A crude RTA estimate with a
    boundary rate ~ v/L gives N ~ L^(1/2) (heuristic).
* A Cattaneo-type closure with a finite memory time is not the low-frequency limit of the bulk three-phonon BTE.
  This bears directly on paper 1.
* The mesh-transfer test compares quantities dominated by different modes on each mesh (M2a).

**Fix.** Establish the infrared behaviour before any further identifiability work (test T1), then rewrite S8,
section 1 and section 4 accordingly.

---

## 5. Major issues

### M1. S5, S6 and "all large changes … are artefacts" are contradicted by the campaign's own transfer data
**Claims attacked.**
* Section 1: "All large changes found (factors 1.3–6) are artefacts of the sampling on which they were optimised."
* S5: "Every single-mesh witness tested on a finer mesh loses or reverses its memory change".
* S6: "On a mesh unused in the optimisation, the largest tau change of any witness is 1.18 x".

**Why they fail.**
* **E19's own table.** The orb 0.1 max z witness (5x5x3 ratio 1.370) gives tau_z 1.347 on 7x7x5, a mesh unused in
  its optimisation. It keeps 94 % of its excess.
* **The review's transfers** (`review/m775/r8_mesh_transfer_newwitness_cons.json`; review 3.3 and 4.2):
  * per-cluster max z continuation witness: 1.436 -> 1.315;
  * cal1nz single-mesh max x witness: tau_z 1.330 on 7x7x5, while its objective axis reverses to 0.872;
  * cal1nz single-mesh max z witness: 1.820 -> 1.184.
* The two cal1nz cases belong to the calibrated class that the headline equates with "first-principles precision".
* The section-1 summary of E19 quotes only the calibrated rows (0.83 / 1.08 / 0.71 / 0.53). It omits the rows that
  contradict the narrative.
* The 9x9x6 transfer of the one witness that survived the first transfer is missing (E21 note).

**Fix.**
* Restate S5 as follows. The calibrated maxima lose or reverse their effect on 7x7x5. The per-cluster c-axis maxima
  keep 72–94 % of their tau_z excess (1.370 -> 1.347 and 1.436 -> 1.315). All violate the declared tolerance on the
  new mesh.
* Correct S6 to the actual maxima on unused meshes: 1.35 on 7x7x5, 1.18 on 9x9x6.
* Restore the omitted rows in section 1.
* Complete the missing 9x9x6 transfers (T3).
* Drop the word "artefact" unless a test able to establish it is passed (M2).

### M2. The mesh-transfer test cannot discriminate "sampling artefact" from "physical ambiguity"
**Claims attacked.**
* S5 (status NUMERICALLY DEMONSTRATED).
* E21's conclusion "Two-mesh constraints do not produce sampling-robust witnesses".
* Section 4: "they disappear on a mesh not used to construct them, where they also break the data".

**Why it is not a fair test of the stronger claim.**
* **(a) The objective is mesh-specific.**
  * N is dominated by the lowest-|q| modes of each mesh (F1). These are different physical modes on each mesh:
    f_min is 3.70, 2.30 and 1.93 THz. The meshes are not nested: 5x5x3 and 9x9x6 share only three Gamma-line
    points, and 7x7x5 shares only Gamma with either.
  * Every large witness effect sits in exactly these modes (`referee_anatomy.json`):
    * cal1 max x on 5x5x3 gains dN/N_ref = +1.54 from the 3–5 THz band, the lowest band present;
    * cal1nz joint max z on 7x7x5 gains +1.77 from its four lowest modes;
    * orb 0.1 joint max z on 7x7x5 gains +0.70 from its four lowest modes.
  * A tensor tuned to the lowest modes of mesh A is then judged on the lowest modes of mesh B, which it never saw.
    Loss of the effect is the generic expectation whatever the physics.
* **(b) Boundary witnesses against ratio noise larger than the tolerance.**
  * The witnesses saturate the constraints by construction (lifetimes 0.0998, kappa 1.0098 / 0.9902).
  * For a benign alternative tensor, the same ratios move between meshes by more than the tolerance. ALAMODE's
    kappa_zz ratio is 1.033 / 1.014 / 0.994 on 5x5x3 / 7x7x5 / 9x9x6, and its tau_z ratio is 1.017 / 1.006 /
    0.932.
  * A point on the boundary of a +-1 % box leaves it under shifts of this size even if it is "physically
    equivalent". No optimisation margin was used and no allowance was made for ratio noise.
  * The per-mesh relative tolerance is the only coherent choice, given that the reference kappa itself moves by
    24 % (214 to 266 W/(m K)). But it makes the constraint set itself mesh-dependent, and the test must account for
    that.
* **(c) The criterion rejects the real second workflow.**
  * ALAMODE violates the +-10 % lifetime tolerance on every mesh: maximum 18 %, 27 % and 28 %, with 2.4 %, 2.3 % and
    1.5 % of degenerate groups outside +-10 %. It also violates kappa_zz on two of three meshes.
  * By the criterion of S5/S6, the second first-principles tensor would be an "artefact".
  * On 9x9x6, the E20 orb 0.1 joint max x witness (tau_x 1.182) keeps kappa within 0.6 %. Its largest lifetime
    deviation is 21 % (2.7 % of groups outside +-10 %). It is about as data-consistent there as ALAMODE (28 %,
    1.5 %).
* **(d) No positive or null control.** The power of the test, i.e. its probability of retaining a genuinely robust
  ambiguity, was never measured. E5 shows that untuned perturbations lower tau and kappa on average. The null
  expectation for a transferred ratio is therefore below 1, not 1. On 9x9x6, values of 1.18 are above the null.
* **(e) The three meshes are not a convergence sequence.**
  * The reference is non-monotone: kappa_zz is 202 -> 272 -> 257.
  * Witness ratios are non-monotone: cal1 max x is 2.345 -> 0.830 -> 0.985.
  * All three meshes use a sigma far below the frequency spacing of neighbouring q points (M7).
  * "Finer" does not mean "closer to the limit" here.

**What E19–E21 do establish.** They falsify "the 5x5x3 witnesses are properties of AlN", which is what the withdrawal
needed. They cannot support "all large changes are artefacts" or any robust bound.

**Fix.** Add controls (T4). Optimise with margins and validate with tolerances widened by the benign-tensor ratio
noise (T5). Downgrade S5 to "mesh-specific", not "artefact".

### M3. "No robust change above ~20 %" has no defined meaning as stated
**Claims attacked.**
* Section 1: "No fc3 change within first-principles precision has been found that moves tau_mem by more than about
  20 % …".
* S6.
* Section 4: "no demonstrated fc3-driven ambiguity of tau_mem beyond about 20 %".

**Why it fails.**
* **(a) The admissible set on the validation mesh is empty apart from the reference.** Every tested tensor violates
  the declared data there. "The largest tau change of any witness" is therefore a maximum over inadmissible
  tensors. It is neither a lower bound on the robust ambiguity (the tensors are not admissible) nor an upper bound
  (the search is not exhaustive).
* **(b) No search ever targeted the property.**
  * E20 maximised tau on 7x7x5 only; a robust objective was never optimised.
  * The ~18 % is an incidental by-product.
  * The per-cluster joint max-z witness keeps the data on both 5x5x3 and 7x7x5 with tau_z 1.18 and 1.69. The
    "about 20 %" threshold sits at the incidental 5x5x3 value of a run that did not aim at it.
* **(c) Budget.**
  * There were four joint runs, each with a single start (u = 0) and 15 outer x 300 inner iterations.
  * All 60 inner L-BFGS-B solves stopped at maxiter (`results/joint/*.json`), while the penalty rose to
    8e5–3e6. The tau plateau reflects stiffening, not stationarity (review 3.1).
  * E22 shows that continuation moves the 5x5x3 maxima by 5–7 %, and its 47.87 ps polish iterate still had
    nit = 1000 at every outer step.
* **(d) Only increases were searched.** No joint minimisation was run. On the meshes in question, the rigorous floor
  (1 - dK) tau_CS allows ratios down to 0.40 / 0.39 (7x7x5, x / z) and 0.36 / 0.31 (9x9x6) (`referee_floor.json`). The
  transferred 5x5x3 minima already reach 0.54–0.72 on 7x7x5 (E19). The lower half of "+-20 %" is unexplored.
* **(e) Dependence on the tolerance.** Section 3 item 3 concedes that the declared tolerance is stricter than
  first-principles precision. Under a tolerance matched to the observed inter-workflow spread (kappa +-3.3 %,
  lifetimes +-28 %):
  * the E20 orb 0.1 max x witness is admissible on all three meshes (tau_x 1.057 / 1.326 / 1.182);
  * the E15 orb 0.1 max z witness is admissible on two (1.370 / 1.347).

  These per-cluster tensors change the dominant on-site and < 2 A clusters by 2.9–5.3 % (anatomy files). That lies
  outside the single-pair calibration but not demonstrably outside first-principles precision (M6).
* **(f) Resolution floor.** The framework's own noise is 5–10 %:
  * inter-convention differences of witness tau ratios, 2.429 vs 2.345 and 4.066 vs 3.721 (review 2.2);
  * benign-tensor mesh noise of tau ratios, up to 8.5 points.

  "About 20 %" is only 2–4 times this floor.

**Fix.** Replace the headline by an absence-of-evidence statement without a number. For example: "No tensor was
found that satisfies the declared data on a mesh not used in its construction; on such meshes the tested witnesses
give tau ratios 0.53–1.35 (7x7x5) and 0.96–1.18 (9x9x6)." A quantitative robust statement requires T3–T6.

### M4. Two mechanisms are conflated: the report underclaims what E19–E21 show and misses a route to a real result
**Claims attacked.**
* S5, S6 and section 4 treat all mesh-dependence as one "artefact".
* E20's "Reading" ("the optimiser exploits whichever triplets the objective mesh samples").

**Referee analysis** (`referee_drag.json`, `referee_drag2.json`; NUMERICAL OBSERVATION).
* **Method.** Over the modes that carry 50 % of N_ref, log(x_w/x_ref) splits into two parts:
  * a diagonal part, -log(r_w/r_ref), which is the RTA expectation;
  * an off-diagonal remainder, the in-scattering or "drag" part.
* **On objective meshes the enhancement is mostly drag.** Lifetimes are clamped to +-10 % there.

  | witness | mesh | log(x_w/x_ref) | diagonal | drag |
  |---|---|---|---|---|
  | cal1 max x | 5x5x3 | 0.675 | 0.050 | 0.625 |
  | cal1nz joint max z | 7x7x5 | 0.706 | 0.062 | 0.644 |
  | orb 0.1 joint max z | 7x7x5 | 0.428 | 0.056 | 0.373 |

  For the lowest modes, x/(b/r) rises from about 1.0 to 1.6–3.2.
* **The effect that transfers intact is diagonal.**
  * E15 orb 0.1 max z on 7x7x5 (1.347): 0.265 = 0.226 diagonal + 0.038 drag. The lowest TA rates are lowered by
    23 %.
  * On 9x9x6, the surviving 1.144 of the orb 0.1 joint max z witness is also diagonal (0.072 - 0.008).

**Interpretation (HEURISTIC).**
* In entropy coordinates, the off-diagonal couplings of a mode of energy eps are of relative order eps/k_B T with
  respect to its source term b, while its diagonal rate scales as eps^2.
* On these meshes eps_min/k_B T = 0.31–0.59. Drag is therefore of order one and can be tuned at the few lowest mesh
  modes. This is a pre-asymptotic, mesh-specific lever.
* The diagonal lever, i.e. the low-frequency acoustic scattering prefactor, survives refinement. It is limited only
  by how tightly the data fix the low-frequency lifetimes.

**CONJECTURE worth pursuing.** Suppose that:
* (i) N is infrared-divergent and K convergent (F1);
* (ii) x_mu/(b_mu/r_mu) -> 1 uniformly in the fc3 class as eps_mu -> 0;
* (iii) the data fix every low-frequency lifetime to +-delta.

Then the continuum-limit tau ratio lies in [(1+delta)^-2/(1+dK), (1-delta)^-2/(1-dK)], which is [0.82, 1.25] for
delta = 10 % and dK = 1 %. This would be an asymptotic identifiability statement of exactly the kind the decision
question asks for. It might also be why the 9x9x6 values cluster near 1.2. With realistic data, which contain no
low-frequency linewidths, the infrared prefactor is not fixed and the answer can differ.

**Also underclaimed.** E20 is a NUMERICALLY DEMONSTRATED negative result for the discretised problem. On 7x7x5, with
the data on two meshes, the cal1nz class admits 1.80 (x) and 3.05 (z), against 1.86 / 1.82 on 5x5x3 alone. The gap
of each discretised operator does not shrink with refinement. The report presents this only as an artefact.

**Fix.** Report the two channels, state the E20 result as a property of each discretised operator, and test the
conjecture (T7).

### M5. The conclusion for paper 1 (section 4) does not follow
**Claims attacked and why each fails.**
* **(a)** "The free-rate memory gap is an artefact of ignoring the amplitude structure".
  * E10 shows that the amplitude structure alone, under a ~10 % Frobenius bound, still admits 14x (x) and 624x (z)
    on the same mesh.
  * What closes the gap on 5x5x3 is the per-cluster plausibility prior.
  * Neither the free-rate witnesses nor the Frobenius witnesses were ever transferred to another mesh. The roles of
    amplitude structure, prior and sampling are therefore not separated.

  UNSUPPORTED.
* **(b)** "a discretisation dependence of tau_mem … larger than any fc3 effect".
  * The statement is unqualified here, although section 1 qualifies it with "robust".
  * As written it is false: E10 gives 624x and E20 gives 3.05x on 7x7x5.
  * It also compares different quantities: a cross-mesh change of the absolute reference value against within-mesh
    ratios. A uniform scaling of fc3 changes the ratio identically on every mesh while the absolute value changes
    by 79 %.
* **(c)** "tau_mem from first principles is plausible". No evidence is offered beyond one same-functional second
  workflow at fixed fc2 on unconverged meshes. The quantity may be undefined in bulk (F1).
* **(d)** "must be quoted with a convergence study … only then the identifiability question on the converged
  operator". This presupposes a finite limit (F1).

**What does follow.**
* The 10 October conclusion is correctly withdrawn.
* The single-mesh ranges are properties of the discretised operators.
* Identifiability is unresolved.
* The second workflow changes tau by <= 7 % on all three meshes, which is a single-sample, fixed-fc2,
  same-functional observation.
* If F1 holds, paper 1 should not quote a bulk tau_mem. It should state the infrared structure, and pose the
  identifiability question for the tau ratio, for the infrared coefficient or for a film-regularised tau_mem.

### M6. "First-principles precision" rests on one workflow pair at a fixed functional, and the data tolerances contradict documented spreads
**Claims attacked.**
* Section 1: "within first-principles precision".
* S9.
* Section 0: "calibrated per-cluster bound eps_o = max(k delta_o, 0.5 %)".

**Why it is under-supported.**
* **Same functional.** The Phonon Olympics benchmark is the source of both tensors (McGaughey et al., J. Appl. Phys.
  (2025), doi:10.1063/5.0289819). It prescribed PBEsol PAW for every team, AlN included.
* **What delta_o therefore measures.** It captures workflow differences (supercell, displacements, fitting, cutoffs)
  at fixed functional. Three sources of first-principles spread are absent:
  * functional dependence, a recognised source of spread in first-principles kappa (Jain & McGaughey, Comput. Mater.
    Sci. (2015), doi:10.1016/j.commatsci.2015.08.014; on force-constant accuracy, Xie, Gu & Bao, Comput. Mater. Sci.
    (2017), doi:10.1016/j.commatsci.2017.06.047);
  * supercell-size dependence of the folded long-range constants;
  * fc2 uncertainty (section 3 item 4).

  delta_o is thus a lower bound on precision for the dominant clusters, where it equals 0.3–0.6 %.
* **Tolerances versus documented spreads.**
  * The same benchmark reports 300 K AlN conductivities of roughly 250–300 W/(m K) in-plane and 230–290 W/(m K)
    cross-plane across packages and solution methods. It attributes the variation to the three-phonon lifetimes
    obtained from the cubic constants.
  * The declared kappa tolerance (+-1 %) is an order of magnitude tighter than that spread.
  * The lifetime tolerance (+-10 %) is tighter than the campaign's own inter-workflow spread (18–28 %).
* **Consequence.** The headline class is ill-defined. Whether the per-cluster eps = 0.1 class (dominant clusters
  changed by 3–5 %) lies within first-principles precision is undecided, not "no", and the robust conclusions
  depend on it (M3e). The headline also does not say which class (cal1, cal1nz, orb) it refers to, although S6
  pools them.

**Fix.** Add a third workflow and a second functional, state which class the headline means, and report every
conclusion as a function of tolerance and precision model (T6).

### M7. No controls on Brillouin-zone integration, and missing literature positioning
**Claims attacked.**
* Section 0 and section 3 item 2: "Three meshes, one smearing".
* S8.
* The use of 7x7x5 and 9x9x6 as validation samplings.

**Why it is a major gap.**
* **sigma is far below the frequency spacing.**
  * sigma is fixed at 0.1 THz (4 sigma cutoff) on meshes of 75–486 q points.
  * Neighbouring acoustic frequencies differ by v Dq/2pi, with in-plane Dq = 0.26–0.47 A^-1 and v = 6–11 km/s.
    That is about 2.5–8 THz, i.e. 25–80 sigma.
  * Energy conservation is therefore sampled sparsely and erratically. This is the mechanism behind the
    non-monotone kappa_zz and behind witnesses that tune "whichever triplets the mesh samples".
* **Standard practice was not followed.** The usual remedies are:
  * adaptive Gaussian smearing (Li, Lindsay, Broido, Stewart & Mingo, PRB 86, 174307 (2012); ShengBTE: Li, Carrete,
    Katcho & Mingo, Comput. Phys. Commun. (2014), doi:10.1016/j.cpc.2014.02.015);
  * the tetrahedron method (Togo, Chaput & Tanaka, PRB 91, 094306 (2015)), with sigma reduced together with Dq.

  No sigma scan and no tetrahedron cross-check exist for kappa, for tau or for any witness ratio.
* **The meshes are far from converged.**
  * The benchmark's AlN calculations used q-grids from 31x31x17 (phono3py) to 53x53x33 (ALAMODE), 30–190 times
    more q points than 9x9x6.
  * The campaign's 9x9x6 kappa (266 / 257) lies 5–7 % below the benchmark's full-solution values (285 / 271) for the
    phono3py dataset from which the campaign's inputs are taken. The 5x5x3 values lie about 25 % below.
  * These figures belong in section 0 to give "coarse" a quantitative meaning.
* **Literature absent from the report.**
  * **Framework.** The relaxon formulation of the conserving BTE (Cepellotti & Marzari, PRX 6, 041013 (2016)), in
    which tau_mem = N/K is the kappa-weighted mean relaxon lifetime. Viscous heat equations (Simoncelli, Marzari &
    Cepellotti, PRX 10, 011019 (2020)).
  * **Infrared physics.** Herring (1954); Pomeau & Résibois (1975); Koh & Cahill (2007); Vermeersch et al. (2015).
  * **Identifiability methodology.** Profile-likelihood structural and practical identifiability (Raue et al.,
    Bioinformatics 25, 1923 (2009)); the campaign's "max/min tau subject to the data" is exactly a set-membership
    profile. Sloppy-model analysis of the data-Jacobian spectrum (Gutenkunst et al., PLoS Comput. Biol. 3, e189
    (2007)), which E9/E14 already compute without naming it.
  * **Inverse problems for phonon relaxation times.** Minnich, PRL 109, 205901 (2012); Forghani, Hadjiconstantinou &
    Péraud, PRB 94, 155439 (2016).
  * The benchmark paper itself, which the report uses only as a data source.

**Fix.** Run T2, add the convergence figures to section 0, and position the work against this literature.

### M8. The data class is synthetic, mesh-relative, and not the data a paper about AlN can use
**Claims attacked.**
* The decision criterion of `00-question.md` and section 0: "Data: all mode-orbit lifetimes … within +-10 %, kappa
  … within +-1 %".
* Section 4's reading for "a conserving closure of AlN".

**Why it is under-supported for the paper-1 question.**
* **What the "data" are.**
  * "Lifetimes" are r = diag C in entropy coordinates (the memory-gap convention), not measured linewidths.
  * The "data" are the reference tensor's own predictions on the same unconverged mesh.
  * The question answered is the self-consistency of discretised models, not determination from measurements.
* **Richer than measurements, tighter than precision.**
  * The class contains every mode lifetime, including the near-Gamma acoustic modes that control tau_mem (F1, M4).
  * Real data are kappa(T), Gamma-point Raman linewidths, a few inelastic-scattering linewidths and thermal
    expansion. With such data the infrared prefactor is not measured directly. Its link to measurable Grüneisen
    parameters and thermal expansion would have to be established.
  * The class is also tighter than first-principles precision (M6).
* **Worst case versus probability.** Worst-case inner bounds over a box can exceed a probabilistic uncertainty by up
  to about sqrt(dimension) for isotropic priors; compare E5 with E10. For paper 1, a data-conditioned probabilistic
  propagation would answer the practical question ("what uncertainty should accompany a quoted tau_mem?") more
  directly than adversarial corners. A linear-Gaussian posterior of log tau, using the existing exact Jacobians and
  a prior calibrated on at least three workflows and two functionals, would do this.

**Fix.** State the data class plainly as synthetic. Add a realistic-data variant (kappa(T) and Gamma-point
linewidths). Add the probabilistic propagation as a complement to the inner bounds.

---

## 6. Minor issues
1. **E19 table, row "orb 0.1 max x (1.308)".** The transferred tensor has a 5x5x3 ratio of 1.287 (E15 zero start),
   not 1.308 (E16 best). *Fix:* correct the label, or transfer the E16 tensor.
2. **"Discretisation error" (section 1).** The change from 5x5x3 to 9x9x6 is a difference between two unconverged
   values, and possibly a cutoff dependence (F1). It is not an error estimate. *Fix:* say "change between meshes".
3. **S10.**
   * Full column rank of a 287 x 218 Jacobian is the generic outcome.
   * The smallest singular value, 0.118, means that a 1 % change in the data allows up to an 8.5 % change in the
     short-range fc3. Practical identifiability is weak.
   * The statement is first-order, at one mesh. The lesson of 09 item 14 applies to S10 as it does to S4.

   *Fix:* qualify the statement and repeat it on 9x9x6.
4. **S12 and section 3 item 5.** The phono3py-LBTE convention test and the kappa(T) exclusion were applied only to
   the 5x5x3 witnesses. The E20 joint witnesses' data consistency is untested in phono3py's convention and in
   temperature. *Fix:* apply both checks to E20.
5. **E22's 47.87 ps.** This is the polish run's best feasible iterate. It sits at the kappa boundary (1.00997), and
   every inner solve stopped at maxiter. The value phono3py verified is 47.514 ps (review 3.2). *Fix:* label 47.87
   as a model value, or verify it with phono3py amplitudes.
6. **Lower end of the range.** Section 1 discusses the minima only on 5x5x3. The transferred minima (0.54–0.72 on
   7x7x5, E19) are not mentioned. *Fix:* add them.
7. **Untransferred witnesses.** The Frobenius-class witnesses (E10, 14x / 624x) and the 3.5 A subspace witnesses
   (E11, up to 2.8x / 3.3x) were never transferred. The latter are smooth and short-range, so they are a natural
   control for M2. *Fix:* transfer them (T4).
8. **Asymmetric presentation of E21 in section 1.** ALAMODE is reported with "kappa within 0.6 %", but not with its
   9x9x6 lifetime deviations (-28 % / +19 %). These are larger than those of the orb 0.1 joint max x witness.
   *Fix:* report both quantities for both tensors.
9. **No uncertainty on any tau ratio.** *Fix:* attach error bars that combine the convention effects (orbit versus
   block averaging, phono3py detuned-event weighting) and the benign-tensor mesh noise.
10. **Reproduction.**
    * The commands depend on about 15 GB of scaled arrays on C:, with hard-coded paths (`joint_mesh.py`:
      RES/MESH).
    * The witness tensors (4.5 MB each) and pp files are the evidence, but they carry no SHA-256 record, unlike the
      inputs.

    *Fix:* archive the tensors and pp files with hashes.
11. **Shell labels.** The anatomy tables (E17, E18, E20) use shell labels that are not cluster diameters for 60
    orbits (review 1.2). *Fix:* add a footnote, or re-bin by consistent diameter.
12. **Unnamed class.** The headline's "within first-principles precision" does not name its class. *Fix:* name it
    (see M6).

## 7. Status of each statement (referee assessment)
| statement | report status | assessment |
|---|---|---|
| Headline: no robust change > ~20 % | asserted | **UNSUPPORTED** as a quantitative statement (M3) |
| Headline: all large changes are artefacts | asserted | **CONTRADICTED** for the per-cluster c-axis witnesses (M1). For the calibrated maxima, mesh-specificity is NUMERICALLY DEMONSTRATED but "artefact" is not established (M2) |
| S1 | NUMERICALLY DEMONSTRATED | NUMERICALLY DEMONSTRATED, independently reproduced here |
| S2 | NUMERICALLY DEMONSTRATED | NUMERICALLY DEMONSTRATED (5x5x3, first order) |
| S3 | PROVED | accepted on the memory-gap audit (not re-audited here) |
| S4 | NUMERICALLY DEMONSTRATED (5x5x3) | NUMERICALLY DEMONSTRATED as budget-limited inner bounds |
| S5 | NUMERICALLY DEMONSTRATED | **CONTRADICTED as worded**; holds for the calibrated maxima only |
| S6 | NUMERICAL OBSERVATION | **CONTRADICTED as worded** (1.35 on 7x7x5). The 1.18 on 9x9x6 refers to inadmissible tensors |
| S7 | UNRESOLVED | UNRESOLVED; possibly ill-posed for the absolute bulk tau_mem (F1) |
| S8 | NUMERICALLY DEMONSTRATED | Non-convergence NUMERICALLY DEMONSTRATED. "Discretisation error" UNSUPPORTED. Infrared divergence PLAUSIBLE |
| S9 | NUMERICALLY DEMONSTRATED | Mapping NUMERICALLY DEMONSTRATED. The calibration excludes functional, supercell and fc2 uncertainty |
| S10 | NUMERICALLY DEMONSTRATED (first order) | First order, one mesh; generic and weakly conditioned |
| S11 | NUMERICAL OBSERVATIONS | accepted |
| S12 | NUMERICALLY DEMONSTRATED | NUMERICALLY DEMONSTRATED for 5x5x3; not applied to E20 |
| Section 4: free-rate gap is an artefact of ignoring the amplitude structure | asserted | **UNSUPPORTED** (E10; never transferred) |
| Section 4: discretisation dependence larger than any fc3 effect | asserted | **CONTRADICTED** as unqualified, and it compares different quantities |
| Section 4: tau_mem from first principles is plausible | asserted | **UNSUPPORTED**; possibly undefined in bulk |

## 8. Decisive tests
* **T1 — existence of the continuum limit** (cheap, memory-light).
  * Compute phono3py RTA linewidths of the acoustic branches on fine q-lines Gamma–A and Gamma–M, with the q2
    integration on >= 31x31x17, using both the tetrahedron method and adaptive smearing. Extract the exponent and
    prefactor of r(omega) per branch.
  * Compute RTA-level N and K on a Gamma-centred sequence (9x9x6, 12x12x8, 15x15x10, 18x18x12, 24x24x16), reducing
    sigma consistently.
  * Test N ~ 1/q_min against saturation, and K -> constant.
  * If N diverges, rewrite S8 and section 4, and reformulate the question as the tau ratio, the infrared coefficient
    or a film-regularised tau_mem.
* **T2 — integration controls at fixed mesh.** On 7x7x5 and 9x9x6, compute kappa, tau and the ratios for the
  reference and for the E20 witnesses with sigma in {0.05, 0.1, 0.2 THz, adaptive} and with the tetrahedron method.
* **T3 — missing transfers.** Run the 9x9x6 transfers of:
  * E15 orb 0.1 max z;
  * review `c1_z_cont02`;
  * the cal1nz single-mesh witnesses.

  The memory failure noted in E21 can probably be avoided with fewer threads or chunked triplets.
* **T4 — power and null distribution of the transfer test.** Transfer each of the following:
  * random perturbations with each witness's per-shell anatomy (>= 20 draws per witness);
  * the E11 short-range witnesses;
  * a role-reversed experiment that takes ALAMODE as the truth and phono3py as the alternative.
* **T5 — robust objective.**
  * Maximise and minimise min over meshes M of log(tau_M/tau_M,ref), or the infrared coefficient.
  * Impose the data on 5x5x3, 7x7x5 and 9x9x6 with margins (+-5 % / +-0.5 %).
  * Validate on a fourth mesh (RTA level is acceptable on 12x12x8), with tolerances widened by the benign-tensor
    ratio noise.
  * Use multistart, and report KKT residuals or run the inner solves to convergence.
* **T6 — tolerance and precision sweep.** Repeat the E20/E21 conclusions for:
  * lifetimes within +-10 / 20 / 30 % and kappa within +-1 / 3 / 5 %;
  * precision models built from at least three workflows and two functionals.
* **T7 — asymptotic conjecture (M4).** Track x_mu r_mu / b_mu of the lowest modes of the witnesses as f_min
  decreases along the T1 sequence. If it tends to 1 uniformly, attempt a proof of the [0.82, 1.25] bound under the
  idealised data class.
* **T8 — paper-1 setting.** Repeat the identifiability question with boundary scattering at the thesis film
  thicknesses, where tau_mem is finite.

## 9. Unresolved questions
1. Does bulk tau_mem (three-phonon, conserving operator) have a finite continuum limit? If not, which regulator is
   physical for paper 1, and how does tau_mem scale with it?
2. Does the drag lever on the lowest modes vanish as f_min -> 0, uniformly in the fc3 classes considered?
3. Are the per-cluster c-axis effects (1.35 on 7x7x5) robust on 9x9x6 and in the limit sigma -> 0?
4. Which static, measurable data constrain the low-frequency acoustic scattering prefactor: Grüneisen parameters,
   thermal expansion, or low-temperature kappa(T)?
5. What is the first-principles precision of fc3 once functional, supercell and fc2 uncertainties are included?
6. How large is the detuned-event convention effect (review 2.2) in the infrared-dominated regime? There the
   permitted detuning (4 sigma = 0.4 THz) is up to 20 % of the mode frequency.

## 10. Recommendation
**Major revision.**
1. **What stands.** The numerical work is careful and reproducible: every number this review re-derived agreed.
   S1–S4 (as 5x5x3 inner bounds), S9 (mapping), S11 and S12 stand.
2. **Why the conclusions do not stand.**
   1. S5, S6 and "all large changes are artefacts" are contradicted by the campaign's own E19 table and by the
      review's transfers.
   2. The "~20 %" robust figure is not a measured quantity. The admissible set is empty, the searches were
      untargeted, single-start and budget-limited, the downward side is unexplored, and the tolerances are
      stricter than first-principles precision.
   3. The transfer test has no demonstrated power, and its criterion rejects the real second workflow.
   4. The physical premise of a converged bulk tau_mem is unexamined. The campaign's own operators are consistent
      with an infrared divergence which, if confirmed, overturns the recommendation for paper 1.
3. **Required headline.** S7 (UNRESOLVED) is the correct overall status and should be the only headline. The
   withdrawn conclusion should stay withdrawn. The present conclusion should be replaced by a statement of three
   things:
   * **Falsified:** the 5x5x3 witnesses as properties of AlN.
   * **Observed:** large mesh-specific ranges on each discretised operator; tau_mem carried by the lowest acoustic
     modes; the second workflow within 7 % on three meshes.
   * **Open:** existence of the limit, robust identifiability, and the realistic-data question.
4. **Path to a stronger paper.** T1 and T7 together could turn the campaign into a substantially stronger result: an
   infrared characterisation of the memory moment with an asymptotic identifiability statement.

## 11. Files of this review
All under `review/peer-review-checks/`. The scripts read the C: mirror results.

| file | content |
|---|---|
| `referee_check.py` | independent minimal operator from the definitions; reference and witness ratios; degeneracy-averaged lifetime statistics |
| `referee_ir.py` | infrared structure: share of N by frequency; low-frequency exponent of r; tau versus 1/f_min |
| `referee_anatomy.py` | spectral location of each witness's dN and dK |
| `referee_drag.py`, `referee_drag2.py` | diagonal versus off-diagonal (drag) decomposition |
| `referee_floor.py` | tau/tau_CS and the rigorous floor per mesh |
| `ref_m553.json`, `ref_m775.json`, `ref_m996.json` | raw outputs of `referee_check.py` |
| `referee_ir.json`, `referee_anatomy.json`, `referee_drag.json`, `referee_drag2.json`, `referee_floor.json` | raw outputs of the other scripts |

**Environment.** Research Python (NumPy, SciPy), 4 threads, peak memory below 2 GB. Run time is about 2 min per
script on 9x9x6.
