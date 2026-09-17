# Audited spectral closure campaign — final report

**Date:** 17 September 2026. **Status:** bounded research cycle completed; quantitative AlN closure and experimental work remain open.

## Exact question

What can the available mode frequencies, velocities, diagonal rates and bulk conductivity determine about a conserving spectral hydrodynamic closure for AlN, and what additional data and mathematical conditions are indispensable?

The work concerns linearized entropy-coordinate kinetics. It distinguishes finite matrix statements, microscopic event constraints, a continuum limit, and the actual film experiment. The existing 300 K Rao spectrum and independent temperature table are not combined into fabricated normal/resistive lifetimes.

## Strongest justified conclusion

**The admissible operator class matters.** In the broad symmetric-positive class, fixed diagonal rates and even a fixed DC tensor allow different or unbounded memory. But in a **fixed finite nonnegative event cone**, fixed current and dissipative space, bounded DC response implies bounded first memory moment. The latter follows from an audited Cauchy–Binet/convex-combination proof; no lower spectral gap is needed for that finite-class bound.

These statements are compatible, not a disagreement to resolve by voting. Neither supplies a useful numerical continuum bound or a unique hydrodynamic model for real AlN. The microscopic operator and its convergence remain missing.

## Classification of results

| Status | What has actually been established |
|---|---|
| **PROVED** | Explicit four- and eight-mode algebra, B1 under its feasibility hypotheses, surrogate pseudoinverse, conditional finite Schur/Jensen identities, and fixed finite event-cone bound. Independent proof audits give complete arguments. |
| **FORMALLY VERIFIED** | None. Symbolic/rational checks are not a proof-assistant formalization. |
| **DERIVED UNDER ASSUMPTIONS** | Slow/fast kinetic reduction, infrared moment thresholds, reflecting-film asymptotics and conditional physical interpretation. Stated regularity, parity, source, gap and boundary assumptions are essential. |
| **NUMERICALLY DEMONSTRATED** | Reproducible discrete 300 K RTA diagnostics; exact-rational/high-precision synthetic checks; both singular-endpoint floating-point failures; Schur test cases and numerical cancellation correction. |
| **LITERATURE SUPPORTED** | Collective collision modes, viscous heat equations, projection memory, event decomposition, and the odd-equivalence limitation of the inspected Chaput/phono3py operator path. Exact read scope is recorded. |
| **SUPPORTED BUT UNPROVED** | A useful small, observable-specific closure for the eventual AlN operator is a plausible route; no accuracy/usefulness certificate for that material has been obtained. |
| **CONJECTURED** | Mesh-converged dressed variables or an informative passive infrared tail for AlN; the four discovery proposals remain conditional research directions. |
| **DISPROVED** | Positive diagonals certify a collective gap; DC agreement certifies memory accuracy; value asymptotics alone imply a pointwise rate-density asymptotic; bounded DC permits unbounded memory in one fixed finite event cone. |
| **UNKNOWN** | Actual AlN continuum memory behavior, N/U trajectory, physical collision/stress coefficients, converged hydrodynamic window and experimental validity. Novelty is not established or claimed. |

## Quantitative evidence and its limits

- Bulk RTA at 300 K: conductivity 309.775 basal and 302.443 along c, in W/(m K); finite-mesh memory 149.30 and 231.66 ps.
- Modes below 2 THz carry 61.61% and 74.42% of the second inverse-rate moment, but only 13.48% and 16.16% of DC conductivity. This is a sensitivity diagnostic, not an extrapolation.
- At 200 MHz, the single pole matched to memory has 6.49% and 14.73% relative complex error against the total discrete RTA response. Relative to the dynamic correction, the same errors are 42.09% and 79.64%. These are not FDTR measurement errors.
- Independent 75-digit response sums agree within 2.363e-16 relative for the stored inputs. This verifies arithmetic, not input physics.
- Exactly grazing nodes create a formal fixed-grid in-plane floor of 40.165 W/(m K). It is not a continuum or atomically thin film prediction, nor a quantified error at any finite thickness.
- Near the eight-mode family's t=1 endpoint, computed spectral memory can err by 77.6% while DC is accurate to machine precision. Near t=0, even very small residuals coexist with large inverse-moment errors. Exact proofs remain valid.
- Second-generation tests compare 600 full and eliminated systems, with high-precision coordinate/pole cases. The report states which pole examples fall outside the causal-half-plane bound.
- Repository validation after the numerical fix: **205 tests passed in 186.54 seconds**. Five new tests compare against independent 100-digit direct subtraction; 23 targeted tests passed.

## What the red team changed

1. Rejected material-specific use of the abstract counterexample. A new fixed-event-cone inequality and stronger memory bound received a separate proof audit.
2. Found an explicit monotone rate with oscillatory derivative that breaks a pointwise density inference. The synthesis adds the missing derivative assumption.
3. Found an additional numerical failure near the other singular endpoint, t=1.
4. Identified cancellation in the reported very small error, duplicate implementations mistaken for independent evidence, and missing provenance. The bulk diagnostic now uses algebraic differences and records versions/input/script hashes.
5. Confirmed that scalar array dimensions do not establish the physical even collision sector. Independent primary-code/equation checks support the restriction on the inspected Python/C path.

The physical/source report retains optional unexecuted normalization/event-test caveats. No exported AlN operator was validated or used.

## Four alternative formulations and why they were not implemented yet

- **Dressed slow variables:** a finite static Schur form can have an infinite induced entropy metric. Require existence, convergence and a useful full streaming residual.
- **Passive infrared continuum:** a zero directional leading coefficient defeats a universal divergent-slope claim; temporal data do not identify spatial streaming. Require nondegenerate tail and separate spatial validation.
- **Positive spectral-measure envelopes:** useful response bounds do not identify memory; exact DC, mass and a finite-shift sample can coexist with unbounded slope in that relaxed class. Require certified envelope widths and feasible constraints.
- **Joint collision/streaming space:** distinct velocities and nonzero energy support can force this space to be the full mode space. Exact reduction may offer no compression; approximation needs a band-specific certificate.

These attacks were completed before substantial implementation. No branch was continued merely to justify prior effort.

## Prior art and access limits

Projection memory, relaxons, variational matrix response and viscous heat equations are established frameworks. The exact read scope and URLs are in the two literature reports and peer-review source notes. The 2020 relaxon erratum corrects a Matthiessen inequality; it does not correct the zero-mode statement. The 2026 preprints remain provisional and do not validate wurtzite AlN here.

The local multi-index search included Crossref, OpenAlex and Semantic Scholar, plus primary web sources and bounded citation tracing. Rate limits, unavailable full texts and uninspected software paths are explicit. There is no exhaustive novelty claim.

## Reproduction

From the repository root, install requirements-dev.txt (including SymPy and its mpmath dependency), then:

~~~powershell
python -B scripts/analyse_aln_response.py
python -B theory/aln/closure/experiments/constructive_collision.py
python -B theory/aln/closure/experiments/infrared_audit.py
python -B theory/aln/closure/experiments/pi_response_bounds.py
python -B theory/aln/closure/experiments/red_proof_checks.py
python -B theory/aln/closure/experiments/red_numerical_audit.py
python -B theory/aln/closure/experiments/second_generation_schur.py
python -B -m pytest tests/ -q -p no:cacheprovider
~~~

The review reports are historical evidence; [errata.md](errata.md) and this synthesis govern corrected conclusions. A prematurely named counterexample script was never executed, and its originating report now explicitly corrects that provenance mistake. It is not included or counted as evidence.

[Note 18](../../../notes/18_Audit_fermeture_spectrale.pdf) retains the preceding notes' typography. Large future datasets/cache and campaign backups use D:/ResearchLab. Extra SSD space does not remove the approximately 28 GiB RAM constraint; one naive dense 24^3-by-12-band matrix alone requires about 205 GiB.

## Why this cycle stops here

The bounded theoretical/audit objective is achieved. The material-specific continuation reaches an identified **data and validated-operator barrier**, not a claim that further physics is impossible. More abstract examples cannot supply the missing AlN collision action, grid series or sample measurements.

The next concrete milestone is a verified physical signed-event action, starting with a small event/reciprocal-pair test and both parity sectors; then one consistent force-constant source, N/U and temperature outputs, and convergence of current/stress responses. No expensive first-principles rerun was performed.

**What is the strongest remaining reason our conclusion could be wrong?** The actual microscopic event geometry, continuum limit and experimental boundaries may impose properties absent from the finite models. Extending any result beyond its stated class, or using a conductivity-equivalent operator for even-sector physics, would invalidate the material inference. Those extensions are expressly not made.

## Final publication validation — 2026-09-17

- Full repository suite: 205 tests passed (186.54 s).
- Portable replays completed: constructive_collision.py, infrared_audit.py, pi_response_bounds.py, red_proof_checks.py, red_numerical_audit.py and second_generation_schur.py.
- The numerical audit uses an explicitly preserved historical pre-correction fixture; the current response generator uses the corrected stable identity.
- Notes 17 and 18 rebuilt successfully: 5 and 6 pages, respectively, with no reported LaTeX warnings. Note 18 retains note 17's preamble and typography.
- PDF contact sheets were rendered, but the image-viewing tool failed in this Windows environment; no completed visual inspection is claimed.
- No material collision operator, continuum grid convergence or experimental validation is inferred from these checks.
