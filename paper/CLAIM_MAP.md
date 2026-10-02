# Manuscript claim-to-evidence map

Working synthesis, 27 September 2026. This maps reviewed local evidence; it does not certify novelty or replace a final manuscript audit. Historical branch reports are subordinate to each campaign's final report and errata.

| Candidate statement | Evidence trail | Accepted scope | Cannot support |
|---|---|---|---|
| Grey kinetic reduction gives an effective closure | [Note 14](../notes/14_Derivation_fermeture_cinetique.pdf), [scope plan](RESEARCH_PLAN.md) | Derivation under its grey/asymptotic assumptions; retain hypotheses in manuscript | Quantitative spectral AlN coefficients without a material reduction |
| Static/diagonal information does not identify arbitrary positive dynamics | [Note 18](../notes/18_Audit_fermeture_spectrale.pdf), [reviewed operator-class report](../theory/aln/closure/13-final-report.md) | Explicit finite matrix classes and counterexamples | Microscopic AlN nonidentifiability from an unrestricted matrix example |
| Fixed finite nonnegative event geometry changes memory bounds | [Closure proof/audit trail](../theory/aln/closure/13-final-report.md) | Fixed current, dissipative space and event cone; finite-class hypotheses | A uniform useful continuum bound or an actual AlN closure |
| A finite exactly resonant Bose-event reference has a positive entropy action | [Note 19](../notes/19_Validation_evenements_collision.pdf), [event final report](../theory/aln/events/13-final-report.md), [implementation](../src/collision_events.py) | Equal weights, distinct modes, declared kinetic flux; compensated detuning and explicit approximation caveats | Broadening as roundoff, physical event selection, repeated-mode extension or absolute material rates |
| Selected material exports can be reproduced and checked | [Note 20](../notes/20_Export_microscopique_AlN.pdf), [pilot archive](../theory/aln/interaction_pilot/README.md) | Pinned inputs and selected coarse outputs; original failure and corrected cutoff retained | Mesh-converged conductivity/collision action, force-constant uncertainty estimate |
| Invariant block sums do not imply closed block dynamics | [Note 21](../notes/21_Degenerescences_et_fermeture.pdf), [reviewed report](../theory/aln/degenerate_action/13-final-report.md) | Finite counterexamples and explicit weighted closure criterion; established mathematical context | Universal scalar-closure failure or a material coherence timescale |
| Tensor, Fock, spectral and population counting factors differ | [Note 22](../notes/22_Orientation_et_comptage.pdf), [reviewed counting report](../theory/aln/channel_counting/13-final-report.md) | Explicit symmetric Hamiltonian, compatible sewing, kinetic/geometric assumptions | Measured factor-two lifetime correction or phono3py defect |
| Geometric populations are not an invariant family of the assumed finite Fock jump process | [C2 derivation](../theory/aln/channel_counting/branches/C2-closure-check.md), [proof review](../theory/aln/channel_counting/branches/red-proof.md) | Exact instantaneous moment obstruction in the declared model; equilibrium exception | Quantitative AlN closure error or disproof of a controlled continuum kinetic limit |
| Existing spectral/FDTR calculations illustrate sensitivity | [Note 16](../notes/16_Analyse_experimentale_FDTR.pdf), [Note 17](../notes/17_AlN_spectral_temperature.pdf), [reviewed RTA limits](../theory/aln/closure/13-final-report.md) | Synthetic design and discrete bulk illustration; keep separate datasets/assumptions | Laboratory validation, converged temperature trajectory or measured hydrodynamic window |
| A finite entropy-variable weak form integrates resonant events with proved local structure | [Note 23](../notes/23_Integration_resonante.pdf), [final report](../theory/aln/resonance_measure/campaign/13-final-report.md), [four reviews](../theory/aln/resonance_measure/campaign/12-peer-review.md) | Finite hypotheses: positive quadratures, full rank, positive entropy variables, exact energy representation, exact resonant nodes; local solution only | Global realizability (finite-time exit shown), rate or measure validation, an AlN operator, novelty |
| Unchanged positive off-shell Bose events cannot conserve energy | Note 23 §2; campaign proof audit | Unchanged event rows and a common Bose rate | A verdict on all broadening or spectral theories |

## Audit admissions, 2 October 2026
The conserving-resonance-measure campaign is closed after four independent hostile reviews; the two rows above are admitted.
The repository audit also established that the formal Fourier-resonance line x = 1.8 depends on τ_c ≈ τ_N (x = 0.8 historical, none conservative, with τ_c kept; notes 07, 12–14),
that isotopes do not explain the Rao/Olympics κ difference (0.1–0.2%; input uncertainty about 15–20%, note 17),
and that the default-cutoff export differs from corrected exports by degenerate-basis redistribution plus a defect-level perturbation of about 4e-8 (interaction_pilot/audit-20261002).

## Manuscript decisions still needed
1. Select one central contribution and check prior art for that precise statement.
2. Trace each final theorem to assumptions and an actual proof; move ancillary audits to supplements.
3. Report numerical errors against the right reference, distinguishing arithmetic, discretization and physical uncertainty.
4. Use one notation for occupations, entropy variables, event weights, operator domains and spectral measures.
5. Include physical unknowns explicitly; do not imply the available notes form a completed experimental article.
