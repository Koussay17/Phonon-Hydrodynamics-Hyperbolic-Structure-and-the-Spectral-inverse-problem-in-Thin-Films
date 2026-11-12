import ResearchProofs.FdtrGauge.RationalIdentity
import ResearchProofs.FdtrGauge.KernelLemma
import ResearchProofs.FdtrGauge.MomentCoefficients
import ResearchProofs.FdtrGauge.MomentMinors
import ResearchProofs.FdtrGauge.LeadingOrder
import ResearchProofs.FdtrGauge.RiccatiHierarchy

/-!
# Campaign `20261111-fdtr-gauge-memory`: machine-checked finite steps of `proofs/B-analytic.md`

* `RationalIdentity`: Theorem 1, Steps 4–6 (rational identity for `g`, fibres, case assembly).
* `KernelLemma`: the `2 × 3` two-load kernel lemma (Theorem 1, Steps 1–3; Theorem 2, Proof A)
  and the algebraic endgame of Theorem 2, Proof A.
* `MomentCoefficients`, `MomentMinors`: the minors `G_K` of the moment Jacobian (Proof C,
  Theorem N2), with the Laurent coefficients taken as given rational functions.
* `LeadingOrder`: Theorem N2, step (3) (determinants `t_K`, index sets, leading coefficient).
* `RiccatiHierarchy`: an independent characterisation of the Laurent coefficients.
-/
