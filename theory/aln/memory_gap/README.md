# Memory gap — identifiability of the closure memory time (paper 1), October 2026

Question: do static transport data (mode lifetimes, collision invariants, DC conductivity tensor) determine the
memory time τ_mem = −K′(0)/K(0) of a conserving phonon closure? Start with `13-final-report.md`; reviews in
`review/` (proof audit, counterexample review).

Declared class: arbitrary non-negative rates on the fixed energy-allowed three-phonon event set (crystal and
time-reversal symmetric, selection-rule zeros kept), no microscopic amplitude model.
- Lower end identified: τ_mem ≥ K/|b|² (proved), nearly attained.
- Upper end unidentified within the class: certified witnesses exceed the reference by factors of 10²–10⁵
  (AlN 5×5×3, 300 K: basal τ ≥ 7,038.17 ps, c axis τ ≥ 4,322,183.47 ps, reference 34.25 / 27.45 ps).
- Under a factor-F amplitude prior the gap closes: τ ≤ F/λ_min(C_ref) = 81.50·F ps (proved).
- Lemma: K ≥ K_RTA/s_max for any event operator with r = diag C (proved, sharp).
Conclusion: the information that pins τ_mem is the microscopic amplitude structure, not static transport data.
Large arrays (event sets, operators) are kept outside the repository.
