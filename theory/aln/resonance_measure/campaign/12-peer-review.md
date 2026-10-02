# Hostile review — completed 2 October 2026

All four independent reviews of the frozen candidate are complete (frozen input: REVIEW_CANDIDATE.md, note23-draft.tex, review-input-hashes.json). The reviewers did not read one another's reports.

| Review | Report | Verdict |
|---|---|---|
| Proof audit | branches/red-proof.md (28 Sep) | The finite local theorem is VALID under its hypotheses; 21 exact checks pass. Required: a trace hypothesis for the coarea identity, and a density bound for the A2 tube rate. |
| Counterexample | branches/red-counterexamples.md (28 Sep) | The narrow finite/local theorem is not refuted. Found: frozen C2 returns NaN for valid states; the domain is left in finite time (T ≤ 1.11e-7); the strong forms of the A2/B2 alternatives fail; the integrability wording is imprecise. |
| Peer | branches/red-peer.md (28 Sep) | No fatal issue for the restricted theorem, but promotion to a convergent AlN method is rejected. Major issues: invariants do not constrain the rates; excluded sectors; representation; coincident labels versus the continuum surface. |
| Numerical (final) | branches/red-numerics-final.md (2 Oct; supersedes the unfinished red-numerics.md checkpoint) | FIX VERIFIED: log-space F and Λ are backward stable and finite on [1e-12,1e3]^3 (30,021 inputs, Arb and mpmath references); the harness passes; no reported number is wrong. One confirmed test-design defect (the vacuous revised C2 thermal-identity assertion) and eight fragilities. |

## Resolutions by the principal investigator (2 Oct 2026)
- **NaN in reaction_data.** Fixed by log-space evaluation with a branch on Re(d), which keeps the complex step analytic. The Decimal regression now uses the exact float inputs, requires zero flux at resonance, and sets the flux tolerance at the conditioning bound 8uκ. The reviewer's independent case (Φ=I2, α=(1,2)) is reproduced to all digits.
- **python -O.** Every script, including D-limits-check, and the runner refuse it before creating any output. The runner also removes PYTHONOPTIMIZE from the child environment.
- **Vacuous C2 assertion.** Replaced by two checks: the a-priori distributivity bound |d−βΔ| ≤ γ_n Σ|b||α|, and |E−PS−T| ≤ 16u|E| with T summed by fsum. A permanent negative control (α=1.3βe, 30% violation) must fail, and does.
- **errstate masking.** The selected outputs must be finite; otherwise FloatingPointError is raised.
- **Complex-step cancellation in ∂Λ for 1e-5 ≤ |d| ≤ 1e-2.** The exprel series now extends to |z| < 1e-2 with 11 terms (truncation < 3e-28).
- **Note 23 revised:**
  - coarea trace hypotheses;
  - drift order σ² for a bounded density and σ^{3/2} for a fold, with constant 0.4300 (independently verified as 0.43001999366);
  - explicit finite-time domain-exit example, independently integrated (δ reaches 0 at t ≈ 9.46e-8);
  - entropy-decreasing interpolation example, verified (−3.5877701352e-3);
  - invariants do not constrain the rates;
  - coincident labels versus the continuum tube;
  - log-space evaluation.
- **Accepted as documented limitations, not fixed:**
  - complex-step underflow for |F| ≲ 1e-280;
  - the harness direct-product Jacobian for ξ ≪ 1 (harness states have ξ ≥ 0.69);
  - genuine floating-point underflow of Λ for ξ ≳ 372 at d=0 (the true value is below the double range);
  - quadrature-identity checks, which are consistency checks, not independent evidence.

Post-fix reproduction: 14/14 checks (repository runner, 2 Oct 2026).
