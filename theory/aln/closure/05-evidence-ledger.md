# Evidence ledger

| ID | Evidence | Scope |
|---|---|---|
| E0 | Starting clean commit 74a4826 and input hashes in environment.json | Repository baseline |
| E1 | 18 existing spectral tests passed on campaign start | Code regression |
| E2 | Exact graph/Walsh derivations plus rational, symbolic and precision-refined solves | Finite operator class |
| E3 | red-proof.md and event-cone-proof-audit.md | Independent mathematical audit; no formal proof assistant |
| E4 | red_numerical_results.json | Exact-rational and 50/80/110-digit checks; both singular endpoints |
| E5 | Independent 75-digit AlN sums, maximum response discrepancy 2.363e-16 | Stored binary64 input arithmetic, not physical accuracy |
| E6 | infrared_benchmarks.json | Synthetic radial/angular benchmarks, not AlN mesh convergence |
| E7 | second_generation_results.json: 600 Schur cases, high-precision coordinate/pole examples | Synthetic checks; report distinguishes causal domain |
| E8 | Five new independent high-precision tests; 23 targeted tests passed | Cancellation-free diagnostic |
| E9 | Full suite: 205 passed in 186.54 s | Repository test run after numerical changes |
| E10 | L2 and independent operator-parity-check source inspection | Exact version/code path; no AlN operator run |

Source-rate, mesh, material and boundary uncertainties are not replaced by E5. Unexecuted proposed checks are not evidence; see errata.md.

## Final publication validation — 2026-09-17

- Full repository suite: 205 tests passed (186.54 s).
- Portable replays completed: constructive_collision.py, infrared_audit.py, pi_response_bounds.py, red_proof_checks.py, red_numerical_audit.py and second_generation_schur.py.
- The numerical audit uses an explicitly preserved historical pre-correction fixture; the current response generator uses the corrected stable identity.
- Notes 17 and 18 rebuilt successfully: 5 and 6 pages, respectively, with no reported LaTeX warnings. Note 18 retains note 17's preamble and typography.
- PDF contact sheets were rendered, but the image-viewing tool failed in this Windows environment; no completed visual inspection is claimed.
- No material collision operator, continuum grid convergence or experimental validation is inferred from these checks.
