# Reviewed errata and provenance

The initial A–D and L reports are preserved as research history. This file and 13-final-report.md supersede their unqualified or provisional readings.

1. **C density asymptotic:** value asymptotics plus smoothness on q>0 do not imply the stated pointwise pushed-forward density. Add uniform derivative control and monotonicity, or supply a separate cumulative-measure proof. The explicit counterexample and proof are in red-proof.md. Moment thresholds and the exact q^2 example survive.
2. **Broad class versus event cone:** the eight-mode result is not a microscopic AlN realization. R2 excludes its entire small-t path from a fixed finite event cone; separately audited R3 bounds memory whenever DC is bounded in that fixed cone. The latter does not give a uniform continuum bound or uniqueness.
3. **R2 zero case:** if R=0, the positive cross term is impossible; do not write 1/R. This is handled in event-cone-proof-audit.md.
4. **Numerical intervals:** the original D sweep's good errors apply only to sampled points up to t=0.75, not all t>=1e-3. Red numerical review finds a severe additional failure near t=1.
5. **Bulk errors:** initial tiny errors were subtraction-limited. The current script uses variance/algebraic differences and distinguishes error relative to total response from error relative to dynamic correction. Duplicate scripts were never independent physical validation.
6. **Provenance:** current bulk and infrared outputs include versions and hashes. Frozen reports may name original local paths/interpreters; the reproducible repository scripts and actual JSON metadata control.
7. **Unexecuted artifact:** red-counterexample.md initially named a planned red_counterexample_checks.py/json as if executed. Its appended correction explicitly retracts that statement. No such artifact is included or counted. Event-cone mathematical proof rests on the independently audited derivation.
8. **Jensen domain:** nonzero current is required. Bounds order positive-real-Laplace responses, not complex harmonic components.
9. **Relaxon erratum:** it corrects a Matthiessen inequality, not a zero mode. L2 resolves the earlier uncertainty.
10. **Software scope:** the inspected phono3py Python/C scalar reducible path implements an odd-equivalent conductivity operator. This is not a statement about all backends/versions. Full physical even action and a validated exported AlN operator remain unavailable here.
11. **Discovery claims:** tail divergence requires nonzero directional weight; dressed-state existence, spatial identification and practical compression are unproved for AlN. See discovery-attacks.md.
12. **Proof status:** no theorem is described as machine-formally verified. SymPy checks are supplementary.
