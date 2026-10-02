# Independent numerical hostile review — checkpoint

Date: 28 September 2026. Reviewer scope: A/B/C/C2 and the PI reproduction runner in campaign 20260927-164023-conserving-resonance-measure. D is excluded from self-approval. No other red reports read. No material computations or repository edits authorized for this review.

## Status

IN PROGRESS. The frozen REVIEW_CANDIDATE.md and reproduce_resonance_measure.py have been read. Earlier access attempts were interrupted; the resumed reads succeeded. No independent numerical claim has yet been accepted.

## Checks to complete

1. Derivative checks and cancellation near zero affinity, including the recorded C2 assertion failure and tolerance change.
2. Physical versus surrogate energy, root coverage, regular-root conditioning, and quadrature error.
3. Independent small arithmetic/reproduction checks and audit of thresholds/provenance. The runner copies Python inputs into a new output directory, so pre-existing JSON is not directly reused, but its checks still depend on candidate assertions and their scope.

Final findings will replace this checkpoint. Unexecuted work remains explicitly unverified.
