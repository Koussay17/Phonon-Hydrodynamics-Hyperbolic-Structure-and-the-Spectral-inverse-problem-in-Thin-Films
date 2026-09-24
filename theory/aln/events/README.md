# Finite physical-event validation reference

Start with [the synthesis](13-final-report.md), [review resolutions](12-peer-review.md),
and [note 19](../../../notes/19_Validation_evenements_collision.pdf).

The bounded reference-validation step is complete. No actual AlN event exporter,
absolute event-rate normalization, conserving material integration or experimental
validation is claimed. The force-constant download manifest records inputs only.

## Two deliberately different code versions
- Current implementation: src/collision_events.py at the repository root.
- review-inputs/: the original BROKEN candidate, preserved byte-for-byte so the
  adversarial counterexamples remain reproducible. Do not use it as production code.

See [commands and interpretation](07-experiments.md). Historical review scripts
import the frozen version, while tests/test_collision_events.py tests the current
corrected version. The post-fix replay uses the corrected code in an isolated
copy and is labelled separately; it is not an additional independent derivation.

Four independent first passes, literature, second-generation comparison and four
isolated hostile reviews are included. Failed numerical regimes were preserved.
The final synthesis and review resolutions supersede preliminary claims in the
historical reports. Check validation.json for actual final execution results.

The archived inputs/experiments retain exact bytes for hash checks. Manifests
describe working-tree bytes; regenerating outputs may change environment/path
metadata. External source payloads stay on D and are not redistributed here.

The scientific manuscript is still in planning: [scope and completion criteria](../../../paper/RESEARCH_PLAN.md).
