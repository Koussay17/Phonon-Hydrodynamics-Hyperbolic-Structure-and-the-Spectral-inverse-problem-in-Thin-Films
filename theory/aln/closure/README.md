# Audited spectral closure campaign

Start with [the final report](13-final-report.md), [reviewed errata](errata.md), and [note 18](../../../notes/18_Audit_fermeture_spectrale.pdf).

The bounded research/audit cycle is complete. Quantitative real-AlN closure, a normal/Umklapp temperature trajectory and experimental validation remain open. The final report explicitly separates nine evidence categories.

## Contents

- Numbered ledgers: question, assumptions, evidence, disagreements, failures, proof status and remaining work.
- reports/: four independent first passes, literature, four independent red-team reviews, second-generation/discovery reports and their attacks, event-cone proof audit and operator parity check.
- experiments/: reproducible finite-matrix, infrared, symbolic and independent numerical checks, including failed numerical regimes.
- Input/environment/storage manifests and checkpoint records.

Initial reports are preserved as history. The final report and errata supersede provisional wording. No result is credited to an unexecuted proposed script.

## Reproduce

From the repository root, install requirements-dev.txt, then run the commands listed in the final report. SymPy's mpmath dependency supplies arbitrary precision. The bulk diagnostic is scripts/analyse_aln_response.py. Tests include independent 100-digit subtraction checks.

Portable experiment copies resolve the repository/data paths relative to themselves. The infrared generator additionally records provenance; its scientific sums are unchanged. Numerical output metadata and input hashes identify the actual run. artifact-manifest.json lists the snapshot hashes.

## Scope

Graph/matrix examples are not alternative AlN force-constant realizations. The fixed finite event-cone bound is not uniform under mesh refinement. AlN sums use one published 300 K RTA mesh. In-plane reflecting-film suppression and bulk harmonic conductivity are not a cross-plane FDTR signal.

New bulky storage and the campaign mirror use D:/ResearchLab. Existing campaign reports on C: are retained for continuity. No personal files were deleted, no dense full collision allocation was attempted, and no first-principles rerun is claimed.

The experiment fixture pi_bulk_response.json intentionally preserves the pre-correction output so the numerical reviewer can reproduce the originally detected cancellation. The current corrected output is ../bulk_response.json relative to this directory. This fixture is not a second independent implementation.
