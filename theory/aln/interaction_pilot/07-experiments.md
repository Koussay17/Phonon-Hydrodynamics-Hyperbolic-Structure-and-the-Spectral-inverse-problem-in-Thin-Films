# Reproducible experiments
Run scripts/run_aln_interaction_pilot.py from the repository in the pinned isolated environment, using a fresh output directory and hash-matching input directory. See README.md here.
The runner executes the original pilot, reciprocal/width controls, direct reconstruction, harmonic toggle controls, gauge-aligned matrices, independent array review and integer orbit test. The additional cutoff control holds Lambda fixed while increasing only G_cutoff.
Keep the original failing 1e-10 THz duplicate-frequency threshold. Do not relabel original data as corrected.
Most informative next experiment: rerun phonons AND pp/gamma with an explicitly enlarged reciprocal cutoff, compare all duplicate representatives and target rates, then examine joint mesh/width refinement. Do not infer convergence from this one-pair matrix diagnostic.

Update, 2 Oct 2026: the enlarged-cutoff rerun was performed on 25 Sep (continuation-export-report.md). Remaining next experiments: a degenerate-block comparison of the default export against factor 1.5 (see 06-disagreements.md), then joint mesh/width refinement.

Update, 2 Oct 2026: the default-vs-corrected block comparison is done (audit-20261002/README.md). Remaining: joint mesh/width refinement.
