# Real AlN interaction export pilot

Read [the final report](13-final-report.md). **The original export fails one declared periodicity check.** A fixed-Lambda cutoff experiment controls the discrepancy for one matrix pair; it does not correct the original exported interactions.

## Reproduction

Use an isolated Python 3.14 environment (recorded: 3.14.7). Install requirements-pilot.txt there, not into the main project environment:

    python -m pip install -r theory/aln/interaction_pilot/requirements-pilot.txt
    python -B scripts/run_aln_interaction_pilot.py --input-dir D:/path/to/AlN_inputs --output-dir D:/path/to/new-pilot

The output directory must not exist. POSCAR, BORN, fc2.hdf5 and fc3.hdf5 must match pilot-run.json byte counts and SHA256 hashes. Input source: [pinned Phonon Olympics directory](https://github.com/McGaughey-Lab/Phonon-Olympics/tree/0640f07735059be9717a7565c2a0f22dc0da7a17/Aluminum%20Nitride/phono3py/AlN_kappa_input_files). Software source provenance is in branches/C-unit-source-trace.md.

The runner records each command and preserves the original failed BZ_duplicate_frequencies validation. A zero runner exit means the experiment executed as recorded, not that every physical validation passed. Large HDF5/NPZ outputs stay in the selected scratch directory.

## Interpretation

The mesh is 3x3x3, at 300 K with selected external points. Full interaction arrays at a selected point are not a full crystal collision operator. Gaussian widths are 0.05, 0.1 and 0.2 THz. The direct accumulator shares upstream interactions, phonons and constants; its independent check covers accumulation only.

First-pass reports and reviews remain historical records. The final report and later cutoff review distinguish subsequent evidence from what reviewers initially knew. No novelty is claimed.

## Cutoff continuation (25 September 2026)

After the base runner completes, execute export_cutoff_comparison.py from its output directory using the same isolated Python. It creates three new export subdirectories (existing directories are rejected), repeats the original inspection and direct reconstruction for cutoff factors 1.5 and 2, and compares with factor 1. Lambda stays fixed.

    python -B D:/path/to/new-pilot/export_cutoff_comparison.py

Read the continuation report before interpreting raw pp differences. Original records are retained; the new exports have separate hashes and paths.

The independent saved-array diagnostic review_export_cutoff.py reads run paths from the locally generated export-cutoff-comparison.json and prints JSON. Run it in the same output workspace after the continuation; the archived comparison paths refer to the original D-drive runs and are provenance, not portable input locations.

    python -B D:/path/to/new-pilot/review_export_cutoff.py

The publication preserves each generated pilot_run.py beside its metadata because those scripts contain the effective NAC settings. The inherited metadata settings alone do not list the enlarged cutoff.
