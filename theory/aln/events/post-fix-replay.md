# Post-fix integration replay

The unchanged independent C implementation was compared with the corrected live operator. In a D scratch copy only the target file and its expected hash were replaced; the historical second-generation artifacts were not edited.

Corrected operator SHA-256: 77196934d0e3ee209e7121f31655af48e074864d7fff5668cfb6f07b1ea99e6c

Command: python -B D:/ResearchLab/scratch/20260924-event-postfix-replay/experiments/second_import_check.py
Results: post-fix-import-results.json. This is a PI regression replay, not a new independent derivation.
