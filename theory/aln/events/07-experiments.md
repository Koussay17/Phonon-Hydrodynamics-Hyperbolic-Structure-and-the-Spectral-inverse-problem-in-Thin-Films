# Reproduce
From the repository root:
    python -B -m pytest tests/test_collision_events.py -q -p no:cacheprovider
    python -B theory/aln/events/experiments/C_parity_events.py
    python -B theory/aln/events/experiments/D_bose_jacobian.py
    python -B theory/aln/events/experiments/second_import_check.py
    python -B theory/aln/events/experiments/red_counter_checks.py
    python -B theory/aln/events/experiments/red_numerical_api.py
    python -B theory/aln/events/experiments/red_proof_event_checks.py
Historical audits deliberately import the BROKEN frozen review-inputs version
so the original failures remain reproducible. They do not test current src/.
Current regressions test src/collision_events.py. Post-fix replay procedure
and outputs are separately labelled in post-fix-replay.md.
Dependencies: requirements-dev.txt. No AlN input payload is needed for these tests.
