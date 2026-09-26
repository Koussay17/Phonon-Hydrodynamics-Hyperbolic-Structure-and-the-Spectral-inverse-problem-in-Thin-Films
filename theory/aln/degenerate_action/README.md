# Degenerate collision action: audited research checkpoint

Start with [the final report](13-final-report.md) and [claim ledger](CLAIM_LEDGER.md).
This package separates finite mathematical results, synthetic examples and a selected material contraction. It does not provide a complete AlN collision operator.

## Reproduce the finite checks

Use Python 3.14 and an isolated environment containing requirements-analysis.txt.
From the repository root:

    python -B scripts/reproduce_degenerate_action.py --output-dir D:/path/to/new-results

This creates a new directory and runs the exact lumpability example, phase-information example, memory comparison and clustering regressions. The memory test compares the pinned repository event source; both canonical LF and tested Windows CRLF hashes are accepted. Other source versions are not certified.

## Include the selected material contraction

First reproduce the previous interaction_pilot and its cutoff continuation; keep its original directory structure and hash-matching raw outputs on local storage.
Use the prior pinned phono3py environment:

    python -B scripts/reproduce_degenerate_action.py --output-dir D:/path/to/new-results --pilot-root D:/path/to/prior-pilot --material-python D:/path/to/phono3py-env/python.exe

Material output remains in the selected scratch directory. The script rejects existing output directories, verifies prior export hashes, and rejects ambiguous frequency chains. It is specific to the recorded nonzero-frequency triplet 4 and pinned internal APIs, not a general importer.

## Scope and historical artifacts

The original complex_triplet.py and reconstruct_complex_triplet_legacy.py preserve the unsafe clustering rule for the archived failure. Use reconstruct_complex_triplet.py with frequency_blocks.py for new runs.
R-numerical-audit.py inspects the original saved tensor artifact; it requires original local NPZ/input paths and is not part of the portable finite-check runner. Raw force constants and amplitude NPZ files are not redistributed here.

All four first passes were independent. Second-generation reports and four hostile reviews are preserved. The clustering fix was reviewed separately. Exact algebra is not formal verification, and source-consistent contraction is not independent physical normalization. No novelty is claimed.
