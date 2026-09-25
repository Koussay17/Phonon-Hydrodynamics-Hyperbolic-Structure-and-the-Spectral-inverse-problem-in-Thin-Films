# Retained failures
The first inspection script failed only during JSON serialization of a NumPy boolean. Checks were converted to native bool, and the script rerun; no numerical formula was changed.

The first NAC toggle was invalid: explicit born_filename overrides is_nac=False in the pinned loader. Those outputs are preserved as harmonic-controls-invalid-toggle.json and are not evidence of NAC independence. The corrected control omits born_filename when disabled and records actual loaded NAC state.

The original duplicate-frequency validation fails at 6.8181e-8 THz against 1e-10 THz. This remains a failed export check even though a later fixed-Lambda matrix control reduces the tested defect to roundoff. No tolerance was relaxed.
