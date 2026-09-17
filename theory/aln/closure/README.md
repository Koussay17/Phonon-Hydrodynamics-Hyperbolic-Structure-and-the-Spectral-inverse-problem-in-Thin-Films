# Spectral closure campaign - evidence snapshot

IN PROGRESS: first-pass reports are frozen; independent red-team reviews and second-generation work are pending. No quantitative real-AlN hydrodynamic closure is certified.

Campaign: 20260916-210742-aln-spectral-closure, persistent original under ResearchLab/orchestration/campaigns/. Copies here make the relevant evidence available with the repository.

## Reproduction
From the repository root:

```powershell
python -B theory/aln/closure/experiments/constructive_collision.py
python -B theory/aln/closure/experiments/infrared_audit.py
python -B theory/aln/closure/experiments/pi_response_bounds.py
python -B scripts/analyse_aln_response.py
```

Install requirements-dev.txt and mpmath. The copied infrared script resolves the repository root relative to itself; this is the only change from its original campaign version. Outputs are regenerated alongside the scripts. Original results can reflect a different Python interpreter; each result's stated environment governs its numeric provenance. Reports are historical first passes, not substitutes for the eventual reviewed synthesis.

## Interpretation
Graph and matrix examples concern explicitly specified finite operator classes. They are not alternative AlN force-constant realizations. AlN sums use one published 300 K RTA mesh; stable arithmetic does not certify continuum convergence. Film suppression is stationary in-plane conduction. Harmonic bulk conductivity is not a FDTR measurement.
