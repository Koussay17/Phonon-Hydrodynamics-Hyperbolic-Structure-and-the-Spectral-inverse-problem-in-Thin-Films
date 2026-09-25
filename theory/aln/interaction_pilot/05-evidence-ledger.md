# Evidence ledger
- pilot-run.json: pinned inputs and exported array provenance; 3x3x3 mesh, 300 K, selected target q only.
- pilot-inspection.json: detailed gamma sum agrees to 1.48e-16 globally scaled; duplicate BZ frequency check FAILS (6.8181e-8 THz versus 1e-10).
- gamma-reconstruction.json: independent accumulation max absolute error 6.94e-18 THz; error divided by maximum reference linewidth 2.96e-16. This is not maximum branchwise relative error.
- results/B_grid_orbits.json: independent integer group action recovers all eight multiplicities [2,1,4,2,4,2,8,4], sum 27.
- pair-checks.json: reciprocal linewidth differences below 2e-13 globally scaled; large width sensitivity.
- harmonic-controls.json and matrix-control.json: defect present with NAC in both backends; eigenpair residuals near machine precision.
- cutoff-control.json: fixed splitting parameter, reciprocal-cutoff factors 1,1.25,1.5,2; local matrix defects 1.35e-8,1.55e-12,4.94e-16,4.94e-16.
- reproduction-status.json: fresh-directory reproduction retains the declared failed original check. Successful process execution does not mean all validation passed.
Raw HDF5/NPZ remain on D; compact results and executable scripts are published.
