# Existing integration interface audit — first-pass PI inventory

This is code inspection, not a physical convergence conclusion. No new peer report read.

## Repository reference
src/collision_events.py accepts equal-weight, distinct-index, exactly resonant finite events. Compensated detuning is checked against both absolute dimensionless and relative roundoff bounds. Its documented PSD approximation at accepted nonzero mismatch is not the exact nonlinear Bose Jacobian. Repeated modes, irreducible weights and spectral broadening lie outside this API.
Do not expand its resonance tolerance to accommodate integration widths.

## Pinned phono3py 4.5.0
phono3py/phonon3/triplets.py get_triplets_integration_weights returns g with two or three components.
Gaussian Python path evaluates:
g0 = gaussian(f0-f1-f2, sigma)
g1 = gaussian(f0+f1-f2, sigma)
g2 = gaussian(f0-f1+f2, sigma)
Stored g[0]=g0, g[1]=g1-g2; collision-matrix path also stores g[2]=g0+g1+g2.
Therefore array components are signed channel combinations, not uniformly nonnegative unique-event measures.

The Python tetrahedron path uses daughter frequency values at paired vertices with opposite q shifts; it integrates f1+f2, -f1+f2 and f1-f2 via TetrahedronMethod before forming those same combinations.
This computes fixed-external-frequency integration weights. It does not itself expose a globally unique resonant-triple quadrature with energy/occupations evaluated at its points. A direct event importer cannot assume one positive saved mesh weight equals one exactly resonant mesh event.

No source defect is asserted. A redesigned weak-form action may have a different state representation and must document its relation to this existing scalar integration.

## Evidence
Inspected source files and repository reference are hashed in experiments/source-interface-hashes.json. Existing pilot provides a preserved comparison target, not a convergence certificate.
