# 00 — Question: identifiability gap of the memory moment in a fixed event geometry

**Date opened:** 4 October 2026. **Context:** Paper 1 (PRIOR_ART.md, "Results still needed", items 1–2).
Read-only inputs: thesis repository `C:\Users\Koussay\these` (notes 18, 19, 22, 23;
`theory/aln/transport_baseline`, `theory/aln/interaction_pilot`). Nothing in that repository is edited.

## Exact question
Linear, energy-conserving three-phonon kinetics in entropy coordinates,
`C(g) = sum_alpha g_alpha z_alpha z_alpha^T`, `g >= 0`, with a FIXED event set `{z_alpha}`
determined by the harmonic spectrum (momentum and energy selection), `C e = 0`,
current `b` orthogonal to every invariant.

Static data held fixed:
* mode relaxation rates (lifetimes): `diag C(g) = r`;
* the collision invariants: `ker C(g) = span(e)` (energy only, umklapp present);
* the DC response `K(0) = b^T C(g)^+ b = K0`.

Unknown: the nonnegative event rates `g`.

Target: the feasible range `[tau_min, tau_max]` of the first memory moment
`tau_mem = b^T (C^+)^2 b / b^T C^+ b = -K'(0)/K(0)` (the Cattaneo time of the matched single pole),
and in particular whether `tau_max / tau_min` is O(1) ("tight") or orders of magnitude ("large").

## Sub-questions
1. Formulation of the feasible set; tractable outer bounds (Ben-Tal–Teboulle constant `M`,
   `tau <= M^2 K0`; convex/SDP relaxations; Cauchy–Schwarz/Jensen lower bounds) with proofs.
2. Debye-type event sets (isotropic linear dispersion, periodic meshes in d = 1, 2, 3,
   normal + umklapp three-phonon events, exact resonance where possible, otherwise a declared
   tolerance with conserving stoichiometry): exact `M` by enumeration for small sets, certified
   bounds, and the feasible `tau_mem` interval by optimization (inner interval = explicit
   feasible witnesses; outer interval = proved bounds).
3. Scaling of the interval with mesh size N and infrared rate exponent alpha (`r ~ q^alpha`);
   the continuum statement that the memory moment diverges when `2 alpha >= d`.
4. AlN event geometry (phono3py 4.5.0, Phonon Olympics force constants, physical operator
   convention Omega' rather than phono3py's odd-sector-equivalent Omega), 300 K:
   largest feasible mesh; `tau_mem` interval with lifetimes and K(0) fixed; Cattaneo time in ps.

## Decision this campaign informs
Whether Paper 1 can claim a *quantitative* identifiability gap (orders of magnitude) in realistic
event geometries, or only a qualitative non-identifiability with a modest spread.
