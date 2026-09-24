# Frozen review target, 23 September 2026
This is a candidate for hostile review, not an accepted final result.

Implementation: C:/Users/Koussay/these/src/collision_events.py
Tests: C:/Users/Koussay/these/tests/test_collision_events.py
Draft note: C:/Users/Koussay/these/notes/19_Validation_evenements_collision.tex

Candidate claim: for finite, positive-frequency, exactly resonant, distinct-mode
reversible three-phonon events on a common-weight full grid, the stated nonlinear
Bose flux linearizes in entropy coordinates to C=R.T R. A sparse reference
implements that formula with explicit user-supplied rate conventions.
It is not a material exporter. Repeated indices and unequal weights are excluded.
Tolerance admits roundoff mismatches, not physical broadening; residual energy
and thermal detuning are exposed. Extreme occupation factors may underflow.

Strongest additional candidate: matching a homogeneous odd-sector collision
action can miss self-reciprocal event miscounting; streaming exposes the error
at nonzero spatial wavevector (branch C toy examples, not AlN).

Evidence so far: 20 new tests; complete repository suite 225 passed in 192.26 s;
note 19 PDF compiled with zero reported warnings. These are arithmetic/code
checks, not a validation of actual AlN rates or a novelty claim.

Attack the formulas, numerical contract, source attribution, missing invariants,
event counting, or overclaims. Do not defend or improve the claims. Record
severity, exact counterexample/reproduction and what remains unresolved.
