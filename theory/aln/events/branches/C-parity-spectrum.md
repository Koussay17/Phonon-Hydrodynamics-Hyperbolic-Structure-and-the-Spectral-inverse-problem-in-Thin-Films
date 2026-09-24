# C: independent parity blocks and spectral response

Status: bounded independent first pass completed, 2026-09-18. No other new branch report or repository prototype was read. Only the campaign question and common assumption/known-result stubs were consulted.

## 1. Main conclusion

Agreement on the complete homogeneous odd collision block does not, by itself, validate a full-population event implementation. A self-reciprocal exactly resonant event can be invisible in that block. Counting it twice changes even relaxation and finite-spatial-wavevector response while preserving every homogeneous odd-sector linear-response test.

This conclusion has an explicit three-mode example and a seven-mode example containing a counted pair of distinct reciprocal events. It concerns a finite event model, not an AlN calculation or an absolute material rate.

## 2. Definitions, assumptions, and event convention

Use dimensionless positive energies and finite temperature beta = 1. All mode weights are equal to one. The seven modes are ordered
\[
 (1+,2+,3+,1-,2-,3-,0).
\]
Their energies, crystal momenta in grid units, and assigned velocities are
\[
 \epsilon=(1,2,3,1,2,3,2),\quad
 q=(1,2,3,-1,-2,-3,0),\quad
 v=(1,0.7,0.4,-1,-0.7,-0.4,0).
\]
One can regard the momenta as a full seven-point reciprocal grid with reciprocal period 7. The velocity values are an assigned odd streaming operator; no actual AlN dispersion is represented. The final mode is positive-frequency and optical-like at zero wavevector, not a zero-frequency translation.

The counted physical-event list is:
\[
 (1+,2+\leftrightarrow3+),\quad
 (1-,2-\leftrightarrow3-),\quad
 (1+,1-\leftrightarrow0).
\]
Each event is included once, with its reverse population process included in its net flux. The first two are a distinct reciprocal pair. The third is a reciprocal orbit of size one because its incoming pair is unordered. Every listed event has exactly zero energy and momentum mismatch.

For one event \(a+b\leftrightarrow c\), write
\[
 \nu_i=\delta_{ia}+\delta_{ib}-\delta_{ic},\quad
 J=\Gamma[n_an_b(1+n_c)-(1+n_a)(1+n_b)n_c],
 \qquad \dot n=-\nu J.
\]
This defines the convention for the event coefficient Gamma. It is not a conversion from a material matrix element. All numerical Gamma values are one except in the deliberate double-counting comparison.

The three states in each event are distinct. In particular, \(1+\) and \(1-\) are distinct modes. This experiment does not settle the combinatorics of a genuine repeated index \(a=b\).

## 3. Independent entropy-coordinate derivation

At Bose equilibrium,
\[
 n_i^0=(e^{\epsilon_i}-1)^{-1},\quad
 w_i=n_i^0(1+n_i^0),\quad
 y_i=\delta n_i/\sqrt{w_i}.
\]
Resonance implies detailed balance. If
\(F=n_a^0n_b^0(1+n_c^0)\), differentiating the net flux gives
\[
 \delta J=\Gamma F\sum_i\nu_i\,\delta n_i/w_i.
\]
Consequently
\[
 \dot y=-Cy,\qquad
 C_{\rm event}=\Gamma F\,dd^T,\qquad d_i=\nu_i/\sqrt{w_i}.
\]
The rank-one matrix is positive semidefinite. Its action conserves the entropy-coordinate energy and momentum vectors
\[
 e_i=\epsilon_i\sqrt{w_i},\qquad p_i=q_i\sqrt{w_i},
\]
because \(d^Te=\nu^T\epsilon=0\) and \(d^Tp=\nu^Tq=0\).

These statements are derivations under the stated event convention. The numerical validation below separately differentiates the unexpanded nonlinear Bose flux.

## 4. Parity blocks and the hidden event

Let R exchange each plus/minus pair and fix mode 0. The normalized parity basis contains
\[
 |a,e\rangle=(|a+\rangle+|a-\rangle)/\sqrt2,\quad
 |a,o\rangle=(|a+\rangle-|a-\rangle)/\sqrt2,
\]
and the extra even vector \(|0\rangle\).

The counted reciprocal pair alone gives identical three-by-three even and odd collision blocks. That equality is a special property of this disjoint pair, not a general certificate for a collision implementation.

The self-reciprocal event has
\[
 d_{\rm self}=(1/\sqrt{w_1},0,0,1/\sqrt{w_1},0,0,-1/\sqrt{w_0})^T.
\]
It is entirely even. Define
\[
 C_{\rm good}=C_{\rm pair}+F_0d_{\rm self}d_{\rm self}^T,\qquad
 C_{\rm double}=C_{\rm pair}+2F_0d_{\rm self}d_{\rm self}^T.
\]
The complete odd blocks are identical exactly; the even blocks differ by a nonzero rank-one matrix.

Both operators are conservative, reciprocal, and positive semidefinite. Thus those structural checks also cannot identify the specified event count by themselves. For a fixed event list, the nonlinear full-population derivative distinguishes them. Equivalently, as two models with different event coefficients, their odd data do not identify the full operator.

## 5. Nonlinear and structural validation

The script computes a central-difference full entropy Jacobian directly from the nonlinear population collision action. It also checks a generic direction at 80 decimal digits.

| Check | Executed result |
|---|---:|
| Nonlinear equilibrium residual | \(9.61\times10^{-17}\) |
| Relative energy collision residual | \(7.48\times10^{-17}\) |
| Relative momentum collision residual | \(7.35\times10^{-17}\) |
| Reciprocal collision residual and streaming parity residual | exactly zero in stored arithmetic |
| Odd-block change from double counting | \(3.79\times10^{-35}\), analytically zero |
| Even-block change norm | 3.01487154198 |
| Full nonlinear-Jacobian relative error, step \(10^{-2}\) | \(1.09\times10^{-14}\) |
| Full nonlinear-Jacobian relative error, step \(10^{-8}\) | \(1.39\times10^{-9}\) |
| 80-digit nonlinear directional derivative, step \(10^{-20}\) | relative error \(4.30\times10^{-62}\) |
| Double-counted action error on the tested odd perturbation | \(2.33\times10^{-13}\) |
| Double-counted action error on the tested even perturbation | 0.9866883026 |

The decreasing finite-difference step eventually amplifies cancellation; it does not monotonically improve this calculation. For these distinct-mode events, the cubic terms cancel in the net Bose flux, leaving a quadratic polynomial. Central differentiation therefore has no ordinary quadratic truncation error in a fixed direction; roundoff is the relevant observed limitation.

The positive full collision eigenvalues are approximately
\[
 2.2555198183,\quad2.3721973063,\quad3.1315490300.
\]
There are four zero eigenvalues up to roundoff. The even block has two positive eigenvalues and the odd block has one. This sparse event model has accidental invariants in addition to energy and momentum. No finite-DC transport claim is made.

## 6. Finite-spatial-wavevector experiment

Use the direct linear resolvent
\[
 {\cal R}(z,k)=(zI+C+ikV)^{-1},\qquad z=0.25,\qquad V=\operatorname{diag}(v).
\]
This can be read as a Laplace response to initial data. No spatial boundary law or time-stepping approximation enters.

Parity gives \(RCR=C\) and \(RVR=-V\). Streaming therefore couples the two sectors. For an odd input, exact even-sector elimination gives
\[
 S_o=z+C_o+k^2V_{oe}(z+C_e)^{-1}V_{eo}.
\]
Identical \(C_o\) does not fix this operator when \(C_e\) differs.

For current \(j=Ve\) and energy vector e, compare
\[
 H_j=j^T{\cal R}j,\qquad H_e=e^T{\cal R}e.
\]

| k | Relative current-response change from double counting | Relative energy-response change |
|---|---:|---:|
| 0 | \(3.41\times10^{-16}\) | \(1.77\times10^{-16}\) |
| 0.03 | \(8.91\times10^{-5}\) | \(6.67\times10^{-7}\) |
| 0.1 | \(8.87\times10^{-4}\) | \(7.32\times10^{-5}\) |
| 0.3 | 0.00411092 | 0.00292475 |
| 1 | 0.00681672 | 0.04580047 |

At k = 1, the correct and double-counted energy responses are 2.59420771625 and 2.47539178079. The effect is resolved in a well-conditioned calculation, not inferred from a near-pole solve.

Independent comparison paths were:

1. Direct full seven-dimensional solve.
2. Exact parity-block elimination for the odd response.
3. The rank-one resolvent identity, with
\[
 {\cal R}_{\rm double}-{\cal R}_{\rm good}
 =-\frac{F_0{\cal R}_{\rm good}d\,d^T{\cal R}_{\rm good}}
 {1+F_0d^T{\cal R}_{\rm good}d}.
\]

The maximum parity-response relative discrepancy is \(6.25\times10^{-16}\). The maximum rank-one-update relative discrepancy is \(2.28\times10^{-14}\). Normalized direct equation residuals are at most \(7.46\times10^{-17}\); full-matrix condition numbers range from about 5.14 to 13.53.

The current-response difference scales as \(k^2\), and the energy-response difference as \(k^4\), at fixed z. The final 80-digit refinement powers are 1.999990479 and 3.999990479. At \(k=10^{-4}\), direct high-precision differences are
\[
 \Delta H_j=1.0435069522\times10^{-8},\qquad
 \Delta H_e=-1.6696111236\times10^{-15}.
\]
Their rank-one-update relative discrepancies are below \(1.41\times10^{-72}\) and \(1.38\times10^{-65}\), respectively.

The two powers are linked by conservation, rather than being fully independent evidence. Since \(Ce=0\), \(j=Ve\), and \(e^Tj=0\),
\[
 H_e=\frac{\|e\|^2}{z}-\frac{k^2}{z^2}H_j,\qquad
 \Delta H_e=-\frac{k^2}{z^2}\Delta H_j.
\]
This explains why an energy-response test can miss the error at lower spatial order than a current-response test.

## 7. Minimal three-mode example

Keep only modes \((+q,-q,0)\), energies \((1,1,2)\), velocities \((1,-1,0)\), and the event \(+q+(-q)\leftrightarrow0\). This is the smallest construction used here with a nontrivial odd mode and distinct incoming states; no universal minimality theorem is claimed.

Its odd collision block is identically zero for every event strength. If \(g=\Gamma F_0\), define
\[
 a=2g/w_q,\quad b=g/w_0.
\]
The exact current response is
\[
 H_j(z,k)=\frac{2w_q}
 {z+k^2(z+b)/[z(z+a+b)]}.
\]
At z = 0.25 and k = 0.3, doubling Gamma changes the response from 3.56730083479 to 3.59330364470, a relative change of 0.0072892114. Direct solves and the scalar formula agree to \(1.25\times10^{-16}\).

Thus the seven-mode result is not dependent on numerical complexity or on the additional reciprocal pair.

## 8. Limitations and next test

Established within this event convention: a complete homogeneous odd-sector check can fail to detect an erroneous reciprocal-orbit count, and finite-k streaming can expose the error.

Not established: material matrix elements, selection-rule permission for the chosen channels in AlN, physical absolute rates, repeated-index statistical factors, nonuniform grid weights, energy broadening, continuum limits, or a material transport coefficient. Generic even-sector normalization errors can be tested by the same principle, but this calculation specifically constructs a reciprocal-orbit fixed-point example.

The next discriminating implementation test is to include both a non-self-reciprocal orbit and a self-reciprocal orbit in a full-population nonlinear derivative test, then separately test generic even and odd perturbations. A genuine repeated-index event should be a separate test with its combinatorial convention stated explicitly.

## 9. Reproducibility and status

Artifacts owned by this branch:

- experiments/C_parity_events.py
- experiments/C_parity_results.json
- branches/C-parity-spectrum.md

Executed command:

~~~powershell
py 'C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260917-232509-aln-physical-events\experiments\C_parity_events.py'
~~~

Environment: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.3.0. The run exited with code 0. All ten recorded assertion checks passed. No repository code was imported or changed. The initial design write hit a usage-limit rejection; after the explicit resume, the report, script, and computation were completed using the authorized elevated shell.

DERIVED UNDER ASSUMPTIONS: entropy-coordinate rank-one event action, exact odd invisibility, parity elimination, and conservation response identity.
OBSERVED NUMERICALLY: the residuals, response differences, and refinement results reported above.
FAILED VALIDATION STRATEGY: homogeneous odd-sector agreement alone.
UNRESOLVED: independent audit and real-material event normalization/selection rules. No novelty claim is made.

