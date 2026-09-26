# Second-generation retained resolvent and memory check

Completed bounded four-mode calculation. Uses only C's exact first-pass model, common campaign assumptions, and the existing repository operator. No other branch or PI result was used. No repository or first-pass file was edited.

## 1. Question and model

Test whether exact elimination reproduces the retained response, whether the equal-rate case closes without memory, and how conserved zero modes affect error interpretation.

Use the first-pass model with positive energies (log6,log6,log2,log3), equilibrium variances
\[
 W=\operatorname{diag}(6/25,6/25,2,3/4),
\]
events p1 <-> a+b and p2 <-> a+b, common equilibrium factor 3/5, and prefactors (2,2) or (3,1). The original population generator is L=MW^{-1}; the symmetric entropy-coordinate generator is C=W^{-1/2}MW^{-1/2}. This is a stipulated finite scalar-population model, not a full quantum kinetic generator.

Retain the normalized parent sum and the two daughters:
\[
 E=\big((1,1,0,0)^T/\sqrt2,\ (0,0,1,0)^T,\ (0,0,0,1)^T\big),
\]
and eliminate q=(1,-1,0,0)^T/sqrt2. Write
\[
 A=E^TCE,\quad K=q^TCE,\quad D=q^TCq.
\]
Both models have the same rank-one A, with nonzero eigenvalue 47/5, and D=5. Equal rates give K=0. Unequal rates give
\[
 K=(5/2,-\sqrt{3/2},-2),\quad \|K\|^2=47/4,\quad K^TK=(5/4)A.
\]

## 2. Exact resolvent and memory formulas tested

For Re z>0,
\[
 R_P(z)=E^T(zI+C)^{-1}E
       =[zI+A-K^T(z+D)^{-1}K]^{-1}.
\]
The hidden-source response is
\[
 E^T(zI+C)^{-1}q=-R_P(z)K^T/(z+D).
\]
The time equation for retained a and hidden initial value b(0) is
\[
 \dot a(t)+Aa(t)-\int_0^t K^Te^{-5(t-s)}K a(s)\,ds
       =-K^Te^{-5t}b(0).
\]
The initial-state term is part of the exact elimination.

For equal rates, K=0 and the projected generator A is an exact closure. For unequal rates, the Schur operator is
\[
 S(z)=zI+A-\frac{(5/4)A}{z+5}.
\]
The comparison called “local” in the script is the original projected model (zI+A)^{-1}, which drops memory. It is not a separately corrected static Schur model.

The script also checks the population-total transformation:
\[
 B(zI+L)^{-1}H
 =\Sigma^{1/2}R_P(z)\Sigma^{-1/2},
\]
with the first-pass B, H, and Sigma=BWB^T.

## 3. Zero modes are kept explicitly

Both full models have exactly two zero modes: the entropy-weighted forms of the conserved population combinations
\[
 (1,1,1,0)^T,\qquad (1,1,0,1)^T.
\]
Energy is one linear combination. These two directions are retained; there is no hidden zero mode because D=5.

Computed positive full eigenvalues:

| Model | Positive eigenvalues |
|---|---|
| Equal rates | 5, 9.4 |
| Unequal rates | 3.12691762912, 11.27308237088 |

Numerically computed zero eigenvalues are of order 10^-16. Inverse calculations use Re z>0 throughout. Neither C nor A is inverted at z=0, and no artificial regularization or deletion of an invariant is used.

Let r be the unit relaxing eigenvector of A. The exact unequal-rate scalar response is
\[
 H_r(z)=r^TR_P(z)r
       =\frac{z+5}{z^2+(72/5)z+141/4}.
\]
The projected local response is H_0(z)=1/(z+47/5). The other two retained responses are exactly 1/z.

## 4. Executed resolvent checks

Twenty double-precision cases use both models, z=sigma+i omega with sigma in {0.001,0.1,1}, omega in {0,0.4,3}, and z=10^-6. Full dense inversion, block Schur inversion, hidden-source response, and the population-coordinate transformation are evaluated separately.

| Diagnostic | Maximum |
|---|---:|
| Full-retained versus Schur relative 2-norm error | 3.65e-10 |
| Normalized full equation residual | 9.59e-17 |
| Hidden-source absolute discrepancy | 9.31e-11 |
| Population-coordinate relative discrepancy | 3.54e-10 |
| Full-matrix condition number | 1.13e7 |

The largest forward discrepancies occur near the conserved poles. They must not be equated with the much smaller equation residual.

Six additional 80-decimal-digit calculations use both models at positive real z=10^-3, 10^-6, 10^-10. The maximum full-retained versus Schur relative Frobenius discrepancy is 1.69e-71. The independently evaluated relaxing-channel scalar formula agrees to 1.06e-71 or better. Conserved residues z R_P P_0=P_0 agree to 3.58e-71, where P_0 projects onto the two retained invariants.

### Conserved modes can mask a closure error

For unequal rates, as z decreases along the positive real axis,
\[
 \frac{|H_r-H_0|}{|H_r|}\longrightarrow\frac14.
\]
Yet the relative full retained resolvent error tends to zero, because its norm is dominated by 1/z:

| Positive real z | Relaxing-channel relative local error | Full retained relative spectral error |
|---|---:|---:|
| 10^-3 | 0.2499234224 | 3.54427e-5 |
| 10^-6 | 0.2499999234 | 3.54610e-8 |
| 10^-10 | 0.249999999992 | 3.54610e-12 |

The full spectral-error column uses the exact orthogonal decomposition into the two conserved directions and the one relaxing direction. It is not a floating-point cancellation artifact. Equal-rate local closure agrees to the high-precision arithmetic floor.

## 5. Time-domain memory and initial data

Use y(0)=Er+q, so both retained and hidden coordinates are initially present. Since this is a linearized calculation, this direction may be scaled arbitrarily small. Evaluate direct matrix exponentials at t=0.05,0.2,1, and independently integrate the memory convolution using adaptive vector quadrature with absolute and relative tolerances 10^-12.

The maximum normalized memory-equation residual is 8.01e-17; the largest quadrature error estimate is 2.85e-14. This checks the convolution identity against the known full trajectory; it is not a separate time-stepping solution of a memory equation.

For unequal rates, omitting the hidden-initial-state term leaves absolute equation residuals:

| t | Residual after omission |
|---|---:|
| 0.05 | 2.66959458563 |
| 0.2 | 1.26102719163 |
| 1 | 0.023096518671 |

They match the norm of K^T exp(-5t)b(0). Equal rates have no such term because K=0.

No timestep or spatial discretization is used. The high-precision and algebraic checks distinguish model error from finite-precision loss; no convergence claim for a material grid follows.

## 6. Repository action and resonance-rounding caveat

The optional check imported C:/Users/Koussay/these/src/collision_events.py read-only, using the parent-first events (0,2,3),(1,2,3) and the two prefactor pairs.

Source SHA-256 before and after:
~~~
77196934d0e3ee209e7121f31655af48e074864d7fff5668cfb6f07b1ea99e6c
~~~

The repository entropy action differs from the exact-rational model by 2.99e-16 relatively for equal rates and 9.15e-17 for unequal rates.

The energies passed were separately evaluated binary64 logarithms. Although log6=log2+log3 mathematically, their stored values have exact binary mismatch
\[
 \delta=-1.1102230246251565404\ldots\times10^{-16}
\]
for both events. The operator accepts this under its default 10^-12 roundoff tolerance.

This is not exact resonance of the stored numbers. At 80 digits, using those exact binary inputs, the true nonlinear Bose Jacobian has relative Frobenius asymmetry 1.29e-16 and 1.31e-16 in the two cases; the Bose base-state RHS norms are 4.21e-16 and 4.32e-16. The repository's symmetric positive factor construction is therefore an approximation to that stored-input nonlinear Jacobian, while matching the ideal resonant toy model at ordinary rounding accuracy.

No uniform response-error bound follows from accepting a resonance tolerance, especially near conserved poles. The present comparison supports only these moderate-energy finite cases. Genuine repeated-index support, material normalization, and export completeness were not tested.

## 7. Artifacts, status, and limits

Owned artifacts:

- experiments/second_memory_check.py
- experiments/second_memory_check.json
- branches/second-memory-check.md

Reproduce:
~~~powershell
py -B D:\ResearchLab\orchestration\campaigns\20260925-113533-aln-degenerate-collision-action\experiments\second_memory_check.py
~~~

Environment: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.3.0. Script hash is recorded in the JSON. Run exited with code 0; all eleven recorded checks passed. Bytecode writes were disabled. The preliminary save was interrupted by a usage-limit rejection; after the explicit resume the report, script, and calculation were completed.

Strongest supported conclusion: exact elimination restores the retained response through memory and hidden-initial-state information. Equal-rate closure is the positive control. Unequal-rate block compression remains inaccurate in its relaxing channel even when a conserved-mode-dominated relative norm appears small.

This is a synthetic scalar event calculation. No material calculation, AlN transport result, quantum generator, or material export is validated. Independent root derivation and subsequent audit are separate from this report.


## 8. Portable reproduction and finalized script

The script now accepts --repo-root and --output-dir. The same bounded computation was rerun using both options on 2026-09-26; all eleven checks passed again, with unchanged numerical results. The JSON records the resolved paths, interpreter executable, and finalized source hash.

Final script SHA-256:
~~~
3f3b83eecac99a2556475ba60ee3a07769ebe19f744c7df13c015895cf1915cb
~~~

Interpreter used:
~~~
C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe
~~~

Portable requirements: Python with NumPy, SciPy, and mpmath; tested versions are Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, and mpmath 1.3.0. SymPy, the first-pass script, external data, and network access are not required. Provide a writable output directory and a repository root containing src/collision_events.py.

Example command, with paths replaced for the target system:
~~~text
python -B second_memory_check.py --repo-root "/path/to/repository" --output-dir "/path/to/results"
~~~

The complete eleven-check reproduction requires the repository source hash pinned in Section 6. A different hash is reported and the repository comparison is skipped; a missing source file is an error. The default repository root remains the original local checkout, so another machine should explicitly supply --repo-root. The default output directory is the directory containing the script. Other dependency versions were not tested.

The finalization changed CLI/provenance handling only; no scientific model, sampling range, test threshold, or repository implementation was changed.

