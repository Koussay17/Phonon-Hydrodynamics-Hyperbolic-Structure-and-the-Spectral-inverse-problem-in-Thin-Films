# Second-generation D: exact elimination, static error, and pole conditioning

Status: completed bounded experiment, 2026-09-17. Computational evidence, with accompanying elementary derivations; this branch does not claim a proof audit or independent numerical-analysis review.

Read: checkpoint-01.md, branches/A-projection.md, and the research protocol. No red-team report was read and no reviewer was contacted. Frozen artifacts were not modified.

## 1. Question and scope

Test whether exact Schur elimination reproduces the full response with noncommuting streaming, whether A5 bounds the static operator error, and whether a small absolute operator error guarantees a good slow response. Separate model error from floating-point error near a pole.

The model is one spatial Fourier mode:
\[
 \dot y=-(C+ikV)y+s,\qquad C=C^T\ge0,\quad V=V^T.
\]
All numerical parameters are dimensionless relative to a fast collision scale. The spatial setting is periodic or infinite, so no kinetic wall law is involved. Resolvents include arbitrary forcing and initial data through the right-hand side. No time integrator, spatial discretization, or continuum limit is used.

These are abstract finite collision operators. They do not identify AlN rates, viscosity, microscopic scattering events, or a physical hydrodynamic window.

## 2. Independent calculation and implementation

Partition the full linear system into retained and eliminated coordinates:
\[
 (z+A+iU)a+Yw=f_P,\qquad Xa+(z+D+iW)w=f_Q.
\]
Solving the second equation and substituting into the first gives
\[
 S=z+A+iU-YR_QX,\quad
 a=S^{-1}(f_P-YR_Qf_Q),\quad
 w=R_Q(f_Q-Xa).
\]
At fixed k, both off-diagonal blocks contain +i streaming. Thus Y is generally not the adjoint of X.

The program constructs the full matrix in a mode basis with diagonal velocities and solves all right-hand sides simultaneously. A separate block calculation reconstructs the full inverse, including the fast-sector source term. This tests the entire response map rather than one chosen forcing.

For A5, the error is evaluated through
\[
 S-S_0=YR_Q(z+iW)D^{-1}X.
\]
This avoids cancellation when subtracting almost equal Schur matrices. The program also records the discrepancy from direct subtraction. The tiny cubic errors below are evaluated with this identity; double-precision subtraction does not independently resolve a difference of order \(10^{-19}\).

## 3. Parameter sweep and exact-elimination results

Seed: 20260917. Five eight-dimensional systems retain three coordinates. Each has an exact invariant before roundoff, a constant positive energy-like vector in the mode basis, and eliminated eigenvalues
\[
 (0.5,0.7,1.1,1.6,2).
\]
The collision matrix is constructed as
\[
 C_{\rm block}=
 \begin{pmatrix}
 H+K^TD^{-1}K&K^T\\K&D
 \end{pmatrix},\qquad
 H=\operatorname{diag}(0,0.025,0.08),
\]
with the first column of K zero. The random coupling is scaled by 0, 0.01, 0.2, 1, or 5. An orthogonal rotation places the invariant along the constant vector. Mode velocities are
\[
 (-1.7,-1.1,-0.6,-0.2,0.2,0.6,1.1,1.7).
\]
The collision and velocity matrices do not commute. These generic models do not additionally impose phonon collision parity or microscopic event realizability.

The 600 cases use
\[
 \Re z\in\{10^{-7},10^{-3},0.1,1\},\quad
 \Im z\in\{-2,-0.3,0,0.4,1.5\},\quad
 k\in\{0,10^{-4},0.02,0.2,1,3\}.
\]
All norms are spectral 2-norms. The normalized reconstruction residual is
\[
 \frac{\|M R_{\rm reconstructed}-I\|}
 {\|M\|\|R_{\rm reconstructed}\|+1},\qquad M=zI+C+ikV.
\]

| Diagnostic | Executed result |
|---|---:|
| Maximum normalized reconstruction residual | \(5.13\times10^{-16}\) |
| Maximum normalized direct-solve residual | \(1.96\times10^{-16}\) |
| Maximum projected full-versus-Schur relative discrepancy, 590 cases with condition number below \(10^6\) | \(6.90\times10^{-12}\) |
| Maximum discrepancy across all 600 cases | \(2.51\times10^{-8}\) |
| Full-matrix condition number in that worst case | \(2.10\times10^9\) |
| Maximum A5 actual-error/bound ratio | 0.5202674411 |
| Maximum A6 actual-error/bound ratio | 0.9491041431 |
| Maximum fast-source correction/bound ratio | 0.6714499088 |
| Minimum static-generator Hermitian eigenvalue | \(-4.16\times10^{-17}\) |

The largest response discrepancy occurs at strong coupling, k = 0, z = \(10^{-7}\), with a normalized reconstruction residual of \(7.63\times10^{-17}\). Its forward discrepancy is consistent with conditioning amplification; it is not evidence against the exact identity.

A3's maximum computed ratio is \(1+1.8\times10^{-15}\). The maximum computed \(\sigma\|S^{-1}\|\) and \(\sigma\|S_0^{-1}\|\) is \(1+4.13\times10^{-10}\). These tiny violations must not be described as exact compliance: the implemented rotated collision matrices and solves have finite precision. The invariant residuals and complete collision spectra are preserved in the JSON output.

A5 is also attained exactly in the scalar calibration D = K = 1, B = W = 0, z = 0.3:
\[
 |S-S_0|=\frac{0.3}{1.3}=0.23076923076923078.
\]
Thus the absence of saturation in the random sweep does not imply the general bound can simply be reduced by its observed maximum ratio.

## 4. Controlled K = 0 scaling

With K = 0, take \(z=(0.3+0.4i)k\), and decrease k from \(10^{-1}\) to \(10^{-6}\) in half-decade steps. The error behaves as
\[
 \|S-S_0\|\simeq0.4935997\,k^3.
\]
The final measured convergence power is 3.000000469. At \(k=10^{-6}\), the identity-based operator error is \(4.9360\times10^{-19}\), while A5 is \(5.0139\times10^{-18}\).

This supports the expected extra small factor in the aligned case. It is a finite-dimensional asymptotic calculation and does not establish a uniform continuum gap.

## 5. Decisive poorly aligned coordinate with noncommuting streaming

Use
\[
 C_\varepsilon=\begin{pmatrix}1+\varepsilon&1\\1&1\end{pmatrix},
 \quad V=\begin{pmatrix}1&0\\0&-1\end{pmatrix},
 \quad E=(1,0)^T,\qquad\varepsilon>0.
\]
C is positive definite, and \(\|[C,V]\|=2\). An independent zero collision block can be appended if an exact invariant is desired. This two-coordinate demonstration by itself is not a microscopic phonon model.

Here D = K = 1, B = 0, U = k, W = -k. Therefore
\[
 S=z+1+\varepsilon+ik-\frac1{1+z-ik},\qquad
 S_0=z+\varepsilon+ik.
\]
The full decay eigenvalues are roots of
\[
 \lambda^2-(2+\varepsilon)\lambda+
       (\varepsilon+k^2-i\varepsilon k)=0.
\]
The code compares a stable quadratic-root evaluation with an independent arbitrary-precision dense eigensolver. It also compares dense resolvent solves with the scalar Schur response. Precision is 90 decimal digits; epsilon runs from \(10^{-2}\) to \(10^{-12}\).

At k = 0 the true slow decay is
\[
 \lambda_{\rm slow}=\varepsilon/2+O(\varepsilon^2),
\]
whereas static elimination predicts epsilon.

More strongly, with \(k=\varepsilon\),
\[
 \lambda_{\rm slow}
 =\varepsilon/2+3\varepsilon^2/8-i\varepsilon^2/2+O(\varepsilon^4),
 \qquad \lambda_0=\varepsilon+i\varepsilon.
\]
Thus the static approximation has both a factor-two decay error and an incorrect first-order propagation frequency.

| epsilon, with k = epsilon | Static/actual decay rate | Static complex-eigenvalue relative error | Causal static-response relative error |
|---|---:|---:|---:|
| \(10^{-2}\) | 1.985111878 | 2.224889563 | 0.675318741 |
| \(10^{-6}\) | 1.999998500 | 2.236066859 | 0.676027724 |
| \(10^{-12}\) | 1.9999999999985 | 2.236067977499 | 0.676027791337 |

The response column uses \(z=0.1\varepsilon>0\), \(f_P=1\), \(f_Q=0\); it lies inside the causal half-plane. At epsilon = \(10^{-12}\),
\[
 \|S-S_0\|=1.00499\times10^{-12},\qquad
 \text{A5}=1.10000\times10^{-12},
\]
but the response error is still 67.6%. In fact,
\[
 \lim_{\varepsilon\to0}
 \frac{|S_0^{-1}-S^{-1}|}{|S^{-1}|}
 =\sqrt{\frac{101}{221}}=0.6760277913\ldots.
\]
The controlling perturbation quantity \(|S_0^{-1}|\ |S-S_0|\) approaches the same nonzero constant. There is no contradiction with A5: A5 controls an absolute operator error, and the inverse grows like \(1/\varepsilon\).

The largest dense-eigenvalue versus quadratic-root relative difference in the 90-digit run is \(1.67\times10^{-79}\); the largest dense versus Schur response difference is \(3.62\times10^{-80}\). The persistent approximation error is not explained by ordinary floating-point loss.

The mechanism is visible without numerical fitting:
\[
 S=\varepsilon+2z+O((|z|+|k|)^2).
\]
The neglected fast component changes the leading time-derivative coefficient and cancels the leading streaming term. Using \(\varepsilon/2\) instead of the static eigenvalue reduces the relative eigenvalue error at \(k=\varepsilon=10^{-12}\) to \(1.25\times10^{-12}\). This is a correction for this example, not a validated general replacement closure.

## 6. Pole error versus roundoff

A separate aligned two-mode model uses
\[
 C=\operatorname{diag}(0,1),\quad
 V=\begin{pmatrix}0&1\\1&0\end{pmatrix}
\]
in the collision basis. In a mode basis the same model has diagonal velocities and
\(C=\frac12\begin{pmatrix}1&-1\\-1&1\end{pmatrix}\), so its invariant is positive and collision/streaming do not commute.

For \(k=0.01\),
\[
 S=z+\frac{k^2}{1+z},\quad S_0=z+k^2,\quad
 \lambda_{\rm slow}=\frac{1-\sqrt{1-4k^2}}2
 =0.0001000100020005001442\ldots.
\]
The static rate is 0.0001 and the pole displacement is \(1.00020005\times10^{-8}\).

Compare responses at \(z=-\lambda_{\rm slow}+\delta\), using 110-digit dense solves as reference for the exact binary z and k supplied to the floating-point computation. The float64 errors include matrix assembly as well as solution arithmetic. Actual rounded distances to the pole are recorded.

| Intended delta | Static response relative error | Direct float64 response relative error |
|---|---:|---:|
| \(10^{-6}\) | 0.0100020 | \(8.22\times10^{-15}\) |
| \(10^{-8}\) | 4999.25 | \(1.21\times10^{-12}\) |
| \(10^{-12}\) | 1.00010 | \(7.77\times10^{-9}\) |
| \(10^{-16}\) | 1.00000001 | \(1.34\times10^{-4}\) |
| \(10^{-20}\) | 1.000000000001 | 0.253248 |

At delta near \(10^{-8}\), the frequency is close to the static pole: the approximation error is enormous while the direct computation remains accurate. Closer to the exact pole, the static approximation misses the large exact response. Still closer, even the exact full-system float64 calculation loses accuracy despite a tiny normalized residual.

These listed z values are negative. They are outside A5's stated \(\Re z\ge0\) domain, and are not counterexamples to that theorem. A scalar resolvent estimate remains algebraically valid for \(z>-1\), but positive-half-plane inverse bounds cannot be carried over unchanged. The poor-coordinate experiment above separately establishes nonvanishing response error within the stated half-plane.

The largest dense-versus-Schur discrepancy at 110 digits is \(3.73\times10^{-96}\). SciPy emitted three expected ill-conditioning warnings for the closest cases. Float64 SVD condition estimates plateau in the most singular cases and should not be interpreted as accurate condition numbers there. No second precision-refinement run was performed; the independent dense/scalar high-precision agreement is the available check.

As a causal aligned control, at \(z=k^2\) the relative static-response error decreases from 0.004950495 at k = 0.1 to approximately \(5.0\times10^{-9}\) at k = \(10^{-4}\). Failure is therefore not inevitable for every local closure.

## 7. Executed checks and limits

Executed in the saved run:

- Full direct inverse versus reconstructed Schur inverse, for all right-hand sides in 600 cases.
- Normalized equation residuals and matrix condition estimates.
- Collision spectra and invariant residuals.
- A3, A5, A6, source-correction bounds, static accretivity, and inverse-coercivity diagnostics.
- K = 0 asymptotic refinement across eleven k values.
- Stable quadratic eigenvalues versus independent 90-digit dense eigenvalues.
- Direct 90-digit resolvent versus scalar Schur elimination for the poor coordinate.
- Float64 direct/scalar responses versus 110-digit dense/scalar responses near a pole.
- A causal aligned control and the deliberately wrong-sign static rate from replacing Y by the adjoint of X.

These are quantitative comparisons recorded by the script; the original executed script does not use assertion-based pass/fail tests. The sampled checks are evidence, not a proof of the candidate bounds for all admissible matrices.

No timestep, mesh, finite-domain, or boundary refinement is applicable to these direct finite-matrix calculations. No material data, continuum convergence, or physical sample model has been tested.

## 8. Strongest justified conclusion and next discriminating test

Exact Schur elimination and the stated operator-error structure survive this bounded computational test. Static elimination can nevertheless fail at leading dynamical order when the retained coordinate has order-one coupling to the eliminated sector. Noncommuting streaming can add a wrong leading propagation frequency. A vanishing absolute A5 bound alone does not certify a vanishing relative response error.

For a proposed material closure, the next discriminating calculation is to evaluate the actual slow-subspace residual K, the time-derivative correction \(K^TD^{-2}K\), the streaming correction, and the observable-specific inverse perturbation quantity over the intended frequency and wave-vector range. That requires collision action beyond diagonal lifetimes.

## 9. Reproducibility and ledger

Files:

- experiments/second_generation_schur.py
- experiments/second_generation_results.json (all individual cases)
- experiments/second_generation_console.txt (summary, asymptotic tables, and high-precision cases)

Executed command from the campaign directory:

~~~powershell
py 'C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260916-210742-aln-spectral-closure\experiments\second_generation_schur.py'
~~~

Environment: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1, mpmath 1.3.0. The completed run exited with code 0.

- DERIVED: scalar characteristic polynomials, joint poor-coordinate expansion, and its nonzero causal relative-response-error limit.
- OBSERVED NUMERICALLY: the quantitative comparisons above.
- FAILED: static closure of the poorly aligned coordinate; ordinary double precision extremely close to a pole.
- INFRASTRUCTURE: the initial script write hit a usage-limit rejection; the resumed authorized write and execution succeeded. The default shell had a sandbox ACL initialization failure, so the run used the approved elevated shell. A later optional script-patch attempt hit the same sandbox ACL issue; the report describes only the already executed script.
- UNRESOLVED: general proof audit, independent numerical-analysis audit, real-AlN collision action and slow-space alignment, continuum limits, and boundary data.
