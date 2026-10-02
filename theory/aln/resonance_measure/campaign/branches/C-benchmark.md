# C: analytically solvable resonance-measure benchmark

**Independent first pass complete, 2026-09-28; campaign audit pending.** Read only the question and assumptions. No new peer reports, material data, or material computation. The experiment separates exact-resonance integration, finite Gaussian width, mesh resolution, and invariant tests.

## 1. Model, representation and analytic reference

Use \(q\in[-\pi,\pi)\) with normalized measure \(dq/(2\pi)\), parent momentum zero, and daughter momenta \(q,-q\). Momentum conservation is exact. In fixed dimensionless energy units,
\[
E_p=2+d,\quad E_a(q)=1+0.3\cos q,\quad E_b(k)=1+0.7\cos k,
\qquad \Delta(q)=d-\cos q.
\]
All tested mode energies are strictly positive. The dispersions and the final positive weight
\[
w(q)=1+0.2\cos(2q)\ge0.8
\]
are reciprocal and smooth. The weak-form test affinity is \(\phi(q)=1+0.4\sin q\). Set \(\beta=1\), with Bose occupation \(n(E)=1/(\exp(\beta E)-1)\).

This is a scalar kinetic-integration benchmark under an assumed weak-coupling population model. Energies and occupations are evaluated from the analytic functions at every quadrature node; no mesh interpolation, surrogate energy, or occupation-projection closure is used. The parent is fixed, so this is not a complete material collision operator.

For \(|d|<1\), there are two roots \(q_\pm=\pm\arccos d\), with
\(s=|\Delta'(q_\pm)|=\sqrt{1-d^2}\). Changing variable locally from \(q\) to \(\Delta\) gives
\[
I[f]=\int {dq\over2\pi}w(q)\delta(\Delta(q))f(q)
=\sum_{\pm}{w(q_\pm)f(q_\pm)\over2\pi s}.
\]
Writing \(W=0.8+0.4d^2\), independently derived targets are
\[
I[1]={W\over\pi s},\qquad
I[\phi^2]={W(1+0.16s^2)\over\pi s}.
\]
The root quadrature has positive weights and exactly resonant event energies. Accordingly it preserves the event energy invariant and the Bose forward/reverse identity, and its quadratic weak form is nonnegative.

**Regular case:** \(d=0.6,\ s=0.8\):
\[
I[1]=1.18/\pi=0.37560566569687304,\quad
I[\phi^2]=0.41406768586423287.
\]
The root implementation matched both targets to the displayed floating-point precision. Evaluated energy-affinity and equilibrium-flux defects were zero in this run.

**Decisive negative control:** retain these exact nodes and positive weights but omit \(1/|\Delta'|\). Both integrals are underestimated by **20%**, while all the same energy, equilibrium and positivity tests still pass. This explicitly falsifies using invariant preservation alone as evidence for continuum convergence.

## 2. Gaussian integration and separate error sources

Use the standard-deviation convention
\[
\delta_\sigma(z)={e^{-z^2/(2\sigma^2)}\over\sqrt{2\pi}\sigma}.
\]
A shifted periodic trapezoidal rule uses
\(q_j=-\pi+(j+\eta)h,\ h=2\pi/N\), and weight \(\delta_\sigma(\Delta_j)/N\).
Tested \(N=32,64,\ldots,2048\), offsets \(\eta=0,0.37\), and six widths from 0.2 to 0.00625.

A separate adaptive integration of the even angular integrand on \([0,\pi]\) supplies the continuum Gaussian reference. Subinterval breakpoints are determined from \(\cos q=d+k\sigma\); no uniform mesh is used in that reference. Quadrature tolerances were absolute \(2\,10^{-13}\), relative \(2\,10^{-11}\), with reported estimator outputs saved. These estimators are not rigorous bounds.

At fixed \(\sigma=0.025,\eta=0.37\):

| \(N\) | \(\sigma/(hs)\) | Relative mesh error against continuum Gaussian mass |
|---:|---:|---:|
| 32 | 0.159 | \(+1.92\,10^{-1}\) |
| 64 | 0.318 | \(+1.74\,10^{-1}\) |
| 128 | 0.637 | \(-5.09\,10^{-4}\) |
| 256 | 1.273 | \(+1.26\,10^{-13}\) |
| 512 | 2.546 | \(-7.77\,10^{-16}\) |

Across all six regular widths and both offsets, the largest mass discrepancy for \(N=2048\) against the adaptive Gaussian reference was \(8.89\,10^{-16}\). This observed agreement concerns the **broadened** integral, not the exact resonance measure.

The broadened mass is the Gaussian convolution of
\[
\rho(t)={0.8+0.4t^2\over\pi\sqrt{1-t^2}}\quad(-1<t<1).
\]
For fixed interior \(d\), its local small-width expansion gives
\[
I_\sigma[1]-I[1]
={\sigma^2\over2}\rho''(d)+O(\sigma^4),\qquad
{\rho''(d)\over2}={0.8+d^2\over\pi(1-d^2)^{5/2}}.
\]
At \(d=0.6\), the predicted coefficient is 1.126829431070548; the measured coefficient at \(\sigma=0.00625\) was 1.127380151702084.

| \(\sigma\) | Relative continuum Gaussian mass bias |
|---:|---:|
| 0.05 | \(7.7510\,10^{-3}\) |
| 0.025 | \(1.8899\,10^{-3}\) |
| 0.0125 | \(4.6967\,10^{-4}\) |
| 0.00625 | \(1.1725\,10^{-4}\) |

The shrinking-width bias follows the independently derived quadratic behavior once the Gaussian is resolved. The condition involving \(hs\) is a local diagnostic for these regular roots, not a universal width prescription.

## 3. Energy and Bose-equilibrium defects

Define the exact-energy quadratic defect
\[
Q_\sigma[E]=\int {dq\over2\pi}w\delta_\sigma(\Delta)\Delta^2.
\]
It is positive for a finite-width Gaussian. Thus positive quadrature weights give a PSD weak form, but true energy is not in its kernel.

For the assumed nonlinear reversible event, let
\[
F_\beta=n_p(1+n_a)(1+n_b)-(1+n_p)n_an_b.
\]
The stable identity used in the code is
\[
F_\beta=(1+n_p)n_an_b\,\operatorname{expm1}(-\beta\Delta).
\]
It follows directly that \(F_\beta\) has the opposite sign to \(\Delta\). Hence the physical energy drift of the broadened event model is
\[
\dot{\mathcal E}_\sigma
=-\int {dq\over2\pi}w\delta_\sigma(\Delta)\Delta F_\beta>0
\]
for every positive width in this example. This is an equilibrium heating defect with no sign-cancellation loophole. Nonlinear Bose equilibrium is not stationary.

For positive forward and reverse factors \(A,R\), the event entropy production is proportional to \((A-R)\log(A/R)\ge0\). At the tested Bose distribution \(\log(A/R)=-\beta\Delta\), so that entropy production equals \(\beta\dot{\mathcal E}_\sigma\). Entropy production can therefore remain nonnegative while energy conservation and Bose stationarity fail. A symmetric weak form should not be identified with a true equilibrium linearization when the proposed equilibrium drifts.

At the regular root, \(B=n_p(1+n_a)(1+n_b)=0.1527425870300458\). Independent expansions predict
\[
Q_\sigma[E]\sim I[1]\sigma^2,\qquad
\dot{\mathcal E}_\sigma\sim\beta B I[1]\sigma^2,
\]
and the integrated absolute event flux scales as
\(\beta B I[1]\sqrt{2/\pi}\,\sigma\).

| \(\sigma\) | Equilibrium energy drift |
|---:|---:|
| 0.05 | \(1.46036\,10^{-4}\) |
| 0.025 | \(3.60136\,10^{-5}\) |
| 0.0125 | \(8.97392\,10^{-6}\) |
| 0.00625 | \(2.24166\,10^{-6}\) |

At the smallest width, the heating divided by its predicted leading term was 1.0002699703. Also \(Q_\sigma[E]=1.46773\,10^{-5}\), integrated net event flux \(=1.32265\,10^{-6}\), and absolute event flux \(=2.86148\,10^{-4}\). Mesh convergence leaves these finite-width defects intact.

## 4. Empty, critical and unresolved cases

### Empty resonance: \(d=1.1\)

The exact measure is zero; the energy mismatch is at least 0.1. Gaussian masses at widths \(0.1,0.05,0.025\) were approximately \(0.303643,\ 0.0760505,\ 1.98142\,10^{-4}\). All produce positive equilibrium heating.

Since the normalized integral of \(w\) is one,
\[
0<I_\sigma[1]\le
{1\over\sqrt{2\pi}\sigma}\exp[-0.1^2/(2\sigma^2)].
\]
This analytic bound establishes decay of the leakage. At width 0.0125 the adaptive estimate was approximately \(7.6\,10^{-15}\), below the nominal absolute tolerance; no high relative accuracy is claimed for that tiny value.

### Critical resonance: \(d=1\)

Here \(\Delta(q)=1-\cos q\sim q^2/2\), so the regular-root formula is singular and no finite regular delta-pullback measure exists. Rescaling \(q=\sqrt{\sigma}\,t\) gives
\[
I_\sigma[1]\sim {C\over\sqrt{\sigma}},\qquad
C={1.2\,\Gamma(1/4)\over2^{7/4}\pi^{3/2}}
=0.23229240994074474.
\]
The same rescaling gives \(Q_\sigma[E]/(\sigma^2I_\sigma[1])\to1/2\).

At widths \(0.1,0.025,0.00625,0.0015625\), the normalized masses
\(I_\sigma[1]\sqrt{\sigma}/C\) were
\(0.981432,\ 0.995103,\ 0.998760,\ 0.999689\).
The final energy-defect ratio was 0.49968935. The mass itself rose to 5.87476. This is the expected singular limit, not convergence to a finite resonance integral.

### Missed regular roots near criticality

At \(d=0.9999\), the exact roots are \(\pm0.01414225348\), their slope magnitude is 0.01414178207, and the exact mass is 27.00836416. A simple sign-change scan with offset 0.37 missed **both** roots for \(N=16,32,64,128\), returning zero. For \(N\ge256\), it found both, with mass errors below \(3\,10^{-13}\) relative in the tested scans.

A scan finding no sign change does not certify an empty resonance set. Exact invariant tests on an empty output would be vacuous.

### Order of limits on a fixed mesh

At \(d=0,N=64\), the exact continuum mass is \(0.8/\pi\). A mathematically root-aligned mesh has mass asymptotic to
\(1.6/[64\sqrt{2\pi}\sigma]\); a shifted mesh without roots has an exponentially vanishing sum as \(\sigma\to0\).

At width 0.0002, the aligned mass was **49.86778505**, while the shifted computation underflowed to **0**. The exact continuum mass is about 0.254648. This falsifies shrinking width on an unresolved fixed mesh. The aligned analytical limit assumes exact grid-root alignment; floating-point root residuals would intervene at extremely smaller widths.

## 5. Conclusion, scope and reproduction

**Qualified candidate:** complete regular-root enumeration, positive coarea weights, and occupations evaluated at those exact-energy nodes pass the analytic integration and invariant checks in this one-dimensional benchmark. This does not establish a convergent higher-dimensional material operator or a closed discretization on fixed occupation nodes.

**Failed alternatives preserved:** missing Jacobian despite exact invariants; finite-width Bose heating despite positivity; missed pairs of regular roots; fixed-mesh width collapse; treating a critical root as a finite regular measure.

The unresolved next step is a controlled mapping between off-grid resonance nodes and a finite population representation while retaining the correct continuum measure. This branch did not implement such interpolation or run a material mesh campaign.

Artifacts: [script](../experiments/C_resonance_benchmark.py), [results](../experiments/C_resonance_benchmark.json), [pre-run checkpoint](../experiments/C_benchmark_checkpoint.md). The checkpoint allowed an additional odd weight term; before execution it was removed to retain reciprocal symmetry. All results and exact weak-form references above use the final even weight.

Executed:
~~~powershell
py -B 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments\C_resonance_benchmark.py' --output-dir 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments'
~~~

All script assertions passed. Portable requirements: Python >=3.10, NumPy, SciPy. The script accepts --output-dir and otherwise writes beside itself.

Interpreter: C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe. Versions: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1.

Script SHA-256: b3edc1ada2d6e0a8c8368c2ab714b41f1eae698c3d0bd94cf300b4fb5e487949.

