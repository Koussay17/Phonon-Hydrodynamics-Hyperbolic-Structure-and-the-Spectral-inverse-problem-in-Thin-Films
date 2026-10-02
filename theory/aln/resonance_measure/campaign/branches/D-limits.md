# D — independent asymptotics and order of limits

**First pass COMPLETE, 2026-09-28.** Read only the campaign question/assumptions and this branch's checkpoint. No peer results, material inputs or large calculations. The results below are explicitly limited to the stated Gaussian point-sampling and scalar-event models.

## 1. Setting and an exact off-shell obstruction

Fix positive mode energies bounded away from zero, beta > 0, a specified scalar Bose kinetic approximation, and positive quadrature/event coefficients. Momentum has already eliminated one daughter variable modulo the reciprocal lattice. Write Delta = epsilon_p-epsilon_a-epsilon_b and nu = (-1,1,1). For the usual common forward/reverse prefactor k >= 0,

\[
F=k\{n_p(1+n_a)(1+n_b)-(1+n_p)n_an_b\}.
\]

At the Bose occupations of the **original** energies, define R=(1+n_p)n_a n_b > 0. Since n/(1+n)=exp(-beta epsilon),

\[
F=kR(e^{-\beta\Delta}-1),\qquad
\dot E=-\Delta F=kR\Delta(1-e^{-\beta\Delta})\ge0.
\]

Equality holds only if k=0 or Delta=0. Thus every positively weighted off-shell event produces energy drift with the same sign. Positive Gaussian weights on unchanged off-shell tuples cannot cancel this drift by summing more channels. Bose stationarity and exact bare-energy conservation both fail in this class. This is an obstruction for the stated incidence, energies and common-rate flux, not an impossibility theorem for all representations.

A separate positive rank-one surrogate also exposes the invariant defect: if M=sum_e lambda_e nu_e nu_e^T with lambda_e>0, then

\[
\epsilon^TM\epsilon=\sum_e\lambda_e\Delta_e^2.
\]

Positive semidefiniteness therefore does not put energy in the kernel. Such a symmetric surrogate must not be confused with the actual off-equilibrium Jacobian. Enforcing conservation of an interpolated/surrogate energy would establish a different invariant unless its relation to the original energy were controlled.

## 2. Regular resonances: continuum and sampling limits are different

Use the normalized Gaussian delta_sigma(z)=exp(-z²/(2 sigma²))/(sqrt(2 pi) sigma). Sigma has energy units, h momentum units, and v=|grad Delta| energy/momentum units. For a compact smooth regular resonance surface with v bounded away from zero, coarea gives

\[
I_\sigma=\int a(q)\delta_\sigma(\Delta(q))dq
 =\int\delta_\sigma(z)A(z)dz,\quad
A(z)=\int_{\Delta=z}\frac{a}{|\nabla\Delta|}\,dS.
\]

If A is twice continuously differentiable near zero and far-level contributions are controlled, I_sigma=A(0)+(sigma²/2)A''(0)+o(sigma²). Stronger smoothness gives the corresponding higher-order expansion. This describes mollification bias, before discretization. Critical levels, acoustic singularities and unsuitable domain boundaries are excluded.

An exact sampling counterexample is already supplied by Delta(x)=v x on the line. Poisson summation yields

\[
h\sum_{k\in\mathbb Z}\delta_\sigma(vh(k+\theta))
=\frac1{|v|}\left[1+2\sum_{m\ge1}
 e^{-2\pi^2m^2r^2}\cos(2\pi m\theta)\right],\qquad
r=\frac{\sigma}{|v|h}.
\]

Consequences with nonzero weights:

- At fixed h, sigma -> 0 gives zero if the root is off-grid, and divergence proportional to h/sigma if a grid node hits the root.
- Fixed finite r leaves an offset-dependent error as h,sigma -> 0. At r=1/2, normalized values are 1.014383772062229 for theta=0 and 0.985616238638923 for theta=1/2.
- Uniform convergence over offsets for this model requires and is obtained by r -> infinity. At theta=0 the first Fourier harmonic bounds the positive error below; the full series tends uniformly to zero as r grows.

This is a normal-sampling result, not a universal necessity for every multidimensional grid, orientation or integration rule. Arithmetic averaging along a surface, interpolation or different quadratures can alter the requirement. A robust local scale to inspect is sigma/(h |grad Delta|); merely decreasing sigma is not a resolution test. On an empty set with min|Delta|>0, the vanishing Gaussian integral is instead the correct limit.

## 3. Periodic solvable benchmark

For Delta(q)=cos(q) on [-pi,pi), unit weight, the exact target is 2. The exact broadened integral is

\[
I_\sigma=\frac{\sqrt{2\pi}}\sigma
 e^{-1/(4\sigma^2)}I_0(1/(4\sigma^2))
 =2+\sigma^2+O(\sigma^4).
\]

The script evaluates the scaled Bessel function to avoid overflow and compares it with a periodic mesh. At N=1024, h=2 pi/N:

| Width protocol | Root-aligned mesh | Half-step mesh | Exact broadened integral |
|---|---:|---:|---:|
| sigma=sqrt(h) | 2.006222886119902 | 2.006222886119902 | 2.006222886119902 |
| sigma=h/2 | 2.028778815604739 | 1.971240031017115 | 2.000009412587573 |
| sigma=h^(3/2) | 10.18591635788130 | 2.895231697447828e-8 | 2.000000231014883 |

The first row's sampling error is 4.44e-16, but its finite-width bias is still 0.006222886. The second retains aliasing despite an almost unbiased continuum mollifier. The third misses the resonance measure catastrophically. Saved sweeps contain N=32 through 1024. The independent real-space and reciprocal-space Gaussian-comb sums agree within 4.45e-16.

## 4. Equilibrium-drift asymptotics and a false convergence certificate

For the nonlinear flux in section 1, let

A_R(z)=integral over Delta=z of kappa R/|grad Delta|,

where kappa excludes delta_sigma. With a smooth regular level, Gaussian symmetry and fixed beta,

\[
\dot E_\sigma=\int\delta_\sigma(z)z(1-e^{-\beta z})A_R(z)dz
=\beta\sigma^2 A_R(0)+o(\sigma^2).
\]

The absolute integrated flux is generally O(sigma); cancellation can make signed observables smaller. These estimates require a suitable finite-moment mollifier and regular level density; they do not hold universally for every broadening profile or critical root. Beta sigma is the relevant thermal width ratio.

A positive-energy periodic example uses epsilon_p=5, epsilon_a(q)=2-cos(q), epsilon_b(-q)=3 and beta=1, so Delta=cos(q). The daughters carry q and -q; no material dispersion is asserted. Here A_R(0)=0.016512965295032284. At N=8192, sigma=0.2 to 0.0125, the energy drift stays positive. The final observed halving order is 2.0017231; drift/[beta sigma² A_R(0)]=1.000397584 at sigma=0.0125. This is a small asymptotic check, not a full operator convergence test.

**Crucial counterexample:** with sigma=h^(3/2), N=1024, equilibrium energy drift is only 3.15e-34 on the aligned mesh and 2.25e-15 on the half-step mesh. Yet their unit-weight resonance measures are respectively 10.1859 and 2.90e-8 instead of 2. Vanishing drift can arise because only exact on-shell nodes survive, or because all sampled collision weights vanish. It does not certify the correct measure or relaxation rate.

## 5. Critical roots require a separate limit

For the isolated one-dimensional quadratic root Delta(x)=x²,

\[
\int_{\mathbb R}\delta_\sigma(x^2)dx
=C\sigma^{-1/2},\quad
C=\frac{\Gamma(1/4)}{2^{5/4}\sqrt\pi}
=0.8600399873245197.
\]

There is no finite regular coarea value. The relevant width in momentum is sqrt(sigma), so the local resolution ratio is sigma/h². On h=1/256 with sigma=h^(3/2), the two offsets' normalized values sqrt(sigma) I_h,sigma differ from C by +/-2.06e-10. With sigma=h²/4, they remain 0.7984198817 and 0.9678828981; with sigma=h³, they are 6.3830764864 and numerical zero. The latter zero includes exponential underflow in an already severely underresolved calculation. The finite test interval's tails are negligible at the tested widths.

More generally, a nonzero smooth local weight at an isolated radial minimum Delta approximately c|q-q0|^m gives a local contribution scaling as sigma^(d/m-1), with spatial width (sigma/c)^(1/m). This statement is for that homogeneous model: a two-dimensional quadratic saddle instead has a logarithmic level-density singularity. Vanishing coupling weights can change the exponent. Even when a critical regularization can be sampled accurately, a finite limiting rate need not exist. For the one-sided quadratic example, equilibrium energy drift scales as sigma^(3/2), not the regular sigma² law.

## Classification and reproducibility

The regular sqrt(h) sequence resolves sampling but retains a known mollification bias. Fixed-ratio and too-narrow sequences are underresolved in the stated models. The critical test resolves a divergent family after normalization; it does not produce a finite resonance measure. The equilibrium obstruction is algebraic and independent of roundoff.

Saved artifacts: `D-limits-check.py` and `D-limits-check.json`. Run with `D:\ResearchLab\envs\aln-phono3py-4.5.0\Scripts\python.exe`; Python/package versions and the script hash are recorded. No input data files are needed. Numerical kernels use expm1 for Bose flux differences and a scaled Bessel function. No solver iteration or timestep is involved. No universal coupled-scaling theorem, conserving material discretization, or AlN conclusion is claimed.
