# C2: finite entropy weak-form verification

**Complete, 2026-09-28; ready for independent audit.** Read A/B/D/L after the first passes were complete. This standalone check implements the proposed finite equations without importing a peer script or production implementation. First-pass C artifacts are preserved. No time integration, material run, or additional model was added.

## 1. Finite model and independent differentiation

The reusable EntropyWeakForm class accepts positive volume weights \(m_\ell\), a full-column-rank volume evaluation matrix \(\Phi\), positive event weights \(\omega_r\), and parent/daughter evaluation rows \(P,A,B\). Set \(b=P-A-B\). One common reconstruction defines
\[
\xi=\Phi\alpha,\quad n=(e^\xi-1)^{-1},\quad
U=\Phi^T(mn),\quad
M=\Phi^T\operatorname{diag}(mn(1+n))\Phi.
\]
At event nodes the same coefficient vector is evaluated by \(P,A,B\), not by interpolating the volume occupations. All tested volume and reaction entropy variables are positive.

With \(a_r=b_r\alpha\) and \(R_r=(1+n_p)n_a n_b\), use
\[
F_r=R_r\operatorname{expm1}(-a_r),\quad
\Lambda_r=R_r\,{\operatorname{expm1}(-a_r)\over-a_r}>0,
\]
including \(\Lambda_r=R_r\) at zero affinity. Then
\[
G(\alpha)=\dot U=-b^T(\omega F)=K(\alpha)\alpha,\quad
K=b^T\operatorname{diag}(\omega\Lambda)b,\quad
\dot\alpha=-M^{-1}G.
\]
The implementation uses a short local series for expm1(z)/z near zero.

Complex-step differentiation independently checks \(DU=-M\). For the collision derivative, the differentiated function forms the direct Bose products \(n_p(1+n_a)(1+n_b)-(1+n_p)n_an_b\), rather than differentiating the assembled matrix \(K\). At a resonant Bose state \(\alpha_*=\beta e\), it checks
\[
DG(\alpha_*)=K_*,\qquad
D\dot\alpha(\alpha_*)=-M_*^{-1}K_*.
\]
**Away from equilibrium, \(K(\alpha)\) is not generally \(DG(\alpha)\).** Their relative differences were 0.681 in the two-node state and 0.0329 in the periodic state. Treating the nonlinear mobility matrix as its Jacobian would therefore be incorrect.

The volume entropy is evaluated as
\[
S_Q=\sum_\ell m_\ell
\{\xi_\ell n_\ell-\log(1-e^{-\xi_\ell})\},
\]
equivalent to the Bose entropy. Its directional derivative along the coefficient equation is checked independently against
\[
\dot S_Q=\sum_r\omega_r\Lambda_r a_r^2\ge0.
\]
For a represented energy \(E_Q=e^TU\), the corresponding check is
\[
\dot E_Q=-\sum_r\omega_r F_r(b_re).
\]
These are quadrature-defined functionals; this experiment does not identify them with exact continuum volume integrals.

## 2. B's two-node exact reference

Use \(\Phi=I_2,\ m=(1,1),\ P=(1/2,1/2)\), both daughter evaluations \((1,0)\), \(\omega=1\), \(e=(1,3)\), and \(\beta=\log2\). Thus \(b=(-3/2,1/2)\) and \(be=0\). The analytic targets are
\[
M_*=\operatorname{diag}(2,8/49),\quad \Lambda_*=4/3,\quad
K_*=\begin{pmatrix}3&-1\\-1&1/3\end{pmatrix},
\]
with eigenvalues \(0,85/24\) for \(C=M_*^{-1/2}K_*M_*^{-1/2}\).
The implementation reproduces all these targets. The capacity condition number is 12.25.

This is a finite interpolation example with fractional weak redistribution, not a physical reaction between only two mode labels. Event normalization is fixed to the stated unit weight.

| Independent numerical check | Two-node case | Periodic case |
|---|---:|---:|
| Relative \(DU+M\) error at equilibrium | 0 | \(5.42\,10^{-16}\) |
| Relative \(DU+M\) error away from equilibrium | \(3.68\,10^{-17}\) | \(2.82\,10^{-16}\) |
| Relative \(DG-K_*\) error at equilibrium | \(1.64\,10^{-16}\) | \(3.03\,10^{-16}\) |
| Relative coefficient-vector-field Jacobian error | \(1.54\,10^{-16}\) | \(3.10\,10^{-16}\) |
| Equilibrium RHS norm | 0 | 0 |
| Relative entropy directional-derivative identity error | \(8.11\,10^{-16}\) | \(5.68\,10^{-15}\) |

For the two-node off-equilibrium state \(\alpha=(0.8,1.7)\), entropy production was 0.13686908974888448. Event energy drift was zero; its independently differentiated value was \(-1.12\,10^{-16}\).

### Deliberate occupation-interpolation failure

At the exact Bose nodal values, linearly interpolating occupations gives parent \(4/7\) instead of \(1/3\). The measured false equilibrium bracket is \(5/7\), and the weak RHS is \((15/14,-5/14)\), although the energy drift remains zero.

The failure is stronger than a stationary-state mismatch. At
\[
\alpha=\beta(1,3)+10^{-3}(-3/2,1/2),
\]
occupation interpolation gives affinity 0.0025 and flux 0.7160959794. Evolving the same retained volume moments with that incorrect flux gives
\[
\dot S_Q=-0.00179023994854948<0,\qquad \dot E_Q=0.
\]
The negative entropy derivative was checked by differentiating the volume entropy along the resulting coefficient velocity. Therefore conserving this energy does not repair the inconsistent occupation reconstruction.

## 3. Physical energy versus a conserved surrogate

Keep the resonant surrogate vector \(e_h=(1,3)\), but assign represented physical energies \(e=(1,2.8)\) or \((1,3.2)\). The event's physical mismatch is respectively \(-0.1\) or \(+0.1\), while \(be_h=0\).

At the physical Bose coefficients \(\alpha=\beta e\):

| Physical mismatch | Surrogate energy drift | Physical energy drift |
|---:|---:|---:|
| -0.1 | 0 | 0.00980435927589948 |
| +0.1 | 0 | 0.00873397180711749 |

Both agree with \(\beta\sum\omega\Lambda(be)^2>0\). At the surrogate Bose coefficients the finite RHS vanishes, but the volume occupations differ from the physical Bose values by up to 0.0248 or 0.0208. Off-equilibrium physical energy derivatives were also checked directly; they have either sign, while the surrogate remains conserved.

These cases confirm the support distinction in B/D. They do not assume that an interpolated dispersion becomes the physical energy merely because the algebra conserves it.

## 4. A's periodic regular example and residual sensitivity

Use A's three branches \(E_P=5.6,\ E_A(q)=3+\cos q,\ E_B=2\), with branchwise basis \(\{1,\cos q,\sin q\}\). The 9 coefficients contain the physical energy exactly. There are 64 normalized volume nodes per branch and 16 parent nodes, each with two analytic daughter roots \(\pm\arccos0.6\). The second daughter momentum is \(p-k\); periodic basis evaluation enforces the momentum convention.

The positive event weights sum to
\[
1/(\pi\sqrt{1-0.6^2})=0.3978873577297384.
\]
The measured sum differs by less than \(10^{-15}\) relative. Set \(\beta=0.7\) and use a fixed small coefficient perturbation specified in the script. Its smallest sampled \(\xi\) is 1.34765.

The capacity condition number is 42.07. The normalized collision matrix has rank 6: three null eigenvalues lie at roundoff scale, and the remaining eigenvalues are approximately 0.01683 (twice), 0.08268, 0.57287 (twice), and 0.80528. The normalized energy-kernel residual is \(3.04\,10^{-17}\). Extra null directions are retained; no unique thermal equilibrium is inferred.

Off-equilibrium entropy production is \(3.69675765803087\,10^{-5}\). Event energy drift is zero, and its independent directional derivative is \(-4.38\,10^{-19}\). Centered differences of \(U\) provide another check: capacity errors decrease from \(5.34\,10^{-7}\) at step \(10^{-3}\) to \(6.19\,10^{-11}\) at \(10^{-5}\), then increase to \(1.96\,10^{-10}\) at \(10^{-6}\) as subtraction error intervenes.

For a residual test only, move daughter nodes to
\(k=\pm\arccos(0.6-\rho)\), keeping the nominal positive coarea weights fixed. These are intentionally off-shell nodes with \(\Delta=\rho\); this is not an alternative resonance quadrature. Use \(\rho=\pm10^{-2},\ldots,\pm10^{-6}\) and the physical Bose state.

- Energy drift is positive for both signs.
- At \(\rho=\pm10^{-2}\), drift is approximately \(8.10\,10^{-7}\) and \(8.17\,10^{-7}\).
- At \(\rho=\pm10^{-6}\), drift is approximately \(8.13821\,10^{-15}\).
- Dividing by the independently calculated leading coefficient times \(\rho^2\) gives 0.999999589 and 1.000000411 at the smallest residual.
- The moment RHS norm divided by \(|\rho|\) approaches 0.0149176: vector-field drift is first order even though thermal energy drift is second order.
- The eventwise Cauchy-Schwarz energy-drift bound passes for every case.

### Initial failed assertion and its resolution

The first execution stopped on an excessively strict relative comparison of thermal heating with its positive-square formula at small residual. Finite arithmetic evaluates \(b(\beta e)\) and \(\beta(be)\) differently. Rather than silently discard this failure, the final output records their difference and the associated absolute error allowance:
\[
\sum_r\omega_r\Lambda_r|\Delta_r|\,
|a_r-\beta\Delta_r|,
\]
plus a stated floating-point summation allowance. This follows by subtracting the two algebraic heating formulas.

For \(\rho=10^{-6}\), the largest affinity distributivity residual was \(7.11\,10^{-16}\); the heating-formula discrepancy was \(8.26079\,10^{-24}\), covered by the recorded allowance \(8.26090\,10^{-24}\). The rerun passed all checks. This diagnoses arithmetic conditioning; it is not a formal floating-point proof.

## 5. Stability, agreement and limits

The stable flux and mobility evaluation was compared with 80-digit standard-library Decimal arithmetic for affinities \(0,\pm10^{-16},\pm10^{-12},\pm10^{-8},\pm0.1\). Maximum observed relative flux error was \(2.21\,10^{-16}\); the mobility rounded identically to the reference in these cases. Direct subtraction of \(R e^{-a}-R\) had relative errors 0.480 and 1.0 at \(a=\pm10^{-16}\). Stable evaluation cannot recover an affinity already lost when forming \(b\alpha\).

The experiment agrees with A/B's finite identities and D's off-shell heating obstruction. L supplies relevant prior-art context; no new literature attribution or novelty claim is made here. Failed approaches explicitly preserved are occupation interpolation, conflating physical and surrogate invariants, identifying \(K\) with a general nonlinear Jacobian, and an inadequately conditioned residual assertion.

**Scope:** local finite identities and bounded arithmetic only. No time integrator, global admissible-domain preservation, continuum volume-quadrature convergence, spectral/dynamical convergence, or material validation is established. No further calculation is needed for this bounded check; hostile review can use the saved exact references and negative controls.

Artifacts: [script](../experiments/C2_weak_form.py), [results](../experiments/C2_weak_form.json), [design checkpoint](../experiments/C2_weak_form_checkpoint.md).

Executed command:
~~~powershell
py -B 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments\C2_weak_form.py' --output-dir 'D:\ResearchLab\orchestration\campaigns\20260927-164023-conserving-resonance-measure\experiments'
~~~
Portable requirements: Python >=3.10 and NumPy only; the other imports are standard library. Actual interpreter: C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe. Versions: Python 3.14.7, NumPy 2.5.3.

Final executed script SHA-256: 933e7ae55004274b2c33c6fcc709879bc1f49bb1bd7b3a5cf425d69e7c136e5c.


### Preserved first-failure values

The initial failed case was rho=-0.0001: observed heating 8.1378800375857028e-11 versus positive-square value 8.13788003763733e-11. Their difference 5.1627702450672954e-22 exceeded the original allowance 8.13788003763733e-23. The measured affinity multiplication-order residual was 4.4408920985006262e-16; its explicit diagnostic allowance is 5.174342788440831e-22.

The [failure record](../experiments/C2_failure_record.md) preserves the original predicate and traceback, the unchanged numerical values identifying the first failing case, and exactly what changed. The initial traceback did not print case values; those are identified from the unchanged calculation in the final JSON. No additional experiment was performed.

