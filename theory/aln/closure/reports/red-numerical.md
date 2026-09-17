# Independent red-team numerical audit

Final, 2026-09-17. Targets: `D-constructive.md`, `C-infrared.md`, `experiments/pi_bulk_response.py`, and the repository's `scripts/analyse_aln_response.py`. No other red-team reports were read, and no other red team's conclusions were used. The original targets were not modified.

## Verdict

The strongest objections concern the reliability and scope of numerical evidence, rather than a newly found arithmetic contradiction in the exact finite-operator construction:

1. **High severity if numerical validation is extrapolated across 0<t<1:** the eight-mode computation also fails near **t=1**, an endpoint absent from D's reported sweep. The DC value can be accurate to machine precision while the spectral memory has **77.6% error**. This is a floating-point artifact involving an exactly current-invisible slow mode.
2. **High severity for physical interpretation:** neither the reproducible discrete RTA memory nor agreement of bulk-response implementations establishes the continuum AlN memory, a collision spectral gap, or an FDTR response. A single physical q mesh cannot settle those questions.
3. **Moderate severity for validation:** the two bulk scripts are duplicate computations, and D's dense and spectral paths share the same rounded physical matrix. Agreement within either pair is weaker evidence than independent construction and arithmetic.
4. **Low severity but measurable:** subtraction used to report the 1 kHz memory-pole error loses accuracy in the error itself. The response is accurate; the many displayed digits of that tiny error are not.
5. **Moderate reproducibility defect:** C's JSON files and both bulk-response outputs omit executable/package/input-fingerprint metadata. Launcher names alone do not identify an environment.

These findings do not turn an exact, explicitly restricted operator example into a material-specific counterexample. Realizability by an AlN event network remains unestablished, as D acknowledges.

## Reproducible evidence

Owned artifacts:

- `experiments/red_numerical_audit.py`
- `experiments/red_numerical_results.json`
- this report

Run:

```powershell
py -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260916-210742-aln-spectral-closure\experiments\red_numerical_audit.py
```

The completed audit run used Python **3.14.7**, NumPy **2.5.3**, mpmath **1.3.0**, with executable `C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe`. The exact input NPZ SHA-256 was:

`b6712ecbbbfd2388fecce89557c6d6138309c8fe8bce148a5b9449fa404bbb87`

The JSON also records the script hash. It was saved successfully after both the conditioning tests and the independent physical-data calculation. Source functions were inspected; the independent physical-data calculation did not import the repository's heat-capacity or moment implementation.

## 1. D: residuals and accurate DC do not certify memory

### 1.1 Reproduced small-t failure

The lifted operator C+ee^T has condition number asymptotic to 4/t^2 as t tends to zero. A stable linear solver can therefore return a small backward error while the desired moments have large forward error. This is ill-conditioning, not evidence that the underlying positive exact operator is dynamically unstable.

Using the target's physical-matrix assembly independently reproduced:

| t | DC absolute error | Memory relative error, dense | Memory relative error, spectral | Dense normwise backward error |
|---|---:|---:|---:|---:|
| 1e-3 | 1.832e-11 | 5.494e-11 | 4.685e-11 | 5.98e-17 |
| 1e-6 | 1.693e-4 | 5.079e-4 | 3.586e-5 | 1.05e-16 |
| 1e-8 | 0.28164 | 73.45% | 238.04% | 7.69e-17 |
| 1e-9 | 0.48559 | 99.84% | 909.43% | 5.59e-17 |

At t=1e-9 the projected eigensolver even returns a negative least eigenvalue, -2.770e-19, and negative DC and memory. The exact family remains PSD. That apparent instability is an artifact of an unresolved eigenvalue, approximately 5e-19.

The p=0.3 resolvent remains accurate to about 1.7e-16 absolute in these runs. Its accuracy gives no certification of the near-singular DC slope.

Entrywise rational matrix construction followed by one rounding changes the physical matrix by only about 1.1e-16 in spectral norm at t=1e-8, yet changes the spectral memory error from 238% to 390%. Thus the dense and spectral methods must not be treated as independent validation of the common matrix assembly.

### 1.2 Additional failure near t=1

At the binary64 input obtained from `0.9999999999999999`, the F block has an eigenvalue near 1e-16. In exact arithmetic that block is orthogonal to b, so the exact memory remains approximately **2**. The physical-coordinate computations give:

| Method | DC | Memory | Memory relative error |
|---|---:|---:|---:|
| Dense solve, target assembly | 0.9999999999999999 | 2.000635824 | 0.03179% |
| Projected eigensolver, target assembly | 0.9999999999999993 | 3.551993212 | 77.60% |

The dense unlifted relative residual is only **2.48e-16** and its normwise backward error is **4.87e-17**. The exact condition number is of order 2e16. In the JSON, `exact_condition` uses the stated decimal parameter; this is distinguished from the rounded binary64 parameter used in the target assembly near this endpoint.

Mechanism: roundoff can give a nominally invisible eigenmode a spurious overlap of order machine precision. Its contribution to DC is then negligible, while division by the squared tiny eigenvalue makes its contribution to the memory numerator substantial. A current-invisible exact mode can therefore contaminate a computed memory even when DC and the solve residual look satisfactory.

D's reported numerical error bounds concern its sampled t values, whose maximum is 0.75. They are not bounds for every t>=1e-3 or for the entire interval 0<t<1. The exact family formula is not contradicted by this endpoint failure.

### 1.3 Precision is a parameter, not a certificate

Independent fraction arithmetic constructed each physical entry and used rational Gaussian elimination on the lifted eight-by-eight matrix. At t=1/2, 1e-3, 1e-8, and 1e-20, it gave zero exact residual, unit diagonal, energy conservation, DC=1 and memory=1+t^-2.

Independent physical-coordinate mpmath solves gave:

| t | Decimal digits | DC absolute error | Memory relative error |
|---|---:|---:|---:|
| 1e-20 | 50 | 6.286e-12 | 1.886e-11 |
| 1e-20 | 80 | 2.861e-42 | 8.584e-42 |
| 1e-20 | 110 | 7.106e-73 | 2.132e-72 |
| 1e-40 | 50 | Numerically singular | Not available |
| 1e-40 | 80 | 0.03221 | 9.557% |
| 1e-40 | 110 | 2.344e-32 | 7.031e-32 |

At t=1e-40, the inaccurate 80-digit solution still has a relative residual of **7.97e-43**. D's 80-digit evidence is adequate for its tested t=1e-20 case; it is not uniform evidence as t approaches zero. The exact rational checks establish sampled identities, while the general theorem still requires its mathematical proof audit.

An additional exact sensitivity test changes Walsh A11 by delta and A33 by -delta. Physical diagonals, energy conservation, b, and reversal parity are preserved. At t=1e-8, delta=t^2=1e-16 preserves PSD but changes DC from 1 to 2/3. Delta=-0.49t^2 also preserves PSD and changes DC to 25.5. This quantifies why an O(machine-precision) perturbation of O(1) matrix entries can destroy a moment estimate. It is not an equal-DC counterexample and is not presented as one.

## 2. Bulk RTA arithmetic and error metrics

A separate 75-digit computation evaluated oscillator capacities through the sinh formula using exact SI h and k_B, and summed a/(r+i omega) directly. It used the stored binary64 inputs as the input data; it does not certify their physical accuracy.

| Quantity | Basal | c axis |
|---|---:|---:|
| DC conductivity, W m^-1 K^-1 | 309.7751567849200 | 302.4426160153510 |
| Memory, ps | 149.3040059204473 | 231.6551512578447 |
| Static time, ps | 24.0121410380425 | 20.5811712875260 |
| Second-moment fraction below 2 THz | 61.607511% | 74.421495% |

The maximum relative discrepancy between the stored bulk complex responses and this computation, across both directions and all eight supplied frequencies, was **2.363e-16**. Roundoff in the response sums is therefore not the important uncertainty for this dataset.

### Error-of-error cancellation

At 1 kHz the independent memory-pole relative errors are **2.4048774982894e-12** basal and **8.8927779240247e-12** c axis. The reported errors differ by **2.489e-5** and **2.435e-5** relatively. Subtraction of nearly equal conductivities causes this loss of significant digits.

This was independently checked using the identity, with k_j=a_j tau_j and T=sum(k_j tau_j)/sum(k_j),

    K(s) - K0/(1+s T)
      = s^2/(1+s T)^2 * sum_j k_j (tau_j-T)^2/(1+s tau_j).

It avoids subtracting two O(K0) quantities. Its binary64 evaluation agrees with the high-precision error to about machine precision. This audit identity is evidence of the cancellation defect; the target scripts were not changed.

### Total-response error can hide dynamic error

At 100 MHz the memory-pole error relative to the total response is **2.071%** basal and **5.749%** c axis. Relative to the dynamic correction K(i omega)-K(0), the same errors are **23.81%** and **50.23%**. At low frequency the static-pole error in that correction tends to **83.92%** and **91.12%**, despite a small total-response error. Therefore a small total error cannot by itself support a claim of accurate dynamics.

C uses the conjugate harmonic convention from the two bulk scripts. Their imaginary signs are consistent with that choice; error magnitudes are unaffected. Neither convention supplies spatial streaming, interfaces, or an FDTR observation model.

## 3. C: continuum claims remain undecidable from this mesh

The NPZ has 793 irreducible q representatives, 12 branches and multiplicity sum 13824=24^3. These are not multiple physical mesh resolutions.

Independent summation gives an exact-zero-v_z basal DC fraction **0.1296600236013592**, second-moment fraction **0.3338086546214306**, and fixed-mesh diffuse-film floor **40.16545413984746 W m^-1 K^-1**. The floor follows because exactly grazing nodes retain suppression 1. It is a property of this quadrature, not a continuum thin-film prediction. It also does not measure the error at any particular finite thickness: a continuum grazing neighborhood carries real heat.

The large low-frequency second-moment fractions above reproduce C's sensitivity objection. Agreement of DC or of the direct finite sums supplies no convergence order, omitted-cell bound, or evidence that the continuum second moment is finite. Smooth continuations below the first sampled radius can agree with all stored samples while having different acoustic exponents and different second-moment convergence behavior. This is an inference limitation, not a numerical determination of the actual AlN exponent.

The inspected synthetic radial and angular benchmarks test their stated models. They cannot establish that the actual unresolved AlN dispersion and rates satisfy the synthetic assumptions. They were not independently rerun in this bounded audit; no new synthetic-convergence claim is made here.

## 4. Reproducibility and numerical classification

The two bulk JSON files are exactly equal as parsed data. Their source computations are duplicates apart from paths, so this is replication of execution, not independent verification. This audit's separate arithmetic supplies a stronger cross-check of the discrete sums.

D's JSON records Python 3.14.7 / NumPy 2.5.3 / SymPy 1.14.0 / mpmath 1.3.0. C's two JSON files and the bulk JSON omit environment metadata. In the elevated audit shell, both py and python resolved to the same Python 3.14.7 executable. This observation does not establish what an unrecorded earlier invocation used; a launcher name is not version evidence.

- **Consistency:** sampled exact D identities and the stated discrete bulk formula passed independent arithmetic checks. Material-specific collision and boundary consistency remain outside this evidence.
- **Stability:** direct bulk sums show no relevant roundoff instability. D's physical-matrix moments lose forward accuracy near singular endpoints despite tiny backward errors.
- **Convergence:** D's sampled t=1e-20 calculation converges under precision refinement. Actual AlN q-mesh, angular, and physical-model convergence are **inconclusive**.
- **Conditioning:** D is ill-conditioned as t approaches either 0 or 1. Near t=1, exact current decoupling does not automatically survive finite arithmetic.
- **Error:** D failures are **floating-point limited**. The tiny 1 kHz reported reduction error is cancellation-limited in its last digits. Actual AlN uncertainty is not bounded by the machine-level agreement of finite sums.
- **Under-resolution:** the finite-node grazing floor demonstrates a nonuniform fixed-mesh thin-film limit; this audit does not quantify the finite-thickness error on the actual grid.
- **Not applicable to these direct calculations:** time integration, CFL, PDE spatial discretization, iterative-solver tolerance, numerical dispersion and dissipation, artificial outer boundaries, finite-element coercivity, and nonlinear spectral aliasing. No such numerical methods are being validated by these scripts.

Unexecuted checks: no additional physical q meshes, event-resolved collision operator, scattering broadening sequence, actual FDTR boundary problem, or physical input-uncertainty propagation was available or run. No full rerun of C's synthetic quadrature suite was performed. None of these missing checks is replaced by the high-precision finite sums.
