# D: independent asymptotic and high-precision numerical limits

**First pass complete, 2026-09-23.** Computations were completed on 2026-09-18; final report persistence resumed after a tool-usage interruption. No peer reports, repository prototype, or peer conclusions were read. No repository files were edited. The saved results contain 12 precision cases, 24 directional-difference steps, and 56 binary64 small-energy cases.

## 1. Definitions, assumptions, and scope

Let x_i=E_i/(k_B T)>0, nbar_i=1/expm1(x_i), d_i=nbar_i(1+nbar_i), and R=diag(sqrt(d_i)). Use frozen linear entropy coordinates y=R^-1(n-nbar); the full nonlinear entropy is not being claimed quadratic.

For a declared decay event p -> a+b, accumulate its incidence vector s_p=-1 and +1 for each daughter occurrence. The stipulated mean-population model is

    dn/dt = gamma s F(n),
    F = n_p(1+n_a)(1+n_b) - (1+n_p)n_a n_b,
    gamma >= 0.

The inverse process is already included in F. A spatially reversed event is a separate channel. Gamma=1 in this experiment sets an arbitrary event scale and is not an absolute AlN rate.

Assumptions: equal full-grid weights, positive finite temperature and energies, fixed gamma during differentiation, exact resonance when asserting stationarity, homogeneous populations, and no streaming or boundary terms. This branch does not test momentum conservation without actual wavevectors.

For repeated daughters a=b, the tested polynomial uses the same n_a twice and s_a=2. **This nonlinear closure is an assumption.** Exact same-mode quantum transition factors involve factorial occupation moments, which general mean populations alone do not determine. Thermal geometric statistics give factors proportional to 2 n_a^2 and 2(1+n_a)^2, but a microscopic symmetry factor and a closure away from equilibrium must be specified separately. Differentiating the stipulated polynomial cannot settle that microscopic counting question.

## 2. Independent derivative

The cubic terms cancel exactly:

    F = n_p + n_p n_a + n_p n_b - n_a n_b.

For distinct indices,

    g_p=1+n_a+n_b,  g_a=n_p-n_b,  g_b=n_p-n_a,
    C = -gamma R^-1 s g^T R,

where dy/dt equals a possible base-state drift minus C y at first order. Repeated indices accumulate their gradient contributions.

At exact resonance x_p=x_a+x_b, the common equilibrium flux is

    Q = nbar_p(1+nbar_a)(1+nbar_b)
      = (1+nbar_p)nbar_a nbar_b.

Direct substitution gives g_i=-Q s_i/d_i and therefore

    C = gamma Q (R^-1 s)(R^-1 s)^T.

This is symmetric positive semidefinite. The entropy-coordinate energy vector e_i=x_i sqrt(d_i) obeys e^T C=0 because x^T s=0. These are consequences of the declared finite event model, not a proof of material realizability.

The experiment differentiated the **unexpanded nonlinear Bose products** and compared them with the independent quadratic gradient and the resonant outer product. It did not merely compare two outer-product implementations.

## 3. Differentiation convergence and precision

### A misleading convergence test

The net flux is quadratic, so central differences in the frozen entropy coordinates have **zero algebraic truncation error**. They cannot demonstrate an O(h^2) truncation regime here. For distinct daughters a one-coordinate forward difference is also exact because the polynomial is affine in each coordinate separately.

A genuine first-order test therefore varied several populations together. Directions in relative populations were (1,0.3,-0.4) for x=(1.8,0.7,1.1) and (1,0.3) for repeated x=(1.4,0.7). All perturbed populations were positive. At 90 decimal digits:

| Relative step eta | Distinct forward error | Repeated forward error |
|---|---:|---:|
| 1e-1 | 1.39418e-2 | 1.82089e-2 |
| 1e-3 | 1.39418e-4 | 1.82089e-4 |
| 1e-6 | 1.39418e-7 | 1.82089e-7 |
| 1e-9 | 1.39418e-10 | 1.82089e-10 |
| 1e-12 | 1.39418e-13 | 1.82089e-13 |

Observed forward order is 1 throughout. Central errors instead rise from approximately 2e-90 to 2e-79 as the step shrinks: arithmetic cancellation is amplified by division by the step.

### Precision refinement

Coordinate central differences used h_j=eta*nbar_j/sqrt(d_j), eta=1e-6. Thus the two populations are nbar_j*(1 +/- eta). Relative Frobenius errors against the analytic derivative were:

| Energies | 50 digits | 80 digits | 120 digits |
|---|---:|---:|---:|
| (1.8,0.7,1.1) | 5.46e-46 | 8.38e-76 | 5.53e-116 |
| (3e-30,1e-30,2e-30) | 2.83e-16 | 2.22e-46 | 1.66e-86 |
| (1000,400,600) | 6.48e-46 | 6.73e-76 | 4.01e-117 |
| Repeated (1.4,0.7) | 9.61e-46 | 7.01e-76 | 3.32e-116 |

Maximum entrywise relative errors are also recorded in the JSON. The small-energy product evaluation loses additional digits as expected from cancellation; raising precision resolves it.

## 4. Small energies: a zero residual can accompany a zero, wrong Jacobian

At x=(3 epsilon,epsilon,2 epsilon), binary64 Bose populations were evaluated with expm1. With relative population step eta=1e-4:

| epsilon | Raw-product Jacobian relative error | Quadratic-expression error |
|---|---:|---:|
| 1 | 1.01e-12 | 7.12e-13 |
| 1e-6 | 3.11e-7 | 7.86e-13 |
| 1e-10 | 5.29e-3 | 6.25e-13 |
| 1e-14 | 2.779e1 | 8.64e-13 |
| 1e-18 | 1.0: all entries zero | 9.89e-13 |
| 1e-30 | 1.0: all entries zero | 9.31e-13 |

The raw-product equilibrium flux was numerically zero in these displayed runs. **Zero equilibrium residual does not certify the Jacobian.** At large occupation the '+1' terms can disappear in binary64 products, and cancellation destroys derivative information.

For distinct resonant daughters, the one nonzero event eigenvalue satisfies

    lambda/gamma = 1+2(nbar_a+nbar_b-nbar_p).

For x=(3 epsilon,epsilon,2 epsilon), epsilon*lambda/gamma tends to 7/3. The 120-digit calculation gives relative leading-asymptotic errors 2.14e-10 at epsilon=1e-2 and 2.14e-26 at 1e-6. For repeated x=(2 epsilon,epsilon), the limiting coefficient is 4, with respective errors 4.17e-6 and 4.17e-14.

These divergences hold at fixed illustrative gamma. Actual acoustic matrix elements may vanish with energy; no physical rate divergence is inferred.

## 5. Large energies and zero-energy modes

Population positivity requires a symmetric coordinate perturbation to satisfy

    |h_j| < nbar_j/sqrt(d_j) = exp(-x_j/2).

| Energies | Maximum common coordinate-step bound | Numerical issue |
|---|---:|---|
| (20,8,12) | 4.5400e-5 | Binary64 outer product finite |
| (200,80,120) | 3.7201e-44 | Fixed h=1e-6 gives negative populations |
| (1000,400,600) | 7.1246e-218 | Parent population underflows in binary64; outer product becomes nonfinite |

At x_p=1000, the true parent population is approximately 5.07596e-435. High precision nevertheless gives C_pp approximately 1 and lambda/gamma approximately 1. This is an intermediate-scaling failure, not a divergent collision eigenvalue. The entropy transformation can be numerically singular while the nonzero event eigenvalue remains well behaved.

A derivative test can agree on the polynomial continuation through negative populations. That agreement does not validate its physical nonlinear domain.

At **exact x=0**, nbar is undefined; direct evaluation raises ZeroDivisionError. Setting it to zero or a finite cap changes the model. For x=(1+epsilon,epsilon,1), the tested limit is epsilon*lambda/gamma -> 2 from epsilon=1e-2 to 1e-30. Excluding the exact singular point does not give a uniform bound nearby. The physical coefficient scaling and a defined limiting representation remain necessary.

## 6. Detuning changes the equilibrium problem

Let delta=x_p-x_a-x_b and B=(1+nbar_p)nbar_a*nbar_b. Then exactly

    F(nbar)/B = exp(-delta)-1.

The Bose state for the actual energies is not stationary if delta differs from zero. For x=(1.8+delta,0.7,1.1), differentiation at 100 decimal digits gives:

| delta | ||C-C^T||/||C|| | Minimum eigenvalue of (C+C^T)/2 | Normalized left-energy residual |
|---|---:|---:|---:|
| 1e-1 | 1.1715e-1 | -6.2410e-3 | 2.3936e-2 |
| 1e-3 | 1.2041e-3 | -6.4793e-7 | 2.5015e-4 |
| 1e-6 | 1.2044e-6 | -6.4818e-13 | 2.5026e-7 |

The negative eigenvalues are well above 100-digit arithmetic noise. Their sign has an independent analytic check: for C=u v^T, the symmetric part has possibly nonzero eigenvalues (u.v +/- ||u|| ||v||)/2. Nonparallel u and v make the smaller one negative.

**Do not call this automatically a spectral instability.** The tested rank-one matrix can retain a positive nonzero eigenvalue despite an indefinite symmetric part. Nor does a negative quadratic form about a nonstationary base alone disprove a full nonlinear entropy law. The supported finding is loss of the claimed stationary, conserving, symmetric PSD linearization.

Base-state energy drift satisfies x^T s F=-delta F and is O(delta^2), whereas stationarity and symmetry defects are O(delta). At delta=1e-6, the energy drift divided by B is 9.999995e-13 while the symmetry defect is 1.2044e-6. A tiny base-state energy drift can therefore conceal a much larger Jacobian defect.

This counterexample concerns applying the unchanged nonlinear event rule to detuned energies while retaining exact-resonance claims. It is not a conclusion about every broadening scheme; a particular scheme needs its own consistency and refinement analysis.

## 7. Repeated indices, counting, and reciprocity

For event 0 -> 1+1 at x=(1.4,0.7), the incidence is (-1,2) and the differentiated Jacobian is

    [ 2.9728677272689  -2.7996325694343 ]
    [-2.7996325694343   2.6364921829326 ].

Its daughter diagonal is 4Q/d_a. Replacing the repeated incidence by (-1,1) in the outer product produces a 49.88% relative matrix error and normalized energy residual 0.31013. Thus overwriting rather than accumulating an index can yield a symmetric PSD matrix that still violates energy conservation.

For distinct daughters, the orderings (p,a,b) and (p,b,a) have identical flux and incidence. Counting both with identical gamma doubles the action. Positivity and conservation cannot detect that factor-of-two normalization error. The exporter's counting convention must determine whether both are intended.

For a disjoint six-mode spatial-reversal pair, reversal swaps two identical event triples. The parity commutator is zero; the normalized energy residual is 1.18e-91 at 90 digits. A selected test vector gives dissipation 4.62950705347. Nonnegative dissipation for all vectors follows from the stated outer-product argument, not from this one sample.

This pair has rank 2 and kernel dimension 4. It validates an action on six populations, not a full material collision network with energy as its only invariant.

## 8. Artifacts and reliability classification

Owned artifacts:

- `experiments/D_bose_jacobian.py`
- `experiments/D_bose_jacobian_results.json`
- `branches/D-numerical-limits.md`

Reproduce with:

```powershell
py -X utf8 -B C:\Users\Koussay\ResearchLab\orchestration\campaigns\20260917-232509-aln-physical-events\experiments\D_bose_jacobian.py
```

JSON metadata records Python 3.14.7, NumPy 2.5.3, mpmath 1.3.0, executable, source hash, precision, and step sizes. No material input data are hidden in this experiment.

Retained coding failure: the first run used a negative index on an mpmath eigenvalue matrix; it returned zero and caused division by zero in a diagnostic ratio. An explicit final-row index corrected it. The JSON records the failure, and the completed rerun succeeded. No physical conclusion used the failed ratio.

- **Consistency:** the declared resonant event derivative passes independent nonlinear and analytic checks. Repeated-mode microscopic closure and event normalization remain assumptions.
- **Convergence:** directional forward differences show order 1; the Jacobian checks converge under precision refinement on all tested cases. Central differences have no truncation-error order to demonstrate here.
- **Stability and conditioning:** raw small-energy products are cancellation-sensitive; rare-mode entropy scaling is positivity-sensitive and can underflow. These arithmetic failures are not physical spectral instabilities.
- **Floating-point limited:** displayed binary64 raw-product failures and naive high-energy outer-product formation.
- **Structurally inconsistent with the resonant equilibrium form:** tested detuned events; high precision does not remove the defect.
- **Singular/excluded:** exact-zero energy. Its material resolution is **inconclusive**.
- **Not assessed:** actual AlN coefficients, microscopic combinatorics, momentum conservation with actual wavevectors, event completeness, mesh convergence, broadening convergence, force-constant convergence, and material-rate accuracy.

No time integrator, PDE grid, finite-domain boundary closure, iterative solver, or nonlinear spectral discretization is used. CFL, timestep error, spatial dispersion/dissipation, artificial-boundary effects, and aliasing are not validated by these calculations. Differentiation-step convergence is not timestep convergence.

## 9. First-pass conclusion

The declared exactly resonant positive-event model has a symmetric PSD entropy Jacobian, and independent nonlinear differentiation verifies it in the tested positive-energy regimes when arithmetic precision and perturbation scales resolve the calculation. This does not settle a microscopic repeated-mode closure or absolute counting factor.

The concrete failure cases are: a numerically zero yet wrong Jacobian; physically inadmissible fixed entropy steps; high-energy intermediate underflow; undefined zero modes; lost conservation from repeated-index overwrites; and loss of the resonant equilibrium identities for unchanged detuned event rules. No result here is an AlN calculation, an absolute material rate, or a novelty claim.
