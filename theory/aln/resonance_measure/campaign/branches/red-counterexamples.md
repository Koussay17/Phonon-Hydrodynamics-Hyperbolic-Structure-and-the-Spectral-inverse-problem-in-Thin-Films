# Independent counterexample audit — conserving resonance measure

**COMPLETE — bounded independent review, 28 September 2026.** Read the frozen candidate/draft, synthesis, A/B/C/D, A2/B2/C2 and the relevant C2 source/failure record. No other red reports read; no material run or repository edit.

**Verdict:** the narrow finite/local theorem was not refuted. There is an actual numerical failure of the frozen C2 helper on data satisfying its mathematical domain, and explicit counterexamples to extending the finite structure into global realizability or autonomous reduced kinetics. The latter are limitations already excluded by the candidate, not contradictions of its stated local conclusion.

## 1. Actual frozen-code counterexample: finite answer becomes NaN

**MEDIUM — numerical range gap, not a failure of the real-arithmetic theorem or of an existing saved test case.** Use unit weights,
Phi=I2, e=(1,1), alpha=(1,2), parent row P=(1,0), and daughter rows A=B=(-399.5,400). Then b=(800,-800), b.e=0; physical leg energies are (1,0.5,0.5); entropy evaluations are (1,400.5,400.5). Every required evaluation is positive. The volume capacity condition number is only 5.08616.

Frozen C2 reaction_data returns NaN for both flux and mobility, and its default rhs is all NaN. The reverse Bose product underflows to zero while expm1(-affinity), with affinity=-800, overflows. The independently evaluated direct Bose products give finite RHS (-465.58136549546117,+465.58136549546117). A 100-digit Decimal calculation gives

    F = 0.581976706869326424385002005109...
    Lambda = F/800 = 0.000727470883586658030481252506...

This reproduces overflow/invalid-multiply warnings without an ill-conditioned volume capacity or inconsistent energy. Large extrapolating evaluation rows are allowed by the finite hypotheses. The saved small-affinity tests do not cover this regime; “stable near zero” cannot be promoted to stability throughout Omega.

## 2. Exact finite weak form can leave its admissible domain in finite time

**HIGH against a global interpretation; explicitly permitted by the theorem's local scope.** Set Phi=I3, unit weights, e=(1,1,1), P=(2,-1,0), A=B=(1/2,-1,1). Then b=(1,1,-2), leg energies are (1,1/2,1/2), and resonance is exact. Start at alpha=(1.0005,2,3), giving parent entropy delta=P.alpha=0.001 and daughter entropy 1.50025.

Write h(x)=1/[n(x)(1+n(x))]=exp(x)+exp(-x)-2. The actual equation is alpha_dot=M^-1 b F. Consequently,

    delta_dot = D(alpha) F,
    D(alpha)=2h(alpha1)-h(alpha2).

On the box alpha1 in [1,1.002], alpha2 in [2,2.01], alpha3 in [2.98,3], while 0<delta<=0.001, D<-3.3426 and delta F>1.5728. Hence d(delta^2)/dt<-9 and the positive-domain endpoint occurs by T<=1.11112e-7 in the declared unit-weight time units. A bootstrap closes the box: integrating d alpha_i/d delta gives absolute coordinate changes bounded by (0.000364,0.001866,0.012091). These are smaller than the available margins in the direction of motion.

Thus all volume occupations remain finite and positive, with conserved quadrature energy and increasing finite quadrature entropy, while the off-grid parent occupation diverges as delta approaches zero at a finite time. This is a genuine resonant finite weak-form example. It demonstrates why volume conservation/entropy cannot supply the missing global off-grid realizability conclusion. The proof uses differential inequalities; no time integration was credited.

## 3. Required attacks on A2/B2 alternatives

**Projected entropy mobility.** Let e=(1,1,1), g_i=log(1+1/n_i), v=(sqrt(n1),1,-1-sqrt(n1)), K0=vv^T. This is smooth and PSD in the positive interior and already energy tangent, so projection leaves dot n=v(v.g). Starting at (1e-8,2,1), on the relevant box L=v.g<-0.2792 and |L|<0.6933. Hence d sqrt(n1)/dt=L/2 and n1 reaches zero by t<=0.0008, while the other coordinates stay positive. Energy, entropy increase and thermal stationarity all hold until then. Strict-positive entropy coordinates end in finite time; a boundary continuation is additional information. This does not establish negative occupations, nor refute B2's explicitly local algebra.

**Dual/constrained entropy step.** For only p<->a+b with energies (2,1,1), D=n_a-n_b is an exact invariant. At n=(1,2,1), maximizing Bose entropy minus ||n-n_old||^2/(2 dt) at fixed energy has initial velocity P_e grad S, giving dot D=log(3/4)=-0.287682... rather than zero. Positivity, energy and entropy checks therefore do not certify the desired kinetics or its additional invariants. B2 already acknowledges this missing derivation.

**Probabilistic/conditional-expectation reduction.** In the reversible repeated-decay model with c=1, the sector 2m+n=4 contains (0,4),(1,2),(2,0), with parent-mean drifts 12,-8,-4. The point mass at (1,2) and the equal mixture of the two outer states both have means (m,n)=(1,2), but parent-mean derivatives -8 and +4. Thus even same-sector reversible laws cannot yield an autonomous equation of these retained means. This attacks closure, not the valid full Markov process; the stated memory/lumpability qualification is essential.

**Enlarged spectral energy.** A normalized positive Lorentzian rho_sigma(x) proportional to sigma/[(x-epsilon)^2+sigma^2] on x>0 has rho_sigma(0)>0 for every sigma>0. Since n_beta(x)~1/(beta x), its reconstructed occupation integral diverges logarithmically, and its constant-test capacity integral rho_sigma n_beta(1+n_beta) dx also diverges. Spectral-energy integrals can remain infrared finite. Yet the delta-width limit has the finite occupation n_beta(epsilon): concentration/normalization of spectral weights does not justify this unbounded-observable limit. For epsilon=beta=1, sigma=0.1, the logarithmic coefficient is 0.03254845. This is a conditional failure of a common broadened ansatz, not a claim that an unspecified physical spectral function is Lorentzian. It blocks an asserted finite population/zero-width limit without actual infrared control.

**Joint geometry/symmetry/Hamiltonian directions.** A regular joint fold such as Delta=p-k^2 still has a fixed-parent density proportional to p^(-1/2); its joint integrability does not give a bounded pointwise parent rate. A symmetry quotient requires invariant states/observables; energy or passive basis covariance alone does not close discarded variables. A finite Hamiltonian or positive finite-time sinc-squared kernel does not supply irreversible Markov dynamics. No counterexample to the explicitly qualified chart/symmetry statements was found.

## 4. Measure assumptions and failed attacks

**LOW — integrability wording.** Ambient L1 integrability alone does not supply a trace on a resonance surface: Delta(x,y)=x and f=|x|^(-1/2) on (-1,1)^2 give an ambient-integrable function with divergent delta regularization. Even representatives equal almost everywhere can differ on x=0. Note23's introductory “integrable f” needs to be understood with the declared trace/weighted surface integrability.

Likewise, a smooth transverse geometry with merely integrable surface density need not give A2's O(ell^d) diagonal-tube rate. In d=1, Delta=p and daughter diagonal p=2a, with surface density |2a|^(-1/2), give tube mass proportional to sqrt(ell). Zero diagonal mass survives. A's standing continuous-kernel regular assumptions exclude this singular density, so this is only a warning against an integrability-only reading, not a counterexample to that stronger setting.

No sign, logarithmic-mean regularity, rank/nullity or omitted derivative-of-K defect was found inside the finite hypotheses. K(alpha)e=0 throughout Omega removes that derivative term at beta e. Opposite off-shell detunings cannot cancel the stated positive thermal energy drift. Critical/vacuum limits are excluded rather than solved. The saved missing-root/wrong-Jacobian controls expose their advertised failures; they do not certify a material atlas. No new defect in those saved root measures is alleged.

C2's retained original heating-tolerance failure and its correction were inspected. The measured affinity multiplication-order allowance follows by subtracting the two computed expressions; it is a conditioning diagnostic using shared quantities, not an independent floating-point error proof. No further flaw in that corrected bounded claim was established.

## Reproduction and closure

Executed py -B experiments/red_counter_resonance_probe.py with exit code 0; all counterexample assertions passed. Script and JSON are saved in experiments/ with that basename. It imports the frozen C2 helper without modifying it; inspected SHA-256:
933e7ae55004274b2c33c6fcc709879bc1f49bb1bd7b3a5cf425d69e7c136e5c.

The actual NaN failure was executed, Decimal checked its finite answer, exact Fraction arithmetic checked the probabilistic counterexample, and elementary arithmetic checked the analytic domain-exit bounds. The spectral/trace limits are analytic arguments, not sampled material evidence. The domain-exit derivation has not had an independent proof audit.

All work required for this bounded review is complete. Global/continuum convergence, physical spectral weights, a material collision operator and novelty remain unestablished; none was attempted or credited.

