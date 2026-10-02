# C2: finite entropy weak-form check

Status: design checkpoint saved; experiment pending. Authorized second-generation reading completed: A-coarea.md, B-weak-form.md, D-limits.md and L-prior-art.md. First-pass C artifacts preserved.

Own a standalone finite experiment with positive volume weights m and evaluation matrix Phi; reaction evaluations P,A,B give row b=P-A-B. Reconstruct xi=Phi alpha and n=1/expm1(xi), retained U=Phi^T(m n), capacity M=Phi^T diag(m n(1+n)) Phi. At reactions use stable F=R expm1(-b alpha), R=(1+n_p)n_a n_b, Lambda=-F/(b alpha) with its continuous value. Then G=dot U=-b^T(omega F)=K alpha, K=b^T diag(omega Lambda)b, dot alpha=-M^-1 G.

Planned bounded checks:
- Independently differentiate U and direct Bose-product G by complex steps, checking dU/dalpha=-M and dG/dalpha=K at resonant Bose equilibrium.
- Differentiate volume energy and Bose entropy along dot alpha; compare energy/entropy event identities off equilibrium. K is not assumed to equal the nonlinear Jacobian away from equilibrium.
- B's exact two-node example: Phi=I, P=(1/2,1/2), both daughters (1,0), b=(-3/2,1/2), energy e=(1,3), beta=log 2; M=diag(2,8/49), Lambda=4/3, normalized generator eigenvalues 0 and 85/24.
- Deliberately interpolate occupations to reproduce the false equilibrium bracket 5/7; test a nearby state for entropy failure.
- Compare physical and surrogate energy vectors, and prescribed small residuals in A's periodic regular model.
- A's inexpensive branchwise {1,cos,sin} test, with positive volume and coarea weights; no dynamic convergence or uniqueness claim.
- Stable exprel evaluation near zero affinity, checked against 80-digit decimal arithmetic.

Portable script will accept --output-dir. Only experiments/C2_* and this report are edited. No time evolution is required to answer the bounded local-identity question; no material runs or global domain-invariance claim.

