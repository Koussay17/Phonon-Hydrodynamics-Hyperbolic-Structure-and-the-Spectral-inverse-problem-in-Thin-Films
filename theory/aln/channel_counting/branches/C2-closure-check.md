# C2: repeated-daughter closure, entropy and moment tangency

**Completed 2026-09-27; independent exact checks, awaiting campaign audit.** Read the A-D first-pass reports as authorized. Preserved C's first pass. No material run, repository import, or independent re-audit of D's source mapping.

## 1. Assumed kinetic model and agreement of conventions

A/B/C agree after setting \(V=T/6\): the repeated-daughter monomial coefficient is \(g=3V=T/2\). Its Fock-state forward and reverse occupation factors are \(X(Y+1)(Y+2)\) and \((X+1)Y(Y-1)\). Geometric averaging supplies a further factor two, so the population-event coefficient is \(\kappa=2c\), where \(c>0\) is the finite transition multiplier proportional to \(|g|^2\).

For this test, **assume** a classical Markov jump model for \(p\rightleftarrows d+d\):
\[
(X,Y)\longrightarrow(X-1,Y+2),\quad
w_+=cX(Y+1)(Y+2),
\]
\[
(X,Y)\longrightarrow(X+1,Y-2),\quad
w_-=c(X+1)Y(Y-1).
\]
Forbidden jumps have zero rate. This is an additional kinetic model, not the coherent evolution of a finite cubic Hamiltonian or a numerical assignment to \(\delta(0)\).

Exact resonance means \(E_p=2E_d\); set \(E_d=1\) for counting. The conserved quantity is \(2X+Y\). Each fixed-energy sector is finite, and geometric initial distributions have all polynomial moments finite. The instantaneous generator calculations below therefore do not require a time integrator or occupation-boundary approximation.

## 2. Exact mean hierarchy and geometric flux

Write \(a=\langle X\rangle\), \(b=\langle Y\rangle\), and let \(J=\langle w_+-w_-\rangle\). Direct expansion, with no statistical closure, gives
\[
J=c\left[4\langle XY\rangle+2a-\langle Y(Y-1)\rangle\right],
\qquad
(\dot a,\dot b)=(-J,2J).
\]
Both an intermode correlation and a daughter factorial moment appear.

At an independent product-geometric distribution,
\[
\langle XY\rangle=ab,\qquad
\langle Y(Y-1)\rangle=2b^2,
\]
and hence
\[
J=2c[a(1+2b)-b^2]
=\kappa\{a(1+b)^2-(1+a)b^2\},\qquad \kappa=2c.
\]
This identity is exact for the instantaneous flux of that initial distribution. Continuing it as a closed mean equation is a separate assumption.

**Mean-only counterexample (exact finite distributions, \(c=1\)):**

| Initial law | \((a,b)\) | \(\langle Y(Y-1)\rangle\) | Forward rate | Reverse rate | Net \(J\) |
|---|---|---:|---:|---:|---:|
| \(X=1,Y=1\) | (1,1) | 0 | 6 | 0 | 6 |
| \(X=1\), \(Y=0,2\) with probabilities 1/2 each | (1,1) | 1 | 7 | 2 | 5 |
| \(X=1\), geometric \(Y\) of mean 1 | (1,1) | 2 | 8 | 4 | 4 |

The parent is fixed, so no parent-daughter correlation ambiguity is needed to produce the discrepancy.

## 3. Closed Jacobian, entropy PSD and the energy kernel

For the geometric mean equation, the loss matrix is minus its flow Jacobian:
\[
L=\kappa
\begin{pmatrix}
1+2b & -2(b-a)\\
-2(1+2b)&4(b-a)
\end{pmatrix}.
\]
At a positive equilibrium, \(a=b^2/(1+2b)\). Define
\[
r=(1,-2)^T,\quad
W=\operatorname{diag}(a(1+a),\,b(1+b)),\quad
B=a(1+b)^2.
\]
The symbolic calculation verifies exactly
\[
L=\kappa B\,rr^TW^{-1},\qquad
C=W^{-1/2}LW^{1/2}
=\kappa B\,(W^{-1/2}r)(W^{-1/2}r)^T.
\]
Thus the susceptibility-normalized matrix \(C\) is symmetric positive semidefinite, of rank one for \(b,c>0\). Equivalently, the quadratic departure from equilibrium entropy,
\(\mathcal H_2=\tfrac12\delta n^TW^{-1}\delta n\), satisfies
\[
\dot{\mathcal H}_2
=-\kappa B(r^TW^{-1}\delta n)^2\le0.
\]
For \(e=(2,1)^T\), the exact kernel identities are
\[
e^TL=0,\qquad LWe=0,\qquad C W^{1/2}e=0.
\]
The nonzero eigenvalue is
\[
\lambda=2c\,{1+8b+8b^2\over1+2b}>0.
\]

An exact sample at \(a=1/3,b=1,c=1\) gives
\[
L=\begin{pmatrix}6&-8/3\\-12&16/3\end{pmatrix},\quad
C=\begin{pmatrix}6&-4\sqrt2\\-4\sqrt2&16/3\end{pmatrix},
\]
with eigenvalues \(0,34/3\) and symmetric energy-kernel vector \((4/3,\sqrt2)^T\).

**Comparison with D:** D's stated inverse-linewidth comparators are
\(\kappa(1+2b)\) for the parent and \(2\kappa(b-a)\) for the daughter. The independently derived loss diagonals give exact ratios **(1,2)**. The sample above gives event diagonals \((6,16/3)\) and comparators \((6,8/3)\). This confirms D's algebraic ratio within the declared closure; this branch has not independently validated that source comparator or identified an error in a self-energy calculation.

## 4. Product-geometric statistics are not an invariant family

Let the initial law be independent geometric, and differentiate moments using the full jump generator. The exact identities are
\[
{d\over dt}\langle X(X-1)\rangle=-4aJ,
\]
\[
{d\over dt}\langle Y(Y-1)\rangle=2(6b+1)J,\qquad
{d\over dt}\langle XY\rangle=(4a-3b)J.
\]
A geometric daughter marginal would instead require
\[
{d\over dt}(2b^2)=4b\dot b=8bJ.
\]
Consequently the instantaneous departure from the geometric factorial relation is
\[
\boxed{
{d\over dt}\left[\langle Y(Y-1)\rangle-2b^2\right]_{t=0}
=2(2b+1)J.
}
\]
For \(b\ge0\), every initial product-geometric state with \(J\ne0\) fails this tangency condition. The covariance derivative is also
\[
{d\over dt}\left[\langle XY\rangle-ab\right]_{t=0}
=2(a-b)J.
\]
The covariance check alone would miss \(a=b\); the daughter factorial check does not. The parent factorial relation happens to have zero first-derivative defect and cannot certify closure.

At \(a=1/3,b=1/2,c=1\):
- \(J=5/6\);
- actual daughter factorial-moment derivative \(=20/3\);
- geometric tangent prediction \(=10/3\);
- tangency defect \(=10/3\);
- covariance derivative \(=-5/18\).

For \(J=0\), \(a=b^2/(1+2b)\), or equivalently \(q_p=q_d^2\) for geometric ratios \(q_p=a/(1+a)\), \(q_d=b/(1+b)\). A forward edge and its destination's reverse edge have equal rates, and the initial/destination probabilities also agree. These equilibrium product distributions are stationary.

Thus the stationary thermal curve survives, but the two-parameter product-geometric family does not. The tangency defect is already first order for generic small perturbations around equilibrium: the closed two-mean Jacobian is not, by this argument, an exact autonomous description of the full jump-process response. A projection, fast rethermalization, or another closure justification would be needed.

This resolves an open assumption rather than a coefficient disagreement: D explicitly called geometric closure an ansatz. The entropy PSD and conservation checks above are valid properties of that ansatz; they do not establish its dynamical exactness.

## 5. Independent exact-rational check and reproduction

The symbolic route obtains geometric raw moments by differentiating the generating series with \(q\,d/dq\). A separate Python Fraction calculation enumerates initial Fock states directly for \(a=1/3,b=1/2,c=1\). It sums generator changes for five observables and the flux, evaluating all destination states without clipping transitions at the summation boundary.

Partial initial probabilities remain unnormalized. Exact rational expectations converge to the hand-derived targets
\[
J=5/6,\quad \dot a=-5/6,\quad \dot b=5/3,\quad
{d\over dt}\langle X(X-1)\rangle=-10/9,
\]
\[
{d\over dt}\langle Y(Y-1)\rangle=20/3,\qquad
{d\over dt}\langle XY\rangle=-5/36.
\]

| Maximum initial occupations \((X,Y)\) | Omitted probability | Maximum absolute error across the six expectations |
|---|---:|---:|
| (8,12) | \(4.44\,10^{-6}\) | \(3.09\,10^{-3}\) |
| (16,24) | \(5.94\,10^{-11}\) | \(1.39\,10^{-7}\) |
| (32,48) | \(1.36\,10^{-20}\) | \(1.16\,10^{-16}\) |

The errors are exact rational truncation errors displayed in decimal, not floating-point agreement thresholds. All symbolic identities and script assertions passed. No time stepping, fitted coefficients, plotting, or material parameters were used.

Artifacts: [script](../experiments/C2_closure_check.py), [results](../experiments/C2_closure_check.json).

Executed:
~~~powershell
py -B 'D:\ResearchLab\orchestration\campaigns\20260926-105330-aln-channel-counting\experiments\C2_closure_check.py'
~~~
Portable requirements: Python >=3.10 and SymPy; remaining imports are standard library. The optional --output-dir PATH flag overrides the script-directory output location.

Actual interpreter: C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe. Versions: Python 3.14.7, SymPy 1.14.0.

Script SHA-256: 14535a5311af3f07653b42d4d991794e3ee21684556088d0ccb8713acda754fd.

**Limits:** conditional jump model and scalar geometric closure only; no material validation, source-normalization re-audit, linewidth definition change, continuum-limit derivation, or novelty claim. The next substantive question is what physical projection or relaxation mechanism could control this closure error; no such additional investigation was performed here.

