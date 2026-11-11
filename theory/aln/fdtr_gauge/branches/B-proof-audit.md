# B-proof-audit — line-by-line audit of `proofs/B-analytic.md`

Campaign `20261111-fdtr-gauge-memory`, Workstream B. Audited object: `proofs/B-analytic.md` (11 Nov 2026,
19:49) and its scripts `derivations/B-analytic/`. Model: note 7 (Assumption M of the proof; B1–B3 of
`02-assumptions.md`). Independent recomputations: `derivations/B-proof-audit/` (§1). Each step is classified
VALID / VALID UNDER EXTRA ASSUMPTION / INCOMPLETE / UNJUSTIFIED / FALSE / UNKNOWN. Defects D1–D9 are in §3 with
their minimal repairs.

---

## 0. Verdicts

| Statement | Verdict | Defects |
|---|---|---|
| Lemma 0 (principal branches on $\operatorname{Re}p\ge0$) | VALID | — |
| Lemma 1 (impedance) | VALID under Assumption M; one unstated non-vanishing | D1 |
| Lemma 2 (single-valued quadrupole, essential singularity at $-1/\tau_\ell$) | VALID | — |
| Lemma 3 (s-plane, identity principle, second sheet) | VALID; (iv) needs a one-line joint-analyticity remark | D3 |
| **Theorem 1** (global, $R=0$; fibres on $\Delta$, $D_2$) | **PROVED** (Proof A complete; Proof B complete) | — |
| **Theorem 2** (rank table, $\Sigma=\Delta$) | **PROVED** (A and B complete; C proves the part off $\Delta$ only) | D1 (one sentence), D4 (status wording) |
| Corollary 2.1 | PROVED | — |
| Lemma T ($e^h\notin\mathcal K$) | VALID (classical; hypotheses hold at every use) | D8 (presentation) |
| Lemma LW | VALID | — |
| **Theorem 1′** (global, with $R$) | **PROVED** | — |
| **Theorem 2′** (rank with $R$) | **PROVED** | — |
| **Theorem N2** ($N=2$, rank 4 at every $\theta\notin\Delta$) | **PROVED WITH REPAIRS** (minor, both supplied) | D1, D2 |
| **Proposition TW** (certified two-frequency twin) | **PROVED** (certified; independently re-certified) | D5 (presentation) |
| Proposition NR, Remark (iii) | PROVED WITH REPAIRS (D1); "$N=4$ everywhere" rests on Theorem 2′ (proved) | D1 |
| **Theorem G** (gauge, breaking) | **PROVED** | D9 (wording) |
| Corollary 7.1 | PROVED, with $T$ analytic near $i\mathbb R_+$ and $T_{21}Z+T_{22}\ne0$ at the data | — |
| Proposition L, structural part | INCOMPLETE | D7 |
| Counterexample (free in-plane conductivity) | VALID | — |

No step is FALSE. In logical order the **earliest defect is D1**: Lemma 1 uses, without saying so, that
$BL+1\neq0$ and $1-rE\neq0$ on $\{\operatorname{Re}p\ge0\}\setminus\{0\}$. This is true; a three-line proof
is given in §3. Theorem 2 (analyticity into $C([\omega_a,\omega_b])$), Lemma F in N2/NR (connectedness of
$(0,\infty)^N$) and the log-sensitivities of Definition 4 ($Z\neq0$) depend on it. No conclusion changes.

**Cross-branch finding (§4).** The B4 verdict of `branches/B-counterexamples.md` gives the rank-3 set with
resistance as $\{\tau_R=\tau_\ell,\ b=b_s,\ R=0\}$ only. Theorem 2′ gives rank 3 on all of
$\Delta\cap\{R=0\}$, including $b\neq b_s$. Function-level checks in both proof branches confirm this. The
counterexample branch never sampled that stratum, so its summary should be amended.

---

## 1. Independent recomputations (`derivations/B-proof-audit/`)

**(a) `01_riccati_moments.py` (independent route to the moments).** The impedance $W(\xi)=T/\varphi$ at
Liouville depth $\xi=z/\sqrt a$ obeys $W'=b\,p\,W^2-g/b$, $W(\xi_1)=L$, $Z=W(0)$. This follows from
$\partial_\xi\varphi=-b\,p\,T$ and $\partial_\xi T=-(g/b)\varphi$. With $p=s^2$, $W=1/(b_ss)+V$, the $c_k$ come
from the order-by-order polynomial recursion $v_k'=-g_k/b+\delta_{k0}b/b_s^2+(2b/b_s)v_{k-1}+b\sum_{i+j=k-2}v_iv_j$,
$v_k(\xi_1)=R\,\delta_{k0}$, $c_k=v_k(0)$. This route never uses the closed-form quadrupole or `laurent.py`.
Results:
- $c_0,\dots,c_6$ agree exactly with `07_moment_jacobian.out` (symbolic difference 0 for each); $c_7$ is new.
- $G_{0123}\equiv0$ and $G_{0124}=-2\xi_1^4(b-b_s)(b+b_s)(\tau_R-\tau_\ell)/(b^4b_s^3)$, which is the
  document's $2\xi_1^4(b^2-b_s^2)(\tau_\ell-\tau_R)/(b^4b_s^3)$. Also
  $G_{0234}|_{b=b_s}=-2\xi_1^4(\tau_R-\tau_\ell)^2/b_s^5$ (difference 0).
- At $b=b_s$: $\nabla c_0=(0,-2\xi_1/b_s^2,0,0)$ and $\nabla c_1=(0,2\xi_1^2/b_s^2,0,0)$ are parallel. Of the
  seven index sets with $|K|\le9$, only $\{0,2,3,4\}$ survives, as step (3) of N2 claims.
- $t_{0124}$ and $t_{0234}$ agree with the document (difference 0). $\#\{K:|K|=m\}\le(m+1)^3$ holds for $m\le39$.
- On $\Delta$ the moment matrix $[\nabla c_0;\dots;\nabla c_7]$ has rank 3 ($b\ne b_s$) and rank 2 ($D_2$) at
  rational points.
- Resistance: $c_0=R+\xi_1/b-\xi_1b/b_s^2$ and $G^R_{01234}=-2\xi_1^5(\tau_R-\tau_\ell)Q/(3b^3b_s^4)$, both
  identical to the document. The kernel vectors of Theorem 2′, $(0,0,1,1,0)$ on $\Delta$,
  $(bb_s^2/(b^2-b_s^2),0,0,0,1)$ on $\Delta\cap\{R=0\}$ and $(1,0,0,0,0)$ on $D_2\times\{0\}$, annihilate
  $\nabla c_k$ for $k\le7$ symbolically. The moment ranks on the five strata are 4, 3, 3, 5, 5, as in the table.
- Step ($\star$) of Theorem 2′: the equivalence with $b+kg\equiv kb^2(Rq+1/b_s)^2$ is confirmed symbolically.
  Its odd part is $-2Rb^2kq/b_s$.

**(b) `02_twin_N2_poles.py`.**
- Twin. The Jacobian is computed by the quotient rule on $(B\cosh A+S\sinh A)/(B(B\sinh A+S\cosh A))$. This is
  not the log-derivative formula of `09_*.py`, and it was validated against numerical differentiation of the
  matrix product. The Krawczyk test, run on a box of radius $10^{-15}$ centred on the **printed** 20-digit
  $\theta'$, gives $K(X)\subset\operatorname{int}X$ with $|K_{\rm mid}-X_{\rm mid}|\le4.3\cdot10^{-21}$.
  $\det\mathcal J_2(\theta')=[0.003148421043\pm3\cdot10^{-13}]$, and $Z(1.7i)$ differs as stated.
- N2. Derivatives are taken by numerical differentiation of the tanh form. For $\lambda=3$ and $\lambda=1/2$,
  at six points (interior, $b=b_s$, $\tau_\ell=0$, $\tau_R=0$, $b=b_s$ with $\tau_R=0$, $\tau_\ell>\tau_R$),
  the ratio $\det\mathcal J_2/(\sigma^{m_0}t_KG_K)$ at $\sigma=10^{-2},10^{-3},10^{-4}$ tends to 1 with error
  $\propto\sigma$. Example: 0.8928, 0.98875, 0.998869 (interior, $\lambda=3$).
- Poles. Sampling 20 000 random points of $\{\operatorname{Re}p\ge0\}\setminus\{0\}$, including the boundary
  strata and $R>0$, gives $\max|rE|=0.99957<1$ and $\min\operatorname{Re}(BL)/|BL|=1.6\cdot10^{-6}>0$.

**(c) `03_gauge_misc.py`.**
- The two-load minors and their common zeros are reproduced.
- Step 5: the CAS `solve` returns only the generic branch $\tau'=\tau$. The $\delta=0$ branch comes from the
  hand argument of the proof, which is correct, so CAS output alone would not certify this step.
- The $D_2$ derivative formula is confirmed by numerical differentiation (relative error $1.5\cdot10^{-20}$).
- Theorem G: $\Phi\circ g_\mu=\Phi$. The kernel $(e,\lambda,-C,2\Lambda,0)$ and $\operatorname{rank}D\Phi=4$
  (also at $\Lambda=0$) are confirmed. The derivatives of the six breaking constraints along the orbit are
  $e,\lambda,-C,2a,2\Lambda,2a$, all nonzero.
- $\mathcal C'=\mathcal S/2$ is confirmed as a series identity.

**(d) `04_resonance_R0_stratum.py`.** Independent implementation, derivatives by numerical differentiation,
41 frequencies, 50 digits.
- On $\Delta\cap\{R=0,b\ne b_s\}$ (points with $b=2$ and $b=0.4$) the rank is 3: the singular values
  $s_4/s_1,s_5/s_1$ are of order $10^{-52}$, and both predicted kernel vectors are annihilated to
  $\le5\cdot10^{-52}$.
- Controls: rank 4 on $\Delta$ with $R=0.3$, and rank 5 off $\Delta$ with $R=0$.

---

## 2. Audit trail

### 2.1 Lemma 0
- L0.1 $\arg(1+\tau p)\in[0,\phi)$ for $\phi\in(0,\pi/2]$, $\tau\ge0$ (symmetrically for $\phi<0$): VALID.
- L0.2 $\arg g\in(-|\phi|,|\phi|)$, as a difference of two numbers in $[0,|\phi|)$: VALID.
- L0.3 $\arg(pg),\arg(p/g)\in(0,2\phi)\subset(0,\pi)$ for $\phi>0$, an open interval, so the cut is avoided;
  the sum of principal arguments is the principal argument: VALID.
- L0.4 $AB=\xi_1bp$ and $A/B=(\xi_1/b)g$ on the closed half-plane minus 0 (moduli and arguments): VALID.
  Extending them to $D$ by the identity theorem is implicit: VALID.
- L0.5 $\operatorname{Re}A,\operatorname{Re}B,\operatorname{Re}S>0$, since the root arguments lie strictly
  inside $(-\pi/2,\pi/2)$: VALID.
- L0.6 $D$ is open and connected: VALID.

### 2.2 Lemma 1
- L1.1 $\lambda_{\rm eff}=\lambda/g$, $\sigma^2=pg/a$, $(\sigma e)^2=A^2$, $(\lambda_{\rm eff}\sigma)^2=B^2$,
  product $=e\rho c\,p=AB$, hence the same sign (rederived): VALID.
- L1.2 Transfer matrix $M$, $\det M=1$, even under $(A,B)\to(-A,-B)$ (rederived from
  $T=\alpha\cosh\sigma z+\beta\sinh\sigma z$): VALID.
- L1.3 Substrate load $\varphi_1=ST_1$ on $D$ ($\operatorname{Re}\sqrt p>0$); the contact adds $R\varphi_1$: VALID.
- L1.4 Möbius form: VALID.
- L1.5 Reflection form: $r=(BL-1)/(BL+1)$ needs $BL+1\ne0$, and finiteness of $Z$ needs $1-rE\ne0$. Neither is
  stated: **INCOMPLETE (D1)**. Both are true on $\{\operatorname{Re}p\ge0\}\setminus\{0\}$.
- L1.6 In 1D the GK term is eliminated through the energy balance, so no extra flux condition is required:
  VALID within Assumption M (acknowledged in §11 of the proof).

### 2.3 Lemma 2
- L2.1 Entries $\mathcal C(z)$, $(\xi_1/b)g\mathcal S(z)$, $\xi_1bp\,\mathcal S(z)$ and $C^2-XY=1$: VALID.
- L2.2 Holomorphic on $\mathbb C\setminus\{-1/\tau_\ell\}$ (entire if $\tau_\ell=0$ or $\tau_R=\tau_\ell$);
  value at $-1/\tau_R$: VALID.
- L2.3 Essential singularity. $z$ has a simple pole at $-1/\tau_\ell$ exactly when $\tau_\ell>0$ and
  $\tau_R\ne\tau_\ell$, because the numerator $p(1+\tau_Rp)$ is nonzero there iff $\tau_R\ne\tau_\ell$.
  $\mathcal C$ and $\mathcal S$ are transcendental entire, and $z$ maps a punctured neighbourhood biholomorphically
  onto a neighbourhood of $\infty$. VALID.
- Remark: near the singularity $C$ takes every value infinitely often, so the poles of $F$ (zeros of $YL+C$)
  may accumulate at $s=\pm i\tau_\ell^{-1/2}$. They are still discrete in $U_\theta$, and that is all that is
  used later.

### 2.4 Lemma 3 and the "substrate second sheet"
- L3.1 (i) $F$ is meromorphic on $U_\theta$ because $YL+C\to1$; $F-1/(b_ss)=R+(X-YL^2)/(YL+C)$ is holomorphic
  at 0; the value of $c_0$ is correct: VALID.
- L3.2 (ii) $Z=F(\sqrt p)$ on $D$: VALID.
- L3.3 (iii) Identity principle: VALID. The following details were checked.
  - $U_\theta\cap U_{\theta'}$ is $\mathbb C$ minus at most four points, all on the imaginary axis, so it is
    connected.
  - The data points $e^{i\pi/4}\sqrt{\omega_n}$ lie on the open ray $\arg s=\pi/4$, so the limit point
    $s_*\in U$ for every $\omega_*\in[0,\infty)$.
  - $H$ has no pole at a limit of its zeros. At $s_*=0$, $H$ is holomorphic because the principal part
    $1/(b_ss)$ does not depend on $\hat\theta$ ($R$ only enters $c_0$).
  - Poles of $H$ are discrete in $U$, and removing them keeps $U$ connected.
- L3.4 Second sheet: VALID.
  - Is the continuation legitimate? The proof never continues the data. $F$ is one explicit meromorphic
    function on the connected, $s\mapsto-s$ symmetric set $U$, and $F(-s)$ is its value at a point of $U$.
  - The accumulation hypothesis is used exactly once, in the identity theorem at $s_*$, and it suffices.
    Data accumulating only at $\omega=\infty$ would not suffice for this argument.
  - The essential singularities, and the possible pole accumulation next to them, lie outside $U$ and do not
    disconnect it.
  - A physical continuation along a loop around $p=0$ would pass through $\operatorname{Re}p<0$, but it is
    never used logically.
- L3.5 (iv) Identifying $D_cc_k$ with the Taylor coefficients of $D_cF$ needs joint holomorphy of
  $F-1/(b_ss)$ in $(s,\hat\theta)$ on a $\hat\theta$-uniform disc. This is asserted ("real-analytic in
  $\hat\theta$"), not argued: VALID, with the short justification D3 to be added.

### 2.5 Theorem 1, Proof A
- T1.0 (c)⇒(b): on $\Delta$, $g\equiv1$; on $D_2$, $r\equiv0$ and $Z=1/S$. (b)⇒(a) is trivial. VALID.
- T1.1 Step 1, the two linear relations $\mathsf Lv=0$ (rederived; $C,X,Y$ even in $s$): VALID.
- T1.2 Step 2: the minors and common zeros (rederived, `03`(1)); the identity
  $F-1/S=X(S^2-B^2)/(S(Y+CS))$; $X(0)\ne0$; $F\equiv1/S\iff g\equiv b^2/b_s^2\iff\theta\in D_2$. VALID.
- T1.3 Step 3: VALID.
  - $\mathsf L$ has rank 2 off a set discrete in $U$, so its kernel is a line, and $v\ne0$ by $C^2-XY=1$.
  - $v'=\lambda v$ gives $\lambda^2=1$. $\lambda$ is continuous on the connected set $U\setminus\Sigma_0$
    (locally it is a quotient of components), so it is constant, and $\lambda=1$ because $C,C'\to1$.
  - $C,X,Y$ have no poles in $U$.
- T1.4 Step 4: $\mathcal C'(0)=\tfrac12$ gives local injectivity, so $z=z'$; $zB^2$ and $z/B^2$ with
  $g(0)=1$ and positivity give $\xi_1=\xi_1'$, $b=b'$, $g=g'$. VALID.
- T1.5 Step 5: rational identity and the case split $\delta\ne0$ / $\delta=0$: VALID. The proof's hand
  argument covers both branches, whereas the CAS `solve` alone misses $\delta=0$.
- T1.6 Step 6: VALID.
- **Verdict: PROVED.** The fibres on $\Delta\setminus D_2$ (the half-line in $t$) and on $D_2$ ($\xi$ and $t$
  free) are exactly as stated. The matched case $b=b_s$ off $\Delta$ is identifiable.

### 2.6 Theorem 2
- T2.1 "Real-analytic into $C([\omega_a,\omega_b];\mathbb C)$": needs $\omega_a>0$ and no poles of $Z$ on
  $i\mathbb R_+$. **INCOMPLETE (D1)**; true.
- T2.2 Boundary strata: the derivative is that of the analytic extension in $\tau$ (and $R$). For fixed
  $p=i\omega$, $\arg(pg)$ stays inside $(-\pi,\pi)$ for small negative $\tau$. VALID.
- T2.3 Proof A, $\theta\notin D_2$: VALID.
  - Differentiating $\mathsf Lv=0$ at fixed $F,\tilde F$ gives $D_cv=\mu v$, and $\mu=0$ from $C^2-XY=1$.
  - Hence $D_cz\equiv0$ (since $\mathcal S(z)\to1$) and $D_cB^2\equiv0$. Evaluating at $p=0$ gives
    $c_1=c_2=0$; then $c_3=c_4$ and $c_3(\tau_\ell-\tau_R)=0$.
  - On $\Delta\setminus D_2$, $(0,0,1,1)$ is in the kernel because $Z$ is constant along $\Delta$, including
    $t<0$ in the extension.
- T2.4 Proof A, $\theta\in D_2$: VALID.
  - The derivative formula was rederived by hand: $X/(S(Y+CS))=(1-e^{-2A})/(2B^3)$ on $D_2$. It was also
    checked numerically (`03`(3)).
  - $1$ and $p/(1+tp)$ are linearly independent on any set with an accumulation point, so the rank is 2.
- T2.5 Proof C: VALID for rank 4 off $\Delta$ at every point, boundary strata included.
  - The kernel inclusion follows from Lemma 3(iv).
  - The minors $G_{0124}$ and $G_{0234}|_{b=b_s}$ were rederived by the independent Riccati route (§1a).
  - Proof C gives no information on $\Delta$ (D4).
- T2.6 Real and complex ranks coincide (real on $p>0$, then the identity theorem): VALID.
- **Verdict: PROVED** ($\Sigma=\Delta$; rank 4/3/2 with the stated kernels), with D1 added for T2.1.
- Corollary 2.1: the restrictions (intersecting fibres and kernels with $\{c_4=0\}$, $\{c_3=0\}$,
  $\{c_3=c_4=0\}$) are correct, including the corners. VALID.

### 2.7 Lemma T and Lemma LW
- LT.0 $\mathcal K$ is a field closed under $d/dp$. It is the relative algebraic closure of $\mathbb C(p)$ in
  the field of meromorphic functions on the connected $D_+$. If $P(p,f)=0$ is minimal then
  $f'=-P_p/P_f\in\mathbb C(p,f)$ (characteristic 0). Not stated, but standard: VALID.
- LT.1 $A,B,S,L,r\in\mathcal K$: they are single-valued on $D_+$ by Lemma 0, their squares are rational, and
  $BL+1\not\equiv0$. VALID.
- LT.2 Analytic proof: VALID (classical); precise wording in D8.
  - $\operatorname{tr.deg}\mathbb C(p,h,y)=1$, so there is a compact connected $X$.
  - The germs are identified via a primitive element.
  - $y=e^h$ holds on the connected $X'$.
  - At a pole of $h$, $e^h$ grows like $\exp(c|t|^{-m})$, while $y$ grows only polynomially.
- LT.3 Algebraic variant: VALID, and complete once $y'/y=h'$ is transferred to $\mathcal M(X)$ (it holds on an
  open set). In a uniformizer at a pole of $h$ of order $m\ge1$, $dh/dt$ has a pole of order $m+1\ge2$, while
  $(dy/dt)/y$ has at most a simple pole.
- LT.4 **Hypotheses at every use**: VALID.
  - Exponent differences used: $2A,4A,2A',2(A+A')$ and $2(A'-A)$.
    - Each is single-valued and holomorphic on $D_+$ and algebraic.
    - Each is non-constant: positive on $p>0$ and tending to 0 as $p\to0^+$. $2(A'-A)$ is non-constant by
      the case hypothesis.
  - Coefficients used: $B\pm B'$, $(B+B')r$, $(B'-B)rr'$, $D_cB$, $D_cr$, $D_cA$.
    - All lie in $\mathcal K$. Real-parameter derivatives of elements of $\mathcal K$ stay in $\mathcal K$
      here: $D_cB=B(c_2/b-\tfrac12D_c\log g)$, $D_cA=A(c_1/\xi_1+\tfrac12D_c\log g)$ with
      $D_c\log g\in\mathbb C(p)$, and $D_cr=r_BD_cB+r_Lc_5$.
- LW.1 Minimal relation ($m\ge2$), normalization, strictly shorter derived relation, then
  $(a_2e^{f_2-f_1})'=0$, so $e^{f_2-f_1}=\kappa/a_2\in\mathcal K$, which contradicts Lemma T: VALID.

### 2.8 Theorem 1′ (Proof B)
- 1′.1 $Z\equiv Z'$ on $D_+$ (Lemma 3(iii)): VALID.
- 1′.2 $|rE|<1$ on $p>0$, so the denominators are $\not\equiv0$: VALID.
- 1′.3 The four-term relation (rederived by cross-multiplication): VALID.
- 1′.4 If $A-A'$ is non-constant, all six exponent differences are non-constant: VALID.
- 1′.5 LW gives $B=B'$ and $r\equiv r'\equiv0$ (using $B+B'>0$ on $p>0$): VALID.
- 1′.6 If $A-A'$ is constant, the constant is 0, and the three-term relation with LW gives $B=B'$, $r=r'$: VALID.
- 1′.7 From $(A,B,r)$ to $(\xi_1,b,g,R)$: $r$ is a Möbius function of $BL$, so $L=L'$ and $R=R'$. VALID.
- 1′.8 If $r\equiv0$: $(1+Rb_s\sqrt p)^2\in\mathbb C(p)$ forces $R=0$ (since $\sqrt p\notin\mathbb C(p)$),
  then $\theta\in D_2$. VALID.
- 1′.9 The listed fibres do give equal responses: VALID.
- **Verdict: PROVED.**

### 2.9 Theorem 2′ (Proof B)
- 2′.1 $D_cZ$ and multiplication by $B(1-Q)^2$; $D_cQ=E(D_cr-2rD_cA)$: VALID.
- 2′.2 LW with exponents $0,-2A,-4A$: VALID.
- 2′.3 $D_cB\equiv0$ gives $c_2=0$, $c_3=c_4$, $c_3(\tau_\ell-\tau_R)=0$: VALID.
- 2′.4 ($\star$): $\partial_Rr=2B/(BL+1)^2$ (rederived); equivalence with
  $b+kg\equiv kb^2(Rq+1/b_s)^2$, confirmed symbolically in `01`(5). The step needs $B^2L^2\not\equiv1$, which
  holds since $r\not\equiv0$ and $\operatorname{Re}BL>0$. VALID.
- 2′.5 Parity argument, so $R=0$; then $g$ constant, so $\Delta$; then the value of $k$ (no solution if
  $b=b_s$): VALID.
- 2′.6 Every step is an equivalence, so the listed vectors really are in the kernel (sufficiency): VALID.
  Confirmed at the moment level for all $c_k$, $k\le7$, and at the function level (`04`).
- **Verdict: PROVED.** The rank is 5 off $\Delta$ (all $R\ge0$, all $b$, both boundary strata). On $\Delta$
  it is 4 if $R>0$, and 3 if $R=0$, for $b\ne b_s$ and for $b=b_s$.

### 2.10 Theorem N2
- N2.1 Step (1): the expansion along $s=e^{i\pi/4}\sigma$ and real $\nabla c_k$ are VALID. Termwise
  differentiation needs D3.
- N2.2 Step (2): Cauchy–Binet over an infinite index set is used before absolute convergence is established.
  **INCOMPLETE (D2)**: the bounds of step (4) supply it for $\mu\sigma<\rho'$.
- N2.3 $t_{0124}$, $t_{0234}$ (rederived): VALID.
- N2.4 Step (3): only $\{0,1,2,3\}$ has $|K|=6$; only $\{0,1,2,4\}$ has $|K|=7$; the seven sets with
  $|K|\le9$ are as listed; the parallelism at $b=b_s$ holds (rederived). VALID.
- N2.5 Step (4): $|\det\mathcal T_K|\le24(\mu\sigma)^{|K|}$ by the Leibniz expansion; $|G_K|\le24M^4\rho'^{-|K|}$;
  the count is $\le(m+1)^3$; hence the tail is $\le576M^4\sum_{m>m_0}(m+1)^3(\mu\sigma/\rho')^m
  =O(\sigma^{m_0+1})$. VALID.
- N2.6 The Lemma F conclusion on all of $(0,\infty)^2$ requires $\det\mathcal J_2(\theta;\cdot,\cdot)$ to be
  real-analytic on that connected set, i.e. no poles of $Z$ or $\partial_\theta Z$ on $i\mathbb R_+$.
  **INCOMPLETE (D1)**; true. Without D1 the generic statement would follow only on the connected component
  containing small pairs. The "More precisely" sentence does not depend on D1.
- N2.7 The independent numerical cross-check ($\lambda=3$ and $1/2$, six points) is consistent.
- **Verdict: PROVED WITH REPAIRS** (D1, D2; both minor and supplied in §3). The CAS inputs ($G_K$, $t_K$) are
  now confirmed by two derivation routes and two systems.

### 2.11 Proposition TW
- TW.1 The Krawczyk operator encloses the interval Jacobian (formula verified twice), the data is enclosed as a
  ball, and $K(X)\subset\operatorname{int}X$. This gives a unique root in the box for every datum in the data
  ball. VALID.
- TW.2 $\theta'\ne\theta^*$, $\theta'\notin\Delta$, $\theta'\in\operatorname{int}\Theta$: VALID.
- TW.3 Regularity: $K(X)\subset\operatorname{int}X$ implies that every matrix of the interval Jacobian is
  nonsingular (standard Krawczyk theory). Independently, $\det\mathcal J_2(\theta')\approx3.148\cdot10^{-3}$
  is certified. VALID.
- TW.4 Persistence on an open set of $(\theta^*,\omega_1,\omega_2)$ follows from the *implicit* function
  theorem at a regular root, plus openness of $\theta'\neq\theta^*$, $\theta'\notin\Delta$ and
  $\theta'\in\operatorname{int}\Theta$. VALID (wording, D5). The open set contains points with $b\ne b_s$.
- TW.5 The third-frequency difference was reproduced: VALID.
- **Verdict: PROVED.** It is certified twice; the second certificate is independent (different Jacobian
  formula, box centred on the printed digits). Quantifier check: local identifiability is generic (N2), but
  global identifiability fails on a nonempty open set, so it is not generic in the full-measure sense.

### 2.12 Proposition NR and Remark (iii)
- NR.1 $N\ge3$ is necessary ($2N\ge5$): VALID.
- NR.2 Generic $N=3$ from the certified points plus joint Lemma F: VALID, with D1.
  - Each open stratum is connected, e.g. the open interior contains both $\tau_R>\tau_\ell$ and
    $\tau_R<\tau_\ell$, because the minors are analytic across $\Delta$.
  - Joint genericity in $(\hat\theta,\omega)$ gives, by Fubini, generic $\omega$ for a.e. $\hat\theta$.
- NR.3 Remark (iii) and "$N=4$ at every $\hat\theta\notin\Delta$": VALID. The ratio
  $\partial_{\tau_R}Z/\partial_{\tau_\ell}Z=-(1+i\omega\tau_\ell)/(1+i\omega\tau_R)\notin\mathbb R$ was
  rederived and does not depend on $R$; the step uses Theorem 2/2′ for a trivial total kernel.

### 2.13 Theorem G
- G.1 (i): VALID.
- G.2 (ii): $\Phi\circ g_\mu=\Phi$; $\operatorname{rank}D\Phi=4$ everywhere, including $\Lambda=0$; the fibres
  are single orbits (choose $e$). VALID (`03`(4)).
- G.3 (iii): the kernel is $(e,\lambda,-C,2\Lambda,0)$, and the fibre is the orbit by Theorem 1. VALID.
- G.4 (iv): the two- and three-dimensional fibres on $\Phi^{-1}(\Delta\setminus D_2)$ and $\Phi^{-1}(D_2)$:
  VALID.
- G.5 (v): VALID (wording, D9).
  - Each constraint is transverse to the orbits; the derivatives along the generator are
    $e,\lambda,-C,2a,2\Lambda_0,2a\ne0$.
  - Each constraint meets an orbit at most once, and exactly once except that $\{\Lambda=\Lambda_0>0\}$ misses
    the orbits with $\Lambda=0$.
  - So $\Phi|_{\mathcal P_*}$ is a diffeomorphism onto its image, and Theorems 1–2 give injectivity and rank 4
    off $\Phi^{-1}(\Delta)$. The fibre curve on $\mathcal P_v$ was rechecked.
- G.6 (vi): $\partial\tau_\ell/\partial\ell=0$ at $\ell=0$: VALID.
- **Verdict: PROVED.** Of the breaking mechanisms in the brief: known $e$, $\lambda$, $\rho c$, $a$,
  $\Lambda>0$, and velocity with the Debye relation are proved. Lateral spreading is demoted (counterexample
  valid; isotropic case only structural, D7). The interior source is structural only (Remark S).

### 2.14 Other items
- Corollary 7.1: VALID, provided $T$ is analytic near $i\mathbb R_+$ and $T_{21}Z+T_{22}\ne0$ at the data. The
  real $2\times2$ block of multiplication by $w\ne0$ has determinant $|w|^2>0$.
- Proposition L, structural part: **INCOMPLETE (D7)**. The counterexample with free in-plane conductivity is
  VALID (invariance of $\xi_1$, $b_z$, $a_r$ rechecked).
- Remark S: structural and labelled; identifiability is not claimed.
- Remark "real-axis asymptotics": not load-bearing. The step "$Z$ determines $g/b^2$ exactly" is asserted, not
  proved.
- The numerical observations and Conjecture N3 are labelled as such and were not audited as proofs.

---

## 3. Defects and minimal repairs

**D1 (earliest; Lemma 1; used in Theorem 2 T2.1, N2.6, NR.2, Definition 4).** Insert after Lemma 0:

> **LEMMA 0′.** For $\operatorname{Re}p\ge0$, $p\ne0$, $\theta\in\Theta$, $R\ge0$: $\operatorname{Re}(BL)>0$.
> Hence $BL+1\ne0$, $|r|<1$, $|E|<1$, $|rE|<1$, $1\pm rE\ne0$. So $Z$ and all its parameter derivatives are
> holomorphic, and $Z\ne0$, on a neighbourhood of $\{\operatorname{Re}p\ge0\}\setminus\{0\}$.
> *Proof.* By Lemma 0, $\arg B-\arg S=-\tfrac12\arg g\in(-\pi/4,\pi/4)$, so $\operatorname{Re}(B/S)>0$;
> also $\operatorname{Re}(RB)=R\operatorname{Re}B\ge0$. Then $|BL-1|^2<|BL+1|^2\iff\operatorname{Re}(BL)>0$,
> and $|E|=e^{-2\operatorname{Re}A}<1$. Finally
> $BL\sinh A+\cosh A=\tfrac12e^A(BL+1)(1-rE)$. $\square$

The omission is INCOMPLETE and repairable; no conclusion changes. Numerically consistent: over 20 000 points
$\max|rE|=0.99957$ (`02`(C)).

**D2 (N2 step (2)).** Move the bounds of step (4) before step (2) and state: for $0<\sigma<\rho'/\mu$,
$\sum_K|\det\mathcal T_K||G_K|<\infty$. Hence the multilinear expansion of $\det(\mathcal T\mathcal G)$ converges
absolutely and regroups into $\sum_K\det\mathcal T_K\,G_K$; terms with repeated indices cancel by the usual
Cauchy–Binet argument.

**D3 (Lemma 3(iv); N2 step (1)).** Add: $F-1/(b_ss)=R+(X-YL^2)/(YL+C)$ is jointly holomorphic in
$(s,\hat\theta)$ on $\{|s|\le\rho'\}\times N$, with $N$ a complex neighbourhood of $\hat\theta_0$, because
$YL+C\to1$ uniformly. So $c_k=(2\pi i)^{-1}\oint(F-1/(b_ss))s^{-k-1}ds$ may be differentiated under the
integral. This gives $\nabla c_k$ as the Taylor coefficients of $\nabla F$, and the Cauchy bound
$|\nabla c_k|\le M\rho'^{-k}$.

**D4 (status table, §0 and §11).** "Theorem 2: PROVED (three independent proofs)" holds only for rank 4 off
$\Delta$. The values 3 and 2 on $\Delta\setminus D_2$ and $D_2$ have two proofs (A, B). Proof C gives only
lower bounds, and on $\Delta$ only at sampled points.

**D5 (Proposition TW, presentation).**
- (a) The 20 printed digits are followed by "$(\pm10^{-40})$". The radius-$10^{-40}$ box is centred on an
  85-digit midpoint that is not printed. Either print that midpoint, or state "a unique root lies within
  $10^{-15}$ of the printed values" (certified here; the root is in fact within $5\cdot10^{-21}$ of them).
- (b) Replace "inverse function theorem" by "implicit function theorem".
- (c) Cite the Krawczyk regularity statement, or report $\det\mathcal J_2(\theta')\approx3.148\cdot10^{-3}$.

**D6 (§8 "Check", reporting only).** The quoted ratios 0.918, 0.9915, 0.99915 belong to the interior and
$\tau=0$ points. At $b=b_s$ the output (`08`(b)) gives 0.958, 0.9958, 0.99958 ($m_0=9$), and at the AlN-like
point 0.891, 0.9886, 0.99886. The sentence should be corrected; the conclusion (error $\propto\sigma$) stands.

**D7 (Proposition L).** "The gauge is not a symmetry of $H$" needs $H$ to depend non-trivially on $a$. Repair:
compute the $O(k^2)$ term of $Z(i\omega,k)$ and show its $a$-derivative is $\not\equiv0$. Alternatively,
certify $H(g_\mu P)\ne H(P)$ at one $(P,\mu,\omega)$ in ball arithmetic. Until then the status should read
"structural, dependence on $a$ asserted".

**D8 (Lemma T, presentation).** State the germ identification precisely.
- Take a primitive element $w$ of $\mathbb C(p,h,y)/\mathbb C(p)$ with irreducible minimal polynomial $P$, and
  let $X$ be the compact connected Riemann surface of $P$ (Forster §8).
- Choose $p_0\in D_+$ not over a branch point. The germ of $w$ at $p_0$ is the germ of some $W\in\mathcal M(X)$,
  and $h,y$ are rational in $(\pi,W)$.
- Lead with the algebraic variant.
- The statement is the one-function special case of Ax (1971), already cited in §11.

**D9 (Theorem G(v), wording).** "Injective" should read: injective on $\mathcal P_*\setminus\Phi^{-1}(\Delta)$,
and no point of $\mathcal P_*\cap\Phi^{-1}(\Delta)$ shares the response of a point off it. Both follow from
Theorem 1.

---

## 4. Disagreement with `branches/B-counterexamples.md` (for `06-disagreements.md`)

**The disagreement.** The counterexample branch (§0, row B4, and "Verdict B4" in §5) states the resistance
degenerate set as "rank 4 on $\tau_R=\tau_\ell$; rank 3 on $\tau_R=\tau_\ell$, $b=b_s$, $R=0$". This misses the
stratum $\Delta\cap\{R=0,\ b\ne b_s\}$, where the rank is **3**, with kernel
$\mathrm{span}\{(0,0,1,1,0),\ (bb_s^2/(b^2-b_s^2),0,0,0,1)\}$.

**Evidence that the rank there is 3.**
- Theorem 2′ (Proof B, audited above).
- The symbolic moment identities of §1(a).
- Function-level SVD at 50 digits at $b=2$ and $b=0.4$: deficient singular values of order $10^{-52}$ (`04`).
- `derivations/B-analytic/02` at 60 digits.
- The thin-slab mechanism: since $p(1/S)^2=1/b_s^2$, a slab $d\xi$ at the interface shifts the load by the
  constant $d\xi(1/b-b/b_s^2)$, which mimics a resistance.

**Why the branch missed it.** The stratum list of `b4_resistance.py` (lines 143–151) contains no resonance
point with $R=0$ and $b\ne b_s$. It samples the resonance with $b=2b_s$ only at $R=0.3$, and with $b=b_s$ at
$R=0.3$ and $R=0$.

**Resolution (by proof and computation, not by count).** Theorem 2′ stands, and Verdict B4 should be amended.

---

## 5. Formal-verification candidates

These are short, central steps, suitable for Lean 4 + Mathlib (`field_simp`, `ring`, `linear_combination`):
1. From the stated $c_0,\dots,c_4$: $G_{0123}\equiv0$, the formula for $G_{0124}$, and $G_{0234}|_{b=b_s}$.
   These are the only CAS inputs of Theorem N2 and Proof C. Formalising that the $c_k$ *are* the Laurent
   coefficients is harder; the polynomial Riccati recursion of `01` is the more formalisable route.
2. Step 5 of Theorem 1 (the rational identity with its two branches).
3. The $2\times3$ kernel lemma of Step 3: a rank-2 matrix, two kernel vectors, and the normalization
   $C^2-XY=1$ give equality up to sign.
4. Lemma LW in a differential-field setting, with Lemma T as a hypothesis or via its algebraic variant.

---

## 6. Overall verdict

The analytic branch is sound. **Theorems 1, 2, 1′, 2′ and G are PROVED**, and **Proposition TW is PROVED**
by certified computation, re-certified independently. **Theorem N2 is PROVED WITH REPAIRS**: two minor gaps,
the no-pole Lemma 0′ (D1) and the absolute convergence behind the infinite Cauchy–Binet expansion (D2), both
repaired above.

The specific points of concern check out:
- Branch choices: Lemma 0 holds on the closed half-plane, including the boundary strata $\tau=0$.
- Single-valuedness: Lemma 2 holds, with an essential singularity only at $-1/\tau_\ell$.
- Second sheet: legitimate. It is an identity of one explicit meromorphic function on a connected domain,
  and the accumulation hypothesis is used correctly, at $s_*=e^{i\pi/4}\sqrt{\omega_*}$, including
  $\omega_*=0$.
- Lemma T: its hypotheses (single-valued, algebraic, non-constant exponents; coefficients in $\mathcal K$) hold
  at every use in Theorems 1′ and 2′.
- CAS identities: every CAS identity on which a proof rests was rederived. The Laurent coefficients and minors
  were obtained by a different route (Riccati recursion), and the twin Jacobian by a different formula.

Everything proved is conditional on Assumption M (note 7's model: local temperature convention, no extra GK
flux condition, 1D) and on exact, noise-free data. The second-sheet argument is a uniqueness proof, not an
inversion method. The remaining open items are the proof's own: Conjecture N3, the identifiability part of
Proposition L, Remark S, and D7. The counterexample branch's B4 verdict needs the correction of §4.
