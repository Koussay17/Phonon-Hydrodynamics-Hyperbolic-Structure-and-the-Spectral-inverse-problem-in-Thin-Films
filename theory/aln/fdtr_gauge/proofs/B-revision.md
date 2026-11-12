# B-revision — Interface-temperature conventions (objection B-2) and the transducer configuration (objection B-4)

Campaign `20261111-fdtr-gauge-memory`, revision branch R-B, 12 November 2026. Answers objections B-2 and B-4 of
`12-peer-review.md`. Builds on `proofs/B-analytic.md` (with its addendum and the audit `branches/B-proof-audit.md`);
lemma numbers of that file are quoted as "B-analytic, Lemma 3" etc. Scripts and raw outputs:
`derivations/B-revision/` (§9). Units of all nondimensional examples: $b_s=1$ and the time unit of the frequencies.

Notation. The parameter vector is $\theta=(\xi_1,b,\tau_R,\tau_\ell)$ as in B-analytic; to avoid a clash, the interface
temperature of the referee and of the twin ("theta") is written $\vartheta$ here.

---

## 0. Summary

### 0.1 Conventions and what the twin does

On each face of the GK film the jump condition reads $\vartheta_{\rm above}-\vartheta_{\rm below}=R\,q$ with
$\vartheta=T-c_n\,\partial_zq$ on the GK side ($\vartheta=T$ on Fourier sides). At $k=0$ energy conservation gives
$\partial_zq=-C\,p\,T$, so $\vartheta=(1+\tau_np)\,T$ with the **interface time** $\tau_n=c_nC$. The dimensionless form natural
to $\vartheta$ is
$$
c_n=\kappa\,\frac{\ell^2}{\lambda}\quad\Longleftrightarrow\quad \tau_n=\varepsilon\,\tau_\ell,\qquad
\varepsilon=\frac{\kappa}{1+\alpha}=\frac{c_n\lambda}{(1+\alpha)\ell^2},\qquad \tau_\ell=\frac{(1+\alpha)\ell^2C}{\lambda}.
$$
local: $\varepsilon=0$; entropic (twin): $\varepsilon=1$; kinetic (Beardo et al. 2020, eq. A6, $c_n=\gamma^{-1}(\beta-\chi_{nn})$):
$\varepsilon\in\mathbb R$, between $0.03$ and $0.23$ for their Si coefficients (81–416 K; §1.3); free: $c_n$ unknown, $\tau_n$ an extra
parameter. **The twin reads $T(0)$ at the heated face** (Proposition O, proved from the code and reproduced to $1.5\cdot10^{-15}$),
and applies the layer's convention at **both** faces of a GK layer.

### 0.2 Classification, bare film on a Fourier half-space, read-out $T(0)$ (Theorems C1, C2)

$\Delta=\{\tau_R=\tau_\ell\}$, $D_2=\Delta\cap\{b=b_s\}$, $\Pi=\{\tau_R=\tau_\ell=0\}$ (Fourier corner).
Kernels in the coordinates $(\xi_1,b,\tau_R,\tau_\ell[,\tau_n][,R])$; $k_s:=bb_s^2/(b^2-b_s^2)$ (thin-slab ratio); "$R$ known" means
$R$ is a known value ($=0$ or $>0$), "$R$ unknown" adds it to the parameters.

| convention | stratum | $R=0$ known | $R>0$ known | $R$ unknown, at $R>0$ | $R$ unknown, at $R=0$ |
|---|---|---|---|---|---|
| local $\varepsilon=0$ (4/5 par.) | off $\Delta$ | 4 | 4 | 5 | 5 |
| | $\Delta\setminus D_2$ | 3: $(0,0,1,1)$ | 3: $(0,0,1,1)$ | 4: $(0,0,1,1,0)$ | 3: $(0,0,1,1,0)$, $(k_s,0,0,0,1)$ |
| | $D_2$ | 2: $(0,0,1,1)$, $e_\xi$ | 3: $(0,0,1,1)$ | 4: $(0,0,1,1,0)$ | 3: $(0,0,1,1,0)$, $e_\xi$ |
| tied $\varepsilon\neq0$ (entropic, kinetic; 4/5 par.) | off $\Pi$ | 4 | 4 | 5 | 5 |
| | $\Pi$, $b\ne b_s$ | 4 | 4 | 5 | 4: $(k_s,0,0,0,1)$ |
| | $\Pi$, $b=b_s$ | 3: $e_\xi$ | 4 | 5 | 4: $e_\xi$ |
| free $c_n$ (5/6 par.) | off $\Delta$ | 5 | 5 | 6 | 6 |
| | $\Delta$, $\tau_n\ne0$ | 4: $(0,0,1,1,0)$ | 4: $(0,0,1,1,0)$ | 5: $(0,0,1,1,0,0)$ | 5: $(0,0,1,1,0,0)$ |
| | $\Delta$, $\tau_n=0$, $b\ne b_s$ | 4: $(0,0,1,1,0)$ | 4: $(0,0,1,1,0)$ | 5: $(0,0,1,1,0,0)$ | 4: $(0,0,1,1,0,0)$, $(k_s,0,0,0,0,1)$ |
| | $\Delta$, $\tau_n=0$, $b=b_s$ | 3: $(0,0,1,1,0)$, $e_\xi$ | 4: $(0,0,1,1,0)$ | 5: $(0,0,1,1,0,0)$ | 4: $(0,0,1,1,0,0)$, $e_\xi$ |

Global fibres (Theorem C1): with any fixed $\varepsilon\neq0$ the map is **injective** except on
$F_0=\{b=b_s,\ \tau_R=\tau_\ell=0,\ R=0\}$, where $Z=1/S$ and $\xi_1$ is free. With $c_n$ free, $\tau_n$ and $R$ are always
identified and the resonance line $\{(\xi_1,b,t,t,\tau_n,R):t\ge0\}$ is again a fibre.

**Answers to B-2.** (1) $\Delta$ is degenerate only for $\kappa=0$: for every fixed $\varepsilon\neq0$ the map is rank-full and
injective on $\Delta\setminus\Pi$ (and on all of $\Delta$ if $R>0$). (2) There is no new degenerate set for $\kappa\neq0$, for any value of $\varepsilon$ (exact real decision, §3.3): the
degenerate set shrinks to the Fourier corner with perfect contact, $\Pi\cap\{R=0\}$ (and $b=b_s$ if $R$ is known), where all
conventions coincide because $\tau_n=\varepsilon\tau_\ell=0$. (3) A free
$c_n$ is identifiable (globally and locally) and is never confounded with an interface resistance ($R$ and $\tau_n$ are the
$s^0$ and $s^1$ coefficients of the load $L/(1+\tau_np)$); its price is the return of the resonance degeneracy, because
$\tau_n$ then no longer carries $\tau_\ell$. (4) The resonance direction $(0,0,1,1)$ is in a kernel only if the coefficients at **both**
film faces (or interface and read-out) vanish; with any coefficient tied to $\tau_\ell$ at either face it disappears (§4). The
thin-slab confound on $\Delta\cap\{R=0\}$ ($R$ unknown) persists whenever the film/substrate coefficient vanishes, whatever the
top face. Mixed conventions (entropic read-out or top face, local or negative-kinetic bottom face, unknown $R$) produce isolated
exotic strata of codimension 3 or 4 (Proposition C4).

### 0.3 Transducer configuration (B-4): known Fourier transducer, unknown $R_1,R_2$, film $\theta$ — 6 unknowns (Theorems T1–T4)

| convention (same $c_n$ at both film faces, as in the twin) | stratum | rank | kernel in $(\xi_1,b,\tau_R,\tau_\ell,R_1,R_2)$ |
|---|---|---|---|
| local | off $\Delta$ | 6 | 0 |
| | $\Delta$, $R_2>0$ | 5 | $(0,0,1,1,0,0)$ |
| | $\Delta$, $R_2=0$, $b\neq b_s$ | 4 | $(0,0,1,1,0,0)$, $(k_s,0,0,0,0,1)$ |
| | $\Delta$, $R_2=0$, $b=b_s$ | 4 | $(0,0,1,1,0,0)$, $(1,0,0,0,0,0)$ |
| tied $\varepsilon\ne0$ (entropic $\varepsilon=1$) | off $\Pi\cap\{R_2=0\}$ | 6 | 0 |
| | $\Pi\cap\{R_2=0\}$, $b\ne b_s$ | 5 | $(k_s,0,0,0,0,1)$ |
| | $\Pi\cap\{R_2=0\}$, $b=b_s$ | 4 | $(1,0,0,0,0,0)$, $(0,0,1-2\varepsilon,1,0,0)$ |
| free $c_n$ (7 unknowns: $\tau_n$ after $\tau_\ell$) | off $\Delta$ / $\Delta$ / $\Delta\cap\{\tau_n=0,R_2=0\}$ | 7 / 6 / 5 | 0 / $(0,0,1,1,0,0,0)$ / $+(k_s,0,0,0,0,0,1)$ or $+e_\xi$ ($b=b_s$) |
| | $\Pi\cap\{\tau_n=0,R_2=0,b=b_s\}$ | 4 | $e_\xi$, $(0,0,1,1,0,0,0)$, $(0,0,-2,0,1,0,0)$ |

$R_1$ is never in a kernel. Global (Theorem T2): local convention — injective off $\Delta$, fibre $t\mapsto(\xi_1,b,t,t,R_1,R_2)$
on $\Delta$ ($\xi_1$ also free on $D_2\cap\{R_2=0\}$); tied $\varepsilon\neq0$ — injective except on $\{b=b_s\}\cap\Pi\cap\{R_2=0\}$
($\xi_1$ free). Frequencies (Theorem T4): $N\ge3$ is necessary; three generic frequencies give local rank 6 at generic points
of every stratum checked (21 certified stack points, both conventions, including the Au/AlN/sapphire values of A); five
frequencies suffice at **every** full-rank point; three suffice at every point off an explicit hypersurface (local:
$\Delta\cup\{P=0\}$; entropic: $\{Q=0\}$, which does not contain $\Delta$), by the leading low-frequency term of $\det\mathcal J_3$, and
at every sampled full-rank point (Conjecture T4′ for the hypersurface). Three frequencies do **not** give global uniqueness:
three certified twins under the local convention (Proposition TW3); none was found under the entropic convention, and none
with four frequencies (Conjecture T4″).

### 0.4 Statement status

| Statement | Status |
|---|---|
| Lemma R ($k=0$ reduction), Proposition O (twin observable and face rows) | PROVED (elementary; code inspection) and NUMERICALLY CONFIRMED to $1.5\cdot10^{-15}$ |
| Lemma 3′ (meromorphy, identity principle) | PROVED (as B-analytic, Lemma 3) |
| Lemma K (kernel reduction to three polynomial identities) | PROVED, elementary (exponential versus algebraic decay along a ray) |
| Lemma G (global reduction) | PROVED via the functional Lindemann–Weierstrass lemma (B-analytic, Lemma LW, audited) |
| Theorem C1 (fibres, read-out $T(0)$, all conventions) | PROVED |
| Theorem C2 (rank table, read-out $T(0)$) | PROVED (hand case analysis) and CAS-DECIDED (independent real quantifier elimination); 50-digit check on 62 points |
| Proposition C4 ($\vartheta_e$ read-out), C5 (film block with $\vartheta$ at both faces) | PROVED (case analysis) and CAS-DECIDED; exotic strata verified exactly and at 50 digits |
| Theorems T1, T2, T3 (transducer, local and every tied $\varepsilon$, free $c_n$) | PROVED; 35 stack points at 50 digits (+7 with different face coefficients, +3 free-$c_n$ corner points); 21 stack points certified in ball arithmetic (§5) |
| Theorem T4: necessity $N\ge3$; generic $N=3$; $N=5$ everywhere; $N=3$ off an explicit hypersurface | PROVED (generic $N=3$: certified points + real analyticity; (iv): CAS minors in two systems) |
| Three frequencies: local but not global identifiability (stack, local convention; Proposition TW3) | PROVED by certified computation (Krawczyk, three twins) |
| Four frequencies: global uniqueness at generic points (Conjecture T4″) | CONJECTURED (no twin found in 6 quadruple searches; dimension count) |
| Low-frequency coefficients; leading time combination in the stack (§6) | DERIVED (CAS, two systems for $c_{-1}..c_2$) |

**Weakest step.** All statements are conditional on the 1D ($k=0$) model of Lemma R with a constant $c_n$ (Assumption M′,
§1.1), on exact data, and — for the global statements only — on Lemma LW, whose proof rests on the classical fact that the
exponential of a non-constant algebraic function is not algebraic (B-analytic, Lemma T, audited). The local rank statements do
not need it (Lemma K has an elementary proof). Among the computations, the decision of the degenerate sets for the general
$\varepsilon$ models was done by real quantifier elimination in pieces (the full systems exceeded the time limit of the symbolic
engine); the pieces and the hand case analysis agree, but the decomposition itself is part of the proof (§3.3, §4).

---

## 1. Setting

### 1.1 Model

**ASSUMPTION M′ (replaces Assumption M of B-analytic).** As Assumption M (1D, constant coefficients, zero initial state, no
interior source, unit flux absorbed at the heated face, semi-infinite Fourier substrate of known effusivity $b_s$), except that
at each face of the GK film the temperature entering the jump condition is $\vartheta=T-c_n\partial_zq$ on the GK side, with a
constant $c_n\in\mathbb R$ (the same at both faces unless stated otherwise), and the jump is $\vartheta_{\rm above}-
\vartheta_{\rm below}=R\,q$, $R\ge0$. The read-out at the heated face is $T(0)$ or $\vartheta_o(0)=(1+\tau_op)T(0)$ (§1.2).
Film law and notation as B-analytic §1: $\rho c\,pT+\partial_zq=0$, $(1+\tau_Rp)q=-\lambda\partial_zT+\Lambda\partial_z^2q$,
$\Lambda=(1+\alpha)\ell^2$, $g=(1+\tau_Rp)/(1+\tau_\ell p)$, $A=\xi_1\sqrt{pg}$, $B=b\sqrt{p/g}$, $S=b_s\sqrt p$,
$M=\begin{pmatrix}\cosh A&\sinh A/B\\B\sinh A&\cosh A\end{pmatrix}$, $\det M=1$.

**DEFINITION 1′.** $N(p)=1+\tau_np$ (bottom face), $N_0(p)=1+\tau_op$ (top face or read-out), $L=R+1/S$,
$\mathrm{Mob}(M;x)=(M_{11}x+M_{12})/(M_{21}x+M_{22})$. Tied conventions: $\tau_n=\varepsilon\tau_\ell$, $\tau_o=\varepsilon_0\tau_\ell$ with
known constants; free conventions: $\tau_n$ (and/or $\tau_o$) unknown real parameters. The **film block**
$$
\Phi:=N_0\,\mathrm{Mob}(M;L/N)=\frac1{\tilde B}\,\frac{1+rE}{1-rE},\qquad \tilde B=\frac B{N_0},\quad u=\frac{BL}{N},\quad r=\frac{u-1}{u+1},\quad E=e^{-2A}.
$$
(The second form follows from B-analytic, Lemma 1, with load $L/N$, multiplied by $N_0$.)

**LEMMA R ($k=0$ reduction).** Under Assumption M′, $\vartheta=(1+\tau_np)T$ at a face where the coefficient is $c_n$, with
$\tau_n=c_nC$. The film relation in the variables $(\vartheta,q)$ is $(\vartheta_0,q_0)^{\sf T}=\mathrm{diag}(N_t,1)\,M\,
\mathrm{diag}(N_b^{-1},1)(\vartheta_1,q_1)^{\sf T}$, $N_t,N_b$ the factors of the top and bottom faces. Hence
(a) bare film, read-out $T(0)$: $Z=\mathrm{Mob}(M;L/N)$ ($N=N_b$), explicitly
$$
Z=\frac{L\cosh A+N\sinh A/B}{BL\sinh A+N\cosh A}=\frac1B\,\frac{1+r_ne^{-2A}}{1-r_ne^{-2A}},\qquad
r_n=\frac{BL-(1+\tau_np)}{BL+(1+\tau_np)},\quad L=R+\frac1{b_s\sqrt p};
$$
(b) bare film, read-out $(1+\tau_op)T(0)$: $Z=N_0\,\mathrm{Mob}(M;L/N)=\Phi$;
(c) known Fourier two-port $T_A$ ($\det T_A=1$) above a contact $R_1$: $Z_{\rm tot}=\mathrm{Mob}(T_A;\,W)$,
$W=R_1+\Phi$ with $N_0=N_t$, $L=L_2=R_2+1/S$; with the same coefficient at both faces ($N_t=N_b=N$, the twin)
$W=R_1+\mathrm{Mob}(NMN^{-1};L_2)=R_1+B_\vartheta^{-1}(1+r_ne^{-2A})/(1-r_ne^{-2A})$, $B_\vartheta=B/N$: the film is a
$\vartheta$-quadrupole with the same $A$ and with $B_\vartheta^2=b^2p/(gN^2)$; entropic: $B_\vartheta^2=b^2p/((1+\tau_Rp)(1+\tau_\ell p))$,
symmetric in $(\tau_R,\tau_\ell)$, while $A^2=\xi_1^2p(1+\tau_Rp)/(1+\tau_\ell p)$ is not.
*Proof.* $\partial_zq=-CpT$ gives $\vartheta=T+c_nCpT$. With $q_0=1$ and $\vartheta_1=T_s+Rq_1=Lq_1$ ($T_s=q_1/S$, B-analytic
Lemma 1) one has $T_1=Lq_1/N_b$ and $T_0/q_0=\mathrm{Mob}(M;T_1/q_1)$: (a); (b) multiplies by $N_0$. In (c) the transducer bottom
temperature is $T_A=R_1q+\vartheta_{\rm film}(0)=R_1q+N_tT_{\rm film}(0)$, so $T_A/q=R_1+N_t\,\mathrm{Mob}(M;L_2/N_b)$, and the
transducer acts by its known quadrupole. $\square$

At $k=0$ the extra GK conditions of the twin (slip) involve only the solenoidal (S) family, which has $T=0$ and $Q_z=0$ there,
decouples and is unforced; the longitudinal family is the 1D problem above (twin docstring, "At $k=0$ the two families
decouple").

### 1.2 Which temperature the twin observes

**PROPOSITION O.** In `src/twin_fdtr.py` (commit as read on 12 November):
(i) `kernel`, `gaussian_response`, `gaussian_response_adaptive` and `plane_response` have `observable="T"` by default, and
`synthetic_measurement` calls `gaussian_response` without an observable argument; the returned `impedance` is assembled from
the $T$-component `comp[0]` of the modes at the top face (`solve_kernel`, lines 513–516, 528). The read-out of every synthetic
measurement is therefore $T(0)$.
(ii) The option `observable="theta_e"` returns $T-(\ell^2/\lambda)(Q_z'+\alpha D)$ at the top face (line 526), i.e.
$(1+\tau_\ell p)T(0)$ at $k=0$, whatever the declared interface convention; with a Fourier transducer it equals $T(0)$.
(iii) The jump rows use `"theta"` with the layer's `GKBoundary` on the GK side of **both** faces of a GK layer (lines 489–492,
430–432: $\vartheta=T+c_n(CpT+kQ_r)-c_tkQ_r$), with `coefficients` (lines 221–228): $c_n=0$ (local), $(1+\alpha)\ell^2/\lambda$
(entropic), supplied (kinetic); $c_n=0$ for Fourier and Cattaneo media.
*Proof.* (i)–(iii) are read off the code. Numerically (`b2_01_twin_observable.py`): for local, entropic and four kinetic
coefficients ($\kappa=0.3,-0.5,2$ and the Beardo Si values), $\alpha\in\{1/3,2\}$, $(\tau_R,\tau_\ell)\in\{(25,10),(25,25),(0,20)\}$ ps,
$R\in\{0,2\cdot10^{-8}\}$ m²K/W, 49 frequencies from 10 kHz to 100 GHz, the twin's $k=0$ kernel equals Lemma R(a) to
$\le1.15\cdot10^{-15}$ (relative), its `theta_e` output equals $(1+\tau_\ell p)$ times it to $\le1.12\cdot10^{-15}$, and the Au/AlN/sapphire
stack equals Lemma R(c) with $N_t=N_b$ to $\le1.48\cdot10^{-15}$. The alternative readings ($\vartheta$ read at the heated face;
convention applied at the bottom face only) differ by at least $0.53$ and $5.6\cdot10^{-5}$ respectively whenever $c_n\ne0$. $\square$

*Consequence for Workstream A.* Cattaneo films ($\ell=0$) have $c_n=0$ in every convention (iii), so all Cattaneo results of A
are convention-free; A's GK variants used the default `GKBoundary()` (local).

### 1.3 Size of the kinetic coefficient

From Beardo et al., Phys. Rev. B 101, 075303 (2020), eq. (A6) and Table I (Si; $\gamma^{-1}=2.71\gamma_0^{-1}$, $\alpha=1/3$):
$c_n=\gamma^{-1}(\beta-\chi_{xx})$ ($x$ the interface normal, so $\chi_{nn}=\chi_{xx}$; the nonlocal term is weighted by
$\gamma_0$, not by the transducer-side $\gamma_{0\Gamma}$, as stated under Table I) and $\varepsilon=c_n\lambda/((1+\alpha)\ell^2)$ give
$\varepsilon\approx0.23,\ 0.034,\ 0.083,\ 0.067$ at 81, 154, 311, 416 K. At 311 K: $c_n=2.54\cdot10^{-17}$ m³K/W, $\tau_n=43$ ps against
$\tau_\ell=515$ ps. These are Si–Cr/Au values; no AlN value is available. The sign of $\beta-\chi_{nn}$ is not fixed in general.

---

## 2. Analytic tools

**LEMMA 3′ (meromorphy, identity principle).** For each model of Lemma R, $F(s):=Z(s^2)$, obtained from the principal-root
values on the ray $\arg s=\pi/4$ by analytic continuation, is meromorphic on $U_\theta=\mathbb C\setminus\{\pm i\tau_\ell^{-1/2}\}$
($\mathbb C$ if $\tau_\ell=0$), and $F(s)-1/(b_ss)$ is holomorphic at $s=0$. If two parameter points give equal data on a set of
frequencies with an accumulation point $\omega_*\in[0,\infty)$ at which both responses are finite, then $F\equiv F'$ on
$U_\theta\cap U_{\theta'}$; if $D_cZ=0$ on such a set, $D_cF\equiv0$ on $U_\theta$.
*Proof.* $C,X,Y$ (B-analytic, Lemma 2), the entries of $T_A$ and $N,N_0$ are even and holomorphic on $U_\theta$; $L$ is
meromorphic. Near $s=0$, $YL\to0$, $NC\to1$, so $F=N_0(CL+NX)/(YL+NC)=1/(b_ss)+O(1)$, and the transducer Möbius map
preserves the principal part ($T_{A,21}W=O(s)$). The rest is B-analytic, Lemma 3(iii)–(iv). $\square$

*Remark (poles on the imaginary axis).* For $\varepsilon=0$ there are none (B-analytic, addendum D1). For the entropic convention
($\vartheta$ at both faces with $\tau_n=\tau_\ell$, or read-out $T(0)$ with $\tau_n=\tau_\ell$) there are none either: at $k=0$,
$q=-\lambda\vartheta'/(1+\tau_Rp)$ with $\vartheta=(1+\tau_\ell p)T$, and $\int\overline\vartheta\,(Cp T+q')\,dz=0$ gives, at $p=i\omega$,
$\mathrm{Re}\,\overline{q(0)}\vartheta(0)=C\omega^2\tau_\ell\!\int|T|^2+\int|q|^2/\lambda+R|q(e)|^2+|q(e)|^2\,\mathrm{Re}(1/S)$, which excludes a free
mode ($q(0)=0\Rightarrow q\equiv0\Rightarrow T\equiv0$). For other $\varepsilon$ the interface term
$(1-\varepsilon)\tau_\ell\,\mathrm{Re}(p\,\overline{q(e)}T(e))$ has no sign and isolated poles on the axis are not excluded; data sets
must then avoid them (they are locally finite).

**LEMMA K (kernel reduction).** Let $c$ be a real direction in the parameters of a film block $\Phi$ (any $\tau_o,\tau_n$, tied
or free; $R$ known or unknown) and write $\delta\tau_o=D_c\tau_o$, $\delta\tau_n=D_c\tau_n$, $c_R=D_cR$,
$\kappa=c_{\xi}N_0-\xi_1\,\delta\tau_o\,p$. Then $D_c\Phi\equiv0$ if and only if
$$
\text{(i)}\ \ c_b=0,\quad D_c\log g=-\frac{2\,\delta\tau_o\,p}{N_0};\qquad
\text{(O)}\ \ N\,\delta\tau_o-N_0\,\delta\tau_n=2bR\,\kappa;\qquad
\text{(E)}\ \ c_R\,NN_0=\kappa\Big(\frac b{b_s^2}-bR^2p-\frac{gN^2}b\Big),
$$
identities between rational functions of $p$ (polynomial identities after multiplication by $(1+\tau_Rp)(1+\tau_\ell p)N_0$ and
$b(1+\tau_\ell p)$). For the stack, $D_cW\equiv0$ with $W=R_1+\Phi$ iff either $c_{R_1}=0$ and (i), (O), (E), or $r\equiv0$ at
the point, $c_{R_1}=D_c\tilde B/\tilde B^2$ and $D_cu=0$.

*Proof.* (⇐) Direct differentiation (below). (⇒) Let $V$ be the open sector $\{s\ne0:-\delta_0<\arg s<\pi/2\}$, $\delta_0>0$
small. $V$ is simply connected and $g(s^2)$ has neither zeros nor poles in $V$ (they lie at $s=\pm i\tau_R^{-1/2},\pm i\tau_\ell^{-1/2}$), so
$\sqrt g$, $A=\xi_1s\sqrt g$ and $B=bs/\sqrt g$ are analytic and single-valued on $V$ (continued from the data ray), and $\tilde B$,
$u$, $r$, $\Phi$ are meromorphic there (poles only at zeros of $N$, $N_0$, $u+1$); the identities $AB=\xi_1bp$ and $A/B=\xi_1g/b$
hold on the data ray (B-analytic, Lemma 0) and hence on $V$. By Lemma 3′, $D_c\Phi\equiv0$ on $V$. Put
$Q=rE$; $D_cQ=E\varphi$, $\varphi:=D_cr-2r\,D_cA$. Differentiating $\Phi=\tilde B^{-1}(1+Q)/(1-Q)$ and multiplying by
$\tilde B^2(1-Q)^2$:
$$
-D_c\tilde B\;+\;E\,(2\tilde B\varphi)\;+\;E^2\,(r^2D_c\tilde B)\;\equiv\;0\quad\text{on }V. \tag{K.1}
$$
$D_c\tilde B,\ \tilde B\varphi,\ r^2D_c\tilde B$ are algebraic functions of $s$ (Puiseux expansions at every point and at $\infty$, hence
at most power growth along any ray). Choose a path in $V$ along which $\mathrm{Re}A\to+\infty$: if $\tau_R>0$ or $\tau_\ell=0$, the
positive real $s$-axis, where $A\sim\kappa_0s$ or $\kappa_0s^2$, $\kappa_0>0$; if $\tau_R=0<\tau_\ell$, the path
$s=i\tau_\ell^{-1/2}+\rho e^{i\pi/4}$, $\rho\downarrow0$, which lies in $V$ and along which $A^2=\xi_1^2p/(1+\tau_\ell p)\sim
c\,\rho^{-1}e^{i\pi/4}$, $c>0$, so $|A|\to\infty$ with $\arg A\to\pi/8$ or $\pi/8+\pi$; in the second case replace
$(A,B,r,E)$ by $(-A,-B,1/r,1/E)$, which leaves $\Phi$, (K.1) and the conclusions below invariant (direct check). Along the path
$|E|=e^{-2\mathrm{Re}A}$ decays faster than any power, so (K.1) gives $D_c\tilde B=o(\text{any power})$; an algebraic function with
this property has a vanishing Puiseux expansion, hence vanishes on $V$. Then $2\tilde BE\varphi\equiv0$ and, $E$ having no zeros,
$\varphi\equiv0$: $D_c\tilde B\equiv0$ and $D_cr\equiv2rD_cA$.
Translation: $\log\tilde B^2=2\log b+\log p-\log g-2\log N_0$ gives (i) (its value at $p=0$ gives $c_b=0$). With
$r=(u-1)/(u+1)$, $D_cr=2rD_cA$ is $D_cu=(u^2-1)D_cA$. Using (i), $D_c\log B=\delta\tau_op/N_0$,
$D_cA=A(c_\xi/\xi_1-\delta\tau_op/N_0)$, $D_c\log L=c_R/L$, $D_c\log N=\delta\tau_np/N$, $uA=\xi_1bpL/N$, $A/u=\xi_1gN/(bL)$:
$$
pL\,(N\delta\tau_o-N_0\delta\tau_n)+c_R\,NN_0=\Big(bpL^2-\frac{gN^2}b\Big)\kappa .
$$
With $L=R+1/(b_sq)$, $q=\sqrt p$: $pL=Rq^2+q/b_s$, $pL^2=R^2q^2+2Rq/b_s+1/b_s^2$, and all other factors are even in $q$; the odd
part is (O), the even part minus $Rq^2\times$(O) is (E). Each step is reversible, which gives (⇐).
Stack: $D_cW=c_{R_1}+D_c\Phi$; multiplying by $\tilde B^2(1-Q)^2/\tilde B$ as above gives
$c_{R_1}\tilde B(1-Q)^2-(D_c\tilde B/\tilde B)(1-Q^2)+2D_cQ\equiv0$; the same asymptotic argument gives
$c_{R_1}\tilde B^2=D_c\tilde B$ and then $\varphi=c_{R_1}\tilde Br(1-rE)$, i.e. $c_{R_1}\tilde Br^2E$ algebraic, so $c_{R_1}r\equiv0$.
If $r\not\equiv0$: $c_{R_1}=0$ and the bare conditions. If $r\equiv0$ at the point: $Q\equiv0$, $D_cQ=ED_cr$, and the argument gives
$c_{R_1}=D_c\tilde B/\tilde B^2$, $D_cr=0$, i.e. $D_cu=0$. $\square$

*Remark.* This is an elementary replacement for Lemma LW in all rank statements (here and in B-analytic Theorem 2′; referee
item B-m4 for the rank part).

**LEMMA G (global reduction).** Let $\Phi,\Phi'$ be film blocks at two parameter points and $\rho\in\mathbb R$. If
$\rho+\Phi\equiv\Phi'$ on $U_\theta\cap U_{\theta'}$, then either (a) $A=A'$, $\tilde B=\tilde B'$, $r=r'$ and $\rho=0$; or
(b) $r\equiv r'\equiv0$ and $\rho+1/\tilde B=1/\tilde B'$.
*Proof.* Clear the denominators $(1-rE)(1-r'E')\tilde B\tilde B'$ (not identically zero: at $p\to0^+$, $E\to1$,
$r\to(b-b_s)/(b+b_s)$):
$(\rho\tilde B\tilde B'+\tilde B'-\tilde B)+rE(\tilde B+\tilde B'-\rho\tilde B\tilde B')-r'E'(\tilde B+\tilde B'+\rho\tilde B\tilde B')
+rr'EE'(\rho\tilde B\tilde B'-\tilde B'+\tilde B)\equiv0$, with coefficients algebraic and meromorphic on $\mathrm{Re}\,p>0$. If
$A-A'$ is non-constant the four exponents are pairwise distinct modulo constants and Lemma LW (B-analytic §6) makes all
coefficients vanish: $\rho\tilde B\tilde B'=\tilde B-\tilde B'$, then $2\tilde B'r=0$, $2\tilde Br'=0$: case (b). If $A-A'$ is constant it
is $\lim_{p\to0^+}(A-A')=0$; the exponents $1,E,E^2$ give $\rho\tilde B\tilde B'=\tilde B-\tilde B'$, $r\tilde B'=r'\tilde B$ and
$rr'(\tilde B-\tilde B')=0$; either $\tilde B=\tilde B'$ (then $\rho=0$, $r=r'$: case (a)) or $rr'=0$, whence $r=r'=0$: case (b). $\square$

---

## 3. Bare film, read-out $T(0)$ (the twin's observable)

### 3.1 Global

**THEOREM C1 (fibres).** Read-out $T(0)$ ($\tau_o=0$), bottom coefficient $\tau_n$ local ($0$), tied ($\varepsilon\tau_\ell$,
$\varepsilon\ne0$ fixed) or free; parameters $(\theta,R)$ [and $\tau_n$ if free], $R$ known or unknown. Data equal on a set with
an accumulation point (Lemma 3′) iff
(a) $r\not\equiv0$ and $(\xi_1,b,g,\tau_n,R)=(\xi_1',b',g',\tau_n',R')$, or
(b) $r\equiv0\equiv r'$; and $r\equiv0$ iff $R=0$, $b=b_s$, $\tau_n=0$, $\tau_R=\tau_\ell$.
Hence the fibres are: local — those of B-analytic Theorems 1, 1′; tied $\varepsilon\ne0$ — $\{\hat\theta\}$, except
$\{(\xi,b_s,0,0,0):\xi>0\}$ on $F_0=\{b=b_s,\tau_R=\tau_\ell=0,R=0\}$; free — $\{\hat\theta\}$ off $\Delta$,
$\{(\xi_1,b,t,t,\tau_n,R):t\ge0\}$ on $\Delta\setminus\{b=b_s,\tau_n=0,R=0\}$, $\{(\xi,b_s,t,t,0,0)\}$ on $\Delta\cap\{b=b_s,\tau_n=0,R=0\}$.
*Proof.* Lemma G with $\rho=0$, $\tilde B=B$. (a): $A^2=A'^2$ gives $\xi_1^2g=\xi_1'^2g'$ and, at $p=0$, $\xi_1=\xi_1'$, $g=g'$;
$B^2=B'^2$ gives $b=b'$; $r=r'$ gives $L/N=L'/N'$, i.e. $(R+1/(b_sq))(1+\tau_n'p)=(R'+1/(b_sq))(1+\tau_np)$, whose odd part gives
$\tau_n=\tau_n'$ and even part $R=R'$. $g=g'$ is B-analytic Theorem 1, step 5. With $\tau_n=\varepsilon\tau_\ell$, $\varepsilon\ne0$,
$\tau_n=\tau_n'$ fixes $\tau_\ell$ on $\Delta$. (b): $u\equiv1$ means $(b/b_s)(1+Rb_sq)=\sqrt g\,(1+\tau_np)$; squaring, the odd term
forces $R=0$, $p=0$ forces $b=b_s$, and $(1+\tau_Rp)(1+\tau_np)^2\equiv1+\tau_\ell p$ forces (degrees) $\tau_n=0$ and $\tau_R=\tau_\ell$;
conversely these give $u\equiv1$. In case (b) $Z=1/S$. $\square$

### 3.2 Local rank

**THEOREM C2 (rank, read-out $T(0)$).** For any frequency set as in Lemma 3′, the rank of $DZ(\hat\theta)$ is given by the table
of §0.2, with the kernels listed there; at boundary points one-sided derivatives (analytic extension) are meant.
*Proof.* Lemma K with $\tau_o=0$: $\delta\tau_o=0$, $\kappa=c_\xi$, (i) is $c_b=0$ and
$c_{\tau_R}p/(1+\tau_Rp)=c_{\tau_\ell}p/(1+\tau_\ell p)$, i.e. $c_{\tau_R}=c_{\tau_\ell}=0$ off $\Delta$ and $c_{\tau_R}=c_{\tau_\ell}=:c_t$ on $\Delta$;
(O) is $-\delta\tau_n=2bRc_\xi$; (E) is $c_RN=c_\xi(b/b_s^2-bR^2p-gN^2/b)$.
*Tied, off $\Delta$:* $\delta\tau_n=\varepsilon c_{\tau_\ell}=0$, so $Rc_\xi=0$. If $c_\xi=0$ then $c_R=0$. If $c_\xi\ne0$ then $R=0$ and
$gN^2=b^2/b_s^2-(bc_R/c_\xi)N$ is a polynomial of degree $\le1$; $gN^2=(1+\tau_Rp)(1+\varepsilon\tau_\ell p)^2/(1+\tau_\ell p)$ is a
polynomial only if $\tau_\ell=0$ (then $N=1$, the right side is constant and $\tau_R=0$: on $\Delta$) or $\varepsilon=1$; for $\varepsilon=1$
it is $(1+\tau_Rp)(1+\tau_\ell p)$, of degree $\le1$ only if $\tau_R=0$, and evaluation at $p=-1/\tau_\ell$ gives $0=b^2/b_s^2$.
Kernel $0$. *Tied, on $\Delta$* ($g=1$, $N=1+\varepsilon tp$): the $p^2$ coefficient of (E) is $-c_\xi\varepsilon^2t^2/b$. If $t>0$:
$c_\xi=0$, (O) gives $\varepsilon c_t=0$, so $c_t=0$, and $c_R=0$. If $t=0$ (corner $\Pi$): (E) gives $c_\xi bR^2=0$ ($p^1$) and
$c_R=c_\xi(b^2-b_s^2)/(bb_s^2)$ ($p^0$), (O) gives $\varepsilon c_t=-2bRc_\xi$; for $R>0$ everything vanishes; for $R=0$, $c_t=0$ and
the kernel is $\mathbb R(1,0,0,0,(b^2-b_s^2)/(bb_s^2))$ (unknown $R$), $\{0\}$ ($R$ known, $b\ne b_s$), $\mathbb Re_{\xi}$ ($R$ known, $b=b_s$).
*Free:* (O) is $-c_{\tau_n}=2bRc_\xi$. Off $\Delta$, $c_\xi\ne0$ forces $gN^2=b^2/b_s^2-b^2R^2p-(bc_R/c_\xi)N$ to be a polynomial of
degree $\le1$: either $\tau_\ell=0$, $\tau_n=0$ and then $\tau_R=-b^2R^2\le0$, i.e. $\tau_R=0$ (on $\Delta$), or $\tau_n=\tau_\ell>0$,
$\tau_R=0$, and evaluation at $p=-1/\tau_\ell$ gives $0=b^2/b_s^2+b^2R^2/\tau_\ell$; so the kernel is $0$. On $\Delta$, $c_t$ is free (the direction $(0,0,1,1)$ leaves
$g$, hence $\Phi$, unchanged), and the $p^2$ coefficient $-c_\xi\tau_n^2/b$ of (E) gives $c_\xi=0$ unless $\tau_n=0$; for $\tau_n=0$,
$c_\xi\ne0$ needs $R=0$ ($p^1$) and then $c_R=c_\xi(b^2-b_s^2)/(bb_s^2)$. *Local:* B-analytic Theorems 2, 2′ (recovered by the same
lines with $\varepsilon=0$). $\square$

### 3.3 Independent decision of the degenerate sets

The three identities of Lemma K were expanded in powers of $p$ (`kernel_system.py`), giving a matrix $M(\hat\theta)$ with
$\mathrm{rank}\,DZ=n-\dim\ker M$. (1) Exact ranks and kernels on 15 candidate strata for 14 models (`b2_02_kernel_system_run.out`,
SymPy) reproduce every entry of the tables of §0.2–0.3. (2) The **real** degenerate sets $\{\hat\theta:\mathrm{rank}\,M<n\}$ were
decided by quantifier elimination over the reals (Wolfram Language `Reduce` on all maximal minors, domain $b>0$, $\tau\ge0$,
$R\ge0$, $\varepsilon\in\mathbb R$; record in `b2_03_reduce_record.md`): read-out $T(0)$ — local, $R$ known or unknown:
$\{\tau_\ell=\tau_R\}$; tied, $R$ known: $\{\varepsilon=0,\tau_\ell=\tau_R\}\cup\{b=b_s,\tau_R=\tau_\ell=0,R=0\}$; tied, $R$ unknown:
$\{\varepsilon=0,\tau_\ell=\tau_R\}\cup\{\tau_R=\tau_\ell=0,R=0\}$; free, $R$ known or unknown: $\{\tau_\ell=\tau_R\}$, with rank $\le n-2$
exactly on $\{\tau_\ell=\tau_R,\tau_n=0,R=0\}$ (and $b=b_s$ when $R$ is known). In particular no special value of $\varepsilon$ produces a
degeneracy. (3) 50-digit function-level Jacobians (41 frequencies) on 62 bare-film points (all conventions, $\varepsilon\in\{1,0.084,
-0.5,2.5\}$) agree: retained singular values $\ge4.7\cdot10^{-6}$, deficient ones $\le7.3\cdot10^{-30}$, predicted kernels annihilated
to $\le8\cdot10^{-33}$ (`b2_04_numeric_rank.json`). (4) Ball arithmetic certifies nonzero $4\times4$ two-frequency and $5\times5$
three-frequency minors on $\Delta$, $D_2$, $\Pi$ ($R>0$), $\{\tau_R=0\}$, $\{\tau_\ell=0\}$ for the entropic convention and on $\Delta$ for
$\varepsilon=0.084$ (`b2_05_certified.json`): in particular **two generic frequencies give local rank 4 on the resonance** under the
entropic convention.

### 3.4 Corollary (answers to B-2)

**COROLLARY C3.** Read-out $T(0)$. (1) For $\kappa=0$ the degenerate set is $\Delta$ (B-analytic Theorem 2); for every fixed
$\kappa\neq0$ it is $\Pi\cap\{R=0\}$ (with $b=b_s$ if $R$ is known), a codimension-$\ge3$ set on which all conventions coincide
because $\tau_n=\varepsilon\tau_\ell=0$. (2) No value of $\kappa$ creates a new degenerate set. (3) With $c_n$ free, $\tau_n=c_nC$ is
identifiable globally (Theorem C1) and locally (it never enters a kernel), and is not confounded with $R$; only $\tau_n$, not $c_n$,
is identifiable: the gauge of B-analytic Theorem G extends by $c_n\mapsto\mu c_n$ (and $c_n=\kappa\ell^2/\lambda$ scales so).
(4) The conditioning of the lifted resonance is approximately proportional to $|\varepsilon|$ for $|\varepsilon|\lesssim1$: on $\Delta$
the smallest relative singular value is
$0.030$ ($\varepsilon=1$), $0.0027$ ($\varepsilon=0.084$), $0.016$ ($\varepsilon=-0.5$), $0.060$ ($\varepsilon=2.5$) at the test point of
`b2_04` (local: $4\cdot10^{-33}$). For Beardo-type coefficients the resonance is lifted but nearly degenerate.
*Proof.* Theorems C1, C2 and §3.3. $\square$

---

## 4. Other read-outs and face conventions

**PROPOSITION C4 (read-out $\vartheta_e(0)=(1+\tau_\ell p)T(0)$, the twin option; bottom coefficient $\varepsilon\tau_\ell$).**
Off $\Delta$ the rank is full in every case. On $\Delta$ ($\tau_R=\tau_\ell=t$):
(a) $R$ known: degenerate only on $\{b=b_s,R=0\}$ intersected with $\{\varepsilon=0\}\cup\{t=0\}$, kernel $\mathbb Re_\xi$.
(b) $R$ unknown, $\varepsilon\ne0,1$, $t>0$, $R=0$: full; $\varepsilon=1$, $t>0$: full; corner $t=0$, $R=0$: thin-slab kernel as in §0.2.
(c) $R$ unknown, $\varepsilon=0$: degenerate on $\Delta\cap\{R=0\}$ (kernel $(k_s,0,0,0,1)$ or $e_\xi$, **without** the direction
$(0,0,1,1)$) and on the exotic set
$E_0=\{t=2bR\xi_1,\ (b_s^2-b^2)\,t=b^2b_s^2R^2\}=\{R=2\xi_1(b_s^2-b^2)/(bb_s^2),\ t=4\xi_1^2(b_s^2-b^2)/b_s^2,\ 0<b<b_s\}$,
kernel $(1,0,-t/\xi_1,t/\xi_1,(b^2-b_s^2)/(bb_s^2))$, rank 4 of 5.
(d) $R$ unknown, $\varepsilon<0$: in addition the exotic set
$E_\varepsilon=\{b=b_s\sqrt{1-\varepsilon},\ t=4|\varepsilon|\xi_1^2/(1+|\varepsilon|),\ R=2|\varepsilon|\xi_1/b\}$,
kernel $(1,0,-t/\xi_1,t/\xi_1,-\varepsilon/b)$, rank 4 of 5. For $0<\varepsilon$ no exotic set exists.
*Proof.* Lemma K with $\tau_o=\tau_\ell$, $\delta\tau_o=c_{\tau_\ell}$: (i) becomes $c_{\tau_R}p/(1+\tau_Rp)+c_{\tau_\ell}p/(1+\tau_\ell p)\equiv0$, so
$c_{\tau_R}=c_{\tau_\ell}=0$ off $\Delta$ (then (O), (E) reduce to the $T(0)$ residual system, which has no solution off $\Delta$) and
$c_{\tau_R}=-c_{\tau_\ell}=:-c$ on $\Delta$. On $\Delta$, with $\kappa=c_\xi+(c_\xi t-\xi_1c)p$, (O) splits into $c(1-\varepsilon)=2bRc_\xi$ and
$bR(c_\xi t-\xi_1c)=0$; for $R>0$, $c=c_\xi t/\xi_1$, and (E) gives three equations in $(c_\xi,c_R)$ whose compatibility with
$c_\xi\ne0$ is $b^2=b_s^2(1-\varepsilon)$, $b^2R^2=\varepsilon(\varepsilon-1)t$, $(1-\varepsilon)t=2bR\xi_1$ for $\varepsilon\ne0$ (solvable iff
$\varepsilon<0$), and $(b_s^2-b^2)t=b^2b_s^2R^2$, $t=2bR\xi_1$ for $\varepsilon=0$. The case $R=0$ is the same computation with
$c(1-\varepsilon)=0$. The exotic kernels satisfy (i), (O), (E) identically (checked symbolically) and are annihilated to
$\le3\cdot10^{-33}$ by 50-digit Jacobians at $(\xi_1,b,t,R)=(1,0.6,2.56,2.133)$ ($\varepsilon=0$) and $(1,\sqrt2,2,\sqrt2)$
($\varepsilon=-1$); a point of $\Delta\cap\{R>0\}$ off $E_0$ is rank-full. Real quantifier elimination: off $\Delta$ and $R$ known — no
degeneracy; on $\Delta$ the Gröbner basis of the minors contains $b\varepsilon^2t^2(1-\varepsilon)(1-\varepsilon-b^2)$, and the three
branches $\varepsilon=0$, $\varepsilon=1$, $b^2=1-\varepsilon$ ($b_s=1$) give exactly the sets above. $\square$

**PROPOSITION C5 (film block with $\vartheta$ at both faces, same coefficient: $\tau_o=\tau_n=\varepsilon\tau_\ell$).** For every fixed
$\varepsilon\ne0$ the block $\Phi=N\,\mathrm{Mob}(M;L/N)=\mathrm{Mob}(N M N^{-1};L)$ has full rank (4 with $R$ known, 5 with $R$
unknown) except on $\Pi\cap\{R=0\}$, where the kernel is $\mathbb R(k_s,0,0,0,1)$ ($R$ unknown, $b\ne b_s$; full rank for $R$ known),
and $\mathrm{span}\{e_\xi,(0,0,1-2\varepsilon,1[,0])\}$ ($b=b_s$). With $\tau_n$ free (coordinates $(\xi_1,b,\tau_R,\tau_\ell,\tau_n,R)$):
degenerate set $\Delta$, kernel $(0,0,1,1,0,0)$, plus $(k_s,0,0,0,0,1)$ on $\Delta\cap\{\tau_n=0,R=0\}$ ($e_\xi$ if $b=b_s$), and at the
corner $\Pi\cap\{\tau_n=0,R=0,b=b_s\}$ the kernel is $\mathrm{span}\{e_\xi,(0,0,1,1,0,0),(0,0,-2,0,1,0)\}$ (there (i) leaves only
$c_{\tau_R}-c_{\tau_\ell}+2c_{\tau_n}=0$). Globally (Lemma G with $\rho=0$, $\tilde B^2=b^2p/(gN^2)$): $(\xi_1,b,g,\tau_n,R)$ are
identified when $r\not\equiv0$, so for tied $\varepsilon\ne0$ the block is injective except on $\{b=b_s\}\cap\Pi\cap\{R=0\}$.
*Proof.* (i) reads $c_b=0$, $c_{\tau_R}p/(1+\tau_Rp)-c_{\tau_\ell}p/(1+\tau_\ell p)+2\varepsilon c_{\tau_\ell}p/(1+\varepsilon\tau_\ell p)\equiv0$. The functions
$p/(1+\tau p)$ with distinct $\tau$ are independent, so $c_{\tau_R}=c_{\tau_\ell}=0$ except where two of $\tau_R,\tau_\ell,\varepsilon\tau_\ell$ coincide
with a cancelling coefficient: $\Pi$; $\Delta$ with $\varepsilon\in\{0,1\}$; $\{\tau_\ell=0,\varepsilon=1/2\}$ (exact locus of the $3\times3$
block, quantifier elimination). Off this locus the residual system in $(c_\xi,c_R)$ — the same as in Theorem C2 — is degenerate only
on $\Pi\cap\{R=0\}$ and on $\Delta\cap\{R=0,\varepsilon=0\}$. On the locus: $\Delta$ with $\varepsilon=1$ is degenerate only at $t=0$, $R=0$;
$\{\tau_\ell=0,\varepsilon=1/2\}$ only at $\tau_R=0$, $R=0$ (both by quantifier elimination and by hand: $R>0$ forces $\kappa\equiv0$, i.e.
$c_\xi=0$ and $\varepsilon c_{\tau_\ell}=0$; $R=0$ is settled by the top coefficient of (E)). At the corner, (i) gives $c_{\tau_R}=(1-2\varepsilon)c_{\tau_\ell}$
and (O), (E) give the stated kernels. Free $\tau_n$: the $4\times4$ block of (i) has determinant $4(\tau_\ell-\tau_n)(\tau_\ell-\tau_R)(\tau_n-\tau_R)/b$;
the strata $\tau_n=\tau_R$, $\tau_n=\tau_\ell$ off $\Delta$ are non-degenerate (quantifier elimination), and $\Delta$ is degenerate. $\square$

*Different coefficients at the two faces* ($\tau_o=\varepsilon_0\tau_\ell$, $\tau_n=\varepsilon\tau_\ell$, $\varepsilon_0\ne\varepsilon$; e.g. Au/AlN and
AlN/sapphire kinetic values). The same case analysis gives: degenerate set $\{\varepsilon_0=\varepsilon=0\}\cap\Delta$, the corner
$\Pi\cap\{R=0\}$ (kernel $(k_s,0,0,0,1)$, or $e_\xi$ at $b=b_s$), $\Delta\cap\{R=0,\varepsilon=0\}$, and, only for $\varepsilon_0=1$, the
exotic sets $E_0$ ($\varepsilon=0$) and $E_\varepsilon$ ($\varepsilon<0$) of Proposition C4 (the $(i)$-block forces
$c_{\tau_R}=c_{\tau_\ell}=0$ on $\Delta$ unless $\varepsilon_0\in\{0,1\}$; off $\Delta$ it forces them to vanish except on
$\{\tau_\ell=0,\varepsilon_0=1/2\}$, which is non-degenerate). The last three sets concern unknown $R$. Generic $(\varepsilon_0,\varepsilon)$:
SymPy strata (`eps0_eps_Runknown`) agree. Stack, 50 digits (`b2_04b_mixed_faces.json`): rank 5 on $E_0$ (entropic top, local
bottom), on $E_{-1}$ (entropic top, $\varepsilon=-1$ bottom), on $\Delta\cap\{R_2=0\}$ (entropic top, local bottom) and on
$\Pi\cap\{R_2=0\}$ ($(\varepsilon_0,\varepsilon)=(0.3,0.7)$, $b=2$ and $b=b_s$); rank 6 at the other points tested (on $\Delta$ off $E_0$;
$\Delta\cap\{R_2=0\}$ for $(0.3,0.7)$); predicted kernels annihilated to $\le7\cdot10^{-34}$. This case was decided by hand and on strata (not by a single
quantifier elimination).

---

## 5. B-4: known Fourier transducer, unknown $R_1$ and $R_2$

Parameters $\hat\theta=(\xi_1,b,\tau_R,\tau_\ell,R_1,R_2)\in(0,\infty)^2\times[0,\infty)^4$; transducer (Fourier layer, or any known
two-port with $\det T_A=1$) known; film convention as in the twin: the same $c_n$ at both film faces ($\tau_o=\tau_n$), local or tied
($\varepsilon\ne0$) or free.

**THEOREM T1 (reduction).** $Z_{\rm tot}=\mathrm{Mob}(T_A;W)$, $W=R_1+\Phi$ (Lemma R(c)); $Z_{\rm tot}\equiv Z_{\rm tot}'$ iff $W\equiv W'$,
and $D_cZ_{\rm tot}=D_cW/(T_{A,21}W+T_{A,22})^2$, so ranks (at every frequency, and at function level) of $Z_{\rm tot}$ and $W$
coincide. *Proof.* Möbius maps with $\det=1$ are injective, and $T_{A,21}W+T_{A,22}\not\equiv0$. $\square$

**THEOREM T2 (global).** $Z_{\rm tot}(\hat\theta)\equiv Z_{\rm tot}(\hat\theta')$ iff $R_1=R_1'$ and the film blocks coincide. Hence
(local) injective off $\Delta$; fibre $\{(\xi_1,b,t,t,R_1,R_2):t\ge0\}$ on $\Delta\setminus(D_2\cap\{R_2=0\})$ and
$\{(\xi,b_s,t,t,R_1,0):\xi>0,t\ge0\}$ on $D_2\cap\{R_2=0\}$; (tied $\varepsilon\ne0$, entropic included) injective except on
$\{b=b_s,\ \tau_R=\tau_\ell=0,\ R_2=0\}$, fibre $\{(\xi,b_s,0,0,R_1,0)\}$; (free $c_n$) as the local case with $\tau_n$ identified.
*Proof.* Lemma G with $\rho=R_1-R_1'$. Case (a) gives $\rho=0$. In case (b), $1/\tilde B-1/\tilde B'=\rho$ with $1/\tilde B^2=gN^2/(b^2p)$
rational: if $\rho\ne0$ then $2\rho/\tilde B'=1/\tilde B^2-1/\tilde B'^2-\rho^2$ is rational, so $g'N'^2/p$ is a square in $\mathbb C(p)$,
i.e. $p(1+\tau_R'p)(1+\tau_\ell'p)$ is a square in $\mathbb C[p]$, impossible ($p$ is a simple factor). So $R_1=R_1'$ in all cases,
and the film-block fibres are those of Theorem C1 with $\tilde B=B/N$: $\tilde B^2=b^2p/(gN^2)$ and $A^2=\xi_1^2pg$ give $b$, $\xi_1$, $g$,
then $N$; with $\tau_n=\varepsilon\tau_\ell$, $\varepsilon\neq0$, $N$ gives $\tau_\ell$. $r\equiv0$ iff $R_2=0$, $b=b_s$, $gN^2\equiv1$,
i.e. $\tau_n=0$, $g=1$ (degrees). $\square$

**THEOREM T3 (local rank).** The rank of $DZ_{\rm tot}(\hat\theta)$ is $1+\mathrm{rank}$ of the film block with unknown $R_2$, and the
table of §0.3 holds: local — 6 off $\Delta$, 5 on $\Delta\cap\{R_2>0\}$, 4 on $\Delta\cap\{R_2=0\}$; tied $\varepsilon\ne0$ (entropic
$\varepsilon=1$, every kinetic value) — 6 except on $\Pi\cap\{R_2=0\}$ (5 for $b\ne b_s$, 4 for $b=b_s$); free $c_n$ — 7, 6 on $\Delta$, 5 on
$\Delta\cap\{\tau_n=0,R_2=0\}$, 4 on $\Pi\cap\{\tau_n=0,R_2=0,b=b_s\}$. $R_1$ never enters a kernel.
*Proof.* Lemma K (stack part). If $r\not\equiv0$, $c_{R_1}=0$ and the film-block conditions apply: B-analytic Theorem 2′ (local),
Proposition C5 (tied, free). If $r\equiv0$ ($R_2=0$, $b=b_s$, $\tau_n=0$, $g=1$, $\tau_R=\tau_\ell=t$; $t=0$ when $\varepsilon\ne0$):
$c_{R_1}=D_c\tilde B/\tilde B^2=(2b_sq)^{-1}\big(2c_b/b_s-(c_{\tau_R}-c_{\tau_\ell})\,p/(1+tp)-2\,\delta\tau_n\,p\big)$ is a constant only if it
vanishes, so $c_{R_1}=0$, $c_b=0$ and $(c_{\tau_R}-c_{\tau_\ell})p/(1+tp)+2\delta\tau_np\equiv0$: tied, $c_{\tau_R}=(1-2\varepsilon)c_{\tau_\ell}$
($\varepsilon=0$: $c_{\tau_R}=c_{\tau_\ell}$); free, $c_{\tau_R}=c_{\tau_\ell}$, $c_{\tau_n}=0$ if $t>0$ and $c_{\tau_R}-c_{\tau_\ell}+2c_{\tau_n}=0$ if $t=0$;
$D_cu=0$ gives $c_{R_2}=0$; $c_\xi$ is free. Kernel $\mathrm{span}\{e_\xi,(0,0,1-2\varepsilon,1,0,0)\}$ (rank 4), resp. the free-$c_n$
kernels of Table 0.3 (`b2_04c_free_corner.json`: rank 4 of 7 at the corner with $b=b_s$, kernels annihilated to $4\cdot10^{-30}$). $\square$

*Checks.* 50 digits, 35 stack points (both conventions, $\varepsilon\in\{0.5,-0.5,0.084\}$, free $c_n$; all strata of the table):
all ranks as predicted, kernels annihilated to $\le6\cdot10^{-30}$ (`b2_04_numeric_rank.json`); this reproduces and extends review
script R3 (rank 6 generic, 5 on $\Delta$, 4 on $\Delta\cap\{R_2=0\}$ for the local convention). Elementary partial second proof
(low-frequency moments, no Lemma K): with $R_1$ only in $c_0$, rank $\ge1+\mathrm{rank}[\nabla c_1;\dots;\nabla c_5]$, and
$$
G'_{12345}\big|_{\rm local,\,R_2=0}=\frac{2(\tau_\ell-\tau_R)^2\xi_1^8\,[\,45(\tau_R-b^2\tau_\ell)+3(b^2-1)\xi_1^2\,]}{135\,b^2},\qquad
G'_{12345}\big|_{{\rm entropic},\,\Delta,\,R_2=0}=-\frac{4}{45}(3+b^2)\,\xi_1^{10}t^2\qquad(b_s=1),
$$
so rank 6 holds by elementary means off $\Delta\cup\{45(\tau_R-b^2\tau_\ell)+3(b^2-1)\xi_1^2=0\}$ (local, $R_2=0$) and on all of
$\Delta\cap\{R_2=0,t>0\}$ (entropic). The general-$R_2$ minors are recorded in `b2_03_reduce_record.md`.

**THEOREM T4 (finitely many frequencies).** Let $\mathcal J_N(\hat\theta;\omega)\in\mathbb R^{2N\times6}$ be the real Jacobian at $N$
frequencies. (i) $N\ge3$ is necessary for rank 6. (ii) For each stratum in Table 0.3 on which the rank is 6 and on which a
point was certified (local: both components $\tau_R>\tau_\ell$ and $\tau_R<\tau_\ell$ of the interior, $b=b_s$, $\tau_\ell=0$, $\tau_R=0$,
$R_1=0$, $R_2=0$; entropic: interior, $b=b_s$, $\tau_\ell=0$, $\tau_R=0$, $R_1=0$, $R_2=0$, $\Delta\cap\{R_2>0\}$, $\Delta\cap\{R_2=0\}$,
$D_2\cap\{R_2=0\}$, $\Pi\cap\{R_2>0\}$; and the Au/AlN/sapphire values of A in both conventions), $\mathrm{rank}\,\mathcal J_3=6$ for all
$(\hat\theta,\omega)$ in that stratum outside a closed nowhere dense null set.
(iii) At **every** $\hat\theta$ of rank 6, $\mathrm{rank}\,\mathcal J_5(\hat\theta;\omega)=6$ for $\omega$ outside a closed nowhere dense null set.
(iv) Let $G'_{12345}(\hat\theta)=\det[\nabla c_1;\dots;\nabla c_5]$ (Laurent coefficients of the film block, columns
$\xi_1,b,\tau_R,\tau_\ell,R_2$); local: $G'_{12345}=2(\tau_\ell-\tau_R)\xi_1^7P/(135b^2)$, entropic: $G'_{12345}=-2\xi_1^7Q/(2025b^2)$ ($b_s=1$;
polynomials $P,Q$ in `b2_03_reduce_record.md`; $Q\ne0$ on $\Delta\cap\{R_2=0,t>0\}$). At **every** $\hat\theta$ with
$G'_{12345}(\hat\theta)\ne0$ — local: off $\Delta\cup\{P=0\}$; entropic: off $\{Q=0\}$, a proper hypersurface that does not contain
$\Delta$ — and for all $\lambda,\mu>0$ distinct and different from 1,
$$
\det\mathcal J_3(\hat\theta;\sigma^2,\lambda^2\sigma^2,\mu^2\sigma^2)=\pm\sigma^{15}\big[t(\lambda,\mu)\,G'_{12345}(\hat\theta)+O(\sigma)\big],\qquad
t(\lambda,\mu)=\frac{(\lambda-1)^2\lambda(1+\lambda^2)(\lambda-\mu)^2(\mu-1)^2\mu(1+\mu^2)(\lambda^2+\mu^2)}{2\sqrt2},
$$
so three generic frequencies give local rank 6 at every such point.
*Proof.* (i) $2N\ge6$. (ii) Ball arithmetic (256 bits, forward-mode derivatives, frequencies $(1/2,4,17/10)$, and
$(1/10,1,10)$ in units $e^2/a=3.19$ ns for the Au/AlN/sapphire points) certifies $\det\mathcal J_3\ne0$ at one point of each listed
stratum (`b2_05_certified.json`, 21 stack points; e.g. local interior $1.5846\cdot10^{-6}\pm3\cdot10^{-16}$, Au/AlN/sapphire
Cattaneo $5.07\cdot10^{-10}\pm1.3\cdot10^{-21}$); $\det\mathcal J_3$ is real-analytic in $(\hat\theta,\omega)$ on each connected stratum (B-analytic
Lemma F, D1); conclude. (iii) The per-frequency rank is 2 for generic $\omega$: the ratio of the $R_2$ and $R_1$ columns is
$\partial_{R_2}W=1/(YL_2/N+C)^2=1-2\xi_1b\,e^{i\pi/4}\sqrt\omega/b_s+O(\omega)$, not real for small $\omega>0$, hence for generic
$\omega$. If $K$ is the common kernel after $j$ frequencies and $K\ne0$, the set of $\omega$ with $K\subset\ker\mathcal J_1(\omega)$ is a
proper analytic subset (function-level rank 6), so a generic next frequency lowers $\dim K$ by at least one: $4\to3\to2\to1\to0$.
(iv) As B-analytic Theorem N2 with repairs D2–D3: the rows of $\mathcal J_3$ for $Z_{\rm tot}$ and for $W$ differ by invertible real
$2\times2$ blocks (Theorem T1), $W-1/(b_ss)=\sum_kc_ks^k$ with $c_0\ni R_1$ only, the series converges on a disc uniform in
$\hat\theta$ near the point and may be differentiated termwise; Cauchy–Binet over 6-element index sets $K$ gives
$\det\mathcal J_3=\sum_K\sigma^{\Sigma K}t_K(\lambda,\mu)G_K$, absolutely convergent for small $\sigma$; the minimal weight $\Sigma K=15$ is
attained only by $K=\{0,\dots,5\}$, for which $G_K=\pm G'_{12345}$ (expansion along the $R_1$ column) and $t_K=t(\lambda,\mu)$; the
remainder is $O(\sigma^{16})$ by Cauchy estimates. $t(\lambda,\mu)$ and the local $G'_{12345}$ were computed in two systems
(Wolfram Language, `b2_03_reduce_record.md`; SymPy, `b2_07b_stack_minor.out`; identical); the entropic $G'_{12345}$ in general
in one system, its restriction to $\Delta\cap\{R_2=0\}$ in both. $\square$

### 5.4 Three frequencies at every point; uniqueness with three or four frequencies (numerical)

**NUMERICAL OBSERVATION (three frequencies, local rank; `b4_09_three_frequency_scan.json`).** For 40 random points on each of eight
strata (interior, $b=b_s$, $\tau_\ell=0$, $\tau_R=0$, $R_1=0$, $R_2=0$, near-resonant $\tau_\ell=1.01\tau_R$, resonant with $R_2>0$) and both
conventions, the best inverse condition number of $\mathcal J_3$ over 286 frequency triples in $[10^{-2},10^2]^3$ is at least
$3.5\cdot10^{-5}$ (local; median $\ge6.9\cdot10^{-4}$) and $1.3\cdot10^{-4}$ (entropic) on every full-rank stratum, against
$7.7\cdot10^{-10}$ (the differentiation floor) on the local resonance, where the rank is 5. This supports
**CONJECTURE T4′**: three generic frequencies give local rank 6 at every point of rank 6, including the hypersurface
$\{G'_{12345}=0\}$ left open by Theorem T4(iv) (a proof would need the next Cauchy–Binet terms on that hypersurface, as in step 3
of B-analytic Theorem N2).

**NUMERICAL OBSERVATION (global uniqueness at finitely many frequencies; `b4_06_twin_search.json`).** Multistart root
collection for $Z_{\rm tot}(\omega_j;\hat\theta)=Z_{\rm tot}(\omega_j;\hat\theta^*)$ (Levenberg–Marquardt in log parameters, 150 starts per
triple, 80 per quadruple, roots accepted at relative residual $<10^{-11}$), targets generic $(0.8,2,0.5,0.1,0.3,0.2)$, $b=b_s$
$(0.8,1,0.5,0.1,0.3,0.2)$ and AlN-like $(1,2.07,0.03,0.012,0.05,0.1)$; triples $(1/2,4,1.7)$, $(0.2,1.1,6)$; quadruple
$(1/2,4,1.7,0.2)$. Local convention: second exact roots for the generic target at both triples and for the $b=b_s$ target at
$(0.2,1.1,6)$; none for the AlN-like target; with the quadruple, only the target. Entropic convention: only the target in all
nine cases.

**PROPOSITION TW3 (three frequencies do not give global uniqueness; certified).** Local convention, transducer
$(\xi_A,b_A)=(0.2,3)$, $b_s=1$, $\hat\theta^*=(0.8,2,0.5,0.1,0.3,0.2)$, $\omega=(1/2,4,17/10)$. There is
$$
\hat\theta'=(0.985850394267,\ 1.81985674482,\ 0.380174095623,\ 0.0542762976219,\ 0.274726741561,\ 0.291556496799)\ (\pm10^{-40})
$$
with $\hat\theta'\notin\Delta$ and $Z_{\rm tot}(i\omega_j;\hat\theta')=Z_{\rm tot}(i\omega_j;\hat\theta^*)$, $j=1,2,3$: after Newton refinement at
80 digits, the Krawczyk operator (320-bit balls, forward-mode Jacobian over the box) maps the box of radius $10^{-40}$ into a
box of radius $1.4\cdot10^{-74}$ about its centre, so the box contains exactly one root, and it is regular (`b4_07_certify_twin.json`;
40-digit values there). By the implicit function theorem the twin persists for $(\hat\theta^*,\omega)$ in an open set. Two further
twins are certified in the same way, for the same target at $\omega=(0.2,1.1,6)$ and for $\hat\theta^*=(0.8,1,0.5,0.1,0.3,0.2)$ at
$(0.2,1.1,6)$. At $\omega=0.9$ the twins differ from the target by $1.2$–$5.3\cdot10^{-4}$ (relative), as Theorem T2 requires.
Hence, as for two frequencies in the bare film (B-analytic Proposition TW), **three frequencies give generic local but not global
identifiability** of the six-parameter stack.

**CONJECTURE T4″ (four frequencies, global).** For $\hat\theta$ outside a proper analytic subset and generic
$(\omega_1,\dots,\omega_4)$, the four-frequency data determine $\hat\theta$. HEURISTIC: the pairs with equal data at $N$ frequencies
form a set of expected dimension $12-2N$, whose projection is a proper subset iff $N\ge4$; consistent with the searches above.

**Function level** (`b2_06_global_search.json`, 60 frequencies in $[10^{-2},10^2]$, 100 starts): the local controls reproduce the
resonance continuum (exact twins at misfit $\le2\cdot10^{-16}$ along $\Delta$, bare film and stack); for the entropic and kinetic
($\varepsilon=0.084$, $-0.5$) conventions no second minimum below misfit $1.7\cdot10^{-2}$ (bare film on $\Delta$, $D_2$, generic) or
$1.2\cdot10^{-3}$ (stack on $\Delta$, generic) was found, as Theorems C1 and T2 predict. (In the stack case with $R_2=0$ the
"distant" exact minima differ from the target only in $\log R_2$ with $R_2\to0$, i.e. they are the target.)

---

## 6. Low-frequency structure and practical consequences

**PROPOSITION LF (low-frequency coefficients; two CAS).** For the film block $\Phi=N_0\,\mathrm{Mob}(M;L/N)$,
$\Phi=\sum_{k\ge-1}c_ks^k$, $p=s^2$:
$$
c_{-1}=\frac1{b_s},\quad c_0=R+\frac{\xi_1(b_s^2-b^2)}{bb_s^2},\quad
c_1=\frac{\xi_1^2(b^2-b_s^2)}{b_s^3}-\frac{2bR\xi_1}{b_s}-\frac{\tau_n-\tau_o}{b_s},
$$
$$
c_2\big|_{R=0}=\frac{\xi_1}{b}\Big[(\tau_R-\tau_\ell)+\frac{2b^2}{b_s^2}\tau_n-\frac{b^2-b_s^2}{b_s^2}\tau_o\Big]-\frac{\xi_1^3(b^2-b_s^2)(3b^2-b_s^2)}{3bb_s^4}.
$$
Consequences. (a) Read-out $T(0)$ ($\tau_o=0$): the interface time enters already at order $s^1$ ($c_1\ni-\tau_n/b_s$), one half-order
before the bulk combination $\tau_R-\tau_\ell$. (b) Stack ($\tau_o=\tau_n=\varepsilon\tau_\ell$): the $s^1$ terms cancel, and the leading
time combination is
$$
\delta_\varepsilon=\tau_R-\tau_\ell+\varepsilon\Big(1+\frac{b^2}{b_s^2}\Big)\tau_\ell,
$$
i.e. $\tau_R-\tau_\ell$ (local), $\tau_R+(b/b_s)^2\tau_\ell$ (entropic; $(b/b_s)^2\approx4.29$ for AlN on sapphire). This is a statement
about $\omega e^2/a\ll1$ ($e^2/a=3.2$ ns for 500 nm AlN, i.e. $f\ll50$ MHz); it does not describe what the band constrains best
(below). (c) At first order in $\omega\tau$ (any $\xi_1\sqrt\omega$) the response depends on the times through $\tau_R-\tau_\ell$ and
$\tau_n$ (and $\tau_o$), whose derivative directions are independent (Theorems C2, T3 at the corner $\Pi$; for the stack except on
$\Pi\cap\{b=b_s,R_2=0\}$): for $\varepsilon\ne0$, $\tau_R$ and $\tau_\ell$ separate structurally at **first** order in $\omega\tau$, whereas
for $\varepsilon=0$ they separate only at second order (referee B-3, item 2, is a property of the local convention).

**NUMERICAL OBSERVATION (practical conditioning in the FDTR band; `b4_08_practical.json`).** Au 80 nm / $R_1=10^{-8}$ / AlN 500 nm
($\lambda=188.6$, $C=2.41\cdot10^6$, $\tau_R=25$ ps) / $R_2=2\cdot10^{-8}$ m²K/W / sapphire; 80 frequencies 10 kHz–200 MHz; 1 % amplitude,
$\sigma_\phi$ phase noise; 1D; parameters $(\ln\lambda,\ln C,\tau_R,\tau_\ell,\ln R_1,\ln R_2)$ plus instrument gain, phase and delay.

| $\tau_\ell$ | convention | $\sigma_\phi$ | best-determined combination, $\sigma$ | $\sigma(\tau_\ell)$ | $\sigma(\tau_R)/\tau_R$ | $s_6/s_1$ (no instrument) |
|---|---|---|---|---|---|---|
| 0 | local | 0.01° | $\tau_R-1.00\,\tau_\ell$, 49 ps | 3.3 ns | 130 | $6.4\cdot10^{-3}$ |
| 0 | entropic | 0.01° | $\tau_R-0.73\,\tau_\ell$, 49 ps | 0.28 ns | 8.2 | $9.1\cdot10^{-3}$ |
| 10 ps | local / entropic | 0.01° | $\tau_R-1.00\,\tau_\ell$ / $\tau_R-0.72\,\tau_\ell$, 50 / 49 ps | 5.2 / 0.26 ns | 206 / 7.5 | $4.5\cdot10^{-3}$ / $9.5\cdot10^{-3}$ |
| 25 ps (on $\Delta$) | local / entropic | 0.01° | $\tau_R-\tau_\ell$ / $\tau_R-0.72\,\tau_\ell$, 51 / 50 ps | (rank-deficient) / 0.23 ns | $7\cdot10^6$ / 6.7 | $5\cdot10^{-9}$ / $1.0\cdot10^{-2}$ |
| any | either | 0.1° | same combinations, 350 ps | ≥ 1.4 ns | ≥ 40 | — |

With $\tau_\ell=0$ known (Cattaneo), the same set reproduces review script R8 (A's N3, 1D): $\sigma(\ln\tau_R)=1.96$ (0.01°), 14 (0.1°).
At fixed parameters the convention itself changes the stack phase by at most $0.019°$ ($\tau_\ell=10$ ps) and $0.047°$ ($\tau_\ell=25$ ps),
at 16 MHz; the resonant GK film differs from the Fourier film by $0.047°$ (entropic) and $4\cdot10^{-15}$° (local). Reading: the entropic
convention lifts the resonance degeneracy and improves the second time direction by a factor of 10–20, but in the band only one
time combination is accessible at the 50 ps (0.01°) or 350 ps (0.1°) level, against values of 10–25 ps; rank 4 (or 6) is still
practically rank 3 (or 5). The best-determined combination moves from $\tau_R-\tau_\ell$ (local) to about $\tau_R-0.7\,\tau_\ell$
(entropic), which is the convention dependence relevant to objection A-2; the Cattaneo projection used in A is convention-free.

---

## 7. Answers to the referee

**B-2.** (i) Theorems 1, 2, 1′, 2′ of B-analytic are restated as conditional on $\kappa=0$ at the interface **and** at the read-out.
(ii) Classification for $c_n=\kappa\ell^2/\lambda$, every $\kappa\neq0$ (entropic $\kappa=1+\alpha$, kinetic): Theorems C1, C2, Propositions C4,
C5, Corollary C3. (iii) $c_n$ free: identifiable, not confounded with $R$, resonance degeneracy restored (Theorem C1, C2). (iv) Which
convention follows from the conserving closure: not decided here; the entropic one is the one with a sign-definite 1D dissipation
identity (Remark after Lemma 3′), the kinetic one is fixed by a boundary distribution function (Beardo et al., App. B) and is small
($0<\varepsilon<0.25$ for Si). (v) "Fourier resonance" is a degeneracy of the local convention, not of the GK law; it should not be
presented as a physical degeneracy of measurements.
**B-4.** Proved (Theorems T1–T4) for the local and every tied convention and for a free $c_n$; the referee's numbers (rank 6/5/4) are
the local column of Table 0.3; under the entropic convention the stack is rank-full on $\Delta$ and injective off
$\{b=b_s\}\cap\Pi\cap\{R_2=0\}$. Minimal number of frequencies: three are necessary and generically sufficient for local
identifiability (at every point off an explicit hypersurface), but not for global uniqueness (certified twins); four are
conjectured to suffice globally. In the FDTR band the practical conditioning (§6) is that of a rank-5 problem: one time
combination is accessible at 50 ps (0.01°) / 350 ps (0.1°), which connects B to the N2/N3 numbers of A (R8 reproduced).

---

## 8. Status, weakest steps, open items

| Item | Status |
|---|---|
| Lemma R, Proposition O | PROVED; NUMERICALLY CONFIRMED ($\le1.5\cdot10^{-15}$) |
| Lemma K (all rank statements) | PROVED (elementary) |
| Lemma G (all global statements) | PROVED (via Lemma LW/Lemma T, audited in B-analytic) |
| Theorems C1, C2; Corollary C3 | PROVED; CAS-decided; 50-digit and certified checks |
| Propositions C4, C5 | PROVED; CAS-decided in pieces; exotic kernels exact |
| Different face coefficients | DERIVED (hand case analysis and strata), not decided by a single elimination |
| Theorems T1–T3 | PROVED; 50-digit and certified checks |
| Theorem T4 (i)–(iv) | PROVED ((ii) for the certified strata; (iv) off $\{G'_{12345}=0\}$) |
| Three frequencies at every full-rank point, including $\{G'_{12345}=0\}$ | CONJECTURED (numerical scan, §5.4) |
| Proposition TW3 (three-frequency twins, local convention) | PROVED by certified computation |
| Entropic convention: no three-frequency twin found (9 searches) | NUMERICAL OBSERVATION |
| Conjecture T4″ (four frequencies, global) | CONJECTURED |
| Proposition LF, $\delta_\varepsilon$ | DERIVED (two CAS for $c_{-1}..c_2$) |

**Weakest steps.** (1) Model: $k=0$, constant $c_n$; the S-family decoupling at $k=0$ is used (twin docstring and Lemma R); finite
beams are not covered. (2) Global statements use Lemma LW (classical; audited). (3) The exact degenerate sets for the general-$\varepsilon$
models rest on a decomposition into strata (the $(i)$-block locus, the residual system, and the coincidence strata) each decided by
real quantifier elimination; the decomposition is proved in the text, but the complete systems were not eliminated in one piece.
(4) The exotic strata of Proposition C4 are found for the read-out $\vartheta_e$ with local or negative-kinetic interfaces only; they
do not occur in the twin's configurations (read-out $T$; same coefficient at both faces).

**Open.** Three-frequency local identifiability on the hypersurface $\{G'_{12345}=0\}$ (next Cauchy–Binet terms); global
uniqueness with $N=4$ (dimension count $12-8<6$ suggests generic uniqueness; Sontag's $2r+1$ bound would give 13 real
measurements, i.e. 7 frequencies, if its hypotheses hold, which were not checked); which interface temperature the conserving
closure implies; finite beams ($k>0$) under the non-local conventions (the S family no longer decouples).

---

## 9. Verification record (`derivations/B-revision/`)

| Script / output | Content | Result |
|---|---|---|
| `b2_01_twin_observable.py/.json` | twin defaults and code checks; $k=0$ kernel vs Lemma R for 6 conventions × 2 α × 3 time pairs × 2 R; stack vs Lemma R(c) | observable $T(0)$; agreement $\le1.5\cdot10^{-15}$; alternatives excluded ($\ge5.6\cdot10^{-5}$) |
| `kernel_system.py`, `b2_02_kernel_system_run.py/.out`, `b2_02_matrices.wl` | polynomial systems of Lemma K for 14 models; exact ranks and kernels on 15 strata | all table entries reproduced |
| `b2_03_reduce_record.md` (batch version `b2_03w_reduce.wl`) | real quantifier elimination of the degenerate sets (Wolfram Language), inputs and outputs; low-frequency coefficients and moment minors | §3.3, §4, §6 |
| `models_mp.py`, `b2_04_numeric_rank.py/.json/.log` | 50-digit Jacobians, 97 points (62 bare film, 35 stack) | all consistent; retained $\ge4.7\cdot10^{-6}$, deficient $\le7.3\cdot10^{-30}$ |
| `b2_04b_mixed_faces.py/.json` | 50-digit stack Jacobians with different coefficients at the two film faces, 7 points | all consistent |
| `b2_04c_free_corner.py/.json` | 50-digit stack Jacobians, free $c_n$, Fourier corner and $\Delta\cap\{\tau_n=0,R_2=0\}$, 3 points | all consistent |
| `b2_05_certified.py/.json` | ball-arithmetic minors (256 bits, forward-mode derivatives), 31 points (10 bare film, 21 stack) | all certified nonzero |
| `b2_06_global_search.py/.json/.log` | function-level multistart for twins (60 frequencies), 10 cases | local controls: continuum found; entropic/kinetic: no twin (§5.4) |
| `b2_07_laurent_check.py/.out` | SymPy recomputation of $c_{-1}..c_2$ and $G_{0123}$ | identical to the Wolfram Language results |
| `b2_07b_stack_minor.py/.out` | SymPy recomputation of $G'_{12345}$ (local in full; entropic on $\Delta\cap\{R_2=0\}$) and $t(\lambda,\mu)$ | identical to the Wolfram Language results |
| `ad_arb.py`, `b4_07_certify_twin.py/.json` | Newton refinement (80 digits) and Krawczyk certification of the three-frequency twins | 3 of 3 certified (box $10^{-40}$ mapped into $\le3.6\cdot10^{-74}$) |
| `b4_06_twin_search.py/.json/.log` | roots of the 3- and 4-frequency stack equations, 18 cases | local: twins in 3 of 6 triple cases; entropic: none; quadruples: none |
| `b4_08_practical.py/.json` | Fisher conditioning, Au/AlN/sapphire, both conventions | §6 |
| `b4_09_three_frequency_scan.py/.json` | best three-frequency conditioning on random points of every stratum | §5.4 |

---

## Addendum after the proof audit (12 November; `branches/B-revision-audit.md`, scripts a01–a10)

Status after audit: Lemma K, C2, C5, T1, T3, T4 (i)–(iv), TW3 and Proposition O **PROVED** (T2 proved given the audited
Lemma LW); Lemma K replaces Lemma LW in every rank statement. One error and several repairs:

- **D-1 (error, C4).** At ε = 1 on Π ∩ {b = b_s, R = 0} the kernel is span{e_ξ, (0,0,−1,1[,0])}: rank 2 of 4 (R known) and
  3 of 5 (R unknown), not ℝe_ξ with ranks 3 and 4. The step "c(1−ε) = 0 ⇒ c = 0" fails at ε = 1, and at t = 0, b = b_s the
  remaining p¹ condition ξ₁c(b² − b_s²) = 0 is vacuous. Confirmed by hand, by 60-digit Jacobians (both vectors annihilated to
  8e-40) and by exact elimination on the record's matrices. The degenerate sets are unchanged. This configuration (read-out
  θ_e with the entropic convention) is the twin's optional observable, not its T(0) default, and does not affect the
  transducer results. Tables §0.4 and §8 corrected accordingly; consistent with C5 at ε = 1.
- **D-2.** The interface term is nonnegative for 0 ≤ ε ≤ 1, and for ε ≤ 0 the identity in T is sign-definite: poles on iℝ₊
  are excluded for every ε ≤ 1 (bare film) and 0 ≤ ε ≤ 1 (both-face stack), covering the Beardo values 0.03–0.23.
- **D-3.** The statement that only the entropic convention has a sign-definite 1D dissipation identity is withdrawn: at
  k = 0 the local convention satisfies Re Z_T > 0 (minimum Re Z_T/|Z_T| = 1.0e-5 over 4·10⁵ samples).
- **D-4.** Real analyticity in T1/T4(ii) requires the stack to be pole-free on iℝ₊; supplied by a layered energy identity.
- **D-5.** Off Δ, {b = b_s}, {R₁ = 0}, {R₂ = 0} have two components (τ_R ≷ τ_ℓ); only τ_R > τ_ℓ was certified for (ii); the other
  components are covered by (iv).
- **D-6.** Sign typo in T2(b): 2ρ/B̃′ = 1/B̃′² + ρ² − 1/B̃²; conclusion unchanged.
- **D-7.** The substitution (A, B, r, E) → (−A, −B, 1/r, 1/E) in Lemma K needs r ≢ 0 (holds on {τ_R = 0 < τ_ℓ}); the two forms of
  (K.1) differ by the factor r²E².
- **D-9, D-10.** The elimination record omits the t = 0 branch of C4; the piecewise eliminations decide the degenerate sets
  {rank < n} but not kernel dimensions (rank ≤ n − 2) — weakest step (3) amended. Exact scans at ε ∈ {0, ½, 1, −1, 2, ⅓} found
  no other anomaly.
- **D-11.** TW3 printed digits carry "±1e-40" for the stored 40-digit values.
