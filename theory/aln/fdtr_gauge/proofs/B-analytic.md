# B-analytic — Rank, fibres and gauge of the 1D non-Fourier front-face response

Campaign `20261111-fdtr-gauge-memory`, Workstream B (claims B1, B2, B3 and extensions of
`00-question.md`; assumptions B1–B3 of `02-assumptions.md`; model of note 7). Analytic branch.
Scripts and raw outputs: `derivations/B-analytic/` (§10). Units in all numerical examples:
$\xi_1^2$-type time unit and $b_s=1$ unless stated.

---

## 0. Summary

| Item | Statement | Status |
|---|---|---|
| Impedance (Lemma 1) | $Z=\dfrac{B\cosh A+S\sinh A}{B\,(B\sinh A+S\cosh A)}=\dfrac1B\,\dfrac{B+S\tanh A}{S+B\tanh A}$, $S=b_s\sqrt p$. **Erratum to the brief:** the brief's $\frac{B\sinh A+S\cosh A}{B(B\cosh A+S\sinh A)}$ equals $1/(B^2Z)$; it tends to $S/B^2$ instead of $1/S$ as $e\to0$. The code (`src/forward_model.py::response`) and everything below use the correct form; nothing here was built on the swapped one. | PROVED, CAS-checked |
| Analytic structure (Lemma 2) | The film quadrupole is a single-valued function of $p$: **no** branch points at $-1/\tau_R$, $-1/\tau_\ell$; regular at $-1/\tau_R$; isolated essential singularity at $-1/\tau_\ell$ ($\tau_\ell>0$, $\tau_R\neq\tau_\ell$). The only branching of $Z$ is the substrate root $\sqrt p$. | PROVED |
| B2 global (Thm 1) | $Z_\theta=Z_{\theta'}$ on a frequency set with an accumulation point in $[0,\infty)$ iff $\theta'\in\mathfrak F(\theta)$: $\{\theta\}$ off $\Delta=\{\tau_R=\tau_\ell\}$; the half-line $\{(\xi_1,b,t,t):t\ge0\}$ on $\Delta\setminus D_2$; the set $\{(\xi,b_s,t,t):\xi>0,t\ge0\}$ on $D_2=\Delta\cap\{b=b_s\}$. The matched case $b=b_s$ **off** $\Delta$ is identifiable. | PROVED (two independent proofs) |
| B1 local (Thm 2) | $\Sigma=\Delta$ exactly. Rank 4 on $\Theta\setminus\Delta$ (including $b=b_s$, $\tau_R=0$, $\tau_\ell=0$); rank 3 on $\Delta\setminus D_2$, kernel $\mathbb R(0,0,1,1)$; rank 2 on $D_2$, kernel $\mathrm{span}\{(1,0,0,0),(0,0,1,1)\}$. | PROVED (three independent proofs) |
| Resistance (Thms 1′, 2′) | Parameters $(\xi_1,b,\tau_R,\tau_\ell,R)$: injective and of rank 5 off $\Delta$ (all $R\ge0$, all $b$, both boundary strata); on $\Delta$: rank 4 if $R>0$; rank 3 if $R=0$, with extra kernel $(\frac{bb_s^2}{b^2-b_s^2},0,0,0,1)$ ($b\neq b_s$) or $(1,0,0,0,0)$ ($b=b_s$). | PROVED (one proof + algebraic and numerical checks) |
| Finite $N$ (§8) | $N=1$ never suffices. $N=2$ generic frequencies give rank 4 at **every** $\theta\notin\Delta$ (Thm N2). $N=2$ does **not** give global uniqueness: certified twin pair. Resistance model: $N\ge3$ necessary; $N=3$ generic at generic points; $N=4$ at every point off $\Delta$. | PROVED; twin CERTIFIED (interval Newton) |
| B3 gauge (Thm G) | Response depends on $(e,\lambda,\rho c,\ell,\tau_R)$ only through $\theta$; off $\Phi^{-1}(\Delta)$ the non-identifiable set is exactly the orbit of $e\to\mu e,\ \lambda\to\mu\lambda,\ \rho c\to\rho c/\mu,\ \ell\to\mu\ell$ (dimension $5-4=1$, not 2). Broken (injective and full rank off $\Delta$) by known $e$, $\lambda$, $\rho c$, $a$ or $\ell>0$, and by a known velocity with $a=v^2\tau_R/3$ ($\tau_R>0$). None breaks the resonance line. | PROVED |
| Lateral diffusion | Breaks the gauge only if the in-plane conductivity is tied to the cross-plane one; with a free in-plane conductivity an extended gauge leaves the full 2D response invariant (exact counterexample). Identifiability from the beam-averaged response: not proved. | structural PROVED; identifiability CONJECTURE |

---

## 1. Setting and conventions

**ASSUMPTION M (model; B1 of `02-assumptions.md`).** One homogeneous film $0<z<e$ on a
semi-infinite Fourier substrate; 1D transfer; constant coefficients; zero initial state; no interior
source; unit flux density absorbed at $z=0$; response $T(0)$; "local" temperature convention at the
interface: continuity of $\varphi$, and continuity of $T$ (perfect contact) or a jump $R\varphi$ (contact
resistance). No additional boundary condition is imposed on the GK flux. Film law in the Laplace
variable $p$:
$$
\rho c\,p\,T+\partial_z\varphi=0,\qquad (1+\tau_Rp)\,\varphi=-\lambda\,\partial_zT+\Lambda\,\partial_z^2\varphi ,
$$
$\Lambda=(1+\alpha)\ell^2$ ($\alpha=2$ historical GK, $\alpha=1/3$ conserving closure, $\Lambda=0$ Cattaneo),
$a=\lambda/\rho c$, $\tau_\ell=\Lambda/a$. Substrate: Fourier, effusivity $b_s>0$ known.

**DEFINITION 1 (parameter sets).**
$\Theta=\{\theta=(\xi_1,b,\tau_R,\tau_\ell):\ \xi_1>0,\ b>0,\ \tau_R\ge0,\ \tau_\ell\ge0\}$,
$\Delta=\{\theta\in\Theta:\tau_R=\tau_\ell\}$ (Fourier resonance, including the Fourier point $\tau_R=\tau_\ell=0$),
$D_2=\{\theta\in\Delta:b=b_s\}$. With a contact resistance: $\hat\Theta=\Theta\times[0,\infty)\ni\hat\theta=(\theta,R)$.

**DEFINITION 2 (functions and branches).** $\sqrt{\cdot}$ is the principal root on $\mathbb C\setminus(-\infty,0]$;
$$
g=\frac{1+\tau_Rp}{1+\tau_\ell p},\quad A=\xi_1\sqrt{pg},\quad B=b\sqrt{p/g},\quad S=b_s\sqrt p,\quad
L=R+\frac1S,\quad E=e^{-2A},\quad r=\frac{BL-1}{BL+1}\ \Big(=\frac{B-S}{B+S}\text{ if }R=0\Big).
$$
Data: $Z(i\omega)$ for $\omega$ in a set $\Omega\subset(0,\infty)$ (convention $e^{+i\omega t}$: Laplace transform at $p=i\omega$).
"Rank" means the dimension of the image of the real derivative $c\mapsto D_cZ=\sum_kc_k\partial_{\theta_k}Z$,
$c\in\mathbb R^4$ (or $\mathbb R^5$), as a map into functions on $\Omega$.

**LEMMA 0 (branches on the closed right half-plane).** For $p\neq0$, $\operatorname{Re}p\ge0$, $\tau_R,\tau_\ell\ge0$:
$pg$, $p/g$, $p\in\mathbb C\setminus(-\infty,0]$. Hence $A,B,S$ are holomorphic on an open connected set
$D\supset\{\operatorname{Re}p\ge0\}\setminus\{0\}$, have positive real parts on $\{\operatorname{Re}p\ge0\}\setminus\{0\}$, and
$$
AB=\xi_1b\,p,\qquad A/B=(\xi_1/b)\,g .
$$
*Proof.* Let $\phi=\arg p\in[-\pi/2,\pi/2]$. For $\tau>0$ and $\phi\in(0,\pi/2]$, adding $1$ to $\tau p$ gives
$\arg(1+\tau p)\in(0,\phi)$; symmetrically for $\phi<0$; for $\tau=0$ the argument is $0$. So
$\arg g\in(-|\phi|,|\phi|)$ when $\phi\neq0$ and $g>0$ when $\phi=0$. Then $\arg(pg)=\phi+\arg g$ and
$\arg(p/g)=\phi-\arg g$ lie in $(0,2\phi)\subset(0,\pi)$ for $\phi>0$, in $(2\phi,0)$ for $\phi<0$, and are $0$ for
$\phi=0$. Principal roots of such numbers have arguments in $(-\pi/2,\pi/2)$, half of the above; the
arguments of $\sqrt{pg}$ and $\sqrt{p/g}$ add to $\phi$ and subtract to $\arg g$, while the moduli multiply to
$|p|$ and divide to $|g|$. Openness and connectedness: $1+\tau_\ell p\neq0$ and the cut is avoided in a disc
around each such $p$; the union with the open right half-plane is connected. $\square$

---

## 2. The front-face impedance

**LEMMA 1 (quadrupole and impedance).** Let $\lambda_{\rm eff}=(\lambda+\Lambda\rho c\,p)/(1+\tau_Rp)=\lambda/g$ and $\sigma$ any root of
$\sigma^2=\rho c\,p/\lambda_{\rm eff}$. Then $(\sigma e,\lambda_{\rm eff}\sigma)=\pm(A,B)$ (same sign), the film obeys
$$
\begin{pmatrix}T_0\\\varphi_0\end{pmatrix}=M\begin{pmatrix}T_1\\\varphi_1\end{pmatrix},\qquad
M=\begin{pmatrix}\cosh A&\sinh A/B\\ B\sinh A&\cosh A\end{pmatrix},\qquad \det M=1,
$$
independently of the sign, and with $T_1=L\varphi_1$, $\varphi_0=1$,
$$
Z=\frac{L\cosh A+\sinh A/B}{BL\sinh A+\cosh A}=\frac1B\,\frac{1+rE}{1-rE};\qquad
R=0:\quad Z=\frac{B\cosh A+S\sinh A}{B\,(B\sinh A+S\cosh A)}=\frac1B\,\frac{B+S\tanh A}{S+B\tanh A}.
$$
Limits: $\xi_1\to0$ gives $Z\to L$ (bare substrate plus contact); $\operatorname{Re}A\to\infty$ gives $Z\to1/B$ (film half-space).
*Proof.* $\partial_z\varphi=-\rho c\,pT$ gives $\Lambda\partial_z^2\varphi=-\Lambda\rho c\,p\,\partial_zT$, so $(1+\tau_Rp)\varphi=-(\lambda+\Lambda\rho c\,p)\partial_zT$, i.e.
$\varphi=-\lambda_{\rm eff}\partial_zT$ and $\partial_z^2T=\sigma^2T$ with $\sigma^2=(p/a)g$. Hence $(\sigma e)^2=\xi_1^2pg=A^2$,
$(\lambda_{\rm eff}\sigma)^2=b^2p/g=B^2$ (using $\lambda/\sqrt a=b$), and $(\sigma e)(\lambda_{\rm eff}\sigma)=e\rho c\,p=\xi_1b\,p=AB$, so the signs agree.
With $T=\alpha\cosh\sigma z+\beta\sinh\sigma z$: $(T_1,\varphi_1)^{\mathsf T}=M^{-1}(T_0,\varphi_0)^{\mathsf T}$,
$M^{-1}=\begin{pmatrix}\cosh A&-\sinh A/B\\-B\sinh A&\cosh A\end{pmatrix}$; $M$ is invariant under $(A,B)\to(-A,-B)$. In the
substrate the bounded solution is $T_1e^{-\sqrt{p/a_s}(z-e)}$ with $\operatorname{Re}\sqrt{p/a_s}>0$, so $\varphi_1=b_s\sqrt p\,T_1=ST_1$
(principal root; Lemma 0); a contact adds $R\varphi_1$ to the film-side temperature. The Möbius form is
$T_0/\varphi_0=(M_{11}L+M_{12})/(M_{21}L+M_{22})$; the reflection form follows from
$\cosh A=\tfrac12e^A(1+E)$, $\sinh A=\tfrac12e^A(1-E)$. $\square$

*Checks* (`01_symbolic_identities.out`, items 1, 8): the matrix product equals both closed forms
(residual $<10^{-180}$ at random points); $\lim_{\xi_1\to0}Z=1/(b_s\sqrt p)$; the brief's expression has limit
$b_s(1+\tau_Rp)/(b^2\sqrt p(1+\tau_\ell p))=S/B^2$ and satisfies $Z_{\rm brief}Z=1/B^2$. For $\tau_\ell=0$ the film
half-space phase is $-\pi/4+\tfrac12\arctan\omega\tau_R$ (note 7; Camacho de la Rosa et al.).

---

## 3. Analytic structure

Let $\mathcal C(w)=\sum_{n\ge0}w^n/(2n)!$, $\mathcal S(w)=\sum_{n\ge0}w^n/(2n+1)!$ (entire), so $\cosh A=\mathcal C(A^2)$,
$\sinh A/A=\mathcal S(A^2)$ for all $A\in\mathbb C$. Put $z(p)=A^2=\xi_1^2\,p\,g(p)$.

**LEMMA 2 (single-valued film quadrupole).** On $D$,
$$
M=\begin{pmatrix}C&X\\Y&C\end{pmatrix},\qquad C=\mathcal C(z),\quad X=\frac{\xi_1}{b}\,g\,\mathcal S(z),\quad Y=\xi_1b\,p\,\mathcal S(z),\qquad C^2-XY=1,
$$
and the right-hand sides are holomorphic on $\mathbb C\setminus\{-1/\tau_\ell\}$ (on $\mathbb C$ if $\tau_\ell=0$). Consequently $M$
has no branch point; $M(-1/\tau_R)=\begin{pmatrix}1&0\\-\xi_1b/\tau_R&1\end{pmatrix}$ is a regular value ($\tau_R>0$, $\tau_R\neq\tau_\ell$); if
$\tau_\ell>0$ and $\tau_R\neq\tau_\ell$, $z$ has a simple pole at $-1/\tau_\ell$ and $M$ has an isolated essential singularity
there (a transcendental entire function composed with a pole).
*Proof.* $\sinh A/B=(A/B)(\sinh A/A)$, $B\sinh A=(AB)(\sinh A/A)$ and Lemma 0;
$C^2-XY=\cosh^2A-z\mathcal S(z)^2=\cosh^2A-\sinh^2A$. $\square$

So the branch points at $-1/\tau_R$, $-1/\tau_\ell$ suggested in the brief are absent from the response:
the film enters only through even functions of $A$ at fixed $AB$. The only branching of $Z$ is the
substrate root. Numerically (`03_two_load_recovery.out` (a)): after continuing all roots along loops
around $0$, around $-1/\tau_R$, and around both, $M$ returns to itself to $10^{-50}$ while $S$ changes sign
exactly for the loops around $0$.

**LEMMA 3 ($s$-plane; identity principle).** Put $p=s^2$, $L(s)=R+1/(b_ss)$ and
$$
F_{\hat\theta}(s)=\frac{C(s^2)L(s)+X(s^2)}{Y(s^2)L(s)+C(s^2)},\qquad s\in U_\theta=\mathbb C\setminus\{\pm i\tau_\ell^{-1/2}\}\quad(U_\theta=\mathbb C\text{ if }\tau_\ell=0).
$$
(i) $F_{\hat\theta}$ is meromorphic on $U_\theta$, and $F_{\hat\theta}(s)=\frac1{b_ss}+\sum_{k\ge0}c_k(\hat\theta)s^k$ on a punctured disc
$0<|s|<\rho(\hat\theta)$; the principal part does not depend on $\hat\theta$.
(ii) $Z_{\hat\theta}(p)=F_{\hat\theta}(\sqrt p)$ on $D$; in particular $Z_{\hat\theta}(i\omega)=F_{\hat\theta}(e^{i\pi/4}\sqrt\omega)$.
(iii) If $Z_{\hat\theta}(i\omega)=Z_{\hat\theta'}(i\omega)$ for $\omega$ in a set with an accumulation point $\omega_*\in[0,\infty)$, then
$F_{\hat\theta}\equiv F_{\hat\theta'}$ on $U_\theta\cap U_{\theta'}$, so also $F_{\hat\theta}(-s)\equiv F_{\hat\theta'}(-s)$ and $Z_{\hat\theta}\equiv Z_{\hat\theta'}$ on $D$.
(iv) Likewise, if $D_cZ_{\hat\theta}(i\omega)=0$ on such a set, then $D_cF_{\hat\theta}\equiv0$ on $U_\theta$, and the Laurent
coefficients satisfy $D_cc_k(\hat\theta)=0$ for all $k$.
*Proof.* (i) $C,X,Y$ are holomorphic in $s$ on $U_\theta$ (Lemma 2) and $L$ is meromorphic. Near $s=0$,
$F-L=(X-YL^2)/(YL+C)$ with $YL=\xi_1b\,s\,\mathcal S(z)(Rs+1/b_s)$ and $YL^2=\xi_1b\,\mathcal S(z)(Rs+1/b_s)^2$ holomorphic, and
$YL+C\to1$; so $F-L$, hence $F-1/(b_ss)=R+(F-L)$, is holomorphic at $0$, with value
$c_0=R+\xi_1/b-\xi_1b/b_s^2$ there. (ii) Lemma 1,
Lemma 2 and $S=b_s\sqrt p$. (iii) $U_\theta\cap U_{\theta'}$ is $\mathbb C$ minus at most four points, hence connected;
$H=F_{\hat\theta}-F_{\hat\theta'}$ is meromorphic there and holomorphic at $0$ by (i). It vanishes at the points
$e^{i\pi/4}\sqrt\omega$, which accumulate at $s_*=e^{i\pi/4}\sqrt{\omega_*}$; $H$ cannot have a pole at $s_*$ (it vanishes on a
sequence tending to it), so it is holomorphic near $s_*$ and $H\equiv0$. (iv) The formula for $F$ is
real-analytic in $\hat\theta$ (also at $\tau=0$, $R=0$: $z$, $g$ are rational in $\tau$), $D_cF$ is meromorphic with no
principal part at $0$, and the argument of (iii) applies; its Laurent coefficients are $D_cc_k$. $\square$

$F(-s)$ is the response continued once around $p=0$ ("second sheet" of the substrate root):
physically, the impedance seen through the film when the substrate impedance $1/S$ is replaced by $-1/S$.

---

## 4. Global identifiability (perfect contact)

**THEOREM 1.** Let $R=0$ and $\theta,\theta'\in\Theta$. The following are equivalent:
(a) $Z_\theta(i\omega)=Z_{\theta'}(i\omega)$ for all $\omega$ in a set with an accumulation point in $[0,\infty)$;
(b) $Z_\theta(i\omega)=Z_{\theta'}(i\omega)$ for all $\omega>0$;
(c) $\theta'\in\mathfrak F(\theta)$, where
$$
\mathfrak F(\theta)=\{\theta\}\ (\theta\notin\Delta),\qquad \{(\xi_1,b,t,t):t\ge0\}\ (\theta\in\Delta\setminus D_2),\qquad \{(\xi,b_s,t,t):\xi>0,t\ge0\}\ (\theta\in D_2).
$$

*Proof A (substrate second sheet).* (c)⇒(b): on $\Delta$, $g\equiv1$, so $A=\xi_1\sqrt p$, $B=b\sqrt p$ do not depend on $t$
(Fourier resonance); on $D_2$ moreover $B=S$, $r\equiv0$ and $Z=1/S$ for every $\xi$. (b)⇒(a) is trivial.
(a)⇒(c):
*Step 1.* By Lemma 3, $F:=F_\theta\equiv F_{\theta'}$ on $U=U_\theta\cap U_{\theta'}$; put $\tilde F(s)=F(-s)$, $S=b_ss$. For
$v=(C,X,Y)$ of either point, the identities $F(Y+CS)=C+XS$ and $\tilde F(Y-CS)=C-XS$ read $\mathsf L(s)v(s)=0$,
$$
\mathsf L=\begin{pmatrix}FS-1&-S&F\\-\tilde FS-1&S&\tilde F\end{pmatrix}.
$$
*Step 2 (rank of $\mathsf L$).* The $2\times2$ minors are $S(S(F-\tilde F)-2)$, $2F\tilde FS+F-\tilde F$, $-S(F+\tilde F)$ (CAS item 3);
they vanish simultaneously iff $F=1/S$, $\tilde F=-1/S$; as meromorphic functions they all vanish
identically iff $F\equiv1/(b_ss)$. Now
$$
F-\frac1S=\frac{X\,(S^2-B^2)}{S\,(Y+CS)},\qquad B^2=\frac YX=\frac{b^2p}{g},
$$
and $X\to\xi_1/b\neq0$ as $s\to0$; so $F\equiv1/S$ iff $b^2/g\equiv b_s^2$ iff $g$ is constant, i.e. $g\equiv g(0)=1$ and $b=b_s$,
i.e. $\theta\in D_2$.
*Step 3 ($\theta\notin D_2$).* $\mathsf L(s)$ has rank 2 off a discrete set $\Sigma_0\subset U$; its kernel is a line and
$v(s)\neq0$ (because $C^2-XY=1$). So $v'=\lambda v$, and $1=C'^2-X'Y'=\lambda^2$. The meromorphic function $\lambda$
($=C'/C$ off the zeros of $C$) takes values in $\{\pm1\}$ on a connected set, hence is constant; $C,C'\to1$ as
$s\to0$, so $\lambda=1$ and $M_\theta\equiv M_{\theta'}$.
*Step 4 (from $M$ to $\theta$).* $\mathcal C'(0)=\tfrac12$, so $\mathcal C$ is invertible near $0$; $z,z'\to0$ as $p\to0$ and $\mathcal C(z)=\mathcal C(z')$,
hence $z=z'$ near $0$ and everywhere. Also $B^2=Y/X=B'^2$. From $zB^2=\xi_1^2b^2p^2$ and $z/B^2=(\xi_1/b)^2g^2$ with
$g(0)=1$: $\xi_1b=\xi_1'b'$ and $\xi_1/b=\xi_1'/b'$, so $\xi_1=\xi_1'$, $b=b'$, $g^2=g'^2$, $g=g'$.
*Step 5 (rational identity).* $g=g'$ iff $(1+\tau_Rp)(1+\tau_\ell'p)\equiv(1+\tau_R'p)(1+\tau_\ell p)$ iff
$\tau_R-\tau_\ell=\tau_R'-\tau_\ell'=:\delta$ and $\tau_R\tau_\ell'=\tau_R'\tau_\ell$; substituting gives $\delta(\tau_\ell'-\tau_\ell)=0$. Either $\delta\neq0$ and
$(\tau_R',\tau_\ell')=(\tau_R,\tau_\ell)$, or $\delta=0$ and both points lie on $\Delta$.
*Step 6 ($\theta\in D_2$).* Then $F\equiv1/S$, so $\theta'\in D_2$ by Step 2. $\square$

(`03_two_load_recovery.out` (b)–(c): $(C,X,Y)$ and $\theta$ are recovered from $F(s)$, $F(-s)$ alone to
$10^{-49}$, including $b=b_s$ off $\Delta$ and $\tau_R=0$.)

*Scope.* Theorems 1–2 (and 1′–2′) are structural, noise-free statements. The second-sheet values $F(-s)$
are reached by analytic continuation from exact data; with noisy data this continuation is severely
ill-conditioned, so Proof A is a uniqueness argument, not an inversion algorithm.

---

## 5. Local rank and the exact degenerate set

**THEOREM 2.** For any frequency set with an accumulation point in $[0,\infty)$ (e.g. a band), the map
$\theta\mapsto Z_\theta$ is real-analytic (into $C([\omega_a,\omega_b];\mathbb C)$ with the sup norm, uniformly on compacts), and
$$
\operatorname{rank}DZ(\theta)=4\ \ (\theta\in\Theta\setminus\Delta);\qquad 3,\ \ker=\mathbb R(0,0,1,1)\ \ (\theta\in\Delta\setminus D_2);\qquad
2,\ \ker=\mathrm{span}\{(1,0,0,0),(0,0,1,1)\}\ \ (\theta\in D_2).
$$
At boundary points ($\tau_R=0$ or $\tau_\ell=0$) $DZ$ is the derivative of the local real-analytic extension of
$\theta\mapsto Z_\theta(i\omega)$ (for fixed $\omega$, $1+i\omega\tau\ne0$ for real $\tau$ and the roots stay off the cut near
$\tau\ge0$), i.e. one-sided partial derivatives. Thus **$\Sigma=\Delta$ exactly**; none of the candidate strata
$\{b=b_s\}$, $\{\tau_R=0\}$, $\{\tau_\ell=0\}$ drops rank off $\Delta$. Real and complex ranks coincide (the
derivatives are real on $p>0$).

*Proof A.* Let $c\in\mathbb R^4$ with $D_cZ_\theta=0$ on the frequency set; by Lemma 3, $D_cF\equiv0$, so $D_c\tilde F\equiv0$.
Case $\theta\notin D_2$: $\mathsf L$ depends on $\theta$ only through $F,\tilde F$, so differentiating $\mathsf Lv=0$ gives $\mathsf L\,D_cv=0$,
$D_cv=\mu v$; differentiating $C^2-XY=1$ gives $0=2CD_cC-XD_cY-YD_cX=2\mu(C^2-XY)=2\mu$. So $D_cM\equiv0$;
$D_cC=\tfrac12\mathcal S(z)D_cz$ with $\mathcal S(z)\to1$ gives $D_cz\equiv0$, and $D_c(Y/X)=0$ gives $D_cB^2\equiv0$:
$$
\frac{2c_1}{\xi_1}+D_c\log g\equiv0,\qquad\frac{2c_2}{b}-D_c\log g\equiv0,\qquad D_c\log g=c_3\frac{p}{1+\tau_Rp}-c_4\frac{p}{1+\tau_\ell p}.
$$
$D_c\log g$ vanishes at $p=0$, so $c_1=c_2=0$ and $c_3(1+\tau_\ell p)\equiv c_4(1+\tau_Rp)$: $c_3=c_4$ and $c_3(\tau_\ell-\tau_R)=0$.
Off $\Delta$, $c=0$; on $\Delta\setminus D_2$, $c\in\mathbb R(0,0,1,1)$, which is in the kernel since $\partial_{\tau_R}g+\partial_{\tau_\ell}g=0$ at
$\tau_R=\tau_\ell$ (CAS item 5).
Case $\theta\in D_2$: from $F-1/S=X(S^2-B^2)/(S(Y+CS))$ and $S^2=B^2$ at $\theta$,
$$
D_cZ=-\frac{X\,D_c(B^2)}{S(Y+CS)}=-\frac{1-e^{-2A}}{2B}\Big(\frac{2c_2}{b}-(c_3-c_4)\frac{p}{1+tp}\Big)\qquad(\tau_R=\tau_\ell=t,\ b=b_s)
$$
(CAS item 7). The prefactor is $\not\equiv0$ and $1$, $p/(1+tp)$ are independent: image of dimension 2,
kernel $\{c_2=0,\ c_3=c_4\}$. $\square$

*Proof C (low-frequency moments; purely algebraic).* By Lemma 3(iv), $D_cZ\equiv0$ implies $D_cc_k=0$ for all
$k$, so $\operatorname{rank}DZ(\theta)\ge\operatorname{rank}[\nabla c_0;\dots;\nabla c_K]$ for every $K$. The coefficients are (SymPy and
Wolfram Language agree; `07_moment_jacobian.out`, `07w_wolfram_crosscheck.out`)
$$
c_0=\frac{\xi_1(b_s^2-b^2)}{b\,b_s^2},\quad c_1=\frac{\xi_1^2(b^2-b_s^2)}{b_s^3},\quad
c_2=\frac{\xi_1(\tau_R-\tau_\ell)}{b}-\frac{\xi_1^3(3b^4-4b^2b_s^2+b_s^4)}{3b\,b_s^4},\quad
c_3=\frac{\xi_1^2\,[\,3b_s^4(\tau_\ell-\tau_R)+(3b^4-5b^2b_s^2+2b_s^4)\xi_1^2\,]}{3b_s^5},
$$
and $c_4,c_5,c_6$ in the output file ($c_0,\dots,c_3$ depend on the times only through $\tau_R-\tau_\ell$; $c_4$ is the
first coefficient that separates them). With $G_K=\det[\nabla c_k]_{k\in K}$:
$$
G_{0123}\equiv0,\qquad G_{0124}=\frac{2\xi_1^4(b^2-b_s^2)(\tau_\ell-\tau_R)}{b^4b_s^3},\qquad G_{0234}\big|_{b=b_s}=-\frac{2\xi_1^4(\tau_R-\tau_\ell)^2}{b_s^5}.
$$
Off $\Delta$, $G_{0124}\neq0$ if $b\neq b_s$ and $G_{0234}\neq0$ if $b=b_s$: rank 4 everywhere off $\Delta$, uniformly over the
boundary strata (the $c_k$ are polynomial in $\tau_R,\tau_\ell$). $\square$

**COROLLARY 2.1 (sub-models).** (i) Cattaneo, $\tau_\ell=0$, parameters $(\xi_1,b,\tau_R)$: injective and of rank 3
everywhere except at the corner $\{b=b_s,\tau_R=0\}$, where $\xi_1$ is free and the rank is 2. (ii) The same for
$\tau_R=0$, parameters $(\xi_1,b,\tau_\ell)$. (iii) Fourier $(\xi_1,b)$: injective and of rank 2 iff $b\ne b_s$.
*Proof.* Intersect the fibres of Theorem 1 and the kernels of Theorem 2 with $\{c_4=0\}$, $\{c_3=0\}$,
$\{c_3=c_4=0\}$. $\square$

**NUMERICAL OBSERVATION (conditioning, not rank).** Near $\Delta$ the smallest relative singular value of the
function-level Jacobian (61 log-spaced frequencies, $10^{-3}$–$10^3$, log-scaled parameters, $\theta=(1,2,\tfrac12,\tfrac12(1+\delta))$) is
$0.0318\,\delta$ for $\delta=10^{-1},\dots,10^{-5}$ (`06_two_frequency_scan.out` (C)): linear in the distance to $\Delta$.
Identifiability off $\Delta$ is practically weak for $|\delta|\ll1$. At $b=b_s$ off $\Delta$ the reflection coefficient
$r=(1-\sqrt g)/(1+\sqrt g)$ is small where $g\approx1$, and the low-frequency information on $(\tau_R,\tau_\ell)$ enters
at second order in $\tau_R-\tau_\ell$ ($G_{0234}$). Rank statements say nothing about these magnitudes.

---

## 6. Second proof by exponential transcendence (covers $R\ge0$)

Let $\mathcal K$ be the field of functions meromorphic on $D_+=\{\operatorname{Re}p>0\}$ and algebraic over $\mathbb C(p)$; it is closed
under $d/dp$. $A,B,S,L,r\in\mathcal K$ ($A^2,B^2,S^2\in\mathbb C(p)$).

**LEMMA T (classical).** If $h\in\mathcal K$ is non-constant, $e^h\notin\mathcal K$.
*Proof.* Suppose $y=e^h\in\mathcal K$. The field $\mathbb C(p,h,y)$ has transcendence degree 1, hence is the function
field of a compact connected Riemann surface $X$, and the germs of $h$, $y$ at a point of $D_+$ are
germs of meromorphic functions on $X$ (standard construction of the Riemann surface of an algebraic
function, e.g. Forster, *Lectures on Riemann Surfaces*, §8). $h$ is non-constant on compact $X$, so it has a
pole $x_0$. On the connected set $X'=X\setminus(\text{poles of }h,y)$, $y=e^h$ by the identity theorem. In a chart
$t$ at $x_0$, $h=\kappa t^{-m}(1+O(t))$, $\kappa\neq0$, $m\ge1$, and $|y|\le K|t|^{-N}$; on the rays where $\kappa t^{-m}>0$,
$|e^h|\ge\exp(\tfrac12|\kappa||t|^{-m})$: contradiction. (Algebraic variant: $y'/y=h'$, so $dh/dt=(dy/dt)/y$ in a uniformizer $t$
at a pole of $h$ of order $m$; the left side has a pole of order $m+1\ge2$, the right side at most a simple pole.) $\square$

**LEMMA LW (functional Lindemann–Weierstrass, simplest form).** If $f_1,\dots,f_m\in\mathcal K$ with $f_i-f_j$
non-constant for $i\neq j$, and $a_i\in\mathcal K$ satisfy $\sum_ia_ie^{f_i}\equiv0$, then all $a_i=0$.
*Proof.* Take a nontrivial relation with the least number of nonzero coefficients; $m\ge2$; normalize
$a_1=1$. Differentiating and subtracting $f_1'$ times the relation gives $\sum_{i\ge2}(a_i'+a_i(f_i'-f_1'))e^{f_i}=0$,
which is shorter, so every coefficient vanishes; then $(a_2e^{f_2-f_1})'=0$, $e^{f_2-f_1}=\kappa/a_2\in\mathcal K$, contradicting
Lemma T. $\square$

**THEOREM 1′ (global, with resistance).** For $\hat\theta,\hat\theta'\in\hat\Theta$: $Z_{\hat\theta}=Z_{\hat\theta'}$ on a frequency set with an
accumulation point in $[0,\infty)$ iff $\hat\theta'\in\hat{\mathfrak F}(\hat\theta)$, where $\hat{\mathfrak F}(\hat\theta)=\{\hat\theta\}$ if $\tau_R\ne\tau_\ell$;
$\{(\xi_1,b,t,t,R):t\ge0\}$ if $\tau_R=\tau_\ell$ and $(b,R)\ne(b_s,0)$; $\{(\xi,b_s,t,t,0):\xi>0,t\ge0\}$ if $\tau_R=\tau_\ell$, $b=b_s$, $R=0$.
For $R=R'=0$ this is Theorem 1 (second proof).
*Proof B.* By Lemma 3, $Z\equiv Z'$ on $D_+$. With $Z=B^{-1}(1+rE)/(1-rE)$ and $|rE|<1$ on $p>0$ (so the
denominators are $\not\equiv0$), clearing denominators:
$$
(B'-B)+(B'+B)\,r\,e^{-2A}-(B'+B)\,r'\,e^{-2A'}-(B'-B)\,rr'\,e^{-2(A+A')}\equiv0 .
$$
$A$, $A'$, $A+A'$ are non-constant (positive on $p>0$, tending to $0$ as $p\to0^+$). If $A-A'$ is non-constant,
the four exponents are pairwise distinct modulo constants and Lemma LW gives $B=B'$, $r\equiv r'\equiv0$. If
$A-A'$ is constant, the constant is $\lim_{p\to0^+}(A-A')=0$; the relation becomes
$(B'-B)(1-rr'e^{-4A})+(B'+B)(r-r')e^{-2A}\equiv0$ and Lemma LW gives $B=B'$, $r=r'$.
If $A=A'$, $B=B'$, $r=r'$: as in Steps 4–5 of Theorem 1, $\xi_1=\xi_1'$, $b=b'$, $g=g'$; with $B=B'$, $r=r'$ gives $L=L'$, i.e.
$R=R'$. If $r\equiv0$: $BL\equiv1$, so $B^2=S^2/(1+RS)^2$, i.e. $(1+Rb_s\sqrt p)^2=b_s^2p/B^2\in\mathbb C(p)$, which forces $R=0$
(the cross term $2Rb_s\sqrt p\notin\mathbb C(p)$); then $B\equiv S$, i.e. $\theta\in D_2$, and the same for $\hat\theta'$. The listed points
give equal responses. $\square$

**THEOREM 2′ (local rank, with resistance).** In coordinates $(\xi_1,b,\tau_R,\tau_\ell,R)$:
$$
\operatorname{rank}=\begin{cases}5,&\tau_R\neq\tau_\ell,\\
4,&\tau_R=\tau_\ell,\ R>0,\qquad\ker=\mathbb R(0,0,1,1,0),\\
3,&\tau_R=\tau_\ell,\ R=0,\ b\neq b_s,\qquad\ker=\mathrm{span}\{(0,0,1,1,0),\ (\tfrac{bb_s^2}{b^2-b_s^2},0,0,0,1)\},\\
3,&\tau_R=\tau_\ell,\ R=0,\ b=b_s,\qquad\ker=\mathrm{span}\{(0,0,1,1,0),\ (1,0,0,0,0)\}.\end{cases}
$$
Restricted to $R=0$ and the first four coordinates this is Theorem 2 (Proof B).
*Proof B.* Let $D_cZ\equiv0$ on $D_+$, $Q=rE$. Multiplying $D_cZ=-\frac{D_cB}{B^2}\frac{1+Q}{1-Q}+\frac{2D_cQ}{B(1-Q)^2}$ by $B(1-Q)^2$ and
using $D_cQ=E(D_cr-2rD_cA)$:
$$
-\frac{D_cB}{B}+2\,(D_cr-2r\,D_cA)\,e^{-2A}+\frac{D_cB}{B}\,r^2\,e^{-4A}\equiv0 .
$$
Lemma LW (exponents $0,-2A,-4A$) gives $D_cB\equiv0$ and $D_cr\equiv2rD_cA$. From $D_c\log B=c_2/b-\tfrac12D_c\log g\equiv0$ and
$D_c\log g(0)=0$: $c_2=0$, $c_3=c_4$, $c_3(\tau_\ell-\tau_R)=0$. Then $D_cg=0$, $D_cA=c_1\sqrt{pg}$, $D_cr=c_5\partial_Rr$, and
$$
c_5\,\partial_Rr=2c_1\,r\sqrt{pg},\qquad \partial_Rr=\frac{2B}{(BL+1)^2}=r\,\frac{2B}{B^2L^2-1}\ \ (r\neq0). \tag{$\star$}
$$
If $r\equiv0$ ($R=0$, $\theta\in D_2$): $c_5=0$ and $c_1$ is free. If $r\not\equiv0$: $c_1=0$ forces $c_5=0$; if $c_1\ne0$ then $c_5\ne0$, and
with $k=c_1/c_5$, $q=\sqrt p$, $(\star)$ is equivalent to (CAS item 8)
$$
b+k\,g(q^2)\equiv k\,b^2\Big(Rq+\frac1{b_s}\Big)^2 .
$$
The left side is even in $q$, so the odd term $2kb^2Rq/b_s$ forces $R=0$; then $g=b^2/b_s^2-b/k$ is constant, so
$\tau_R=\tau_\ell$, $g=1$ and $k=bb_s^2/(b^2-b_s^2)$ (no solution if $b=b_s$). Collecting the cases gives the table. $\square$

*Remark (the $R=0$ resonance degeneracy).* On a Fourier film at $R=0$, a thin slab $d\xi$ of film material
at the interface shifts the load $1/S$ by $d\xi(1/b-b/b_s^2)$ because $p(1/S)^2=1/b_s^2$ is constant for a
semi-infinite Fourier substrate; a resistance shifts it by $R$. Hence $\partial_RZ=\frac{bb_s^2}{b_s^2-b^2}\partial_{\xi_1}Z$ at $R=0$ on
$\Delta$ (CAS item 8: ratio $-10/21$ at $b=5/2$, $b_s=1$, for all $p$). Off $\Delta$ the slab shifts the load by
$d\xi(g/b-b/b_s^2)$, not constant. The degeneracy is first-order only: Theorem 1′ gives injectivity.

*Remark (real-axis asymptotics, partial third route).* For $\tau_R>0$ or $\tau_R=\tau_\ell=0$, $A(p)\ge c\sqrt p$ as
$p\to+\infty$, so $Z-1/B=O(e^{-2c\sqrt p})$ and $Z$ determines the rational function $p/B^2=g/b^2$ exactly, hence
$(b,\tau_R,\tau_\ell)$ off $\Delta$, then $\xi_1$ and $R$ from $rE=(BZ-1)/(BZ+1)$ (exponential versus algebraic growth in $\sqrt p$).
This route fails on $\{\tau_R=0<\tau_\ell\}$, where $A\to\xi_1/\sqrt{\tau_\ell}$; Proofs A and B do not need it.

---

## 7. Interface resistance and transducer — answer to the extension question

Rank five holds at every point off the resonance, including $R=0$, $b=b_s$ and both boundary
strata (Theorem 2′), and the map is injective there (Theorem 1′). On the resonance the resistance is
identifiable when $R>0$; at $R=0$ it is first-order confounded with the thickness (Fourier
coincidence above). All eleven five-parameter strata tested numerically agree with the table; the
predicted kernel vectors are annihilated to $10^{-62}$ (`02_numerical_rank.out`).

*Algebraic cross-check (moments, `08_low_frequency_two_point.out` (c)).* With $R$ the Laurent coefficients
start with $c_0=R+\xi_1(b_s^2-b^2)/(bb_s^2)$, and the $5\times5$ moment minor
$$
G^R_{01234}=-\frac{2\xi_1^5(\tau_R-\tau_\ell)\,Q}{3b^3b_s^4},\qquad Q=3b_s^2(\tau_R-\tau_\ell)+2Rb^3\xi_1-R^2b^2b_s^2 ,
$$
equals $-2\xi_1^5(\tau_R-\tau_\ell)^2/(b^3b_s^2)\neq0$ at $R=0$ for every $b$ (including $b=b_s$): an independent algebraic proof
of rank 5 at $R=0$ off $\Delta$. For $R>0$ it vanishes on the hypersurface $Q=0$ (one positive root $R$ when
$\tau_R>\tau_\ell$; all minors containing rows $0,1,2,3$ carry the factor $Q$). On $Q=0$ (`11_resistance_Q0_check.out`), the
minors $G^R_{01245}$, $G^R_{01345}$, $G^R_{01246}$ reduce to $2R^2\xi_1^4(2b\xi_1-Rb_s^2)$ times cofactors $P_1$ (the same for the
first two) and $P_2$; the common factor is nonzero off $\Delta$, because $Q=0$ and $R=2b\xi_1/b_s^2$ together give $Q=3b_s^2(\tau_R-\tau_\ell)=0$.
At four exact rational points on $Q=0$ the rank of $[\nabla c_0;\dots;\nabla c_7]$ is 5 (that of $[\nabla c_0;\dots;\nabla c_4]$ is 4). The
residual set $\{Q=P_1=P_2=0\}$ was not analysed, so the algebraic route is complete only for $R=0$; for $R>0$
Theorem 2′ rests on Proof B. At sample points the moment rank is 4 on $\Delta$ with $R>0$, 3 on $\Delta$ with $R=0$
($b\neq b_s$ and $b=b_s$), and 5 off $\Delta$, as in Theorem 2′.

**COROLLARY 7.1 (known transducer).** If a two-port $T$ with $\det T=1$ and known entries precedes the
film (metal transducer and its contact, all known), the measured impedance is
$Z_{\rm tot}=(T_{11}Z+T_{12})/(T_{21}Z+T_{22})$, a Möbius map with known coefficients, injective on functions;
Theorems 1, 1′, 2, 2′ and N2 transfer verbatim (derivatives are multiplied by the nonvanishing
$1/(T_{21}Z+T_{22})^2$; for N2 this multiplies each frequency's complex row by a nonzero number). An unknown
transducer or top resistance is not covered here.

---

## 8. Finitely many frequencies

**DEFINITION 4.** For $\omega=(\omega_1,\dots,\omega_N)\in(0,\infty)^N$, $\mathcal J_N(\theta;\omega)\in\mathbb R^{2N\times4}$ has rows $\operatorname{Re}\partial_\theta Z(i\omega_j)$,
$\operatorname{Im}\partial_\theta Z(i\omega_j)$ (amplitude and phase: multiplying a frequency's complex row by $1/Z$ is an
invertible real $2\times2$ operation and does not change the rank).

**LEMMA F (standard).** A real-analytic function on a connected open subset of $\mathbb R^n$, not identically
zero, has a closed, nowhere dense, Lebesgue-null zero set (e.g. B. Mityagin, *The zero set of a real
analytic function*, arXiv:1512.07276). Hence if a maximal minor of $\mathcal J_N(\theta;\cdot)$ is nonzero at one
$\omega$, then $\operatorname{rank}\mathcal J_N(\theta;\omega)=4$ for all $\omega$ outside such a set ("generic frequencies"); the same holds
jointly in $(\theta,\omega)$ on each connected stratum of $\Theta$.

**PROPOSITION N1.** $N=1$ never gives rank 4; on $\Delta$, $\operatorname{rank}\mathcal J_N\le3$ for every $N$. (Trivial / Theorem 2.)

**THEOREM N2 (two frequencies suffice locally, everywhere off the resonance).** For every $\theta\in\Theta\setminus\Delta$
(all strata, including $b=b_s$, $\tau_R=0$, $\tau_\ell=0$), $\operatorname{rank}\mathcal J_2(\theta;\omega_1,\omega_2)=4$ for all $(\omega_1,\omega_2)$ outside a
closed, nowhere dense, null subset of $(0,\infty)^2$. More precisely, for every $\lambda\in(0,1)\cup(1,\infty)$ there is
$\sigma_0(\theta,\lambda)>0$ such that $\det\mathcal J_2(\theta;\sigma^2,\lambda^2\sigma^2)\neq0$ for $0<\sigma<\sigma_0$.
*Proof.* (1) By Lemma 3(i), $F_\theta(s)-\frac1{b_ss}=\sum_{k\ge0}c_ks^k$ for $|s|<\rho(\theta)$, the $c_k$ analytic in $\theta$ and the series
differentiable termwise in $\theta$; the principal part does not depend on $\theta$. With $s=e^{i\pi/4}\sigma$ ($p=i\sigma^2$):
$\operatorname{Re}\partial_\theta Z=\sum_k\cos(k\pi/4)\sigma^k\nabla c_k$, $\operatorname{Im}\partial_\theta Z=\sum_k\sin(k\pi/4)\sigma^k\nabla c_k$.
(2) Cauchy–Binet: $\mathcal J_2=\mathcal T\mathcal G$ with $\mathcal G$ the $\infty\times4$ matrix of rows $\nabla c_k$ and $\mathcal T$ the $4\times\infty$ matrix with entries
$\cos(k\pi/4)\sigma^k$, $\sin(k\pi/4)\sigma^k$, $\cos(k\pi/4)(\lambda\sigma)^k$, $\sin(k\pi/4)(\lambda\sigma)^k$; so
$\det\mathcal J_2=\sum_{K=\{k_1<\dots<k_4\}}\sigma^{|K|}t_K(\lambda)G_K$, $|K|=\sum k_j$, with (SymPy and Wolfram Language agree)
$$
t_{0124}(\lambda)=\frac{\lambda(\lambda-1)^2(\lambda+1)(\lambda^2+1)}{\sqrt2},\qquad t_{0234}(\lambda)=\frac{\lambda^2(\lambda-1)^2(\lambda+1)(\lambda^2+1)}{\sqrt2}.
$$
(3) Leading order. $|K|\ge6$, with equality only for $K=\{0,1,2,3\}$, and $G_{0123}\equiv0$. If $b\ne b_s$: the only $K$ with
$|K|=7$ is $\{0,1,2,4\}$ and $t_{0124}(\lambda)G_{0124}\neq0$. If $b=b_s$: $\nabla c_0\parallel\nabla c_1$ (both are multiples of $e_b$), so
$G_K=0$ whenever $\{0,1\}\subset K$; the $K$ with $|K|\le9$ are $\{0,1,2,3\},\{0,1,2,4\},\{0,1,2,5\},\{0,1,3,4\},\{0,1,2,6\},
\{0,1,3,5\},\{0,2,3,4\}$, all but the last containing $\{0,1\}$; so the expansion starts with
$\sigma^9t_{0234}(\lambda)G_{0234}$, $G_{0234}|_{b=b_s}=-2\xi_1^4(\tau_R-\tau_\ell)^2/b_s^5\neq0$.
(4) Remainder. Cauchy estimates give $|\nabla c_k|\le M\rho'^{-k}$ for some $0<\rho'<\rho$; with $\mu=\max(1,\lambda)$,
$|t_K\sigma^{|K|}|\le24(\mu\sigma)^{|K|}$, $|G_K|\le24M^4\rho'^{-|K|}$ and $\#\{K:|K|=m\}\le(m+1)^3$, so
$\sum_{|K|>m_0}|\sigma^{|K|}t_KG_K|\le576M^4\sum_{m>m_0}(m+1)^3(\mu\sigma/\rho')^m=O(\sigma^{m_0+1})$. Hence
$\det\mathcal J_2=\sigma^{m_0}\big(t_{K_0}(\lambda)G_{K_0}+O(\sigma)\big)$ with $(m_0,K_0)=(7,\{0,1,2,4\})$ or $(9,\{0,2,3,4\})$, nonzero
for small $\sigma$. Lemma F concludes. $\square$

*Check* (`08_low_frequency_two_point.out` (b), 80 digits, $\lambda=2$): the ratio $\det\mathcal J_2/(\sigma^{m_0}t_{K_0}G_{K_0})$ equals
$0.918,\ 0.9915,\ 0.99915$ at $\sigma=10^{-2},10^{-3},10^{-4}$ (error $\propto\sigma$) at the interior point, at $b=b_s$, at $\tau_\ell=0$, at $\tau_R=0$, at
$b=b_s$ with $\tau_\ell=0$ or $\tau_R=0$, and at an AlN-like point; the $G_K$ with $|K|\le9$ printed there show the vanishing
pattern used in step (3). Double-precision scan (`06_two_frequency_scan.out` (B)): over 150 random interior
points and 60 points on each of five boundary/matched strata, the best two-frequency inverse condition
number on a $31\times31$ grid of pairs in $[10^{-2},10^2]^2$ is at least $5\cdot10^{-4}$ (smallest at a near-resonant point,
$\tau_R/\tau_\ell=0.988$); on the resonance it is $6\cdot10^{-17}$.

*Remarks.* (i) The leading coefficient is $\propto(\tau_R-\tau_\ell)$ for $b\neq b_s$ and $\propto(\tau_R-\tau_\ell)^2$ for $b=b_s$: the two
low-frequency determinant degenerates as $\theta\to\Delta$, faster in the matched case. (ii) The proof uses
small $\omega_1,\omega_2$ as a device; it says nothing about noise. (iii) An elementary alternative gives $N=3$ at
every $\theta\notin\Delta$: at each $\omega$ with $W(i\omega)\neq0$, $W=A\partial_AZ-B\partial_BZ$, the per-frequency rank is 2 because
$\partial_{\tau_R}Z/\partial_{\tau_\ell}Z=-(1+i\omega\tau_\ell)/(1+i\omega\tau_R)\notin\mathbb R$ (CAS item 5); since the kernels have trivial total
intersection, each new frequency can be chosen to lower the common kernel dimension, $2\to1\to0$.

**PROPOSITION NR (interface resistance).** $N\ge3$ is necessary ($2N\ge5$). The $5\times5$ minors of $\mathcal J_3$ are
certified nonzero (Arb ball arithmetic, `04_certified_finite_frequency.out`) at six points covering the
interior, $R=0$, $b=b_s$, $\tau_R=0$, $\tau_\ell=0$ and an AlN-like point; by Lemma F, $N=3$ generic frequencies give
rank 5 at generic points of each of these strata. At **every** $\hat\theta$ off $\Delta$, $N=4$ generic frequencies give
rank 5 (kernel dropping $3\to2\to1\to0$ as in Remark (iii); the ratio $\partial_{\tau_R}Z/\partial_{\tau_\ell}Z$ is unchanged by $R$).

**PROPOSITION TW (two frequencies do not give uniqueness; certified).** Let $\theta^*=(1,1,7/10,3/10)$ (units
$\xi_1^2$-time and $b_s=1$; $b=b_s$, $\theta^*\notin\Delta$) and $\omega=(1/2,4)$. There is
$$
\theta'=(0.38213304403876405621,\ 0.99455364692717951823,\ 1.58320216151666891086,\ 0.77319048590574559013)\ (\pm10^{-40})
$$
with $\theta'\neq\theta^*$, $\theta'\notin\Delta$, and $Z_{\theta'}(i\omega_j)=Z_{\theta^*}(i\omega_j)$, $j=1,2$: the Krawczyk test in 300-bit ball
arithmetic proves a unique root of the four real equations in the box of radius $10^{-40}$
(`09_certified_two_frequency_twin.out`). $\theta'$ is a regular root, so by the inverse function theorem the twin
persists for all $(\theta^*,\omega_1,\omega_2)$ in an open set: two frequencies give generic **local** but not
**global** identifiability. Consistently with Theorem 1, $Z_{\theta'}(1.7\,i)-Z_{\theta^*}(1.7\,i)\approx-0.0147-0.0125\,i$ (certified
nonzero).

**NUMERICAL OBSERVATION (twins are common for $N=2$, absent for $N=3$ in the cases tried).** Multistart root
collection (`10a_two_frequency_roots_partial.out`, `10_three_frequency_global.out`; 300–800 wide log-normal
starts per case; roots accepted at rms $<10^{-11}$): with two frequencies, exact second roots exist for the
$b=b_s$ target at both pairs tried, $(1/2,4)$ (the certified $\theta'$) and $(0.2,1.1)$ (twin $\approx(0.0990,0.9663,3.946,1.032)$),
and for $\theta^*=(0.8,0.5,0.2,0.9)$ at $(1/2,4)$ (twin $\approx(0.8405,0.5114,0.2037,0.8646)$); none was found for the
interior and AlN-like targets. With three frequencies, $(1/2,4,1.7)$ and $(0.2,1.1,6)$, only the target was found
for each of the four targets. This is evidence about the explored region, not a uniqueness proof.

**CONJECTURE N3 (global, three frequencies).** For $\theta$ outside a proper analytic subset of $\Theta\setminus\Delta$ and
generic $(\omega_1,\omega_2,\omega_3)$, the three-frequency data determine $\theta$. HEURISTIC: the set of pairs $(\theta,\theta')$ with
equal data has expected dimension $8-2N$, whose projection has dimension $<4$ iff $N\ge3$. A uniform
statement with $2\cdot4+1=9$ generic real measurements is the content of E. D. Sontag, *For differential
equations with $r$ parameters, $2r+1$ experiments are enough for identification*, J. Nonlinear Sci.
(DOI 10.1007/s00332-002-0506-0; arXiv:math/0111135), whose hypotheses were not checked here.

---

## 9. The gauge (claim B3)

**DEFINITION 3.** $P=(e,\lambda,C,\Lambda,\tau_R)\in\mathcal P=(0,\infty)^3\times[0,\infty)^2$, $C=\rho c$, $\Lambda=(1+\alpha)\ell^2$, $a=\lambda/C$. (The
list $(e,a,\lambda,\rho c,\ell,\tau_R)$ of the brief has five independent entries because $a=\lambda/\rho c$.)
$$
\Phi(P)=\Big(e\sqrt{C/\lambda},\ \sqrt{\lambda C},\ \tau_R,\ \Lambda C/\lambda\Big),\qquad
g_\mu(e,\lambda,C,\Lambda,\tau_R)=(\mu e,\mu\lambda,C/\mu,\mu^2\Lambda,\tau_R),\ \mu>0 .
$$

**THEOREM G.** (i) The front-face response at $P$ is $Z_{\Phi(P)}$ (Lemma 1).
(ii) $\Phi$ is a real-analytic surjective submersion $\mathcal P\to\Theta$, and $\Phi(P)=\Phi(P')$ iff $P'=g_\mu P$ for some $\mu>0$.
(iii) If $\Phi(P)\notin\Delta$: $\{P':Z_{\Phi(P')}=Z_{\Phi(P)}\}=\{g_\mu P:\mu>0\}$ and $\ker D(Z\circ\Phi)(P)=\mathbb R(e,\lambda,-C,2\Lambda,0)$ — in the
variables $(e,a,\lambda,\ell)$ this is note 7's identity $e\partial_e+2a\partial_a+\lambda\partial_\lambda+\ell\partial_\ell=0$. Dimension count: $5-4=1$.
(iv) If $\Phi(P)\in\Delta\setminus D_2$, the fibre is the 2-dimensional set $\{g_\mu(e,\lambda,C,t\lambda/C,t):\mu>0,\ t\ge0\}$ and the rank
is 3; on $D_2$ the fibre is 3-dimensional ($\xi_1$ free) and the rank is 2.
(v) (Breaking.) Let $\mathcal P_*$ be one of $\{e=e_0\}$, $\{\lambda=\lambda_0\}$, $\{C=C_0\}$, $\{a=a_0\}$, $\{\Lambda=\Lambda_0>0\}$, or
$\mathcal P_v=\{\lambda/C=v_0^2\tau_R/3,\ \tau_R>0\}$ (known velocity with the kinetic relation $a=v^2\tau_R/3$). Then
$Z\circ\Phi|_{\mathcal P_*}$ is injective and has derivative of rank 4 at every $P\in\mathcal P_*$ with $\Phi(P)\notin\Delta$ (on
$\{\Lambda=\Lambda_0>0\}$ one has $\tau_\ell=\Lambda_0/a>0$ automatically). On $\Phi^{-1}(\Delta)$ each restriction keeps a one-parameter fibre
(two-parameter on $D_2$): **none of these data breaks the Fourier resonance**. On $\mathcal P_v$ the fibre over
$(\xi_1,b,t,t)$ is the curve $t\mapsto(\xi_1v_0\sqrt{t/3},\ bv_0\sqrt{t/3},\ b/(v_0\sqrt{t/3}),\ v_0^2t^2/3,\ t)$.
(vi) (Parametrisation.) In the coordinate $\ell$ instead of $\Lambda$, $\partial\tau_\ell/\partial\ell=0$ at $\ell=0$: Cattaneo points acquire an
artificial extra kernel direction $\partial_\ell$. $\Lambda\propto\ell^2$ is the regular coordinate.
*Proof.* (i) Lemma 1. (ii) In the coordinates $(\log e,\log\lambda,\log C,\Lambda,\tau_R)$: $\log\xi_1=\log e+\frac12\log C-\frac12\log\lambda$,
$\log b=\frac12(\log\lambda+\log C)$, $\tau_\ell=\Lambda C/\lambda$; the Jacobian has rank 4 everywhere. Given $\theta$, a choice of $e>0$
determines $C/\lambda=\xi_1^2/e^2$ and $\lambda C=b^2$, hence $(\lambda,C)$, and $\Lambda=\tau_\ell\lambda/C$; the choices form one orbit. (iii)–(iv)
Theorems 1–2 composed with (ii); $\ker D\Phi$ is spanned by $\frac{d}{d\mu}g_\mu|_{\mu=1}$, and $DZ$ is injective off $\Delta$. (v) Along
an orbit $e,\lambda,C,a,\Lambda$ scale as $\mu,\mu,\mu^{-1},\mu^2,\mu^2$, so each constraint meets each orbit once; explicit inverses,
e.g. on $\mathcal P_v$: $a=v_0^2\tau_R/3$, $\lambda=b\sqrt a$, $C=b/\sqrt a$, $e=\xi_1\sqrt a$, $\Lambda=\tau_\ell a$; on $\{\Lambda=\Lambda_0\}$: $a=\Lambda_0/\tau_\ell$. Hence
$\Phi|_{\mathcal P_*}$ is a diffeomorphism onto its image and Theorems 1–2 apply. On $\Delta$, $t=\tau_R=\tau_\ell$ is not determined by
$Z$; on $\mathcal P_v$ (resp. $\{\Lambda=\Lambda_0\}$) it moves $a$. (vi) $\tau_\ell=(1+\alpha)\ell^2C/\lambda$. $\square$

**PROPOSITION L (lateral spreading; structural).** ASSUMPTIONS: axisymmetric 2D, Hankel wavenumber
$k$, Gaussian pump/probe of known radii (known weights $W(k)$ in $H(\omega)=\int_0^\infty Z(i\omega,k)W(k)k\,dk$), Fourier
substrate with known $(a_s,b_s)$, isotropic film, and only longitudinal flux modes in the film, so that
$\sigma^2=(p/a)g+k^2$ (exact for Fourier and Cattaneo; for GK an extra modelling assumption: transverse modes
with $\partial_z^2q_T=((1+\tau_Rp)/\ell^2+k^2)q_T$ need further boundary conditions). Then
$A_k=\xi_1\sqrt{pg+ak^2}$, $B_k=(b/g)\sqrt{pg+ak^2}$, so $H$ depends on $P$ through $(\xi_1,b,\tau_R,\tau_\ell,a)$, and $a\to\mu^2a$ under
$g_\mu$: the gauge is not a symmetry of $H$. Whether $(\xi_1,b,\tau_R,\tau_\ell,a)$ is identifiable from $H$ is **not proved**
(CONJECTURE L; an analogue of Lemma 3 for the $k$-averaged response would be needed).
**COUNTEREXAMPLE (free in-plane conductivity).** Anisotropic Fourier film with unknown $\lambda_z,\lambda_r$:
$\sigma^2=(Cp+\lambda_rk^2)/\lambda_z$, so $\sigma e=\xi_1\sqrt{p+a_rk^2}$, $\lambda_z\sigma=b_z\sqrt{p+a_rk^2}$, with $\xi_1=e\sqrt{C/\lambda_z}$, $b_z=\sqrt{\lambda_zC}$,
$a_r=\lambda_r/C$ all invariant under $(e,\lambda_z,\lambda_r,C)\to(\mu e,\mu\lambda_z,\lambda_r/\mu,C/\mu)$. The full 2D response is invariant
for every beam size: lateral diffusion does **not** break the gauge unless the in-plane conductivity is
tied to the cross-plane one (isotropy, or a known ratio $\lambda_r/\lambda_z$, which scales as $\mu^{-2}$). This demotes the
unconditional form of the B3 breaking claim for lateral wavenumbers.

**REMARK S (interior source; structural).** A volumetric source $\delta^{-1}f(z/\delta)$, written in the Liouville
variable $\xi=z/\sqrt a$, has profile $f(\xi\sqrt a/\delta)$: with $\delta$ known, $\sqrt a$ enters and the gauge is not a symmetry;
with $\delta$ unknown, the extended gauge $\delta\to\mu\delta$ restores invariance. Identifiability with known $\delta$: not proved.

---

## 10. Verification record (`derivations/B-analytic/`)

| Script / output | Content | Result |
|---|---|---|
| `01_symbolic_identities` | 26 identities (impedance forms, brief's expression, single-valued $M$, two-load system and minors, log-derivatives, $\partial_{\tau_R}Z+\partial_{\tau_\ell}Z=0$ on $\Delta$, Laurent $c_0,c_1$, $D_2$ derivative, resistance algebra, thin-layer ratio) | all pass; residuals $<10^{-180}$ at random rational points |
| `02_numerical_rank` | analytic vs. numerical-differentiation Jacobian; singular values over 61 frequencies on 12 four-parameter and 11 five-parameter strata | discrepancy $4\cdot10^{-59}$; every rank as predicted; deficient singular values $<10^{-61}$, predicted kernels annihilated to $<10^{-61}$ |
| `03_two_load_recovery` | monodromy of $M$ and $S$; recovery of $(C,X,Y)$ and $\theta$ from $F(s),F(-s)$ | $M$ single-valued to $10^{-50}$; recovery to $10^{-49}$ incl. $b=b_s$, $\tau_R=0$, $\tau_\ell=0$ |
| `04_certified_finite_frequency` | Arb (256-bit) determinants: $N=2$, 4 parameters, 9 points; $N=3$, 5 parameters, 6 points | all certified nonzero (radii $\le3\cdot10^{-75}$); resonance control ball contains 0 |
| `05_injectivity_search` | 400-start least squares (log-parameters) for twins, 60 frequencies $10^{-2}$–$10^2$; 600-start root collection for $N=2$ | no function-level twin off $\Delta$: the best local minima farther than $10^{-3}$ (log-distance) from the target have rms misfit between $6\cdot10^{-6}$ (target with $R=10^{-3}$, minimum at $R\to0$) and $0.13$, none tending to 0; continuum along $\Delta$ found on the resonance; for $N=2$ an exact second root at $b=b_s$ |
| `06_two_frequency_scan` | two-frequency scan over random $\theta$ and the boundary strata; conditioning near $\Delta$ | best pair inverse condition $\ge5\cdot10^{-4}$ at all 450 points ($6\cdot10^{-17}$ on $\Delta$); $\sigma_{\min}/\sigma_{\max}=0.0318\,\delta$ near $\Delta$ |
| `07_moment_jacobian`, `07w_wolfram_crosscheck` | Laurent coefficients $c_0..c_6$, minors $G_K$ (SymPy; Wolfram Language for $c_{-1}..c_3$, $G_{0123}$, $G_{0124}$, $G_{0234}$, $t_K$) | identical in both systems |
| `08_low_frequency_two_point` (with `laurent.py`) | $t_K(\lambda)$, $G_K$ for $|K|\le9$, high-precision check of the leading term of Theorem N2, resistance moment minors | vanishing pattern as in step (3); leading-term ratio $\to1$ with error $\propto\sigma$ on 7 points; resistance minors as in §7 |
| `09_certified_two_frequency_twin` | Krawczyk certification of the $N=2$ twin | passes (box $10^{-40}$) |
| `10_three_frequency_global` (+ `10a_two_frequency_roots_partial.out`) | multistart root collection for $N=2,3$ | $N=2$: second exact roots in 3 of 7 (target, pair) cases; $N=3$: only the target in all 8 cases |
| `11_resistance_Q0_check` | resistance model on the hypersurface $Q=0$ where $G^R_{01234}$ vanishes | see §7 |

---

## 11. Status, fragile steps, open items

| Statement | Status |
|---|---|
| Lemma 0, 1, 2, 3 | PROVED |
| Theorem 1 (global, $R=0$) | PROVED twice (Proof A elementary; Proof B via Lemma T/LW) — pending audit |
| Theorem 2 (rank, $\Sigma=\Delta$, ranks 3/2) | PROVED three times (A, B, C) — pending audit |
| Theorems 1′, 2′ (resistance) | PROVED (Proof B) — pending audit; numerically cross-checked on 11 strata; algebraic (moment) proof of the rank only at $R=0$ and off $\{Q=P_1=P_2=0\}$ |
| Corollaries 2.1, 7.1 | PROVED |
| Theorem N2 ($N=2$, every $\theta\notin\Delta$) | PROVED (relies on CAS-computed $G_K$, $t_K$, agreed by two systems) — pending audit |
| Proposition NR | PROVED (generic $N=3$ via certified evaluations + Lemma F; $N=4$ everywhere) |
| Proposition TW (twin for $N=2$) | PROVED by certified computation (Krawczyk) |
| Conjecture N3 (global, $N=3$) | CONJECTURED; consistent with multistart root collection (8 target/triple cases, one root each) |
| Twins for $N=2$ beyond the certified one | NUMERICAL OBSERVATION (3 of 7 target/pair cases) |
| Theorem G (gauge, breaking) | PROVED |
| Proposition L, Remark S | structural part PROVED; identifiability CONJECTURED; counterexample PROVED |
| Numerical observation (conditioning $\propto\delta$) | NUMERICALLY SUPPORTED |

**Earliest fragile step.** In logical order: Lemma 1 is model-dependent (Assumption M: local
temperature convention and no extra GK flux condition); all results are conditional on it. Among the
mathematical steps, the least self-contained is Lemma T (identification of algebraic germs with
meromorphic functions on a compact Riemann surface), used only in Proof B; Proof A avoids it for $R=0$,
so for $R>0$ (Theorems 1′, 2′) Lemma T is load-bearing. Theorem N2 rests on symbolic minors computed
by two independent systems, not by hand.

**What would close the remaining gaps.** (1) Audit of Proofs A/B/C and of Theorem N2's remainder
estimate. (2) A hand or formal verification of the polynomial identities $G_{0123}\equiv0$, $G_{0124}$, $G_{0234}|_{b=b_s}$
from the stated $c_0..c_4$ (finite polynomial algebra; suitable for Lean `ring`/`field_simp`). (3) For
the resistance model, a second proof independent of Lemma T (e.g. extending the two-load argument:
for each trial $R$ the two loads $R\pm1/S$ determine a symmetric $M_R$, and one would show $M_R$ has film
structure only for the true $R$). (4) For lateral spreading: an identity principle for the
$k$-averaged response and a decision on the GK transverse-mode boundary conditions. (5) Conjecture N3.

**Literature pointers (for the novelty search; not verified here unless stated).** Front-face
photothermal non-identifiability under Fourier (Krapez & Rigollet, arXiv:1708.07362; Krapez 2023);
Fourier resonance (Both et al. 2016; Kovács 2018; Hennessy & Myers 2021); identification of a
symmetric reciprocal two-port from input impedances under two terminations ("two-load" methods for
transfer matrices in acoustics/network theory); unphysical sheet of a radiation condition; functional
Lindemann–Weierstrass / Ax, *On Schanuel's conjectures*, Ann. Math. 93 (1971); structural
identifiability (Bellman & Åström 1970) and finite-experiment results (Sontag, verified to exist:
DOI 10.1007/s00332-002-0506-0); inverse impedance problems for layered media (Krein strings,
Gel'fand–Levitan, Sondhi–Gopinath). Alternative names: structural/generic identifiability, Fisher
information rank, transfer-function identifiability.

**Formalisation candidates.** Step 5 of Theorem 1 (rational identity); the $2\times3$ kernel lemma of Step 3;
the polynomial minor identities of Proof C / Theorem N2. Lemma LW is formalisable in a differential
field setting; Lemma T requires more infrastructure.

---

## Addendum after the proof audit (11 November; audit in `branches/B-proof-audit.md`)

Status after audit: Theorems 1, 2, 1′, 2′, Corollary 2.1, Theorem G and Proposition TW **PROVED**; Theorem N2 and Proposition NR
**PROVED WITH REPAIRS** (below); Proposition L, structural part: **GAP** (D7). All results remain conditional on Assumption M and
exact, noise-free data. Repairs, in the order of the argument:

- **D1 (Lemma 1: Z has no zeros or poles on {Re p ≥ 0} \ {0}).** By Lemma 0, arg B − arg S = −½ arg g ∈ (−π/4, π/4), so
  Re(B/S) > 0, and Re(RB) ≥ 0; hence Re(BL) > 0, |r| < 1, |E| = |e^{−2A}| < 1 and |rE| < 1. Since
  BL sinh A + cosh A = ½ e^{A}(BL + 1)(1 − rE), it is nonzero; the same holds for the numerator. This is used for the real
  analyticity of θ ↦ Z into C([ω_a, ω_b]) (Theorem 2), for the genericity conclusions of N2 and NR on all of (0, ∞)^N
  (Lemma F requires a connected domain of analyticity), and for Z ≠ 0 in the log-sensitivities of Definition 4.
- **D2 (N2, step 2).** The Cauchy–Binet expansion over infinitely many index sets is justified by absolute convergence for
  μσ < ρ′; the bounds of step 4 give it and are moved before step 2.
- **D3 (Lemma 3(iv), N2 step 1).** The Laurent series in s converges on a disc that does not shrink as θ varies:
  F − 1/(b_s s) = R + (X − YL²)/(YL + C) with YL + C → 1 locally uniformly in θ. Term-by-term θ-differentiation is therefore
  legitimate.
- **D4.** Off Δ the rank has three proofs (A, B, C); on Δ (ranks 3 and 2) two proofs (A, B). Wording corrected accordingly.
- **D5 (TW).** The certified box (radius 10⁻⁴⁰) is centred on an 85-digit value stored in `derivations/B-analytic/`; the printed
  20 digits lie within 5·10⁻²¹ of the root (independent re-certification, radius 10⁻¹⁵, det J₂(θ′) = 3.148·10⁻³, so the root is
  regular). "Inverse function theorem" reads "implicit function theorem".
- **D6 (§8, Check).** The quoted ratios 0.918, 0.9915, 0.99915 are replaced by the computed values: 0.958, 0.9958, 0.99958
  (b = b_s) and 0.891, 0.9886, 0.99886 (AlN-like point). The conclusion (error ∝ σ) is unchanged.
- **D7 (Proposition L) — open.** "The gauge is not a symmetry of H" requires that the beam-averaged response actually depends on a;
  this is not shown. Proposition L is downgraded to a structural remark (CONJECTURE) until the k² term is computed or a
  certified evaluation is given.
- **D8 (Lemma T).** The identification of germs of h and y with meromorphic functions on a compact Riemann surface is stated via
  a primitive element of the algebraic function field and a base point not over a branch point; the algebraic version of the
  argument is given first.
- **D9 (Theorem G(v)).** "Injective" means injective off Φ⁻¹(Δ).
- **CAS caveat.** In Theorem 1, step 5, SymPy's `solve` omits the δ = 0 branch; the hand case split of the proof covers it.
