# AB-A — Asymptotic regimes and scaling (workstreams B and A)

Campaign `20261111-fdtr-gauge-memory`, branch AB-A, 11 November 2026. Independent branch: written without access to the
other branches' outputs. Scripts and raw results: `derivations/AB-asymptotic/` (`b2`–`b5` for B, `a1`–`a3` for A,
`x1`, `x2` cross-checks; JSON in `results/`). Labels: EXACT, ASYMPTOTIC (order and uniformity stated), FORMAL, HEURISTIC;
claim statuses as in the campaign ledger.

---

## 0. Conventions and a correction to the brief

$p=i\omega$, $\exp(+i\omega t)$, front-face flux and temperature, per unit area (notes/conventions.md). For the film,

$$
G(p)=\frac{1+\tau_R p}{1+\tau_\ell p}=g^2,\qquad A=\xi_1\sqrt p\,g,\qquad B=\frac{b\sqrt p}{g},\qquad S=b_s\sqrt p ,
$$
$$
Z=\frac1B\,\frac{B+S\tanh A}{S+B\tanh A}.
$$

This is `zin` of `src/forward_model.response` with $R=0$ (impedance of the film on the load $1/S$). **The brief writes
$(1/B)(B\tanh A+S)/(B+S\tanh A)$; that is the reciprocal ratio. Its thin-film limit is $S/B^2$ rather than $1/S$, so it
is not the front-face impedance.** All results below use the form above.

Notation: $R_f=\xi_1/b=e/\lambda$ (film resistance), $C_f=\xi_1 b=\rho c\,e$ (film heat capacity),
$\delta=\tau_R-\tau_\ell$, $\beta=b_s/b$, $\Gamma=(B-S)/(B+S)$, $E=e^{-2A}$, $Z_s=1/S$.

**EXACT (1D model).** The film is a Fourier film with the complex conductivity
$\lambda_{\rm eff}(p)=\lambda(1+\tau_\ell p)/(1+\tau_R p)=\lambda/G$ at fixed $C$ and $e$: $A=e\sqrt{pC/\lambda_{\rm eff}}$ and
$B=\sqrt{\lambda_{\rm eff}Cp}$, with principal branches consistent on $p=i\omega$. Every statement about the times
below follows from this substitution.

---

# Part B — identifiability by expansions

## B1. Dimensionless groups

The quantities are $\omega,\xi_1,b,b_s,\tau_R,\tau_\ell,Z$, with two independent dimensions (time and impedance).
By Buckingham there are five groups, and $Z/R_f=\mathcal F(\Omega;\varepsilon_R,\varepsilon_\ell,\beta)$:

| group | definition | meaning | AlN on sapphire, 300 K |
|---|---|---|---|
| $\Omega$ | $\omega\xi_1^2=\omega e^2/a$ | (thickness / thermal penetration depth)$^2$ | band max $2\pi\cdot200$ MHz: 0.035 (30 nm), 0.23 (100 nm), 1.5 (300 nm), 3.8 (500 nm) |
| $\varepsilon_R$ | $\tau_R/\xi_1^2=a\tau_R/e^2$ | (Cattaneo length $\sqrt{a\tau_R}$ / thickness)$^2$ = Vernotte number$^2$ | 0.19, 0.059, 0.016, 0.009 |
| $\varepsilon_\ell$ | $\tau_\ell/\xi_1^2=(1+\alpha)\ell^2/e^2$ | (nonlocal length / thickness)$^2$ | not fixed by $\tau_{\rm mem}$ (nuisance) |
| $\beta$ | $b_s/b$ | effusivity contrast | 0.75, 0.58, 0.49, 0.47 ($\Gamma_0=(1-\beta)/(1+\beta)\approx0.14$–$0.36$) |
| $\omega\tau_R$ | $\Omega\,\varepsilon_R$ | relaxation number | band max 0.0067, 0.014, 0.025, 0.035 |

Inputs: $\lambda_z\approx78,130,180,200$ W/mK (E26 film $\kappa_z$ ratios 0.26, 0.43, 0.61, about 0.69 of bulk),
$\tau_R\approx5.3,11,19.5,27.5$ ps (E26/E30), $C=2.41\times10^6$, sapphire $b_s=10\,298$ W s$^{1/2}$ m$^{-2}$K$^{-1}$
(`b3` results, `films`).

**Orderings.** (O1) $\varepsilon_R,\varepsilon_\ell\ll1$ is the realistic ordering for every campaign film. As $\omega$
increases, the regimes are (i) $\Omega\ll1$, (ii) $1\ll\Omega\ll\varepsilon^{-1}$, (iii) $\omega\tau=O(1)$ with $\Omega\gg1$,
(iv) $\omega\tau\gg1$.
(O2) $\varepsilon_R\gtrsim1$ holds for films thinner than $\sqrt{a\tau_R}\approx50$ nm. A regime (i′) then exists, a thin
film with $\omega\tau=O(1)$, in which the film is a lumped element. In both orderings the FDTR band (10 kHz–200 MHz) has
$\omega\tau\le0.035$, so only regime (i) and the (i)/(ii) crossover are sampled.

## B2. Exact reflection representation

**LEMMA B2 (EXACT).** On $p=i\omega$, $\omega>0$: $|\Gamma|<1$ and $\operatorname{Re}A>0$, and
$$
Z=\frac1B\,\frac{1+\Gamma E}{1-\Gamma E}=\frac1B\Big[1+2\sum_{n\ge1}\Gamma^nE^n\Big],\qquad
\Gamma=\frac{b-b_s g}{b+b_s g},
$$
with absolute, geometric convergence of ratio $|\Gamma|e^{-2\operatorname{Re}A}$.

*Proof.* Write $\tanh A=(1-E)/(1+E)$ and simplify. $|b-b_sg|<|b+b_sg|$ is equivalent to $\operatorname{Re}g>0$, and
$\arg g=\tfrac12[\arctan\omega\tau_R-\arctan\omega\tau_\ell]\in(-\pi/4,\pi/4)$, so $\arg A=\pi/4+\arg g\in(0,\pi/2)$.
Numerically: over 2000 random samples, $\max|\Gamma|=0.99991$ and $\min\operatorname{Re}A/|A|=1.2\times10^{-7}$ (Cattaneo
wave limit); partial sums converge geometrically (`b3`, R1).

The derivative with respect to $\xi_1$ follows exactly. **(EXACT)**
$$\partial\ln Z/\partial\ln\xi_1=-4A\Gamma E/(1-\Gamma^2E^2)$$
(checked against the analytic Jacobian to $2\times10^{-10}$, `b5`). $\xi_1$ is seen only through the reflections; if
$\Gamma\equiv0$, $Z$ does not depend on $\xi_1$.

## B3. Regime map and the combinations each regime fixes

### (i) $\Omega\ll1$ (thin film, substrate-dominated): convergent Puiseux series

**LEMMA B3 (EXACT expansion, CAS + 50-digit check).** $Z-Z_s$ is analytic in $q=\sqrt p$ at $q=0$, jointly analytic in
$\theta$, and $Z=Z_s+\sum_{k\ge0}L_k(\theta)\,p^{k/2}$ with
$$
\begin{aligned}
L_0&=R_f-\frac{C_f}{b_s^2}=\xi_1\Big(\frac1b-\frac{b}{b_s^2}\Big),\qquad
L_1=\frac{\xi_1^2}{b_s}\Big(\frac{b^2}{b_s^2}-1\Big),\\
L_2&=R_f\,\delta+\frac{\xi_1^3}{b_s}K(\beta),\qquad K(\beta)=\frac{4}{3\beta}-\frac1{\beta^3}-\frac\beta3\quad(K(1)=0),\\
L_3&=-\frac{\xi_1^2}{b_s}\,\delta+\xi_1^4\Big(\frac{b^4}{b_s^5}-\frac{5b^2}{3b_s^3}+\frac{2}{3b_s}\Big),\\
L_4&=-R_f\,\tau_\ell\,\delta+\delta\,\xi_1^3\Big(\frac{4b}{3b_s^2}-\frac{2}{3b}\Big)+\xi_1^5\Big(\frac{2}{15b}-\frac{17b}{15b_s^2}+\frac{2b^3}{b_s^4}-\frac{b^5}{b_s^6}\Big).
\end{aligned}
$$
($L_5$ is in `results/b2_lowfreq_series.json`.)

The structure is verified symbolically: $L_0,L_1$ do not depend on the times; $L_2$ and $L_3$ depend on them only through
$\delta$; $\tau_\ell$ first enters at $q^4$ through $-R_f\tau_\ell\delta$, and $L_4+R_f\tau_\ell\delta$ depends on the
times only through $\delta$. Numerically, the remainder after $K$ terms decays with slope $K/2$ in $\omega$ for
$K=0,\dots,6$ (generic, Cattaneo, $\tau\gg\xi_1^2$). At $b=b_s$ the first nonzero term is $L_2$, giving slope 1
(`b2`). Physically, $L_0$ is the classical series resistance plus shunt capacitance ($Z\simeq Z_s+R_f-pC_fZ_s^2$), and
the times enter through the series resistance $e/\lambda_{\rm eff}=R_fG(p)$.

**Fixed by (i), order by order:**
- $O(1)$: $R_f-C_f/b_s^2$.
- $O(\sqrt p)$: $C_f(b^2-b_s^2)\xi_1/b_s^3$. Together with $O(1)$ this gives $(\xi_1,b)$ explicitly when $b\neq b_s$ (B5).
- $O(p)$: $\delta=\tau_R-\tau_\ell$ through $R_f\delta$. It is the only time combination at first order.
- $O(p^2)$: $\tau_\ell$ (equivalently $\tau_R+\tau_\ell$) through $-R_f\tau_\ell\delta$, which vanishes when $\delta=0$.

### Uniform small-$\omega\tau$ statement (all $\Omega$)

**PROPOSITION B3′ (ASYMPTOTIC, uniform in $\Omega$).** Write $s_\lambda=\partial\ln Z_F/\partial\ln\lambda$. From
$\ln\lambda_{\rm eff}=\ln\lambda-\delta p+\tfrac12\delta(\tau_R+\tau_\ell)p^2+O(p^3\tau^3)$,
$$
\ln Z=\ln Z_F(\xi_1,b)-\delta\,p\,s_\lambda+\Big[\tfrac12\delta(\tau_R+\tau_\ell)\,s_\lambda+\tfrac12\delta^2\,\partial_{\ln\lambda}s_\lambda\Big]p^2+O(|p\tau|^3).
$$
At relative order $\omega\tau$ the times enter only through $\delta$. The sum $\tau_R+\tau_\ell$ enters at relative order
$(\omega\tau)^2$, multiplied by $\delta$. The coefficients are bounded functions of $\Omega$, so the statement is uniform
across regimes (i), (ii) and their crossover. For FDTR it means that only $\tau_R-\tau_\ell$ is accessible at first order
(see A7, nuisance set N4).

### (ii) $1\ll\Omega\ll\varepsilon^{-1}$ (thick film, Fourier-like)

$Z=\dfrac{g}{b\sqrt p}\big[1+2\Gamma E+O(E^2)\big]$ with $\Gamma=\Gamma_0+O(\omega\tau)$ and $\Gamma_0=(b-b_s)/(b+b_s)$.

- Phase: $\varphi=-\pi/4+\tfrac12\omega\delta-\tfrac16\omega^3(\tau_R^3-\tau_\ell^3)+\dots$
- Amplitude: $\ln|Z|=-\ln b-\tfrac12\ln\omega+\tfrac14\omega^2(\tau_R^2-\tau_\ell^2)+\dots$

This regime fixes $b$ at $O(1)$, $\delta$ at $O(\omega\tau)$ in the phase, $\tau_R+\tau_\ell$ at $O((\omega\tau)^2)$ in the
amplitude, and $\xi_1$ only through the exponentially small reflection. **ASYMPTOTIC:**
$|\partial\ln Z/\partial\ln\xi_1|\simeq4|\Gamma||A|e^{-2\operatorname{Re}A}$.

*Non-uniformity.* $2\operatorname{Re}A=\xi_1\sqrt{2\omega}\,[1-\omega\delta/2+O((\omega\tau)^2)]$. The relaxation
correction to the decay exponent is $O(1)$ once $\omega\gtrsim\omega_\times=(\sqrt2/\xi_1\delta)^{2/3}$, even though
$\omega\tau\ll1$; the Fourier exponent $\xi_1\sqrt{2\omega}$ is then wrong by a factor $e^{O(1)}$. Numerically, the
ratio of the exact term to the Fourier asymptote goes from 1.01 to 12 as $\omega\delta$ goes from 0.006 to 0.24
(`b5`). For 500 nm AlN, $\omega_\times/2\pi\approx1.5$ GHz, outside the band.

### (iii) $\omega\tau=O(1)$, $\Omega\gg1$ (semi-infinite film)

$Z=g/(b\sqrt p)+O(e^{-2\operatorname{Re}A})$, with $2\operatorname{Re}A\gtrsim\sqrt{2/\varepsilon}\approx15$ (500 nm), so
$e^{-15}\approx3\times10^{-7}$.

- $\varphi=-\pi/4+\tfrac12[\arctan\omega\tau_R-\arctan\omega\tau_\ell]$ fixes $(\tau_R,\tau_\ell)$ when $\delta\ne0$; the
  amplitude adds $b$.
- $\xi_1$ is invisible to all algebraic orders. This is note 7's Cattaneo half-space ($\tau_\ell=0$).

### (iv) $\omega\tau\gg1$

**GK, $\tau_\ell>0$ (ASYMPTOTIC, error $O(\omega^{-2})$ plus exponentially small terms).**
$$
Z=\frac{\sqrt{\tau_R/\tau_\ell}}{b\sqrt p}\Big[1+\frac{1}{2p}\Big(\frac1{\tau_R}-\frac1{\tau_\ell}\Big)+O(p^{-2})\Big]+O\big(e^{-\xi_\infty\sqrt{2\omega}}\big).
$$
The film is diffusive, not wave-like, with effective $b_\infty=b\sqrt{\tau_\ell/\tau_R}$,
$\xi_\infty=\xi_1\sqrt{\tau_R/\tau_\ell}$ and $a_\infty=a\tau_\ell/\tau_R$.

- The large-$p$ coefficients of $\ln g$ are $\tfrac12(\alpha-\beta')$ and $-\tfrac14(\alpha^2-\beta'^2)$, with
  $\alpha=1/\tau_R$, $\beta'=1/\tau_\ell$. With the prefactor, the Jacobian determinant with respect to $(b,\alpha,\beta')$
  is $(\alpha-\beta')/4b$ (CAS, `b3` G3). So (iv) fixes $(b,\tau_R,\tau_\ell)$ iff $\tau_R\ne\tau_\ell$.
- Reflections: $\Gamma\to\Gamma_\infty=(b_\infty-b_s)/(b_\infty+b_s)$ with $\Gamma-\Gamma_\infty=O(\omega^{-1})$, and
  $2\operatorname{Re}A=\xi_\infty\sqrt{2\omega}+\xi_\infty(1/\tau_R-1/\tau_\ell)/\sqrt{2\omega}+O(\omega^{-3/2})$.
  Numerical checks: remainder slope $-2.000$; the correction ratio to the prediction goes to 1 as $\omega^{-1}$; slope of
  $|\Gamma-\Gamma_\infty|$ is $-1.000$ (`b3` G1–G2).
- Hence $\xi_1$ is carried by $2\Gamma_\infty e^{-2A}$, of relative size
  $2|\Gamma_\infty|\exp(-\xi_1\sqrt{\tau_R/\tau_\ell}\sqrt{2\omega})$. It vanishes at leading order on the
  high-frequency matching surface $b\sqrt{\tau_\ell/\tau_R}=b_s$, which is a different surface from $b=b_s$.

**Cattaneo, $\tau_\ell=0$ (ASYMPTOTIC, error $O(\omega^{-1/2})$, not small in $\omega$ otherwise).**
$A=t_e p+t_e/2\tau_R-t_e/(8\tau_R^2p)+\dots$, with $t_e=\xi_1\sqrt{\tau_R}=e/v$ and $v=\sqrt{a/\tau_R}$.
$B\to b/\sqrt{\tau_R}=\rho c\,v$ (wave impedance). $\Gamma=-1+2B/(b_s\sqrt p)+O(1/p)$.
$$
Z(i\omega)=\frac{\sqrt{\tau_R}}{b}\tanh\Big(\frac{t_e}{2\tau_R}+i\omega t_e\Big)+O(\omega^{-1/2}).
$$
The film is wave-like: $\operatorname{Re}A\to t_e/2\tau_R=1/(2\sqrt{\varepsilon_R})$ is finite, so the echoes do not die
out with frequency.

- The response is almost periodic in $\omega$: period $\pi/t_e$, mean $\sqrt{\tau_R}/b$, first harmonic
  $2e^{-t_e/\tau_R}\sqrt{\tau_R}/b$.
- These three numbers give $(\xi_1,b,\tau_R)$ explicitly. A reconstruction from a window at $\omega_0$ converges as
  $\omega_0^{-1/2}$ (13 %, 4.6 %, 1.5 %, 0.5 % at $\omega_0\xi_1^2=10^2,\dots,10^5$; `b3` C1–C2).
- Echo amplitude per round trip, $e^{-1/\sqrt{\varepsilon_R}}$: 0.10 (30 nm), 0.016 (100 nm), $4\times10^{-4}$ (300 nm),
  $3\times10^{-5}$ (500 nm).
- In this regime $|\partial\ln Z/\partial\ln\xi_1|\sim4\omega t_e e^{-t_e/\tau_R}$ grows linearly in $\omega$.

*Non-uniformity between $\tau_\ell\to0$ and $\omega\to\infty$.* The limits do not commute. The distinguished limit is
$\omega\tau_\ell=O(1)$ with $\omega\tau_R\gg1$, where
$\operatorname{Re}A\simeq t_e/2\tau_R+t_e\omega(1+x^2)^{-1/4}\sin(\tfrac12\arctan x)$, $x=\omega\tau_\ell$. Wave echoes
occupy $1/\tau_R\ll\omega\ll1/\tau_\ell$ and are damped like $t_e\omega^2\tau_\ell/2$. Beyond that the response is
diffusive with $\xi_\infty$.

### Summary: regime map

| regime | condition | leading form | fixed at leading orders (relative order) | $\xi_1$ enters through |
|---|---|---|---|---|
| (i) | $\Omega\ll1$, $\omega\tau\ll1$ | $Z_s+R_fG-C_f/b_s^2+O(\sqrt p)$ | $R_f-C_f/b_s^2$ ($1$); $C_f(1-\beta^{-2})$ ($\sqrt\Omega$); $\delta$ ($\omega\tau$); $\tau_\ell\delta$ ($(\omega\tau)^2$) | $L_0,L_1$ (algebraic); at $b=b_s$ only via $L_2=\xi_1\delta/b_s$ |
| (i′) (O2 only) | $|A|\ll1$, $\omega\tau=O(1)$ | $Z_s+R_fG(p)-C_f/b_s^2$ | $R_f,\tau_R,\tau_\ell$ (shape of $G$), $C_f$: all four if $\delta\ne0$ | algebraic |
| (ii) | $1\ll\Omega\ll\varepsilon^{-1}$ | $g/(b\sqrt p)$ | $b$ ($1$); $\delta$ (phase, $\omega\tau$); $\tau_R+\tau_\ell$ (amplitude, $(\omega\tau)^2$) | $2\Gamma_0e^{-\xi_1\sqrt{2\omega}}$, exponentially small |
| (iii) | $\omega\tau\sim1$, $\Omega\gg1$ | $g/(b\sqrt p)$ | $b,\tau_R,\tau_\ell$ ($\delta\ne0$) | exponentially small ($\lesssim e^{-\sqrt{2/\varepsilon}}$) |
| (iv) GK | $\omega\tau\gg1$, $\tau_\ell>0$ | $\sqrt{\tau_R/\tau_\ell}/(b\sqrt p)$ | $b\sqrt{\tau_\ell/\tau_R}$ ($1$); $\tau_R^{-1}-\tau_\ell^{-1}$ ($1/\omega\tau$); $\tau_R^{-2}-\tau_\ell^{-2}$ ($1/(\omega\tau)^2$) | $2\Gamma_\infty e^{-\xi_\infty\sqrt{2\omega}}$ |
| (iv) Cattaneo | $\omega\tau_R\gg1$, $\tau_\ell=0$ | $(\sqrt{\tau_R}/b)\tanh(t_e/2\tau_R+i\omega t_e)$ | $b/\sqrt{\tau_R}$, $t_e$, $t_e/\tau_R$: all three | echo phase $2\omega t_e$, amplitude $e^{-t_e/\tau_R}$ |

## B4. Rank four from the expansions (independent route)

**THEOREM B4 (DERIVED; proof below, CAS coefficients; proof audit pending).** Let $b_s>0$ be known and
$\theta=(\xi_1,b,\tau_R,\tau_\ell)\in(0,\infty)^2\times[0,\infty)^2$, with one-sided derivatives on the boundary. The four
functions $\partial_{\theta_i}Z(i\omega;\theta)$, $\omega>0$, are linearly independent (over $\mathbb C$, hence over
$\mathbb R$) **iff $\tau_R\ne\tau_\ell$**. On $\{\tau_R=\tau_\ell,\ b\ne b_s\}$ the rank is 3. On
$\{\tau_R=\tau_\ell,\ b=b_s\}$ it is 2. The exceptional set is exactly $\Sigma=\{\tau_R=\tau_\ell\}$. **$b=b_s$ is not a
rank-deficient set when $\tau_R\ne\tau_\ell$, and neither are the boundaries $\tau_\ell=0$ (Cattaneo) or $\tau_R=0$.**

*Proof.*

1. **Reduction to Puiseux coefficients.** Suppose $\sum_ic_i\partial_iZ=0$ on the ray $q=\sqrt{i\omega}$. Since $Z_s$
   does not depend on $\theta$, this function equals $\sum_ic_i\partial_i(Z-Z_s)$, which is analytic in $q$ in a disc
   around 0 (B3). The ray contains a segment of that disc, so the function vanishes identically and all its Taylor
   coefficients vanish: $\sum_ic_i\partial_iL_k=0$ for every $k$. Differentiation commutes with the coefficients because
   the series converges with analytic dependence on $\theta$.
2. **Off the resonance.** CAS gives
   $\det\partial(L_0,L_1,L_2,L_4)/\partial(\xi_1,b,\tau_R,\tau_\ell)=-2\xi_1^4(b^2-b_s^2)(\tau_R-\tau_\ell)/(b^4b_s^3)$,
   which is the triangular product of the hand derivation. At $b=b_s$, where this vanishes, CAS gives
   $\det\partial(L_0,L_2,L_3,L_4)/\partial(b,\xi_1,\tau_R,\tau_\ell)=2\xi_1^4(\tau_R-\tau_\ell)^2/b_s^5$. One of the two
   determinants is nonzero whenever $\delta\ne0$, so $c=0$.
3. **On the resonance.** $G\equiv1$, so $(\partial_{\tau_R}+\partial_{\tau_\ell})Z\equiv0$ exactly. More generally
   $(\partial_{\tau_R}+\partial_{\tau_\ell})Z=-\partial_GZ\,p^2\delta/(1+\tau_\ell p)^2$. With $b=b_s$ as well,
   $\Gamma\equiv0$ and $Z\equiv Z_s$, so $\partial_{\xi_1}Z\equiv0$ too.
4. **Lower bounds on the resonance.** The rank of $\partial(L_0,\dots,L_5)/\partial\theta$ is 3 on $\{\delta=0\}$ and 2
   on $\{\delta=0,b=b_s\}$ (CAS), so the upper bounds of step 3 are attained. ∎

A second route that does not use the low frequencies needs $\tau_\ell>0$: (iv) gives $(b,\tau_R,\tau_\ell)$ when
$\delta\ne0$, then $L_0$ gives $\xi_1$ when $b\ne b_s$, or $L_2$ gives it ($\partial_{\xi_1}L_2=\delta/b$) when
$b=b_s$.

**Numerical confirmation (`b4`, independent of the CAS).** I took the SVD of the column-normalised real Jacobian
$\partial\ln Z/\partial(\ln\xi_1,\ln b,\tau_R,\tau_\ell)$ on 561 frequencies with $\Omega\in[10^{-6},10^8]$ (units
$\xi_1=b_s=1$), using analytic derivatives checked against finite differences to $10^{-9}$ and against 40-digit
evaluations.

| point | smallest normalised singular value |
|---|---|
| generic GK ($b=2.13$, $\tau_R=0.01$, $\tau_\ell=0.004$) | 0.112 |
| $b=b_s$, same times | 0.113 (unchanged) |
| $\tau_R=0$ | 0.37 |
| $\tau_\ell>\tau_R$ | 0.112 |
| resonance | $1.3\times10^{-16}$ (rank 3) |
| resonance with $b=b_s$ | with fixed scaling $\sigma=(1.41,0.43,1.3\times10^{-16},8.9\times10^{-20})$ (rank 2; $\xi_1$ column $\equiv0$) |

As $\delta\to0$, $\sigma_{\min}\propto|\delta|^{1.000}$ for $b\ne b_s$. Along $b=b_s$ two singular values vanish as
$|\delta|^{0.97}$ and $|\delta|^{1.00}$, and $\|\partial_{\ln\xi_1}\ln Z\|\propto|\delta|^{1.00}$. Near $b=b_s$ with
$\delta$ fixed, $\|\partial_{\ln\xi_1}\ln Z\|\simeq c_1|b-b_s|+c_2|\delta|/\xi_1^2$: weak, but not dependent.

The Cattaneo stratum needs care. Normalised $\sigma_{\min}$ falls from 0.53 to 0.10 to $1.1\times10^{-3}$ as the grid
extends to $\Omega=10^4,10^6,10^8$. The cause is the echo phase $2\omega t_e$: its derivative grows linearly and
dominates both the $\ln\xi_1$ and the $\tau_R$ columns ($\cos\to1.00000$), while the complementary combinations
$t_e/\tau_R$ and $b/\sqrt{\tau_R}$ carry bounded weight. This is a weighting effect; the rank is 4.

## B5. Global injectivity and explicit fibres (constructive)

**PROPOSITION B5 (DERIVED; CAS coefficients).** With $b_s$ known, $Z$ on any interval $\omega\in(0,\epsilon)$
determines $\theta$ uniquely on $\{\tau_R\ne\tau_\ell\}$, by the following explicit steps.

- **$L_0\ne0$ (that is, $b\ne b_s$):**
  - $C_f=-b_sL_1/L_0$ and $R_f=L_0-L_1/(b_sL_0)$, hence $\xi_1=\sqrt{R_fC_f}$ and $b=\sqrt{C_f/R_f}$;
  - $\delta=[L_2-\xi_1^3K(\beta)/b_s]/R_f$;
  - $\tau_\ell=-[L_4-F_4(\xi_1,b,\delta)]/(R_f\delta)$, where $F_4=L_4+R_f\tau_\ell\delta$ depends on the times only
    through $\delta$.
- **$L_0=0$ (that is, $b=b_s$):**
  - $\xi_1=-L_3/L_2$;
  - $\delta=b_sL_2/\xi_1$;
  - $\tau_\ell$ from $L_4$ as above.

Tests with 20 random $\theta$ each, at 50 digits: maximum relative error $1.2\times10^{-45}$ for $b\ne b_s$ and
$6.5\times10^{-48}$ for $b=b_s$ (`b2`, `b5`).

Since $\delta$ is a function of the data, the fibre of a resonant point is resonant:
- $\{(\xi_1,b,\tau,\tau):\tau\ge0\}$ if $b\ne b_s$;
- $\{(\xi_1',b_s,\tau,\tau):\xi_1'>0,\tau\ge0\}$ if $b=b_s$, in which case $Z\equiv Z_s$.

This answers B2 of the question, with $b_s$ known, at function level.

*Caveat (essential).* This identifiability holds in exact arithmetic. The fourth combination enters at relative order
$(\omega\tau)^2$ and the third at $\omega\tau$. In the FDTR band, $\omega\tau\le0.035$, so only $(\xi_1,b)$ and at best
$\delta$ are practically accessible (Part A, N3 and N4).

## B6. Where the asymptotic route fails, and how $\xi_1$ is carried at high frequency

- **$\tau_R=\tau_\ell$:** genuine. It is an exact invariance (Fourier resonance), the rank is 3, and the near-resonance
  conditioning is $\propto|\delta|$.
- **$b=b_s$:** not degenerate when $\delta\ne0$. The leading low-frequency $\xi_1$-information ($L_0$, $L_1$) vanishes,
  and $\xi_1$ is then carried by $L_2=\xi_1\delta/b_s$ and $L_3=-\xi_1^2\delta/b_s$. Equivalently,
  $\Gamma=(1-g)/(1+g)\simeq-\delta p/4$ at low frequency, so $|\partial\ln Z/\partial\ln\xi_1|\simeq\xi_1\delta\,\omega^{3/2}$
  (observed slope 1.44, tending to 3/2). At high frequency, $\Gamma\to\Gamma_\infty=(1-\sqrt{\tau_R/\tau_\ell})/(1+\sqrt{\tau_R/\tau_\ell})\ne0$.
- **Exponentially small reflections in (ii)–(iv):** $|\partial\ln Z/\partial\ln\xi_1|=4|A\Gamma E|/|1-\Gamma^2E^2|$ exactly.
  - (ii): $\simeq4|\Gamma_0|\xi_1\sqrt\omega\,e^{-\xi_1\sqrt{2\omega}}$, using the corrected exponent of B3(ii).
  - (iv) GK: $\simeq4|\Gamma_\infty|\xi_\infty\sqrt\omega\,e^{-\xi_\infty\sqrt{2\omega}}$, with
    $|\Gamma-\Gamma_\infty|=O(\omega^{-1})$.
  - (iv) Cattaneo: $\simeq4\omega t_e e^{-t_e/\tau_R}/|1-e^{-2t_e/\tau_R}e^{-4i\omega t_e}|$, which grows like $\omega$.
  - For 500 nm AlN at 200 MHz, $2\operatorname{Re}A\approx2.75$, so the reflection term has relative weight
    $2\Gamma_0e^{-2.75}\approx0.05$. It is not negligible there: the band top is in the (i)/(ii) crossover, not in (ii).
- **Boundary $\tau_\ell=0$ with $\omega\to\infty$:** non-uniform (B3, iv).

---

# Part A — scaling of the FDTR phase sensitivity to $\tau_R$

1D (large-beam) limit of the stack Au (80 nm, 150 W/mK, $2.49\times10^6$) / $R_1$ / AlN (500 nm, $\lambda_z$, $2.41\times10^6$,
Cattaneo $\tau_R$) / $R_2$ / sapphire (35, $3.03\times10^6$). Base values: $\lambda_z=200$, $\tau_R=27.5$ ps,
$R_1=10^{-8}$, $R_2=2\times10^{-8}$ m$^2$K/W. Implementation checks:
- The 1D stack equals the twin kernel (`src/twin_fdtr.kernel` at $k=0$, a modal linear system) to $10^{-15}$ for
  Cattaneo, GK (local temperature) and Fourier, and $\partial\varphi/\partial\ln\tau_R$ agrees to $2\times10^{-11}$ deg
  (`x1`).
- Derivatives use forward-mode dual numbers, exact to rounding. Finite differences of $\ln Z$ are pure roundoff at low
  frequency, where $|\partial\ln Z/\partial\ln\tau_R|\sim10^{-9}$.

## A1. Scales

| symbol | definition | base value |
|---|---|---|
| $R_f$ | $e/\lambda_z$ | $2.5\times10^{-9}$ m$^2$K/W ($R_\Sigma=R_1+R_f+R_2=3.25\times10^{-8}$) |
| $C_f$, $C_A$ | $\rho c\,e$ of AlN, Au | 1.205, 0.199 J m$^{-2}$K$^{-1}$ |
| $f_*$ | $1/[2\pi(b_sR_\Sigma)^2]$: substrate ↔ resistances | 1.4 MHz |
| $f_s$, $f_2$ | $1/[2\pi(b_sR_2)^2]$; $1/(2\pi C_fR_2)$: film capacitance shunts $R_2$ | 3.75, 6.6 MHz |
| $f_a$ | $1/(2\pi C_fR_1)$ | 13.2 MHz |
| $f_{\rm film}$ | $a/(2\pi e^2)$: AlN thermally thick above | 52.8 MHz |
| $f_b$ | $(1+C_A/C_f)/(2\pi C_AR_1)$: transducer capacitance shunts $R_1$ | 93.1 MHz |
| $\omega\tau_R$ at 200 MHz | | 0.035 |

The scales are not separated (3.75–93 MHz). The closed forms below are therefore sub-regime asymptotes and an
order-one-accurate peak law, not a uniform expansion.

## A2. Exact identities

**(EXACT)** Since $\lambda\to\lambda_{\rm eff}$,
$$
\frac{\partial\ln Z}{\partial\ln\tau_R}=-\frac{\tau_Rp}{1+\tau_Rp}\;s_\lambda\big(p;\lambda_{\rm eff}(p)\big),\qquad
s_\lambda=\frac{\partial\ln Z}{\partial\ln\lambda}\Big|_{C,e},\qquad
\partial_{\tau_\ell}Z=-G(p)\,\partial_{\tau_R}Z .
$$
Checked to $10^{-13}$ (`a1`, I1).

**(EXACT, reciprocity of the 1D two-port)** Let $T(z)=q(z)/q(0)$ for unit front flux, and
$\langle T^2\rangle=e^{-1}\int_{\rm film}T^2dz$. Then
$$
\partial_{\tau_R}Z=\frac{pR_f}{1+\tau_\ell p}\langle T^2\rangle,\qquad\partial_{R_i}Z=T_i^2,\qquad\partial_{\ln\lambda}Z=-R_f\langle T^2\rangle\ (\tau=0),
$$
where $T_1$ and $T_2$ are taken at the AlN top and bottom. Checked to $10^{-13}$ and $10^{-10}$ (I2).

## A3. Leading order: phase change ∝ $\omega\tau_R w$ times a film weight

**(ASYMPTOTIC; error $O(w^2)$ and $O((\omega\tau)^2)$ absolute)**
$$
\delta\varphi(\omega)=w\,\omega\tau_R\,W(\omega),\qquad W(\omega)=-\frac{\partial\ln|Z|}{\partial\ln\lambda_z}=R_f\operatorname{Re}\frac{\langle T^2\rangle}{Z},
\qquad \delta\ln|Z|=w\,\omega\tau_R\,\frac{\partial\varphi}{\partial\ln\lambda_z}.
$$
The phase sensitivity to $\tau_R$ is $\omega\tau_R$ times the amplitude sensitivity to the film conductivity, and the
roles of phase and amplitude exchange; the time acts as an imaginary conductivity perturbation $-i\omega\tau_R$.

Accuracy:
- At $w=0.08$ the deviation from linearity in $w$ is $6.5\times10^{-3}$ relative at 200 MHz.
- $|\text{first order}-\text{leading order}|/\max|\cdot|$ scales as $(\omega\tau)^{2.25}$ over the fitted range, and at
  fixed frequency $(\cdot)/\tau^2\to$ const as $\tau\to0$, which is the predicted $O(\tau^2)$. At 200 MHz the error is
  0.027 of the maximum.

$W$ for the base stack: 0.014 (100 kHz), 0.033 (1 MHz), 0.048 (10 MHz), 0.051 (30 MHz), 0.018 (100 MHz), and $-0.0016$
(200 MHz). $\delta\varphi$ changes sign at 178 MHz. The amplitude signature is at most $7\times10^{-4}w$, negligible
against 1 % noise.

**Interface form (ASYMPTOTIC, thin AlN, error $O(|A_f|^2)$, observed slope 1.000).** With
$\langle T^2\rangle=(T_1^2+T_1T_2+T_2^2)/3$,
$$
\frac{\partial\varphi}{\partial\tau_R}\simeq\frac{\omega R_f}{3}\Big[\frac{\partial\ln|Z|}{\partial R_1}+\frac{\partial\ln|Z|}{\partial R_2}+\operatorname{Re}\frac{T_1T_2}{Z}\Big].
$$
So the $\tau_R$ phase signature is $\omega$ times the *amplitude* sensitivities to the interface resistances, not their
phase sensitivities (whitened cosines: $-0.42$ with $\partial\varphi/\partial R_1$, $-0.18$ with
$\partial\varphi/\partial R_2$, $+0.91$ with $\omega\,\partial\ln|Z|/\partial R_1$). The approximation error is 3 % at
10 MHz, 25 % at 100 MHz and 51 % at 200 MHz.

## A4. Closed forms

1. **Substrate-dominated, $f\ll f_*$ (ASYMPTOTIC, relative error $O(b_sR_\Sigma\sqrt\omega)$):**
   $\delta\varphi=w\,\tau_RR_fb_s\,\omega^{3/2}/\sqrt2$. Exact/closed = 0.999 at 1 kHz and 0.90 at 300 kHz. The signal
   is negligible there ($<10^{-5}$ deg for $w=0.05$).
2. **Transducer-limited band (ASYMPTOTIC for a thin AlN film fully shunted by $C_f$, with $R_f\ll R_1$):**
   $\delta\varphi\simeq K\,\dfrac{\omega^3(\omega_a+\omega_b)}{(\omega^2+\omega_a^2)(\omega^2+\omega_b^2)}$, with
   $K=\dfrac{w\tau_RR_f}{3C_AR_1^2}$.
   - For $\omega_a\ll\omega\ll\omega_b$: $\delta\varphi\simeq\omega\,w\tau_RR_f/[3R_1(1+C_A/C_f)]$, which is
     **delay-like**.
   - For $\omega\gg\omega_b$: $\delta\ln Z\to K$ (real), which is **gain-like**, and $\delta\varphi\to0$.
3. **Peak law (HEURISTIC prefactor):**
   $$
   |\delta\varphi|_{\max}=\kappa\,\frac{w\,\tau_R\,R_f}{6\,C_A\,R_1^2},\qquad \kappa=0.50\text{–}0.79,\quad f_{\rm peak}\approx(0.54\text{–}1.05)\,f_b .
   $$
   $\kappa=1$ is the thin-film limit. $\kappa$ was fitted on 14 scenarios ($\lambda_z$, $\tau_R$, $R_1\in[5\times10^{-9},2\times10^{-8}]$,
   $R_2$, Au 40–120 nm, 100/300/500 nm films). It is below 1 because the AlN is in its own crossover at the peak
   ($\Omega\approx1.2$). Practical form:
   $|\delta\varphi|_{\max}\approx0.024^\circ\,w\,(\tau_R/30\ \text{ps})(e/500\ \text{nm})(200/\lambda_z)(10^{-8}/R_1)^2(0.2/C_A)$,
   accurate to ±20 % while $f_b$ lies in the band. $R_2$ barely matters ($\kappa$ 0.63–0.71 for $R_2=1$–$5\times10^{-8}$),
   because $C_f$ shunts it near the peak.
4. **Lumped-transducer form, valid through the AlN crossover (Au thermally thin, $\omega e_{\rm Au}^2/a_{\rm Au}\le0.13$):**
   $\delta\varphi/w=\omega\tau_R\operatorname{Re}\{\zeta_f/[X(1+i\omega C_AX)]\}$, with $X=R_1+Z_f$ and
   $\zeta_f=-\partial Z_f/\partial\ln\lambda$. It gives the peak to 10 % and the curve to 21 % of the maximum (`a3`).

## A5. Size of the effect versus phase noise (exact 1D)

Maximum over 10 kHz–200 MHz of $|\partial\varphi/\partial\ln\tau_R|$, in degrees per unit $w$ (`a1`, S1):

| scenario | max deg/w | $f$ at max |
|---|---|---|
| 500 nm, $\lambda_z=150$; $\tau_R=20/27.5/35$ ps | 0.019 / 0.027 / 0.034 | 59 MHz |
| 500 nm, $\lambda_z=200$; same | 0.016 / 0.022 / 0.028 | 63 MHz |
| 500 nm, $\lambda_z=250$; same | 0.014 / 0.019 / 0.024 | 66 MHz |
| 300 nm, $\lambda_z=180$ (150), $\tau_R=20$ ps | 0.012 (0.014) | 83 MHz |
| 100 nm, $\lambda_z=130$ (100), $\tau_R=11$ ps | 0.0032 (0.0041) | 151 MHz |
| 500 nm base with $R_1=R_2=2\times10^{-8}$ | 0.0063 | 36 MHz |
| 500 nm base with $R_1=10^{-9}$ | 0.36 | 200 MHz (band edge) |
| 500 nm base with $R_1=R_2=0$ | 0.62 | 200 MHz |
| bare AlN on sapphire (no Au, no $R$) | 1.13 | 200 MHz (half-space $\omega\tau/2$) |

For $w=0.03$–$0.08$ at 500 nm (base), $|\Delta\varphi|_{\max}=6.6\times10^{-4}$°–$1.8\times10^{-3}$°, and
$4.1\times10^{-4}$°–$2.7\times10^{-3}$° over the $\lambda_z$ and $\tau_R$ ranges. That is **4–27 % of a 0.01° noise
level and 0.4–2.7 % of 0.1°**. At 100 nm, $(1$–$3)\times10^{-4}$°.

The dominant suppressor is the Au/AlN contact. Relative to bare AlN, the Au layer alone costs a factor of 2, and
$R_1=10^{-8}$ costs a further 28: the $C_AR_1$ low-pass sits inside the band.

## A6. Structural question: is $\partial\varphi/\partial\tau_R$ collinear with a delay?

**Leading form.** $\partial\varphi/\partial\ln\tau_R=\omega\tau_RW(\omega)$, against $\partial\varphi/\partial t_d=-\omega$.

**PROPOSITION A6 (EXACT identity at leading order; phase data, equal noise).** Let
$\mu_k=\omega_k^2/\sum_j\omega_j^2$ on the sampled frequencies. Then
$$
\cos^2\angle(\omega W,\omega)=\frac{E_\mu(W)^2}{E_\mu(W^2)},\qquad
1-\cos^2=\frac{\operatorname{Var}_\mu W}{E_\mu W^2}.
$$
This is the fraction of $\tau_R$-information that survives a free delay. With the offset free as well, it becomes
$\min_{a,c}E_\mu[(W-a-c/\omega)^2]/E_\mu W^2$.

Values for the 80 log-spaced frequencies (leading order / exact signature; `a3`):

| stack | surviving a free delay | surviving delay and offset |
|---|---|---|
| 500 nm base | 0.657 / 0.667 ($\cos=0.58$) | 0.565 / 0.574 |
| 300 nm | 0.34 | 0.30 |
| 100 nm | 0.088 | 0.073 |
| 500 nm, $R_1=R_2=0$ | 0.0072 | 0.0056 |
| bare AlN | 0.024 | 0.020 |

**Verdict on the delay.** The $\tau_R$ signature is a pure delay exactly when the film weight $W$ is flat under the
$\omega^2$-weighted measure, as in the Cattaneo half-space: $\varphi=-\pi/4+\tfrac12\arctan\omega\tau$, so
$\delta\varphi=\omega\tau w/2$, an advance $t_{\rm eq}=-\tau_Rw/2$. It is also nearly a delay in the resistive plateau
$f_a\ll f\ll f_b$.

In the declared Au/AlN/sapphire stack with $R_1=10^{-8}$, the transducer corner $f_b=93$ MHz, the film crossover
$f_{\rm film}=53$ MHz and the sign change of $W$ at about 180 MHz all lie inside the band. The signature is therefore
**not** nearly collinear with the delay: the delay removes 33 % of the information ($w^*\times1.22$), and N1 removes
43 % ($\times1.32$). The equivalent delay of the projected part is $-0.27$ ps per unit $w$, that is
$-0.014$ ps at $w=0.05$, three orders of magnitude below any realistic delay calibration (note 16 prior: 20 ps).

There is a trade-off. Interfaces that de-correlate the signature from the delay also suppress it by about 30×.
Removing them ($R_1=R_2=0$) restores the size but makes it 99.3 % delay-collinear. Even then $w^*(N1)=0.094$ at 0.01°
(A7).

**Stronger degeneracy: the thermal inductance (ASYMPTOTIC, resistive regime).** At leading order
$\delta Z_\tau=i\omega\,w\tau_RR_f\langle T^2\rangle$, a distributed series "thermal inductance" $L_\tau=\tau_RR_f$.
Where the circuit around the film is resistive ($T\simeq$ const and $\theta\simeq R_{\rm below}$), every heat-capacity
perturbation and the delay also give $\delta Z\propto i\omega\times$(real constant):
- $\delta Z_{C_j}=-i\omega\,\delta C_j\langle Z_{\rm below}^2\rangle$;
- $\delta Z_{t_d}=-i\omega t_dZ$.

At this order only the net coefficient of $i\omega$ is identifiable, $\tau_RR_f\langle T^2\rangle-\sum_jC_j\langle Z_{{\rm below},j}^2\rangle-t_dR_\Sigma$.
The degeneracy is broken at relative order $\omega/\omega_{\rm corner}$ by the different spatial weights ($T^2$ against
$\theta^2$).

Whitened single-column cosines with the $\tau_R$ column (amplitude and phase, $\sigma_\varphi=0.01°$):

| nuisance | $\tau_\ell$ | $C_{\rm Au}$ | $\lambda_{\rm AlN}$ | $\lambda_{\rm Au}$ | $C_{\rm AlN}$ | delay | $\varphi_0$ | $R_1$ | $R_2$ | gain |
|---|---|---|---|---|---|---|---|---|---|---|
| cosine | $-0.9994$ | $-0.92$ | 0.71 | $-0.67$ | 0.60 | $-0.58$ | 0.53 | $-0.42$ | $-0.18$ | 0.01 |

The $\tau_\ell$ entry follows from the exact relation $\partial_{\tau_\ell}Z=-G\,\partial_{\tau_R}Z$: the angle is of
order $\omega\tau$, consistent with B3′.

Information retained with the $\tau_R$ column for nested sets: N1 0.57; N0+$R_{1,2}$ 0.69; N2 0.22; N0+$C_{\rm Au}$
0.15; N1+$C_{\rm Au}$ 0.034; N1+$\lambda$ 0.048; N2+$C$ 0.075; N2+$\lambda$ 0.038; N3 $5.9\times10^{-3}$; N0+$\tau_\ell$
$1.2\times10^{-3}$; N4 $8.7\times10^{-9}$. The normalised Fisher condition number reaches $1.5\times10^9$ for N4; an
independent QR projection reproduces every retention to 5 digits.

## A7. Consequence for detectability (1D Fisher, scaling level)

$w^*(1\sigma)$ is the relative change of $\tau_R$ resolved at 1σ, with 80 frequencies from 10 kHz to 200 MHz and 1 %
amplitude noise. $2\sigma$ is twice these values. N5 (beam radii) has no 1D analogue; $\lambda_{\rm Au}$ and
$C_{\rm Au}$ are included in N4.

| scenario | $\sigma_\varphi$ | N0 | N1 | N2 | N3 | N0+$\tau_\ell$ | N4 |
|---|---|---|---|---|---|---|---|
| 500 nm, $\lambda_z=200$, 27.5 ps | 0.01° | 0.14 | 0.18 | 0.29 | 1.8 | 4.0 | $1.5\times10^3$ |
|  | 0.1° | 1.3 | 1.7 | 2.4 | 13 | 39 | $1.3\times10^4$ |
| 500 nm, $\lambda_z=150$, 35 ps (most favourable) | 0.01° | 0.090 | 0.11 | 0.17 | 1.3 | 1.9 | 970 |
|  | 0.1° | 0.87 | 1.05 | 1.4 | 9.7 | 19 | $8.3\times10^3$ |
| 500 nm, $\lambda_z=250$, 20 ps | 0.01° | 0.22 | 0.31 | 0.50 | 2.7 | 9.2 | $3.0\times10^3$ |
| 300 nm, $\lambda_z=180$, 20 ps | 0.01° | 0.24 | 0.44 | 0.62 | 2.4 | 13 | $4.0\times10^3$ |
| 100 nm, $\lambda_z=130$, 11 ps | 0.01° | 0.99 | 3.7 | 4.8 | 12 | 167 | $4\times10^5$ |
| 500 nm, $R_1=R_2=0$ | 0.01° | 0.0072 | 0.094 | 0.12 | 12 | 1.3 | $2.7\times10^3$ |

Phase-only results are identical to amplitude-plus-phase for N0–N3 to the digits shown.

Against the admissible half-widths of 0.03–0.05 (orb 0.1, 30–300 nm) and up to 0.08 (calibrated classes, 100 nm),
$w^*$ is above them in every declared case **including N0 at 0.01°**. The margins $w^*/w$ are:

| nuisance set, noise | margin |
|---|---|
| N0, 0.01° | 1.1–7 |
| N1, 0.01° | $\ge1.4$ |
| N2, 0.01° | $\ge2.2$ |
| N ⊇ N2, 0.1° (the falsification point of $H_A$) | $\ge17$ at 500 nm, $\ge60$ at 300 nm |

The weakest case is the most favourable 500 nm scenario ($\lambda_z=150$, $\tau_R=35$ ps) at N0: $w^*=0.090$ against
the largest half-width 0.08, marginal at 1σ and not resolved at 2σ (0.18). The only clearly resolvable configuration in
the table is the counterfactual perfect-contact stack at N0 ($w^*=0.0072$). It is no longer resolved at N1
($w^*=0.094$).

**Conclusion of the scaling analysis (1D):** consistent with $H_A$, which is not falsified at the scaling level. The
three reasons are:
1. The interface $R_1$ caps the signal at about $0.02°\,w$.
2. The signal is degenerate at leading order with the transducer and film heat capacities and, partly, with the delay
   (the thermal-inductance structure).
3. $\tau_R$ is inseparable from $\tau_\ell$ at first order (only $\tau_R-\tau_\ell$ is accessible).

## A8. Robustness and limits of the 1D estimate

The finite-beam twin (`gaussian_response`, $w_p=10$ µm, $w_r=8$ µm, isotropic Cattaneo AlN) agrees with 1D where the
information lies: 3D/1D = 0.94–0.99 from 40 to 200 MHz and 0.82–0.88 from 10 to 30 MHz (`x2`). The 3D peak is
$0.021°/w$, against $0.022°/w$ in 1D.

Below 5 MHz the 3D signature exceeds 1D by up to 15× at 100 kHz. The twin's Cattaneo law relaxes the in-plane
conductivity too, and radial spreading in the AlN then matters. The absolute size is negligible ($\le2\times10^{-4}°/w$
at 100 kHz, $\le10^{-3}°/w$ at 1 MHz). Its shape is not delay-like, but in 3D it is confounded with in-plane
$\lambda_r$ and the beam radii, which the full-twin Fisher analysis must settle. Two further caveats:
- the in-plane memory time of the film is not $\tau_{\rm mem,z}$ (assumption A2);
- the 1D nuisance correlations at low frequency are not those of 3D, so the 3D $w^*$ for N2–N5 should be taken from the
  full analysis. The margins above make the scaling verdict insensitive to factors of 2.

---

## Weakest steps and open points

1. **Part A is 1D.** The verdict rests on exact 1D numerics, cross-checked against the twin kernel, with large margins.
   It does not include finite-beam nuisances (N5) or the low-frequency radial component of an isotropic Cattaneo law.
   The admissible half-width at 500 nm is extrapolated, not computed.
2. **The peak law has a fitted prefactor.** In $|\delta\varphi|_{\max}=\kappa\,w\tau_RR_f/(6C_AR_1^2)$, $\kappa=0.5$–0.8
   is fitted, not derived. The thin-film closed form ($\kappa=1$) fails at the peak because the AlN is in its own
   crossover there. This is the weakest analytic step; the verdict uses the exact 1D values, not this law.
3. **Part B is function-level.** Theorem B4 and Proposition B5 rest on CAS-computed coefficients up to $L_5$
   (independently confirmed by 50-digit remainder slopes and exact inversions) and on a standard analyticity argument;
   they still need the proof audit. They establish identifiability in exact arithmetic from any low-frequency interval.
   The practical information about the fourth combination is $O((\omega\tau)^2)$ in relative terms, which is far below
   any FDTR noise.

## Claims for the ledger

| ID | claim | label | status |
|---|---|---|---|
| AB-1 | Reflection representation (B2) and $\partial\ln Z/\partial\ln\xi_1=-4A\Gamma E/(1-\Gamma^2E^2)$ | EXACT | PROVED (elementary) |
| AB-2 | Low-frequency coefficients $L_0$–$L_5$ and their time structure | EXACT | DERIVED (CAS + 50-digit check) |
| AB-3 | Rank 4 iff $\tau_R\ne\tau_\ell$; rank 3 / 2 on $\Sigma$ without / with $b=b_s$; $b=b_s$ not degenerate | THEOREM | DERIVED UNDER ASSUMPTIONS (B1–B3), proof audit pending; NUMERICALLY SUPPORTED |
| AB-4 | Injectivity off $\Sigma$, explicit fibres on $\Sigma$, explicit inversion | PROPOSITION | DERIVED UNDER ASSUMPTIONS; NUMERICALLY SUPPORTED |
| AB-5 | Regime map (i)–(iv), Cattaneo almost-periodic high-frequency form, GK effective half-space, non-uniformities | ASYMPTOTIC | NUMERICALLY SUPPORTED (slopes) |
| AB-6 | Times enter only via $\delta$ at first order, uniformly in $\Omega$ | ASYMPTOTIC | DERIVED |
| AB-7 | $\partial\ln Z/\partial\ln\tau_R=-[\tau_Rp/(1+\tau_Rp)]\,s_\lambda(\lambda_{\rm eff})$; reciprocity forms | EXACT | PROVED (1D), checked to $10^{-13}$ |
| AB-8 | $\delta\varphi=w\omega\tau_RW$, $W=-\partial\ln|Z|/\partial\ln\lambda$; closed forms A4 | ASYMPTOTIC / HEURISTIC (peak $\kappa$) | NUMERICALLY SUPPORTED |
| AB-9 | Delay collinearity $=1-\operatorname{Var}_\mu W/E_\mu W^2$; not collinear in the declared stack (cos 0.58) | EXACT identity / numerical | NUMERICALLY SUPPORTED (1D) |
| AB-10 | Thermal-inductance degeneracy with $C_{\rm Au}$, $C_{\rm AlN}$, $\lambda$, delay; $\tau_\ell$ cosine $-0.9994$ | ASYMPTOTIC / numerical | NUMERICALLY SUPPORTED (1D) |
| AB-11 | 1D $w^*$ exceeds the admissible half-widths for all declared N, including N0 at 0.01° (marginal in the most favourable N0 case: 0.090 against 0.08); margin $\ge17$ at N ⊇ N2, 0.1° | numerical | NUMERICALLY SUPPORTED (1D); 3D confirmation pending |

## Files

All scripts are in `derivations/AB-asymptotic/` and use the research Python with 2 threads. `run_all.py` reruns
everything and writes `results/run_all.log`.

| script | content |
|---|---|
| `ab_common.py` | Film and stack models; overflow-safe tanh/sech; analytic Jacobian of the film; forward-mode dual numbers for the stack |
| `b2_lowfreq_series.py` | $L_0$–$L_5$ (SymPy), structure checks, coefficient determinants and ranks, 50-digit remainder slopes, explicit inversion ($b\ne b_s$) |
| `b3_highfreq_reflections.py` | Lemma B2 sampling; GK and Cattaneo high-frequency forms; reflection asymptotics; Cattaneo reconstruction; film table |
| `b4_rank_strata.py` | SVD on the strata; $\delta$- and $b$-scans; fixed-scaling column norms; 40-digit cross-check |
| `b5_xi1_sensitivity_and_bs_inversion.py` | Exact $\xi_1$-sensitivity; regime-(ii) non-uniformity; matched-stratum inversion |
| `a1_phase_sensitivity.py` | Exact identities I1–I3; leading order and its errors; closed forms; scenario table |
| `a2_nuisance_collinearity.py` | Delay/offset geometry; whitened cosines; 1D Fisher $w^*$ and retentions (N0–N4 and partial sets) |
| `a3_delay_variance_and_lumped.py` | Variance identity (A6); lumped-transducer form |
| `x1_twin_crosscheck.py` | 1D stack against `twin_fdtr.kernel(k=0)` |
| `x2_beam_vs_1d.py` | Finite-beam twin against the 1D $\tau_R$ signature |
