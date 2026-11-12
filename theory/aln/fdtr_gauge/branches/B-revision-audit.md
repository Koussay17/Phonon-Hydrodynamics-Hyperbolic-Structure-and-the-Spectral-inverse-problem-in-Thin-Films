# B-revision-audit — line-by-line audit of `proofs/B-revision.md`

Campaign `20261111-fdtr-gauge-memory`, revision branch R-B, 12 November 2026. Audited object: `proofs/B-revision.md`
(12 Nov, 02:28; objections B-2 and B-4 of `12-peer-review.md`) with its scripts and outputs in
`derivations/B-revision/`, including the elimination record `b2_03_reduce_record.md`. Lemma LW and Lemma T are taken as
audited in `branches/B-proof-audit.md` and are not re-audited. Model code: `C:\Users\Koussay\these\src\twin_fdtr.py`
(last commit touching the file: `663771d`; no local modification). Independent recomputations:
`derivations/B-revision-audit/` (§1). Each step is classified VALID / VALID UNDER EXTRA ASSUMPTION / INCOMPLETE /
UNJUSTIFIED / FALSE / UNKNOWN; defects D-1 … D-11 are collected in §3 with their minimal repairs.

---

## 0. Verdicts

| Item (priority order of the brief) | Verdict | Defects |
|---|---|---|
| **Lemma K** (kernel reduction; real-axis path, and the path into $s_0=i\tau_\ell^{-1/2}$ when $\tau_R=0<\tau_\ell$) | **VALID** (elementary; replaces Lemma LW in every rank statement) | D-7 (presentation) |
| **Theorem C2** and table §0.2 (read-out $T(0)$; local, tied $\varepsilon\ne0$, free $c_n$) | **PROVED** (hand analysis complete and re-derived; every entry confirmed) | — |
| **Proposition C4** ($\vartheta_e$ read-out, bottom coefficient $\varepsilon\tau_\ell$) | degenerate **sets** PROVED; exotic strata $E_0$, $E_\varepsilon$ PROVED; **ERROR** in the kernel/rank at $\varepsilon=1$ on $\Pi\cap\{b=b_s,R=0\}$ | **D-1** |
| **Proposition C5** (film block with $\vartheta$ at both faces; tied and free $\tau_n$) | **PROVED** (every coincidence stratum re-derived by hand) | — |
| Different coefficients at the two faces | DERIVED, as labelled; consistent with an exact special-value scan | — |
| **Theorem T1** (Möbius reduction preserves rank) | **VALID**; per-frequency form needs the stack to be pole-free on $i\mathbb R_+$ | D-4 |
| **Theorem T2** (global, via Lemma G) | **PROVED** given Lemma LW | D-6 (sign typo) |
| **Lemma G** | **VALID** given Lemma LW | — |
| **Theorem T3** and table §0.3 (local, every tied $\varepsilon$, free $c_n$) | **PROVED** (both cases of the stack part of Lemma K re-derived; 11 stack points at 60 digits, a03 and a10) | — |
| **Theorem T4** (i) | VALID | — |
| T4 (ii) | VALID for the certified connected components | D-4, D-5 |
| T4 (iii) | VALID | — |
| T4 (iv) ($\sigma^{15}t(\lambda,\mu)G'_{12345}$; excluded sets $\{P=0\}$, $\{Q=0\}$) | **PROVED**; $t$, $P$, $Q$ re-derived by an independent (Riccati) route, exact agreement | — |
| **Proposition TW3** (three-frequency twins, local convention) | **PROVED**; all three twins re-certified independently | D-11 (presentation) |
| **Proposition O** (twin reads $T(0)$; convention at both faces of a GK layer) | **PROVED** (code read line by line) and independently confirmed to $1.3\cdot10^{-15}$ | — |
| Lemma R, Lemma 3′ | VALID | D-8 (wording) |
| Remark after Lemma 3′ (poles on $i\mathbb R_+$) | VALID for $\varepsilon\in\{0,1\}$; the sentence on other $\varepsilon$ is **FALSE** for $0<\varepsilon<1$ (conservatively so) | D-2 |
| Piecewise exact eliminations for general $\varepsilon$ (§3.3, §4, weakest step 3) | **VALID** for the degenerate sets $\{\operatorname{rank}<n\}$ (decomposition exhaustive); **INCOMPLETE** for kernel dimensions on degenerate strata — the gap that contains D-1 | D-9, D-10 |
| Proposition LF, $\delta_\varepsilon$ | VALID (independent route) | — |
| §7, answer B-2 (iv) ("the entropic one is the one with a sign-definite 1D dissipation identity") | **FALSE** as stated (interpretive) | D-3 |

**Earliest failure.** In logical order the first false sentence is the last clause of the Remark after Lemma 3′ (D-2),
but it is a conservative non-claim with no downstream use. The earliest failure that changes a stated result is in the
proof of Proposition C4: the step "The case $R=0$ is the same computation with $c(1-\varepsilon)=0$" yields $c=0$ only for
$\varepsilon\ne1$; the case $\varepsilon=1$ is treated for $t>0$ only, and at $t=0$, $b=b_s$ the remaining constraint ($p^1$
coefficient of (E), $\xi_1c\,(b^2-b_s^2)=0$) is vacuous as well, so a second kernel direction survives (D-1). Nothing in the transducer results (T1–T4, TW3) or in the
twin's own observable ($T(0)$, Theorem C2) depends on it.

---

## 1. Independent recomputations (`derivations/B-revision-audit/`)

All scripts single-threaded, < 0.2 GB. None imports `kernel_system.py`, `models_mp.py` or `ad_arb.py` (a10 reuses the
audit's own a03 model).

| Script / output | Method (independent of the branch by …) | Result |
|---|---|---|
| `a01_prop_O_twin.py/.json` | twin `kernel(f, 0, stack)` vs the audit's **reflection-form** closed forms; $\alpha\in\{1/3,2\}$; $(\tau_R,\tau_\ell)\in\{(25,10),(25,25),(0,20),(30,0)\}$ ps; local, entropic, kinetic $c_n=+0.4,-0.7$ (units $\ell^2/\lambda$); $R\in\{0,3\cdot10^{-8}\}$; 32 frequencies 10 kHz–100 GHz | bare $T(0)$: $1.34\cdot10^{-15}$; `theta_e` $=(1+\tau_\ell p)\times$ it: $1.33\cdot10^{-15}$; Au/AlN/sapphire stack with $N_t=N_b$: $1.12\cdot10^{-15}$. Excluded: $\vartheta$ read at the heated face ($\ge0.64$), convention at bottom face only ($\ge7.5\cdot10^{-5}$), at top face only ($\ge1.5\cdot10^{-5}$). Slip coefficient $C_s\in\{0,1,\infty\}$ changes $T(0)$ at $k=0$ by $\le1.3\cdot10^{-16}$ (S family decoupled). Defaults `observable="T"`; `synthetic_measurement` passes none; Cattaneo layers get $(c_n,c_t)=(0,0)$ in all three conventions |
| `a02_lemmaK_symbolic.py/.out` | SymPy, own chain-rule derivation on base quantities | (K.1); invariance of $\Phi$ and of (K.1) (factor $r^2E^2$) under $(A,B,r,E)\to(-A,-B,1/r,1/E)$; stack form and its $E^0,E^1,E^2$ coefficients; $D\tilde B\equiv0$ under (i); combined identity $\Leftrightarrow D_cu=(u^2-1)D_cA$; odd part $=(s/b_s)$(O), even part $-Rs^2$(O) $=$ (E); $\rho A^2\to\xi_1^2e^{i\pi/4}/(2\tau_\ell^{3/2})$ on the special path; $u\equiv1$ impossible on $\{\tau_R=0<\tau_\ell\}$ — all differences $0$ |
| `a03_rank_numeric.py/.json` | 60 digits, single-valued $(C,X,Y)$ form (B-analytic Lemma 2), 25 frequencies $10^{-2}$–$10^2$, central differences $h=10^{-22}$ | 17 points (C2, C4, C5, T3): ranks and kernels as derived by hand in §2; deficient singular values $\le1.5\cdot10^{-38}$, kernels annihilated to $\le2.2\cdot10^{-39}$. **C4 at $\varepsilon=1$, $\Pi$, $b=b_s$, $R=0$: rank 2/4 ($R$ known) and 3/5 ($R$ unknown); $e_\xi$ and $(0,0,-1,1[,0])$ annihilated to $8\cdot10^{-40}$**; control $\varepsilon=1/2$: rank 4/5 as stated |
| `a04_moments_riccati.py/.out` | Laurent coefficients from the **Riccati equation** $w_T'=bp\,w_T^2-g/b$, $w_T(\xi_1)=L/N$, $\Phi=N_0w_T(0)$ (no quadrupole, no cosh/sinh series) | $c_0,c_1,c_2|_{R=0}$ of Prop. LF and $\delta_\varepsilon$: exact; stack $G'_{12345}$ local $=2(\tau_\ell-\tau_R)\xi^7P/(135b^2)$ and entropic $=-2\xi^7Q/(2025b^2)$ with $P,Q$ **exactly** those of the record (40 and 49 terms); $G'|_{{\rm ent},\Delta,R_2=0}=-\frac4{45}(b^2+3)\xi^{10}t^2$; $t(\lambda,\mu)$ exact |
| `a05_tw3_recertify.py/.json` | own forward-mode jets over Arb balls, reflection form (film) and tanh form (transducer), exact rational target and frequencies; Newton at 400 bits; Krawczyk at 300 bits on boxes of radius $10^{-35}$ | all three twins certified: $|K(X)-\tilde x|\le1.2\cdot10^{-64}$; $\det J(X)=-1.494\cdot10^{-7}$, $2.30\cdot10^{-8}$, $1.57\cdot10^{-8}$ (certified, regular roots); $\tau_R-\tau_\ell>0$ on $X$; mismatch at $\omega=9/10$ certified $5.264\cdot10^{-4}$, $1.194\cdot10^{-4}$, $3.262\cdot10^{-4}$; roots agree with `b4_07` to $\le4.2\cdot10^{-40}$ |
| `a06_special_eps_scan.py/.out` | own builder of (i), (O), (E); exact ranks on 18 strata for $\varepsilon\in\{$symbolic$,0,\tfrac12,1,-1,2,\tfrac13\}$ (read-outs $T$, $\vartheta_e$, $\vartheta$ at both faces; $R$ known/unknown) and mixed faces | the only value $\varepsilon\ne0$ at which a kernel exceeds its generic-$\varepsilon$ dimension: **$\vartheta_e$ read-out, $\varepsilon=1$, $\Pi\cap\{R=0,b=b_s\}$**, kernel span$\{e_\xi,(0,0,-1,1[,0])\}$; all $\varepsilon=0$ deviations are the local degeneracies listed in the tables; mixed faces: only the thin-slab and $e_\xi$ kernels stated in §4 |
| `a07_certified_minors.py/.json` | own jets, reflection/tanh form, 256 bits | 7 of the 21 certified stack determinants of `b2_05` (local generic, $R_2=0$; entropic $\Delta\cap\{R_2=0\}$, $D_2\cap\{R_2=0\}$, $\Pi\cap\{R_2>0\}$; both Au/AlN/sapphire points) reproduced to 10 digits, all certified nonzero |
| `a08_passivity_1d.py/.json` | $4\cdot10^5$ random points, double precision | bare film: $\operatorname{Re}Z_T>0$ for $\varepsilon\in\{0,-\tfrac12\}$; $\operatorname{Re}Z_T<0$ found for $\varepsilon\in\{0.3,1,2\}$; $\operatorname{Re}[(1+\tau_\ell p)Z_T]>0$ for $\varepsilon\in\{0,0.3,1,2\}$, $<0$ found for $-\tfrac12$. Stack (same $\varepsilon$ both faces): $\operatorname{Re}Z_{\rm tot}>0$ for $\varepsilon\in\{0,0.3,1\}$, $<0$ found for $\varepsilon\in\{-\tfrac12,2\}$ |
| `a09_reduce_recheck.wl/.out` | exact real elimination (Wolfram Language, remote kernel) on the record's own matrices | `T_eps_Rknown`: the elided object is $\mathrm{Root}[t_\ell^2+2t_\ell t_R-3t_R^2+(2t_\ell-2t_R)x^2+x^4,4]$, $=0$ on $t_\ell=t_R$: the record's reading is right. $\vartheta_e$, tied: $\{\operatorname{rank}\le n-2\}$ $=\{\varepsilon=1,t=0\}$ on $\{\Delta,b=1,R=0\}$ ($R$ known), $\{b=1,\varepsilon=1\}$ on $\{\Pi,R=0\}$ ($R$ unknown), empty on $\{\Delta,R=0,t>0\}$; null spaces as in a03 |
| `a10_free_cn_stack.py/.json` | as a03, seven unknowns | free-$c_n$ stack: rank 7 off $\Delta$ incl. $\tau_n=\tau_\ell$, $\tau_n=\tau_R$, $\tau_n=\tau_R$ with $\tau_\ell=0$, $\tau_n<0$; rank 6 on $\Delta$ with $\tau_n=t$, $R_2=0$; rank 5 on $\Delta\cap\{\tau_n=0,R_2=0\}$, $b\ne b_s$ |

---

## 2. Line-by-line audit

### 2.1 Setting, Lemma R, Proposition O (§1)

1. *Assumption M′, Definition 1′.* $\vartheta=T-c_n\partial_zq$ at both faces with the same $\partial_zq_z$ (not $\partial_nq_n$): harmless, since
   $\partial_nq_n=\partial_zq_z$ for either orientation of the normal. VALID.
2. *Film block, reflection form.* $\mathrm{Mob}(M;L/N)=B^{-1}(1+rE)/(1-rE)$ with $r=(u-1)/(u+1)$, $u=BL/N$ (from
   $\cosh A=\tfrac12e^A(1+E)$, $\sinh A=\tfrac12e^A(1-E)$). VALID.
3. *Lemma R.* $\partial_zq=-CpT$ at $k=0$ gives $\vartheta=(1+c_nCp)T$; the bottom jump gives $\vartheta_1=Lq_1$, hence
   $T_1=Lq_1/N_b$; (b) multiplies by $N_0$; (c) $T_A=R_1q+N_tT_{\rm film}(0)$, so $W=R_1+N_t\,\mathrm{Mob}(M;L_2/N_b)$, and with
   $N_t=N_b=N$, $W=R_1+\mathrm{Mob}(\mathrm{diag}(N,1)M\,\mathrm{diag}(N,1)^{-1};L_2)$. VALID.
4. *S-family decoupling at $k=0$.* In `_modes`, P modes have $Q_r=k\,y=0$, $D_r=-\gamma_P k\,y=0$; S modes have $T=0$,
   $Q_z=k/\gamma_S=0$; the `theta` row is $T+c_n(CpT+kQ_r)-c_tkQ_r$, which vanishes on S modes at $k=0$. The S block
   enters only the homogeneous slip rows. VALID; confirmed by the slip-independence test of a01.
5. *Proposition O (i).* Defaults `observable="T"` (`kernel`, `gaussian_response`, `gaussian_response_adaptive`,
   `plane_response`); `synthetic_measurement` (lines 910–911) passes no observable; `impedance` is built from `comp[0]`
   (lines 513–516) and returned (528). VALID.
6. *(ii).* Line 526: $\vartheta_e=T-(\ell^2/\lambda)(Q_z'+\alpha D)$; at $k=0$, $D=Q_z'=-CpT$, so $\vartheta_e=(1+\tau_\ell p)T(0)$,
   independent of `GKBoundary`; $\vartheta_e=T(0)$ for a Fourier/Cattaneo top layer. VALID.
7. *(iii).* Lines 489–492: row $\vartheta_j({\rm bottom})-R\,Q_{z,j}({\rm bottom})-\vartheta_{j+1}({\rm top})=0$, with
   $\vartheta$ evaluated by `bcs[j]` resp. `bcs[j+1]`; a GK layer therefore uses its own coefficient at its top face (row of
   the interface above it) and at its bottom face. `coefficients` (221–228) returns $(0,0)$ for non-GK media and for
   `local`, $(1+\alpha)\ell^2/\lambda$ for `entropic`, the supplied value for `kinetic`. VALID. Numerical confirmation: a01.
8. *Consequence for Workstream A* (Cattaneo results convention-free). VALID (`law != "gk"` returns $(0,0)$ even for a
   supplied kinetic $c_n$). That A's GK variants used the default boundary was not re-checked (A's scripts not in scope).
9. *§1.3 (Beardo et al. coefficients).* Arithmetic consistent ($c_n=2.71/1.068\cdot10^9\times10^{-8}=2.54\cdot10^{-17}$ m³K/W;
   $\tau_n=c_nC\approx42$ ps; $\varepsilon=43/515=0.083$). The reading of eq. (A6) and Table I was not re-checked against the
   source: UNKNOWN (outside the priority list).

### 2.2 Lemma 3′ and the Remark on poles (§2)

10. *Lemma 3′.* $F=N_0(CL+NX)/(YL+NC)$ is meromorphic on $U_\theta$; near $s=0$, $YL\to0$, $NC\to1$, so $F-1/(b_ss)$ is
    holomorphic; the transducer map preserves the principal part since $T_{A,21}W=O(s)$, $T_{A,12}=O(1)$. Identity
    principle as in B-analytic Lemma 3(iii)–(iv). VALID; wording D-8.
11. *Remark, $\varepsilon=0$ and $\varepsilon=1$.* For $\varepsilon=1$: with $q=-\lambda\vartheta'/(1+\tau_Rp)$ the identity
    $\operatorname{Re}\overline{q(0)}\vartheta(0)=C\omega^2\tau_\ell\!\int|T|^2+\int|q|^2/\lambda+R|q(e)|^2+|q(e)|^2\operatorname{Re}(1/S)$ is
    correct (re-derived) and excludes free modes. VALID.
12. *Remark, "for other $\varepsilon$ the interface term $(1-\varepsilon)\tau_\ell\operatorname{Re}(p\,\overline{q(e)}T(e))$ has no sign".*
    With $T(e)=Lq(e)/N$, the term is $(1-\varepsilon)\tau_\ell\operatorname{Re}(pL/N)|q(e)|^2$ and, at $p=i\omega$,
    $\operatorname{Re}(pL/N)=\dfrac{R\,\omega^2\varepsilon\tau_\ell+\omega^{1/2}(\omega\varepsilon\tau_\ell+1)/\sqrt2}{1+\omega^2\varepsilon^2\tau_\ell^2}\ge0$ for $\varepsilon\ge0$.
    So the term is nonnegative for $0\le\varepsilon\le1$. FALSE for $0<\varepsilon<1$ (D-2). Conservative: it only withholds a
    conclusion that in fact holds.

### 2.3 Lemma K (§2) — central

13. *Domain $V=\{-\delta_0<\arg s<\pi/2\}$.* Zeros/poles of $g(s^2)$ lie at $\pm i\tau_R^{-1/2}$, $\pm i\tau_\ell^{-1/2}$, on $\arg s=\pm\pi/2$;
    $V$ is a simply connected sector containing the data ray $\arg s=\pi/4$ and the positive real axis, contained in
    $U_\theta\setminus\{0\}$. $\sqrt g$ has a single-valued branch on $V$ continued from the principal root on the data ray
    (where $|\arg g|<\pi/2$ by Lemma 0). $AB=\xi_1bp$, $A/B=\xi_1g/b$ hold identically for $A=\xi_1s\sqrt g$, $B=bs/\sqrt g$.
    $D_c\Phi\equiv0$ on $V$ by Lemma 3′. VALID. (The list "poles only at zeros of $N$, $N_0$, $u+1$" omits the zeros of
    $1-rE$; immaterial — D-7.)
14. *(K.1).* $\tilde B^2(1-Q)^2D_c\Phi=-D_c\tilde B+2\tilde B\varphi E+r^2D_c\tilde B\,E^2$, $\varphi=D_cr-2rD_cA$. VALID (a02 (1)).
15. *Coefficients algebraic.* $D_c\tilde B$, $\tilde B\varphi$, $r^2D_c\tilde B$ lie in $\mathbb C(s)(\sqrt{g(s^2)})$, a field of degree $\le2$
    over $\mathbb C(s)$; hence Puiseux expansions at every point and at $\infty$, and at most power growth along any path
    ending at a point or at $\infty$. VALID. (An even shorter argument avoiding Puiseux series is given in §4.)
16. *Path 1 ($\tau_R>0$ or $\tau_\ell=0$): positive real $s$-axis.* There $g>0$, the continued $\sqrt g$ is the positive root, and
    $A\sim\xi_1\sqrt{\tau_R/\tau_\ell}\,s$ ($\tau_\ell>0$), $\xi_1\sqrt{\tau_R}\,s^2$ ($\tau_\ell=0<\tau_R$), $\xi_1s$ ($\tau_R=\tau_\ell=0$). VALID.
17. *Path 2 ($\tau_R=0<\tau_\ell$): $s=i\tau_\ell^{-1/2}+\rho e^{i\pi/4}$.* $\operatorname{Re}s=\rho/\sqrt2>0$ and $\arg s<\pi/2$: in $V$.
    $1+\tau_\ell s^2=2\tau_\ell^{1/2}\rho e^{3i\pi/4}+O(\rho^2)$, $\rho A^2\to\xi_1^2e^{i\pi/4}/(2\tau_\ell^{3/2})$ (a02 (6)): $|A|\to\infty$,
    $\arg A\to\pi/8$ or $\pi/8+\pi$, and $|E|=e^{-2\operatorname{Re}A}$ is $o(\rho^N)$ for every $N$ in the first case. VALID.
18. *Second case of path 2.* The substitution $(A,B,r,E)\to(-A,-B,1/r,1/E)$ leaves $\Phi$ invariant and maps (K.1) to
    (K.1)$/(r^2E^2)$ (a02 (2)); it requires $r\not\equiv0$, which holds on $\{\tau_R=0<\tau_\ell\}$ ($u\equiv1$ would force $R=0$,
    $\tau_n=0$, $\tau_\ell=0$; a02 (7)) but is not stated. Equivalently: divide (K.1) by $E^2$ and use $E^{-1}\to0$. VALID;
    presentation D-7.
19. *Conclusion.* $D_c\tilde B=E(\dots)+E^2(\dots)$ is algebraic and $o(\text{any power})$ along the path, hence $\equiv0$; then
    $2\tilde B\varphi E\equiv0$ and $\varphi\equiv0$. VALID.
20. *Translation, (i).* $D_c\log\tilde B^2=2c_b/b-D_c\log g-2\delta\tau_op/N_0\equiv0$; value at $p=0$ gives $c_b=0$. VALID (a02 (4)).
21. *Translation, (O) and (E).* $D_cr=2rD_cA\iff D_cu=(u^2-1)D_cA$; with (i), $D_c\log B=\delta\tau_op/N_0$,
    $D_cA=A\kappa/(\xi_1N_0)$, $uA=\xi_1bpL/N$, $A/u=\xi_1gN/(bL)$, giving
    $pL(N\delta\tau_o-N_0\delta\tau_n)+c_RNN_0=(bpL^2-gN^2/b)\kappa$. Verified as an identity (a02 (4)). The odd/even split in
    $s$ is legitimate (all other factors are rational in $s^2$; an identity $E_1(s^2)+sO_1(s^2)\equiv0$ forces $E_1\equiv O_1\equiv0$);
    odd part $=(s/b_s)\cdot$(O), even part $-Rs^2\cdot$(O) $=$ (E) (a02 (5)). VALID.
22. *Reversibility (⇐).* Each step is an equivalence given $\tilde B,u,L,u+1\not\equiv0$. VALID.
23. *Stack part.* $\tilde B(1-Q)^2D_cW=c_{R_1}\tilde B(1-Q)^2-(D_c\tilde B/\tilde B)(1-Q^2)+2D_cQ$ (a02 (3)); $E^0$:
    $c_{R_1}\tilde B^2=D_c\tilde B$; then $\varphi-c_{R_1}\tilde Br=-c_{R_1}\tilde Br^2E$, both sides forced to vanish; $r\not\equiv0$ gives
    $c_{R_1}=0$ and the bare conditions; $r\equiv0$ gives $c_{R_1}=D_c\tilde B/\tilde B^2$, $D_cu=0$. VALID.
24. *Remark (replacement of Lemma LW in rank statements).* VALID: the argument uses only growth of elements of
    $\mathbb C(s)(\sqrt g)$ against $e^{-2\operatorname{Re}A}$; Lemma T is not needed for any rank statement.

**Lemma K: VALID.** It is elementary and its algebra is confirmed exactly by an independent derivation.

### 2.4 Lemma G, Theorem C1 (§2, §3.1)

25. *Lemma G, cleared relation.* Coefficients of $1,E,E',EE'$: $(\rho\tilde B\tilde B'+\tilde B'-\tilde B)$, $r(\tilde B+\tilde B'-\rho\tilde B\tilde B')$,
    $-r'(\tilde B+\tilde B'+\rho\tilde B\tilde B')$, $rr'(\rho\tilde B\tilde B'-\tilde B'+\tilde B)$ — re-expanded by hand. Denominator
    $\not\equiv0$ ($E\to1$, $r\to(b-b_s)/(b+b_s)$ as $p\to0^+$). VALID.
26. *$A-A'$ non-constant.* Exponents $0,-2A,-2A',-2(A+A')$ pairwise distinct modulo constants (all six differences are
    non-constant: $A,A',A\pm A'$ vanish at $p=0$ and $A,A',A+A'>0$ on $p>0$); LW gives $r\equiv r'\equiv0$ and
    $\rho\tilde B\tilde B'=\tilde B-\tilde B'$. VALID.
27. *$A-A'$ constant.* $=0$; exponents $0,-2A,-4A$; $\rho\tilde B\tilde B'=\tilde B-\tilde B'$, $r\tilde B'=r'\tilde B$,
    $rr'(\tilde B-\tilde B')=0$; integral domain of meromorphic functions gives (a) or (b). VALID.
28. *Theorem C1 (a).* $A^2=A'^2$ at $p=0$ gives $\xi_1=\xi_1'$, $g=g'$; $B^2=B'^2$ gives $b=b'$; $L/N=L'/N'$, odd part
    $\tau_n=\tau_n'$, even part $R=R'$; $g=g'$ by B-analytic Theorem 1 step 5; tied $\varepsilon\ne0$: $\tau_n=\tau_n'$ fixes $\tau_\ell$.
    VALID. *(b)* $u\equiv1$ squared: odd term forces $R=0$, $p=0$ forces $b=b_s$, degrees force $\tau_n=0$ (also for $\tau_n<0$)
    and $\tau_R=\tau_\ell$; then $Z=1/S$ for both points. VALID. Fibre list. VALID.

### 2.5 Theorem C2 and table §0.2 (§3.2)

29. *Reduction with $\tau_o=0$.* (i): $c_b=0$, $c_{\tau_R}(1+\tau_\ell p)\equiv c_{\tau_\ell}(1+\tau_Rp)$; (O): $-\delta\tau_n=2bRc_\xi$;
    (E): $c_RN=c_\xi(b/b_s^2-bR^2p-gN^2/b)$. VALID.
30. *Tied, off $\Delta$.* $Rc_\xi=0$; for $c_\xi\ne0$, $gN^2=(1+\tau_Rp)(1+\varepsilon\tau_\ell p)^2/(1+\tau_\ell p)$ must be a polynomial of
    degree $\le1$: divisibility forces $\tau_\ell=0$ (then $\tau_R=0$, on $\Delta$) or $\varepsilon=1$ (then $\tau_R=0$, and evaluation at
    $p=-1/\tau_\ell$ gives $0=b^2/b_s^2$). VALID.
31. *Tied, on $\Delta$, $t>0$.* $p^2$ coefficient of (E) $-c_\xi\varepsilon^2t^2/b$; then (O) $\varepsilon c_t=0$, (E) $c_R=0$. VALID.
    *Corner $t=0$.* $p^1$: $c_\xi bR^2=0$; $p^0$: $c_R=c_\xi(b^2-b_s^2)/(bb_s^2)$; (O): $\varepsilon c_t=-2bRc_\xi$. For $\varepsilon\ne0$
    no special value can enlarge the kernel (the coefficient of $c_t$ is $\varepsilon$). VALID.
32. *Free, off $\Delta$.* Polynomial degree argument: $\tau_\ell=0\Rightarrow\tau_n=0\Rightarrow\tau_R=-b^2R^2\le0\Rightarrow$ on $\Delta$;
    $\tau_\ell>0\Rightarrow\tau_n=\tau_\ell$, $\tau_R=0$, evaluation at $-1/\tau_\ell$ contradicts. VALID (covers $\tau_n<0$).
    *Free, on $\Delta$.* $p^2$: $c_\xi\tau_n^2=0$; for $\tau_n=0$: $R=0$ and the thin-slab ratio. VALID.
33. *Local.* B-analytic Theorems 2, 2′ (audited). VALID.
34. *Table §0.2.* All 36 entries re-derived from 29–33; a03 (4 points), a06 (all strata, special $\varepsilon$ values: no
    anomaly at $\varepsilon\ne0$), and the re-run elimination for `T_eps_Rknown` (a09) agree. **PROVED.**
35. *§3.3 (2).* The elimination outputs for read-out $T(0)$ are full-system decisions; the one elided object is resolved
    (a09; D-9). VALID.
36. *Corollary C3 (1)–(3).* Consequences of C1, C2. VALID. (4) and the conditioning numbers are numerical observations
    ($0.030/1$, $0.0027/0.084$, $0.016/0.5$, $0.060/2.5$ are $0.024$–$0.032$; "approximately proportional" is fair).

### 2.6 Proposition C4 ($\vartheta_e$ read-out) (§4)

37. *(i)-block.* $c_{\tau_R}p/(1+\tau_Rp)+c_{\tau_\ell}p/(1+\tau_\ell p)\equiv0$: off $\Delta$ $c_{\tau_R}=c_{\tau_\ell}=0$; on $\Delta$ $c_{\tau_R}=-c_{\tau_\ell}=:-c$.
    VALID. Off $\Delta$ the residual system is that of item 30: full rank. VALID.
38. *On $\Delta$.* $\kappa=c_\xi+(c_\xi t-\xi_1c)p$; (O): $c(1-\varepsilon)=2bRc_\xi$ ($p^0$), $bR(c_\xi t-\xi_1c)=0$ ($p^1$). VALID.
39. *$R>0$, $t>0$.* $c=c_\xi t/\xi_1$, $\kappa=c_\xi$; (E) $p^0,p^1,p^2$ give, for $c_\xi\ne0$ and $\varepsilon\ne0$: $b^2=b_s^2(1-\varepsilon)$,
    $b^2R^2=\varepsilon(\varepsilon-1)t$, $(1-\varepsilon)t=2bR\xi_1$, solvable iff $\varepsilon<0$, giving $E_\varepsilon$ with
    $t=4|\varepsilon|\xi_1^2/(1+|\varepsilon|)$, $R=2|\varepsilon|\xi_1/b$; for $\varepsilon=0$: $(b_s^2-b^2)t=b^2b_s^2R^2$, $t=2bR\xi_1$, i.e. $E_0$. Kernels as
    stated. VALID (re-derived; a03 at $E_0$ and $E_{-1}$: rank 4/5; control point on $\Delta\cap\{R>0\}$ off $E_0$: rank 5/5).
    *$R>0$, $t=0$:* $p^1$ of (O) gives $c=0$, then $c_\xi=0$: full rank. VALID.
40. *$R=0$, $\varepsilon\ne1$.* $c=0$, $\kappa=c_\xi N_0$, $c_RN=c_\xi(b/b_s^2-N^2/b)$: full for $\varepsilon\ne0$, $t>0$; thin slab for
    $\varepsilon=0$ or $t=0$. VALID.
41. *$R=0$, $\varepsilon=1$, $t>0$.* $p^3$ of (E) forces $\kappa=c_\xi$; then $p^0,p^2$ force $c_\xi=0$: full. VALID.
42. ***$R=0$, $\varepsilon=1$, $t=0$.*** Here (O) reads $0=0$, so $c$ is **not** forced to vanish ("the same computation with
    $c(1-\varepsilon)=0$" does not give $c=0$ at $\varepsilon=1$). (E) at $t=0$: $c_R=(c_\xi-\xi_1cp)(b^2-b_s^2)/(bb_s^2)$; its $p^1$
    coefficient $\xi_1c(b^2-b_s^2)$ vanishes identically when $b=b_s$. Hence at $\{\varepsilon=1,\ \tau_R=\tau_\ell=0,\ b=b_s,\ R=0\}$
    the kernel is $\mathrm{span}\{e_\xi,(0,0,-1,1)\}$ ($R$ known: rank **2** of 4, not 3) and $\mathrm{span}\{e_\xi,(0,0,-1,1,0)\}$
    ($R$ unknown: rank **3** of 5, not 4). Direct first-order check: at $r\equiv0$ the block is $N_0/B$ to first order, which
    depends on $\tau_R+\tau_\ell$ only, and $D_cu=(1-\varepsilon)c\,p\,u$ vanishes at $\varepsilon=1$. Confirmed three ways: hand
    (above), 60-digit Jacobian (a03: deficient singular values $6.7\cdot10^{-40}$, both kernel vectors annihilated to
    $8\cdot10^{-40}$), exact elimination of $\{\operatorname{rank}\le n-2\}$ on the record's matrices (a09). Statements (a)
    ("kernel $\mathbb Re_\xi$") and (b) ("corner $t=0$, $R=0$: thin-slab kernel as in §0.2") are **FALSE** at this point (D-1).
    For $\varepsilon=1$ and $b\ne b_s$ the corner kernel is the thin-slab line, as stated.
43. *Gröbner branch record.* The first basis element $b\varepsilon^2t^2(1-\varepsilon)(1-\varepsilon-b^2)$ has the factor $t^2$; the
    record lists the branches $\varepsilon=0$, $\varepsilon=1$, $b^2=1-\varepsilon$ only. The $t=0$ branch is treated in the text by hand
    (item 42), which is where the error sits (D-9).
44. *Status.* The degenerate **sets** of C4 (a)–(d) are correct for every $\varepsilon$ (the point of item 42 lies in the stated
    set). The rank/kernel at $\varepsilon=1$ on $\Pi\cap\{b=b_s,R=0\}$ is wrong. Note that C4 at $\varepsilon=1$ is the same function as
    C5 at $\varepsilon=1$, whose corner kernel $\mathrm{span}\{e_\xi,(0,0,1-2\varepsilon,1)\}$ is stated correctly: the two propositions
    contradict each other at that point.

### 2.7 Proposition C5 and different coefficients (§4)

45. *(i)-block.* Coefficients $c_{\tau_R}$, $-c_{\tau_\ell}$, $2\varepsilon c_{\tau_\ell}$ of $p/(1+\tau p)$ at $\tau_R$, $\tau_\ell$, $\varepsilon\tau_\ell$; these functions are
    linearly independent for distinct $\tau$ (including $\tau=0$). Case enumeration (all distinct; $\tau_R=\tau_\ell\ne\varepsilon\tau_\ell$;
    $\tau_R=\varepsilon\tau_\ell\ne\tau_\ell$; $\tau_\ell=\varepsilon\tau_\ell\ne\tau_R$; all equal) gives exactly the locus
    $\Pi\cup(\Delta\cap\{\varepsilon\in\{0,1\}\})\cup\{\tau_\ell=0<\tau_R,\ \varepsilon=\tfrac12\}$, as in the record. VALID.
46. *Residual system* (= item 30): no degeneracy off the locus for $\varepsilon\ne0$. VALID.
47. *On the locus.* $R>0$: (O) is $2bR\kappa\equiv0$, so $c_\xi=0$, $\varepsilon c_{\tau_\ell}=0$, then $c=0$ — full rank on the whole locus,
    uniformly. $R=0$: $\Delta\cap\{\varepsilon=1,t>0\}$ full ($p^3$ then $p^0,p^2$ of (E)); $\{\tau_\ell=0<\tau_R,\varepsilon=\tfrac12\}$ full
    ($p^2$, $p^1$ of (E)); corner: $c_{\tau_R}=(1-2\varepsilon)c_{\tau_\ell}$, $p^1$ of (E) $\varepsilon c_{\tau_\ell}(b^2-b_s^2)=0$, giving the thin-slab
    line ($b\ne b_s$) or $\mathrm{span}\{e_\xi,(0,0,1-2\varepsilon,1,0)\}$ ($b=b_s$). All re-derived by hand. VALID; a03, a06 (no anomaly
    at $\varepsilon\in\{\tfrac12,1,-1,2,\tfrac13\}$).
48. *Free $\tau_n$.* (i)-block determinant $\propto(\tau_\ell-\tau_n)(\tau_\ell-\tau_R)(\tau_n-\tau_R)$; on $\tau_n=\tau_R\ne\tau_\ell$ and
    $\tau_n=\tau_\ell\ne\tau_R$ (including $\tau=0$ sub-cases) re-derived by hand: (E) evaluated at $p=-1/\tau_\ell$ and $p=-1/\tau_R$
    forces $\kappa\equiv0$ or $c_{\tau_n}=0$, then $c=0$; on $\Delta$ the stated kernels, including the three-dimensional corner
    kernel. VALID; a10 (stack, seven unknowns) confirms rank 7 on both coincidence strata and for $\tau_n<0$.
49. *Global part of C5.* Lemma G with $\tilde B^2=b^2p/(gN^2)$: $A^2$ gives $\xi_1$, $g$; $\tilde B^2$ then $b$ and $N$. VALID.
50. *Different coefficients at the two faces.* Labelled DERIVED; the (O) equation at the corner, $(\varepsilon_0-\varepsilon)c_{\tau_\ell}=2bR\kappa$,
    explains why the extra corner direction occurs only for $\varepsilon_0=\varepsilon$. The a06 scan of $(\varepsilon_0,\varepsilon)$ pairs on
    $\Delta$, $\Delta\cap\{R=0\}$, $\Pi\cap\{R=0\}$, $\Pi\cap\{R=0,b=b_s\}$, $\{\tau_\ell=0\}$ finds only the kernels stated. Consistent.

### 2.8 Transducer: Theorems T1–T3 (§5)

51. *T1.* $\frac{d}{dW}\mathrm{Mob}(T_A;W)=\det T_A/(T_{A,21}W+T_{A,22})^2$; Möbius maps with $\det\ne0$ are injective. Function
    level: VALID. Per frequency: the complex row is multiplied by a nonzero number (invertible real $2\times2$ block) wherever
    $T_{A,21}W+T_{A,22}\ne0$, i.e. wherever $Z_{\rm tot}$ is finite; that the stack is pole-free on $i\mathbb R_+$ for the local and
    entropic conventions is true but is not shown in the text (D-4).
52. *T2.* Case (a) of Lemma G gives $\rho=0$. Case (b): sign slip in "$2\rho/\tilde B'=1/\tilde B^2-1/\tilde B'^2-\rho^2$" (correct:
    $2\rho/\tilde B'=1/\tilde B'^2+\rho^2-1/\tilde B^2$, D-6); the conclusion ($p(1+\tau_R'p)(1+\tau_\ell'p)$ would be a square, impossible)
    is unaffected. Film-block fibres: $A^2\to\xi_1,g$; $\tilde B^2\to b,N$; $r\to R_2$; tied $\varepsilon\ne0$: $N\to\tau_\ell$. $r\equiv0$ iff
    $R_2=0$, $b=b_s$, $\tau_n=0$, $g=1$. VALID; **PROVED** given Lemma LW.
53. *T3, $r\not\equiv0$.* Kernel of $DW$ = $\{c_{R_1}=0\}\times$ kernel of the film block with $R_2$ unknown, so rank $=1+$rank of
    the block. VALID.
54. *T3, $r\equiv0$.* $c_{R_1}=(2b_ss)^{-1}\big(2c_b/b_s-(c_{\tau_R}-c_{\tau_\ell})p/(1+tp)-2\delta\tau_np\big)$ re-derived; odd in $s$, hence
    $c_{R_1}=0$, $c_b=0$; then $c_{\tau_R}=(1-2\varepsilon)c_{\tau_\ell}$ (tied, $t=0$), $c_{\tau_R}=c_{\tau_\ell}$ (local), free as stated; $D_cu=0\iff c_{R_2}=0$
    (the $\delta\tau_n$ terms cancel); $c_\xi$ free. VALID.
55. *Table §0.3.* All rows re-derived; a03 (4 stack points incl. $\varepsilon=-\tfrac12$ corner, kernel $(0,0,2,1,0,0)$) and a10 (7 free-$c_n$
    points) agree. "$R_1$ never enters a kernel": VALID. **PROVED.**
56. *Elementary partial second proof (moments).* $\operatorname{rank}DZ\ge1+\operatorname{rank}[\nabla'c_1;\dots;\nabla'c_5]$ is correct (the $R_1$
    column of the moment matrix is $e_0$). Both restrictions of $G'_{12345}$ reproduced (a04). VALID.

### 2.9 Theorem T4 (§5)

57. *(i)* $2N\ge6$. VALID.
58. *(ii)* Ball-arithmetic certificates: 7 of 21 re-certified independently (a07), identical to 10 digits. Real analyticity
    of $\det\mathcal J_3$ in $(\hat\theta,\omega)$ on a connected stratum requires $Z_{\rm tot}$ to be pole-free on $i\mathbb R_+$; the text cites
    B-analytic D1, which concerns the bare film (D-4). For the local convention the listed strata $\{b=b_s\}$, $\{R_1=0\}$,
    $\{R_2=0\}$ minus $\Delta$ have two components ($\tau_R\gtrless\tau_\ell$), and only $\tau_R>\tau_\ell$ was certified, so Lemma F covers
    only that component (D-5). Neither changes the conclusion: (iv) covers both components generically (e.g.
    $P|_{b=b_s,R_2=0}=-45\xi_1(\tau_\ell-\tau_R)^2\ne0$). VALID for the certified components.
59. *(iii)* Ratio of the $R_2$ and $R_1$ columns $\partial_{R_2}W=1/(YL_2/N+C)^2=1-2\xi_1be^{i\pi/4}\sqrt\omega/b_s+O(\omega)$: correct
    ($YL_2/N\sim\xi_1bs/b_s$); non-real for small $\omega$, hence per-frequency rank 2 off a discrete set. Kernel-dropping:
    for $c\in K\setminus0$, $\{\omega:\mathcal J_1(\omega)c=0\}$ is discrete because $D_cZ_{\rm tot}\not\equiv0$; $4\to3\to2\to1\to0$. Lemma F
    on the sum of squares of the $6\times6$ minors gives the closed nowhere dense null exceptional set (polar hyperplanes, if
    any for other $\varepsilon$, are themselves closed and null). VALID.
60. *(iv)* Rows of $\mathcal J_3$ for $W$: $\cos(k\pi/4)(\lambda_j\sigma)^k\nabla c_k$, $\sin(k\pi/4)(\lambda_j\sigma)^k\nabla c_k$; Cauchy–Binet over
    6-sets with absolute convergence and termwise differentiation (uniform disc, as D2–D3 of the B-analytic audit: here
    $YL+NC\to1$ locally uniformly); minimal weight 15 only for $K=\{0,\dots,5\}$; $G_K=\pm G'_{12345}$ by expansion along the
    $R_1$ column ($\partial_{R_1}c_k=\delta_{k0}$); remainder $O(\sigma^{16})$. $t(\lambda,\mu)$ exact (a04, $=\sqrt2\lambda\mu(\dots)/4$), nonzero
    for distinct $\lambda,\mu\ne1$. Local $G'=2(\tau_\ell-\tau_R)\xi^7P/(135b^2)$ and entropic $G'=-2\xi^7Q/(2025b^2)$ reproduced
    **exactly** by the Riccati route (a04), which makes three systems for the local and two independent derivations for the
    entropic polynomial. $Q|_{\Delta,R_2=0}=90b^2t^2\xi^3(b^2+3)\ne0$ ($t>0$), so $\{Q=0\}\not\supset\Delta$; on $\Delta\cap\{R_2>0\}$, $Q$
    does vanish on a hypersurface (the bracket of $G'|_\Delta$ has the negative terms $-R_2b^4\xi^5$, $-45R_2\tau_\ell^2\xi$), which the
    text leaves to Conjecture T4′, correctly. Excluded sets as stated (local $\Delta\cup\{P=0\}$, entropic $\{Q=0\}$; $b_s=1$ is
    no loss of generality: $b\to b/b_s$, $R\to Rb_s$, $b_A\to b_A/b_s$, $Z\to Zb_s$). **PROVED.**

### 2.10 Proposition TW3 and the numerical sections (§5.4)

61. *TW3.* Krawczyk logic in `b4_07` correct ($K(X)\subset\operatorname{int}X$ tested componentwise; $Y$ a real midpoint inverse;
    exact decimal balls for target and frequencies). Independent re-certification of all three twins (a05): different
    formulation, own jets, radius $10^{-35}$, 300 bits; contraction to $10^{-64}$; $\det J(X)$ certified nonzero (regular
    roots, so the implicit-function persistence is justified); $\tau_R-\tau_\ell>0$ certified on the box (off $\Delta$); mismatch at
    $\omega=9/10$ certified nonzero, values identical to the record. The printed 12-digit $\hat\theta'$ is labelled "$(\pm10^{-40})$",
    which refers to the 40-digit centres stored in `b4_07_certify_twin.json` (D-11). **PROVED.**
62. *"Three frequencies give generic local but not global identifiability."* Local part from T4(ii)/(iv); non-uniqueness
    on an open set of $(\hat\theta^*,\omega)$ under the local convention. For the entropic convention no twin is known; the text
    says so. VALID as stated.
63. *Numerical observations, Conjectures T4′, T4″.* Correctly labelled; not proofs.

### 2.11 §6–§8

64. *Proposition LF.* $c_{-1}$, $c_0$, $c_1$, $c_2|_{R=0}$ reproduced exactly by the Riccati route for general $\tau_n,\tau_o$ (a04);
    $\delta_\varepsilon=\tau_R-\tau_\ell+\varepsilon(1+b^2/b_s^2)\tau_\ell$ exact. Consequence (c): first-order dependence on
    $(\tau_R-\tau_\ell,\tau_n,\tau_o)$ is immediate; independence of the two time directions at $\Pi$ follows from T3 except on
    $\Pi\cap\{b=b_s,R_2=0\}$, as stated. VALID.
65. *Practical conditioning (§6 table).* Numerical; not audited here (numerical-analysis gate).
66. *§7, B-2 (iv).* "The entropic one is the one with a sign-definite 1D dissipation identity." At $k=0$ the film obeys
    $q=-(\lambda/g)T'$, and
    $\overline{T(0)}q(0)=Cp\!\int|T|^2+\overline{T(e)}q(e)+(\bar g/\lambda)\!\int|q|^2$; with the local jump $T(e)=Lq(e)$,
    $\operatorname{Re}Z_T=\operatorname{Re}L\,|q(e)|^2+(\operatorname{Re}g/\lambda)\!\int|q|^2>0$ ($\operatorname{Re}g=(1+\omega^2\tau_R\tau_\ell)/(1+\omega^2\tau_\ell^2)>0$). The
    same identity is sign-definite for every $\varepsilon\le0$ ($\operatorname{Re}(L/N)>0$), and the $\vartheta$-identity for every
    $0\le\varepsilon\le1$ (item 12). a08: $\min\operatorname{Re}Z_T/|Z_T|=1.0\cdot10^{-5}>0$ for $\varepsilon=0$ over $4\cdot10^5$ samples. The statement
    is therefore FALSE as a characterization (D-3); the local/entropic distinction of the twin's documentation concerns the
    per-wavenumber identity at $k>0$.
67. *§8 status table and "weakest steps".* "Propositions C4, C5 — PROVED" must be qualified by D-1; weakest step (3) should
    add that kernel dimensions on degenerate strata were not decided by elimination (D-10).

### 2.12 The piecewise eliminations for general $\varepsilon$ (§3.3, §4, weakest step 3)

68. *Structure.* For a direction $c$, (i) involves only $(c_b,c_{\tau_R},c_{\tau_\ell}[,c_{\tau_n}])$. Where the (i)-block has full column rank,
    (i) forces these to vanish, $\delta\tau_o=\delta\tau_n=0$, $\kappa=c_\xi N_0$, and (O), (E) reduce to the residual system in
    $(c_\xi,c_R)$, identical to that of read-out $T(0)$. Where it is deficient, the full system is analysed on each component
    of the deficiency locus. Hence $\{\operatorname{rank}<n\}=(\text{residual locus}\setminus L_{\rm i})\cup(\text{full locus}\cap L_{\rm i})$
    — an exhaustive decomposition. VALID.
69. *Coverage.* $T(0)$ read-out: full systems eliminated in one piece (no decomposition needed). $\vartheta_e$: $L_{\rm i}=\Delta$;
    residual off $\Delta$: empty; on $\Delta$: reduced $3\times3$ minors, branches $\varepsilon=0$, $\varepsilon=1$, $b^2=1-\varepsilon$, and $t=0$ (by
    hand). $\vartheta$ at both faces: $L_{\rm i}$ has the four components of item 45, each settled (elimination and hand);
    free $\tau_n$: $L_{\rm i}=\{\tau_n=\tau_\ell\}\cup\{\tau_n=\tau_R\}\cup\Delta$, each settled. Every branch was re-derived by hand in
    items 37–48. The degenerate **sets** are proved. VALID.
70. *What the decomposition does not decide.* No step decides $\{\operatorname{rank}\le n-2\}$ for the tied-$\varepsilon$ models; the exact
    strata check (`b2_02`) keeps $\varepsilon$ symbolic, i.e. generic. Kernel dimensions on degenerate strata therefore rest on the
    hand analysis alone, and the one special value where it fails (item 42) passed through both checks. INCOMPLETE (D-10);
    the special-value scan a06 and the elimination a09 close it: no other anomaly.

---

## 3. Defects and repairs

| # | Location | Class | Defect | Minimal repair |
|---|---|---|---|---|
| **D-1** | Prop. C4 (a), (b) and its proof; §0.4, §8 status | **FALSE** (localized) | At $\varepsilon=1$ on $\Pi\cap\{b=b_s,R=0\}$ the kernel is $\mathrm{span}\{e_\xi,(0,0,-1,1)\}$ ($R$ known, rank 2 of 4) resp. $\mathrm{span}\{e_\xi,(0,0,-1,1,0)\}$ ($R$ unknown, rank 3 of 5), not $\mathbb Re_\xi$ (rank 3 / 4). Proof step: "the same computation with $c(1-\varepsilon)=0$" assumes $\varepsilon\ne1$ | (a): "kernel $\mathbb Re_\xi$, and $\mathrm{span}\{e_\xi,(0,0,-1,1)\}$ if $\varepsilon=1$, $t=0$". (b): "corner $t=0$, $R=0$: thin-slab kernel; for $b=b_s$ and $\varepsilon=1$ in addition $(0,0,-1,1,0)$". Add to the proof: "for $\varepsilon=1$, (O) is void; at $t=0$ the $p^1$ coefficient of (E) is $\xi_1c(b^2-b_s^2)$". Note the coincidence with C5 at $\varepsilon=1$. Status: degenerate sets PROVED; kernel table PROVED after correction |
| D-2 | Remark after Lemma 3′ | FALSE (conservative) | "For other $\varepsilon$ the interface term has no sign": it is $\ge0$ for $0\le\varepsilon\le1$ | Replace by: "the term is nonnegative for $0\le\varepsilon\le1$ ($\operatorname{Re}(pL/N)\ge0$ for $\varepsilon\ge0$); for $\varepsilon\le0$ the identity in $T$ is sign-definite ($\operatorname{Re}(L/N)>0$); hence $Z$ has no poles on $i\mathbb R_+$ for every $\varepsilon\le1$ (bare film), and the both-face block and stack have none for $0\le\varepsilon\le1$". This also removes the data-set caveat for the Beardo values $0.03\le\varepsilon\le0.23$ |
| D-3 | §7, B-2 (iv) | FALSE (interpretive) | Sign-definiteness of a 1D identity does not single out the entropic convention | "At $k=0$ the local convention is passive in $T(0)$ and the entropic one in $\vartheta_e(0)$; sign-definite 1D identities exist for every $\varepsilon\le1$; the distinction drawn in the twin's documentation concerns the identity at $k>0$" |
| D-4 | T1 (per frequency), T4 (ii) (analyticity on a connected domain) | INCOMPLETE | Pole-freeness of the **stack** on $i\mathbb R_+$ is used, but only the bare film (D1) and the film block (Remark) are covered | One line: layered identity — local: $T$-identity in each layer plus $R_1|q|^2$; entropic: $\vartheta$-identity in the film plus the Fourier identity in the transducer; gives $\operatorname{Re}Z_{\rm tot}>0$, hence no free mode (a08 consistent: $\min\operatorname{Re}Z_{\rm tot}/|Z_{\rm tot}|>0$ for $\varepsilon\in\{0,0.3,1\}$) |
| D-5 | T4 (ii), local strata list | INCOMPLETE | $\{b=b_s\}$, $\{R_1=0\}$, $\{R_2=0\}$ minus $\Delta$ have two components; one certified | Restrict the wording to the certified component, or certify one point with $\tau_R<\tau_\ell$ on each, or cite (iv) |
| D-6 | T2, proof of case (b) | typo | Sign in $2\rho/\tilde B'=\dots$ | $2\rho/\tilde B'=1/\tilde B'^2+\rho^2-1/\tilde B^2$; or note directly that $r\equiv r'\equiv0$ forces $\tilde B=\tilde B'=S$, so $\rho=0$ |
| D-7 | Lemma K proof | presentation | Substitution $(A,B,r,E)\to(-A,-B,1/r,1/E)$ needs $r\not\equiv0$ (true on $\{\tau_R=0<\tau_\ell\}$); relating factor is $r^2E^2$; pole list incomplete | State "$r\not\equiv0$ on $\{\tau_R=0<\tau_\ell\}$" or divide (K.1) by $E^2$ |
| D-8 | Lemma 3′ | wording | "accumulation point at which both responses are finite" literally excludes $\omega_*=0$ | "at which $F-F'$ (resp. $D_cF$) has no pole; this holds at $\omega_*=0$ by the common principal part" |
| D-9 | `b2_03_reduce_record.md` | record | (a) `Root[...]` elided in the quoted `T_eps_Rknown` output; (b) the $t=0$ factor of the C4 Gröbner element is not listed as a branch | Quote $\mathrm{Root}[t_\ell^2+2t_\ell t_R-3t_R^2+(2t_\ell-2t_R)x^2+x^4,4]$ ($=0$ on $\Delta$; a09); list the $t=0$ branch with its kernel (D-1) |
| D-10 | §3.3, §4, weakest step (3) | INCOMPLETE | Eliminations decide $\{\operatorname{rank}<n\}$ only; kernel dimensions at special $\varepsilon$ unchecked | Add the $\{\operatorname{rank}\le n-2\}$ elimination on each degenerate set (a09 does it for C4) or the special-value scan (a06) |
| D-11 | TW3 | presentation | "$(\pm10^{-40})$" after 12 printed digits | Print the 40-digit values or say "centres in `b4_07_certify_twin.json`, radius $10^{-40}$" |

No other step of §1–§6 is FALSE or UNJUSTIFIED.

---

## 4. Optional strengthening of Lemma K, and formal verification

The growth step of Lemma K (item 15) can be made Puiseux-free. All coefficients lie in $K=\mathbb C(s)(\sqrt g)$. If $g$ is a
square in $\mathbb C(s)$, $K=\mathbb C(s)$ and a rational function that is $o(\text{every power})$ along the path is zero. Otherwise
write $f=a+b\sqrt g$, $\sigma f=a-b\sqrt g$ ($a,b\in\mathbb C(s)$). Along either path $\sigma f$ grows at most polynomially, so
$\mathrm N(f)=f\,\sigma f=a^2-b^2g\in\mathbb C(s)$ is $o(\text{every power})$, hence $\mathrm N(f)\equiv0$; $g$ not a square forces
$b\equiv0$, then $a\equiv0$. This removes the only analytic ingredient beyond rational functions and makes the step
directly formalisable.

*Formal verification.* Not required for the conclusions: the algebra of Lemma K is checked exactly by two independent
derivations (record and a02), and the asymptotic step is three lines. If a machine check is wanted, the natural
candidate is the norm argument above together with the finite statement of D-1 (a $4\times4$ / $5\times5$ rational
kernel computation), both small enough for the `formal-verifier` workflow.

---

## 5. Overall verdict

Lemma K is valid and elementary; it supplies the rank statements without Lemma LW, as claimed (referee item B-m4, rank
part). Theorem C2 with table §0.2, Proposition C5, Theorems T1–T3 with table §0.3, Theorem T4 (i)–(iv) and Proposition TW3
are **proved** (T2 given the audited Lemma LW), up to the minor repairs D-4 to D-8 and D-11. Proposition O is proved from
the code and confirmed independently to $1.3\cdot10^{-15}$: the twin reads $T(0)$ and applies the layer's convention at
both faces of a GK layer. Proposition C4 is correct for its degenerate sets and exotic strata but states the wrong
kernel at one exceptional point ($\varepsilon=1$, $\Pi\cap\{b=b_s,R=0\}$; D-1), the only error found in the rank tables; it does
not touch the twin's observable, the transducer results, or the answers to B-2 (1)–(3) and B-4. The piecewise
eliminations cover every case for the degenerate sets; they do not decide kernel dimensions, which is where D-1 sits
(D-10). The interpretive sentence of §7, B-2 (iv) should be corrected (D-3), and the Remark after Lemma 3′ can be
strengthened (D-2): no poles on $i\mathbb R_+$ for every $\varepsilon\le1$ in the bare film.
