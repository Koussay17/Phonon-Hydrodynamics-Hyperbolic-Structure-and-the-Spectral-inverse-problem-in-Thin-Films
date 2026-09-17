# RED peer review — source and calculation notes

Date: 2026-09-17. These notes belong only to the independent peer reviewer. No other red-team reports were read.

## Primary-source read scope

1. Fugallo, Lazzeri, Paulatto, Mauri, PRB 88, 045430 (2013), DOI 10.1103/PhysRevB.88.045430. Actual arXiv v2 HTML: https://arxiv.org/html/1212.0470v2 . Appendix A, Eqs. (31)-(34), explicitly decomposes the collision matrix into nonnegative probabilities times rank-one PSD event matrices. Triplet vector (1,1,-1); isotope pair vector (1,-1). Appendix B also displays the three-phonon energy delta and isotope equal-energy delta. The review's finite-library inequality is its own derivation, not a result claimed by Fugallo et al.
2. Herbert Spohn, Energy Current Correlations For Weakly Anharmonic Lattices, arXiv:0706.0815v1, 6 June 2007. Actual PDF: https://arxiv.org/pdf/0706.0815 . Read targeted Sections 3-4, especially (3.23), (3.24), (4.1), (4.2), and spectral-gap discussion. The model is weakly anharmonic scalar lattice dynamics, not AlN. It supplies an independent primary event quadratic form and relevant historical transport structure.
3. Simoncelli, Marzari, Cepellotti, Generalization of Fourier's Law into Viscous Heat Equations, PRX 10, 011019 (2020), DOI 10.1103/PhysRevX.10.011019. Actual arXiv v2 HTML: https://arxiv.org/html/1906.09743v2 . Read Sections I-II and Appendix A.1, especially the explicit normalized Bose energy vector Eq. (30). That formula, rather than nearby shorthand statements of proportionality, is used for the equal-weight argument.

## Access limitations

Opening https://arxiv.org/html/2008.07596v2 and https://www.nature.com/articles/s41467-021-23618-7 returned internal errors. A Nature search result described the usual resonance rule, but the review uses the inspected primary event formulas instead. The search is bounded; no absence or novelty inference follows. No local multi-index query was run in this bounded peer pass because the frozen L report already covered that infrastructure and the PI requested no expansion beyond the verdict.

## Independent symbolic check

Run with PowerShell using this Python source through `py -`. It does not modify campaign artifacts.

```python
import sympy as s
x, p = s.symbols('x p', positive=True)
H = s.Matrix(8, 8, lambda i, j:
             (-1)**((int(i) & int(j)).bit_count()) / s.sqrt(8))
A = s.diag(0, x*x, 2, 2-x*x, 1, 1, 1, 1)
A[1,2] = A[2,1] = x
A[4,7] = A[7,4] = -x
C = s.simplify(H*A*H.T)
K = (p+x*x)/(p*p+(2+x*x)*p+x*x)
print('diagonal:', [s.simplify(C[i,i]) for i in range(8)])
print('energy:', s.simplify(C*s.ones(8,1)))
print('C04:', s.factor(C[0,4]))
print('response difference:', s.factor(K-1/(p+2)))
print('negative zero-frequency derivative:',
      s.simplify(-s.diff(K,p).subs(p,0)))
```

Observed output: eight unit diagonals; zero energy residual; C04=x/2; response difference x²/[(p+2)(p²+p x²+2p+x²)]; derivative 1+x^-2. The first attempt omitted int() conversions and raised AttributeError because a SymPy Zero lacks bit_count. The corrected attempt completed.

## Hand checks and their logical status

- For a fixed finite event list, define R as in report (R1). Each term obeys |(z.v)(b.v)| <= R(z.v)²; terms with z.v=0 vanish. Nonnegative weighting and summing proves the inequality. It excludes D's entire near-zero path, not every individual matrix and not every possible unbounded-memory family. Independent proof audit is still required by campaign protocol.
- The Bose factor f(x)=x/[2sinh(x/2)] has strictly negative derivative for x>0 because tanh(x/2)<x/2. Equal positive entries at equal weights and finite temperature imply equal positive frequencies. Classical equipartition is explicitly outside this inference.
- The harmonic bound follows by factoring the exact quadratic into (p+lambda_-)(p+lambda_+) with lambda_+>=2. At p=i omega, the three denominator factors have moduli at least sqrt(omega²+4), |omega|, sqrt(omega²+4), respectively. It does not require numerical fitting.

## Tool/runtime scope

Windows sandbox creation failed with helper_sandbox_lock_failed / SetNamedSecurityInfoW access denied. Authorized reads, in-memory checks, and writes to reviewer-owned files therefore used approved require_escalated calls. No remote messages, repository mutations, or publication actions were performed by this reviewer.
