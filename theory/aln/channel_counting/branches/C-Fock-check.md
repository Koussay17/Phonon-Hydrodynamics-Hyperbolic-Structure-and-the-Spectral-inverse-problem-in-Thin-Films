# C: independent Fock-space channel-counting check

**Completed independent first pass.** Read only the campaign question and assumptions; no peer report, material code, or supplied coefficient was used. Results below concern an explicitly defined toy Hamiltonian. The mapping to pinned phono3py amplitudes remains unvalidated by this branch.

## 1. Convention and hypothesis

Assume positive mode energies with \(E_{\bar\lambda}=E_\lambda\), a fixed reciprocal sewing gauge, and
\[
X_\lambda=a_\lambda+a^\dagger_{\bar\lambda},\qquad
H_3={1\over3!}\sum_{\lambda\mu\nu\ {\rm ordered}}
T_{\lambda\mu\nu}X_\lambda X_\mu X_\nu.
\]
Here \(T\) has energy units, is fully permutation symmetric, and satisfies
\(T_{\bar\lambda\bar\mu\bar\nu}=T_{\lambda\mu\nu}^*\).
These are assumptions, not asserted properties of an uninspected export.

For \(p\to a+b\), the all-incoming tuple is \((p,\bar a,\bar b)\). The ordered sum visits each ordered index tuple once. Distinct daughters yield six tuples; a repeated daughter \(d\) yields three. Consequently, the collected coefficients multiplying the indicated normal-ordered monomials are
\[
H_{\rm decay}=g_{pab}a_pa_a^\dagger a_b^\dagger+{\rm h.c.},
\quad g_{pab}=T_{p\bar a\bar b},
\]
\[
H_{\rm repeat}=g_{pdd}a_p(a_d^\dagger)^2+{\rm h.c.},
\quad g_{pdd}={1\over2}T_{p\bar d\bar d}.
\]
Positive exact resonance excludes the parent being either daughter or its reciprocal. The numerical examples additionally use daughter modes whose reciprocal labels are distinct.

Reaction reversal gives the adjoint monomial on the same mode labels. Wavevector reversal gives \(\bar p\to\bar a+\bar b\), with conjugated coefficient in this gauge. These are different operations. The two all-incoming reciprocal tensor orbits already generate the forward, reverse, reciprocal-forward, and reciprocal-reverse monomials on expansion; adding them again duplicates terms. Self-reciprocal orbits must be enumerated only once, but were **not** tested here.

More generally, under \(H_3=\alpha\sum_{\rm ordered}VXXX\), the collected coefficients are \(6\alpha V\) and \(3\alpha V\). No absolute factor is portable between conventions without identifying \(\alpha\) and \(V\).

## 2. Independent finite Fock calculation

Direct ladder algebra predicts forward amplitudes
\[
M_{\rm distinct}=g_{pab}\sqrt{N_p(N_a+1)(N_b+1)},\qquad
M_{\rm repeat}=g_{pdd}\sqrt{N_p(N_d+1)(N_d+2)}.
\]
Reverse amplitudes contain \(g^*\sqrt{(N_p+1)N_aN_b}\) or
\(g^*\sqrt{(N_p+1)N_d(N_d-1)}\), evaluated in the reverse initial state.

The script independently constructs sparse oscillator annihilation matrices, forms every unique ordered \(XXX\) product, and projects its matrix entries onto exact equal harmonic energies. It compares this result with the collected resonant monomials and with individual ladder amplitudes. Projection is an algebraic rotating-wave check; it introduces no irreversible dynamics.

Toy tensor value: \(T=1+2i\). Energies are arbitrary common units:

| Test | Modes and energies | Occupation cutoffs | Matrix dimension |
|---|---|---|---:|
| Distinct | \(p,a,b,\bar p,\bar a,\bar b\): \(3,1,2,3,1,2\) | 0 through 2 each | 729 |
| Repeated | \(p,d,\bar p,\bar d\): \(2,1,2,1\) | 0 through 3 each | 256 |

Eight distinct and four repeated transitions use \(N_p=1,2\) and daughter occupations 0,1. All selected ladder paths remain within the cutoffs.

| Check | Distinct | Repeated |
|---|---:|---:|
| Unique ordered tensor tuples per orbit | 6 | 3 |
| Relative Frobenius discrepancy: projected vs collected Hamiltonian | \(2.14\,10^{-17}\) | \(1.47\,10^{-16}\) |
| Maximum sampled amplitude discrepancy | \(9.93\,10^{-16}\) | \(9.93\,10^{-16}\) |
| Maximum reverse or reciprocal conjugation discrepancy | \(9.93\,10^{-16}\) | \(9.93\,10^{-16}\) |
| Absolute raw-Hamiltonian Hermiticity discrepancy | 0 | \(1.38\,10^{-14}\) |
| Energy-commutator residual after projection | 0 | 0 |
| \(|M|^2\), one parent and daughter vacuum | 5 | 2.5 |

Each forward transition changes total particle number by +1 and harmonic energy by exactly zero. Absolute Hermiticity discrepancies reflect floating-point multiplication order; matrix comparison errors are at roundoff scale. This is no assertion about long-time convergence of a truncated bosonic Hamiltonian.

**Failed counting hypothesis / negative control:** enumerating six list permutations for a repeated tuple duplicates each of its three unique index tuples. Measured amplitude ratio: 2.0000000000000004; squared-amplitude ratio: 4.000000000000002. This coherent Hamiltonian duplication differs from duplicating an already formed kinetic event, which doubles its contribution.

## 3. Thermal factorial moments and kinetic interpretation

A finite Hamiltonian does not imply an irreversible population rate. A separate weak-coupling, secular kinetic limit with spectral integration/coarse graining and a specified population closure is required.

For independent geometric occupations,
\[
P(N)=(1-r)r^N,\quad n={r\over1-r},\quad
\langle N(N-1)\rangle=2n^2,\quad
\langle(N+1)(N+2)\rangle=2(1+n)^2.
\]
Exact rational partial sums, left unnormalized with their omitted mass recorded, verify convergence:

| \(r\) | Cutoffs examined | Largest final relative moment error | Final omitted mass |
|---|---|---:|---:|
| 1/2 | 4, 8, 16, 32, 64 | \(5.82\,10^{-17}\) | \(2.71\,10^{-20}\) |
| 3/4 | 8, 16, 32, 64, 128 | \(7.34\,10^{-14}\) | \(7.64\,10^{-17}\) |

**Failed mean-only closure:** a Fock state \(N=1\) has annihilation/creation factorial moments 0 and 6. A geometric state of the same mean has 2 and 8. Thus repeated-mode rates cannot be inferred from mean occupation alone.

Let \(\mathcal G=(2\pi/\hbar)\delta(E_p-E_a-E_b)\) denote the formal golden-rule spectral factor, not a numerical evaluation of \(\delta(0)\) in a finite system. Under the stated independent geometric closure:
\[
J_{\rm distinct}=\mathcal G|g_{pab}|^2
[n_p(1+n_a)(1+n_b)-(1+n_p)n_an_b],
\]
\[
J_{\rm repeat}=2\mathcal G|g_{pdd}|^2
[n_p(1+n_d)^2-(1+n_p)n_d^2].
\]
Accordingly the unordered **population-event coefficients**, with the displayed Bose brackets factored out, are
\(\kappa_{\rm distinct}=\mathcal G|T|^2\) and
\(\kappa_{\rm repeat}=\mathcal G|T|^2/2\).
For the \(\alpha V\) convention they are \(36\mathcal G|\alpha V|^2\) and \(18\mathcal G|\alpha V|^2\).

The repeated daughter still gains **two particles per net event**, \(\dot n_d=2J\), while \(\dot n_p=-J\). This stoichiometric two is distinct from both tensor counting and the factorial moment.

At equilibrium the Bose factor is \(B=n_p(1+n_a)(1+n_b)\), or \(B=n_p(1+n_d)^2\) after the repeated factorial two has been included in \(\kappa\). With \(\psi_i=\delta n_i/[n_i(1+n_i)]\), linearization gives
\(\delta J=\kappa B(\psi_p-\psi_a-\psi_b)\), with \(2\psi_d\) for repeated daughters. Thus an equilibrium linearized event weight is \(\kappa B\), not a bare Hamiltonian coefficient or a linewidth. Exact-rational detailed-balance examples produced equal forward/reverse factors \(3/5\) (distinct) and \(8/3\) (repeated including the factorial two).

## 4. Limits and reproducibility

**Established within the convention:** ordered tuple multiplicities, Fock amplitudes, their reciprocal/adjoint conjugations for the tested disjoint orbits, and the extra geometric factorial moment. **Unresolved here:** exported tensor normalization and phases, self-reciprocal orbit handling, the continuum kinetic limit, absolute rate/THz conversions, half-linewidth mapping, and all material validation. No novelty claim is made.

If permutation symmetry is unavailable, the coherent amplitude must first be formed from the relevant ordered tensor entries; summing their squared magnitudes is generally different. Squared export data alone do not establish relative phases.

The next discriminating step is a source-level identification of the exported amplitude's Hamiltonian convention and enumeration contract. No further toy or material run is proposed by this branch.

Artifacts: [script](../experiments/C_fock_check.py), [full numerical results](../experiments/C_fock_check.json).

Executed command:
~~~powershell
py -B 'D:\ResearchLab\orchestration\campaigns\20260926-105330-aln-channel-counting\experiments\C_fock_check.py'
~~~
All script assertions passed. Portable reproduction requires Python >=3.10, NumPy, SciPy; no repository import or material files. The optional flag --output-dir PATH overrides the default script directory.

Actual interpreter: C:\Users\Koussay\AppData\Local\Python\pythoncore-3.14-64\python.exe; Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1.

Executed script SHA-256: bc706066c8111c3b45740ff8f3a31106da5fe68aba5bc3e965a06382082dc59d.

