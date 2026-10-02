# C: low-dimensional resonance-measure benchmark

Status: independent benchmark definition saved; script and numerical checks pending. Read only this campaign's question and assumptions. No peer reports or material calculation.

## Definition and independently derived references

Use a one-dimensional periodic Brillouin zone q in [-pi,pi), normalized measure dq/(2pi), fixed parent momentum zero, and daughter momenta q and -q. In dimensionless energy units,
E_p=2+d, E_a(q)=1+0.3 cos(q), E_b(k)=1+0.7 cos(k).
Thus all energies stay positive in the cases below, momentum is exact, and
Delta(q)=E_p-E_a(q)-E_b(-q)=d-cos(q).

Choose positive interaction/test weight w(q)=1+0.2 cos(2q)+0.1 sin(q), bounded below by 0.7, and entropy-affinity test phi(q)=1+0.4 sin(q). Occupations, when tested, use exact analytic energies at quadrature nodes and beta=1. Scalar weak-coupling population kinetics is assumed, not derived from a finite Hamiltonian. The weight need not obey reciprocal equality; the even-weight reciprocal completion gives the same scalar mass and thermal diagnostics here.

For |d|<1 the roots are +/-acos(d), with slope magnitude s=sqrt(1-d^2). Write W=0.8+0.4d^2. A direct delta-function change of variables gives
I[1]=W/(pi*s),
I[phi^2]=[W(1+0.16s^2)+0.08s^2]/(pi*s).
At the regular test d=0.6, s=0.8 and I[1]=1.18/pi. The root weights are w(q_root)/(2pi*|sin(q_root)|), each positive.

A negative control omits 1/|sin(q_root)|. It is still positive and exactly on shell, but underestimates both regular targets by 20%. Conservation and equilibrium stationarity alone therefore cannot certify this integration rule.

Compare root quadrature with the normalized Gaussian delta_sigma(Delta)=exp[-Delta^2/(2sigma^2)]/(sqrt(2pi)*sigma), sampled on shifted uniform periodic meshes. An independent adaptive angular quadrature supplies the continuum broadened reference, separating finite-sigma error from mesh error.

## Planned diagnostics and failure cases

Measure total mass, weak-form value, energy weak-form defect integral[w delta_sigma Delta^2], equilibrium net event flux, absolute event flux, and physical energy drift. At exact Bose occupations, F=n_p(1+n_a)(1+n_b)-(1+n_p)n_a n_b has sign opposite Delta, so the broadened energy drift -integral[w delta_sigma Delta F] is nonnegative and generally positive. A symmetric positive weak form does not establish nonlinear detailed balance.

Regular d=0.6: vary mesh and width independently and record sigma/(h*s). Near-critical d=1-1e-4: expose missed pairs of roots in a sign-change scan on a coarse shifted mesh. Empty d=1.1: exact measure is zero but Gaussian leakage is positive.

At critical d=1, Delta~q^2/2 and no finite regular-root weight exists. The Gaussian mass diverges as C/sqrt(sigma), where
C=1.2 Gamma(1/4)/(2^(7/4) pi^(3/2)).
This follows by q=sqrt(sigma)*t in the Gaussian integral, not by fitting. The energy-defect ratio I[Delta^2]/(sigma^2 I[1]) tends to 1/2.

For d=0 on a fixed 64-point mesh, an aligned mesh has exact roots and its Gaussian mass grows as 1/sigma; a shifted mesh without exact roots tends to zero. This tests the noncommuting mesh/width limits.

No claims about material convergence, critical-surface regularity in higher dimension, or an occupation-interpolation scheme are made. Owned artifacts: this report and experiments/C_*; script will support --output-dir.

