"""Bounded exact-rational attack of B2 sectors; no material inputs."""
from fractions import Fraction as F

p, a, c = F(1,2), F(1), F(1)
states = [(1,0),(0,2)]
fluxes = [c*(m*(n+1)*(n+2)-(m+1)*n*(n-1)) for m,n in states]
exact_flux = sum(fluxes)/2
geometric_flux = 2*c*(p*(1+2*a)-a*a)
assert exact_flux == 0 and geometric_flux == 1

# Bose equilibrium: parent ratio 1/4, daughter ratio 1/2.
# Enumerate complete energy sectors, without clipping any within-sector state.
# phi=m-1/3-(4/9)(n-1) has zero mean and zero covariance with mean energy.
mass = mean = covariance = variance = plateau = F(0)
for ell in range(121):
    pairs = [(m,ell-2*m) for m in range(ell//2+1)]
    prob = F(3,8)*F(1,2)**ell
    values = [F(m)-F(1,3)-F(4,9)*(F(n)-1) for m,n in pairs]
    conditional = sum(values)/len(values)
    mass += prob*len(values)
    mean += prob*sum(values)
    covariance += prob*sum(values)*(F(ell)-F(5,3))
    variance += prob*sum(v*v for v in values)
    plateau += prob*len(values)*conditional**2
assert abs(variance-F(68,81)) < F(1,10**28)
assert abs(plateau-F(34,729)) < F(1,10**28)
assert abs(covariance) < F(1,10**28)
assert abs(plateau/variance-F(1,18)) < F(1,10**28)
print("stationary_sector_exact_flux =", exact_flux)
print("same_means_geometric_flux =", geometric_flux)
print("equilibrium_variance_limit = 68/81")
print("exact_sector_plateau_limit = 34/729")
print("normalized_plateau_limit = 1/18")
print("enumerated_energy_sectors = 0..120")
print("omitted_probability =", float(1-mass))
print("variance_abs_error =", float(abs(variance-F(68,81))))
print("plateau_abs_error =", float(abs(plateau-F(34,729))))
print("energy_covariance_abs =", float(abs(covariance)))

