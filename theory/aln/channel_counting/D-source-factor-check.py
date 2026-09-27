"""D independent units/counting check; source kernels and arithmetic only.
No force constants, material structures, peer outputs, or solvers are loaded.
"""
from pathlib import Path
from types import SimpleNamespace
import hashlib, importlib.metadata, itertools, json, sys
import numpy as np
from phonopy.physical_units import get_physical_units
from phono3py.phonon3.imag_self_energy import ImagSelfEnergy

ROOT=Path(__file__).resolve().parent
u=get_physical_units()
Nq=27
hbar_ps=u.Hbar*1e12
h_thz=u.THzToEv
# Literal installed expressions, alongside a distinct oscillator-length route.
cpp=(u.Hbar*u.EV)**3/36/8*u.EV**2/u.Angstrom**6/(2*np.pi*u.THz)**3/u.AMU**3/Nq/u.EV**2
cgamma=18*np.pi/(u.Hbar*u.EV)**2/(2*np.pi*u.THz)**2*u.EV**2
freq=np.array([3.,1.,2.]);mass=np.array([27.,14.,27.]);fc_numeric=1.23
p_source=cpp*fc_numeric**2/(np.prod(mass)*np.prod(freq))
oscillator_A=np.sqrt(u.Hbar*u.EV/(2*mass*u.AMU*(2*np.pi*freq*u.THz)))/u.Angstrom
v_length=fc_numeric*np.prod(oscillator_A)/(6*np.sqrt(Nq))
p_length=v_length**2


def source_gamma(frequencies, observed, pair, p_value, temperature):
    # Gamma arithmetic at Gamma-like synthetic labels; D_f is a formal unit
    # on-shell spectral-density value, not evaluation of a finite delta(0).
    f=np.asarray(frequencies,float);nb=len(f)
    pp=np.zeros((1,1,nb,nb))
    for j,k in set(itertools.permutations(pair)):
        pp[0,0,j,k]=p_value
    g=np.zeros((2,1,1,nb,nb))
    for j,k in itertools.product(range(nb),repeat=2):
        g[0,0,0,j,k]=float(f[observed]-f[j]-f[k]==0.)
        g[1,0,0,j,k]=float(f[observed]+f[j]-f[k]==0.)-float(f[observed]-f[j]+f[k]==0.)
    fake=SimpleNamespace(_temperature=temperature,_frequencies=f[None,:],
        _triplets_at_q=np.array([[0,0,0]]),_weights_at_q=np.array([1]),
        _pp_strength=pp,_g=g,_imag_self_energy=np.zeros(1),
        _cutoff_frequency=0.,_unit_conversion=cgamma)
    ImagSelfEnergy._ise_thm_with_band_indices(fake)
    return float(fake._imag_self_energy[0])

T=300.;P=1e-12;Df=1.
checks={}
for repeated in (False,True):
    frequencies=np.array([2.,1.]) if repeated else np.array([3.,1.,2.])
    n=1/np.expm1(h_thz*frequencies/(u.KB*T));var=n*(1+n)
    if repeated:
        kappa=18*P/hbar_ps**2*Df
        q=n[0]*(1+n[1])**2
        incidence=np.array([-1.,2.])
        pairs=[(1,1),(0,1)]
    else:
        kappa=36*P/hbar_ps**2*Df
        q=n[0]*(1+n[1])*(1+n[2])
        incidence=np.array([-1.,1.,1.])
        pairs=[(1,2),(0,2),(0,1)]
    gammas=np.array([source_gamma(frequencies,i,pair,P,T) for i,pair in enumerate(pairs)])
    event_diag=kappa*q*incidence**2/var
    inverse_lifetime=4*np.pi*gammas
    record={'frequencies_THz':frequencies.tolist(),'temperature_K':T,'P_eV_squared':P,
        'formal_spectral_density_per_THz':Df,'equilibrium_n':n.tolist(),
        'event_kappa_per_ps':float(kappa),'equilibrium_forward_factor_Q':float(q),
        'stoichiometry':incidence.tolist(),'source_gamma_THz':gammas.tolist(),
        'source_inverse_lifetime_per_ps':inverse_lifetime.tolist(),
        'scalar_geometric_closure_diagonal_per_ps':event_diag.tolist(),
        'event_diagonal_over_inverse_lifetime':(event_diag/inverse_lifetime).tolist()}
    if repeated:
        counts=np.arange(1001,dtype=float)
        ratio=n[1]/(1+n[1]);prob=(1-ratio)*ratio**counts
        factorial_minus=float(np.dot(prob,counts*(counts-1)))
        factorial_plus=float(np.dot(prob,(counts+1)*(counts+2)))
        record['thermal_factorial_check']={'cutoff_occupation':1000,
            'geometric_tail_probability':float(ratio**1001),
            'falling_second_moment':factorial_minus,'exact_falling_second_moment':float(2*n[1]**2),
            'raising_second_moment':factorial_plus,'exact_raising_second_moment':float(2*(1+n[1])**2),
            'maximum_relative_error':float(max(abs(factorial_minus/(2*n[1]**2)-1),abs(factorial_plus/(2*(1+n[1])**2)-1)))}
    checks['repeated' if repeated else 'distinct']=record

site=Path(sys.executable).parent.parent/'Lib/site-packages'
files=['phono3py/phonon3/interaction.py','phono3py/phonon3/imag_self_energy.py',
       'phono3py/phonon3/reciprocal_to_normal.py','phono3py/phonon3/triplets.py',
       'phono3py/conductivity/utils.py','phonopy/physical_units.py']
result={'scope':'Independent source/dimensional first pass; no material run, no peers',
 'versions':{p:importlib.metadata.version(p) for p in ('phono3py','phonopy','numpy')},
 'interpreter':sys.executable,'constants':{'hbar_eV_ps':hbar_ps,'h_THz_eV_per_THz':h_thz,'N_q':Nq},
 'unit_checks':{'pp_source_conversion_N27':float(cpp),'gamma_source_conversion':float(cgamma),
  'gamma_simplified_conversion':float(18*np.pi/h_thz**2),
  'source_pp_example_eV_squared':float(p_source),'oscillator_length_pp_example_eV_squared':float(p_length),
  'pp_independent_relative_difference':float(abs(p_source/p_length-1)),
  'distinct_rate_8pi_Cgamma':float(8*np.pi*cgamma),'distinct_rate_36_over_hbar_ps_squared':float(36/hbar_ps**2),
  'repeated_rate_4pi_Cgamma':float(4*np.pi*cgamma),'repeated_rate_18_over_hbar_ps_squared':float(18/hbar_ps**2)},
 'conditional_channel_checks':checks,
 'source_hashes':{f:hashlib.sha256((site/f).read_bytes()).hexdigest() for f in files},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'D-source-factor-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
