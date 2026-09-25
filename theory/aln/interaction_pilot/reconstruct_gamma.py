"""Reconstruct gamma from pp without invoking a collision accumulation kernel."""
from pathlib import Path
import json,h5py,numpy as np
from phonopy.physical_units import get_physical_units
ROOT=Path(__file__).resolve().parent
RUN=ROOT/"runs/gp1-m333"
def read(name):
 with h5py.File(RUN/name) as f:return {k:f[k][()] for k in f}
ph=read("phonon-m333.hdf5");pp=read("pp-m333-g1-s0.1.hdf5")
ref=read("kappa-m333-g1-s0.1.hdf5")
u=get_physical_units();temp=float(ref["temperature"][0]);sigma=float(ref["sigma"])
freq=ph["frequency"];t=pp["triplet"];p=pp["pp"];weight=pp["weight"]
f0=freq[t[:,0],:,None,None];f1=freq[t[:,1],None,:,None];f2=freq[t[:,2],None,None,:]
cutoff=1e-4
def occupation(f):
 out=np.zeros_like(f);mask=f>cutoff
 out[mask]=1/np.expm1(u.THzToEv*f[mask]/(u.KB*temp))
 return out
n1=occupation(f1);n2=occupation(f2)
gaussian=lambda d:np.exp(-0.5*(d/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
kernel=((n1+n2+1)*gaussian(f0-f1-f2)
        +(n1-n2)*(gaussian(f0+f1-f2)-gaussian(f0-f1+f2)))
kernel=np.where((f1>cutoff)&(f2>cutoff),kernel,0.)
conversion=18*np.pi/(u.Hbar*2*np.pi*u.THz)**2
actual=conversion*np.einsum("t,tijk,tijk->i",weight,p,kernel)
expected=ref["gamma"][0]
scaled=float(np.max(abs(actual-expected))/np.max(abs(expected)))
result={"method":"Direct exported-array Gaussian/Bose sum, constants only imported; no phono3py accumulation routine.",
 "scope":"Same stipulated quadrature and force constants; not independent material physics.",
 "C_gamma":conversion,"hbar_eV_s":u.Hbar,"THzToEv":u.THzToEv,"KB_eV_K":u.KB,
 "cutoff_THz":cutoff,"sigma_THz":sigma,"temperature_K":temp,
 "reconstructed_gamma_THz":actual.tolist(),"exported_gamma_THz":expected.tolist(),
 "max_abs_error_THz":float(np.max(abs(actual-expected))),"max_scaled_error":scaled,
 "incorrect_extra_division_by_Nq_scaled_error":float(np.max(abs(actual/27-expected))/np.max(abs(expected))),
 "lifetime_inverse_s":(4*np.pi*1e12*actual).tolist(),
 "checks":{"gamma_matches":scaled<1e-12,"extra_Nq_is_wrong":np.max(abs(actual/27-expected))/np.max(abs(expected))>.9}}
result["checks"]={k:bool(v) for k,v in result["checks"].items()}
(ROOT/"gamma-reconstruction.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
print(json.dumps(result,indent=2))
assert all(result["checks"].values())
