"""Reproducible SYNTHETIC FDTR design study; never substitutes for laboratory data.

Run: python scripts/analyse_experiment.py
Outputs: theory/experimental_results.json and figures/06_experimental_design.png
"""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from fdtr import Medium, Film, gaussian_response, plane_response
from inversion import information_from_jacobian, NonIdentifiableError

NAMES=["film_lam","R_film_substrate","R_metal_film","metal_thickness",
       "film_thickness","film_rho_c","spot_scale","log_gain","phase_deg","delay_ns"]
# Coordinates 0:7 are log-ratios to nominal; calibration coordinates are additive.
PRIOR_SD={2:.20,3:.05,4:.03,5:.05,6:.03,7:.02,8:.10,9:.02}
SIGMA_LOG_AMPLITUDE=.01
SIGMA_PHASE=np.deg2rad(.1)
FREQ=np.geomspace(1e4,2e8,60)

def model(x,lam=60.,plane=False):
    films=[
        Film(Medium(150.,2.49e6),80e-9*np.exp(x[3]),1e-8*np.exp(x[2])),
        Film(Medium(lam*np.exp(x[0]),2.41e6*np.exp(x[5])),
             500e-9*np.exp(x[4]),2e-8*np.exp(x[1]))]
    substrate=Medium(35.,3.03e6)
    radii=np.exp(x[6])*np.array([10e-6,8e-6])
    fun=plane_response if plane else gaussian_response
    z=fun(FREQ,films,substrate,*radii)
    return z*np.exp(x[7]+1j*(np.deg2rad(x[8])-2*np.pi*FREQ*x[9]*1e-9))

def residual(z,observed):
    ratio=z/observed
    return np.r_[np.log(np.abs(ratio))/SIGMA_LOG_AMPLITUDE,
                  np.angle(ratio)/SIGMA_PHASE]

def jacobian(x,lam=60.,step=2e-4):
    baseline=model(x,lam)
    cols=[]
    for i in range(len(x)):
        dx=np.zeros_like(x);dx[i]=step
        cols.append((residual(model(x+dx,lam),baseline)
                    -residual(model(x-dx,lam),baseline))/(2*step))
    return np.array(cols).T

def design(jac,indices,priors=False,phase_only=False):
    data=jac[len(FREQ):] if phase_only else jac
    j=data[:,indices]
    augmented=j
    if priors:
        prior=np.zeros((len(PRIOR_SD),len(indices)))
        for row,(idx,sd) in enumerate(PRIOR_SD.items()):
            if idx in indices: prior[row,indices.index(idx)]=1/sd
        augmented=np.vstack((j,prior))
    norms=np.linalg.norm(augmented,axis=0)
    scaled=augmented/np.where(norms>0,norms,1)
    singular=np.linalg.svd(scaled,compute_uv=False)
    result={"parameters":[NAMES[i] for i in indices],
            "phase_only":phase_only,"calibration_priors":priors,
            "column_scaled_singular_values":singular.tolist()}
    try:
        _,cov=information_from_jacobian(augmented)
        sd=np.sqrt(np.diag(cov))
        result.update(identifiable=True,coordinate_sd=sd.tolist(),
                      correlation=(cov/np.outer(sd,sd)).tolist(),
                      film_lam_relative_sd=float(sd[indices.index(0)]))
    except NonIdentifiableError:
        result.update(identifiable=False,film_lam_relative_sd=None)
    return result

def run():
    x=np.zeros(len(NAMES))
    scenarios=[]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for lam,color in [(60.,"#245d89"),(321.,"#a04926")]:
        z=model(x,lam)
        plane=model(x,lam,True)
        j=jacobian(x,lam)
        jfine=jacobian(x,lam,1e-4)
        step_error=np.linalg.norm(j-jfine)/np.linalg.norm(jfine)
        cases={
          "thermal_only":design(j,[0,1]),
          "thermal_only_phase":design(j,[0,1],phase_only=True),
          "all_free":design(j,list(range(10))),
          "independent_calibrations":design(j,list(range(10)),priors=True),
          "calibrations_phase_only":design(j,list(range(10)),priors=True,phase_only=True)}
        phase_error=np.rad2deg(np.angle(plane/z))
        amp_error=100*(np.abs(plane/z)-1)
        # Anisotropy is not varied here: film is assumed isotropic throughout.
        scenarios.append({"film_lam_W_mK":lam,"designs":cases,
                          "derivative_step_relative_change":float(step_error),
                          "max_1d_phase_error_deg":float(max(abs(phase_error))),
                          "max_1d_amplitude_error_percent":float(max(abs(amp_error))),
                          "phase_error_deg_at_frequencies":phase_error.tolist(),
                          "amplitude_error_percent_at_frequencies":amp_error.tolist()})
        axes[0,0].semilogx(FREQ,np.rad2deg(np.angle(z)),color=color,label=f"{lam:g} W/(m K)")
        axes[0,1].semilogx(FREQ,phase_error,color=color,label=f"{lam:g}")
        axes[1,0].semilogx(FREQ,j[len(FREQ):,0]*SIGMA_PHASE*180/np.pi,color=color,label=f"{lam:g}")
    axes[0,0].set(title="Fourier axisymetrique : phase",ylabel="Phase (degres)")
    axes[0,0].legend()
    axes[0,1].set(title="Erreur de l'approximation 1D",ylabel="Phase 1D - axisymetrique (degres)")
    axes[0,1].axhline(.1,color="gray",ls=":",label="Bruit suppose : 0,1 deg")
    axes[0,1].axhline(-.1,color="gray",ls=":")
    axes[0,1].legend()
    axes[1,0].set(title="Sensibilite a la conductivite du film",ylabel="d phase / d ln(lambda) (degres)")
    axes[1,0].legend()

    rng=np.random.default_rng(20260916)
    truth=model(x)
    observed=truth*np.exp(rng.normal(0,SIGMA_LOG_AMPLITUDE,len(FREQ))
                         +1j*rng.normal(0,SIGMA_PHASE,len(FREQ)))
    def objective(y):
        return np.r_[residual(model(y),observed),
                     [y[i]/sd for i,sd in PRIOR_SD.items()]]
    lower=np.r_[np.full(7,-np.log(3.)), -.3,-2.,-1.]
    upper=-lower
    fits=[]
    for amplitude in [0.,.08,.16]:
        initial=np.zeros(10)
        initial[:7]=rng.normal(0,amplitude,7)
        fit=least_squares(objective,initial,bounds=(lower,upper),
                           diff_step=1e-4,xtol=1e-9,ftol=1e-9,gtol=1e-7,max_nfev=180)
        fits.append(fit)
    best=min(fits,key=lambda r:r.cost)
    zfit=model(best.x)
    rdata=residual(zfit,observed)
    fit_j=jacobian(best.x)
    fit_uncertainty=design(fit_j,list(range(10)),priors=True)
    axes[1,1].semilogx(FREQ,rdata[:len(FREQ)],".",label="log amplitude")
    axes[1,1].semilogx(FREQ,rdata[len(FREQ):],".",label="phase")
    axes[1,1].axhline(0,color="black",lw=.7)
    axes[1,1].set(title="Ajustement synthetique avec calibrations",ylabel="Residus / bruit")
    axes[1,1].legend()
    for ax in axes.flat:
        ax.set_xlabel("Frequence (Hz)");ax.grid(alpha=.2)
    fig.suptitle("SCENARIOS SYNTHETIQUES - aucune donnee experimentale",fontsize=13)
    (ROOT/"figures").mkdir(exist_ok=True)
    fig.savefig(ROOT/"figures/06_experimental_design.png",dpi=180)
    plt.close(fig)
    result={
      "status":"SYNTHETIC ONLY; no laboratory validation",
      "seed":20260916,"frequency_Hz":FREQ.tolist(),
      "assumptions":{"metal_lam":150.,"metal_rho_c":2.49e6,"metal_thickness_m":80e-9,
                     "film_thickness_m":500e-9,"film_rho_c":2.41e6,
                     "R_metal_film_m2K_W":1e-8,"R_film_substrate_m2K_W":2e-8,
                     "substrate_lam":35.,"substrate_rho_c":3.03e6,
                     "pump_radius_m":10e-6,"probe_radius_m":8e-6,
                     "sigma_log_amplitude":SIGMA_LOG_AMPLITUDE,
                     "sigma_phase_deg":.1,"all_properties_isotropic":True,
                     "surface_absorption":True,"noise_independent_across_frequency":True},
      "coordinates":NAMES,
      "coordinate_definition":"0:7 log ratios to nominal; log_gain additive; phase_deg additive; delay_ns additive",
      "independent_calibration_sd":{NAMES[i]:sd for i,sd in PRIOR_SD.items()},
      "scenarios":scenarios,
      "synthetic_fit":{"optimizer_success":bool(best.success),
        "all_start_success":[bool(r.success) for r in fits],
        "start_objectives":[float(2*r.cost) for r in fits],
        "max_coordinate_spread_between_starts":float(np.max(np.ptp([r.x for r in fits],axis=0))),
        "coordinates":best.x.tolist(),
        "film_lam_W_mK":float(60*np.exp(best.x[0])),
        "R_film_substrate_m2K_W":float(2e-8*np.exp(best.x[1])),
        "data_chi_squared":float(rdata@rdata),
        "total_objective":float(2*best.cost),
        "at_bound":bool(np.any(np.minimum(best.x-lower,upper-best.x)<1e-4)),
        "local_uncertainty":fit_uncertainty,
        "limitation":"Same-model synthetic fit tests estimator operation; not model adequacy or global uniqueness."}}
    (ROOT/"theory").mkdir(exist_ok=True)
    (ROOT/"theory/experimental_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for case in scenarios:
        print("lambda",case["film_lam_W_mK"],"max 1D phase error",case["max_1d_phase_error_deg"])
        for name,d in case["designs"].items(): print(name,d["film_lam_relative_sd"])
    print("synthetic fit",result["synthetic_fit"]["film_lam_W_mK"],
          "chi2",result["synthetic_fit"]["data_chi_squared"],
          "starts",result["synthetic_fit"]["start_objectives"],flush=True)

if __name__=="__main__":
    run()
