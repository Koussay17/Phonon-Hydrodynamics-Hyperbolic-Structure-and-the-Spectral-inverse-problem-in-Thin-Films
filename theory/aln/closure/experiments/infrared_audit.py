"""Independent low-frequency and grazing-angle diagnostics.

Read-only with respect to the original repository. Writes only infrared_*
campaign artifacts. The single AlN mesh is NOT treated as a convergence study.
"""
from pathlib import Path
import json
import hashlib
import platform
import scipy
import math
import sys

import mpmath as mp
import numpy as np
from scipy.constants import Boltzmann as KB, hbar
from scipy.integrate import quad
from scipy.special import roots_legendre

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from src.spectral import inplane_suppression, mode_heat_capacity


def positive_sum(values):
    return math.fsum(np.ravel(values).tolist())


def diffuse_s(t):
    """Independent stable diffuse suppression; argument t=L/(|v_n| tau)."""
    if math.isinf(t):
        return 1.0
    if t < 1e-4:
        # Integral representation: S(t)=integral_0^1 [1-exp(-t*x)] dx.
        return t*(0.5+t*(-1/6+t*(1/24+t*(-1/120+t/720))))
    return 1 + math.expm1(-t)/t


def shell_ratio_adaptive(epsilon):
    # Splitting at mu=epsilon resolves the grazing boundary layer explicitly.
    fun = lambda mu: 1.5*(1-mu*mu)*diffuse_s(epsilon/mu) if mu else 1.5
    split = min(epsilon, 1.0)
    left, el = quad(fun, 0, split, epsabs=1e-15, epsrel=1e-11, limit=200)
    right, er = quad(fun, split, 1, epsabs=1e-15, epsrel=1e-11, limit=200)
    return left+right, el+er


def shell_ratio_mp(epsilon):
    # Exact independently transformed exponential-integral expression.
    e = mp.mpf(str(epsilon))
    return 1-mp.mpf(3)/(2*e)*(mp.mpf(1)/4-mp.expint(3,e)+mp.expint(5,e))


def analyze_data():
    with np.load(ROOT/"theory/aln/rao_300K_modes.npz", allow_pickle=False) as z:
        w, v, rate = z["omega_rad_s"], z["velocity_m_s"], z["total_rate_s"]
        deg, vol = z["degeneracy"], float(z["cell_volume_m3"])
        q = z["q_fractional"]
        rate_anh, rate_iso = z["anharmonic_rate_s"], z["isotope_rate_s"]
    freq = w/(2*np.pi*1e12)
    x = hbar*w/(KB*300)
    # Independent heat-capacity evaluation (hyperbolic-sine expression).
    heat = np.full_like(x, KB)
    np.divide(x, 2*np.sinh(x/2), out=heat, where=x>0)
    heat[x>0] = KB*heat[x>0]**2
    cap = heat*deg[:,None]/deg.sum()/vol
    cap_ref = mode_heat_capacity(w, 300)*deg[:,None]/deg.sum()/vol
    tau = np.divide(1.0,rate,out=np.zeros_like(rate),where=rate>0)
    speed = np.linalg.norm(v,axis=2)
    vz = np.abs(v[:,:,2])
    mu = np.divide(vz,speed,out=np.zeros_like(vz),where=speed>0)
    lamz = vz*tau
    velocity_sq = {"basal":(v[:,:,0]**2+v[:,:,1]**2)/2,
                   "c_axis":v[:,:,2]**2}
    result = {
        "scope":"RTA diagnostics on one 24^3 grid, not mesh convergence or identification of continuum exponents",
        "temperature_K":300,
        "q_shape":list(q.shape),
        "q_degeneracy_sum":int(deg.sum()),
        "min_positive_frequency_THz":float(freq[freq>0].min()),
        "capacity_independent_max_relative_difference":float(np.max(np.abs(cap-cap_ref)/cap_ref)),
        "total_rate_component_max_relative_residual":float(np.max(np.abs(rate-rate_anh-rate_iso)[rate>0]/rate[rate>0])),
        "directions":{},
    }
    for key,v2 in velocity_sq.items():
        a = cap*v2
        moments = [a*tau**m for m in range(4)]
        sums = [positive_sum(m) for m in moments]
        active = a>0
        direct64 = [float(np.sum(m)) for m in moments]
        row = {
            "M0":sums[0], "M1_kappa_W_mK":sums[1],
            "M2":sums[2], "M3":sums[3],
            "tau_static_ps":sums[1]/sums[0]*1e12,
            "tau_memory_ps":sums[2]/sums[1]*1e12,
            "r_min_active_s-1":float(rate[active].min()),
            "r_max_active_s-1":float(rate[active].max()),
            "rate_dynamic_range_active":float(rate[active].max()/rate[active].min()),
            "f_min_relaxation_MHz":float(rate[active].min()/(2*np.pi*1e6)),
            "tau_max_active_ps":float(tau[active].max()*1e12),
            "max_numpy_vs_fsum_relative_difference":max(abs(a-b)/a for a,b in zip(sums,direct64)),
            "cumulative_frequency":[],
            "lifetime_tail":[],
            "mode_concentration":{},
            "positive_laplace_slope":[],
            "harmonic_response":[],
        }
        for cutoff in (0.5,0.75,1.0,1.5,2.0,3.0,5.0):
            mask=(freq>0)&(freq<=cutoff)
            retained=[positive_sum(m[~mask]) for m in moments]
            row["cumulative_frequency"].append({
                "cutoff_THz":cutoff,
                "irreducible_modes":int(mask.sum()),
                "M1_fraction":positive_sum(moments[1][mask])/sums[1],
                "M2_fraction":positive_sum(moments[2][mask])/sums[2],
                "M3_fraction":positive_sum(moments[3][mask])/sums[3],
                "tau_memory_after_deleting_resolved_low_f_ps":retained[2]/retained[1]*1e12,
            })
        for cutoff_ps in (100,300,1000,3000):
            mask=tau*1e12>cutoff_ps
            row["lifetime_tail"].append({
                "tau_cutoff_ps":cutoff_ps,
                "M1_fraction":positive_sum(moments[1][mask])/sums[1],
                "M2_fraction":positive_sum(moments[2][mask])/sums[2],
            })
        for m in (1,2,3):
            vals=np.sort(moments[m].ravel())[::-1]
            row["mode_concentration"][f"M{m}"]={
                "max_irreducible_mode_fraction":float(vals[0]/sums[m]),
                "top_10_irreducible_modes_fraction":positive_sum(vals[:10])/sums[m],
                "participation_effective_irreducible_count":sums[m]**2/positive_sum(vals*vals),
            }
        for ratio in (1e-4,1e-3,1e-2,0.1,1,10):
            p=float(rate[active].min())*ratio
            # Cancellation-free secant: [kappa(0)-kappa(p)]/p.
            slope=positive_sum(moments[2]/(1+p*tau))
            row["positive_laplace_slope"].append({
                "p_over_min_rate":ratio,"p_s-1":p,
                "secant_over_M2":slope/sums[2],
                "relative_first_order_remainder_to_pM2":1-slope/sums[2],
            })
        for f_hz in (1e5,1e6,1e7,1e8,1e9):
            omega=2*np.pi*f_hz
            h=np.sum(moments[1]/(1-1j*omega*tau))/sums[1]
            onepole=1/(1-1j*omega*sums[2]/sums[1])
            row["harmonic_response"].append({
                "modulation_frequency_Hz":f_hz,
                "normalized_RTA_real":float(h.real),"normalized_RTA_imag":float(h.imag),
                "relative_onepole_memory_fit_error":float(abs(onepole-h)/abs(h)),
            })
        if key=="basal":
            row["grazing_velocity_thresholds"]=[]
            for threshold in (0,1e-14,1e-12,1e-10,1e-8,1e-6,1,10,100):
                mask=vz<=threshold
                row["grazing_velocity_thresholds"].append({
                    "abs_vz_le_m_s":threshold,"irreducible_modes":int(np.sum(mask&active)),
                    "M1_fraction":positive_sum(moments[1][mask])/sums[1],
                    "M2_fraction":positive_sum(moments[2][mask])/sums[2],
                })
            row["grazing_angle_thresholds"]=[]
            for threshold in (1e-12,1e-4,1e-3,1e-2,0.1):
                mask=mu<=threshold
                row["grazing_angle_thresholds"].append({
                    "abs_cos_velocity_angle_le":threshold,
                    "M1_fraction":positive_sum(moments[1][mask])/sums[1]})
            row["film"]=[]
            for thickness_nm in (500,100,30,10,3,1,0.1,0.01,0.001):
                length=thickness_nm*1e-9
                suppression=inplane_suppression(lamz/length,0)
                individual=np.array([diffuse_s(length/l) if l else 1 for l in lamz.ravel()]).reshape(lamz.shape)
                filmw=moments[1]*suppression
                row["film"].append({
                    "thickness_nm":thickness_nm,
                    "kappa_diffuse_W_mK":positive_sum(filmw),
                    "relative_bulk":positive_sum(filmw)/sums[1],
                    "effectively_grazing_fraction_of_film":positive_sum(filmw[vz<1e-10])/positive_sum(filmw),
                    "max_abs_suppression_difference_independent_eval":float(np.max(np.abs(suppression-individual))),
                })
        result["directions"][key]=row
    # Restrict to first three branch indices near Gamma; every fit is a resolved-
    # window descriptive regression, NOT evidence of a limiting exponent.
    result["acoustic_rate_frequency_window_fits"]=[]
    for branch in range(3):
        for upper in (1.5,2.0,3.0,4.0):
            mask=(freq[:,branch]>0)&(freq[:,branch]<=upper)&(rate[:,branch]>0)
            if np.sum(mask)<3:
                continue
            xf=np.log(freq[mask,branch]); yr=np.log(rate[mask,branch])
            # q degeneracy weights represent the available angular sampling.
            alpha, intercept=np.polyfit(xf,yr,1,w=np.sqrt(deg[mask]))
            residual=yr-(intercept+alpha*xf)
            rms=np.sqrt(np.average(residual**2,weights=deg[mask]))
            result["acoustic_rate_frequency_window_fits"].append({
                "branch_index":branch+1,"upper_frequency_THz":upper,
                "irreducible_points":int(mask.sum()),
                "fitted_rate_vs_frequency_exponent":float(alpha),
                "weighted_log_residual_RMS":float(rms)})
    return result


def benchmarks():
    mp.mp.dps=70
    result={"scope":"Independent synthetic asymptotic/quadrature verification, not AlN convergence",
            "radial_d3_alpha2":[],"dynamic_radial_quadrature":[],
            "angular_reference":[],"angular_fixed_rule_errors":[],
            "mp_dps":mp.mp.dps}
    for n in (16,32,64,128,256,512,1024):
        midpoint=(np.arange(n)+0.5)/n
        right=(np.arange(n)+1)/n
        result["radial_d3_alpha2"].append({
            "N":n,"cutoff_kappa_exact":1-1/n,"cutoff_tau_memory_exact":n,
            "midpoint_kappa":float(np.mean(np.ones(n))),
            "midpoint_tau_memory":float(np.mean(1/midpoint**2)),
            "right_endpoint_tau_memory":float(np.mean(1/right**2)),
            "midpoint_tau_memory_over_N":float(np.mean(1/midpoint**2)/n),
            "right_tau_memory_over_N":float(np.mean(1/right**2)/n),
        })
        for p in (1e-2,1e-4,1e-6,1e-8):
            exact=1-math.sqrt(p)*math.atan(1/math.sqrt(p))
            approx=float(np.mean(midpoint**2/(midpoint**2+p)))
            result["dynamic_radial_quadrature"].append({
                "N":n,"p":p,"p_over_min_rate":p/midpoint[0]**2,
                "kappa_exact":exact,"kappa_midpoint":approx,
                "relative_total_response_error":abs(approx-exact)/exact,
                "relative_dynamic_deficit_error":abs(approx-exact)/(1-exact)})
    for epsilon in (1e-1,1e-2,1e-3,1e-4,1e-6,1e-8):
        reference=float(shell_ratio_mp(epsilon))
        adaptive,estimate=shell_ratio_adaptive(epsilon)
        leading=0.75*epsilon*math.log(1/epsilon)
        next_order=0.75*epsilon*(math.log(1/epsilon)+1-float(mp.euler))
        result["angular_reference"].append({
            "epsilon_L_over_mfp":epsilon,"reference_70digit":reference,
            "adaptive_double":adaptive,"adaptive_reported_abs_error":estimate,
            "relative_adaptive_error":abs(adaptive-reference)/reference,
            "ratio_to_leading_log":reference/leading,
            "relative_two_term_asymptotic_error":abs(reference-next_order)/reference})
        for n in (16,32,64,128,256,512):
            nodes,weights=roots_legendre(n)
            mus=(nodes+1)/2
            gauss=float(np.sum(weights/2*1.5*(1-mus**2)*inplane_suppression(mus/epsilon,0)))
            equi=np.linspace(0,1,n+1)
            trapezoid=float(np.trapezoid(1.5*(1-equi**2)*inplane_suppression(equi/epsilon,0),equi))
            result["angular_fixed_rule_errors"].append({
                "epsilon_L_over_mfp":epsilon,"N":n,
                "gauss_relative_error":(gauss-reference)/reference,
                "trapezoid_relative_error":(trapezoid-reference)/reference,
                "trapezoid_zero_node_floor":0.75/n,
                "gauss_minimum_mu":float(mus.min())})
    return result


def main():
    data,checks=analyze_data(),benchmarks()
    source=ROOT/"theory/aln/rao_300K_modes.npz"
    provenance={"python":platform.python_version(),"numpy":np.__version__,
                "scipy":scipy.__version__,"mpmath":mp.__version__,
                "input_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
                "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    data["provenance"]=provenance
    checks["provenance"]=provenance
    (OUT/"infrared_diagnostics.json").write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    (OUT/"infrared_benchmarks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"minimum_frequency_THz":data["min_positive_frequency_THz"],
                      "directions":{k:{p:v[p] for p in ("M1_kappa_W_mK","tau_memory_ps","r_min_active_s-1","tau_max_active_ps")}
                                    for k,v in data["directions"].items()},
                      "max_angular_independent_relative_error":max(v["relative_adaptive_error"] for v in checks["angular_reference"])},indent=2))


if __name__=="__main__":
    main()
