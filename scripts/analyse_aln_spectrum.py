"""Reproduce spectral AlN diagnostics; no network and no extrapolated lifetimes."""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.spectral import mode_heat_capacity, rta_moments, inplane_suppression

DATA = ROOT/"theory/aln"


def main():
    with np.load(DATA/"rao_300K_modes.npz", allow_pickle=False) as z:
        w, v, r = z["omega_rad_s"], z["velocity_m_s"], z["total_rate_s"]
        degeneracy, vol = z["degeneracy"], float(z["cell_volume_m3"])
    weight = degeneracy[:, None]/degeneracy.sum()/vol
    c = mode_heat_capacity(w, 300)*weight
    tau = np.divide(1., r, out=np.zeros_like(r), where=r > 0)
    vz = abs(v[:, :, 2])
    v2 = {"basal": (v[:,:,0]**2+v[:,:,1]**2)/2, "c_axis": v[:,:,2]**2}
    m = {key: rta_moments(c, val, r) for key,val in v2.items()}
    lamz = vz*tau
    manifest = json.loads((DATA/"sources.json").read_text(encoding="utf-8"))
    refs = manifest["rao"]["references"]
    np.testing.assert_allclose(c.sum(), refs["C_J_m3K"], rtol=1e-6)
    for key, idx in (("basal",0),("c_axis",2)):
        np.testing.assert_allclose(m[key]["kappa"], refs["kappa_RTA_W_mK"][idx][idx], rtol=3e-5)
    result = {"temperature_K":300, "model":"Rao dataset; diagonal RTA postprocessing",
              "C_J_m3K":float(c.sum()), "cell_volume_m3":vol,
              "normal_umklapp_split_available":False,
              "minimum_positive_frequency_THz":float(w[w>0].min()/2/np.pi/1e12),
              "basal_grazing_weight_fraction_abs_vz_lt_1e_10":float(m["basal"]["contributions"][vz<1e-10].sum()/m["basal"]["kappa"]),
              "gamma_zero_modes":int(np.sum(w==0)),
              "gamma_capacity_fraction":float(c[w==0].sum()/c.sum()),
              "directions":{}, "branches":[],"thickness_diagnostics":[]}
    for key, val in m.items():
        contrib = val["contributions"]
        order = np.argsort(lamz.ravel())
        accum = np.cumsum(contrib.ravel()[order])/contrib.sum()
        quant = [float(lamz.ravel()[order][np.searchsorted(accum,p)]) for p in (.1,.5,.9)]
        result["directions"][key] = {k:value for k,value in val.items() if k!="contributions"}
        result["directions"][key]["normal_mfp_quantiles_m_10_50_90"] = quant
        result["directions"][key]["branch_4_to_12_fraction"] = float(contrib[:,3:].sum()/contrib.sum())
    for b in range(w.shape[1]):
        result["branches"].append({"index":b+1,"frequency_min_THz":float(w[:,b].min()/2/np.pi/1e12),
            "frequency_max_THz":float(w[:,b].max()/2/np.pi/1e12),
            "C_J_m3K":float(c[:,b].sum()),
            "kappa_basal_W_mK":float(m["basal"]["contributions"][:,b].sum()),
            "kappa_c_W_mK":float(m["c_axis"]["contributions"][:,b].sum())})
    for d in (30e-9,500e-9,5e-6,20e-6,500e-6):
        row={"thickness_m":d}
        for key,val in m.items():
            wt=val["contributions"]
            for cutoff in (.1,1.):
                row[key+"_weight_Knz_gt_"+str(cutoff)]=float(wt[lamz/d>cutoff].sum()/wt.sum())
        for p in (0.,.5,1.):
            row["inplane_kappa_p_"+str(p)]=float(np.sum(m["basal"]["contributions"]*inplane_suppression(lamz/d,p)))
        result["thickness_diagnostics"].append(row)
    # Temperature affects only Bose weights here. Do not reuse 300 K rates.
    temperatures=np.array([30,50,75,100,150,200,250,300,400,600,800,1000])
    capacity=[]
    for t in temperatures:
        ct=mode_heat_capacity(w,t)*weight
        capacity.append({"T_K":int(t),"C_J_m3K":float(ct.sum()),
                         "branch_4_to_12_fraction":float(ct[:,3:].sum()/ct.sum())})
    result["harmonic_capacity_trajectory"]=capacity
    source=json.loads((DATA/"olympics_temperature.json").read_text(encoding="utf-8"))
    table=np.array([row for row in source["rows"] if row[0]>=100],float)
    result["temperature_transport_scope"]="Independent published bulk kappa(T); no modal rate trajectory inferred."
    result["temperature_display_range_K"]=[100,1000]
    result["high_T_fit_exponent_RTA_c_400_1000"] = float(-np.polyfit(
        np.log(table[table[:,0]>=400,0]),np.log(table[table[:,0]>=400,3]),1)[0])
    (DATA/"spectral_results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

    fig,axes=plt.subplots(2,2,figsize=(10,7.6),layout="constrained")
    ax=axes[0,0]
    for b in range(12):
        mask=(w[:,b]>0)&(r[:,b]>0)
        ax.scatter(w[mask,b]/2/np.pi/1e12,tau[mask,b]*1e12,s=3,alpha=.4,
                   label=str(b+1),color=plt.cm.tab20(b))
    ax.set(xlabel="Fréquence phononique (THz)",ylabel="Durée RTA (ps)",yscale="log",
           title="12 indices de branche, 300 K")
    ax.legend(title="Branche",ncol=4,fontsize=7,title_fontsize=8)
    ax=axes[0,1]
    for key,label in (("basal","Plan basal"),("c_axis","Axe c")):
        wt=m[key]["contributions"].ravel()
        ix=np.argsort(lamz.ravel())
        ax.plot(lamz.ravel()[ix]*1e6,np.cumsum(wt[ix])/wt.sum(),label=label)
    ax.set(xscale="log",xlim=(1e-3,1e2),xlabel="Longueur projetée |v_z| τ (µm)",
           ylabel="Fraction cumulée de κ RTA",title="Poids de transport, pas moyenne grise")
    ax.legend()
    ax=axes[1,0]
    d=np.geomspace(30e-9,500e-6,100)
    for p in (0,.5,1):
        vals=[np.sum(m["basal"]["contributions"]*inplane_suppression(lamz/di,p))/m["basal"]["kappa"] for di in d]
        ax.plot(d*1e6,vals,label=f"p = {p}")
    ax.set(xscale="log",xlabel="Épaisseur (µm)",ylabel="κ parallèle / κ basal massif",
           title="Surfaces réfléchissantes : RTA stationnaire")
    ax.legend()
    ax=axes[1,1]
    b=np.arange(1,13)
    ax.bar(b-.18,c.sum(axis=0)/c.sum(),width=.36,label="Capacité thermique")
    ax.bar(b+.18,m["c_axis"]["contributions"].sum(axis=0)/m["c_axis"]["kappa"],width=.36,label="Conductivité axe c")
    ax.set(xticks=b,xlabel="Indice de branche",ylabel="Fraction",title="Les pondérations ne sont pas interchangeables")
    ax.legend(fontsize=8)
    for ax in axes.ravel():ax.grid(alpha=.2)
    fig.suptitle("AlN à 300 K - données Rao et al. (2025), traitement spectral RTA",fontsize=12)
    fig.savefig(ROOT/"figures/07_aln_spectral.png",dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout="constrained")
    for col,label,style in ((1,"RTA basal","--"),(3,"RTA axe c","--"),(4,"Itératif basal","-"),(6,"Itératif axe c","-")):
        axes[0].loglog(table[:,0],table[:,col],style,marker=".",label=label)
    axes[0].set(xlabel="Température (K)",ylabel="κ massif (W m⁻¹ K⁻¹)",
                title="Phonon Olympics : table publiée")
    axes[0].legend(fontsize=8)
    axes[1].plot(temperatures,[a["C_J_m3K"]/1e6 for a in capacity],marker="o",label="12 branches")
    axes[1].plot(temperatures,[a["C_J_m3K"]*(1-a["branch_4_to_12_fraction"])/1e6 for a in capacity],marker=".",label="Indices 1 à 3")
    axes[1].set(xlabel="Température (K)",ylabel="C harmonique (MJ m⁻³ K⁻¹)",
                title="Spectre Rao : poids de Bose recalculés")
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.savefig(ROOT/"figures/08_aln_temperature.png",dpi=180)
    plt.close(fig)
    print(json.dumps({k:v for k,v in result.items() if k not in ("branches","harmonic_capacity_trajectory")},indent=2))


if __name__=="__main__":
    main()
