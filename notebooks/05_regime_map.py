"""Schematic collision ordering, with paired illustrative AlN mode readings.

Source: Ma, Li and Luo (2014), Fig. 7; restricted acoustic slice of the
angular-frequency axis marked omega (THz), 2 to 10. Readings are approximate
(factor about two), not digitized data and not a full-spectrum average.
Fourier resonance: Kovacs (2018), arXiv:1804.05225.
The formal historical-GK x=1.8 line does not satisfy tau_N << tau_R.
The conserving closure has a different conversion; see note 14.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

BLIND_X = 1.8
ALN_VELOCITY = 6000.
ALN_TAU_N = (4e-10, 1e-8)
ALN_X_RANGE = (10., 200.)
# branch, source-axis coordinate, tau_U seconds, tau_N seconds
MODES = [("LA", 2, 7e-7, 1e-8), ("LA", 10, 1e-8, 4e-10),
         ("TA", 2, 4e-7, 2e-9), ("TA", 10, 4e-9, 4e-10)]
ALN_THICKNESSES = [(500e-9,"500 nm","#b03030"),(5e-6,"5 µm","#8a6a20"),
                   (20e-6,"20 µm","#1f6f8b"),(500e-6,"500 µm","#20591f")]

def aln_y_range(thickness, velocity=ALN_VELOCITY, tau_n=ALN_TAU_N):
    return thickness/velocity/tau_n[1], thickness/velocity/tau_n[0]

def regime(x,y):
    """Order labels only; candidate hydrodynamics needs strong separation."""
    if x <= 0 or y <= 0:
        raise ValueError("positive collision-time ratios required")
    if y < min(1,x):
        return "boundary-dominated"
    if x < min(1,y):
        return "resistive-dominated"
    if 1 < y < x:
        return "hydrodynamic candidate"
    if y > x > 1:
        return "resistive before boundary"
    return "crossover"

def window_width_decades(x):
    """Width of weak ordering 1<y<x, not a validated hydrodynamic window."""
    return np.log10(np.maximum(x,1.))

def trajectory(temperatures,tau_n,tau_r,thickness,velocity):
    """Requires supplied, justified effective collision laws."""
    t=np.asarray(temperatures,dtype=float)
    tn=np.array([tau_n(ti) for ti in t])
    tr=np.array([tau_r(ti) for ti in t])
    if np.any(tn<=0) or np.any(tr<=0) or thickness<=0 or velocity<=0:
        raise ValueError("positive times, thickness and velocity required")
    return tr/tn,thickness/velocity/tn

def main():
    fig,ax=plt.subplots(figsize=(7.4,6.4))
    xs=np.logspace(-1,5,400)
    ax.fill_between(xs,1e-3,np.minimum(xs,1),color="#c8d8e4",alpha=.75)
    ax.fill_between(xs,np.minimum(xs,1),1e4,color="#e8dcc8",alpha=.75)
    wx=np.logspace(0,5,300)
    ax.fill_between(wx,1,wx,color="#c9e0c9",alpha=.85)
    ax.plot(xs,xs,color=".3",lw=1.3)
    ax.axhline(1,color=".3",lw=1.3)
    ax.axvline(BLIND_X,color="#b03030",lw=1.5,ls="--")
    ax.annotate("x = 1,8 : GK historique extrapolé\nsans séparation forte",
                (BLIND_X,2e-3),xytext=(6,3),textcoords="offset points",
                fontsize=8,color="#b03030")
    ax.text(.3,2e3,"résistif",fontsize=11,color="#8a6a20",ha="center")
    ax.text(2e3,.035,"frontières\ndominantes",fontsize=10,color="#2a5a80",ha="center")
    ax.text(2e3,15,"candidat hydrodynamique\nexiger 1 ≪ y ≪ x",
            fontsize=9,color="#20591f",ha="center",rotation=24)
    for d,label,color in ALN_THICKNESSES:
        for branch,marker in (("LA","o"),("TA","s")):
            modes=[m for m in MODES if m[0]==branch]
            x=np.array([m[2]/m[3] for m in modes])
            y=np.array([d/ALN_VELOCITY/m[3] for m in modes])
            ax.plot(x,y,marker=marker,color=color,lw=1.4,ms=5,
                    label=label if branch=="LA" else None)
        print(label,[(m[0],m[1],regime(m[2]/m[3],d/ALN_VELOCITY/m[3])) for m in MODES])
    ax.set(xscale="log",yscale="log",xlim=(.1,1e5),ylim=(1e-3,1e4),
           xlabel=r"$x=\tau_R/\tau_N$",ylabel=r"$y=\tau_B/\tau_N=d/(v\tau_N)$")
    ax.set_title("Hiérarchies de collision : points AlN illustratifs à 300 K",fontsize=11)
    ax.legend(title="Épaisseur ; ○ LA, □ TA",fontsize=8,title_fontsize=8,loc="upper right")
    ax.grid(alpha=.2,which="both")
    fig.text(.5,.012,"Lectures appariées, facteur ≈ 2 ; segments guides ; spectre incomplet",
             ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.03,1,1))
    fig.savefig(Path(__file__).resolve().parents[1]/"figures/05_regime_map.png",dpi=180)
    plt.show()

if __name__=="__main__":
    main()
