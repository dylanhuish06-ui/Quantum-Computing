"""Generate corrected 51-interaction Chapter 5 figures as PNG and SVG.

Usage from repository root:
    python scripts/generate_corrected_figures.py

This is a separate corrected figure series, not a replacement for the original
validated five-curve milestone figures. Needs numpy, scipy, matplotlib.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "corrected_figures"
OMEGA = 2*np.pi*1e6
DELTA_MAX = 2*np.pi*50e6
T = 1e-6
VS = np.arange(51, dtype=float)
X = np.linspace(0., 1., 1401)
TIMES = X*T
KET_00 = np.array([1., 0., 0., 0.], dtype=complex)


def hbar(delta, v_mhz):
    c = OMEGA/2
    v = 2*np.pi*1e6*v_mhz
    return np.array([[delta,c,c,0],[c,0,0,c],
                     [c,0,0,c],[0,c,c,v-delta]], dtype=complex)


def static_states(v):
    ev, vec = np.linalg.eigh(hbar(0., v))
    amp = vec.conj().T @ KET_00
    return (np.exp(-1j*np.outer(TIMES,ev))*amp) @ vec.T


def sweep_states(v):
    def derivative(x, state):
        return -1j*T*(hbar(DELTA_MAX*(2*x-1)**3,v) @ state)
    sol = solve_ivp(derivative,(0.,1.),KET_00,method="DOP853",
                    t_eval=X,rtol=2e-10,atol=3e-12,max_step=.001)
    if not sol.success:
        raise RuntimeError(f"V={v}: {sol.message}")
    return sol.y.T


def measurements(s):
    return np.abs(s[...,0])**2, np.abs((s[...,1]+s[...,2])/np.sqrt(2))**2


def physical_check(s):
    norm_error = float(np.max(abs(np.sum(abs(s)**2,axis=-1)-1)))
    symmetry_error = float(np.max(abs(s[...,1]-s[...,2])))
    if norm_error>2e-8 or symmetry_error>2e-8:
        raise ValueError(f"Unphysical output: norm={norm_error}, symmetry={symmetry_error}")
    return norm_error,symmetry_error


def draw(states, figure_id, include_detuning=False):
    plt.rcParams.update({"font.family":"DejaVu Sans","mathtext.fontset":"stix",
                         "font.size":10,"svg.fonttype":"none"})
    rows = 3 if include_detuning else 2
    fig = plt.figure(figsize=(9.2,8.1 if include_detuning else 6.7),layout="constrained")
    gs = fig.add_gridspec(rows,2,width_ratios=(30,1),hspace=.07,wspace=.07)
    axes = [fig.add_subplot(gs[i,0]) for i in range(rows)]
    cmap=matplotlib.colormaps["Reds"]
    norm=Normalize(0,50)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_subplot(gs[:,1]))
    cb.set_label(r"$V/(2\pi)$ (MHz)",rotation=270,labelpad=19)
    j=1 if include_detuning else 0
    if include_detuning:
        axes[0].plot(X,50*(2*X-1)**3,color="#222222",lw=1.9)
        axes[0].set_ylabel(r"$\Delta/(2\pi)$ (MHz)")
        axes[0].set_ylim(-54,54)
    ground,bell=measurements(states)
    for i,v in enumerate(VS):
        colour="#07599a" if v==0 else cmap(norm(v))
        width=2.15 if v==0 else 1.0
        transparency=1.0 if v==0 else .55
        for ax,value in zip(axes[j:],[ground[i],bell[i]]):
            ax.plot(X,value,color=colour,lw=width,alpha=transparency,
                    zorder=6 if v==0 else 2)
    axes[j].set_ylabel(r"$P_{00}$")
    axes[j+1].set_ylabel(r"$P_{\Psi^+}$")
    for ax in axes[j:]:
        ax.set_ylim(-.035,1.035)
        ax.set_yticks(np.arange(0,1.01,.25))
    for ax in axes:
        ax.set_xlim(0,1)
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(alpha=.15,linewidth=.6)
    for ax in axes[:-1]:ax.tick_params(labelbottom=False)
    axes[-1].set_xlabel(r"Time $t$ ($\mu$s)")
    axes[-1].plot([],[],color="#07599a",lw=2.15,label="V=0 MHz")
    axes[-1].plot([],[],color=cmap(norm(30)),lw=1.5,label="V>0 MHz; 1-MHz steps")
    axes[-1].legend(loc="upper left",fontsize=8,framealpha=.92)
    fig.suptitle(f"Figure {figure_id} | 51 interaction strengths",fontsize=13)
    name=f"rydberg_figure_{figure_id.replace('.','_')}_51_interactions"
    fig.savefig(OUT/(name+".png"),dpi=230,bbox_inches="tight")
    fig.savefig(OUT/(name+".svg"),bbox_inches="tight")
    plt.close(fig)
    return ground,bell


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    static=[];sweep=[];max_norm=0.;max_sym=0.
    for v in VS:
        a,b=static_states(v),sweep_states(v)
        for states in (a,b):
            n,s=physical_check(states)
            max_norm=max(max_norm,n)
            max_sym=max(max_sym,s)
        static.append(a)
        sweep.append(b)
    static=np.asarray(static);sweep=np.asarray(sweep)
    draw(static,"5.1")
    _, bell=draw(sweep,"5.2",include_detuning=True)
    diagnostics={
        "omega_over_2pi_MHz":1,"delta_max_over_2pi_MHz":50,"duration_us":1,
        "interaction_mhz":VS.tolist(),"basis":["00","0r","r0","rr"],
        "bell_probability":"abs((c_0r+c_r0)/sqrt(2))**2",
        "cubic_detuning":"50 MHz * (2t/T-1)**3",
        "solver":"DOP853, rtol=2e-10, atol=3e-12, max_step=0.001 in t/T",
        "max_norm_error":max_norm,"max_exchange_amplitude_error":max_sym,
        "final_P_Psi_plus":{str(v):float(bell[v,-1]) for v in (0,5,10,20,30,40,50)}
    }
    (OUT/"diagnostics.json").write_text(json.dumps(diagnostics,indent=2)+"\n")
    print(json.dumps(diagnostics,indent=2))


if __name__=="__main__":
    main()
