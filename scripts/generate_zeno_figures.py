"""Validate archived two-atom research and export the publication figure suite.

Run: python scripts/generate_zeno_figures.py --all
Only --output-dir is written. No untracked research module is required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

os.environ.setdefault('MPLCONFIGDIR', '/tmp/zeno-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
from PIL import Image, __version__ as pillow_version
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from quantum_project.hamiltonian import rydberg_hamiltonian

DEFAULT = ROOT / 'results/zeno_figures'
NAVY, TEAL, AMBER, RED, GREY = '#172B4D', '#087E8B', '#D99527', '#C0504D', '#667085'
COLORS = [NAVY, TEAL, AMBER, RED]
P = np.diag([1., 0, 0, 0]).astype(complex)
Q = np.eye(4) - P
R = np.diag([0., 0, 0, 1]).astype(complex)
CMAP = LinearSegmentedColormap.from_list('zeno', ['#F5F8FA', '#86BFC5', TEAL, NAVY])
BASE_NOTE = r'$\Omega/2\pi=1\,\mathrm{MHz}$  ·  $\Delta=\phi=0$  ·  $\rho_0=|00\rangle\langle00|$'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def spectral(h, time):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * e * time)) @ v.conj().T


def channel(u, projectors=(P, Q)):
    """Column-major vectorisation of the complete non-selective Kraus map."""
    return sum(np.kron((p @ u).conj(), p @ u) for p in projectors)


def final_density(h, total, count, projectors=(P, Q), spectral_reference=False):
    if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or count < 0:
        raise ValueError('count must be a nonnegative integer')
    if not np.isfinite(total) or total < 0:
        raise ValueError('total must be finite and nonnegative')
    unitary = spectral if spectral_reference else lambda h, t: expm(-1j * h * t)
    if count == 0:
        u = unitary(h, total)
        return u @ P @ u.conj().T
    c = channel(unitary(h, total / count), projectors)
    return (np.linalg.matrix_power(c, int(count)) @ P.ravel(order='F')).reshape(4, 4, order='F')


def spectral_trajectory(h, total, count, times):
    """Independent eigen-propagation between cached powers of the channel."""
    e, v = np.linalg.eigh(h)
    def u(t):
        return (v * np.exp(-1j * e * t)) @ v.conj().T
    if count:
        scaled = times / total * count
        nearest = np.rint(scaled)
        events = np.where(abs(scaled - nearest) < 8*np.finfo(float).eps*count,
                          nearest, np.floor(scaled)).astype(int)
        c = channel(u(total / count))
        cache = {int(k): (np.linalg.matrix_power(c, int(k)) @ P.ravel(order='F')).reshape(4, 4, order='F')
                 for k in np.unique(events)}
        return np.array([u(t-k*total/count) @ cache[k] @ u(t-k*total/count).conj().T
                         for t, k in zip(times, events)])
    return np.array([u(t) @ P @ u(t).conj().T for t in times])


def populations(rho):
    return np.stack([rho[..., 0, 0].real,
                     (rho[..., 1, 1]+rho[..., 2, 2]+rho[..., 1, 2]+rho[..., 2, 1]).real/2,
                     rho[..., 3, 3].real], axis=-1)


def validate_inputs(inputs):
    for record in json.loads((inputs / 'provenance.json').read_text()):
        path = inputs.parent / record['snapshot_path']
        if sha(path) != record['sha256']:
            raise ValueError(f'Input fingerprint mismatch: {path}')
    with np.load(inputs / 'zeno_study.npz', allow_pickle=False) as saved:
        z = {k: saved[k] for k in saved.files}
    expected = {'omega_rad_s': 2*np.pi*1e6, 'total_time_s': .5e-6,
                'delta_rad_s': 0., 'phase_rad': 0.}
    for key, value in expected.items():
        if float(z[key]) != value:
            raise ValueError(f'Unsupported configuration: {key}')
    np.testing.assert_array_equal(z['V_over_omega'], [0, 1, 10, 100])
    np.testing.assert_array_equal(z['trace_counts'], [0, 2, 10, 100, 1000])
    np.testing.assert_array_equal(z['measurement_counts'], [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
    np.testing.assert_array_equal(z['questions'], ['ground', 'any_excitation', 'double'])
    np.testing.assert_allclose(z['times_s'], np.linspace(0, .5e-6, 1001), atol=0, rtol=0)
    metrics = {}
    def check(name, actual, expected, tol=3e-11):
        error = float(np.max(abs(np.asarray(actual) - np.asarray(expected))))
        if not np.isfinite(error) or error > tol:
            raise ValueError(f'{name}: {error} exceeds {tol}')
        metrics[name] = max(metrics.get(name, 0), error)
    for name in ('baseline', 'trajectory', 'final'):
        rho = z[name + '_density']
        check(name+'_trace_error', np.trace(rho, axis1=-2, axis2=-1), 1)
        check(name+'_hermiticity_error', rho, rho.conj().swapaxes(-1, -2))
        minimum = float(np.linalg.eigvalsh(rho).min())
        if minimum < -3e-11:
            raise ValueError('Density positivity failure')
        metrics[name+'_minimum_eigenvalue'] = minimum
        check(name+'_population_reconstruction', populations(rho), z[name+'_probabilities'], 2e-14)
    omega, total = float(z['omega_rad_s']), float(z['total_time_s'])
    for i, ratio in enumerate(z['V_over_omega']):
        h = rydberg_hamiltonian(omega, 0, 0, ratio*omega)
        check('hamiltonian_hermiticity', h, h.conj().T, 0)
        # Explicit independent Hamiltonian convention, including arbitrary phase/sign.
        expected_h = omega*np.array([[0,.5,.5,0],[.5,0,0,.5],[.5,0,0,.5],[0,.5,.5,ratio]])
        check('hamiltonian_convention', h, expected_h, 0)
        variance = float(np.vdot(h[:, 0], h[:, 0]).real - h[0, 0].real**2)
        check('relative_variance_error', z['variance_rad2_s2'][i]/variance, 1, 1e-14)
        check('baseline_reference_error', z['baseline_density'][i], spectral_trajectory(h,total,0,z['times_s']))
        for j, count in enumerate(z['trace_counts']):
            check('trajectory_reference_error', z['trajectory_density'][i,j],
                  spectral_trajectory(h,total,int(count),z['times_s']))
        for j, ps in enumerate([(P,Q),(Q,P),(R,np.eye(4)-R)]):
            for k, count in enumerate(z['measurement_counts']):
                check('final_reference_error',z['final_density'][i,j,k],
                      final_density(h,total,int(count),ps,True))
        check('M0_unmeasured_error',z['final_density'][i,:,0],z['baseline_density'][i,-1])
        check('M1_population_error', z['final_probabilities'][i,:,1,0],z['baseline_probabilities'][i,-1,0])
        tau = z['short_times_s']
        exact = np.array([abs(spectral(h,t)[0,0])**2 for t in tau])
        check('short_time_reference_error',z['short_survival'][i],exact)
        check('quadratic_prediction_error',z['short_prediction'][i],1-variance*tau**2,1e-14)
        check('short_residual_error',z['short_residual'][i],z['short_survival'][i]-(1-variance*tau**2),1e-14)
        counts = z['measurement_counts'][1:]
        selective = np.array([abs(expm(-1j*h*total/m)[0,0])**(2*int(m)) for m in counts])
        check('all_success_reference_error',z['exact_all_yes'][i],selective)
        check('exponential_prediction_error',z['exponential_prediction'][i],np.exp(-variance*total**2/counts),1e-14)
        if np.min(z['final_probabilities'][i,0,1:,0]-selective) < -3e-11:
            raise ValueError('All-success probability exceeds non-selective population')
        if z['final_probabilities'][i,0,-1,0] < .999:
            raise ValueError('High-frequency limit failed')
        coeff = counts[-1]*(1-z['final_probabilities'][i,0,-1,0])
        if abs(coeff/(variance*total**2)-1) > .002:
            raise ValueError('Large-M leading loss scaling failed')
    check('outcome_reversal_error',z['final_density'][:,0],z['final_density'][:,1],0)
    for ps in [(P,Q),(R,np.eye(4)-R)]:
        check('projector_completeness',sum(ps),np.eye(4),0)
        for a in ps:
            check('projector_idempotence',a@a,a,0)
            check('projector_hermiticity',a,a.conj().T,0)
        check('projector_orthogonality',ps[0]@ps[1],0,0)
    metrics['tau_Z_seconds'] = float(1/np.sqrt(z['variance_rad2_s2'][0]))
    metrics['asymptotic_loss_coefficient'] = float(z['variance_rad2_s2'][0]*total**2)
    return z, metrics


def rk4_density(h, total, count, steps):
    """Independent integration diagnostic; not the production propagator."""
    dt = total/count/steps
    rho = P.copy()
    def f(r):
        return -1j*(h@r-r@h)
    for _ in range(count):
        for _ in range(steps):
            k1=f(rho); k2=f(rho+dt*k1/2); k3=f(rho+dt*k2/2); k4=f(rho+dt*k3)
            rho = rho+dt*(k1+2*k2+2*k3+k4)/6
        rho=P@rho@P+Q@rho@Q
    return rho


def calculate(z):
    omega, total = float(z['omega_rad_s']), float(z['total_time_s'])
    h = rydberg_hamiltonian(omega, 0, 0, 0)
    times = np.linspace(.01e-6, 1.e-6, 81)
    counts = np.unique(np.rint(np.geomspace(1,1000,81)).astype(int))
    heat = np.array([[final_density(h,t,int(m))[0,0].real for t in times] for m in counts])
    base = np.array([final_density(h,t,0)[0,0].real for t in times])
    ratios = np.r_[0,np.geomspace(.01,100,100)]
    interactions = np.array([[populations(final_density(rydberg_hamiltonian(omega,0,0,v*omega),total,m))
                             for m in [0,10,1000]] for v in ratios])
    # Short-time loss evaluated as outgoing population, avoiding 1-P cancellation.
    short_t = np.geomspace(.01e-9,100e-9,160)
    short_loss = np.array([[np.sum(abs(expm(-1j*rydberg_hamiltonian(omega,0,0,v*omega)*t)[1:,0])**2)
                           for t in short_t] for v in [0,100]])
    psi=expm(-1j*h*.15e-6)[:,0]
    before=np.outer(psi,psi.conj()); after=P@before@P+Q@before@Q
    ch=rydberg_hamiltonian(omega,0,0,omega)
    steps=np.array([1,2,4,8,16,32,64])
    reference=final_density(ch,total,8,spectral_reference=True)
    rk=np.array([np.linalg.norm(rk4_density(ch,total,8,int(n))-reference,'fro') for n in steps])
    split=[]
    for n in steps:
        u=np.linalg.matrix_power(expm(-1j*ch*total/8/n),int(n))
        r=(np.linalg.matrix_power(channel(u),8)@P.ravel(order='F')).reshape(4,4,order='F')
        split.append(np.linalg.norm(r-reference,'fro'))
    data=dict(heat_times_s=times,heat_counts=counts,heat_P00=heat,heat_unmeasured_P00=base,
              interaction_ratios=ratios,interaction_counts=np.array([0,10,1000]),interaction_populations=interactions,
              short_times_s=short_t,short_ratios=np.array([0,100]),short_loss=short_loss,
              coherence_before=before,coherence_after=after,convergence_steps=steps,
              convergence_rk4_error=rk,convergence_expm_error=np.array(split))
    if not np.all(np.isfinite(np.concatenate([a.ravel() for a in data.values()]))):
        raise ValueError('Nonfinite new results')
    order=np.log2(rk[:-2]/rk[1:-1])
    if not np.all((order>3.7)&(order<4.3)):
        raise ValueError(f'RK4 convergence order failure: {order}')
    if np.max(split)>1e-12:
        raise ValueError('Exact propagation subdivision disagreement')
    if not np.allclose(after[1:,1:],before[1:,1:],atol=0,rtol=0) or np.max(abs(after[0,1:]))!=0:
        raise ValueError('Incorrect measurement coherences')
    for name in ['heat_P00','heat_unmeasured_P00','interaction_populations','short_loss']:
        if data[name].min() < -3e-11 or data[name].max()>1+3e-11:
            raise ValueError(f'Invalid probabilities: {name}')
    return data


def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'stix',
        'font.size':10,'axes.labelsize':11,'axes.titlesize':11,'axes.titleweight':'bold',
        'text.color':NAVY,'axes.labelcolor':NAVY,'axes.edgecolor':'#AAB2BD',
        'xtick.color':GREY,'ytick.color':GREY,'axes.spines.top':False,'axes.spines.right':False,
        'axes.linewidth':.7,'lines.linewidth':2.,'legend.frameon':False,'legend.fontsize':9,
        'xtick.labelsize':9,'ytick.labelsize':9,'savefig.facecolor':'white','figure.facecolor':'white',
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'quantum-zeno-figures-v1',
        'axes.unicode_minus':True})


def frame(number,title,subtitle,size=(8.4,4.8),bottom=.16):
    fig=plt.figure(figsize=size)
    fig.text(.07,.965,f'{number:02d}',color=TEAL,fontsize=13,weight='bold',va='top')
    fig.text(.19 if size[0]<5 else .125,.965,title,fontsize=14,weight='bold',va='top')
    fig.text(.07,.897,subtitle,fontsize=9,color=GREY,va='top')
    fig.subplots_adjust(left=.10,right=.96,bottom=bottom,top=.79,wspace=.32,hspace=.42)
    return fig


def footer(fig,text):
    fig.text(.07,.025,text,fontsize=8,color=GREY,va='bottom')


def panel(ax,label,title):
    ax.set_title(f'{label}   {title}',loc='left',pad=10)


def probability_axis(ax):
    ax.set_ylim(-.015,1.035);ax.set_yticks([0,.25,.5,.75,1])


def export(fig,out,stem,number,title,keys,guide,variant='wide'):
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    from matplotlib.text import Text
    suppressed_ticks=set()
    for ax in fig.axes:
        for axis in (ax.xaxis,ax.yaxis):
            drawn=axis._update_ticks()
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if tick not in drawn:
                    suppressed_ticks.update([id(tick.label1),id(tick.label2)])
    for label in fig.findobj(Text):
        if id(label) in suppressed_ticks:
            continue
        if label.get_visible() and label.get_text():
            bounds=label.get_window_extent(renderer)
            if bounds.width and bounds.height and (bounds.x0 < -1 or bounds.y0 < -1 or bounds.x1 > fig.bbox.width+1 or bounds.y1 > fig.bbox.height+1):
                # Tick labels outside current axis limits are not rendered.
                if label.axes is not None and label.get_clip_on():
                    continue
                raise ValueError(f'Label outside canvas in {stem}: {label.get_text()!r}, {bounds}')
    dimensions={}
    for ext in ['pdf','svg','png']:
        target=out/ext/f'{stem}.{ext}'
        metadata={'Creator':'Quantum Zeno figure suite'}
        if ext=='pdf':metadata.update(CreationDate=None,ModDate=None)
        if ext=='svg':metadata.update(Date=None)
        fig.savefig(target,dpi=400,metadata=metadata)
        dimensions[ext]=str(target.relative_to(out))
    plt.close(fig)
    return dict(number=number,title=title,stem=stem,variant=variant,files=dimensions,
                source_script='scripts/generate_zeno_figures.py',dataset_keys=keys,**guide)


def make_figures(z,d,out):
    entries=[]
    def save(fig,stem,n,title,keys,meaning,finding,axes,limitations,use,params='',variant='wide'):
        entries.append(export(fig,out,stem,n,title,keys,dict(physical_meaning=meaning,main_finding=finding,
            axes=axes,limitations=limitations,suggested_poster_use=use,
            parameters=params or 'Omega/(2pi)=1 MHz; Delta=phi=0; initial |00>; T=0.5 us.'),variant))
    # 01: vector schematic; event count includes the endpoint.
    fig=frame(1,'Measurement protocol','Ideal instantaneous binary projections; outcomes discarded unless explicitly post-selected.',(8.4,5.3))
    ax=fig.add_axes([.07,.15,.88,.64]);ax.axis('off');ax.set_xlim(-.02,1.02);ax.set_ylim(0,1)
    boxes=[(.00,.83,.24,r'$\rho_0=|00\rangle\langle00|$'),(.34,.83,.29,r'$U(\tau)=e^{-iH\tau/\hbar}$'),(.73,.83,.26,r'$\mathcal{D}(\rho)=P\rho P+Q\rho Q$')]
    for x,y,w,label in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,.13,boxstyle='round,pad=0.008,rounding_size=0.02',fc='#F0F5F7',ec='#CAD6DC'))
        ax.text(x+w/2,y+.065,label,ha='center',va='center',fontsize=12)
    for x1,x2 in [(.25,.33),(.64,.72)]:ax.annotate('',(x2,.895),(x1,.895),arrowprops=dict(arrowstyle='->',color=NAVY,lw=1.5))
    for y,m,c in [(.65,0,GREY),(.47,4,AMBER),(.29,12,TEAL)]:
        ax.text(.01,y,f'$M={m}$',va='center',fontsize=11,color=c)
        ax.plot([.18,.95],[y,y],color=c,lw=1.5)
        ax.plot(.18,y,'o',color=NAVY,ms=4)
        if m:
            xx=np.linspace(.18,.95,m+1)[1:]
            ax.plot(xx,np.full(m,y),'|',ms=18,mew=2,color=c)
        ax.text(.18,y-.055,'0',ha='center',va='top',fontsize=8)
        ax.text(.95,y-.055,'$T$',ha='center',va='top',fontsize=10)
    ax.annotate('',(.18+.77/4,.55),(.18,.55),arrowprops=dict(arrowstyle='<->',color=AMBER))
    ax.text(.18+.77/8,.575,r'$\tau=T/M$',ha='center',fontsize=11,color=AMBER)
    ax.text(.01,.075,r'$P=|00\rangle\langle00|,\quad Q=I-P$',fontsize=12)
    ax.text(.54,.075,r'$p_P(T)=\mathrm{Tr}[P\rho(T)]$',fontsize=12)
    footer(fig,r'Non-selective: sum all outcomes.  All-success: $S_M=|\langle00|U(T/M)|00\rangle|^{2M}$.  Events include $T$.')
    save(fig,'fig01_measurement_protocol',1,'Measurement protocol',[],
         'Coherent intervals followed by the complete binary Lüders channel.',
         'Increasing M shortens each coherent interval without altering the initial state.',
         'Schematic time increases left to right; event spacings are uniform.',
         'Ideal mathematical projections, no apparatus, ancilla, finite pulse width or detection inefficiency.',
         'Introduce the measurement model; pair with Figure 03.',params='Illustrative M=0,4,12; tau=T/M; general T.')
    # 02: reuse original density archive.
    fig=frame(2,'Population evolution',BASE_NOTE+r'  ·  $V=0$',(8.4,6.3),bottom=.16)
    axes=fig.subplots(2,2,sharex=True,sharey=True)
    for ax,j,letter in zip(axes.flat,[0,1,2,3],'abcd'):
        for k,(label,c,ls) in enumerate(zip([r'$|00\rangle$',r'$|\Psi^+\rangle$',r'$|rr\rangle$'],[NAVY,TEAL,AMBER],['-','--',':'])):
            ax.plot(z['times_s']*1e6,z['trajectory_probabilities'][0,j,:,k],label=label,color=c,ls=ls)
        m=int(z['trace_counts'][j]);panel(ax,letter,'Unmeasured' if m==0 else f'{m} measurements')
        probability_axis(ax);ax.set_xlim(0,.5)
        if m and m<=10:
            for t in np.linspace(0,.5,m+1)[1:]:ax.axvline(t,color=GREY,lw=.6,alpha=.20,zorder=0)
    for ax in axes[-1]:ax.set_xlabel(r'Time $t$ ($\mu$s)')
    for ax in axes[:,0]:ax.set_ylabel('Population')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.53,.055),ncol=3)
    footer(fig,r'$|\Psi^+\rangle=(|0r\rangle+|r0\rangle)/\sqrt{2}$; the antisymmetric state remains dark.  Non-selective measurement.')
    save(fig,'fig02_population_evolution',2,'Population evolution',['archive:times_s','archive:trajectory_probabilities[0,0:4]'],
         'Resolved populations in the symmetric three-state sector.',
         'Measurements of the ground-state sector inhibit excitation even at V=0, where blockade is absent.',
         'Time in microseconds; populations dimensionless, common scales in all panels.',
         'The Bell population is a coherent projector expectation, not one atom’s excitation probability.',
         'Main dynamics panel.',params='Common defaults; V=0; M=0,2,10,100; 1001 time samples.')
    # 03: wide plus single column.
    ms=z['measurement_counts'][1:]; non=z['final_probabilities'][0,0,1:,0];sel=z['exact_all_yes'][0]
    for compact in [False,True]:
        fig=frame(3,'Zeno survival' if compact else 'Frequent measurements suppress transitions',
                  r'$V=0$ · $T=0.5\,\mu\mathrm{s}$ · $\Omega/2\pi=1\,\mathrm{MHz}$',
                  (4.2,4.6) if compact else (8.4,4.8),bottom=.22 if compact else .19)
        if compact:fig.subplots_adjust(left=.17,right=.97,top=.77)
        axes=[fig.add_subplot(111)] if compact else list(fig.subplots(1,2))
        ax=axes[0]
        ax.semilogx(ms,non,'o-',color=TEAL,ms=4,label=r'Final population $p_P$')
        ax.semilogx(ms,sel,'s--',color=NAVY,ms=3,label=r'All-success probability $S_M$')
        fine=np.geomspace(20,1e4,180)
        ax.semilogx(fine,np.exp(-np.pi**2/(2*fine)),':',color=AMBER,lw=2.3,label=r'$e^{-\pi^2/(2M)}$ (large $M$)')
        ax.set_xlabel('Number of measurements $M$');ax.set_ylabel('Probability at $T$');probability_axis(ax)
        ax.set_xlim(1,1e4);ax.legend(loc='lower right',fontsize=8)
        if not compact:
            panel(ax,'a','Population versus complete history')
            b=axes[1];b.loglog(ms,1-non,'o-',color=TEAL,ms=4,label=r'$1-p_P$')
            b.loglog(ms,1-sel,'s--',color=NAVY,ms=3,label=r'$1-S_M$')
            b.loglog(fine,np.pi**2/(2*fine),':',color=AMBER,lw=2.3,label=r'$\pi^2/(2M)$')
            b.set_xlabel('Number of measurements $M$');b.set_ylabel('Loss probability');b.legend(loc='lower left',fontsize=9)
            panel(b,'b','High-survival regime');b.set_xlim(1,1e4)
        footer(fig,('Non-selective population includes return histories.\n'+r'$M=0$: $p_P=0$; $\Delta=\phi=0$; initial $|00\rangle$.') if compact else
               r'Non-selective population includes return histories; $S_M$ excludes them.  $M=0$: $p_P\simeq0$.  $\Delta=\phi=0$; initial $|00\rangle$.')
        save(fig,'fig03_zeno_survival_vs_measurements'+('_single' if compact else ''),3,'Zeno survival curve',
             ['archive:measurement_counts','archive:final_probabilities[0,0]','archive:exact_all_yes[0]'],
             'Final ground population and the joint probability of every ground outcome.',
             'Both losses approach (T/tau_Z)^2/M; return histories make non-selective population larger at finite M.',
             'Logarithmic measurement count; right panel shows logarithmic loss to resolve near-unity probabilities.',
             'Exponential shown only at M>=20; asymptotic, not an exact formula for either finite-M curve. M=0 stated separately.',
             'Flagship result; single-column version for narrow spaces.',params='Common defaults; V=0; counts 1 to 10000.',variant='single' if compact else 'wide')
    # 04: avoid inaccurate subtraction at very short times.
    fig=frame(4,'The quadratic origin of Zeno suppression',BASE_NOTE,(8.4,4.8),bottom=.19)
    axes=fig.subplots(1,2)
    t=d['short_times_s'];quad=z['variance_rad2_s2'][0]*t*t
    for i,(ratio,c,ls) in enumerate(zip([0,100],[TEAL,RED],['-','--'])):
        axes[0].loglog(t*1e9,d['short_loss'][i],color=c,ls=ls,label=fr'Exact, $V/\Omega={ratio}$')
        axes[1].semilogx(t*1e9,d['short_loss'][i]/quad,color=c,ls=ls)
    axes[0].loglog(t*1e9,quad,':',color=AMBER,label=r'$(t/\tau_Z)^2$',lw=2.3)
    axes[0].set_ylabel('Initial-state survival loss $1-P(t)$');axes[0].legend(loc='upper left',fontsize=8)
    axes[1].axhline(1,color=AMBER,ls=':',lw=2.3);axes[1].set_ylabel(r'Exact loss / $(t/\tau_Z)^2$')
    axes[1].set_ylim(.90,1.006)
    for ax in axes:ax.set_xlabel('Evolution time $t$ (ns)');ax.set_xlim(.01,100)
    panel(axes[0],'a','Short-time exponent: 2');panel(axes[1],'b','Validity of the leading term')
    footer(fig,r'$\tau_Z=\sqrt{2}/|\Omega|=225.08\,\mathrm{ns}$ for every $V$.  The quadratic window depends on higher Hamiltonian moments.')
    save(fig,'fig04_short_time_quadratic_law',4,'Short-time quadratic law',
         ['derived:short_times_s','derived:short_loss','archive:variance_rad2_s2'],
         'Survival loss of the initial pure state during a single unmeasured interval.',
         'The same leading coefficient holds at V/Omega=0 and 100; finite-time corrections need not be the same.',
         'Nanoseconds on a logarithmic x axis; log-log loss at left, loss divided by the leading quadratic term at right.',
         'Static finite-dimensional Hermitian H makes P(t) even: next correction is O(t^4). tau_Z alone does not fix the useful expansion window.',
         'Explain the mechanism beneath the main survival curve.',params='Omega/(2pi)=1 MHz; Delta=phi=0; V/Omega=0,100; t=0.01 to 100 ns; initial |00>.')
    # 05: heatmap plus matched baseline; vector mesh is intentionally not rasterised.
    fig=frame(5,'Where transitions are suppressed',BASE_NOTE+r'  ·  $V=0$',(8.4,5.3),bottom=.15)
    fig.subplots_adjust(right=.88)
    gs=fig.add_gridspec(2,1,height_ratios=[4,1.2],hspace=.16)
    ax=fig.add_subplot(gs[0]);b=fig.add_subplot(gs[1],sharex=ax)
    mesh=ax.pcolormesh(d['heat_times_s']*1e6,d['heat_counts'],d['heat_P00'],shading='nearest',cmap=CMAP,vmin=0,vmax=1,rasterized=False)
    mesh.set_edgecolor('face');mesh.set_linewidth(.5)
    ax.set_yscale('log');ax.set_ylabel('Measurements $M$');ax.tick_params(labelbottom=False)
    contour=ax.contour(d['heat_times_s']*1e6,d['heat_counts'],d['heat_P00'],levels=[.95],colors=[AMBER],linewidths=1.7)
    ax.clabel(contour,fmt={.95:r'$p_P=0.95$'},fontsize=9,inline=True)
    cb=fig.colorbar(mesh,ax=[ax,b],pad=.025,fraction=.035);cb.set_label(r'Final ground population $p_P(T)$');cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
    b.plot(d['heat_times_s']*1e6,d['heat_unmeasured_P00'],color=GREY,label='$M=0$')
    b.set_ylabel('$p_P$',labelpad=10);b.set_yticks([0,1]);b.set_ylim(-.05,1.05);b.legend(loc='upper right',fontsize=8)
    b.set_xlabel(r'Total time $T$ ($\mu$s)')
    footer(fig,'The unmeasured control exposes coherent revivals. A high endpoint population alone does not establish Zeno suppression.')
    save(fig,'fig05_zeno_regime_heatmap',5,'Zeno regime heatmap',['derived:heat_times_s','derived:heat_counts','derived:heat_P00','derived:heat_unmeasured_P00'],
         'Non-selective ground-state population over total time and measurement count.',
         'Sustained protection across the time window requires more measurements as T increases.',
         'Microseconds horizontally, logarithmic M vertically, colour gives final ground population. Amber contour marks 0.95.',
         'M is a count, not a rate: rate=M/T. Low-M revivals can also yield large final populations; see the matched M=0 panel.',
         'Show the range of useful parameters.',params='Omega/(2pi)=1 MHz; Delta=phi=V=0; T=0.01 to 1 us; integer M=1 to 1000; initial |00>.')
    # 06 matched controls separate blockade from inhibition.
    fig=frame(6,'Measurement inhibition and interaction blockade',BASE_NOTE+r'  ·  $T=0.5\,\mu\mathrm{s}$',(8.4,4.8),bottom=.22)
    axes=fig.subplots(1,2)
    for j,(m,c,ls) in enumerate(zip([0,10,1000],[GREY,AMBER,TEAL],['--','-.','-'])):
        for ax,k in zip(axes,[0,2]):
            ax.plot(d['interaction_ratios'],d['interaction_populations'][:,j,k],color=c,ls=ls,label='Unmeasured' if m==0 else f'$M={m}$')
    for ax in axes:
        ax.set_xscale('symlog',linthresh=.01,linscale=.4);ax.set_xlim(0,100);probability_axis(ax)
        ax.set_xlabel(r'Interaction strength $V/\Omega$');ax.set_xticks([0,.01,.1,1,10,100]);ax.set_xticklabels(['0','0.01','0.1','1','10','100'])
    axes[0].set_ylabel(r'Ground population $P_{00}(T)$');axes[1].set_ylabel(r'Double excitation $P_{rr}(T)$')
    panel(axes[0],'a','Protection of the initial state');panel(axes[1],'b','Suppression of double excitation')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.065),ncol=3)
    footer(fig,r'Strong $V$ inhibits double excitation without freezing $|00\rangle$. Matched $M=0$ controls separate the mechanisms.')
    save(fig,'fig06_interaction_dependence',6,'Interaction dependence',['derived:interaction_ratios','derived:interaction_counts','derived:interaction_populations'],
         'Matched measured and unmeasured endpoints as interaction strength changes.',
         'Blockade suppresses |rr> while allowing collective single excitation; ground-sector measurements inhibit excitation even at V=0.',
         'V/Omega dimensionless with a linear segment from 0 to 0.01 then logarithmic; populations at fixed T.',
         'Endpoint oscillations with V are coherent dynamics, not error bars or a monotonic blockade threshold.',
         'Distinguish the two suppression mechanisms.',params='Common defaults; V/Omega=0 and 0.01 to 100; M=0,10,1000.')
    # 07 vector matrix cells, a shared absolute-value scale and explicit phases in archive.
    fig=frame(7,'What a non-selective measurement removes',BASE_NOTE+r'  ·  $V=0$ · $t=0.15\,\mu\mathrm{s}$',(8.4,4.9),bottom=.23)
    axes=fig.subplots(1,3);labels=[r'$00$',r'$0r$',r'$r0$',r'$rr$']
    values=[abs(d['coherence_before']),abs(d['coherence_after']),abs(d['coherence_before']-d['coherence_after'])]
    vmax=max(a.max() for a in values)
    for ax,a,letter,title in zip(axes,values,'abc',[r'Before: $|\rho_{ij}|$',r'After: $|\mathcal{D}(\rho)_{ij}|$',r'Removed: $|(\rho-\mathcal{D}\rho)_{ij}|$']):
        im=ax.pcolormesh(np.arange(5),np.arange(5),a,cmap=CMAP,vmin=0,vmax=vmax,edgecolors='white',linewidth=.8,rasterized=False)
        ax.set_aspect('equal');ax.set_ylim(4,0);ax.set_xticks(np.arange(4)+.5,labels);ax.set_yticks(np.arange(4)+.5,labels)
        ax.tick_params(length=0);panel(ax,letter,title);ax.title.set_fontsize(9);ax._left_title.set_fontsize(9)
        for i in range(4):
            for j in range(4):ax.text(j+.5,i+.5,f'{a[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if a[i,j]>.53*vmax else NAVY)
        ax.add_patch(Rectangle((1,1),3,3,fill=False,edgecolor=AMBER,lw=1.5))
    ca=fig.add_axes([.30,.145,.40,.025]);cb=fig.colorbar(im,cax=ca,orientation='horizontal');cb.set_label('Absolute matrix element (dimensionless)',fontsize=9);cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
    footer(fig,r'$P=|00\rangle\langle00|$, $Q=I-P$: only cross-sector coherences vanish.  The outlined $Q$ block remains intact; magnitudes shown.')
    save(fig,'fig07_density_matrix_coherence',7,'Density-matrix coherence',['derived:coherence_before','derived:coherence_after'],
         'Magnitudes of a physical evolved density matrix before/after the binary Lüders measurement, and magnitudes of removed entries.',
         'Ground-to-excited coherences vanish, while all populations and within-Q coherences survive.',
         'Row/column labels are basis indices |00>, |0r>, |r0>, |rr>; all three matrices use one colour scale.',
         'Magnitude plots do not show phases; complete complex matrices are saved. Not a measurement resolving all four basis states.',
         'Explain the measurement superoperator.',params='Omega/(2pi)=1 MHz; Delta=phi=V=0; pure |00> evolved for 0.15 us, then one binary measurement.')
    # 08 genuine independently integrated errors.
    fig=frame(8,'Numerical accuracy and convergence',BASE_NOTE+r'  ·  $V/\Omega=1$ · $M=8$ · $T=0.5\,\mu\mathrm{s}$',(8.4,4.8),bottom=.19)
    axes=fig.subplots(1,2);n=d['convergence_steps'];err=d['convergence_rk4_error']
    axes[0].loglog(n,err,'o-',color=TEAL,label='RK4 density integration',ms=4)
    axes[0].loglog(n,err[0]/n.astype(float)**4,':',color=AMBER,label=r'$n^{-4}$ reference',lw=2.2)
    axes[0].set_ylabel(r'$\|\rho-\rho_{\mathrm{ref}}\|_F$');axes[0].legend(fontsize=8)
    axes[1].plot(n,d['convergence_expm_error']/1e-14,'o-',color=NAVY,ms=4)
    axes[1].set_xscale('log',base=2);axes[1].set_ylabel(r'$\|\rho-\rho_{\mathrm{ref}}\|_F\,/\,10^{-14}$')
    axes[1].set_ylim(bottom=0)
    for ax in axes:ax.set_xlabel('Substeps per measurement interval $n$');ax.set_xticks([1,4,16,64]);ax.set_xticklabels(['1','4','16','64'])
    panel(axes[0],'a','Fourth-order integration check');panel(axes[1],'b','Exact exponential: roundoff floor')
    footer(fig,'Reference: independent spectral unitary + channel powers. Production uses exact static propagation; no integration timestep is fitted.')
    save(fig,'fig08_numerical_convergence',8,'Numerical convergence',['derived:convergence_steps','derived:convergence_rk4_error','derived:convergence_expm_error'],
         'Actual full-density Frobenius errors relative to an independently diagonalised reference.',
         'RK4 approaches fourth order; subdividing the exact exponential only changes accumulated floating-point roundoff.',
         'Substeps per coherent interval; dimensionless Frobenius error. Left log-log, right log x and linear error divided by 1e-14.',
         'RK4 is an additional diagnostic, not the solver used to generate the production figures. Curves are numerical errors, not statistical uncertainties.',
         'Validation inset or methods panel.',params='Common defaults; V/Omega=1; M=8; n=1,2,4,8,16,32,64; spectral reference.')
    # 09 ground vs double question: compare a common observable and subspace occupancy.
    fig=frame(9,'The measured subspace determines the dynamics',BASE_NOTE+r'  ·  $V=0$ · $T=0.5\,\mu\mathrm{s}$',(8.4,4.9),bottom=.23)
    axes=fig.subplots(1,2)
    ground=z['final_probabilities'][0,0,1:];double=z['final_probabilities'][0,2,1:]
    for ax in axes:ax.set_xscale('log');ax.set_xlim(1,1e4);ax.set_xlabel('Number of measurements $M$');probability_axis(ax)
    axes[0].plot(ms,ground[:,0],'o-',color=TEAL,ms=4,label=r'Measure $\{P_{00},I-P_{00}\}$')
    axes[0].plot(ms,double[:,0],'s--',color=RED,ms=3,label=r'Measure $\{P_{rr},I-P_{rr}\}$')
    axes[1].plot(ms,ground[:,0],'o-',color=TEAL,ms=4)
    axes[1].plot(ms,1-double[:,2],'s--',color=RED,ms=3)
    projected_limit=np.cos(np.pi/np.sqrt(2))**2
    axes[0].axhline(projected_limit,color=GREY,ls=':',lw=1)
    axes[0].text(35,projected_limit+.035,r'$\cos^2(\Omega T/\sqrt{2})$',fontsize=10,color=GREY)
    axes[0].set_ylabel(r'Common observable: $P_{00}(T)$');axes[1].set_ylabel('Population in initial measurement sector')
    panel(axes[0],'a','Initial-state population');panel(axes[1],'b','Subspace population')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.07),ncol=2,fontsize=9)
    footer(fig,r'Right: rank-one $P_{00}$ versus rank-three $I-P_{rr}$. Reversing binary outcome labels leaves the non-selective channel unchanged.')
    save(fig,'fig09_measurement_subspace_comparison',9,'Measurement-subspace comparison',
         ['archive:measurement_counts','archive:final_probabilities[0,0]','archive:final_probabilities[0,2]'],
         'Two genuine binary projector choices, with a common ground observable and the population of the sector containing the initial state.',
         'Protecting a multidimensional sector permits motion inside it. Avoiding |rr> does not freeze the initial ground state.',
         'Log M; left ground population for both channels, right rank-one P00 or rank-three I-Prr sector population.',
         'The original any_excitation question only reverses the ground question’s outcomes; it is not an independent dynamics. Rank-three occupancy is not all-success probability.',
         'Explain Zeno dynamics beyond rank-one freezing.',params='Common defaults; V=0; ground vs double-excitation binary projectors; M=1 to 10000.')
    return entries


def contact_sheet(entries,out):
    fig,axes=plt.subplots(5,2,figsize=(12,16))
    fig.subplots_adjust(left=.025,right=.975,bottom=.025,top=.95,hspace=.16,wspace=.07)
    fig.suptitle('Quantum Zeno effect · figure selection sheet',x=.04,ha='left',fontsize=18,color=NAVY)
    for ax,entry in zip(axes.flat,entries):
        with Image.open(out/entry['files']['png']) as thumbnail:
            thumbnail.thumbnail((1600,1000),Image.Resampling.LANCZOS)
            ax.imshow(np.asarray(thumbnail))
        ax.set_title(f"{entry['number']:02d}  {entry['title']} · {entry['variant']}",loc='left',fontsize=9)
        ax.axis('off')
    for ext in ['png','pdf']:
        fig.savefig(out/f'contact_sheet.{ext}',dpi=300,metadata={'Creator':'Quantum Zeno figure suite',**({'CreationDate':None,'ModDate':None} if ext=='pdf' else {})})
    plt.close(fig)


def check_outputs(entries,out):
    checked=[]
    for entry in entries:
        sizes={}
        for ext,rel in entry['files'].items():
            p=out/rel
            if p.stat().st_size<1000:raise ValueError(f'Empty export: {p}')
            if ext=='pdf' and not p.read_bytes().startswith(b'%PDF-'):raise ValueError('Invalid PDF')
            if ext=='svg':
                tree=ET.parse(p);ns={'svg':'http://www.w3.org/2000/svg'}
                if not tree.findall('.//svg:text',ns):raise ValueError('SVG text missing')
                if tree.findall('.//svg:image',ns):raise ValueError('Unexpected raster figure element')
            if ext=='png':
                with Image.open(p) as im:
                    im.verify()
                with Image.open(p) as im:
                    if min(im.info.get('dpi',(0,0)))<399:raise ValueError('Insufficient raster DPI')
                    sizes['pixels']=list(im.size)
            sizes[ext]={'bytes':p.stat().st_size,'sha256':sha(p)}
        entry['export_checks']=sizes
        checked.append(entry['stem'])
    for ext in ['png','pdf']:
        if not (out/f'contact_sheet.{ext}').is_file():raise ValueError('Contact sheet missing')
    return checked


def guide(entries,out,metrics):
    text=r'''# Quantum Zeno figure guide

This directory is a figure library, not a poster or presentation. There are nine numbered figures and ten export sets (Figure 03 also has a single-column version). All numbered figures are exported as vector PDF, SVG with live text, and 400-DPI PNG. The contact sheet is only a selection aid and uses raster thumbnails.

## Reproduction

```bash
python -m pip install -r results/zeno_figures/requirements.txt -e '.[test]'
python scripts/generate_zeno_figures.py --all
python -m pytest tests/test_zeno_figures.py
```

Run from a checkout with the committed inputs. `--output-dir /tmp/zeno-reproduction` provides a separate output directory; `--validate-only` audits archived inputs without generating graphics. The manual GitHub Actions workflow regenerates to its temporary directory and uploads an artifact. No remote publishing, large-N simulation, or historical-data overwrite occurs. NumPy, SciPy, Matplotlib and Pillow (a Matplotlib dependency) are used; exact versions are recorded in `reproducibility_report.json`. Same-environment exports have fixed metadata and stable SVG IDs; different library/font versions can alter binary rendering. Numerical reproducibility is checked with stated absolute tolerances.

## Model and quantities

Basis: |00>, |0r>, |r0>, |rr>. In the repository convention, Hbar = H/hbar is

```
[[Delta, c, c, 0],
 [c*,    0, 0, c],
 [c*,    0, 0, c],
 [0,    c*, c*, V-Delta]],       c = (Omega/2) exp(-i phi).
```

Omega, Delta and V are angular frequencies in rad/s; U(t)=exp(-i Hbar t), with t in seconds. All plotted experiments use Omega/(2pi)=1 MHz, Delta=phi=0 and initial |00>. Unless explicitly varied, T=0.5 us. The exchange-antisymmetric state is dark; the remaining basis is |00>, |Psi+>=(|0r>+|r0>)/sqrt(2), |rr>. The coherent Bell population must not be mistaken for an incoherent sum for arbitrary states, although they coincide in this symmetric study.

For M>0, project at kT/M for k=1,...,M, including T; M=0 is unmeasured. The ground binary question uses P=diag(1,0,0,0), Q=I-P. Its non-selective channel is D(rho)=P rho P+Q rho Q. It preserves both diagonal blocks, including coherences inside Q. Switching the names/order of P and Q changes no non-selective dynamics. The double-excitation question instead uses R=|rr><rr| and I-R.

- **Non-selective final population:** p_P(T)=Tr(P rho_M), summed over all measurement histories, including departures followed by returns.
- **All-success probability:** S_M=|<00|U(T/M)|00>|^(2M), the joint probability of P at every event. This expression applies to the rank-one initial-state projector.
- **Post-selected state:** after a successful rank-one projection the normalized state is |00>, with conditional ground population exactly one. This is not S_M.
- **Subspace population:** for the double question the initial sector is rank-three I-R, with population 1-Prr. This is neither P00 nor the probability of all successful I-R outcomes.

From the actual Hamiltonian, Var(Hbar)=Omega^2/2 in |00>, independent of Delta, phi and V. Thus tau_Z=hbar/sqrt(Var(H))=sqrt(2)/|Omega|=225.079079 ns. For this static finite-dimensional Hermitian model, P(t)=1-(t/tau_Z)^2+O(t^4), stronger than the general requested O(t^3) remainder, because survival is even in t. No universal finite-time validity window follows from tau_Z alone; higher moments matter. Figure 04 calculates loss from the outgoing amplitudes to avoid cancellation in 1-P at very short times.

At fixed T and Hamiltonian, S_M=[1-a/M^2+O(M^-4)]^M with a=(T/tau_Z)^2, so log S_M=-a/M+O(M^-3). Both 1-S_M and 1-p_P have leading term a/M, but the non-selective population has return contributions at O(M^-2). The exponential exp(-a/M) is a large-M approximation, not an exact finite-M non-selective formula. Here a=pi^2/2. It is drawn only for M>=20 at V=0. The short-interval expansion must hold before applying it; a large interaction introduces additional spectral scales.

Interaction blockade detunes |rr> at large V, but permits collective motion between |00> and |Psi+>. Measuring the ground sector suppresses even this excitation. Measuring the double-excitation sector at high M also leaves internal motion: P00 tends to cos^2(Omega T/sqrt(2)) at Delta=0. These ideal numerical projections are not ancilla experiments or finite-strength continuous measurements. No statistical error bars are assigned to deterministic calculations.

For context on the distinction between freezing and evolution inside a multidimensional Zeno subspace, see P. Facchi and S. Pascazio, *Quantum Zeno dynamics: mathematical and physical aspects*, J. Phys. A 41 (2008) 493001, [author manuscript](https://arxiv.org/abs/0903.3297), [DOI](https://doi.org/10.1088/1751-8113/41/49/493001). The numerical results and asymptotic coefficients above were computed directly from this repository's model.

## Research inventory and preservation

The tracked Hamiltonian is `src/quantum_project/hamiltonian.py`; the milestone spectral state solver is `src/quantum_project/static_solver.py`. Existing Zeno implementations are `src/quantum_project/density_matrix.py` and `src/quantum_project/quantum_zeno.py`; studies and audit live in `scripts/quantum_zeno.py` and `scripts/audit_quantum_zeno.py`, with `tests/test_quantum_zeno.py`. These Zeno files and saved datasets were untracked at the start of this figure task. They were read and left unchanged.

`inputs/zeno_study.npz` is a byte-identical snapshot of the pre-existing `results/data/zeno_study.npz`, not a modified historical dataset. `inputs/provenance.json` records its original path and SHA-256, alongside text snapshots of the original generator, density/measurement routines, tests and audit documents. This makes the figure commit self-contained without staging the user's existing research changes. The independent figure generator imports only the tracked Hamiltonian; all saved complex matrices are revalidated against spectral unitaries and 16x16 channel powers. Plotting reads validated archived values rather than rerunning long M=10000 loops. New heatmaps, interaction controls, coherence examples and convergence diagnostics are small four-dimensional calculations saved in `derived_data.npz`.

The separate `symmetry-scaling/` worktree has a Dicke reduction, full-versus-reduced agreement data and nonuniform-interaction leakage data (`symmetry_agreement.npz`, `symmetry_broken.npz`, `symmetry_validation.json`). Inspection found no projection-protected large-N trajectories or validated frequency experiment there. Optional Figure 10 is therefore omitted. Symmetry leakage alone is not evidence of measurement protection. No large-N computation was run.

## Validation and limits

The report records full-complex-matrix reference errors, trace/Hermiticity/positivity, M=0 and M=1, exact projector properties, outcome-reversal equivalence, all-success ordering, short-time predictions, and large-M leading loss. Saved-array checks use a maximum absolute tolerance of 3e-11, allowing accumulated spectral roundoff at M=10000; no densities are renormalized, clipped or projected onto a positive cone. New RK4 convergence is checked against a separate spectral reference. Production propagation is exact for a static Hamiltonian up to floating-point error, so it has no integration-resolution bias to extrapolate. The manifest ties each graphic to specific dataset keys and hashes; file checks verify PNG dimensions/DPI, PDF signatures, parseable SVG with text and no embedded raster elements. Human visual inspection is recorded separately in `quality_control.md` after rendering, not claimed by the automated script.

Palette: navy ground population/all-success (context-specific), teal Bell population/non-selective protection, amber double-excitation/asymptotic reference, muted red alternative interaction/projector, grey unmeasured controls. State identities are consistent across population panels; comparison figures distinguish protocols by colour and line/marker style. Legends explicitly identify each quantity. White backgrounds, STIX mathematical glyphs, vector linework, and a sequential non-rainbow map support printing. Wide figures are 8.4 inches across; the Figure 03 single variant is 4.2 inches. PDFs retain embedded fonts; SVG live text requires compatible DejaVu/STIX fonts in the downstream editor.

## Figure-by-figure interpretation
'''
    for e in entries:
        text+=f"\n### {e['number']:02d} — {e['title']} ({e['variant']})\n\n"
        for key,label in [('physical_meaning','Physical meaning'),('parameters','Parameters'),('main_finding','Main finding'),('axes','Axes'),('limitations','Limitations'),('suggested_poster_use','Suggested poster use')]:
            text+=f"**{label}:** {e[key]}\n\n"
        text+=f"**Source script:** `{e['source_script']}`. **Dataset keys:** "+(', '.join(f'`{k}`' for k in e['dataset_keys']) or 'Vector schematic derived from the stated protocol; no numerical dataset.')+'\n\n'
        text+='**Exports:** '+', '.join(f'[{ext.upper()}]({path})' for ext,path in e['files'].items())+'.\n'
    text+='\nDataset prefix `archive:` refers to `inputs/zeno_study.npz`; `derived:` refers to `derived_data.npz`. Python-style array indices identify the exact subset. Figure numbering intentionally stops at 09 because the optional large-N protection study is unavailable.\n'
    (out/'figure_guide.md').write_text(text)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='Generate the complete suite (default).')
    parser.add_argument('--validate-only',action='store_true')
    parser.add_argument('--input-dir',type=Path,default=DEFAULT/'inputs')
    parser.add_argument('--output-dir',type=Path,default=DEFAULT)
    args=parser.parse_args()
    z,metrics=validate_inputs(args.input_dir.resolve())
    print('Archived inputs independently validated.',flush=True)
    if args.validate_only:
        print(json.dumps(metrics,indent=2));return
    out=args.output_dir.resolve()
    # Refuse historical result directories and the input directory itself.
    if out in [ROOT,ROOT/'results',ROOT/'results/data',ROOT/'results/figures',args.input_dir.resolve()]:
        parser.error('Choose a dedicated figure output directory.')
    for ext in ['pdf','svg','png']:(out/ext).mkdir(parents=True,exist_ok=True)
    if out/'inputs' != args.input_dir.resolve():
        import shutil
        shutil.copytree(args.input_dir,out/'inputs',dirs_exist_ok=True)
    (out/'requirements.txt').write_text((DEFAULT/'requirements.txt').read_text())
    d=calculate(z);np.savez_compressed(out/'derived_data.npz',**d)
    print('Lightweight calculations and convergence checks passed.',flush=True)
    style();entries=make_figures(z,d,out);contact_sheet(entries,out)
    checked=check_outputs(entries,out)
    guide(entries,out,metrics)
    inputs={str(p.relative_to(out)):sha(p) for p in sorted((out/'inputs').iterdir()) if p.is_file()}
    report={'status':'passed','command':'python scripts/generate_zeno_figures.py --all',
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'pillow':pillow_version},
        'source_hashes':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'src/quantum_project/hamiltonian.py',ROOT/'src/quantum_project/validation.py']},
        'input_hashes':inputs,'validation':metrics,
        'rk4_observed_orders':np.log2(d['convergence_rk4_error'][:-2]/d['convergence_rk4_error'][1:-1]).tolist(),
        'maximum_expm_subdivision_error':float(d['convergence_expm_error'].max()),
        'derived_data_sha256':sha(out/'derived_data.npz'),'requirements_sha256':sha(out/'requirements.txt'),'checked_exports':checked,
        'visual_inspection':'See quality_control.md; automation does not assert human inspection.',
        'optional_figure_10':'Omitted: no validated large-N projection-protection dataset found.',
        'numerical_tolerance':3e-11,'dpi':400,'numbered_figures':9,'export_sets':len(entries)}
    write_json(out/'reproducibility_report.json',report)
    write_json(out/'figure_manifest.json',{'schema_version':1,'datasets':{'archive':{'path':'inputs/zeno_study.npz','sha256':sha(out/'inputs/zeno_study.npz')},
        'derived':{'path':'derived_data.npz','sha256':sha(out/'derived_data.npz')}},'figures':entries,
        'contact_sheet':{'files':['contact_sheet.pdf','contact_sheet.png'],'contents':[e['stem'] for e in entries]},
        'source_script_sha256':sha(__file__)})
    print(f'Exported {len(entries)} sets (9 figures) to {out}',flush=True)


if __name__=='__main__':
    main()
