"""Independent physics and artifact contracts for the self-contained figure suite."""
import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('zeno_figures', ROOT/'scripts/generate_zeno_figures.py')
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)


@pytest.mark.parametrize('ratio', [0, 1, 100, -10])
def test_channel_matches_explicit_outcome_branch_sum(ratio):
    h=f.rydberg_hamiltonian(2*np.pi*1e6, -.3e6, .37, ratio*2*np.pi*1e6)
    t=.5e-6; m=3; u=expm(-1j*h*t/m)
    branches=[f.P.copy()]
    for _ in range(m):
        branches=[p@u@r@u.conj().T@p for r in branches for p in [f.P,f.Q]]
    rho=f.final_density(h,t,m)
    np.testing.assert_allclose(rho,sum(branches),atol=2e-14,rtol=0)
    assert abs(np.trace(rho)-1)<2e-14
    assert np.linalg.eigvalsh(rho).min()>-2e-14
    # Sum of Kraus maps is CPTP even for a coherent mixed input.
    psi=np.array([1,1j,2,-3j])/np.sqrt(15)
    mixed=.4*np.outer(psi,psi.conj())+.6*np.eye(4)/4
    output=(f.channel(u)@mixed.ravel(order='F')).reshape(4,4,order='F')
    assert abs(np.trace(output)-1)<2e-14
    assert np.linalg.eigvalsh(output).min()>0


@pytest.mark.parametrize('m',[0,1,2,10,1000])
def test_unmeasured_limit_and_all_success_analytical_v0(m):
    omega=2*np.pi*1e6;t=.5e-6
    h=f.rydberg_hamiltonian(omega,0,0,0)
    rho=f.final_density(h,t,m)
    if m==0:
        assert rho[0,0].real==pytest.approx(np.cos(omega*t/2)**4,abs=1e-14)
    else:
        selective=np.cos(omega*t/(2*m))**(4*m)
        exact=abs(expm(-1j*h*t/m)[0,0])**(2*m)
        assert exact==pytest.approx(selective,abs=2e-12)
        assert rho[0,0].real>=selective-2e-12
        if m==1:assert rho[0,0].real<1e-14


def test_rank_three_protection_permits_internal_motion():
    omega=2*np.pi*1e6;t=.5e-6
    rho=f.final_density(f.rydberg_hamiltonian(omega,0,0,0),t,10000,(f.R,np.eye(4)-f.R))
    assert rho[3,3].real<.001
    assert rho[0,0].real==pytest.approx(np.cos(omega*t/np.sqrt(2))**2,abs=.001)
    assert rho[0,0].real<.5


def test_short_time_coefficient_from_actual_hamiltonian():
    omega=2*np.pi*1e6
    for v in [0,100,-100]:
        h=f.rydberg_hamiltonian(omega,-.2*omega,.7,v*omega)
        variance=np.vdot(h[:,0],h[:,0]).real-h[0,0].real**2
        assert variance==pytest.approx(omega**2/2,rel=2e-15)
        times=np.array([.8,.4,.2,.1])*1e-9
        losses=np.array([sum(abs(expm(-1j*h*t)[1:,0])**2) for t in times])
        residual=abs(losses-variance*times**2)
        orders=np.log2(residual[:-1]/residual[1:])
        assert np.all((orders>3.9)&(orders<4.1))


@pytest.mark.parametrize('m',[-1,1.5,True])
def test_reject_invalid_measurement_counts(m):
    with pytest.raises(ValueError):f.final_density(np.zeros((4,4)),1,m)


def test_convergence_and_coherence_calculations():
    with np.load(f.DEFAULT/'inputs/zeno_study.npz',allow_pickle=False) as z:
        d=f.calculate(z)
    b,a=d['coherence_before'],d['coherence_after']
    np.testing.assert_array_equal(a[1:,1:],b[1:,1:])
    np.testing.assert_array_equal(np.diag(a),np.diag(b))
    assert np.max(abs(a[0,1:]))==0
    assert abs(a[1,2])>.1
    assert abs(np.trace(a)-1)<1e-14
    assert np.linalg.eigvalsh(a).min()>-1e-14
    assert d['convergence_rk4_error'][-1]<1e-10
    assert d['convergence_expm_error'].max()<1e-12


def test_committed_exports_and_numerical_provenance():
    base=f.DEFAULT
    manifest=json.loads((base/'figure_manifest.json').read_text())
    assert manifest['source_script_sha256']==f.sha(ROOT/'scripts/generate_zeno_figures.py')
    entries=manifest['figures']
    assert len(entries)==10
    assert {e['number'] for e in entries}==set(range(1,10))
    assert manifest['contact_sheet']['contents']==[e['stem'] for e in entries]
    assert len(f.check_outputs(copy.deepcopy(entries),base))==10
    for record in manifest['datasets'].values():
        assert f.sha(base/record['path'])==record['sha256']
    for e in manifest['figures']:
        for ext,path in e['files'].items():
            assert f.sha(base/path)==e['export_checks'][ext]['sha256']
    for record in json.loads((base/'inputs/provenance.json').read_text()):
        assert f.sha(base/record['snapshot_path'])==record['sha256']


def test_reject_corrupted_archived_input(tmp_path):
    inputs=tmp_path/'inputs';inputs.mkdir()
    records=json.loads((f.DEFAULT/'inputs/provenance.json').read_text())
    record=next(r for r in records if r['snapshot_path']=='inputs/zeno_study.npz')
    (inputs/'provenance.json').write_text(json.dumps([record]))
    (inputs/'zeno_study.npz').write_bytes(b'corrupted research data')
    with pytest.raises(ValueError,match='Input fingerprint mismatch'):
        f.validate_inputs(inputs)
