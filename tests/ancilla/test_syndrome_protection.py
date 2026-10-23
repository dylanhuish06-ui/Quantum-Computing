from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project.ancilla import FullParameters, full_hamiltonian_reference, dicke_isometry_bulk, evolve
from quantum_project.ancilla.syndrome import (
    symmetric_projector, syndrome_unitary, measure_density,
    selective_evolution, nonselective_evolution, spin_phase_register,
)
from quantum_project.ancilla.protection import (
    spin_squared, gap_penalty, permutation_operator, symmetric_group,
    permutation_average, permutation_cycles, swap_decomposition,
)


def broken(N=4):
    return full_hamiltonian_reference(FullParameters(
        np.ones(N), np.linspace(-.7, .7, N), np.zeros(N),
        .17*(np.ones((N, N))-np.eye(N)))).toarray()


def test_oracle_unitarity_both_kraus_and_partial_trace():
    P = symmetric_projector(4)
    Q = np.eye(16)-P
    U = syndrome_unitary(P)
    np.testing.assert_allclose(U.conj().T@U, np.eye(32), atol=2e-15)
    rng = np.random.default_rng(411)
    psi = rng.normal(size=16)+1j*rng.normal(size=16)
    psi /= np.linalg.norm(psi)
    joint = (U@np.kron(psi, [1, 0])).reshape(16, 2)
    np.testing.assert_allclose(joint[:, 0], P@psi, atol=1e-15)
    np.testing.assert_allclose(joint[:, 1], Q@psi, atol=1e-15)
    rho = np.outer(psi, psi.conj())
    for outcome, A in enumerate((P, Q)):
        block = measure_density(rho, P, outcome)
        np.testing.assert_allclose(block, np.outer(joint[:, outcome], joint[:, outcome].conj()), atol=1e-15)
        np.testing.assert_allclose(np.trace(block), np.linalg.norm(A@psi)**2, atol=1e-15)
    np.testing.assert_allclose(joint@joint.conj().T, measure_density(rho, P), atol=1e-15)


def test_repeated_selective_nonselective_and_all_record_enumeration():
    H, P = broken(3), symmetric_projector(3)
    initial = np.eye(8)[0].astype(complex)
    times = np.linspace(0, 1.4, 4)
    selected = selective_evolution(lambda t: H, initial, times, P, dt=.07)
    density = nonselective_evolution(lambda t: H, np.outer(initial, initial), times, P, dt=.08)
    midpoint = nonselective_evolution(lambda t: H, np.outer(initial, initial), times, P,
                                     dt=.05, method="dense_midpoint")
    np.testing.assert_allclose(density, midpoint, atol=1e-11)
    step = expm(-1j*H*(times[1]-times[0]))
    total = np.zeros((8, 8), complex)
    for record in product((0, 1), repeat=3):
        branch = initial.copy()
        for outcome in record:
            branch = (P if outcome == 0 else np.eye(8)-P)@step@branch
        total += np.outer(branch, branch.conj())
        if record == (0, 0, 0):
            np.testing.assert_allclose(branch, selected.branch[-1], atol=1e-11)
    np.testing.assert_allclose(total, density[-1], atol=1e-11)
    np.testing.assert_allclose(np.trace(density, axis1=1, axis2=2), 1, atol=2e-12)
    assert np.linalg.eigvalsh(density[-1]).min() > -1e-12
    assert 0 < selected.success[-1] < 1
    assert 1-np.trace(P@density[-1]).real > 0


def test_zeno_limit_is_projected_not_bare_and_dt_independent():
    H, P = broken(), symmetric_projector(4)
    B = dicke_isometry_bulk(4, 0).toarray()
    initial = np.eye(16)[0].astype(complex)
    T = 2.
    target = B@expm(-1j*(B.conj().T@H@B)*T)@np.eye(5)[0]
    bare = expm(-1j*H*T)@initial
    errors = []
    for checks in (8, 24, 96):
        times = np.linspace(0, T, checks+1)
        result = selective_evolution(lambda t: H, initial, times, P, dt=.07)
        errors.append(np.linalg.norm(result.conditional[-1]-target))
        target_path = np.tile(target, (len(times), 1))
        stats = result.metrics(target_path, P)
        np.testing.assert_allclose(stats["successful_target_probability"],
                                  result.success*stats["conditional_fidelity"], atol=1e-15)
        assert stats["conditional_fidelity"][-1] > .99
        assert stats["success_probability"][-1] < 1
    assert errors[2] < errors[1] < errors[0]
    assert np.linalg.norm(result.conditional[-1]-bare) > .1
    finer = selective_evolution(lambda t: H, initial, times, P, dt=.02)
    np.testing.assert_allclose(result.branch, finer.branch, atol=1e-11)


def test_time_dependent_measurements_against_independent_dense_midpoint():
    static, P = broken(3), symmetric_projector(3)
    disorder = np.diag(np.linspace(-1, 1, 8))
    H = lambda t: static+np.sin(3*t)*disorder
    initial = np.eye(8)[0].astype(complex)
    times = np.linspace(0, .7, 8)
    ode = selective_evolution(H, initial, times, P, dt=.025)
    coarse = selective_evolution(H, initial, times, P, dt=.02, method="dense_midpoint")
    fine = selective_evolution(H, initial, times, P, dt=.01, method="dense_midpoint")
    assert np.linalg.norm(fine.branch[-1]-ode.branch[-1]) < np.linalg.norm(coarse.branch[-1]-ode.branch[-1])/3.9
    rho = np.outer(initial, initial)
    ode_rho = nonselective_evolution(H, rho, times, P, dt=.02)
    dense_rho = nonselective_evolution(H, rho, times, P, dt=.002, method="dense_midpoint")
    np.testing.assert_allclose(ode_rho, dense_rho, atol=3e-7)


def test_zero_success_is_not_a_fictitious_normalised_state():
    P = np.diag([1., 0.])
    result = selective_evolution(lambda t: np.zeros((2, 2)), [0, 1], [0, .5], P)
    assert result.success[-1] == 0
    assert np.all(np.isnan(result.conditional[-1]))
    stats = result.metrics(np.array([[1, 0], [1, 0]]), P)
    assert stats["successful_target_probability"][-1] == 0
    assert np.isnan(stats["conditional_fidelity"][-1])


def test_symmetry_preserving_limit_has_unit_success_and_unchanged_dynamics():
    N = 4
    P = symmetric_projector(N)
    B = dicke_isometry_bulk(N, 0).toarray()
    initial = (B[:, 1]+1j*B[:, 3])/np.sqrt(2)
    H = lambda t: full_hamiltonian_reference(FullParameters(
        np.full(N, 1+.2*t), np.full(N, -.3+t), np.full(N, .4*t),
        .13*(np.ones((N, N))-np.eye(N))))
    times = np.linspace(0, 1, 13)
    selected = selective_evolution(H, initial, times, P, dt=.04)
    bare = evolve(H, initial, times, dt=.025)
    np.testing.assert_allclose(selected.success, 1, atol=2e-12, rtol=0)
    np.testing.assert_allclose(selected.branch, bare, atol=2e-12, rtol=0)


@pytest.mark.parametrize("N", range(2, 8))
def test_spin_spectrum_gap_and_phase_separation(N):
    j = N/2
    J2 = spin_squared(N).toarray()
    eigenvalues = np.linalg.eigvalsh(J2)
    expected = np.array([s*(s+1) for s in np.arange(N % 2 / 2, j+1, 1)])
    assert max(min(abs(expected-x)) for x in eigenvalues) < 1e-12
    penalty = np.linalg.eigvalsh(gap_penalty(N, .37).toarray())
    np.testing.assert_allclose(penalty[penalty > 1e-10].min(), N*.37, atol=1e-12)
    B = dicke_isometry_bulk(N, 0).toarray()
    np.testing.assert_allclose(J2@B, j*(j+1)*B, atol=2e-14)
    plan = spin_phase_register(N)
    assert len(set(plan["phases"])) == N//2+1
    assert max(plan["phases"]) < 1
    # Exact finite-register Fourier readout for each allowed eigenphase.
    R = plan["phase_denominator"]
    for label, phase in zip(plan["eigenvalue_labels"], plan["phases"]):
        readout = np.fft.fft(np.exp(2j*np.pi*np.arange(R)*phase))/R
        assert abs(readout[label])**2 > 1-1e-12


def test_gap_suppression_for_static_off_resonant_case():
    H, P = broken(), symmetric_projector(4)
    initial = np.eye(16)[0]
    losses = []
    for G in (0, 1, 8):
        psi = expm(-1j*(H+gap_penalty(4, G).toarray())*2)@initial
        losses.append(np.linalg.norm((np.eye(16)-P)@psi)**2)
    assert losses[-1] < losses[0]/100


def test_permutation_average_pairs_and_cycle_convergence():
    N = 4
    pairs = np.zeros((N, N))
    pairs[0, 1] = pairs[1, 0] = 2.4
    H = full_hamiltonian_reference(FullParameters([1]*N, [.3, -.2, .7, -.1], [.2]*N, pairs))
    group = symmetric_group(N)
    average = permutation_average(H, group)
    cyclic = permutation_average(H, [tuple(np.roll(np.arange(N), i)) for i in range(N)])
    swap = permutation_operator([0, 2, 1, 3])
    np.testing.assert_allclose((swap.T@average@swap-average).toarray(), 0, atol=1e-14)
    assert np.linalg.norm((swap.T@cyclic@swap-cyclic).toarray()) > .1
    initial = np.eye(16)[0]
    target = expm(-1j*.7*average.toarray())@initial
    errors = []
    for cycles in (1, 2, 4):
        result = permutation_cycles(lambda t: H, initial, .7, cycles, group,
                                    dt=.2, method="midpoint")
        errors.append(np.linalg.norm(result["states"][-1]-target))
        assert result["swap_count"] == cycles*2*sum(len(swap_decomposition(p)) for p in group)
    assert errors[2] < errors[1] < errors[0]


def test_swap_compilation_and_finite_pulses():
    for p in symmetric_group(4):
        U = np.eye(16)
        for a, b in swap_decomposition(p):
            order = list(range(4))
            order[a], order[b] = order[b], order[a]
            U = permutation_operator(order)@U
        np.testing.assert_array_equal(U, permutation_operator(p).toarray())
    H = broken(2)
    group = symmetric_group(2)
    ideal = permutation_cycles(lambda t: H, np.eye(4)[0], .4, 2, group)
    finite = permutation_cycles(lambda t: H, np.eye(4)[0], .4, 2, group, pulse_duration=1e-5)
    assert abs(finite["times"][-1]-.4) < 1e-14
    assert np.linalg.norm(ideal["states"][-1]-finite["states"][-1]) < 1e-4
    assert abs(np.linalg.norm(finite["states"][-1])-1) < 1e-10
    with pytest.raises(ValueError):
        permutation_cycles(lambda t: H, np.eye(4)[0], .4, 2, group, pulse_duration=1)


@pytest.mark.parametrize("N", range(2, 8))
def test_disorder_analytic_variance_and_short_time_quartic(N):
    rng = np.random.default_rng(71+N)
    epsilon = rng.normal(size=N)
    epsilon -= epsilon.mean()
    variance = np.mean(epsilon**2)
    B = dicke_isometry_bulk(N, 0).toarray()
    # H_delta=-sum epsilon*n; sign does not affect the norm identity.
    D = full_hamiltonian_reference(FullParameters(np.zeros(N), epsilon, np.zeros(N), np.zeros((N, N)))).toarray()
    np.testing.assert_allclose(B.conj().T@D@B, 0, atol=2e-15)
    Q = np.eye(2**N)-B@B.conj().T
    for k in range(N+1):
        np.testing.assert_allclose(np.linalg.norm(Q@D@B[:, k])**2,
                                  variance*k*(N-k)/(N-1), atol=2e-14)
    omega = 1.2
    H = full_hamiltonian_reference(FullParameters(np.full(N, omega), epsilon,
                                                np.zeros(N), np.zeros((N, N)))).toarray()
    ratios = []
    for t in (.02, .01, .005):
        psi = expm(-1j*t*H)@np.eye(2**N)[0]
        loss = np.linalg.norm(Q@psi)**2  # avoid cancellation in 1-||P psi||²
        ratios.append(loss/(N*omega**2*variance*t**4/16))
    assert abs(ratios[-1]-1) < 1e-4
    assert abs(ratios[-1]-1) < abs(ratios[0]-1)
