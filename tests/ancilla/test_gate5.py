import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from quantum_project.ancilla.action import controlled_action
from quantum_project.ancilla.approximation import adaptive_static_grouping
from quantum_project.ancilla.models import FullParameters, full_hamiltonian_reference
from quantum_project.ancilla.propagation import evolve


@pytest.mark.parametrize("duration", [0., .1, -2., 3.])
def test_controlled_action_against_dense(duration):
    rng = np.random.default_rng(73)
    A = rng.normal(size=(12, 12))+1j*rng.normal(size=(12, 12))
    H = A+A.conj().T+17*np.eye(12)
    state = rng.normal(size=12)+1j*rng.normal(size=12)
    for tol in (1e-5, 1e-9, 1e-13):
        actual, info = controlled_action(sparse.csr_matrix(H), state, duration, atol=tol)
        expected = expm(-1j*duration*H)@state
        assert info["truncation_bound"] <= tol
        assert np.linalg.norm(actual-expected) <= info["truncation_bound"]+2e-12
        assert abs(np.linalg.norm(actual)-np.linalg.norm(state)) < tol+2e-12


def test_action_resource_and_input_guards():
    for kwargs in ({"atol": 0}, {"duration": np.nan}):
        with pytest.raises(ValueError):
            controlled_action(np.eye(2), [1, 0], **({"duration": 1}|kwargs))
    with pytest.raises(RuntimeError):
        controlled_action(np.diag([-100, 100]), [1, 0], 1, max_matvecs=2)
    with pytest.raises(ValueError):
        controlled_action([[0, 1], [0, 0]], [1, 0], 1)


def test_generic_disorder_failure_and_adaptive_bound():
    N, T = 5, .7
    H = full_hamiltonian_reference(FullParameters(np.ones(N), [-.8, .7, -.4, .3, 1.1],
                                                 np.zeros(N), np.zeros((N, N))))
    fixed = adaptive_static_grouping(H, N, T, error_budget=1e-4, max_m=0)
    assert not fixed["budget_met"]
    assert fixed["history"][-1]["state_error_bound"] > .1
    for budget in (1., .1, 1e-4):
        result = adaptive_static_grouping(H, N, T, error_budget=budget)
        B = result["isometry"].toarray()
        initial = np.eye(B.shape[1])[1]
        exact = expm(-1j*T*H.toarray())@(B@initial)
        grouped = B@expm(-1j*T*result["hamiltonian"].toarray())@initial
        assert np.linalg.norm(exact-grouped) <= result["history"][-1]["state_error_bound"]+1e-12
        assert result["budget_met"]
    assert result["history"][-1]["dimension"] == 2**N


def test_approximation_guard():
    with pytest.raises(ValueError):
        adaptive_static_grouping(sparse.eye(512), 9, 1)


def test_controlled_midpoint_separates_action_and_time_error():
    X = np.array([[0, 1], [1, 0]], complex)
    Z = np.diag([1., -1.])
    H = lambda t: (1+.3*t)*X+np.sin(t)*Z
    initial = np.array([1., 0.])
    times = np.linspace(0, 1, 11)
    reference = evolve(H, initial, times, dt=.01)
    errors = []
    for dt in (.05, .025, .0125):
        states = evolve(H, initial, times, method="controlled_midpoint", dt=dt, action_atol=1e-12)
        scipy_states = evolve(H, initial, times, method="midpoint", dt=dt)
        assert np.max(abs(states-scipy_states)) < 2e-12
        errors.append(np.linalg.norm(states[-1]-reference[-1]))
    assert 3.9 < errors[1]/errors[2] < 4.1


@pytest.mark.parametrize("bad", [0, np.nan])
def test_action_budget_validation_in_evolve(bad):
    with pytest.raises(ValueError):
        evolve(lambda t: np.eye(2), [1, 0], [0, 1],
               method="controlled_midpoint", action_atol=bad)
