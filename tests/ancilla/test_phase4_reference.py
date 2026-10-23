"""Independent tensor/algebra checks for the imported Phase 4 reference."""

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project.ancilla import FullParameters, full_hamiltonian_reference, evolve
from quantum_project.ancilla.models import dicke_isometry_bulk
from quantum_project.ancilla.syndrome import symmetric_projector


def tensor_local(op, site, N):
    result = np.ones((1, 1))
    for i in range(N):
        result = np.kron(result, op if site == i else np.eye(2))
    return result


@pytest.mark.parametrize("N", range(2, 7))
def test_arbitrary_full_generator_against_independent_tensor_construction(N):
    rng = np.random.default_rng(400+N)
    pairs = np.triu(rng.normal(size=(N, N)), 1)
    pairs += pairs.T
    p = FullParameters(*rng.normal(size=(3, N)), pairs)
    identity = np.eye(2**N)
    n = [tensor_local(np.diag([0., 1.]), i, N) for i in range(N)]
    H = np.zeros((2**N, 2**N), complex)
    for i in range(N):
        raising = tensor_local([[0, 0], [1, 0]], i, N)
        H += p.delta[i]*(identity/2-n[i])
        H += p.omega[i]/2*(np.exp(1j*p.phase[i])*raising + np.exp(-1j*p.phase[i])*raising.T)
        for j in range(i+1, N):
            H += p.pairs[i, j]*n[i]@n[j]
    actual = full_hamiltonian_reference(p).toarray()
    np.testing.assert_allclose(actual, H, atol=3e-15, rtol=0)
    np.testing.assert_allclose(actual, actual.conj().T, atol=0, rtol=0)
    psi = rng.normal(size=2**N)+1j*rng.normal(size=2**N)
    psi /= np.linalg.norm(psi)
    expected = expm(-1j*H*.37)@psi
    for method in ("dop853", "midpoint", "dense_midpoint"):
        result = evolve(lambda t: actual, psi, [0, .37], dt=.1, method=method)
        np.testing.assert_allclose(result[-1], expected, atol=4e-11, rtol=0)


@pytest.mark.parametrize("N", range(2, 9))
def test_projector_commutator_and_invariance(N):
    B = dicke_isometry_bulk(N, 0).toarray()
    P = symmetric_projector(N)
    Q = np.eye(2**N)-P
    np.testing.assert_allclose(B.conj().T@B, np.eye(N+1), atol=3e-15)
    np.testing.assert_allclose(P@P, P, atol=3e-15)
    p = FullParameters(np.ones(N), np.full(N, -.4), np.full(N, .7),
                       .2*(np.ones((N, N))-np.eye(N)))
    H = full_hamiltonian_reference(p).toarray()
    np.testing.assert_allclose(H@P-P@H, 0, atol=3e-15)
    np.testing.assert_allclose(Q@H@B, 0, atol=8e-15)


@pytest.mark.parametrize("method", ["dop853", "midpoint", "dense_midpoint"])
def test_solver_work_budget_and_zero_time_validation(method):
    H = lambda t: np.array([[t, 1], [1, -t]])
    with pytest.raises(RuntimeError, match="resource budget"):
        evolve(H, [1, 0], [0, 1], method=method, max_evaluations=2)
    with pytest.raises(ValueError, match="Hermitian"):
        evolve(lambda t: [[0, 1], [0, 0]], [1, 0], [0], method=method)


def test_reference_guards_before_exponential_allocation():
    with pytest.raises(ValueError, match="N<=12"):
        full_hamiltonian_reference(FullParameters(np.ones(13), np.zeros(13), np.zeros(13), np.zeros((13, 13))))
    with pytest.raises(ValueError, match="N<=8"):
        symmetric_projector(1024)
