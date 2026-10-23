import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project.ancilla import (
    FullParameters, HybridParameters, dicke_isometry_bulk,
    full_hamiltonian_reference, hybrid_hamiltonian, validate_bulk_symmetry,
    evolve, evolve_hybrid,
)
from quantum_project.hamiltonian import rydberg_hamiltonian
from quantum_project.symmetry import reduced_two_atom_hamiltonian
from quantum_project.trotter import evolve_trotter
from quantum_project.ancilla.models import split_bulk_qubit


def profile(N, m, seed):
    rng = np.random.default_rng(seed)
    omega, delta, phase, cross = rng.normal(size=(4, m))
    pairs = rng.normal(size=(m, m))
    pairs = np.triu(pairs, 1)
    pairs += pairs.T
    def parameters(t):
        return HybridParameters(N, m, 1+0.3*np.sin(2*t), -.4+t*t, .2+.6*t,
                                .12*np.cos(t), omega*(1+.2*t), delta+t,
                                phase-.7*t, cross*(1+.1*t), pairs*np.cos(t))
    return parameters


@pytest.mark.parametrize("N,m", [(N, m) for N in range(2, 9) for m in range(min(3, N)+1)])
def test_exact_hybrid_matrix_and_complex_trajectories(N, m):
    parameters = profile(N, m, 841+10*N+m)
    B = dicke_isometry_bulk(N, m).toarray()
    np.testing.assert_allclose(B.conj().T@B, np.eye(B.shape[1]), atol=3e-15)
    for t in (0., .137, .42):
        p = parameters(t)
        H = full_hamiltonian_reference(p.full_parameters()).toarray()
        h = hybrid_hamiltonian(p).toarray()
        assert h.shape == ((N-m+1)*2**m,)*2
        np.testing.assert_allclose(h, h.conj().T, atol=1e-15)
        np.testing.assert_allclose(B.conj().T@H@B, h, atol=5e-14)
        np.testing.assert_allclose(H@B-B@h, 0, atol=5e-14)
        np.testing.assert_allclose((np.eye(2**N)-B@B.conj().T)@H@B, 0, atol=5e-14)
        recovered = validate_bulk_symmetry(p.full_parameters(), m)
        np.testing.assert_allclose(hybrid_hamiltonian(recovered).toarray(), h, atol=1e-14)
    rng = np.random.default_rng(111+N+m)
    random = rng.normal(size=B.shape[1])+1j*rng.normal(size=B.shape[1])
    random /= np.linalg.norm(random)
    for initial in (np.eye(1, B.shape[1], 0).ravel(), np.eye(1, B.shape[1], B.shape[1]-1).ravel(), random):
        times = np.linspace(0, .42, 5)
        reduced = evolve_hybrid(parameters, initial, times, dt=.1)
        full = evolve(lambda t: full_hamiltonian_reference(parameters(t).full_parameters()),
                      B@initial, times, dt=.08)
        assert np.max(np.linalg.norm(full-(B@reduced.T).T, axis=1)) < 1e-9
        assert np.max(abs(np.sum(abs(full)**2, axis=1)-1)) < 1e-10


def test_milestone_offsets_phases_and_original_propagator():
    omega, delta, phase, V = 1.3, -.72, 1.19, 2.1
    p = HybridParameters(2, 0, omega, delta, phase, V)
    H = full_hamiltonian_reference(p.full_parameters()).toarray()
    expected = rydberg_hamiltonian(omega, delta, phase, V)
    np.testing.assert_allclose(H, expected, atol=1e-15)
    np.testing.assert_allclose(np.diag(H), [delta, 0, 0, V-delta])
    np.testing.assert_allclose(H[1, 0], omega/2*np.exp(1j*phase))
    np.testing.assert_allclose(hybrid_hamiltonian(p).toarray(),
                               reduced_two_atom_hamiltonian(omega, delta, phase, V), atol=1e-15)
    initial = np.array([1, 0, 0, 0], dtype=complex)
    times = np.linspace(0, .9, 13)
    np.testing.assert_allclose(evolve(lambda t: H, initial, times),
                              evolve_trotter(initial, times, lambda t: expected), atol=1e-11)
    # Dropping the identity offset gives the wrong raw complex amplitudes.
    shifted = expm(-1j*(H-delta*np.eye(4))*.9)@initial
    assert np.linalg.norm(shifted-expm(-1j*H*.9)@initial) > .1


@pytest.mark.parametrize("kind", ["delta", "omega", "phase", "bulk_pair", "cross"])
def test_reject_bulk_violations(kind):
    p = profile(6, 2, 2)(.3).full_parameters()
    arrays = {key: getattr(p, key).copy() for key in ("omega", "delta", "phase", "pairs")}
    if kind in arrays:
        arrays[kind][0] += .1
    else:
        j = 1 if kind == "bulk_pair" else 4
        arrays["pairs"][0, j] += .1
        arrays["pairs"][j, 0] += .1
    with pytest.raises(ValueError, match="symmetry"):
        validate_bulk_symmetry(FullParameters(**arrays), 2)


@pytest.mark.parametrize("N,m,dimension", [(64, 1, 128), (256, 2, 1020), (1024, 1, 2048), (1024, 2, 4092)])
def test_large_direct_sparse(N, m, dimension, monkeypatch):
    def no_expansion(*args, **kwargs):
        raise AssertionError("large-N expansion attempted")
    monkeypatch.setattr(HybridParameters, "full_parameters", no_expansion)
    p = profile(N, m, 12)(.1)
    H = hybrid_hamiltonian(p)
    assert H.shape == (dimension, dimension)
    assert H.nnz <= (m+3)*dimension
    assert H.data.nbytes < 2**20


def test_invalid_inputs_and_all_defects():
    for args in ((2, 3), (True, 0), (-1, 0)):
        with pytest.raises(ValueError):
            HybridParameters(*args)
    with pytest.raises(ValueError):
        HybridParameters(4, 1)
    with pytest.raises(ValueError):
        HybridParameters(2, 0, phase_bulk=np.nan)
    with pytest.raises(ValueError):
        FullParameters([1, 1], [0, 0], [0, 0], [[0, 1], [2, 0]])
    with pytest.raises(ValueError):
        dicke_isometry_bulk(1024, 1)
    with pytest.raises(ValueError):
        hybrid_hamiltonian(HybridParameters(1000000, 0))
    p = profile(3, 3, 91)(.13)
    np.testing.assert_array_equal(hybrid_hamiltonian(p).toarray(),
                                  full_hamiltonian_reference(p.full_parameters()).toarray())


@pytest.mark.parametrize("N,m", [(3, 0), (6, 2), (6, 5)])
def test_exact_state_transfer_when_splitting_bulk(N, m):
    rng = np.random.default_rng(816+N+m)
    old = rng.normal(size=(N-m+1)*2**m)+1j*rng.normal(size=(N-m+1)*2**m)
    old /= np.linalg.norm(old)
    new = split_bulk_qubit(N, m, old)
    assert abs(np.linalg.norm(new)-1) < 1e-15
    np.testing.assert_allclose(dicke_isometry_bulk(N, m+1)@new,
                              dicke_isometry_bulk(N, m)@old, atol=1e-15)
