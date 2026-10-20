"""Independent small-space audit; no benchmark or recorded-data dependency."""

import itertools
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm

import quantum_project
from quantum_project import KET_00, evolve_static, rydberg_hamiltonian
from quantum_project.symmetry import (
    compare_states, dicke_hamiltonian, dicke_isometry, evolve_uniform,
    exchange_hamiltonian, exchange_transformation, full_hamiltonian,
    local_operator, reduced_two_atom_hamiltonian,
)

SCALE = 2 * np.pi * 1e6  # All energies below are H/hbar in rad/s; times in s.
PARAMETERS = list(itertools.product([1., -2.], [0., .3], [0., .37], [0., 2., -3.]))


def update_maxima(metrics, **values):
    for name, value in values.items():
        metrics[name] = max(metrics.get(name, 0.), float(value))


def state_errors(full, expanded):
    """Align expanded by exp(-i arg(<full|expanded>)) at each output time."""
    overlap = np.einsum("ti,ti->t", full.conj(), expanded)
    aligned = expanded * np.exp(-1j * np.angle(overlap))[:, None]
    return dict(
        fidelity_error=np.max(abs(abs(overlap)**2 - 1)),  # No clipping.
        probability_error=np.max(abs(abs(full)**2 - abs(expanded)**2)),
        aligned_amplitude_error=np.max(abs(full - aligned)),
        expanded_norm_error=np.max(abs(np.linalg.norm(expanded, axis=1) - 1)),
    )


def audit_n2():
    s = 1 / np.sqrt(2)
    # NEW vectors as COLUMNS in |00>, |0r>, |r0>, |rr> coordinates.
    U = np.array([[1, 0, 0, 0], [0, s, s, 0], [0, s, -s, 0], [0, 0, 0, 1]], complex)
    B = U[:, [0, 1, 3]]
    np.testing.assert_array_equal(exchange_transformation(), U)
    np.testing.assert_allclose(U.conj().T @ U, np.eye(4), atol=3e-16, rtol=0)
    metrics = {"unitarity_error": float(np.max(abs(U.conj().T @ U - np.eye(4))))}
    times = np.linspace(0, .5e-6, 51)
    for omega, delta, phase, interaction in PARAMETERS:
        args = (omega*SCALE, delta*SCALE, phase, interaction*SCALE)
        o, d, p, v = args
        c = o * np.exp(-1j*p) / np.sqrt(2)
        expected4 = np.array([[d, c, 0, 0], [c.conjugate(), 0, 0, c],
                              [0, 0, 0, 0], [0, c.conjugate(), 0, v-d]])
        # Construct independently, before comparing any extracted block.
        expected3 = np.array([[d, c, 0], [c.conjugate(), 0, c], [0, c.conjugate(), v-d]])
        H = rydberg_hamiltonian(*args)
        elements = np.array([[np.vdot(U[:, i], H @ U[:, j]) for j in range(4)] for i in range(4)])
        transformed = exchange_hamiltonian(*args)
        for actual in (elements, U.conj().T @ H @ U, transformed):
            np.testing.assert_allclose(actual, expected4, atol=1e-9, rtol=0)
        np.testing.assert_allclose(H @ U, U @ expected4, atol=2e-9, rtol=0)
        np.testing.assert_allclose(transformed[2], 0, atol=1e-9, rtol=0)
        np.testing.assert_allclose(transformed[:, 2], 0, atol=1e-9, rtol=0)
        np.testing.assert_array_equal(full_hamiltonian(2, *args).toarray(), H)
        for actual in (transformed[np.ix_([0, 1, 3], [0, 1, 3])],
                       reduced_two_atom_hamiltonian(*args), dicke_hamiltonian(2, *args).toarray()):
            np.testing.assert_allclose(actual, expected3, atol=4e-9, rtol=0)
        # Rephase the columns so U != U^dagger: the original real symmetric U
        # alone cannot distinguish the correct transformation from its inverse.
        Q = U @ np.diag(np.exp(1j*np.array([.2, -.4, .7, 1.1])))
        direct = np.array([[np.vdot(Q[:, i], H @ Q[:, j]) for j in range(4)] for i in range(4)])
        np.testing.assert_allclose(Q.conj().T @ H @ Q, direct, atol=4e-9, rtol=0)
        assert np.max(abs(Q @ H @ Q.conj().T - direct)) > .1*SCALE
        full = evolve_static(KET_00, times, H)
        # Dense exponential of an independent 3x3 matrix versus milestone eigensolver.
        reduced = np.array([expm(-1j*expected3*t)[:, 0] for t in times])
        expanded = np.column_stack((reduced[:, 0], reduced[:, 1]*s,
                                    reduced[:, 1]*s, reduced[:, 2]))
        np.testing.assert_array_equal(expanded, reduced @ B.T)
        errors = state_errors(full, expanded)
        assert max(errors.values()) < 2e-13
        dark = float(np.max(abs(full @ U[:, 2].conj())**2))
        assert dark < 1e-28
        update_maxima(metrics, **errors, dark_probability=dark,
                      matrix_entry_error_rad_s=np.max(abs(elements-expected4)),
                      tensor_matrix_error_rad_s=np.max(abs(full_hamiltonian(2, *args).toarray()-H)))
    return {"parameter_sets": len(PARAMETERS), "samples": len(times), **metrics}


def audit_small_n():
    metrics = {}
    for N in (2, 3, 4):
        sectors = np.array([i.bit_count() for i in range(2**N)])
        # Independent equal-weight basis: no use of the implementation isometry.
        B = np.column_stack([(sectors == k)/np.sqrt(comb(N, k)) for k in range(N+1)])
        np.testing.assert_array_equal(dicke_isometry(N).toarray(), B)
        np.testing.assert_allclose(B.T @ B, np.eye(N+1), atol=4e-15, rtol=0)
        for omega, delta, phase, interaction in PARAMETERS:
            args = (omega*SCALE, delta*SCALE, phase, interaction*SCALE)
            o, d, p, v = args
            # Explicit bit flips, independent of Kronecker embedding.
            reference = np.diag(d*(N/2-sectors) + v*sectors*(sectors-1)/2).astype(complex)
            for col in range(2**N):
                for bit in range(N):
                    row = col ^ (1 << bit)
                    reference[row, col] = o/2*np.exp((1j if row > col else -1j)*p)
            full = full_hamiltonian(N, *args).toarray()
            reduced = dicke_hamiltonian(N, *args).toarray()
            assert full.shape == (2**N, 2**N) and reduced.shape == (N+1, N+1)
            np.testing.assert_allclose(full/SCALE, reference/SCALE, atol=4e-15, rtol=0)
            np.testing.assert_allclose((full @ B)/SCALE, (B @ reduced)/SCALE, atol=4e-15, rtol=0)
            full_states = np.array([expm(-1j*reference*t)[:, 0] for t in np.linspace(0, .25e-6, 31)])
            reduced_states = evolve_uniform(reduced, .25e-6, 31)
            expanded = reduced_states @ B.T
            errors = state_errors(full_states, expanded)
            assert max(errors.values()) < 2e-12
            checks = compare_states(full_states, reduced_states, N)
            np.testing.assert_allclose(checks["fidelity"],
                                       abs(np.einsum("ti,ti->t", full_states.conj(), expanded))**2,
                                       atol=2e-15, rtol=0)
            assert checks["infidelity"] < 2e-12
            assert checks["max_excitation_error"] < 2e-12
            assert checks["max_dicke_population_error"] < 2e-12
            update_maxima(metrics, **errors, excitation_error=checks["max_excitation_error"],
                          dicke_population_error=checks["max_dicke_population_error"],
                          intertwining_error_over_scale=np.max(abs(full @ B - B @ reduced))/SCALE)
    return {"N": [2, 3, 4], "parameter_sets_per_N": len(PARAMETERS), "samples": 31, **metrics}


def test_imports_use_this_worktree():
    root = Path(__file__).resolve().parents[1]
    assert Path(quantum_project.__file__).resolve() == root / "src/quantum_project/__init__.py"


def test_independent_n2_audit():
    print("N2 recovery audit:", json.dumps(audit_n2(), sort_keys=True))


def test_independent_small_n_audit():
    print("Small-N recovery audit:", json.dumps(audit_small_n(), sort_keys=True))


def test_tensor_site_order_is_most_significant_first():
    # Uniform drives cannot reveal a site reversal, so check a local operator.
    for N in (2, 3, 4):
        for site in range(N):
            actual = local_operator([[0, 0], [0, 1]], site, N).toarray()
            expected = np.diag([(i >> (N-1-site)) & 1 for i in range(2**N)])
            np.testing.assert_array_equal(actual, expected)


def test_global_phase_and_explicit_n2_expansion():
    reduced = np.array([[1, 1j, -2]], complex) / np.sqrt(6)
    expanded = np.column_stack((reduced[:, 0], reduced[:, 1]/np.sqrt(2),
                                reduced[:, 1]/np.sqrt(2), reduced[:, 2]))
    full = expanded * np.exp(1.234j)
    checks = compare_states(full, reduced, 2)
    np.testing.assert_allclose(checks["fidelity"], 1, atol=1e-15, rtol=0)
    assert checks["max_amplitude_error"] > .5  # Legacy metric is explicitly raw.
    assert checks["max_aligned_amplitude_error"] < 1e-15
    assert checks["max_absolute_fidelity_error"] < 1e-15
    assert max(state_errors(full, expanded).values()) < 1e-15
    for name in ("max_excitation_error", "max_ground_probability_error", "max_dicke_population_error"):
        assert checks[name] < 1e-15
