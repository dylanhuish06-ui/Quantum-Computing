"""Static propagation checks, including an independent SciPy expm reference."""

from itertools import product
from unittest.mock import patch

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project import (
    BASIS_STATES, KET_00, PSI_PLUS, evolve_static, normalise_state,
    rydberg_hamiltonian, state_probability,
)
from quantum_project.units import mhz_to_angular_frequency, microseconds_to_seconds


def evolve_expm_reference(psi0, times, Hbar):
    """Validation only: apply a direct matrix exponential at each supplied time.

    Inputs are valid normalised states, real times (s), and Hermitian Hbar
    (rad/s). This deliberately uses no eigenmethod or production solver code.
    """
    return np.array([expm(-1j * Hbar * t) @ psi0 for t in times],
                    dtype=complex).reshape(-1, 4)


def representative_hamiltonian(delta_mhz=0.4, interaction_mhz=5.0, phase=0.73):
    return rydberg_hamiltonian(
        omega=mhz_to_angular_frequency(1.0),
        delta=mhz_to_angular_frequency(delta_mhz),
        phase=phase,
        V=mhz_to_angular_frequency(interaction_mhz),
    )


COMPLEX_STATE = normalise_state([1 + 2j, -3j, 4, -5 + 1j])
TIMES = np.linspace(0, microseconds_to_seconds(4), 81)


@pytest.mark.parametrize("initial", [KET_00, PSI_PLUS, COMPLEX_STATE])
def test_initial_state_shape_and_complex_dtype(initial):
    states = evolve_static(initial, [0], representative_hamiltonian())
    assert states.shape == (1, 4)
    assert np.issubdtype(states.dtype, np.complexfloating)
    np.testing.assert_allclose(states[0], initial, atol=1e-14, rtol=0)


def test_agreement_with_independent_expm_reference():
    # All eight combinations cover zero/nonzero detuning, interaction and phase.
    # Compare complex amplitudes directly, including their global phase.
    maximum_error = 0.0
    maximum_norm_error = 0.0
    for delta, interaction, phase in product([0.0, -0.4], [0.0, 5.0], [0.0, 0.73]):
        matrix = representative_hamiltonian(delta, interaction, phase)
        for initial in [KET_00, PSI_PLUS, COMPLEX_STATE]:
            actual = evolve_static(initial, TIMES, matrix)
            reference = evolve_expm_reference(initial, TIMES, matrix)
            maximum_error = max(maximum_error, float(np.max(np.abs(actual - reference))))
            norms = np.linalg.norm(actual, axis=1)
            maximum_norm_error = max(maximum_norm_error, float(np.max(np.abs(norms - 1))))
            np.testing.assert_allclose(actual, reference, atol=1e-11, rtol=0)
            np.testing.assert_allclose(norms, 1, atol=1e-13, rtol=0)
    print(f"Maximum static/expm amplitude discrepancy: {maximum_error:.6e}")
    print(f"Maximum state-norm deviation: {maximum_norm_error:.6e}")


def test_reference_does_not_use_eigenmethod():
    matrix = representative_hamiltonian()
    with patch("numpy.linalg.eigh", side_effect=AssertionError("unexpected diagonalisation")):
        states = evolve_expm_reference(KET_00, [0, 1e-6], matrix)
    np.testing.assert_allclose(states[0], KET_00, atol=1e-14, rtol=0)


def test_only_one_diagonalisation_for_all_times():
    with patch("numpy.linalg.eigh", wraps=np.linalg.eigh) as diagonalise:
        evolve_static(COMPLEX_STATE, TIMES, representative_hamiltonian())
    assert diagonalise.call_count == 1


@pytest.mark.parametrize("phase", [0, 0.73, np.pi / 2])
def test_exchange_symmetry_from_ground_state(phase):
    states = evolve_static(KET_00, TIMES, representative_hamiltonian(phase=phase))
    np.testing.assert_allclose(states[:, 1], states[:, 2], atol=1e-13, rtol=0)


def test_basis_probabilities_sum_to_one():
    states = evolve_static(COMPLEX_STATE, TIMES, representative_hamiltonian())
    probabilities = np.array([
        [state_probability(state, basis) for basis in BASIS_STATES] for state in states
    ])
    np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=1e-13, rtol=0)


def test_diagonal_hamiltonian_has_correct_phase_sign_and_time_origin():
    eigenvalues = mhz_to_angular_frequency(1) * np.array([0.2, -0.3, 0.7, 1.1])
    times = np.array([0.7e-6, -0.1e-6, 0.7e-6, 0])
    actual = evolve_static(COMPLEX_STATE, times, np.diag(eigenvalues))
    expected = COMPLEX_STATE * np.exp(-1j * times[:, None] * eigenvalues)
    np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)


def test_zero_hamiltonian_preserves_state_despite_degeneracy():
    states = evolve_static(COMPLEX_STATE, TIMES, np.zeros((4, 4)))
    np.testing.assert_allclose(states, np.tile(COMPLEX_STATE, (len(TIMES), 1)),
                               atol=1e-14, rtol=0)


def test_empty_times_return_empty_state_array():
    states = evolve_static(KET_00, [], representative_hamiltonian())
    assert states.shape == (0, 4)
    assert np.issubdtype(states.dtype, np.complexfloating)


def test_inputs_are_not_modified():
    initial, times, matrix = COMPLEX_STATE.copy(), TIMES.copy(), representative_hamiltonian()
    originals = [value.copy() for value in (initial, times, matrix)]
    evolve_static(initial, times, matrix)
    for value, original in zip((initial, times, matrix), originals):
        np.testing.assert_array_equal(value, original)


@pytest.mark.parametrize("invalid", [
    [1, 0], np.ones((4, 1)), np.ones((1, 4)), [2, 0, 0, 0],
    [0, 0, 0, 0], [np.nan, 0, 0, 0], [np.inf, 0, 0, 0],
])
def test_invalid_initial_state(invalid):
    with pytest.raises(ValueError):
        evolve_static(invalid, TIMES, representative_hamiltonian())


@pytest.mark.parametrize("invalid", [
    0.0, [[0, 1e-6]], [np.nan], [np.inf], [-np.inf], [1j], [0j], [True], ["0"],
])
def test_invalid_times(invalid):
    with pytest.raises(ValueError, match="times"):
        evolve_static(KET_00, invalid, representative_hamiltonian())


@pytest.mark.parametrize("invalid", [
    np.eye(3), np.zeros(4), np.zeros((1, 4, 4)), np.full((4, 4), np.nan),
    np.full((4, 4), np.inf), 1j * np.eye(4), np.triu(np.ones((4, 4))),
])
def test_invalid_hamiltonian(invalid):
    with pytest.raises(ValueError):
        evolve_static(KET_00, TIMES, invalid)
