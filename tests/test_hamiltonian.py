import numpy as np
import pytest

from quantum_project import rydberg_hamiltonian
from quantum_project.validation import is_hermitian


@pytest.mark.parametrize("parameters", [
    (0, 0, 0, 0), (2, 3, 0, 7), (-4, -2, 0.73, -5),
    (2 * np.pi * 1e6, -2 * np.pi * 2e6, 1.31, 2 * np.pi * 10e6),
])
def test_shape_dtype_and_hermiticity(parameters):
    matrix = rydberg_hamiltonian(*parameters)
    assert matrix.shape == (4, 4)
    assert np.issubdtype(matrix.dtype, np.complexfloating)
    assert is_hermitian(matrix)
    np.testing.assert_allclose(matrix, matrix.conj().T, atol=1e-12, rtol=0)


def test_zero_phase_matches_equation_in_full():
    expected = np.array([[3, 1, 1, 0], [1, 0, 0, 1],
                         [1, 0, 0, 1], [0, 1, 1, 4]], dtype=complex)
    np.testing.assert_array_equal(rydberg_hamiltonian(2, 3, 0, 7), expected)


def test_nonzero_phase_has_correct_signs_and_factor_of_two():
    # phi = pi/2 gives exp(-i phi) = -i on each upper coupling.
    expected = np.array([[3, -1j, -1j, 0], [1j, 0, 0, -1j],
                         [1j, 0, 0, -1j], [0, 1j, 1j, 4]])
    np.testing.assert_allclose(rydberg_hamiltonian(2, 3, np.pi / 2, 7),
                               expected, atol=1e-14, rtol=0)


def test_no_drive_leaves_detuning_and_interaction_diagonal():
    np.testing.assert_array_equal(rydberg_hamiltonian(0, -3, 0.8, 11),
                                  np.diag([-3, 0, 0, 14]))


def test_antisymmetric_state_is_dark():
    antisymmetric = np.array([0, 1, -1, 0]) / np.sqrt(2)
    np.testing.assert_allclose(rydberg_hamiltonian(2, 3, 0.7, 7) @ antisymmetric,
                               np.zeros(4), atol=1e-14, rtol=0)


@pytest.mark.parametrize("name", ["omega", "delta", "phase", "V"])
@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf, 1j, [1], "1", True])
def test_invalid_parameter(name, invalid):
    parameters = dict(omega=1, delta=0, phase=0, V=2)
    parameters[name] = invalid
    with pytest.raises(ValueError, match="finite real scalar"):
        rydberg_hamiltonian(**parameters)


def test_unrepresentable_diagonal_is_rejected():
    with pytest.raises(ValueError, match="finite numeric range"):
        rydberg_hamiltonian(0, -1e308, 0, 1e308)
