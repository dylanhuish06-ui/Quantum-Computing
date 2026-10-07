import numpy as np
import pytest

from quantum_project import (
    BASIS_LABELS, BASIS_STATES, KET_00, KET_0R, KET_R0, KET_RR, PSI_PLUS,
    fidelity, normalise_state, state_probability,
)
from quantum_project.validation import are_orthonormal, is_normalised


def test_canonical_basis_order_and_orthogonality():
    assert BASIS_LABELS == ("|00>", "|0r>", "|r0>", "|rr>")
    vectors = np.array([KET_00, KET_0R, KET_R0, KET_RR])
    np.testing.assert_array_equal(vectors, np.eye(4))
    np.testing.assert_array_equal(BASIS_STATES, vectors)
    np.testing.assert_array_equal(vectors.conj() @ vectors.T, np.eye(4))
    assert are_orthonormal(vectors)
    assert all(is_normalised(vector) for vector in vectors)


def test_bell_state():
    np.testing.assert_allclose(PSI_PLUS, [0, 1 / np.sqrt(2), 1 / np.sqrt(2), 0])
    assert is_normalised(PSI_PLUS)


@pytest.mark.parametrize("state", [BASIS_STATES, KET_00, KET_0R, KET_R0, KET_RR, PSI_PLUS])
def test_exported_states_are_read_only(state):
    with pytest.raises(ValueError):
        state.flat[0] = 9


@pytest.mark.parametrize("state", [KET_00, KET_0R, KET_R0, KET_RR, PSI_PLUS,
                                         np.array([1, 1j, -1, -1j]) / 2])
def test_fidelity_with_self(state):
    assert fidelity(state, state) == pytest.approx(1)


def test_fidelity_respects_complex_conjugation_and_global_phase():
    state = np.array([1, 1j, 0, 0]) / np.sqrt(2)
    orthogonal = np.array([1, -1j, 0, 0]) / np.sqrt(2)
    assert fidelity(state, orthogonal) == pytest.approx(0, abs=1e-14)
    assert fidelity(state, np.exp(0.37j) * state) == pytest.approx(1)


def test_basis_and_bell_probabilities():
    np.testing.assert_allclose([state_probability(KET_R0, b) for b in BASIS_STATES],
                               [0, 0, 1, 0], atol=1e-14)
    np.testing.assert_allclose([state_probability(PSI_PLUS, b) for b in BASIS_STATES],
                               [0, 0.5, 0.5, 0], atol=1e-14)
    assert state_probability(KET_0R, PSI_PLUS) == pytest.approx(0.5)


def test_probabilities_sum_to_one():
    state = normalise_state([1 + 2j, -3j, 4, -5 + 1j])
    probabilities = [state_probability(state, basis) for basis in BASIS_STATES]
    assert sum(probabilities) == pytest.approx(1)
    assert all(0 <= probability <= 1 for probability in probabilities)


@pytest.mark.parametrize("scale", [1, 1e-300, 1e300])
def test_normalisation_preserves_relative_amplitudes_and_input(scale):
    original = np.array([1, 1j, -1, 0], dtype=complex) * scale
    snapshot = original.copy()
    state = normalise_state(original)
    np.testing.assert_allclose(state, np.array([1, 1j, -1, 0]) / np.sqrt(3))
    np.testing.assert_array_equal(original, snapshot)
    assert is_normalised(state)
    assert not np.shares_memory(original, state)


@pytest.mark.parametrize("invalid", [
    [0, 0, 0, 0], [1, 0], [[1], [0], [0], [0]],
    [np.nan, 0, 0, 0], [np.inf, 0, 0, 0],
])
def test_normalise_rejects_invalid_state(invalid):
    with pytest.raises(ValueError):
        normalise_state(invalid)


@pytest.mark.parametrize("observable", [fidelity, state_probability])
@pytest.mark.parametrize("invalid", [[0, 0, 0, 0], [2, 0, 0, 0],
                                       [1, 0], [np.nan, 0, 0, 0]])
def test_observables_validate_both_inputs(observable, invalid):
    with pytest.raises(ValueError):
        observable(invalid, KET_00)
    with pytest.raises(ValueError):
        observable(KET_00, invalid)
