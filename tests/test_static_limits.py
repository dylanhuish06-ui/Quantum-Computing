import numpy as np
import pytest

from quantum_project import KET_00, PSI_PLUS, evolve_static, rydberg_hamiltonian
from quantum_project.static_limits import blockade_state, noninteracting_state
from quantum_project.units import mhz_to_angular_frequency


OMEGA = mhz_to_angular_frequency(1)
TIMES = np.linspace(0, 1e-6, 4001)


def numerical_state(interaction_mhz):
    hbar = rydberg_hamiltonian(OMEGA, 0, 0, mhz_to_angular_frequency(interaction_mhz))
    return evolve_static(KET_00, TIMES, hbar)


def observables(states):
    return np.column_stack((np.abs(states)**2, np.abs(states @ PSI_PLUS.conj())**2))


def test_noninteracting_solution_matches_complex_amplitudes():
    np.testing.assert_allclose(numerical_state(0), noninteracting_state(TIMES, OMEGA),
                               atol=1e-13, rtol=0)


def test_noninteracting_half_and_full_rabi_cycle():
    states = noninteracting_state([0, 0.5e-6, 1e-6], OMEGA)
    np.testing.assert_allclose(states, [[1, 0, 0, 0], [0, 0, 0, -1], [1, 0, 0, 0]],
                               atol=1e-14, rtol=0)
    quarter_cycle = noninteracting_state([0.25e-6], OMEGA)
    assert abs(quarter_cycle[0] @ PSI_PLUS.conj())**2 == pytest.approx(0.5)


def test_blockade_collective_bell_pulse():
    bell_time = np.pi / (np.sqrt(2) * OMEGA)
    np.testing.assert_allclose(blockade_state([bell_time], OMEGA)[0], -1j * PSI_PLUS,
                               atol=1e-14, rtol=0)


@pytest.mark.parametrize("limit", [noninteracting_state, blockade_state])
def test_limit_states_preserve_normalisation_and_symmetry(limit):
    states = limit(TIMES, OMEGA)
    np.testing.assert_allclose(np.sum(np.abs(states)**2, axis=1), 1, atol=1e-14, rtol=0)
    np.testing.assert_array_equal(states[:, 1], states[:, 2])
    np.testing.assert_array_equal(states[0], KET_00)


def test_small_interaction_converges_quadratically_to_zero_limit():
    expected = observables(noninteracting_state(TIMES, OMEGA))
    errors = [np.max(np.abs(observables(numerical_state(v)) - expected)) for v in (0.01, 0.001)]
    assert errors[1] < errors[0] / 80
    assert errors[1] < 1e-5


def test_finite_blockade_converges_with_suppressed_double_excitation():
    expected = observables(blockade_state(TIMES, OMEGA))
    previous_error = np.inf
    for interaction in (10, 30, 100):
        actual = observables(numerical_state(interaction))
        error = np.max(np.abs(actual - expected))
        assert error < previous_error
        # Generous O((Omega/V)^2) envelopes over this fixed one-microsecond window.
        assert error < 10 / interaction**2
        assert np.max(actual[:, 3]) < 2 / interaction**2
        previous_error = error


def test_population_dynamics_are_even_in_interaction():
    np.testing.assert_allclose(observables(numerical_state(-100)),
                               observables(numerical_state(100)), atol=1e-12, rtol=0)


@pytest.mark.parametrize("interaction", [0, 0.5, 1, 5, 100])
def test_figure_cases_obey_probability_constraints(interaction):
    values = observables(numerical_state(interaction))
    np.testing.assert_allclose(values[:, :4].sum(axis=1), 1, atol=1e-12, rtol=0)
    np.testing.assert_allclose(values[:, 1], values[:, 2], atol=1e-12, rtol=0)
    np.testing.assert_allclose(values[:, 4], values[:, 1] + values[:, 2], atol=1e-12, rtol=0)
    assert np.all(values >= -1e-12)
    assert np.all(values <= 1 + 1e-12)


@pytest.mark.parametrize("limit", [noninteracting_state, blockade_state])
@pytest.mark.parametrize("times", [0, [[0]], [np.nan], [np.inf], [1j]])
def test_limit_time_validation(limit, times):
    with pytest.raises(ValueError, match="times"):
        limit(times, OMEGA)


@pytest.mark.parametrize("limit", [noninteracting_state, blockade_state])
def test_limits_accept_empty_times_and_reject_invalid_frequency(limit):
    assert limit([], OMEGA).shape == (0, 4)
    with pytest.raises(ValueError, match="omega"):
        limit([0], np.nan)
