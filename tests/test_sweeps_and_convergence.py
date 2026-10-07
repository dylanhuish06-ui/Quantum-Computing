import numpy as np
import pytest

from quantum_project import KET_00, evolve_trotter, rydberg_hamiltonian
from quantum_project.convergence import (
    compare_probabilities, largest_passing_timestep, significant_tolerance, sweep_probabilities,
)
from quantum_project.sweeps import cubic_detuning, linear_detuning
from quantum_project.units import mhz_to_angular_frequency


@pytest.mark.parametrize("sweep", [cubic_detuning, linear_detuning])
def test_sweep_endpoints_midpoint_and_symmetry(sweep):
    amplitude, duration = mhz_to_angular_frequency(50), 1e-6
    values = sweep(np.array([0, duration / 2, duration]), amplitude, duration)
    np.testing.assert_allclose(values, [-amplitude, 0, amplitude], rtol=1e-15, atol=0)
    times = np.linspace(0, duration, 101)
    np.testing.assert_allclose(sweep(times, amplitude, duration),
                               -sweep(duration - times, amplitude, duration), atol=1e-6, rtol=1e-14)


def test_profiles_have_distinct_quarter_time_values():
    assert cubic_detuning(0.25, 8, 1) == -1
    assert linear_detuning(0.25, 8, 1) == -4


@pytest.mark.parametrize("sweep", [cubic_detuning, linear_detuning])
@pytest.mark.parametrize("t,amplitude,duration", [(-1, 1, 1), (2, 1, 1), (0, -1, 1),
    (0, 1, 0), (0, 1, -1), (np.nan, 1, 1), (1j, 1, 1), (0, np.inf, 1)])
def test_sweep_invalid_input(sweep, t, amplitude, duration):
    with pytest.raises(ValueError):
        sweep(t, amplitude, duration)


@pytest.mark.parametrize("sweep", [cubic_detuning, linear_detuning])
def test_deterministic_repeatability_and_conservation(sweep):
    times = np.linspace(0, 1e-6, 501)
    omega, amplitude = mhz_to_angular_frequency(1), mhz_to_angular_frequency(50)
    callback = lambda t: rydberg_hamiltonian(omega, sweep(t, amplitude, times[-1]), 0, 10 * omega)
    first = evolve_trotter(KET_00, times, callback)
    second = evolve_trotter(KET_00, times, callback)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_allclose(np.linalg.norm(first, axis=1), 1, atol=1e-12, rtol=0)
    values = sweep_probabilities(first)
    np.testing.assert_allclose(values[:, :4].sum(axis=1), 1, atol=1e-12, rtol=0)
    np.testing.assert_allclose(values[:, 1], values[:, 2], atol=1e-12, rtol=0)
    assert values.min() >= 0 and values.max() <= 1 + 1e-12


@pytest.mark.parametrize("value,tolerance", [(0.9, 5e-5), (0.09, 5e-6), (0.009, 5e-7), (1, 5e-4)])
def test_significant_digits_are_not_fixed_decimal_places(value, tolerance):
    assert significant_tolerance(value) == pytest.approx(tolerance)


def test_comparison_identical_data_and_fidelity_alias():
    reference = np.array([[1, 0, 0, 0, 0], [0.2, 0.4, 0.4, 0, 0.8]])
    result = compare_probabilities(reference, reference)
    assert result["passes"]
    assert result["maximum_error_ratio"] == 0
    changed = reference.copy()
    changed[-1, 4] += 1e-3
    result = compare_probabilities(changed, reference)
    assert not result["passes"]
    assert result["maximum_P_Psi_plus_error"] == result["maximum_bell_fidelity_difference"]


def test_four_significant_digit_rounding_boundary_is_checked():
    reference = np.array([[1, 0, 0, 0, 0], [0.2, 0.4, 0.4, 0, 0.800049]])
    changed = reference.copy()
    changed[-1, 4] = 0.800051
    result = compare_probabilities(changed, reference)
    assert result["maximum_error_ratio"] < 1
    assert not result["scalar_rounding_matches"]
    assert not result["passes"]


def test_largest_passing_step_requires_all_finer_steps_to_pass():
    rows = [{"dt_s": dt, "passes": passed} for dt, passed in [(4, True), (1, True), (3, False), (2, True)]]
    assert largest_passing_timestep(rows) == 2
    assert largest_passing_timestep([{"dt_s": 1, "passes": False}]) is None


@pytest.mark.parametrize("rows", [[], [{"dt_s": 0, "passes": True}],
    [{"dt_s": 1, "passes": True}, {"dt_s": 1, "passes": False}]])
def test_invalid_convergence_steps(rows):
    with pytest.raises(ValueError):
        largest_passing_timestep(rows)


def test_comparison_rejects_mismatched_grids_and_nonfinite_values():
    with pytest.raises(ValueError):
        compare_probabilities(np.zeros((2, 5)), np.zeros((3, 5)))
    with pytest.raises(ValueError):
        compare_probabilities(np.full((2, 5), np.nan), np.zeros((2, 5)))


def test_global_reference_peak_includes_missed_sample_peak_error():
    values = np.array([[1, 0, 0, 0, 0], [0.2, 0.4, 0.4, 0, 0.8]])
    result = compare_probabilities(values, values, reference_peak_fidelity=0.81)
    assert result["peak_fidelity_error"] == pytest.approx(0.01)
    assert not result["passes"]
    with pytest.raises(ValueError):
        compare_probabilities(values, values, reference_peak_fidelity=0.79)
