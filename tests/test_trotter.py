"""Validate constant-H propagation before the time-dependent examples."""

from unittest.mock import Mock

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from quantum_project import (
    KET_00, PSI_PLUS, evolve_static, evolve_trotter, normalise_state, rydberg_hamiltonian,
)
from quantum_project.units import mhz_to_angular_frequency
from quantum_project.static_limits import noninteracting_state


OMEGA = mhz_to_angular_frequency(1)
COMPLEX_STATE = normalise_state([1 + 2j, -3j, 4, -5 + 1j])


def test_constant_hamiltonian_refinement_study():
    maximum_error = 0.0
    maximum_norm_error = 0.0
    print("Constant H: dt (ns), maximum amplitude discrepancy, maximum norm deviation")
    for steps in (8, 32, 128, 512):
        times = np.linspace(0, 1e-6, steps + 1)
        grid_error = grid_norm_error = 0.0
        for delta_mhz, interaction_mhz, phase in ((0, 0, 0), (0.4, 5, 0.73), (-0.4, 100, np.pi / 2)):
            matrix = rydberg_hamiltonian(
                OMEGA, mhz_to_angular_frequency(delta_mhz), phase,
                mhz_to_angular_frequency(interaction_mhz),
            )
            for initial in (KET_00, PSI_PLUS, COMPLEX_STATE):
                actual = evolve_trotter(initial, times, lambda t: matrix)
                expected = evolve_static(initial, times, matrix)
                grid_error = max(grid_error, float(np.max(np.abs(actual - expected))))
                grid_norm_error = max(grid_norm_error, float(np.max(np.abs(np.linalg.norm(actual, axis=1) - 1))))
                np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=0)
                np.testing.assert_allclose(np.linalg.norm(actual, axis=1), 1, atol=1e-12, rtol=0)
                np.testing.assert_array_equal(actual[0], initial)
        print(f"{1e3 / steps:10.6f} {grid_error:.6e} {grid_norm_error:.6e}")
        maximum_error = max(maximum_error, grid_error)
        maximum_norm_error = max(maximum_norm_error, grid_norm_error)
    # No monotonic decrease is required: for constant H there is only roundoff,
    # and refining the grid increases the number of floating-point operations.
    print(f"Maximum constant-H amplitude discrepancy: {maximum_error:.6e}")
    print(f"Maximum constant-H norm deviation: {maximum_norm_error:.6e}")


@pytest.mark.parametrize("start", [0.0, 2e-6, -1e-6])
def test_one_point_grid_returns_initial_state_without_callback(start):
    callback = Mock(side_effect=AssertionError("no interval to propagate"))
    actual = evolve_trotter(COMPLEX_STATE, [start], callback)
    assert actual.shape == (1, 4)
    assert np.issubdtype(actual.dtype, np.complexfloating)
    np.testing.assert_array_equal(actual[0], COMPLEX_STATE)
    assert not np.shares_memory(actual, COMPLEX_STATE)
    callback.assert_not_called()


def test_nonuniform_grid_and_absolute_start_time_with_constant_h():
    times = np.array([2e-6, 2.01e-6, 2.2e-6, 2.7e-6])
    matrix = rydberg_hamiltonian(OMEGA, 0.4 * OMEGA, 0.73, 5 * OMEGA)
    callback = Mock(return_value=matrix)
    actual = evolve_trotter(COMPLEX_STATE, times, callback)
    expected = evolve_static(COMPLEX_STATE, times - times[0], matrix)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=0)
    assert callback.call_count == len(times) - 1
    assert [call.args[0] for call in callback.call_args_list] == list(times[:-1])


def test_left_endpoint_sampling_and_chronological_noncommuting_propagation():
    times = np.array([3e-6, 3.2e-6, 3.65e-6])
    drive = rydberg_hamiltonian(OMEGA, 0, 0, 0)
    delta = 0.7 * OMEGA
    detuning = rydberg_hamiltonian(0, delta, 0, 0)
    callback = Mock(side_effect=lambda t: drive if t < times[1] else detuning)
    actual = evolve_trotter(KET_00, times, callback)
    # Independent closed forms: drive both atoms, THEN apply diagonal phases.
    after_drive = noninteracting_state([times[1] - times[0]], OMEGA)[0]
    diagonal_phases = np.exp(-1j * np.array([delta, 0, 0, -delta]) * (times[2] - times[1]))
    expected = diagonal_phases * after_drive
    np.testing.assert_allclose(actual[1], after_drive, atol=1e-14, rtol=0)
    np.testing.assert_allclose(actual[2], expected, atol=1e-14, rtol=0)
    reverse_order = diagonal_phases[0] * after_drive
    assert np.linalg.norm(expected - reverse_order) > 0.1
    assert [call.args[0] for call in callback.call_args_list] == list(times[:-1])


def test_commuting_linear_detuning_against_exact_integrated_phase():
    duration = 1e-6
    errors = []
    for steps in (20, 40, 80, 160):
        times = np.linspace(0, duration, steps + 1)
        actual = evolve_trotter(
            COMPLEX_STATE, times,
            lambda t: rydberg_hamiltonian(0, OMEGA * t / duration, 0, 0),
        )
        integrated_delta = OMEGA * times**2 / (2 * duration)
        exact = COMPLEX_STATE * np.exp(-1j * integrated_delta[:, None] * np.array([1, 0, 0, -1]))
        errors.append(float(np.max(np.abs(actual - exact))))
        np.testing.assert_allclose(np.linalg.norm(actual, axis=1), 1, atol=1e-13, rtol=0)
    ratios = np.array(errors[:-1]) / errors[1:]
    np.testing.assert_allclose(ratios, 2, atol=0.01, rtol=0)


def test_detuning_ramp_converges_first_order_against_independent_ode():
    duration = 1e-6

    def hbar(t):
        return rydberg_hamiltonian(OMEGA, OMEGA * (t / duration - 0.5), 0.37, 3 * OMEGA)

    # This reference integrates the differential equation, not step exponentials.
    reference = solve_ivp(
        lambda t, psi: -1j * hbar(t) @ psi,
        (0, duration), KET_00.copy(), method="DOP853", dense_output=True,
        rtol=1e-12, atol=1e-14, max_step=duration / 100,
    )
    assert reference.success
    errors = []
    print("Linear detuning ramp: dt (ns), maximum amplitude discrepancy versus DOP853")
    for steps in (20, 40, 80, 160, 320):
        times = np.linspace(0, duration, steps + 1)
        actual = evolve_trotter(KET_00, times, hbar)
        errors.append(float(np.max(np.abs(actual - reference.sol(times).T))))
        np.testing.assert_allclose(np.linalg.norm(actual, axis=1), 1, atol=1e-12, rtol=0)
        print(f"{1e3 / steps:10.6f} {errors[-1]:.6e}")
    ratios = np.array(errors[:-1]) / errors[1:]
    assert np.all((ratios > 1.8) & (ratios < 2.2))
    assert errors[-1] < 0.01


def test_all_hamiltonian_parameters_can_vary_with_time():
    times = np.linspace(0, 1e-6, 51)
    seen = []

    def callback(t):
        x = t / times[-1]
        seen.append(t)
        return rydberg_hamiltonian(OMEGA * (1 + x), OMEGA * x, 0.7 * x, OMEGA * (5 - x))

    states = evolve_trotter(COMPLEX_STATE, times, callback)
    assert seen == list(times[:-1])
    np.testing.assert_allclose(np.linalg.norm(states, axis=1), 1, atol=1e-12, rtol=0)


def test_inputs_are_not_modified():
    initial = COMPLEX_STATE.copy()
    times = np.linspace(0, 1e-6, 11)
    matrix = rydberg_hamiltonian(OMEGA, 0, 0.3, OMEGA)
    snapshots = [value.copy() for value in (initial, times, matrix)]
    evolve_trotter(initial, times, lambda t: matrix)
    for value, snapshot in zip((initial, times, matrix), snapshots):
        np.testing.assert_array_equal(value, snapshot)


@pytest.mark.parametrize("invalid", [
    [1, 0], np.ones((4, 1)), [2, 0, 0, 0], [0, 0, 0, 0],
    [np.nan, 0, 0, 0], [np.inf, 0, 0, 0],
])
def test_invalid_initial_state(invalid):
    with pytest.raises(ValueError):
        evolve_trotter(invalid, [0, 1e-6], lambda t: np.eye(4))


@pytest.mark.parametrize("invalid", [
    [], 0, [[0, 1]], [0, np.nan], [0, np.inf], [0, -np.inf], [0j, 1j],
    [False, True], ["0", "1"], [0, 0], [1, 0], [0, 2, 1], [-1e308, 1e308],
])
def test_invalid_time_grid(invalid):
    callback = Mock(return_value=np.eye(4))
    with pytest.raises(ValueError, match="time_grid"):
        evolve_trotter(KET_00, invalid, callback)
    callback.assert_not_called()


@pytest.mark.parametrize("invalid", [
    np.eye(3), np.zeros(4), np.full((4, 4), np.nan), np.full((4, 4), np.inf),
    1j * np.eye(4), np.triu(np.ones((4, 4))),
])
def test_every_callback_result_is_validated(invalid):
    callback = Mock(side_effect=[np.eye(4), invalid])
    with pytest.raises(ValueError):
        evolve_trotter(KET_00, [0, 1e-6, 2e-6], callback)
    assert callback.call_count == 2


def test_noncallable_hamiltonian_is_rejected():
    with pytest.raises(ValueError, match="callable"):
        evolve_trotter(KET_00, [0, 1], np.eye(4))


def test_unrepresentable_step_exponent_is_rejected():
    with pytest.raises(ValueError, match="finite numeric range"):
        evolve_trotter(KET_00, [0, 2], lambda t: 1e308 * np.eye(4))
