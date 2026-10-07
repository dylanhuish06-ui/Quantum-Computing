"""Constant common laser phase is a gauge convention for this initial state."""

import numpy as np
import pytest

from quantum_project import KET_00, evolve_static, evolve_trotter, rydberg_hamiltonian
from quantum_project.convergence import sweep_probabilities
from quantum_project.sweeps import cubic_detuning, linear_detuning
from quantum_project.units import mhz_to_angular_frequency


@pytest.mark.parametrize("interaction_mhz", [0, 0.5, 1, 5, 100])
@pytest.mark.parametrize("profile", ["static", "cubic", "linear"])
def test_constant_phase_invariance(interaction_mhz, profile):
    times = np.linspace(0, 1e-6, 1025)
    omega = mhz_to_angular_frequency(1)
    interaction = mhz_to_angular_frequency(interaction_mhz)
    amplitude = mhz_to_angular_frequency(50)
    sweep = {"cubic": cubic_detuning, "linear": linear_detuning}.get(profile)

    def callback(phase):
        return lambda t: rydberg_hamiltonian(
            omega, 0 if sweep is None else sweep(t, amplitude, times[-1]), phase, interaction)

    def propagate(phase):
        h = callback(phase)
        return (evolve_static(KET_00, times, h(0)) if sweep is None
                else evolve_trotter(KET_00, times, h))

    zero_states = propagate(0)
    zero_probabilities = sweep_probabilities(zero_states)
    for phase in (0.37, np.pi / 2, -1.2):
        diagonal = np.exp(1j * phase * np.array([0, 1, 1, 2]))
        transform = np.diag(diagonal)
        for t in times[::256]:
            actual = callback(phase)(t)
            expected = transform @ callback(0)(t) @ transform.conj().T
            np.testing.assert_allclose(actual, expected, atol=2e-7, rtol=2e-15)
            np.testing.assert_allclose(actual, actual.conj().T, atol=0, rtol=0)
        states = propagate(phase)
        np.testing.assert_allclose(states, zero_states * diagonal, atol=2e-12, rtol=0)
        probabilities = sweep_probabilities(states)
        np.testing.assert_allclose(probabilities, zero_probabilities, atol=2e-12, rtol=0)
        np.testing.assert_allclose(probabilities[:, :4].sum(axis=1), 1, atol=2e-12, rtol=0)
        np.testing.assert_allclose(probabilities[:, 1], probabilities[:, 2], atol=2e-12, rtol=0)
        errors = np.max(np.abs(probabilities - zero_probabilities), axis=0)
        print(f"phase invariance: {profile} V={interaction_mhz:g} MHz phi={phase:.8g} "
              f"P00_error={errors[0]:.9e} P_Psi_plus_error={errors[4]:.9e}")
