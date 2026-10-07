"""Chapter 5 Eq. (5.4), represented as Hbar = H / ħ in rad/s."""

import numpy as np

from .validation import _real_scalar


def rydberg_hamiltonian(omega, delta, phase, V):
    """Return H/ħ in the basis |00>, |0r>, |r0>, |rr>.

    omega, delta and V are finite real angular frequencies (rad/s).
    phase is a finite real angle (radians). All signs are accepted.
    Evolution uses exp(-i Hbar t) with time t in seconds.
    """
    omega = _real_scalar(omega, "omega")
    delta = _real_scalar(delta, "delta")
    phase = _real_scalar(phase, "phase")
    V = _real_scalar(V, "V")
    coupling = (omega / 2) * np.exp(-1j * phase)
    hbar = np.array([
        [delta, coupling, coupling, 0],
        [coupling.conjugate(), 0, 0, coupling],
        [coupling.conjugate(), 0, 0, coupling],
        [0, coupling.conjugate(), coupling.conjugate(), -delta + V],
    ], dtype=complex)
    if not np.all(np.isfinite(hbar)):
        raise ValueError("parameters produce a Hamiltonian outside the finite numeric range")
    return hbar
