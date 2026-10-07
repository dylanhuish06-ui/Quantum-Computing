"""Derived resonant limits for |00> initial state and zero laser phase.

These closed-form states are independent of the numerical eigensolver.
See docs/static_validation.md for the derivation and validity conditions.
"""

import numpy as np

from .validation import _real_scalar


def _rabi_angles(times, omega):
    omega = _real_scalar(omega, "omega")
    times = np.asarray(times)
    if times.ndim != 1 or times.dtype.kind not in "iuf" or not np.all(np.isfinite(times)):
        raise ValueError("times must be a one-dimensional array of finite real seconds")
    return np.asarray(times, dtype=float) * omega


def noninteracting_state(times, omega):
    """Exact state for V = delta = phase = 0, in the fixed four-state basis.

    times is a 1D array in seconds; omega is a real angular frequency in rad/s.
    Return shape (len(times), 4), starting from |00> at t = 0.
    """
    angle = _rabi_angles(times, omega) / 2
    c, s = np.cos(angle), np.sin(angle)
    return np.column_stack((c**2, -1j * c * s, -1j * c * s, -s**2))


def blockade_state(times, omega):
    """Leading |V/omega| -> infinity limit for delta = phase = 0.

    Starting from |00>, the coupled states are |00> and |Psi+>, with matrix
    element omega/sqrt(2). Finite-V light shifts and leakage are omitted.
    Units and output shape match noninteracting_state.
    """
    angle = _rabi_angles(times, omega) / np.sqrt(2)
    c, s = np.cos(angle), np.sin(angle)
    single = -1j * s / np.sqrt(2)
    return np.column_stack((c, single, single, np.zeros_like(c)))
