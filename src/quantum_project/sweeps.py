"""Detuning sweeps in angular-frequency units (rad/s), with time in seconds."""

import numpy as np

from .validation import _real_scalar


def _sweep_coordinate(t, delta_max, tmax):
    delta_max = _real_scalar(delta_max, "delta_max")
    tmax = _real_scalar(tmax, "tmax")
    if tmax <= 0 or delta_max < 0:
        raise ValueError("tmax must be positive and delta_max must be nonnegative")
    times = np.asarray(t)
    if times.dtype.kind not in "iuf" or not np.all(np.isfinite(times)):
        raise ValueError("t must contain finite real seconds")
    if np.any(times < 0) or np.any(times > tmax):
        raise ValueError("sweep times must lie in [0, tmax]")
    return (times - tmax / 2) / tmax, delta_max


def cubic_detuning(t, delta_max, tmax):
    """Return exactly 8 delta_max ((t - tmax/2) / tmax)^3 in rad/s."""
    coordinate, amplitude = _sweep_coordinate(t, delta_max, tmax)
    return 8 * amplitude * coordinate**3


def linear_detuning(t, delta_max, tmax):
    """Return 2 delta_max (t/tmax - 1/2) in rad/s."""
    coordinate, amplitude = _sweep_coordinate(t, delta_max, tmax)
    return 2 * amplitude * coordinate
