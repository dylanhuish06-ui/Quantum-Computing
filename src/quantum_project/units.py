"""Explicit scalar conversions: ordinary MHz, angular rad/s, and seconds."""

import numpy as np

from .validation import _real_scalar


def mhz_to_angular_frequency(frequency_mhz):
    """Convert ordinary frequency in MHz to angular frequency in rad/s."""
    result = _real_scalar(frequency_mhz, "frequency_mhz") * (2 * np.pi * 1e6)
    if not np.isfinite(result):
        raise ValueError("converted angular frequency is outside the finite numeric range")
    return result


def angular_frequency_to_mhz(angular_frequency):
    """Convert angular frequency in rad/s to ordinary frequency in MHz."""
    return _real_scalar(angular_frequency, "angular_frequency") / (2 * np.pi * 1e6)


def microseconds_to_seconds(time_us):
    """Convert a time or signed time interval in microseconds to seconds."""
    return _real_scalar(time_us, "time_us") * 1e-6
