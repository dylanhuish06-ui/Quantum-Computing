import numpy as np
import pytest

from quantum_project import BASIS_STATES, KET_00, KET_0R
from quantum_project.units import (
    angular_frequency_to_mhz, mhz_to_angular_frequency, microseconds_to_seconds,
)
from quantum_project.validation import are_orthonormal, is_hermitian, is_normalised


def test_validation_detects_valid_shapes_with_invalid_physics():
    matrix = np.eye(4, dtype=complex)
    matrix[0, 1] = 1j
    assert not is_hermitian(matrix)
    assert not is_hermitian(1j * np.eye(4))
    assert not is_normalised(2 * KET_00)
    assert not is_normalised(np.zeros(4))
    assert not is_normalised([1e300, 0, 0, 0])
    assert not are_orthonormal([KET_00, KET_00])
    assert not are_orthonormal(2 * BASIS_STATES)
    assert are_orthonormal([KET_00, 1j * KET_0R])


def test_normalisation_does_not_use_loose_relative_tolerance():
    assert not is_normalised((1 + 1e-8) * KET_00)


@pytest.mark.parametrize("invalid", [np.eye(3), np.zeros(4), np.full((4, 4), np.nan)])
def test_hermitian_validation_rejects_malformed_input(invalid):
    with pytest.raises(ValueError):
        is_hermitian(invalid)


@pytest.mark.parametrize("invalid", [np.eye(4), np.zeros(3), [np.inf, 0, 0, 0]])
def test_normalisation_validation_rejects_malformed_input(invalid):
    with pytest.raises(ValueError):
        is_normalised(invalid)


@pytest.mark.parametrize("invalid", [[], np.zeros(4), np.eye(3),
                                       np.zeros((5, 4)), np.full((4, 4), np.nan)])
def test_orthonormal_validation_rejects_malformed_input(invalid):
    with pytest.raises(ValueError):
        are_orthonormal(invalid)


def test_known_unit_conversions():
    assert mhz_to_angular_frequency(1) == pytest.approx(2 * np.pi * 1e6)
    assert angular_frequency_to_mhz(2 * np.pi * 1e6) == pytest.approx(1)
    assert microseconds_to_seconds(2.5) == pytest.approx(2.5e-6)


@pytest.mark.parametrize("frequency", [0, -3.2, 0.125, 100])
def test_frequency_round_trip(frequency):
    assert angular_frequency_to_mhz(mhz_to_angular_frequency(frequency)) == pytest.approx(frequency)


@pytest.mark.parametrize("convert", [mhz_to_angular_frequency,
                                       angular_frequency_to_mhz, microseconds_to_seconds])
@pytest.mark.parametrize("invalid", [np.nan, np.inf, 1j, [1], "1", True])
def test_conversions_reject_invalid_scalars(convert, invalid):
    with pytest.raises(ValueError):
        convert(invalid)


def test_conversion_rejects_overflow():
    with pytest.raises(ValueError, match="finite numeric range"):
        mhz_to_angular_frequency(1e308)
