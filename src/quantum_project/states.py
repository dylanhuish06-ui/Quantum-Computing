"""Canonical states ordered as |00>, |0r>, |r0>, |rr>."""

import numpy as np

from .validation import _state_vector

BASIS_LABELS = ("|00>", "|0r>", "|r0>", "|rr>")
BASIS_STATES = np.eye(4, dtype=complex)
BASIS_STATES.setflags(write=False)
KET_00, KET_0R, KET_R0, KET_RR = BASIS_STATES
PSI_PLUS = (KET_0R + KET_R0) / np.sqrt(2)
PSI_PLUS.setflags(write=False)


def normalise_state(psi):
    """Return a new unit vector; reject malformed, nonfinite or zero states."""
    state = _state_vector(psi)
    # Scale first so nonzero inputs with very large/small amplitudes work too.
    scale = max(np.max(np.abs(state.real)), np.max(np.abs(state.imag)))
    if scale == 0:
        raise ValueError("cannot normalise the zero state")
    scaled = state / scale
    return scaled / np.linalg.norm(scaled)
