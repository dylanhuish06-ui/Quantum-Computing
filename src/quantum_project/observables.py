"""Pure-state probabilities and fidelities for normalised two-atom states."""

import numpy as np

from .validation import _state_vector, is_normalised


def fidelity(psi, target_state):
    """Return |<target_state|psi>|²; both inputs must already be normalised."""
    state = _state_vector(psi)
    target = _state_vector(target_state)
    if not is_normalised(state) or not is_normalised(target):
        raise ValueError("fidelity requires normalised states; use normalise_state first")
    # Clip only rounding error at the boundaries of the probability interval.
    return float(np.clip(abs(np.vdot(target, state)) ** 2, 0, 1))


def state_probability(psi, basis_state):
    """Return the projection probability onto a normalised basis/target state."""
    return fidelity(psi, basis_state)
