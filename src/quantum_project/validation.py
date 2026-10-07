"""Validation for the fixed four-dimensional two-atom Hilbert space."""

from numbers import Real

import numpy as np

ATOL = 1e-12


def _real_scalar(value, name):
    """Validate a finite real scalar without silently discarding imaginary parts."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be a finite real scalar")
    return value


def _state_vector(psi):
    """Return a finite complex vector of shape (4,), or raise ValueError."""
    state = np.asarray(psi, dtype=complex)
    if state.shape != (4,) or not np.all(np.isfinite(state)):
        raise ValueError("state must be a finite vector of shape (4,)")
    return state


def is_hermitian(hamiltonian):
    """Check a finite 4 × 4 matrix against its conjugate transpose."""
    matrix = np.asarray(hamiltonian, dtype=complex)
    if matrix.shape != (4, 4) or not np.all(np.isfinite(matrix)):
        raise ValueError("Hamiltonian must be a finite matrix of shape (4, 4)")
    return bool(np.allclose(matrix, matrix.conj().T, atol=ATOL, rtol=0))


def is_normalised(psi):
    """Check whether a finite length-four state's norm is one."""
    state = _state_vector(psi)
    # An amplitude larger than one cannot belong to a normalised state;
    # this guard also avoids overflow when checking very large inputs.
    if np.max(np.abs(state)) > 1 + ATOL:
        return False
    return bool(np.isclose(np.vdot(state, state).real, 1, atol=ATOL, rtol=0))


def are_orthonormal(states):
    """Check 1–4 normalised states, supplied as rows, for mutual orthogonality."""
    vectors = np.asarray(states, dtype=complex)
    if (vectors.ndim != 2 or vectors.shape[1] != 4
            or not 1 <= vectors.shape[0] <= 4
            or not np.all(np.isfinite(vectors))):
        raise ValueError("states must be a finite array of shape (n, 4), with 1 <= n <= 4")
    if not all(is_normalised(state) for state in vectors):
        return False
    gram = vectors.conj() @ vectors.T
    return bool(np.allclose(gram, np.eye(len(vectors)), atol=ATOL, rtol=0))
