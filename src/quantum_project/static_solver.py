"""Time-independent evolution in the basis |00>, |0r>, |r0>, |rr>."""

import numpy as np

from .validation import _state_vector, is_hermitian, is_normalised


def evolve_static(psi0, times, Hbar):
    """Evolve a normalised state with a constant Hermitian Hbar = H / ħ.

    Parameters
    ----------
    psi0 : array-like, shape (4,)
        Finite, normalised initial state at t = 0.
    times : array-like, shape (n_times,)
        Finite real times in seconds, relative to t = 0. Negative, repeated,
        unsorted and empty time arrays are supported.
    Hbar : array-like, shape (4, 4)
        Finite Hermitian Hamiltonian in rad/s, in the fixed two-atom basis.

    Returns
    -------
    numpy.ndarray, shape (n_times, 4)
        Complex state amplitudes, one row per requested time. No input or
        output state is silently renormalised.

    Notes
    -----
    numpy.linalg.eigh returns eigenvectors as columns of R. Hence the initial
    eigenbasis coefficients are R† @ psi0. Each coefficient gets the phase
    exp(-i lambda_j t), then R transforms the state back. With states stored
    as rows, this last operation is right multiplication by R.T, not R†.
    Exactly one diagonalisation is performed per call.
    """
    state = _state_vector(psi0)
    if not is_normalised(state):
        raise ValueError("psi0 must be normalised; use normalise_state first")

    time_array = np.asarray(times)
    if time_array.ndim != 1 or time_array.dtype.kind not in "iuf":
        raise ValueError("times must be a one-dimensional array of finite real numbers")
    time_array = np.asarray(time_array, dtype=float)
    if not np.all(np.isfinite(time_array)):
        raise ValueError("times must be a one-dimensional array of finite real numbers")

    matrix = np.asarray(Hbar, dtype=complex)
    if not is_hermitian(matrix):
        raise ValueError("Hbar must be Hermitian")

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    coefficients = eigenvectors.conj().T @ state
    phases = np.exp(-1j * time_array[:, None] * eigenvalues[None, :])
    return (phases * coefficients[None, :]) @ eigenvectors.T
