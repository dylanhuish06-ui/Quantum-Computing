"""Left-endpoint exponential propagation for Chapter 5 Eq. (5.7)."""

import numpy as np
from scipy.linalg import expm

from .validation import _state_vector, is_hermitian, is_normalised


def evolve_trotter(psi0, time_grid, hamiltonian_function):
    """Propagate with exp(-i Hbar(t_n) (t_{n+1} - t_n)) at every interval.

    Parameters
    ----------
    psi0 : array-like, shape (4,)
        Finite normalised state at time_grid[0], in the fixed basis
        |00>, |0r>, |r0>, |rr>.
    time_grid : array-like, shape (n_times,)
        Nonempty, finite, real, strictly increasing times in seconds.
        Nonuniform intervals and nonzero starting times are supported.
    hamiltonian_function : callable
        Takes an absolute time in seconds and returns a finite Hermitian
        (4, 4) Hbar = H / ħ in rad/s. Called once per interval, at its LEFT
        endpoint. Any Hamiltonian parameter may depend on time.

    Returns
    -------
    numpy.ndarray, shape (n_times, 4)
        Complex states at every grid point, including an unchanged copy of
        psi0 in row zero. A one-point grid does not evaluate the callback.

    Notes
    -----
    Each full-matrix exponential uses scipy.linalg.expm. Later propagators
    multiply the state after earlier ones. No splitting into Hamiltonian
    terms, midpoint sampling, interpolation or renormalisation is performed.
    The method is exact for constant Hbar apart from floating-point error;
    for smooth time-varying Hbar it is a first-order global approximation.
    """
    state = _state_vector(psi0)
    if not is_normalised(state):
        raise ValueError("psi0 must be normalised; use normalise_state first")
    times = np.asarray(time_grid)
    if times.ndim != 1 or times.size == 0 or times.dtype.kind not in "iuf":
        raise ValueError("time_grid must be a nonempty one-dimensional array of finite real times")
    times = np.asarray(times, dtype=float)
    if not np.all(np.isfinite(times)):
        raise ValueError("time_grid must contain only finite real times")
    with np.errstate(over="ignore"):
        intervals = np.diff(times)
    if not np.all(np.isfinite(intervals)) or np.any(intervals <= 0):
        raise ValueError("time_grid must have finite, strictly positive intervals")
    if not callable(hamiltonian_function):
        raise ValueError("hamiltonian_function must be callable")

    states = np.empty((len(times), 4), dtype=complex)
    states[0] = state
    for index, dt in enumerate(intervals):
        matrix = np.asarray(hamiltonian_function(float(times[index])), dtype=complex)
        if not is_hermitian(matrix):
            raise ValueError(f"Hamiltonian at step {index} must be Hermitian")
        with np.errstate(over="ignore", invalid="ignore"):
            exponent = -1j * matrix * dt
        if not np.all(np.isfinite(exponent)):
            raise ValueError(f"Hamiltonian times interval at step {index} exceeds finite numeric range")
        states[index + 1] = expm(exponent) @ states[index]
    return states
