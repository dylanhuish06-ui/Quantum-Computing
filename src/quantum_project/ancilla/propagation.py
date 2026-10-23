"""Time-ordered DOP853 and second-order midpoint exponential action.

SciPy expm_multiply uses a scaled truncated Taylor exponential-action method,
not Lanczos. The independent dense expm backend is limited to small matrices.
Integration dt and measurement times are separate parameters. No unitary output
is renormalised. Split discontinuous profiles into separate evolve calls at
each discontinuity; DOP853 output times do not force integration boundaries.
"""

from collections.abc import Callable
from time import monotonic

import numpy as np
from numpy.typing import ArrayLike
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from .models import count, hybrid_hamiltonian


def checked_matrix(value, dimension):
    matrix = sparse.csr_matrix(value, dtype=complex)
    if matrix.shape != (dimension, dimension) or not np.all(np.isfinite(matrix.data)):
        raise ValueError("Hamiltonian must be finite with the state dimension")
    anti = matrix - matrix.conj().T
    if anti.nnz and np.max(abs(anti.data)) > 1e-12*max(1., np.max(abs(matrix.data))):
        raise ValueError("Hamiltonian must be Hermitian")
    return matrix


def time_grid(times):
    values = np.asarray(times)
    if values.ndim != 1 or not values.size or values.dtype.kind not in "iuf":
        raise ValueError("times must be a nonempty real vector")
    values = values.astype(float)
    if not np.all(np.isfinite(values)) or np.any(np.diff(values) <= 0):
        raise ValueError("times must be finite and strictly increasing")
    return values


def evolve(hamiltonian: Callable, initial: ArrayLike, times: ArrayLike, *,
           method: str = "dop853", dt: float = 0.01,
           rtol: float = 2e-12, atol: float = 2e-14,
           require_normalised: bool = True, max_evaluations: int = 200_000,
           max_seconds: float = 120.) -> np.ndarray:
    """Bounded propagation; dt is a maximum integration step, not an output grid.

    Runtime is checked between evaluations/actions, not an interrupt of a SciPy
    kernel. Exponential-action precision is SciPy's machine-precision target;
    it has no exposed Krylov tolerance. Do not report dt as that tolerance.
    """
    count(max_evaluations, "max_evaluations", 1)
    if not np.isfinite(max_seconds) or max_seconds <= 0:
        raise ValueError("max_seconds must be positive and finite")
    started = monotonic()
    times = time_grid(times)
    state = np.asarray(initial, dtype=complex).copy()
    if state.ndim != 1 or not state.size or not np.all(np.isfinite(state)):
        raise ValueError("initial must be a finite state vector")
    if require_normalised and abs(np.vdot(state, state)-1) > 1e-10:
        raise ValueError("initial must have unit norm")
    if method not in ("dop853", "midpoint", "dense_midpoint"):
        raise ValueError("unknown propagation method")
    if not np.isfinite(dt) or dt <= 0 or not (0 < rtol < 1 and 0 < atol < 1):
        raise ValueError("dt and solver tolerances must be positive and finite")
    d = len(state)
    if 3*16*len(times)*d > 256*2**20:
        raise ValueError("trajectory exceeds 256 MiB workspace budget")
    if method == "dense_midpoint" and d > 256:
        raise ValueError("dense exponential limited to dimension <=256")
    evaluations = 0
    def H(t):
        nonlocal evaluations
        evaluations += 1
        if evaluations > max_evaluations or monotonic()-started > max_seconds:
            raise RuntimeError("propagation resource budget exceeded")
        return checked_matrix(hamiltonian(float(t)), d)
    H(times[0])  # Validate even a zero-duration request.
    if len(times) == 1:
        return state[None, :]
    if method == "dop853":
        sol = solve_ivp(lambda t, y: -1j*(H(t) @ y), (times[0], times[-1]),
                        state, t_eval=times, method="DOP853", rtol=rtol,
                        atol=atol, max_step=dt)
        if not sol.success:
            raise RuntimeError(sol.message)
        result = sol.y.T
        if not np.all(np.isfinite(result)):
            raise RuntimeError("nonfinite propagated state")
        return result
    result = np.empty((len(times), d), dtype=complex)
    result[0] = state
    for i, (left, right) in enumerate(zip(times[:-1], times[1:])):
        steps = max(1, int(np.ceil((right-left)/dt - 1e-12)))
        if steps > max_evaluations-evaluations:
            raise RuntimeError("propagation resource budget exceeded")
        step = (right-left)/steps
        for j in range(steps):
            A = -1j*step*H(left+(j+0.5)*step)
            if method == "midpoint":
                state = expm_multiply(A, state, traceA=A.diagonal().sum())
            else:
                state = expm(A.toarray()) @ state
        result[i+1] = state
    if not np.all(np.isfinite(result)):
        raise RuntimeError("nonfinite propagated state")
    return result


def evolve_hybrid(profile, initial, times, **kwargs):
    """profile(t) returns HybridParameters; its bulk symmetry is structural."""
    return evolve(lambda t: hybrid_hamiltonian(profile(t)), initial, times, **kwargs)
