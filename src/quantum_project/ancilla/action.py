"""Sparse exponential action with an explicit Taylor truncation budget.

For Hermitian H, the integral Taylor remainder of exp(-i h H)v is bounded
by ||v|| (|h| ||H||_2)^(p+1)/(p+1)!. Use the induced 1-norm as a bound
on ||H||_2, shift its real trace, and subdivide until |h| ||H||_1 <= 1.
The reported bound excludes floating-point roundoff and time-ordering error.
This is a Taylor method, not a Krylov tolerance for scipy.expm_multiply.
"""

import numpy as np
from scipy import sparse

from .models import count
from .propagation import checked_matrix


def controlled_action(H, state, duration, *, atol=1e-12, max_matvecs=100_000):
    """Return (state, diagnostics), without normalisation or dense conversion."""
    state = np.asarray(state, dtype=complex).copy()
    if state.ndim != 1 or not state.size or not np.all(np.isfinite(state)):
        raise ValueError("state must be a finite vector")
    if not np.isfinite(duration) or not np.isfinite(atol) or atol <= 0:
        raise ValueError("duration must be finite and atol positive")
    count(max_matvecs, "max_matvecs", 1)
    H = checked_matrix(H, state.size)
    shift = float(H.diagonal().real.mean())
    centered = H - sparse.eye(state.size, format="csr")*shift
    radius = float(np.max(np.asarray(abs(centered).sum(axis=0))))
    action_size = abs(duration)*radius
    if not np.isfinite(action_size) or action_size > max_matvecs:
        raise RuntimeError("exponential action exceeds matvec budget")
    steps = max(1, int(np.ceil(action_size)))
    step = duration/steps
    a = abs(step)*radius
    bound, matvecs = 0., 0
    for _ in range(steps):
        term, total = state.copy(), state.copy()
        remainder = float(np.linalg.norm(state))*a
        degree = 0
        while remainder > atol/steps:
            if matvecs >= max_matvecs:
                raise RuntimeError("exponential action exceeds matvec budget")
            degree += 1
            term = (-1j*step/degree)*(centered @ term)
            total += term
            matvecs += 1
            remainder *= a/(degree+1)
        bound += remainder
        state = np.exp(-1j*step*shift)*total
    if not np.all(np.isfinite(state)):
        raise RuntimeError("nonfinite exponential action")
    return state, {"truncation_bound": bound, "matvecs": matvecs,
                   "substeps": steps, "roundoff_included": False}
