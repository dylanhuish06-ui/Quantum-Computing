"""Small-N adaptive grouping reference for generic static asymmetry.

Project H onto bulk-Dicke x distinguished-register spaces. This discards bulk
leakage and is APPROXIMATE until invariance holds. Promote the last bulk atom
until T ||(I-B B†) H B||_2 meets a requested error budget. Duhamel's formula
bounds the complex-state error by this value for every initial state in B.
This diagnostic deliberately uses full space and is guarded at N<=8. It is
not a scalable generic-disorder solver or a claim of efficient compression.
"""

import numpy as np

from .models import count, dicke_isometry_bulk
from .propagation import checked_matrix


def adaptive_static_grouping(H, N, duration, *, error_budget=1e-3, max_m=None):
    count(N, "N", 1)
    if N > 8:
        raise ValueError("adaptive full-space diagnostic limited to N<=8")
    if not np.isfinite(duration) or duration < 0:
        raise ValueError("duration must be finite and nonnegative")
    if not np.isfinite(error_budget) or error_budget <= 0:
        raise ValueError("error_budget must be positive and finite")
    max_m = N if max_m is None else count(max_m, "max_m")
    if max_m > N:
        raise ValueError("max_m cannot exceed N")
    H = checked_matrix(H, 2**N)
    history = []
    for m in range(max_m+1):
        B = dicke_isometry_bulk(N, m)
        reduced = (B.conj().T @ H @ B).tocsr()
        residual = H @ B - B @ reduced
        coupling = float(np.linalg.norm(residual.toarray(), 2))
        bound = duration*coupling
        history.append({"m": m, "dimension": B.shape[1],
                        "residual_operator_norm": coupling, "state_error_bound": bound})
        if bound <= error_budget:
            break
    return {"isometry": B, "hamiltonian": reduced, "history": history,
            "budget_met": bound <= error_budget,
            "exact_within_roundoff": coupling < 1e-13,
            "limitation": "static small-N full-space diagnostic; no large-N guarantee"}
