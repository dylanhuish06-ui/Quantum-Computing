"""Small-system spin-gap and properly ordered permutation control references."""

from itertools import permutations

import numpy as np
from scipy import sparse

from ..symmetry import local_operator
from .models import count, reference_guard
from .propagation import evolve


def spin_squared(N):
    reference_guard(N)
    count(N, "N", 1)
    axes = ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])
    result = sparse.csr_matrix((2**N, 2**N), dtype=complex)
    for axis in axes:
        J = sum((.5*local_operator(axis, i, N) for i in range(N)),
                sparse.csr_matrix(result.shape, dtype=complex))
        result += J@J
    result.eliminate_zeros()
    return result


def gap_penalty(N, G):
    reference_guard(N)
    count(N, "N", 1)
    if not np.isfinite(G) or G < 0:
        raise ValueError("G must be finite and nonnegative")
    j = N/2
    return G*(j*(j+1)*sparse.eye(2**N, format="csr")-spin_squared(N))


def permutation_operator(permutation):
    p = tuple(permutation)
    N = len(p)
    reference_guard(N)
    if sorted(p) != list(range(N)):
        raise ValueError("invalid permutation")
    cols = np.arange(2**N)
    rows = np.zeros(2**N, dtype=int)
    for i, old in enumerate(p):
        rows |= ((cols >> (N-1-old)) & 1) << (N-1-i)
    return sparse.csr_matrix((np.ones(2**N), (rows, cols)), shape=(2**N, 2**N))


def symmetric_group(N):
    count(N, "N", 1)
    if N > 5:
        raise ValueError("explicit S_N averaging limited to N<=5")
    return list(permutations(range(N)))


def permutation_average(H, group):
    operators = [permutation_operator(p) for p in group]
    if not operators:
        raise ValueError("group cannot be empty")
    return sum((U.conj().T@H@U for U in operators),
               sparse.csr_matrix(H.shape, dtype=complex))/len(operators)


def swap_decomposition(permutation):
    """Ordered position swaps taking identity ordering to the given tuple."""
    current = list(range(len(permutation)))
    swaps = []
    for i, value in enumerate(permutation):
        j = current.index(value)
        if i != j:
            swaps.append((i, j))
            current[i], current[j] = current[j], current[i]
    return swaps


def permutation_cycles(H, initial, total_time, cycles, group, *,
                       dt=.01, pulse_duration=0., pulse_substeps=4, method="dop853"):
    """Apply each U, evolve freely, and undo U, in the supplied group order.

    pulse_duration=0 gives instantaneous permutations. For positive duration,
    each swap uses H_control=pi/(2*t_p)*(I-SWAP), while H(t) remains on.
    Total_time includes pulses; free windows shrink accordingly.
    This deliberately unoptimised compilation uses <=N-1 swaps per permutation
    in each direction, and three CNOTs per ideal SWAP as a digital cost model.
    """
    count(cycles, "cycles", 1)
    count(pulse_substeps, "pulse_substeps", 1)
    if (not group or not np.isfinite(total_time) or not np.isfinite(pulse_duration)
            or total_time <= 0 or pulse_duration < 0):
        raise ValueError("invalid cycle parameters")
    operators = [permutation_operator(p) for p in group]
    swaps = [swap_decomposition(p) for p in group]
    swaps_per_cycle = 2*sum(map(len, swaps))
    free = (total_time/cycles-swaps_per_cycle*pulse_duration)/len(group)
    if free <= 0:
        raise ValueError("pulse windows consume all available time")
    state = np.asarray(initial, complex).copy()
    # Validate before any pulses.
    evolve(H, state, [0])
    t = 0.
    histories, times = [state.copy()], [0.]
    for _ in range(cycles):
        for p, U, pairs in zip(group, operators, swaps):
            if pulse_duration == 0:
                state = U@state
                state = evolve(H, state, [t, t+free], method=method, dt=dt,
                               require_normalised=False)[-1]
                state = U.conj().T@state
                t += free
            else:
                def pulse(a, b, state, start):
                    order = list(range(len(p)))
                    order[a], order[b] = order[b], order[a]
                    S = permutation_operator(order)
                    control = np.pi/(2*pulse_duration)*(sparse.eye(S.shape[0])-S)
                    return evolve(lambda s: H(s)+control, state, [start, start+pulse_duration],
                                  dt=min(dt, pulse_duration/pulse_substeps), method=method,
                                  require_normalised=False)[-1]
                for a, b in pairs:
                    state = pulse(a, b, state, t)
                    t += pulse_duration
                state = evolve(H, state, [t, t+free], method=method, dt=dt,
                               require_normalised=False)[-1]
                t += free
                for a, b in reversed(pairs):
                    state = pulse(a, b, state, t)
                    t += pulse_duration
        histories.append(state.copy())
        times.append(t)
    return {"times": np.array(times), "states": np.array(histories),
            "free_segments": cycles*len(group), "swap_count": cycles*swaps_per_cycle,
            "cnot_equivalent_count": 3*cycles*swaps_per_cycle,
            "serial_cnot_depth_upper_bound": 3*cycles*swaps_per_cycle,
            "pulse_duration": pulse_duration,
            "pulse_control_norm": None if pulse_duration == 0 else np.pi/pulse_duration}
