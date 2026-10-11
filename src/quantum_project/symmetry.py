"""Full tensor-product and symmetric Dicke models; Hbar is always in rad/s.

Computational ordering is binary, with the leftmost atom the most significant
bit and |0>, |r> represented by 0, 1. No milestone solver is modified.
"""

from math import comb
from numbers import Integral

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from .hamiltonian import rydberg_hamiltonian
from .validation import _real_scalar


def atom_count(N):
    if isinstance(N, (bool, np.bool_)) or not isinstance(N, Integral) or N < 1:
        raise ValueError("N must be a positive integer")
    return int(N)


def exchange_transformation():
    """Columns: |00>, |Psi+>, |Psi->, |rr> in the milestone basis."""
    s = 1 / np.sqrt(2)
    return np.array([[1, 0, 0, 0], [0, s, s, 0], [0, s, -s, 0],
                     [0, 0, 0, 1]], dtype=complex)


def exchange_hamiltonian(omega, delta, phase, V):
    transform = exchange_transformation()
    return transform.conj().T @ rydberg_hamiltonian(omega, delta, phase, V) @ transform


def reduced_two_atom_hamiltonian(omega, delta, phase, V):
    """Direct symmetric construction, independently checked against U† H U."""
    return dicke_hamiltonian(2, omega, delta, phase, V).toarray()


def local_operator(operator, site, N):
    """Tensor embedding I ⊗ ... ⊗ operator ⊗ ... ⊗ I, in CSR storage."""
    N = atom_count(N)
    if not isinstance(site, Integral) or not 0 <= site < N:
        raise ValueError("site must be an integer in [0,N)")
    operator = np.asarray(operator, dtype=complex)
    if operator.shape != (2, 2) or not np.all(np.isfinite(operator)):
        raise ValueError("local operator must be a finite 2x2 matrix")
    return sparse.kron(sparse.eye(2**site, format="csr"),
                       sparse.kron(sparse.csr_matrix(operator),
                                   sparse.eye(2**(N-site-1), format="csr"), format="csr"),
                       format="csr")


def pair_matrix(N, V):
    N = atom_count(N)
    if np.ndim(V) == 0:
        value = _real_scalar(V, "V")
        return value * (np.ones((N, N)) - np.eye(N))
    values = np.asarray(V)
    if (values.shape != (N, N) or values.dtype.kind not in "iuf"
            or not np.all(np.isfinite(values))):
        raise ValueError("V must be a finite real scalar or NxN pair matrix")
    if not np.allclose(values, values.T, atol=1e-12, rtol=0) or np.any(np.diag(values) != 0):
        raise ValueError("pair matrix must be symmetric with zero diagonal")
    return values.astype(float)


def full_hamiltonian(N, omega, delta, phase, V):
    """Independent tensor-product construction on all 2**N basis states.

    Sum_i [delta/2 (I-2n_i) + omega/2 (e^(-i phi)|0><r| + h.c.)]
    + sum_(i<j) V_ij n_i n_j. Sparse storage never projects onto symmetry sectors.
    """
    N = atom_count(N)
    if memory_estimate(N, 2, "full") > 256 * 2**20:
        raise ValueError("Full Hamiltonian exceeds the conservative 256 MiB allocation limit")
    omega, delta, phase = (_real_scalar(v, label) for v, label in
                            ((omega, "omega"), (delta, "delta"), (phase, "phase")))
    pairs = pair_matrix(N, V)
    coupling = omega / 2 * np.exp(-1j * phase)
    single = np.array([[delta / 2, coupling], [coupling.conjugate(), -delta / 2]])
    matrix = sparse.csr_matrix((2**N, 2**N), dtype=complex)
    occupations = []
    for site in range(N):
        matrix = matrix + local_operator(single, site, N)
        occupations.append(local_operator([[0, 0], [0, 1]], site, N))
    for i in range(N):
        for j in range(i + 1, N):
            if pairs[i, j] != 0:
                matrix = matrix + pairs[i, j] * occupations[i].multiply(occupations[j])
    matrix.eliminate_zeros()
    return matrix.tocsr()


def dicke_hamiltonian(N, omega, delta, phase, V):
    """Independent (N+1)x(N+1) tridiagonal construction, uniform pairs only."""
    N = atom_count(N)
    omega, delta, phase, V = (_real_scalar(v, label) for v, label in
                             ((omega, "omega"), (delta, "delta"), (phase, "phase"), (V, "V")))
    k = np.arange(N + 1, dtype=float)
    diagonal = delta * (N / 2 - k) + V * k * (k - 1) / 2
    raising = omega / 2 * np.sqrt((N - k[:-1]) * (k[:-1] + 1)) * np.exp(1j * phase)
    return sparse.diags([raising, diagonal, raising.conj()], [-1, 0, 1],
                        shape=(N + 1, N + 1), dtype=complex, format="csr")


def dicke_isometry(N):
    """Columns are normalised Dicke states, with exactly 2**N nonzeros."""
    N = atom_count(N)
    if 64 * 2**N > 256 * 2**20:
        raise ValueError("Dicke expansion exceeds the conservative 256 MiB allocation limit")
    rows = np.arange(2**N)
    sectors = np.fromiter((int(i).bit_count() for i in rows), dtype=int, count=len(rows))
    factors = np.array([1 / np.sqrt(float(comb(N, k))) for k in range(N + 1)])
    return sparse.csr_matrix((factors[sectors], (rows, sectors)), shape=(2**N, N + 1))


def evolve_uniform(Hbar, total_time, samples, initial=None):
    """Exponential action on a uniform grid from 0 to T; retain complex states.

    Both benchmark methods use this same sparse algorithm and output grid.
    Creating a ground state, validation and propagation are included in evolution.
    """
    total_time = _real_scalar(total_time, "total_time")
    if total_time < 0:
        raise ValueError("total_time must be nonnegative")
    if isinstance(samples, (bool, np.bool_)) or not isinstance(samples, Integral) or samples < 2:
        raise ValueError("samples must be an integer >= 2")
    matrix = sparse.csr_matrix(Hbar, dtype=complex)
    if matrix.shape[0] != matrix.shape[1] or not np.all(np.isfinite(matrix.data)):
        raise ValueError("Hbar must be finite and square")
    difference = matrix - matrix.conj().T
    if difference.nnz and np.max(abs(difference.data)) > 1e-12:
        raise ValueError("Hbar must be Hermitian")
    state = np.zeros(matrix.shape[0], dtype=complex) if initial is None else np.asarray(initial, dtype=complex)
    if initial is None:
        state[0] = 1
    if (state.shape != (matrix.shape[0],) or not np.all(np.isfinite(state))
            or not np.isclose(np.vdot(state, state), 1, atol=1e-12, rtol=0)):
        raise ValueError("initial state must have the correct shape and unit norm")
    if total_time == 0:
        return np.tile(state, (samples, 1))
    generator = -1j * matrix
    return expm_multiply(generator, state, start=0, stop=total_time, num=samples,
                         traceA=generator.diagonal().sum())


def compare_states(full, reduced, N):
    """Expand Dicke states and compare phases, fidelity, populations and leakage.

    Raw fidelity is retained (not clipped or renormalised). The scalar reported
    infidelity is floored at zero only for roundoff above unit fidelity.
    """
    basis = dicke_isometry(N)
    expanded = (basis @ reduced.T).T
    overlap = np.einsum("ti,ti->t", full.conj(), expanded)
    aligned = expanded * np.exp(-1j * np.angle(overlap))[:, None]
    fidelity = abs(overlap)**2
    projection = (basis.conj().T @ full.T).T
    sectors = np.array([i.bit_count() for i in range(2**N)])
    full_number = abs(full)**2 @ sectors
    reduced_number = abs(reduced)**2 @ np.arange(N + 1)
    full_norm = np.sum(abs(full)**2, axis=1)
    return {
        "fidelity": fidelity,
        "infidelity": max(0., float(1 - fidelity.min())),
        "raw_infidelity": float(1 - fidelity.min()),
        "max_absolute_fidelity_error": float(np.max(abs(1 - fidelity))),
        "max_amplitude_error": float(np.max(abs(full - expanded))),
        "max_aligned_amplitude_error": float(np.max(abs(full - aligned))),
        "max_excitation_error": float(np.max(abs(full_number - reduced_number))),
        "max_ground_probability_error": float(np.max(abs(abs(full[:, 0])**2 - abs(reduced[:, 0])**2))),
        "max_dicke_population_error": float(np.max(abs(abs(projection)**2 - abs(reduced)**2))),
        "max_norm_error": float(max(np.max(abs(full_norm - 1)), np.max(abs(np.sum(abs(reduced)**2, axis=1) - 1)))),
        "max_symmetric_leakage": max(0., float(np.max(full_norm - np.sum(abs(projection)**2, axis=1)))),
    }


def spatial_interactions(N, c6, spacing):
    """1D equally spaced chain, V_ij=C6/(a |i-j|)^6, angular units.

    If a is in micrometres, C6 must be in rad/s * micrometre^6.
    """
    N = atom_count(N)
    c6, spacing = _real_scalar(c6, "c6"), _real_scalar(spacing, "spacing")
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    distance = abs(np.arange(N)[:, None] - np.arange(N)[None, :]) * spacing
    result = np.zeros((N, N))
    nonzero = distance > 0
    result[nonzero] = c6 / distance[nonzero]**6
    if not np.all(np.isfinite(result)):
        raise ValueError("interactions exceed numeric range")
    return result


def csr_bytes(matrix):
    matrix = sparse.csr_matrix(matrix)
    return matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes


def memory_estimate(N, samples, method):
    """Conservative preflight bytes, not measured peak RSS or an exact bound.

    Includes local operators, multiple Hamiltonian/work buffers and trajectory;
    dense matrices are never allocated in the full calculation.
    """
    N = atom_count(N)
    if method not in ("full", "reduced"):
        raise ValueError("method must be full or reduced")
    dimension = 2**N if method == "full" else N + 1
    nnz = (N + 1) * dimension if method == "full" else 3 * dimension
    matrix = 24 * nnz + 8 * (dimension + 1)
    return int(8 * matrix + 3 * 16 * dimension * samples +
               (N * 24 * dimension if method == "full" else 0))
