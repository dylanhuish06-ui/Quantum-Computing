"""Exact bulk Dicke × physical defect models, with explicit reference guards.

All matrices are H/hbar in rad/s, or H/(hbar Omega_ref) when the caller uses
tau=Omega_ref*t. Raising has exp(+i phi); detuning is delta*(1/2-n).
The bulk consists of the first M sites, defects of the last m sites.
"""

from dataclasses import dataclass
from math import comb
from numbers import Integral

import numpy as np
from scipy import sparse
from numpy.typing import ArrayLike


def count(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def real_array(value: ArrayLike, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.shape != shape or array.dtype.kind not in "iuf" or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite and real, with shape {shape}")
    array = array.astype(float, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class FullParameters:
    omega: ArrayLike
    delta: ArrayLike
    phase: ArrayLike
    pairs: ArrayLike

    def __post_init__(self):
        n = len(np.atleast_1d(self.omega))
        for key in ("omega", "delta", "phase"):
            object.__setattr__(self, key, real_array(getattr(self, key), (n,), key))
        pairs = real_array(self.pairs, (n, n), "pairs")
        if not np.array_equal(pairs, pairs.T) or np.any(np.diag(pairs)):
            raise ValueError("pairs must be symmetric with zero diagonal")
        object.__setattr__(self, "pairs", pairs)

    @property
    def N(self):
        return len(self.omega)


@dataclass(frozen=True)
class HybridParameters:
    N: int
    m: int
    omega_bulk: float = 1.
    delta_bulk: float = 0.
    phase_bulk: float = 0.
    V_bulk: float = 0.
    omega_def: ArrayLike = ()
    delta_def: ArrayLike = ()
    phase_def: ArrayLike = ()
    V_cross: ArrayLike = ()
    V_def: ArrayLike = ()

    def __post_init__(self):
        count(self.N, "N", 1)
        count(self.m, "m")
        if self.m > self.N:
            raise ValueError("m cannot exceed N")
        for key in ("omega_bulk", "delta_bulk", "phase_bulk", "V_bulk"):
            object.__setattr__(self, key, float(real_array(getattr(self, key), (), key)))
        for key in ("omega_def", "delta_def", "phase_def", "V_cross"):
            object.__setattr__(self, key, real_array(getattr(self, key), (self.m,), key))
        pairs = np.asarray(self.V_def)
        if self.m == 0 and pairs.size == 0:
            pairs = np.zeros((0, 0))
        pairs = real_array(pairs, (self.m, self.m), "V_def")
        if not np.array_equal(pairs, pairs.T) or np.any(np.diag(pairs)):
            raise ValueError("V_def must be symmetric with zero diagonal")
        object.__setattr__(self, "V_def", pairs)

    @property
    def M(self):
        return self.N - self.m

    @property
    def dimension(self):
        return (self.M + 1) * 2**self.m

    def full_parameters(self):
        """Small-system validation only; never used by the large-N solver."""
        reference_guard(self.N)
        M, N = self.M, self.N
        pairs = np.zeros((N, N))
        pairs[:M, :M] = self.V_bulk * (np.ones((M, M)) - np.eye(M))
        pairs[:M, M:] = self.V_cross
        pairs[M:, :M] = np.broadcast_to(self.V_cross[:, None], (self.m, M))
        pairs[M:, M:] = self.V_def
        return FullParameters(
            np.r_[np.full(M, self.omega_bulk), self.omega_def],
            np.r_[np.full(M, self.delta_bulk), self.delta_def],
            np.r_[np.full(M, self.phase_bulk), self.phase_def], pairs)


def reference_guard(N):
    count(N, "N")
    if N > 12:
        raise ValueError("full-space reference/expansion is limited to N<=12")


def validate_bulk_symmetry(parameters, m, *, atol=0.):
    """Validate one instant and return its exact hybrid parameters.

    Default equality is exact. Nonzero atol is an explicit approximation;
    callers must quantify its model error. Sampled checks cannot prove
    symmetry between samples for an arbitrary time callback.
    """
    p = parameters
    count(m, "m")
    if p.N < 1 or m > p.N or not np.isfinite(atol) or atol < 0:
        raise ValueError("invalid atom count or tolerance")
    M = p.N - m
    def common(values, name):
        if not len(values):
            return 0.
        if np.max(abs(values - values[0])) > atol:
            raise ValueError(f"bulk exchange symmetry violated: {name}")
        return values[0]
    delta = common(p.delta[:M], "detuning")
    # Equality of complex drive coefficients also admits equivalent phase gauges.
    drive = common(p.omega[:M] * np.exp(1j*p.phase[:M]), "drive")
    pairs = p.pairs[:M, :M][np.triu_indices(M, 1)]
    V = common(pairs, "bulk pairs")
    cross = [common(p.pairs[:M, M+a], f"bulk coupling to defect {a}") for a in range(m)]
    return HybridParameters(p.N, m, abs(drive), delta, np.angle(drive), V,
                            p.omega[M:], p.delta[M:], p.phase[M:], cross, p.pairs[M:, M:])


def full_hamiltonian_reference(p: FullParameters) -> sparse.csr_matrix:
    """Independent sparse bit-string construction, retaining all 2**N states."""
    N = p.N
    reference_guard(N)
    d = 2**N
    states = np.arange(d, dtype=np.int64)
    bits = ((states[:, None] >> np.arange(N-1, -1, -1)) & 1).astype(float)
    diagonal = (0.5-bits) @ p.delta
    diagonal += 0.5*np.einsum("bi,ij,bj->b", bits, p.pairs, bits)
    rows, cols, data = [states], [states], [diagonal.astype(complex)]
    for site in range(N):
        low = states[bits[:, site] == 0]
        high = low | (1 << (N-1-site))
        c = 0.5*p.omega[site]*np.exp(1j*p.phase[site])
        rows.extend((high, low))
        cols.extend((low, high))
        data.extend((np.full(len(low), c), np.full(len(low), c.conjugate())))
    result = sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                               shape=(d, d)).tocsr()
    result.eliminate_zeros()
    return result


def hybrid_hamiltonian(p: HybridParameters) -> sparse.csr_matrix:
    """Direct sparse construction. No full N-qubit matrix or isometry is made."""
    M, m, d = p.M, p.m, p.dimension
    # Conservative workspace heuristic, not a measurement of peak RSS.
    if 8*24*d*(m+3) > 256*2**20 or m > 12:
        raise ValueError("hybrid construction exceeds the 256 MiB workspace budget")
    k = np.arange(M+1, dtype=float)
    diag = p.delta_bulk*(M/2-k) + p.V_bulk*k*(k-1)/2
    bulk = sparse.diags(diag, format="csr", dtype=complex)
    if M:
        raise_c = p.omega_bulk/2*np.exp(1j*p.phase_bulk)*np.sqrt((M-k[:-1])*(k[:-1]+1))
        bulk += sparse.diags([raise_c, raise_c.conj()], [-1, 1], shape=(M+1, M+1))
    defect = full_hamiltonian_reference(FullParameters(p.omega_def, p.delta_def, p.phase_def, p.V_def))
    result = sparse.kron(bulk, sparse.eye(2**m), format="csr")
    result += sparse.kron(sparse.eye(M+1), defect, format="csr")
    defect_bits = ((np.arange(2**m)[:, None] >> np.arange(m-1, -1, -1)) & 1)
    result += sparse.diags(np.outer(k, defect_bits @ p.V_cross).ravel(), format="csr")
    result.eliminate_zeros()
    return result.tocsr()


def dicke_isometry_bulk(N: int, m: int) -> sparse.csr_matrix:
    """B= D_(N-m) tensor I_(2**m), guarded small-system expansion."""
    reference_guard(N)
    count(m, "m")
    if m > N:
        raise ValueError("m cannot exceed N")
    M = N-m
    rows = np.arange(2**M)
    sectors = np.array([int(x).bit_count() for x in rows])
    weights = np.array([1/np.sqrt(comb(M, k)) for k in range(M+1)])
    dicke = sparse.csr_matrix((weights[sectors], (rows, sectors)), shape=(2**M, M+1))
    return sparse.kron(dicke, sparse.eye(2**m), format="csr")


def split_bulk_qubit(N, m, state):
    """Exactly promote the last bulk atom to the first distinguished qubit.

    This changes the representation without changing the physical state.
    Parameter reassignment and the decision to split belong to the caller.
    It cannot recover leakage already lost under an approximate grouped model.
    """
    count(N, "N", 1)
    count(m, "m")
    if m >= N:
        raise ValueError("no bulk qubit remains to split")
    M, d = N-m, 2**m
    state = np.asarray(state, complex)
    if state.shape != ((M+1)*d,) or not np.all(np.isfinite(state)):
        raise ValueError("invalid hybrid state")
    if 16*M*2*d > 256*2**20:
        raise ValueError("split state exceeds allocation limit")
    old = state.reshape(M+1, d)
    new = np.zeros((M, 2, d), complex)
    k = np.arange(M)
    new[:, 0, :] = np.sqrt((M-k)/M)[:, None]*old[:-1]
    new[:, 1, :] = np.sqrt((k+1)/M)[:, None]*old[1:]
    return new.ravel()
