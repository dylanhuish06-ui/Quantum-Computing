"""Ideal, nonlocal total-spin syndrome. No compiled physical oracle is implied."""

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from .models import dicke_isometry_bulk, count
from .propagation import checked_matrix, evolve, time_grid


def symmetric_projector(N):
    if N > 8:
        raise ValueError("dense syndrome reference limited to N<=8")
    B = dicke_isometry_bulk(N, 0).toarray()
    return B@B.conj().T


def checked_projector(P):
    P = np.asarray(P, dtype=complex)
    if (P.ndim != 2 or P.shape[0] != P.shape[1] or not np.all(np.isfinite(P))
            or not np.allclose(P, P.conj().T, atol=1e-12, rtol=0)
            or not np.allclose(P@P, P, atol=1e-12, rtol=0)):
        raise ValueError("P must be an orthogonal projector")
    return P


def syndrome_unitary(P):
    P = checked_projector(P)
    Q = np.eye(len(P))-P
    return np.kron(P, np.eye(2))+np.kron(Q, [[0, 1], [1, 0]])


def measure_density(rho, P, outcome=None):
    """Return the unnormalised outcome block, or the trace-preserving channel."""
    P = checked_projector(P)
    rho = np.asarray(rho, complex)
    if rho.shape != P.shape or not np.all(np.isfinite(rho)):
        raise ValueError("density shape must match P")
    Q = np.eye(len(P))-P
    if outcome is None:
        return P@rho@P + Q@rho@Q
    if outcome not in (0, 1):
        raise ValueError("outcome must be 0, 1 or None")
    A = P if outcome == 0 else Q
    return A@rho@A


@dataclass
class SelectiveResult:
    times: np.ndarray
    branch: np.ndarray
    success: np.ndarray

    @property
    def conditional(self):
        result = np.full_like(self.branch, np.nan)
        valid = self.success > 0
        result[valid] = self.branch[valid]/np.sqrt(self.success[valid, None])
        return result

    def metrics(self, target, P):
        """Target must be a unit vector at each reported time; no clipping."""
        target = np.asarray(target, complex)
        if target.shape != self.branch.shape or not np.allclose(
                np.sum(abs(target)**2, axis=1), 1, atol=1e-9, rtol=0):
            raise ValueError("target trajectory must be normalised with matching shape")
        product = abs(np.einsum("ti,ti->t", target.conj(), self.branch))**2
        conditional = np.divide(product, self.success, out=np.full_like(product, np.nan),
                                where=self.success > 0)
        states = self.conditional
        loss = 1-np.einsum("ti,ij,tj->t", states.conj(), P, states).real
        return {"conditional_fidelity": conditional, "success_probability": self.success,
                "successful_target_probability": product, "conditional_loss": loss}


def selective_evolution(H, initial, checks, P, **solver):
    """Retain the unnormalised all-zero branch; check at checks[1:].

    Integration dt is supplied separately via solver. All checks have ideal
    instantaneous oracle/readout/reset. A zero-probability conditional state is
    undefined (NaN), while its successful-target probability is zero.
    """
    P = checked_projector(P)
    times = time_grid(checks)
    state = evolve(H, initial, times[:1])[0]
    if len(state) != len(P):
        raise ValueError("projector dimension must match the state")
    branch = [state.copy()]
    for left, right in zip(times[:-1], times[1:]):
        state = P@evolve(H, state, [left, right], require_normalised=False, **solver)[-1]
        branch.append(state.copy())
    branch = np.array(branch)
    return SelectiveResult(times, branch, np.sum(abs(branch)**2, axis=1))


def nonselective_evolution(H, initial_density, checks, P, *, dt=.01,
                          method="dop853", rtol=2e-12, atol=2e-14):
    """Physical trace-preserving density trajectory, including both outcomes."""
    times, P = time_grid(checks), checked_projector(P)
    rho = np.asarray(initial_density, complex).copy()
    if (rho.shape != P.shape or not np.all(np.isfinite(rho))
            or not np.allclose(rho, rho.conj().T, atol=1e-12, rtol=0)
            or abs(np.trace(rho)-1) > 1e-10 or np.linalg.eigvalsh(rho).min() < -1e-12):
        raise ValueError("initial_density must be a positive unit-trace density matrix")
    if (dt <= 0 or not np.isfinite(dt) or method not in ("dop853", "dense_midpoint")
            or not (0 < rtol < 1 and 0 < atol < 1)):
        raise ValueError("invalid density solver controls")
    d = len(P)
    if d > 64 or 3*16*len(times)*d*d > 256*2**20:
        raise ValueError("density reference exceeds allocation guard")
    result = [rho.copy()]
    for left, right in zip(times[:-1], times[1:]):
        if method == "dop853":
            def rhs(t, vector):
                matrix = checked_matrix(H(t), d).toarray()
                r = vector.reshape(d, d)
                return (-1j*(matrix@r-r@matrix)).ravel()
            sol = solve_ivp(rhs, (left, right), rho.ravel(), method="DOP853",
                            rtol=rtol, atol=atol, max_step=dt)
            if not sol.success:
                raise RuntimeError(sol.message)
            rho = sol.y[:, -1].reshape(d, d)
        else:
            steps = max(1, int(np.ceil((right-left)/dt-1e-12)))
            step = (right-left)/steps
            for j in range(steps):
                U = expm(-1j*step*checked_matrix(H(left+(j+.5)*step), d).toarray())
                rho = U@rho@U.conj().T
        rho = measure_density(rho, P)
        result.append(rho.copy())
    return np.array(result)


def spin_phase_register(N):
    """An injective *ideal* QPE spectrum; does not compile exp(i J² t).

    A=4J² has distinct integer eigenvalues l(l+2), l=N,N-2,... .
    Choose 2**b > N(N+2), so all phases A/2**b are exact b-bit fractions.
    """
    count(N, "N", 1)
    b = (N*(N+2)).bit_length()
    labels = np.arange(N, -1, -2, dtype=np.int64)
    eigenvalues = labels*(labels+2)
    return {"register_qubits": b, "iterative_reusable_qubits": 1,
            "controlled_evolutions_per_check": b,
            "parallel_register_readouts": b, "iterative_resets_per_check": b,
            "phase_denominator": 2**b, "eigenvalue_labels": eigenvalues.tolist(),
            "phases": (eigenvalues/2**b).tolist(),
            "compiled_two_qubit_gates": None, "compiled_depth": None}
