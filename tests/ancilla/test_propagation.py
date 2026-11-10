import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project.ancilla.propagation import evolve


X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1])


def H(t):
    return (1+.3*t)*X + (.4+np.sin(2*t))*Z


def test_time_ordering_and_independent_midpoint_convergence():
    initial = np.array([1, 1j])/np.sqrt(2)
    assert np.linalg.norm(H(0)@H(.7)-H(.7)@H(0)) > .1
    reference = evolve(H, initial, [0, 1.3], dt=.03)[-1]
    tighter = evolve(H, initial, [0, 1.3], dt=.02, rtol=3e-14, atol=1e-15)[-1]
    np.testing.assert_allclose(reference, tighter, atol=1e-11)
    assert np.linalg.norm(reference-expm(-1j*H(0)*1.3)@initial) > .1
    errors, local_errors = [], []
    for steps in (20, 40, 80):
        dt = 1.3/steps
        sparse = evolve(H, initial, [0, 1.3], method="midpoint", dt=dt)
        dense = evolve(H, initial, [0, 1.3], method="dense_midpoint", dt=dt)
        np.testing.assert_allclose(sparse, dense, atol=2e-14)
        assert abs(np.vdot(sparse[-1], sparse[-1])-1) < 2e-13
        errors.append(np.linalg.norm(sparse[-1]-reference))
        exact_local = evolve(H, initial, [0, dt], dt=dt/4)[-1]
        local = evolve(H, initial, [0, dt], method="midpoint", dt=dt)[-1]
        local_errors.append(np.linalg.norm(local-exact_local))
    assert 3.8 < errors[0]/errors[1] < 4.2
    assert 3.8 < errors[1]/errors[2] < 4.2
    assert 7 < local_errors[1]/local_errors[2] < 9


def test_dimensionless_conversion_and_no_renormalisation():
    initial = np.array([1, 0], complex)
    scale = 2*np.pi*1e6
    dimensionless = evolve(H, initial, [0, .2, .7], dt=.04)
    physical = evolve(lambda t: scale*H(scale*t), initial,
                      np.array([0, .2, .7])/scale, dt=.04/scale)
    np.testing.assert_allclose(dimensionless, physical, atol=1e-12)
    unnormalised = evolve(H, .2*initial, [0, .7], require_normalised=False)
    assert abs(np.linalg.norm(unnormalised[-1])-.2) < 1e-12


@pytest.mark.parametrize("kwargs", [{"dt":0}, {"method":"frozen"}, {"rtol":np.nan}])
def test_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        evolve(H, [1, 0], [0, 1], **kwargs)


def test_invalid_states_times_and_hamiltonians():
    for times in ([1, 0], [], [0, np.inf]):
        with pytest.raises(ValueError):
            evolve(H, [1, 0], times)
    with pytest.raises(ValueError):
        evolve(H, [2, 0], [0, 1])
    with pytest.raises(ValueError):
        evolve(lambda t: [[0, 1], [0, 0]], [1, 0], [0, 1])
