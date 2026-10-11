"""Symmetry reduction, tensor construction and numerical equivalence."""

import os
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from quantum_project import KET_00, evolve_static, rydberg_hamiltonian
from quantum_project.symmetry import (
    compare_states, dicke_hamiltonian, dicke_isometry, evolve_uniform,
    exchange_hamiltonian, exchange_transformation, full_hamiltonian,
    memory_estimate, reduced_two_atom_hamiltonian, spatial_interactions,
)

OMEGA = 2 * np.pi * 1e6


def test_explicit_exchange_transformation():
    transform = exchange_transformation()
    expected = np.array([[1, 0, 0, 0], [0, 1/np.sqrt(2), 1/np.sqrt(2), 0],
                         [0, 1/np.sqrt(2), -1/np.sqrt(2), 0], [0, 0, 0, 1]])
    np.testing.assert_array_equal(transform, expected)
    np.testing.assert_allclose(transform.conj().T @ transform, np.eye(4), atol=3e-16)


@pytest.mark.parametrize("parameters", [(1, 0, 0, 0), (1, -.7, .37, 3), (-2, .3, -1, -5)])
def test_n2_dark_state_exact_hamiltonian_and_dynamics(parameters):
    omega, delta, phase, V = parameters
    args = (OMEGA*omega, OMEGA*delta, phase, OMEGA*V)
    original = rydberg_hamiltonian(*args)
    full = full_hamiltonian(2, *args).toarray()
    np.testing.assert_allclose(full, original, atol=1e-9, rtol=0)
    transformed = exchange_hamiltonian(*args)
    np.testing.assert_allclose(transformed[2], 0, atol=1e-9, rtol=0)
    np.testing.assert_allclose(transformed[:, 2], 0, atol=1e-9, rtol=0)
    reduced = reduced_two_atom_hamiltonian(*args)
    np.testing.assert_allclose(reduced, dicke_hamiltonian(2, *args).toarray(), atol=4e-9, rtol=0)
    times = np.linspace(0, .5e-6, 51)
    pure = evolve_static(KET_00, times, original)
    dark = pure @ exchange_transformation()[:, 2].conj()
    assert np.max(abs(dark)**2) < 1e-28
    reduced_states = evolve_uniform(reduced, times[-1], len(times))
    assert compare_states(pure, reduced_states, 2)["max_amplitude_error"] < 2e-13


@pytest.mark.parametrize("N", range(1, 8))
def test_dimensions_normalisation_pair_count_and_projection(N):
    basis = dicke_isometry(N)
    assert basis.shape == (2**N, N+1)
    np.testing.assert_allclose((basis.T @ basis).toarray(), np.eye(N+1), atol=4e-15)
    full = full_hamiltonian(N, 1.3, -.4, .7, 2.1)
    reduced = dicke_hamiltonian(N, 1.3, -.4, .7, 2.1)
    assert full.shape == (2**N, 2**N)
    assert reduced.shape == (N+1, N+1)
    np.testing.assert_allclose((basis.T @ full @ basis).toarray(), reduced.toarray(), atol=1e-12)
    # Stronger than projection alone: verifies no coupling out of symmetric space.
    np.testing.assert_allclose((full @ basis - basis @ reduced).toarray(), 0, atol=1e-12)
    diagonal = full_hamiltonian(N, 0, 0, 0, 2.1).diagonal()
    expected = [2.1 * i.bit_count() * (i.bit_count()-1) / 2 for i in range(2**N)]
    np.testing.assert_allclose(diagonal, expected, atol=1e-13)
    k = np.arange(N+1)
    np.testing.assert_allclose(dicke_hamiltonian(N, 0, .7, 0, 0).diagonal(), .7*(N/2-k))


@pytest.mark.parametrize("N", [1, 2, 3, 4, 5, 6])
@pytest.mark.parametrize("interaction", [0, .2, -2])
def test_full_reduced_agreement(N, interaction):
    args = (OMEGA, .2*OMEGA, .37, interaction*OMEGA)
    full = evolve_uniform(full_hamiltonian(N, *args), .25e-6, 31)
    reduced = evolve_uniform(dicke_hamiltonian(N, *args), .25e-6, 31)
    checks = compare_states(full, reduced, N)
    assert checks["infidelity"] < 2e-12
    assert checks["max_amplitude_error"] < 2e-12
    assert checks["max_excitation_error"] < 2e-12
    assert checks["max_dicke_population_error"] < 2e-12
    assert checks["max_norm_error"] < 2e-12


def test_collective_coupling_and_excitation_support():
    for n in range(2, 6):
        matrix = dicke_hamiltonian(n, 1.7, .3, .41, -.9).toarray()
        basis = dicke_isometry(n).toarray()
        counts = np.array([i.bit_count() for i in range(2**n)])
        for k in range(n + 1):
            np.testing.assert_array_equal(counts[basis[:, k] != 0], k)
            assert np.sum(abs(basis[:, k])**2 * counts) == pytest.approx(k)
            assert matrix[k, k] == pytest.approx(.3*(n/2-k) - .9*k*(k-1)/2)
            if k < n:
                assert matrix[k+1, k] == pytest.approx(1.7/2*np.exp(.41j)*np.sqrt((n-k)*(k+1)))


def test_allocation_guards_before_exponential_arrays():
    with pytest.raises(ValueError, match="allocation limit"):
        full_hamiltonian(30, 1, 0, 0, 1)
    with pytest.raises(ValueError, match="allocation limit"):
        dicke_isometry(30)


def test_small_launcher_demo_uses_own_worktree_without_zeno(tmp_path):
    root = Path(__file__).resolve().parents[1]
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    code = """
import importlib.abc, runpy, sys
class NoZeno(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if 'zeno' in fullname.lower():
            raise AssertionError('Zeno dependency: ' + fullname)
sys.meta_path.insert(0, NoZeno())
runpy.run_path(sys.argv[1], run_name='__main__')
"""
    result = subprocess.run([sys.executable, "-c", code, str(root / "scripts/symmetry_demo.py")],
                            cwd=tmp_path, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "PASS" in result.stdout and "N=5:" in result.stdout
    assert not list(tmp_path.iterdir())


def test_reduced_sizes_continue_when_no_full_calculation_fits(tmp_path):
    import json
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([sys.executable, str(root / "scripts/symmetry_benchmark.py"),
                             "--full-max", "4", "--reduced-extra", "--memory-budget-mib", ".009", "--dry-run"],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    config = json.loads(result.stdout)
    assert memory_estimate(2, 41, "reduced") <= config["budget_bytes"]
    assert config["budget_bytes"] < memory_estimate(2, 41, "full")
    assert config["full_sizes"] == []
    assert config["reduced_sizes"] == [2]
    assert all(r["reason"] for r in config["skips"])
    assert not list(tmp_path.iterdir())


def test_final_timings_rejected_in_codespaces(tmp_path):
    root = Path(__file__).resolve().parents[1]
    env = {**os.environ, "CODESPACES": "true", "CI": "false", "GITHUB_ACTIONS": "false"}
    result = subprocess.run([sys.executable, str(root / "scripts/symmetry_benchmark.py"),
                             "--session-kind", "final", "--machine-label", "test", "--output-dir", str(tmp_path)],
                            env=env, capture_output=True, text=True)
    assert result.returncode != 0 and "physical machine" in result.stderr
    assert not list(tmp_path.iterdir())


def test_independent_dense_exponential_and_nonsymmetric_state():
    matrix = full_hamiltonian(3, OMEGA, .2*OMEGA, .3, .4*OMEGA)
    initial = np.zeros(8, dtype=complex)
    initial[1] = 1  # A full-space state outside the symmetric sector.
    result = evolve_uniform(matrix, .2e-6, 11, initial)
    np.testing.assert_allclose(result[-1], expm(-1j*matrix.toarray()*.2e-6) @ initial, atol=2e-14)
    assert np.linalg.norm(dicke_isometry(3).T @ initial)**2 == pytest.approx(1/3)


def test_spatial_symmetry_breaking_and_mean_pair_projection():
    pairs = spatial_interactions(4, 2*OMEGA, 1)
    assert pairs[0, 1] == pytest.approx(2*OMEGA)
    assert pairs[0, 2] == pytest.approx(2*OMEGA/64)
    assert pairs[0, 3] == pytest.approx(2*OMEGA/729)
    mean = pairs[np.triu_indices(4, 1)].mean()
    basis = dicke_isometry(4)
    full_h = full_hamiltonian(4, OMEGA, 0, 0, pairs)
    reduced_h = dicke_hamiltonian(4, OMEGA, 0, 0, mean)
    np.testing.assert_allclose((basis.T @ full_h @ basis).toarray(), reduced_h.toarray(), atol=1e-8)
    assert np.linalg.norm((full_h @ basis - basis @ reduced_h).toarray()) > OMEGA
    checks = compare_states(evolve_uniform(full_h, .5e-6, 51), evolve_uniform(reduced_h, .5e-6, 51), 4)
    assert checks["infidelity"] > .01
    assert checks["max_symmetric_leakage"] > .01
    # Two atoms have only one pair: a spatial interaction cannot break exchange.
    n2 = spatial_interactions(2, OMEGA, 2)
    np.testing.assert_allclose(full_hamiltonian(2, 1, 2, .3, n2).toarray(),
                               full_hamiltonian(2, 1, 2, .3, OMEGA/64).toarray())


@pytest.mark.parametrize("N", [0, -1, True, 2.5])
def test_invalid_counts(N):
    with pytest.raises(ValueError):
        dicke_hamiltonian(N, 1, 0, 0, 0)


@pytest.mark.parametrize("V", [np.eye(3), np.ones((2, 2)), [[0, 1], [2, 0]], 1j])
def test_invalid_pairs(V):
    with pytest.raises(ValueError):
        full_hamiltonian(2, 1, 0, 0, V)


def test_memory_estimate_grows_and_benchmark_ci_guard(tmp_path):
    assert memory_estimate(12, 41, "full") > memory_estimate(12, 41, "reduced") * 100
    root = Path(__file__).resolve().parents[1]
    command = [sys.executable, str(root / "scripts" / "symmetry_benchmark.py"), "--full-max", "2", "--output-dir", str(tmp_path)]
    env = {**os.environ, "CI": "true", "PYTHONPATH": str(root / "src")}
    forbidden = subprocess.run(command, env=env, capture_output=True, text=True)
    assert forbidden.returncode != 0 and "prohibited" in forbidden.stderr
    allowed = subprocess.run(command + ["--dry-run"], env=env, capture_output=True, text=True)
    assert allowed.returncode == 0
    assert not list(tmp_path.iterdir())


def test_benchmark_summary_medians_iqr_and_memory_are_from_raw_rows():
    path = Path(__file__).resolve().parents[1] / "scripts/symmetry_benchmark.py"
    spec = importlib.util.spec_from_file_location("symmetry_benchmark_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = [{"N": 2, "method": "full", "construction_s": i, "evolution_s": 2*i,
             "total_s": 3*i, "hamiltonian_bytes": 100, "trajectory_bytes": 200,
             "estimated_working_bytes": 1000} for i in range(1, 11)]
    summary = module.summarize(rows)[0]
    assert summary["repeats"] == 10
    assert summary["construction_median_s"] == 5.5
    assert summary["construction_iqr_s"] == 4.5
    assert summary["total_median_s"] == 16.5
    assert summary["hamiltonian_bytes"] == 100


def test_benchmark_memory_guard_skips_large_full_sizes(tmp_path):
    root = Path(__file__).resolve().parents[1]
    command = [sys.executable, str(root / "scripts/symmetry_benchmark.py"), "--full-max", "20",
               "--memory-budget-mib", "8", "--dry-run", "--output-dir", str(tmp_path)]
    result = subprocess.run(command, env={**os.environ, "CI": "true", "PYTHONPATH": str(root / "src")},
                            capture_output=True, text=True)
    assert result.returncode == 0
    import json
    config = json.loads(result.stdout)
    assert config["skipped_full_sizes"]
    assert max(config["full_sizes"]) < 20
    assert 64 in config["reduced_sizes"]
    assert not list(tmp_path.iterdir())


@pytest.fixture
def recorded_session(tmp_path):
    import shutil
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("symmetry_figures_test", root / "scripts/symmetry_figures.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for path in (root / "results/data/symmetry_provisional").glob("*"):
        shutil.copyfile(path, tmp_path / path.name)
    for name in ("symmetry_validation.json", "symmetry_agreement.csv", "symmetry_broken.csv"):
        shutil.copyfile(root / "results/data" / name, tmp_path / name)
    return module, tmp_path


def test_recorded_figures_validate_existing_session(recorded_session):
    module, directory = recorded_session
    metadata, summary, agreement, broken = module.load_recorded_data(directory, directory)
    assert metadata["complete"]
    assert len(summary) == len(metadata["full_sizes"]) + len(metadata["reduced_sizes"])
    assert len(agreement) == len(metadata["full_sizes"])
    assert broken


@pytest.mark.parametrize("corruption", ["session", "summary", "repeat", "agreement"])
def test_recorded_figures_reject_mixed_or_stale_data(recorded_session, corruption):
    import csv
    import json
    module, directory = recorded_session
    if corruption == "session":
        path = directory / "symmetry_validation.json"
        content = json.loads(path.read_text())
        content["benchmark_session_id"] = "another-machine-session"
        path.write_text(json.dumps(content))
    else:
        filename, field, value = {
            "summary": ("symmetry_benchmark_summary.csv", "total_median_s", "999"),
            "repeat": ("symmetry_benchmark_raw.csv", "repeat", "999"),
            "agreement": ("symmetry_agreement.csv", "infidelity", "0.5"),
        }[corruption]
        path = directory / filename
        rows = module.read_rows(path)
        rows[0][field] = value
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
    with pytest.raises(ValueError):
        module.load_recorded_data(directory, directory)
