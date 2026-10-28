"""Sparse large-N sweep and guarded adaptive-disorder exit experiment.

Run PYTHONPATH=src python scripts/ancilla/validate_gate5.py --output-dir PATH.
Preserves existing reports. No full N-qubit expansion at large N. Wall times
are provisional; no speedup measurement is claimed. Results use dimensionless
tau=Omega_ref*t, h=H/(hbar*Omega_ref), Omega_ref=2*pi*1e6 rad/s.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import scipy
from scipy.linalg import expm
from quantum_project.ancilla.models import HybridParameters, FullParameters, full_hamiltonian_reference
from quantum_project.ancilla.propagation import evolve_hybrid
from quantum_project.ancilla.approximation import adaptive_static_grouping


def profile(N):
    """One physically distinguished atom; bulk permutation symmetry is exact."""
    return lambda tau: HybridParameters(
        N, 1, 1+.1*np.sin(tau), -.5+2*tau, .2*tau, .2/N,
        [1.1], [.35+tau], [-.3*tau], [.3/N], [[0.]])


def error(a, b):
    return float(np.max(np.linalg.norm(a-b, axis=1)))


def study(sizes=(64, 256, 1024)):
    rows, trajectories = [], []
    for N in sizes:
        times = np.linspace(0, .4, 21)
        initial = np.zeros(2*N, complex)
        initial[0] = 1
        kwargs = {"method": "controlled_midpoint", "max_seconds": 90.}
        coarse = evolve_hybrid(profile(N), initial, times, dt=.01, action_atol=1e-11, **kwargs)
        medium = evolve_hybrid(profile(N), initial, times, dt=.005, action_atol=1e-11, **kwargs)
        fine = evolve_hybrid(profile(N), initial, times, dt=.0025, action_atol=1e-11, **kwargs)
        tighter = evolve_hybrid(profile(N), initial, times, dt=.0025, action_atol=1e-13, **kwargs)
        scipy_action = evolve_hybrid(profile(N), initial, times, dt=.0025, method="midpoint")
        output_fine = evolve_hybrid(profile(N), initial, np.linspace(0, .4, 41),
                                   dt=.0025, action_atol=1e-11, **kwargs)[::2]
        e1, e2 = error(coarse, medium), error(medium, fine)
        population = abs(fine.reshape(len(times), N, 2))**2
        bulk = np.einsum("tka,k->t", population, np.arange(N))/(N-1)
        defect = population[:, :, 1].sum(axis=1)
        excitation = ((N-1)*bulk+defect)/N
        row = {"N": N, "m": 1, "dimension": 2*N,
               "coarse_medium_l2": e1, "medium_fine_l2": e2,
               "midpoint_order": float(np.log2(e1/e2)),
               "estimated_fine_time_error_l2": e2/3,
               "action_refinement_l2": error(fine, tighter),
               "independent_action_l2": error(fine, scipy_action),
               "output_grid_l2": error(fine, output_fine),
               "norm_drift": float(np.max(abs(population.sum(axis=(1, 2))-1))),
               "final_excitation_fraction": float(excitation[-1]),
               "final_defect_population": float(defect[-1])}
        assert 1.95 < row["midpoint_order"] < 2.05
        assert row["estimated_fine_time_error_l2"] < 1e-4
        for key in ("action_refinement_l2", "independent_action_l2", "output_grid_l2", "norm_drift"):
            assert row[key] < 1e-9, (N, key, row[key])
        assert np.all(np.isfinite(fine))
        assert min(bulk.min(), defect.min()) >= -1e-12
        assert max(bulk.max(), defect.max()) <= 1+1e-10
        rows.append(row)
        trajectories.extend({"N": N, "tau": float(t), "time_seconds": float(t/(2*np.pi*1e6)),
                             "bulk_excitation_fraction": float(b), "defect_population": float(d),
                             "total_excitation_fraction": float(e)}
                            for t, b, d, e in zip(times, bulk, defect, excitation))
    N, T = 5, .7
    H = full_hamiltonian_reference(FullParameters(np.ones(N), [-.8, .7, -.4, .3, 1.1],
                                                 np.zeros(N), np.zeros((N, N))))
    adaptive = []
    for max_m in range(N+1):
        result = adaptive_static_grouping(H, N, T, error_budget=1e-4, max_m=max_m)
        B = result["isometry"].toarray()
        initial = np.eye(B.shape[1])[0]  # Same physical all-ground state at every m.
        exact = expm(-1j*T*H.toarray())@(B@initial)
        approximate = B@expm(-1j*T*result["hamiltonian"].toarray())@initial
        measured = float(np.linalg.norm(exact-approximate))
        last = result["history"][-1]
        assert measured <= last["state_error_bound"]+1e-12
        adaptive.append(last | {"actual_state_error_l2": measured, "budget_met": result["budget_met"]})
    assert not adaptive[0]["budget_met"] and adaptive[-1]["budget_met"]
    return rows, trajectories, adaptive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir
    if output.exists() and any(output.iterdir()):
        parser.error("preserving existing artifacts")
    started = time.monotonic()
    rows, trajectories, adaptive = study()
    output.mkdir(parents=True, exist_ok=True)
    for name, data in (("residuals.csv", rows), ("observables.csv", trajectories), ("adaptive.csv", adaptive)):
        with (output/name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    sources = sorted((ROOT/"src/quantum_project/ancilla").glob("*.py"))+[Path(__file__).resolve()]
    report = {"gate": 5, "passed": True, "utc": datetime.now(timezone.utc).isoformat(),
              "upstream_base": "c6b1cac624ba995966df1249467c7e1367d83b05",
              "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
              "artifact_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()},
              "parameters": {"sizes": [64, 256, 1024], "m": 1, "T_tau": .4, "output_samples": [21, 41],
                             "dt_tau": [.01, .005, .0025], "action_atol": [1e-11, 1e-13],
                             "Omega_ref_rad_per_s": 2*np.pi*1e6, "profile": "profile(N) in this script"},
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "elapsed_seconds_provisional": time.monotonic()-started,
              "rows": rows, "adaptive": adaptive,
              "limitations": ["Taylor truncation bound excludes roundoff and time discretisation",
                              "midpoint error estimate is asymptotic Richardson, not a rigorous bound",
                              "generic-disorder adaptive reference is static and N<=8; worst case reaches full dimension",
                              "large-N trajectories are exact structured hybrid model, without ancilla measurement",
                              "raw complex states regenerated by recipe; compact observable trajectories stored"]}
    (output/"validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
