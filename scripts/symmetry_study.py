"""Untimed physical validation and broken-symmetry data; never collects timings."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from quantum_project import KET_00, evolve_static, rydberg_hamiltonian
from quantum_project.symmetry import (compare_states, dicke_hamiltonian, dicke_isometry,
    evolve_uniform, exchange_transformation, full_hamiltonian, spatial_interactions,
    memory_estimate, reduced_two_atom_hamiltonian)
from symmetry_benchmark import available_memory

ROOT = Path(__file__).resolve().parents[1]
OMEGA = 2 * np.pi * 1e6


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-dir", type=Path, default=ROOT / "results/data/symmetry_provisional")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/data")
    args = parser.parse_args()
    metadata = json.loads((args.benchmark_dir / "symmetry_benchmark_metadata.json").read_text())
    if not metadata.get("complete"):
        raise ValueError("Benchmark session is incomplete")
    config = metadata["parameters"]
    omega = config["omega_rad_s"]
    parameters = (omega, config["delta_over_omega"]*omega, config["phase_rad"], config["V_over_omega"]*omega)
    arrays = {"times_s": np.linspace(0, config["total_time_s"], config["samples"])}
    agreement = []
    for N in metadata["full_sizes"]:
        # Include retained trajectories plus expansion/comparison temporaries.
        retained = sum(a.nbytes for a in arrays.values())
        estimate = retained + memory_estimate(N, config["samples"], "full") + 8 * 16 * 2**N * config["samples"]
        available = available_memory()
        budget = min(256 * 2**20, available // 4 if available is not None else 64 * 2**20)
        if estimate > budget:
            raise MemoryError(f"Validation N={N} needs approximately {estimate} bytes, above budget {budget}; use a smaller session")
        full = evolve_uniform(full_hamiltonian(N, *parameters), config["total_time_s"], config["samples"])
        reduced = evolve_uniform(dicke_hamiltonian(N, *parameters), config["total_time_s"], config["samples"])
        checks = compare_states(full, reduced, N)
        arrays[f"N{N}_fidelity"] = checks.pop("fidelity")
        arrays[f"N{N}_full_states"], arrays[f"N{N}_reduced_states"] = full, reduced
        agreement.append({"N": N, "full_dimension": 2**N, "reduced_dimension": N+1, **checks})
        if any(checks[key] > 1e-10 for key in ("max_absolute_fidelity_error", "max_excitation_error",
                "max_dicke_population_error", "max_ground_probability_error", "max_norm_error")):
            raise AssertionError(f"Symmetry equivalence failed at N={N}: {checks}")
    # Independent milestone N=2 comparison, including antisymmetric population.
    pure = evolve_static(KET_00, arrays["times_s"], rydberg_hamiltonian(*parameters))
    symmetric = evolve_uniform(reduced_two_atom_hamiltonian(*parameters), config["total_time_s"], config["samples"])
    milestone = compare_states(pure, symmetric, 2)
    milestone.pop("fidelity")
    milestone["max_dark_population"] = float(np.max(abs(pure @ exchange_transformation()[:, 2].conj())**2))
    bell = pure @ exchange_transformation()[:, 1].conj()
    milestone["max_P_Psi_plus_error"] = float(np.max(abs(abs(bell)**2 - abs(symmetric[:, 1])**2)))
    milestone["max_P_rr_error"] = float(np.max(abs(abs(pure[:, 3])**2 - abs(symmetric[:, 2])**2)))
    if milestone["max_absolute_fidelity_error"] > 1e-12 or milestone["max_dark_population"] > 1e-25:
        raise AssertionError("Milestone comparison failed")

    broken = []
    broken_arrays = {"times_s": np.linspace(0, .5e-6, 101), "alpha": np.linspace(0, 1, 11)}
    for N in (2, 3, 4, 5, 6):
        # a=1 micrometre; C6=2*Omega micrometre^6; nearest-neighbour V=2*Omega.
        spatial = spatial_interactions(N, 2*OMEGA, 1)
        upper = np.triu_indices(N, 1)
        mean = spatial[upper].mean()
        uniform = mean * (np.ones((N, N)) - np.eye(N))
        reduced = evolve_uniform(dicke_hamiltonian(N, OMEGA, 0, 0, mean), .5e-6, 101)
        broken_arrays[f"N{N}_reduced_states"] = reduced
        for index, alpha in enumerate(broken_arrays["alpha"]):
            pairs = (1-alpha)*uniform + alpha*spatial
            full_h = full_hamiltonian(N, OMEGA, 0, 0, pairs)
            full = evolve_uniform(full_h, .5e-6, 101)
            checks = compare_states(full, reduced, N)
            broken_arrays[f"N{N}_alpha{index}_fidelity"] = checks.pop("fidelity")
            broken_arrays[f"N{N}_alpha{index}_full_states"] = full
            basis = dicke_isometry(N)
            coupling = full_h @ basis - basis @ dicke_hamiltonian(N, OMEGA, 0, 0, mean)
            broken.append({"N": N, "alpha": float(alpha), "mean_V_over_omega": mean/OMEGA,
                           "pair_std_over_abs_mean": float(np.std(pairs[upper])/abs(mean)),
                           "leakage_coupling_frobenius_over_omega": float(np.linalg.norm(coupling.data)/OMEGA), **checks})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.output_dir / "symmetry_agreement.csv", agreement)
    write_rows(args.output_dir / "symmetry_broken.csv", broken)
    np.savez_compressed(args.output_dir / "symmetry_agreement.npz", **arrays)
    np.savez_compressed(args.output_dir / "symmetry_broken.npz", **broken_arrays)
    report = {"exact_parameters": config, "benchmark_session_id": metadata["session_id"],
              "milestone_N2": milestone, "exact_agreement": agreement,
              "broken_symmetry": {"omega_rad_s": OMEGA, "delta_rad_s": 0, "phase_rad": 0,
                                  "spacing_micrometres": 1, "C6_rad_s_micrometre6": 2*OMEGA,
                                  "total_time_s": .5e-6, "samples": 101,
                                  "interpolation": "Vij(alpha)=(1-alpha)*mean(Vspatial)+alpha*Vspatial; mean fixed",
                                  "physical_chain": "alpha=1 only; intermediate alpha are a controlled mathematical interpolation",
                                  "approximation": "uniform pair mean, exact projected Hamiltonian but omits leakage",
                                  "endpoints": [row for row in broken if row["alpha"] == 1]}}
    (args.output_dir / "symmetry_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"milestone_N2": milestone,
                      "max_full_reduced_infidelity": max(row["infidelity"] for row in agreement),
                      "broken_endpoints": [{k:r[k] for k in ("N", "infidelity", "max_symmetric_leakage")} for r in broken if r["alpha"] == 1]}, indent=2))


if __name__ == "__main__":
    main()
