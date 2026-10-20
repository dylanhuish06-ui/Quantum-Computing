"""Small S1 validation: no timings, plots or saved files; N <= 5."""

from pathlib import Path
import sys

# Select this checkout even when a shared environment is editable elsewhere.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from quantum_project import KET_00, evolve_static, rydberg_hamiltonian
from quantum_project.symmetry import (
    compare_states, dicke_hamiltonian, evolve_uniform, exchange_transformation,
    full_hamiltonian, reduced_two_atom_hamiltonian,
)


def main():
    omega = 2 * np.pi * 1e6
    parameters = (omega, .2 * omega, .37, .2 * omega)
    times = np.linspace(0, .25e-6, 41)
    full = evolve_static(KET_00, times, rydberg_hamiltonian(*parameters))
    reduced = evolve_uniform(reduced_two_atom_hamiltonian(*parameters), times[-1], len(times))
    checks = compare_states(full, reduced, 2)
    dark = np.max(abs(full @ exchange_transformation()[:, 2].conj())**2)
    if checks["max_absolute_fidelity_error"] > 1e-12 or checks["max_norm_error"] > 1e-12 or dark > 1e-25:
        raise AssertionError("N=2 milestone comparison failed")
    print("S1: exact symmetry validation (uniform interactions, ground initial state)")
    print(f"N=2 milestone 4x4 vs 3x3: max(1-F)={checks['raw_infidelity']:.3e}; max P_Psi-={dark:.3e}")
    for n in range(2, 6):
        full = evolve_uniform(full_hamiltonian(n, *parameters), times[-1], len(times))
        reduced = evolve_uniform(dicke_hamiltonian(n, *parameters), times[-1], len(times))
        checks = compare_states(full, reduced, n)
        for key in ("max_absolute_fidelity_error", "max_aligned_amplitude_error",
                    "max_excitation_error", "max_ground_probability_error",
                    "max_dicke_population_error", "max_norm_error"):
            if checks[key] > 1e-12:
                raise AssertionError(f"N={n}: {key}={checks[key]}")
        print(f"N={n}: dimensions {2**n} -> {n+1}; max(1-F)={checks['raw_infidelity']:.3e}")
    print("PASS. Small signed discrepancies include floating-point roundoff.")
    print("Explicit timing study: python scripts/symmetry_benchmark.py --help")


if __name__ == "__main__":
    main()
