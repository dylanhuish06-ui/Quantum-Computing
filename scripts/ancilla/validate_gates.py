"""Bounded, reproducible numerical exit checks for Phase 4 gates 1--4.

Run from any directory. Frequencies use h=Hbar/Omega_ref and tau=Omega_ref*t.
Existing output directories are preserved. JSON is written last, only on success.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
import scipy
from scipy.linalg import expm
from quantum_project.ancilla import (
    FullParameters, HybridParameters, dicke_isometry_bulk, evolve,
    full_hamiltonian_reference, hybrid_hamiltonian,
)
from quantum_project.ancilla.protection import (
    gap_penalty, permutation_average, permutation_cycles, symmetric_group,
)
from quantum_project.ancilla.syndrome import selective_evolution, symmetric_projector

OMEGA_REF = 2*np.pi*1e6
SOLVER = {"dt": .04, "rtol": 2e-12, "atol": 2e-14}


def norm_error(states):
    return float(np.max(abs(np.sum(abs(states)**2, axis=1)-1)))


def amplitude_error(a, b):
    return float(np.max(np.linalg.norm(a-b, axis=-1)))


def broken_profile(N):
    """All coefficients are dimensionless; inhomogeneity has zero-mean detuning."""
    def parameters(tau):
        return FullParameters(
            np.full(N, 1+.2*np.sin(2*tau)),
            -.4+.3*tau+(.65+.1*np.sin(tau))*np.linspace(-1, 1, N),
            np.full(N, .23+.3*tau),
            .17*(1+.1*np.cos(tau))*(np.ones((N, N))-np.eye(N)))
    return lambda tau: full_hamiltonian_reference(parameters(tau))


def hybrid_profile(N, m):
    rng = np.random.default_rng(841+10*N+m)
    omega, delta, phase, cross = rng.normal(size=(4, m))
    pairs = np.triu(rng.normal(size=(m, m)), 1)
    pairs += pairs.T
    def parameters(tau):
        return HybridParameters(N, m, 1+.3*np.sin(2*tau), -.4+tau*tau,
                                .2+.6*tau, .12*np.cos(tau), omega*(1+.2*tau),
                                delta+tau, phase-.7*tau, cross*(1+.1*tau), pairs*np.cos(tau))
    return parameters


def gate1():
    rows, arrays = [], {}
    for N in range(2, 7):
        H = broken_profile(N)
        B = dicke_isometry_bulk(N, 0).toarray()
        P = B@B.conj().T
        initial = np.zeros(2**N, complex)
        initial[0] = 1
        times = np.linspace(0, 1.2, 13)
        ode = evolve(H, initial, times, **SOLVER)
        tighter = evolve(H, initial, times, dt=.02, rtol=3e-14, atol=1e-15)
        errors = []
        for dt in (.04, .02, .01):
            mid = evolve(H, initial, times, dt=dt, method="midpoint")
            errors.append(amplitude_error(mid, tighter))
        # Outputs 0.1 apart make dt=.04 become .1/3; use actual refinements
        # .1/3 -> .1/5 -> .1/10; only the last pair has exact factor two.
        dense = evolve(H, initial, times, dt=.01, method="dense_midpoint")
        static = H(0).toarray()
        exact = expm(-1j*1.2*static)@initial
        constant = evolve(lambda t: static, initial, [0, 1.2], **SOLVER)[-1]
        row = {"N": N, "ode_refinement_l2": amplitude_error(ode, tighter),
               "mid_coarse_l2": errors[0], "mid_medium_l2": errors[1], "mid_fine_l2": errors[2],
               "midpoint_order": float(np.log2(errors[1]/errors[2])),
               "dense_sparse_l2": amplitude_error(mid, dense),
               "constant_dense_l2": amplitude_error(constant, exact),
               "norm_drift": norm_error(ode),
               "projector_idempotence_max": float(np.max(abs(P@P-P))),
               "QHP_spectral_norm": float(np.linalg.norm((np.eye(2**N)-P)@static@B, 2))}
        assert row["ode_refinement_l2"] < 1e-10
        assert row["dense_sparse_l2"] < 1e-12
        assert row["constant_dense_l2"] < 1e-10
        assert row["norm_drift"] < 1e-10
        assert 1.95 < row["midpoint_order"] < 2.05
        rows.append(row)
        arrays.update({f"N{N}_times_tau": times, f"N{N}_ode": ode, f"N{N}_midpoint": mid,
                       f"N{N}_tighter": tighter})
    return rows, arrays, {"N": list(range(2, 7)), "T_tau": 1.2,
                          "profile": "broken_profile", "midpoint_dt": [.04, .02, .01]}


def gate2():
    rows, arrays = [], {}
    T = 1.6
    for N in range(2, 7):
        H = broken_profile(N)
        B = dicke_isometry_bulk(N, 0).toarray()
        P = symmetric_projector(N)
        initial = np.zeros(2**N, complex)
        initial[0] = 1
        for checks in (8, 32, 128):
            times = np.linspace(0, T, checks+1)
            target = evolve(lambda t: B.conj().T@H(t)@B, np.eye(N+1)[0], times, **SOLVER)@B.T
            bare = evolve(H, initial, times, **SOLVER)
            result = selective_evolution(H, initial, times, P, **SOLVER)
            # Separate integration refinement from the fixed measurement schedule.
            fine = selective_evolution(H, initial, times, P, dt=.02, rtol=3e-14, atol=1e-15)
            mid = selective_evolution(H, initial, times, P, dt=(T/checks)/4, method="midpoint")
            finer_mid = selective_evolution(H, initial, times, P, dt=(T/checks)/8, method="midpoint")
            errors = [amplitude_error(x.branch, result.branch) for x in (mid, finer_mid)]
            stats = result.metrics(target, P)
            raw_leak = np.sum(abs(bare-bare@P.T)**2, axis=1)
            row = {"N": N, "checks": checks, "interval_tau": T/checks,
                   "conditional_fidelity": float(stats["conditional_fidelity"][-1]),
                   "success_probability": float(result.success[-1]),
                   "successful_target_probability": float(stats["successful_target_probability"][-1]),
                   "bare_target_probability": float(abs(np.vdot(target[-1], bare[-1]))**2),
                   "bare_leakage": float(raw_leak[-1]),
                   "conditional_target_l2": amplitude_error(result.conditional[-1], target[-1]),
                   "ode_refinement_l2": amplitude_error(fine.branch, result.branch),
                   "midpoint_branch_l2": errors[-1], "midpoint_order": float(np.log2(errors[0]/errors[1])),
                   "oracle_calls": checks, "ancilla_preparations": checks,
                   "between_check_reuses": checks-1}
            assert row["ode_refinement_l2"] < 1e-10
            assert 1.95 < row["midpoint_order"] < 2.05
            assert 0 < row["success_probability"] <= 1+1e-12
            assert abs(row["successful_target_probability"]-row["conditional_fidelity"]*row["success_probability"]) < 1e-13
            assert np.max(np.diff(result.success)) < 1e-12
            rows.append(row)
            key = f"N{N}_checks{checks}"
            arrays.update({key+"_times_tau": times, key+"_branch": result.branch,
                           key+"_success": result.success, key+"_target": target, key+"_bare": bare})
        selected = [r for r in rows if r["N"] == N]
        assert selected[-1]["conditional_target_l2"] < selected[0]["conditional_target_l2"]
        assert selected[-1]["success_probability"] > selected[0]["success_probability"]
    return rows, arrays, {"N": list(range(2, 7)), "T_tau": T, "checks": [8, 32, 128],
                          "profile": "broken_profile", "oracle": "ideal global controlled projector",
                          "reset": "prepare |0> per check; successful branch already returns |0>",
                          "compiled_two_qubit_gates": None, "compiled_depth": None}


def gate3():
    rows, arrays = [], {}
    T = 1.6
    for N in range(2, 7):
        H = broken_profile(N)
        B = dicke_isometry_bulk(N, 0).toarray()
        P = B@B.conj().T
        initial = np.zeros(2**N, complex)
        initial[0] = 1
        times = np.linspace(0, T, 33)
        target = evolve(lambda t: B.conj().T@H(t)@B, np.eye(N+1)[0], times, **SOLVER)@B.T
        for G in (0., .5, 2., 8.):
            penalty = gap_penalty(N, G)
            controlled = lambda t: H(t)+penalty
            states = evolve(controlled, initial, times, **SOLVER)
            fine = evolve(controlled, initial, times, dt=.02, rtol=3e-14, atol=1e-15)
            # Independent exponential action at a fine time-discretisation grid.
            mid = evolve(controlled, initial, times, dt=.0025, method="midpoint")
            leakage = np.sum(abs(states-states@P.T)**2, axis=1)
            target_probability = abs(np.einsum("ti,ti->t", target.conj(), states))**2
            row = {"N": N, "G_over_Omega_ref": G, "first_gap_over_Omega_ref": N*G,
                   "final_leakage": float(leakage[-1]), "max_leakage": float(leakage.max()),
                   "target_probability": float(target_probability[-1]), "norm_drift": norm_error(states),
                   "ode_refinement_l2": amplitude_error(states, fine),
                   "midpoint_reference_l2": amplitude_error(states, mid),
                   "penalty_on_Dicke_max": float(np.max(abs(penalty@B)))}
            assert row["norm_drift"] < 1e-10
            assert row["ode_refinement_l2"] < 1e-10
            assert row["midpoint_reference_l2"] < 3e-6
            assert row["penalty_on_Dicke_max"] < 2e-12
            rows.append(row)
            key = f"N{N}_G{G}"
            arrays.update({key+"_times_tau": times, key+"_states": states,
                           key+"_target_probability": target_probability, key+"_leakage": leakage})
        selected = [r for r in rows if r["N"] == N]
        assert selected[-1]["max_leakage"] < selected[0]["max_leakage"]/100
    # Finite-cycle experiment is deliberately small and independently compiled.
    H = broken_profile(3)(.1)
    group = symmetric_group(3)
    average = permutation_average(H, group)
    initial = np.eye(8)[0]
    target = expm(-1j*.7*average.toarray())@initial
    dynamical = []
    for cycles in (1, 2, 4, 8):
        result = permutation_cycles(lambda t: H, initial, .7, cycles, group, method="midpoint")
        dynamical.append({"cycles": cycles, "state_l2_error": amplitude_error(result["states"][-1], target),
                          "swap_count": result["swap_count"], "cnot_equivalent_count": result["cnot_equivalent_count"]})
    assert dynamical[-1]["state_l2_error"] < dynamical[0]["state_l2_error"]
    return rows, arrays, {"N": list(range(2, 7)), "T_tau": T, "G_over_Omega_ref": [0, .5, 2, 8],
                          "profile": "broken_profile", "success_probability": 1,
                          "penalty": "G*(Jmax*(Jmax+1)*I-J2); ideal all-to-all coherent control",
                          "dynamical_static_N3": dynamical,
                          "dynamical_limitations": "instantaneous pulses; finite pulses separately unit tested; no hardware fidelity claim"}


def gate4():
    rows, arrays = [], {}
    for N in range(2, 9):
        for m in range(min(3, N)+1):
            profile = hybrid_profile(N, m)
            B = dicke_isometry_bulk(N, m).toarray()
            d = B.shape[1]
            rng = np.random.default_rng(111+N+m)
            initial = rng.normal(size=d)+1j*rng.normal(size=d)
            initial /= np.linalg.norm(initial)
            times = np.linspace(0, .6, 13)
            full_H = lambda t: full_hamiltonian_reference(profile(t).full_parameters())
            red_H = lambda t: hybrid_hamiltonian(profile(t))
            full = evolve(full_H, B@initial, times, **SOLVER)
            reduced = evolve(red_H, initial, times, **SOLVER)
            expanded = reduced@B.T
            errors = []
            for dt in (.025, .0125, .00625):
                mid = evolve(red_H, initial, times, method="midpoint", dt=dt)
                errors.append(amplitude_error(mid@B.T, full))
            # Richardson extrapolation checks the *complex* phase-sensitive state.
            medium = evolve(red_H, initial, times, method="midpoint", dt=.0125)
            extrapolated = (4*mid-medium)/3
            bits = ((np.arange(2**N)[:, None] >> np.arange(N-1, -1, -1)) & 1)
            pop_full, pop_red = abs(full)**2@bits, abs(expanded)**2@bits
            corr_full = abs(full)**2@(bits[:, 0]*bits[:, -1])
            corr_red = abs(expanded)**2@(bits[:, 0]*bits[:, -1])
            h, H = red_H(.137).toarray(), full_H(.137).toarray()
            row = {"N": N, "m": m, "dimension": d,
                   "compression_max": float(np.max(abs(B.conj().T@H@B-h))),
                   "invariance_max": float(np.max(abs(H@B-B@h))),
                   "complex_trajectory_l2": amplitude_error(full, expanded),
                   "population_max": float(np.max(abs(pop_full-pop_red))),
                   "pair_correlation_max": float(np.max(abs(corr_full-corr_red))),
                   "mid_coarse_l2": errors[0], "mid_medium_l2": errors[1], "mid_fine_l2": errors[2],
                   "midpoint_order": float(np.log2(errors[1]/errors[2])),
                   "extrapolated_l2": amplitude_error(extrapolated@B.T, full),
                   "norm_drift": max(norm_error(full), norm_error(reduced))}
            assert max(row[k] for k in ("compression_max", "invariance_max")) < 1e-12
            assert max(row[k] for k in ("complex_trajectory_l2", "population_max", "pair_correlation_max", "norm_drift")) < 1e-10
            assert 1.95 < row["midpoint_order"] < 2.05
            assert row["extrapolated_l2"] < errors[-1]/100
            rows.append(row)
            key = f"N{N}_m{m}"
            arrays.update({key+"_times_tau": times, key+"_full": full, key+"_reduced": reduced,
                           key+"_population_full": pop_full, key+"_population_reduced": pop_red})
    return rows, arrays, {"N": list(range(2, 9)), "m": [0, 1, 2, 3], "T_tau": .6,
                          "profile": "hybrid_profile", "seed_parameters": "841+10*N+m", "seed_states": "111+N+m",
                          "midpoint_dt": [.025, .0125, .00625], "initial": "random complex unit state in hybrid basis"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=int, choices=(1, 2, 3, 4), required=True)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--max-seconds", type=int, default=240)
    args = parser.parse_args()
    if args.max_seconds <= 0:
        parser.error("max-seconds must be positive")
    output = args.output_dir or ROOT / "results" / "ancilla" / f"gate{args.gate}"
    if output.exists() and any(output.iterdir()):
        parser.error(f"preserving existing artifacts: {output}")
    def timeout(*_):
        raise TimeoutError("gate validation exceeded the session runtime budget")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(args.max_seconds)
    started = time.monotonic()
    rows, arrays, parameters = globals()[f"gate{args.gate}"]()
    output.mkdir(parents=True, exist_ok=True)
    with (output / "residuals.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(output / "trajectories.npz", **arrays)
    paths = sorted((ROOT / "src/quantum_project/ancilla").glob("*.py"))+[Path(__file__).resolve()]
    report = {"schema": 1, "gate": args.gate, "passed": True,
              "utc": datetime.now(timezone.utc).isoformat(),
              "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
              "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in paths},
              "artifact_sha256": {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()},
              "parameters": parameters, "solver": SOLVER,
              "units": {"Omega_ref_rad_per_s": OMEGA_REF, "tau": "Omega_ref*t_seconds", "h": "Hbar/Omega_ref",
                        "state_errors": "dimensionless L2, no phase alignment or renormalisation",
                        "QHP_spectral_norm": "in units of Omega_ref"},
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                              "platform": platform.platform(), "threads": {n: os.environ[n] for n in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}},
              "limits": {"max_seconds": args.max_seconds, "full_N_max": 12, "dense_expm_dimension_max": 256},
              "provisional_codespaces_elapsed_seconds": time.monotonic()-started,
              "rows": rows}
    (output / "validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    signal.alarm(0)
    print(json.dumps({"gate": args.gate, "passed": True, "rows": len(rows), "output": str(output),
                      "elapsed_seconds": report["provisional_codespaces_elapsed_seconds"]}))


if __name__ == "__main__":
    main()
