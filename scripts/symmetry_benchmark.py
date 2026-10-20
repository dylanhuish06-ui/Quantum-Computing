"""Record a controlled benchmark; never run scientific timings in CI.

Imports, warmups, checks, metadata and file writes are
outside timed regions. Plotting belongs exclusively to symmetry_figures.py.
"""

import os
# Set before importing NumPy/SciPy; each benchmark process uses one BLAS thread.
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_variable] = "1"

import argparse
from contextlib import redirect_stdout
import csv
import hashlib
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import platform
import random
import socket
import subprocess
import time
import uuid

import numpy as np
import scipy

from quantum_project.symmetry import csr_bytes, dicke_hamiltonian, evolve_uniform, full_hamiltonian, memory_estimate
from quantum_project.units import mhz_to_angular_frequency

ROOT = Path(__file__).resolve().parents[1]
PARAMETERS = {"omega_rad_s": mhz_to_angular_frequency(1), "delta_over_omega": 0.2,
              "V_over_omega": 0.2, "phase_rad": 0.37, "total_time_s": 0.25e-6, "samples": 41}


def available_memory():
    """Available bytes, considering Linux cgroup limits when present."""
    estimates = []
    if Path("/proc/meminfo").exists():
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                estimates.append(int(line.split()[1]) * 1024)
    try:
        limit = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if limit != "max":
            used = int(Path("/sys/fs/cgroup/memory.current").read_text())
            estimates.append(max(0, int(limit) - used))
    except (OSError, ValueError):
        pass
    return min(estimates) if estimates else None


def physical_core_count():
    """Count exposed Linux socket/core pairs; unavailable is explicitly null."""
    pairs = set()
    for topology in Path("/sys/devices/system/cpu").glob("cpu[0-9]*/topology"):
        try:
            pairs.add(((topology / "physical_package_id").read_text().strip(),
                       (topology / "core_id").read_text().strip()))
        except OSError:
            continue
    return len(pairs) or None


def allocation_plan(full_max, reduced_extra, budget):
    """Keep every requested reduced size even when its full counterpart is unsafe."""
    sizes = sorted(set(range(2, full_max + 1)) | set(reduced_extra))
    full, skipped = [], []
    for n in range(2, full_max + 1):
        estimate = memory_estimate(n, PARAMETERS["samples"], "full")
        if estimate > budget or memory_estimate(n, 2, "full") > 256 * 2**20:
            skipped.append({"N": n, "method": "full", "estimated_working_bytes": estimate,
                            "reason": "exceeds configured/available memory budget or 256 MiB builder limit"})
        else:
            full.append(n)
    reduced = []
    for n in sizes:
        estimate = memory_estimate(n, PARAMETERS["samples"], "reduced")
        if estimate > budget:
            skipped.append({"N": n, "method": "reduced", "estimated_working_bytes": estimate,
                            "reason": "exceeds configured/available memory budget"})
        else:
            reduced.append(n)
    return full, reduced, skipped


def run_trial(N, method, parameters=PARAMETERS):
    builder = full_hamiltonian if method == "full" else dicke_hamiltonian
    omega = parameters["omega_rad_s"]
    arguments = (N, omega, parameters["delta_over_omega"] * omega,
                 parameters["phase_rad"], parameters["V_over_omega"] * omega)
    start = time.perf_counter_ns()
    matrix = builder(*arguments)
    constructed = time.perf_counter_ns()
    states = evolve_uniform(matrix, parameters["total_time_s"], parameters["samples"])
    finished = time.perf_counter_ns()
    # All checks and memory accounting occur after the clocks stop.
    norm_error = float(np.max(abs(np.sum(abs(states)**2, axis=1) - 1)))
    if norm_error > 1e-9:
        raise RuntimeError(f"Propagation norm failure: {norm_error}")
    return {"construction_s": (constructed-start) * 1e-9,
            "evolution_s": (finished-constructed) * 1e-9,
            "total_s": (finished-start) * 1e-9,
            "hamiltonian_bytes": csr_bytes(matrix), "trajectory_bytes": states.nbytes,
            "estimated_working_bytes": memory_estimate(N, parameters["samples"], method),
            "norm_error": norm_error}


def summarize(rows):
    summary = []
    for N, method in sorted({(row["N"], row["method"]) for row in rows}):
        group = [row for row in rows if (row["N"], row["method"]) == (N, method)]
        result = {"N": N, "method": method, "full_dimension": 2**N,
                  "reduced_dimension": N + 1, "repeats": len(group)}
        for phase in ("construction", "evolution", "total"):
            q25, median, q75 = np.percentile([r[f"{phase}_s"] for r in group], [25, 50, 75])
            result.update({f"{phase}_median_s": float(median), f"{phase}_q25_s": float(q25),
                           f"{phase}_q75_s": float(q75), f"{phase}_iqr_s": float(q75-q25)})
        for field in ("hamiltonian_bytes", "trajectory_bytes", "estimated_working_bytes"):
            result[field] = max(r[field] for r in group)
        summary.append(result)
    return summary


def write_rows(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-max", type=int, default=12)
    parser.add_argument("--reduced-extra", type=int, nargs="*", default=[16, 24, 32, 64])
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--memory-budget-mib", type=float, default=256)
    parser.add_argument("--session-kind", choices=["provisional", "final"], default="provisional")
    parser.add_argument("--machine-label", default="")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "data" / "symmetry_provisional")
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--dry-run", action="store_true", help="Validate configuration and list safe sizes, without timing")
    args = parser.parse_args()
    if args.repeats < 10 or args.full_max < 2 or args.full_max > 20 or any(n < 2 or n > 512 for n in args.reduced_extra):
        parser.error("require repeats>=10, 2<=full-max<=20 and 2<=reduced-extra<=512")
    if not np.isfinite(args.memory_budget_mib) or not 0 < args.memory_budget_mib <= 4096:
        parser.error("memory budget must be in (0,4096] MiB")
    if not args.dry_run and any(os.environ.get(key, "").lower() in ("1", "true", "yes") for key in ("CI", "GITHUB_ACTIONS")):
        parser.error("Scientific timing collection is prohibited in CI/GitHub Actions; use --dry-run for correctness checks")
    if args.session_kind == "final" and (not args.machine_label or "provisional" in str(args.output_dir)):
        parser.error("Final runs require --machine-label and a separate non-provisional --output-dir")
    if args.session_kind == "final" and os.environ.get("CODESPACES", "").lower() == "true":
        parser.error("Final scientific timings require a controlled physical machine, not Codespaces")
    available = available_memory()
    budget = min(int(args.memory_budget_mib * 2**20), available // 4 if available is not None else 64 * 2**20)
    full_sizes, sizes, skips = allocation_plan(args.full_max, args.reduced_extra, budget)
    skipped = [r["N"] for r in skips if r["method"] == "full"]
    if not sizes:
        parser.error("Insufficient budget for any reduced calculation")
    if args.dry_run:
        print(json.dumps({"full_sizes": full_sizes, "reduced_sizes": sizes,
                          "skipped_full_sizes": skipped, "skips": skips, "budget_bytes": budget}))
        return
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = [args.output_dir / f"symmetry_benchmark_{name}" for name in ("raw.csv", "summary.csv", "metadata.json")]
    if any(path.exists() for path in paths):
        parser.error("Output session already exists; choose a new directory to preserve raw timings")
    original_affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    pinned = None
    if original_affinity:
        try:
            os.sched_setaffinity(0, {original_affinity[0]})
            pinned = original_affinity[0]
        except OSError:
            pass
    config = io.StringIO()
    with redirect_stdout(config):
        np.show_config()
        scipy.show_config()
    cpu_model = platform.processor()
    if Path("/proc/cpuinfo").exists():
        cpu_model = next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
                          if line.startswith("model name")), cpu_model)
    metadata = {"session_id": str(uuid.uuid4()), "session_kind": args.session_kind,
                "machine_label": args.machine_label or "Codespaces provisional session",
                "started_utc": datetime.now(timezone.utc).isoformat(), "hostname": socket.gethostname(),
                "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                "operating_system": platform.platform(), "cpu_model": cpu_model,
                "logical_cpu_count": os.cpu_count(), "original_affinity": original_affinity, "pinned_cpu": pinned,
                "physical_core_count": physical_core_count(),
                "thread_environment": {k: os.environ[k] for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS")},
                "blas_lapack_configuration": config.getvalue(), "available_memory_at_start": available,
                "memory_budget_bytes": budget, "parameters": PARAMETERS, "repeats": args.repeats,
                "full_sizes": full_sizes, "reduced_sizes": sizes, "skipped_full_sizes": skipped,
                "skips": skips,
                "seed": args.seed, "timer": "time.perf_counter_ns",
                "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
                "source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in (Path(__file__), ROOT / "src/quantum_project/symmetry.py")},
                "order": "Shuffle N each repeat; within each N alternate full/reduced first across repeats, with seeded initial choice",
                "warmup": "One untimed construction+evolution per (N, method) before timed rounds",
                "memory_note": "CSR+trajectory payloads are exact array bytes; estimated working memory is a conservative preflight heuristic, not peak RSS",
                "scope": "Construction and evolution only; no imports, plotting, expansion to full space, file I/O or observable comparisons timed",
                "limitations": "Host sharing and frequency scaling are not controlled by this process; same process/session used for every raw timing"}
    rows = []
    rng = random.Random(args.seed)
    first = {n: rng.randrange(2) for n in sizes}
    unavailable = set()

    def safe_now(n, method, repeat):
        remaining = available_memory()
        estimate = memory_estimate(n, PARAMETERS["samples"], method)
        if remaining is not None and estimate > remaining // 4:
            skips.append({"N": n, "method": method, "repeat": repeat,
                          "estimated_working_bytes": estimate,
                          "reason": "available memory fell below the safety margin"})
            return False
        return True

    try:
        for n in sizes:
            for method in (["full", "reduced"] if n in full_sizes else ["reduced"]):
                if safe_now(n, method, "warmup"):
                    run_trial(n, method)
                else:
                    unavailable.add((n, method))
        print("Warm-up complete; beginning interleaved timed rounds", flush=True)
        for repeat in range(args.repeats):
            order = sizes.copy()
            rng.shuffle(order)
            for n in order:
                methods = ["full", "reduced"] if n in full_sizes else ["reduced"]
                if (repeat + first[n]) % 2:
                    methods.reverse()
                for method in methods:
                    if (n, method) in unavailable or not safe_now(n, method, repeat):
                        continue
                    stamp = datetime.now(timezone.utc).isoformat()
                    measurement = run_trial(n, method)
                    rows.append({"session_id": metadata["session_id"], "session_kind": args.session_kind,
                                 "order": len(rows), "repeat": repeat, "N": n, "method": method,
                                 "full_dimension": 2**n, "reduced_dimension": n+1,
                                 "started_utc": stamp, **measurement})
            print(f"Completed round {repeat + 1}/{args.repeats}", flush=True)
        metadata["complete"] = len(rows) == args.repeats * (len(full_sizes) + len(sizes))
    finally:
        metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if rows:
            write_rows(paths[0], rows)
            write_rows(paths[1], summarize(rows))
        paths[2].write_text(json.dumps(metadata, indent=2) + "\n")
        if original_affinity is not None and pinned is not None:
            os.sched_setaffinity(0, set(original_affinity))
    print(f"Saved {len(rows)} raw timings to {args.output_dir}")


if __name__ == "__main__":
    main()
