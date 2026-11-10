"""Resumable, content-verified Phase 4 experiments.

From the repository: PYTHONPATH=src python -m quantum_project.ancilla.run --resume
Each gate is a bounded foreground subprocess. A validated manifest is published
atomically only after its experiment succeeds. Changed inputs get a new cache
key; damaged outputs are preserved and rebuilt with --resume.
"""
import argparse
import fcntl
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import uuid

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[3]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def environment():
    import matplotlib
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
            "platform": platform.platform(), "threads": 1}


def specification(gate):
    script = ROOT / "scripts/ancilla" / ("validate_gate5.py" if gate == 5 else "validate_gates.py")
    sources = [p for p in sorted((ROOT/"src/quantum_project/ancilla").glob("*.py"))
               if p.name not in ("run.py", "figures.py")]+[script]
    return {"schema": 1, "gate": gate, "environment": environment(),
            "parameters": {"profile": "canonical validated gate recipe", "gate": gate},
            "source_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sources}}


def valid_cache(directory, spec):
    try:
        manifest = json.loads((directory/"manifest.json").read_text())
        if manifest["input_hash"] != digest(spec) or manifest["specification"] != spec:
            return False
        expected = {"validation.json", "residuals.csv"}
        expected |= {"adaptive.csv", "observables.csv"} if spec["gate"] == 5 else {"trajectories.npz"}
        hashes = manifest["artifact_sha256"]
        if not expected.issubset(hashes):
            return False
        for name, checksum in hashes.items():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts or sha256(directory/relative) != checksum:
                return False
        actual = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()}
        if actual != set(hashes)|{"manifest.json"}:
            return False
        report = json.loads((directory/"validation.json").read_text())
        return report["passed"] is True and report["gate"] == spec["gate"]
    except (OSError, ValueError, KeyError, TypeError):
        return False


def gate_command(gate, output):
    script = "validate_gate5.py" if gate == 5 else "validate_gates.py"
    command = [sys.executable, str(ROOT/"scripts/ancilla"/script), "--output-dir", str(output)]
    if gate != 5:
        command += ["--gate", str(gate)]
    return command


def execute_gate(gate, root, *, resume, timeout):
    spec = specification(gate)
    key = digest(spec)
    directory = root/f"gate{gate}-{key[:16]}"
    if directory.exists():
        if resume and valid_cache(directory, spec):
            print(f"Gate {gate}: verified cache, skipped", flush=True)
            return directory, "skipped"
        if not resume:
            raise ValueError(f"existing outputs preserved: {directory}; use --resume")
        directory.rename(root/(directory.name+".invalid-"+uuid.uuid4().hex[:8]))
    pending = root/(directory.name+".pending-"+uuid.uuid4().hex[:8])
    pending.mkdir()
    print(f"Gate {gate}: running", flush=True)
    env = {**os.environ, "PYTHONPATH": str(ROOT/"src"), "MPLBACKEND": "Agg",
           "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    # Logs live outside the experiment folder until success; old gate scripts
    # intentionally reject nonempty output directories.
    log = root/(pending.name+".log")
    with log.open("w") as stream:
        subprocess.run(gate_command(gate, pending), cwd=ROOT, env=env,
                       stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    log.rename(pending/"execution.log")
    report = json.loads((pending/"validation.json").read_text())
    if report.get("passed") is not True or report.get("gate") != gate:
        raise ValueError(f"gate {gate} did not validate")
    manifest = {"input_hash": key, "specification": spec,
                "local_git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "local_git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True),
                "utc": datetime.now(timezone.utc).isoformat(),
                "artifact_sha256": {str(p.relative_to(pending)): sha256(p) for p in pending.rglob("*") if p.is_file()}}
    (pending/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    if not valid_cache(pending, spec):
        raise ValueError("new experiment manifest failed validation")
    pending.rename(directory)
    return directory, "ran"


@contextmanager
def run_lock(root):
    lock = root/".run.lock"
    with lock.open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise FileExistsError(f"another ancilla runner holds {lock}") from error
        try:
            stream.seek(0)
            stream.truncate()
            stream.write(str(os.getpid()))
            stream.flush()
            yield
        finally:
            stream.seek(0)
            stream.truncate()
            stream.flush()
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/ancilla/runs")
    parser.add_argument("--gates", type=int, nargs="+", choices=range(1, 6), default=[1, 2, 3, 4, 5])
    parser.add_argument("--no-figures", action="store_true", help="allow a selected subset of experiments")
    parser.add_argument("--max-seconds-per-gate", type=float, default=300.)
    args = parser.parse_args(argv)
    if not np.isfinite(args.max_seconds_per_gate) or args.max_seconds_per_gate <= 0:
        parser.error("max-seconds-per-gate must be positive and finite")
    if not args.no_figures and set(args.gates) != set(range(1, 6)):
        parser.error("all five gates are required for figures; use --no-figures for a subset")
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    try:
        with run_lock(root):
            directories, actions = {}, {}
            for gate in sorted(set(args.gates)):
                directories[gate], actions[gate] = execute_gate(gate, root, resume=args.resume,
                                                              timeout=args.max_seconds_per_gate)
            figures = None
            if not args.no_figures:
                from .figures import build_figures
                figures = build_figures(directories, root, resume=args.resume)
            summary = {"gates": {str(k): str(v.relative_to(root)) for k, v in directories.items()},
                       "actions": actions, "figures": str(figures.relative_to(root)) if figures else None}
            temporary = root/(".latest-"+uuid.uuid4().hex+".json")
            temporary.write_text(json.dumps(summary, indent=2)+"\n")
            temporary.replace(root/"latest.json")
            print(json.dumps(summary, indent=2))
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        print(f"Ancilla run stopped; completed artifacts preserved: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
