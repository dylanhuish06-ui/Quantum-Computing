"""Cache corruption, interrupted experiments and resume are behavioral contracts."""
import json
from pathlib import Path
import subprocess
from unittest.mock import Mock

import pytest

from quantum_project.ancilla import run


def fake_experiment(command, **kwargs):
    directory = Path(command[command.index("--output-dir")+1])
    gate = int(command[command.index("--gate")+1]) if "--gate" in command else 5
    (directory/"validation.json").write_text(json.dumps({"passed": True, "gate": gate}))
    (directory/"residuals.csv").write_text("error\n0\n")
    if gate == 5:
        (directory/"adaptive.csv").write_text("m\n0\n")
        (directory/"observables.csv").write_text("tau\n0\n")
    else:
        (directory/"trajectories.npz").write_bytes(b"fake artifact for cache test")


@pytest.fixture
def mocked_experiment(monkeypatch):
    spec = {"schema": 1, "gate": 5, "parameter": 1}
    monkeypatch.setattr(run, "specification", lambda gate: spec|{"gate": gate})
    child = Mock(side_effect=fake_experiment)
    monkeypatch.setattr(run.subprocess, "run", child)
    monkeypatch.setattr(run.subprocess, "check_output", lambda *args, **kwargs: "test-source-snapshot\n")
    return spec, child


def test_resume_skips_only_intact_matching_experiment(tmp_path, mocked_experiment):
    spec, child = mocked_experiment
    directory, action = run.execute_gate(5, tmp_path, resume=True, timeout=10)
    assert action == "ran" and run.valid_cache(directory, spec)
    snapshot = {p.name: p.stat().st_mtime_ns for p in directory.iterdir()}
    assert run.execute_gate(5, tmp_path, resume=True, timeout=10) == (directory, "skipped")
    assert child.call_count == 1
    assert snapshot == {p.name: p.stat().st_mtime_ns for p in directory.iterdir()}
    (directory/"observables.csv").write_text("corrupted")
    assert not run.valid_cache(directory, spec)
    assert run.execute_gate(5, tmp_path, resume=True, timeout=10)[1] == "ran"
    assert child.call_count == 2
    assert len(list(tmp_path.glob("*.invalid-*"))) == 1
    assert run.valid_cache(directory, spec)


def test_changed_parameters_get_new_directory(tmp_path, mocked_experiment):
    spec, child = mocked_experiment
    first, _ = run.execute_gate(5, tmp_path, resume=True, timeout=10)
    spec["parameter"] = 2
    second, _ = run.execute_gate(5, tmp_path, resume=True, timeout=10)
    assert first != second and first.exists() and child.call_count == 2
    with pytest.raises(ValueError, match="preserved"):
        run.execute_gate(5, tmp_path, resume=False, timeout=10)


def test_failed_task_preserves_completed_artifacts(tmp_path, mocked_experiment):
    _, child = mocked_experiment
    first, _ = run.execute_gate(1, tmp_path, resume=True, timeout=10)
    child.side_effect = subprocess.TimeoutExpired("experiment", 10)
    with pytest.raises(subprocess.TimeoutExpired):
        run.execute_gate(2, tmp_path, resume=True, timeout=10)
    assert first.exists()
    assert not list(tmp_path.glob("gate2-*/manifest.json"))
    assert list(tmp_path.glob("*.pending-*"))


def test_missing_report_artifacts_never_skipped(tmp_path, mocked_experiment):
    spec, _ = mocked_experiment
    directory, _ = run.execute_gate(5, tmp_path, resume=True, timeout=10)
    manifest = json.loads((directory/"manifest.json").read_text())
    del manifest["artifact_sha256"]["observables.csv"]
    (directory/"manifest.json").write_text(json.dumps(manifest))
    assert not run.valid_cache(directory, spec)


def test_run_lock_and_cli_validation(tmp_path):
    with run.run_lock(tmp_path):
        with pytest.raises(FileExistsError):
            with run.run_lock(tmp_path):
                pass
    assert (tmp_path/".run.lock").read_text() == ""
    with pytest.raises(SystemExit):
        run.main(["--gates", "5"])
    with pytest.raises(SystemExit):
        run.main(["--max-seconds-per-gate", "nan"])


def test_main_records_resume_outcome(tmp_path, mocked_experiment):
    args = ["--gates", "5", "--no-figures", "--resume", "--output-dir", str(tmp_path)]
    assert run.main(args) == 0
    assert run.main(args) == 0
    assert json.loads((tmp_path/"latest.json").read_text())["actions"] == {"5": "skipped"}
