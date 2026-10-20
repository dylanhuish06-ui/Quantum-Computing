"""Exercise launcher orchestration without running scientific simulations."""

from pathlib import Path
from unittest.mock import Mock, call
import sys

import pytest

import run_project as launcher


@pytest.mark.parametrize("command", launcher.COMMANDS)
def test_known_command_dispatches_without_prompt(command, monkeypatch):
    dispatch = Mock(return_value=0)
    monkeypatch.setattr(launcher, "dispatch", dispatch)
    monkeypatch.setattr("builtins.input", Mock(side_effect=AssertionError("prompted")))
    assert launcher.main([command]) == 0
    dispatch.assert_called_once_with(command)


def test_unknown_command_rejected(monkeypatch, capsys):
    dispatch = Mock()
    monkeypatch.setattr(launcher, "dispatch", dispatch)
    with pytest.raises(SystemExit) as error:
        launcher.main(["unknown"])
    assert error.value.code == 2
    assert "invalid choice" in capsys.readouterr().err
    dispatch.assert_not_called()


def test_help(capsys):
    with pytest.raises(SystemExit) as error:
        launcher.main(["--help"])
    assert error.value.code == 0
    assert "milestone" in capsys.readouterr().out


@pytest.mark.parametrize(("command", "arguments"), [
    ("tests", ["-m", "pytest", "-q"]),
    ("symmetry", [str(launcher.ROOT / "scripts/symmetry_demo.py")]),
    ("fig51", [str(launcher.ROOT / "scripts/figure_5_1.py")]),
    ("fig52", [str(launcher.ROOT / "scripts/figure_5_2.py"), "--milestone",
               "--interaction-mhz", "0", "0.5", "1", "5", "100", "--phase", "0"]),
])
def test_commands_use_active_interpreter_and_repository_cwd(
        command, arguments, monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    run = Mock(return_value=Mock(returncode=0))
    monkeypatch.setattr(launcher.subprocess, "run", run)
    assert launcher.main([command]) == 0
    run.assert_called_once_with([sys.executable, *arguments],
                                cwd=Path(launcher.__file__).resolve().parent, check=False)


@pytest.mark.parametrize(("codes", "steps"), [
    ([5], ["tests"]),
    ([0, 6], ["tests", "fig51"]),
    ([0, 0, 7], ["tests", "fig51", "fig52"]),
    ([0, 0, 0], ["tests", "fig51", "fig52"]),
])
def test_milestone_sequence_stops_at_failure(codes, steps, monkeypatch):
    run = Mock(side_effect=[Mock(returncode=code) for code in codes])
    results = Mock()
    monkeypatch.setattr(launcher.subprocess, "run", run)
    monkeypatch.setattr(launcher, "show_results", results)
    assert launcher.main(["milestone"]) == codes[-1]
    assert run.call_args_list == [
        call(launcher.build_command(step), cwd=launcher.ROOT, check=False) for step in steps
    ]
    assert results.call_count == (1 if codes[-1] == 0 else 0)


def test_results_discovery_and_listing(tmp_path, monkeypatch, capsys):
    figures = tmp_path / "results/figures"
    figures.mkdir(parents=True)
    (figures / "nested").mkdir()
    files = [figures / "z.png", figures / "a.pdf", figures / "nested/b.png"]
    for path in files:
        path.write_bytes(b"existing result")
    monkeypatch.setattr(launcher, "ROOT", tmp_path)
    assert launcher.discover_results(figures) == sorted(files)
    assert launcher.main(["results"]) == 0
    output = capsys.readouterr().out
    for path in files:
        assert str(path.relative_to(figures)) in output
        assert path.read_bytes() == b"existing result"
    assert "Directory is absent" in output
    data = tmp_path / "results/data"
    assert not data.exists()
    data.mkdir()
    assert launcher.main(["results"]) == 0
    assert "No files found (empty)" in capsys.readouterr().out


@pytest.mark.parametrize("command", ["zeno"])
def test_extensions_do_not_run_subprocesses(command, monkeypatch, capsys):
    run = Mock()
    monkeypatch.setattr(launcher.subprocess, "run", run)
    assert launcher.main([command]) == 0
    assert "not yet been implemented" in capsys.readouterr().out
    run.assert_not_called()


@pytest.mark.parametrize(("choice", "command"), [("3", "fig52"), ("4", "milestone")])
@pytest.mark.parametrize("answer", ["y", "n"])
def test_expensive_menu_confirmation(choice, command, answer, monkeypatch, capsys):
    dispatch = Mock(return_value=0)
    monkeypatch.setattr(launcher, "dispatch", dispatch)
    monkeypatch.setattr("builtins.input", Mock(side_effect=[choice, "invalid", answer, "0"]))
    assert launcher.main([]) == 0
    assert "131,072-interval" in capsys.readouterr().out
    if answer == "y":
        dispatch.assert_called_once_with(command)
    else:
        dispatch.assert_not_called()


@pytest.mark.parametrize(("choice", "command"), [
    ("1", "tests"), ("2", "fig51"), ("5", "results"), ("6", "zeno"), ("7", "symmetry"),
])
def test_menu_recovers_from_invalid_input(choice, command, monkeypatch, capsys):
    dispatch = Mock(return_value=0)
    monkeypatch.setattr(launcher, "dispatch", dispatch)
    monkeypatch.setattr("builtins.input", Mock(side_effect=["bad", "", choice, "0"]))
    assert launcher.main([]) == 0
    assert "Invalid option" in capsys.readouterr().out
    dispatch.assert_called_once_with(command)


def test_menu_exits_on_subprocess_failure(monkeypatch, capsys):
    monkeypatch.setattr("builtins.input", Mock(side_effect=["1"]))
    monkeypatch.setattr(launcher.subprocess, "run", Mock(return_value=Mock(returncode=8)))
    assert launcher.main([]) == 8
    assert "tests failed (return code 8)" in capsys.readouterr().err


@pytest.mark.parametrize(("error", "status"), [(KeyboardInterrupt(), 130), (EOFError(), 0)])
def test_menu_closed_cleanly(error, status, monkeypatch):
    monkeypatch.setattr("builtins.input", Mock(side_effect=error))
    assert launcher.main([]) == status


@pytest.mark.parametrize(("result", "status"), [
    (KeyboardInterrupt(), 130), (OSError("cannot start"), 1), (Mock(returncode=-2), 130),
])
def test_subprocess_errors_exit_cleanly(result, status, monkeypatch):
    monkeypatch.setattr(launcher.subprocess, "run", Mock(side_effect=[result]))
    assert launcher.main(["tests"]) == status
