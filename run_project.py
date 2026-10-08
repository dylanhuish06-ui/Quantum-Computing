"""Convenient entry point for the existing Chapter 5 scripts and tests."""

import argparse
from pathlib import Path
import shlex
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
COMMANDS = ("tests", "fig51", "fig52", "milestone", "results", "zeno", "symmetry")
MENU = """
PHYS3561 Chapter 5 — Quantum Dynamics Project

1. Run test suite
2. Generate Figure 5.1
3. Generate Figure 5.2 milestone
4. Run complete milestone
5. Show generated result files
6. Quantum Zeno extension
7. Symmetry-scaling extension
0. Exit
"""


def build_command(command):
    """Build a subprocess command using this launcher's Python interpreter."""
    arguments = {
        "tests": ["-m", "pytest", "-q"],
        "fig51": [str(ROOT / "scripts/figure_5_1.py")],
        "fig52": [str(ROOT / "scripts/figure_5_2.py"), "--milestone",
                  "--interaction-mhz", "0", "0.5", "1", "5", "100", "--phase", "0"],
    }
    return [sys.executable, *arguments[command]]


def run_command(command):
    args = build_command(command)
    print(f"Running: {shlex.join(args)}", flush=True)
    result = subprocess.run(args, cwd=ROOT, check=False)
    if result.returncode:
        print(f"{command} failed (return code {result.returncode}).", file=sys.stderr)
    # POSIX signals produce negative subprocess codes; expose a shell exit code.
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


def discover_results(directory):
    """Return sorted files, including any nested output files, without modifying them."""
    return sorted(path for path in directory.rglob("*") if path.is_file())


def show_results():
    for name in ("figures", "data"):
        directory = ROOT / "results" / name
        print(f"\n{directory}/")
        if not directory.is_dir():
            print("  Directory is absent.")
            continue
        files = discover_results(directory)
        if not files:
            print("  No files found (empty).")
        for path in files:
            print(f"  {path.relative_to(directory)}")


def dispatch(command):
    if command == "milestone":
        for step in ("tests", "fig51", "fig52"):
            status = run_command(step)
            if status:
                return status
        print("Milestone complete. Generated outputs:")
        show_results()
        return 0
    if command == "results":
        show_results()
        return 0
    if command in ("zeno", "symmetry"):
        label = "Quantum Zeno" if command == "zeno" else "Symmetry-scaling"
        print(f"{label} extension has not yet been implemented.")
        return 0
    return run_command(command)


def interactive():
    while True:
        print(MENU)
        choice = input("Choose an option: ").strip()
        if choice == "0":
            return 0
        if choice not in tuple(str(number) for number in range(1, 8)):
            print("Invalid option. Enter a number from 0 to 7.")
            continue
        command = COMMANDS[int(choice) - 1]
        if command in ("fig52", "milestone"):
            print("Figure 5.2 uses the validated 131,072-interval production grid "
                  "and may take some time.")
            while True:
                answer = input("Continue? [y/n]: ").strip().lower()
                if answer in ("y", "n"):
                    break
                print("Please enter y or n.")
            if answer == "n":
                print("Cancelled.")
                continue
        status = dispatch(command)
        if status:
            return status


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run Chapter 5 project tasks; omit the command for a menu.",
        epilog="fig52 and milestone run without confirmation in command-line mode. "
               "zeno and symmetry are reserved for future extensions.",
    )
    parser.add_argument("command", nargs="?", choices=COMMANDS,
                        help="project task to run")
    args = parser.parse_args(argv)
    try:
        return dispatch(args.command) if args.command else interactive()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130
    except EOFError:
        print("\nInput closed. Exiting.")
        return 0
    except OSError as error:
        print(f"Launcher error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
