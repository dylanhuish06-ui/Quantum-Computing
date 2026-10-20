# Symmetry-scaling implementation resumption — 2026-10-11

Resumed the existing `extension/symmetry-scaling` worktree at
`/workspaces/Quantum-Computing/symmetry-scaling`, based on launcher commit
`090f5dc`. The mathematical implementation, benchmark tools and launcher
integration were already present as uncommitted work. The parent Quantum Zeno
worktree was left untouched. No commit, push or tag change was made.

The first complete suite reproduced one failure, with 344 passing tests.
`test_reduced_sizes_continue_when_no_full_calculation_fits` requested 0.012 MiB,
which allows the N=2 full calculation (10,688 estimated bytes). Its expectation
of no full sizes was therefore incorrect. The budget is now 0.009 MiB (9,437
bytes): enough for N=2 reduced (7,952 bytes), below both N=2 full and N=3 reduced
(10,496 bytes). Additional assertions verify the intended memory boundary.
The allocation algorithm and its safety limits did not need adjustment.

The launcher demo now checks absolute fidelity error, including deviations
above one, and phase-aligned amplitudes alongside populations and norms.
The independent global-phase test also verifies the new public aligned-amplitude
and absolute-fidelity metrics. Documentation now matches the direct N=2
construction, checkout-local imports, active launcher command and signed S4 plot.

The untimed study regenerated the agreement and broken-symmetry CSV/NPZ files
and validation JSON. The plotting script validated their session linkage and
regenerated all seven S1–S6 plots in PNG/PDF, including the S2 component plot.
The refreshed S4 image was visually inspected. Existing provisional benchmark
raw observations, summaries and metadata were retained; no timings were collected.

Validation used `/workspaces/Quantum-Computing/.venv/bin/python` from this
worktree, without a `PYTHONPATH` override:

| Check | Result |
|---|---|
| `-m pytest -q` | 345 passed in 27.28 seconds |
| `run_project.py symmetry` | N=2 milestone and N=2–5 full/reduced demo passed |
| `scripts/symmetry_benchmark.py --dry-run` | Full N=2–12 and reduced N=2–12,16,24,32,64 fit the configured budget |
| `scripts/symmetry_study.py` | Exact agreement and spatial symmetry-breaking study completed |
| `scripts/symmetry_figures.py` | Recorded-data validation and all figures completed |
| `-m compileall -q src scripts tests run_project.py` | Passed |
| `-m pip check` | No broken requirements found |
| `git diff --check` | Passed |

The suite includes 254 milestone, 34 launcher, 52 symmetry and 5 recovery tests.
Maximum absolute fidelity discrepancy over N=2–12 is 2.504663143554353e-13;
maximum phase-aligned amplitude error is 7.668228917234767e-15. The independent
N=2 milestone dark population is 1.251854463851703e-32. The physical chain still
shows substantial symmetry leakage for N=3–6, as expected for nonuniform pairs.

Matplotlib used a temporary cache because its default configuration directory
was not writable. Pip disabled its unwritable cache. Neither affected results.

The implementation and correctness validation are complete. Final scientific
timings still require the user's controlled physical computer, using the
commands in the [extension report](symmetry_scaling_extension.md). Existing
Codespaces timings and runtime figures remain explicitly provisional.
