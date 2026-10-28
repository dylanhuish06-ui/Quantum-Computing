# Phase 4 — Persistent Progress Ledger

Canonical plan: [MASTER_PLAN.md](MASTER_PLAN.md). Branch: `extension/ancilla-hybrid-phase4`.

## Checkpoint (2026-10-23 UTC)

Verified starting HEAD: `c816a85a01602e2595cdd780955c4fec0ee0498f`.
The commit containing this ledger records the Gate 0 audit. Resolve its SHA with
`git log -1 -- docs/ancilla_phase4/STATUS.md`; a commit cannot contain its own SHA.
Worktree: `/workspaces/Quantum-Computing/ancilla-phase4`; clean at audit start.
Unrelated edits in the original workspace and `ancilla-hybrid` are preserved.
No main changes or merges.

| Gate | State | Evidence |
| --- | --- | --- |
| 0 — Baseline and audit | passed | 345 inherited tests passed; source/conventions reviewed |
| 1 — Operators/reference | passed | 35 reference tests; N=2–6 independent numerical report |
| 2 — Ancilla/Zeno | passed | 6 protocol tests; 15 N/check-count numerical comparisons |
| 3 — Gap/dynamical | passed | 15 tests; 20 gap comparisons and finite-cycle scan |
| 4 — Exact hybrid | passed | 41 hybrid tests; 27 N/m numerical comparisons |
| 5 — Large-N/adaptive | passed | N=64/256/1024 sparse sweep refinements; bounded static adaptive reference and failure example |
| 6 — Figures/CLI | passed | Four PDF/PNG figure sets, hashed resumable CLI, same-session comparisons; 447 tests and end-to-end resume passed |
| 7 — Review | not_started | — |

## Baseline and environment

Inherited symmetry implementation: `1eb445c`, milestone ancestor `cafa485`.
Read `docs/symmetry_scaling_extension.md`, `symmetry.py`, `hamiltonian.py`,
`units.py`, launcher, dependency configuration and tests. Physics: Hbar=H/hbar
in rad/s, seconds, raising phase exp(+i phi), detuning delta*(1/2-n).
Retain identity offsets for complex amplitudes. Leftmost atom is the most
significant bit. Existing 256 MiB guards are workspace heuristics, not peak RSS
bounds. Existing provisional benchmarks are not rerun or overwritten.

Python 3.14.2; NumPy 2.5.3; SciPy 1.18.1; Matplotlib 3.11.2; pytest 9.1.1.
No lint/type checker is configured. Use this worktree's source with the shared
interpreter:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mpl-phase4 PYTHONPATH=src /workspaces/Quantum-Computing/.venv/bin/python -m pytest -q
```

Baseline used the same command without MPLCONFIGDIR/PYTHONPATH (pytest selects
src through pyproject.toml): **345 passed in 15.85 s**. No test failures.
The nonwritable Matplotlib cache warning is addressed by the temporary cache
path above. Codespaces CPU times are provisional.

## Artifacts, residuals and limitations

Gate 1: `results/ancilla/gate1/{validation.json,residuals.csv,trajectories.npz}`.
Reference tests: 19 passed (propagation, milestone, spin/disorder selections),
16 passed (`test_phase4_reference.py`). `validate_gates.py --gate 1` passed.
Maximum ODE refinement L2 error 1.446e-13, dense/sparse exponential difference
4.862e-15, constant-H error 7.557e-15, norm drift 8.904e-14. Midpoint converges
at order 2.00004 or closer. Frequencies h=Hbar/Omega_ref, tau=Omega_ref*t;
Omega_ref=2*pi*1e6 rad/s. Raw source hashes and environmental provenance in JSON.
 Reviewed prototypes
on `extension/ancilla-hybrid-dynamics` (`c8008f2`, `a118bbd`) are a reuse source,
not evidence of passed Phase 4 gates. Their uncommitted scripts/reports/resources
edits are outside this checkpoint. Exact hybrid dynamics require invariant bulk
symmetry and symmetric initial bulk. Ancilla postselection changes the original
dynamics. Nonlocal oracle compilation costs remain unresolved. No high-N
full-space allocations allowed.

## Immediate next action

Begin Gate 7: final scientific/code audit, artifact/link checks and hand-off review. Do not merge to main without explicit approval.

## Session log

- Packaging: plan and ledger only; no scientific results.
- 2026-10-23: fetched the research branch into an isolated worktree; passed Gate 0.

Gate 2 verification: 6 tests passed, 15 deselected in 2.70 s;
`validate_gates.py --gate 2` passed 15 comparisons, 36.52 s (provisional).
Artifacts: `results/ancilla/gate2/{validation.json,residuals.csv,trajectories.npz}`.

Gate 3 verification: 15 tests passed, 6 deselected in 0.95 s;
`validate_gates.py --gate 3` passed 20 gap comparisons plus finite-cycle checks.
Artifacts: `results/ancilla/gate3/{validation.json,residuals.csv,trajectories.npz}`.

Gate 4 verification: 41 tests passed in 25.33 s; `validate_gates.py --gate 4`
passed 27 comparisons in 24.32 s (provisional). Artifacts in `results/ancilla/gate4`.

Full-suite gate 1–4 checkpoint: **429 passed in 37.62 s** with the documented
single-threaded pytest command. No failures. See [VALIDATION_REPORT.md](VALIDATION_REPORT.md)
for models, derivations, measurements and limitations. Hybrid maximum complex-state
L2 error 4.705e-12; compression/invariance maximum element 3.553e-15;
local population error 2.231e-12. N=6 maximum leakage falls from 0.3448744 to
0.00017543 at G/Omega_ref=8. No generic disorder or large-N dynamical claim yet.

Commit whitespace check initially flagged CSV CRLF endings; normalized to LF and
updated artifact hashes without rerunning numerical work. No scientific failure.

## Gate 5 checkpoint — 2026-10-28 UTC

Verified upstream branch HEAD before implementation:
`c6b1cac624ba995966df1249467c7e1367d83b05`. Source obtained through the connected
GitHub API because shell cloning had no private-repository credentials. Local
source snapshot is not an upstream ancestor; a future remote commit should be
directly parented on the verified upstream HEAD. Only Gate 5 implementation,
tests, reports and compact artifacts are prepared for publication.
The commit containing this section is the Gate 5 checkpoint; resolve its SHA
from the branch history rather than treating the older header as current HEAD.

Added independently budgeted sparse Taylor action, integrated as
`method="controlled_midpoint", action_atol=...`; added a small-N static
adaptive bulk grouping with Duhamel residual bound. Retained existing solver
interfaces and defaults. The adaptive diagnostic is N<=8, uses full-space
reference matrices and can require the full state dimension. It is not a
scalable generic-disorder simulator. Large-N sweep atoms are distinguished
physical qubits, not syndrome ancillas.

Numerical command (passed; fresh output directory):
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python scripts/ancilla/validate_gate5.py --output-dir /tmp/ancilla-gate5-final`.
Artifacts copied byte-for-byte to `results/ancilla/gate5/`: validation JSON,
residuals CSV, raw observables CSV and adaptive CSV. See
[GATE5_REPORT.md](GATE5_REPORT.md) for exact parameters and limitations.

N=1024 dimension 2048, observed midpoint order 2.000004; estimated finest
complex-state error 6.5859e-5 (dimensionless L2, asymptotic estimate);
independent exponential-action difference 1.9456e-14; norm drift 2.8422e-14.
Time step, action budget and output sampling independently refined. The N=5
generic-disorder example meets the 1e-4 error budget only at full dimension 32.
No statistical uncertainty or hardware-performance claim.

Full-suite validation command:
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q`.
Final outcome: **439 passed in 32.68 s** (10 new tests; no failures).
Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, pytest 9.2.0. Timing is
provisional. Artifact and implementation source SHA256 checks passed, as did
the tracked-file whitespace check. Gates 6–7 remain not_started.

Publication: user explicitly authorized the push on 2026-10-28. Gate 5 is
published as `c4352023c5bcb5c1dd51550f458802c0cd1047e3`, directly parented on
`c6b1cac624ba995966df1249467c7e1367d83b05`. The branch update succeeded.
The earlier automatic review block is resolved.

## Gate 6 checkpoint — 2026-10-28 UTC

Upstream starting HEAD: `c4352023c5bcb5c1dd51550f458802c0cd1047e3`.
The commit containing this checkpoint is the Gate 6 hand-off; resolve its remote
SHA from branch history. Local reconstructed source history differs from remote
history, as described in [RUNNING_EXPERIMENTS.md](RUNNING_EXPERIMENTS.md).
Unrelated retrieved scripts/data and cache files were preserved locally. Gate 6
implementation, tests and documentation are prepared for publication. The
selected numerical artifacts and figures remain local pending separate upload
authorization: automatic approval review rejected the Gate 2 trajectory NPZ
upload as a distinct numerical payload. No complete Gate 6 data tree was
published. The source-only checkpoint can reproduce all results through the CLI.

Added `python -m quantum_project.ancilla.run --resume` (use `PYTHONPATH=src`
from an uninstalled checkout). Each experiment is a bounded foreground task;
matching parameter/recipe, numerical-source and environment hashes plus artifact
checksums are required before skipping. Damaged results are preserved before
recomputation, pending results are never treated as passed, and a Unix advisory
lock prevents overlapping runs without stale-lock files blocking recovery.

Gate 1–5 experiments recomputed in the current Python 3.12.14 / NumPy 2.3.5 /
SciPy 1.17.0 environment with single-threaded BLAS. The final command was:
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m quantum_project.ancilla.run --resume`.
Selected completed directories are recorded in `results/ancilla/runs/latest.json`.
New Gate 1–4 NPZ complex trajectories and Gate 5 raw observables/adaptive CSVs
are retained with their validation JSON, residuals, execution logs and manifests.

Four publication figure sets were generated in PDF and PNG and visually checked:
`results/ancilla/runs/figures-ac66ceca6fdce52e/{accuracy,protection,large_n,adaptive}`.
The same directory contains `report.md`, `errors.csv`, `comparison.json`,
`figure_data.npz` and checksum manifest. Exact symmetry, exact hybrid asymmetry,
postselection and approximate dynamics have separate interpretations; coherent
gap protection is also labelled. No timing speedup or statistical-error claim.

End-to-end second invocation: all five experiments and figures skipped after
checksum validation; every artifact SHA256 and modification time remained
unchanged. Cache-corruption, changed-input, interrupted-task and lock behavior
are additionally unit tested. Final suite command:
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q`.
Outcome: **447 passed in 33.24 s** (eight Gate 6 tests, no failures).
Local verification logs and machine-readable summary: `results/ancilla/gate6/`.
The artifact paths above describe the verified local checkpoint; they are not
files included in the source-only remote commit.

No new physical claim beyond Gate 5: structured N=1024 evolution is exact only
under bulk invariance; the strict generic-disorder example reaches full space.
Ideal oracle compilation costs remain unresolved. Stage 7 is not_started and
the programme is not yet ready_for_review. Immediate next action: Gate 7 audit.

