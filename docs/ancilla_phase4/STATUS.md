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
| 5 — Large-N/adaptive | not_started | — |
| 6 — Figures/CLI | not_started | — |
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

Begin Gate 5: independently controlled exponential-action precision, sparse N=64/256/1024 sweep convergence, and a small-N residual-monitored asymmetric approximation.

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
