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
| 1 — Operators/reference | not_started | Committed prototypes located; not validated on this branch |
| 2 — Ancilla/Zeno | not_started | — |
| 3 — Gap/dynamical | not_started | — |
| 4 — Exact hybrid | not_started | — |
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

No Phase 4 numerical artifacts or measured residuals yet. Reviewed prototypes
on `extension/ancilla-hybrid-dynamics` (`c8008f2`, `a118bbd`) are a reuse source,
not evidence of passed Phase 4 gates. Their uncommitted scripts/reports/resources
edits are outside this checkpoint. Exact hybrid dynamics require invariant bulk
symmetry and symmetric initial bulk. Ancilla postselection changes the original
dynamics. Nonlocal oracle compilation costs remain unresolved. No high-N
full-space allocations allowed.

## Immediate next action

Import reviewed committed operator/reference prototypes, extend missing
validation/safeguards, and independently verify Gate 1 before advancing.

## Session log

- Packaging: plan and ledger only; no scientific results.
- 2026-10-23: fetched the research branch into an isolated worktree; passed Gate 0.
