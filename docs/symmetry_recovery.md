# S1 recovery audit — 2026-10-10

This is the historical recovery audit. The subsequent launcher integration,
updated diagnostics and current verification are recorded in
[the implementation resumption record](symmetry_resume.md).

The reproduced failure was test collection from the wrong worktree. No failing
mathematical assumption was found. The existing symmetry Hamiltonians and
solvers required no correction. Numerical equivalence was independently checked
before changing the test configuration.

## State and original failure

The initial working directory, `/workspaces/Quantum-Computing`, is on
`extension/quantum-zeno`. S1 lives in its existing nested worktree,
`/workspaces/Quantum-Computing/symmetry-scaling`, on
`extension/symmetry-scaling`, with HEAD `cafa485`. All recovery changes are in
that S1 worktree. No branch switch, reset, commit, tag update, or benchmark was run.

`milestone-v1.0` is unchanged: annotated tag object
`571acfe17c3396ac61b6450e5f016e6987928b01`, targeting
`cafa4850f313174685fa491cc3099ab308fe6c58`.
The original 83 tracked/untracked S1 files were copied to
`/tmp/s1-recovery/before` with a SHA-256 manifest before edits.

The shell's default Python first failed with `No module named pytest`. Activating
the existing parent `.venv` exposed the actual collection failure:

```text
___________________ ERROR collecting tests/test_symmetry.py ____________________
ImportError while importing test module '/workspaces/Quantum-Computing/symmetry-scaling/tests/test_symmetry.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
/home/codespace/.python/current/lib/python3.14/importlib/__init__.py:88: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
tests/test_symmetry.py:14: in <module>
    from quantum_project.symmetry import (
E   ModuleNotFoundError: No module named 'quantum_project.symmetry'
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
=============================== 1 error in 0.91s ===============================
```

There were **no failed assertions or assertion values**: collection stopped
after finding 254 milestone tests and one S1 import error. The isolated S1
module reproduced the same error. Its expanded traceback and the complete
suite output are preserved in
[`original-isolated.log`](../results/data/symmetry_recovery/original-isolated.log)
and [`original-suite.log`](../results/data/symmetry_recovery/original-suite.log).

The shared virtual environment's editable-install `.pth` points to
`/workspaces/Quantum-Computing/src`. Direct inspection confirmed that
`quantum_project.__file__` resolved there. That worktree lacks the S1 module.
Temporarily setting `PYTHONPATH` to S1's `src` produced **298 passing tests**,
including **44 passing S1 tests**, before any edits. No milestone regression or
physics failure could be reproduced. An earlier failure with a different
traceback cannot be diagnosed from the available evidence.

## Minimal repair

The sole configuration change is `pythonpath = ["src"]` under
`[tool.pytest.ini_options]` in `pyproject.toml`. Pytest now selects source from
its own worktree without changing the shared editable installation.

`tests/test_symmetry_recovery.py` adds five independent checks: source location,
N=2 matrices/dynamics, N=2–4 dynamics, tensor site ordering, and global phase.
No existing source, test, tolerance, README, figure, or recorded scientific data
was edited. This report and `results/data/symmetry_recovery/` are new artifacts.
The existing extension report's 298-test count describes its earlier run; the
recovery suite contains 303 tests.

## Basis and independent matrix checks

With s = 1/sqrt(2), the **columns** of U are the new basis vectors, expressed
in the old ordered basis |00>, |0r>, |r0>, |rr>:

```text
U = [[1, 0,  0, 0],       columns: |00>, |Psi+>, |Psi->, |rr>
     [0, s,  s, 0],
     [0, s, -s, 0],
     [0, 0,  0, 1]]
```

Thus psi_old = U psi_new and H_new = U† H_old U. Each matrix entry was also
computed directly as `<U[:,i] | H_old | U[:,j]>`, and H_old U = U H_new was
checked. Because this particular U equals U†, the inverse convention would
accidentally give the same matrix. A separate check rephases the columns by
exp(i*[0.2,-0.4,0.7,1.1]); direct basis-vector matrix elements then agree with
Q† H Q and distinguish it from Q H Q†.

Writing c = Omega exp(-i phi)/sqrt(2), the independently specified matrices are:

```text
H_new/hbar = [[Delta, c, 0, 0],    H_sym/hbar = [[Delta, c, 0],
              [c*,    0, 0, c],                  [c*,    0, c],
              [0,     0, 0, 0],                  [0,    c*, V-Delta]]
              [0,    c*, 0, V-Delta]]
```

The 3x3 reference is constructed directly from its entries before comparing
it with the verified symmetric block. All entries were tested over the full
24-case Cartesian product Omega/scale={1,-2}, Delta/scale={0,0.3},
phi={0,0.37}, V/scale={0,2,-3}, with scale=2*pi*10^6 rad/s.
The transformed antisymmetric row and column vanish to roundoff. The
independent tensor Hamiltonian at N=2 equals the milestone 4x4 matrix exactly,
with **no global energy shift**.

The milestone one-site detuning term is Delta*(1/2-n_i). Summing gives
Delta*(N/2-k); adding the unordered excited pairs gives V*k*(k-1)/2.
The k-to-k+1 matrix element is
Omega/2 * exp(i phi) * sqrt((N-k)*(k+1)); the reverse is its conjugate.
At N=2 these give diagonals Delta,0,V-Delta and exactly the required
collective coupling. No Hz/rad/s conversion, sign, or factor correction is needed.

For N=2–4, independently enumerated bit strings provide normalised Dicke
columns of weight 1/sqrt(binomial(N,k)). Independent bit flips construct a
full Hamiltonian without Kronecker products. Checks cover dimensions 2^N and
N+1, pair counts, H_full B = B H_D, and explicit local occupation operators
whose leftmost site is the most significant bit. The latter check detects
site reversal even when uniform drives would conceal it.

## Independent evolution and phase handling

At N=2, the milestone eigensolver is compared to a dense matrix exponential of
the independent 3x3 matrix over 51 times in [0,0.5 microseconds]. Reduced
states are explicitly expanded as [a0,a1/sqrt(2),a1/sqrt(2),a2]. Fidelity is
computed in the full space without normalisation or clipping, with norm drift
checked separately. For N=2–4, a dense exponential of the independent
bit-flip matrix is compared to the reduced sparse propagator over 31 times in
[0,0.25 microseconds], for the same 24 parameter sets per N.

The recovery audit uses fidelity and probabilities as primary checks. Its
amplitude diagnostic aligns each expanded state by
exp(-i arg(<psi_full|psi_expanded>)). A deliberately imposed phase of 1.234
radians demonstrates that physical errors stay near zero while the existing
`compare_states` raw amplitude diagnostic exceeds 0.5. That existing diagnostic
is retained unchanged; it is not treated as a phase-invariant physical error.

| Maximum observed error | Value |
|---|---:|
| U†U - I, elementwise | 2.220446049250313e-16 |
| N=2 transformed matrix entry (rad/s) | 3.60325190760751e-10 |
| N=2 tensor vs milestone matrix (rad/s) | 0 |
| N=2 abs(F-1), including roundoff above one | 2.6645352591003757e-15 |
| N=2 computational probability | 2.4424906541753444e-15 |
| N=2 phase-aligned amplitude | 2.1138016631127878e-15 |
| N=2 expanded norm minus one | 5.551115123125783e-16 |
| N=2 antisymmetric population from ground state | 2.9890432736889896e-31 |
| N=2–4 abs(F-1) | 8.881784197001252e-15 |
| N=2–4 computational probability | 5.773159728050814e-15 |
| N=2–4 phase-aligned amplitude | 9.929671126925273e-15 |
| N=2–4 total excitation observable | 2.3092638912203256e-14 |
| N=2–4 Dicke projector population | 6.8833827526759706e-15 |

The largest dimensionless discrepancy among the reported state/observable
checks is the excitation-observable error, 2.31e-14. Matrix errors carry
angular-frequency units and are reported separately. Machine-readable metrics
are in [`metrics.json`](../results/data/symmetry_recovery/metrics.json).

Exact comparisons use uniform all-to-all interactions only. The pre-existing
nonuniform-chain test intentionally checks leakage and disagreement. Correctness
tests do not collect timings: benchmark subprocesses either fail the CI guard
before allocation or use `--dry-run`; summary tests use synthetic or recorded
rows. The new audit never exceeds N=4 (16 full amplitudes); unchanged existing
tests have bounded allocations no larger than N=7 (128 amplitudes).

## Final verification and reproducibility

From the S1 worktree, activate `../.venv/bin/activate`, then run
`python -m pytest -rP`. No PYTHONPATH override is needed for pytest after the
repair. Standalone scripts still require the documented `PYTHONPATH=src` when
using this shared virtual environment.

- Final complete S1-worktree suite: **303 passed** (254 unchanged milestone,
  44 unchanged S1, 5 new audit tests).
- Each symmetry module separately: **44 passed** and **5 passed**.
- Existing launcher tests, run in the parent worktree containing the launcher:
  **34 passed**, unchanged. Total distinct tests verified: **337**.
- `python -m compileall -q src scripts tests run_project.py`: extension sources
  compile, with `Can't list 'run_project.py'` because the milestone-based S1
  worktree predates the launcher. The same complete command in the parent
  worktree succeeds silently and covers the real launcher.
- `python -m pip check`: **No broken requirements found**.
- `git diff --check`: clean. No commits made.

Warnings were not suppressed. Matplotlib logged that
`/home/codespace/.config/matplotlib` is not writable and that it created a
temporary cache in `/tmp`. Pip warned that `/home/codespace/.cache/pip` is not
writable and disabled its cache. These environment warnings do not affect
the numerical validation. Full final suite and launcher outputs are preserved
alongside the original tracebacks in `results/data/symmetry_recovery/`.

There is no remaining S1 recovery blocker. The previously documented final
benchmark on the user's physical computer remains a separate pending activity;
this recovery makes no new performance claims.
