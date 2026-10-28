# Reproduce the Phase 4 experiments

From the repository root with the project's Python dependencies installed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m quantum_project.ancilla.run --resume
```

After `python -m pip install -e '.[test]'`, the `PYTHONPATH=src` prefix is
optional. The dedicated module leaves the existing Chapter 5 menu unchanged.
This command runs the bounded canonical Gate 1–5 validation experiments, then
creates the Gate 6 figures and comparison report. It executes in the foreground
and stops on the first failure. It does not implement future development gates
or launch background work.

By default outputs go under `results/ancilla/runs/`. `latest.json` identifies
the five experiment directories and the figure directory for the last completed
invocation. Each experiment has `validation.json`, `residuals.csv`, execution
logs and a checksum manifest; Gates 1–4 also retain complex trajectories in NPZ.
Gate 5 stores raw population trajectories and adaptive residuals in CSV, with a
recipe for regenerating its large complex arrays. The figures directory contains:

- `accuracy.pdf/png`: exact reduction checks and large-N time convergence.
- `protection.pdf/png`: selective success/fidelity and coherent gap suppression.
- `large_n.pdf/png`: bulk and distinguished-site excitation trajectories.
- `adaptive.pdf/png`: measured model error and a static residual upper bound.
- `errors.csv`, `comparison.json`, `figure_data.npz`, `report.md`, `manifest.json`.

## Safe resume behavior

A cache key includes the experiment identity and parameter recipe, numerical
source and validation-script SHA256 hashes, Python/library versions, platform
and thread policy. Numerical parameters are fixed in the canonical recipes;
changing them in a recipe changes its source hash. Figure keys additionally
include plotting code and the exact experiment manifests. Documentation edits
alone do not invalidate numerical results. Plotting edits regenerate figures
without rerunning scientific experiments.

`--resume` skips a directory only when its input specification, required file
inventory, passed gate report and every artifact checksum verify. It never
trusts a directory name alone. A damaged output is retained with `.invalid-...`
and recomputed. Changed inputs create a different directory. A failed experiment
retains its `.pending-...` folder and log; only successfully validated output is
atomically renamed to a completed cache directory. Without `--resume`, existing
outputs cause an error. There is no implicit overwrite or deletion of old runs.

A per-output-root advisory lock rejects overlapping invocations on Unix.
The operating system releases it even after a hard process kill. The persistent
`.run.lock` file contains the PID during a run and is emptied on normal exit;
its mere existence never blocks resume. No stale-lock deletion is required.
The default foreground subprocess timeout is 300 seconds per gate; the
underlying propagation and reference-allocation guards remain in force.

For a bounded subset without figures:

```bash
PYTHONPATH=src python -m quantum_project.ancilla.run --resume --gates 5 --no-figures --output-dir /tmp/ancilla-smoke
```

For a fresh complete reproduction on the same computer:

```bash
PYTHONPATH=src python -m quantum_project.ancilla.run --resume --output-dir /tmp/ancilla-reproduction
```

Repeat the command to verify that all five experiments and the figures are
skipped. Compare `latest.json` and manifest hashes to confirm that validated
files were not modified. Runtime records are provisional; they are not a
controlled benchmark or evidence of a speedup.

## Scientific interpretation

The comparison report distinguishes exact symmetry, exact hybrid structured
asymmetry, postselected protected dynamics and approximate adaptive dynamics.
Coherent gap protection is labelled separately. The large-N distinguished atom
is a physical model degree of freedom, not a syndrome ancilla. The ideal
projector oracle has unresolved hardware compilation costs.

Numerical L2 residuals, physical leakage, success probabilities and model
approximation errors are different quantities. The CSV table labels their
meaning. No sampling or experimental uncertainty is inferred from deterministic
evolution. The finest midpoint error uses asymptotic Richardson estimation;
the Taylor action truncation bound excludes roundoff. The static adaptive bound
is a Duhamel bound for N<=8 diagnostics, not an efficient general-disorder method.

Source Git HEAD and dirty-state records describe the local source checkout
used for computation; source hashes identify the exact numerical code even
when a working change is validated before committing. In the hosted development
environment, local source history was reconstructed through the GitHub API.
Its local snapshot SHA is not the same as the upstream branch history. The
upstream checkpoints are recorded in `STATUS.md`; do not treat a local snapshot
as a remote ancestor or substitute a later commit SHA into historical data.
