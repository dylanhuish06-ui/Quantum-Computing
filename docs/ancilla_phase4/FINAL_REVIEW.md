# Phase 4 final review and hand-off

Status: **ready_for_review**. Gates 0–7 have passed their documented checks.
The research branch is `extension/ancilla-hybrid-phase4`. Gate 6 source was
published as `7c24f053910d5a2beb83bd81d85f70e0014873f4`; its complete numerical
checkpoint was published as `136677adb3a0f076e568cab3b322bbb434dc2e15`.
The commit containing this report is the Gate 7 hand-off. No merge is included.

## What was verified

The final full suite passed **447 tests in 30.49 seconds**. The command was
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q`.
No lint or type checker is configured; Python compilation of the ancilla
package and scripts passed. Timings are provisional environment measurements.

The documented CLI smoke run verified and skipped all five numerical gates
and the figure step. A separate read-only audit verified the scientific source
and environment signatures, required file inventories, SHA256 checksums,
figure-to-experiment provenance, 310 finite raw trajectory arrays, four PDF/PNG
figure pairs and all relative Markdown links in the Phase 4 documentation.
The 44 published Gate 6 artifact Git blob hashes match the local files.
Evidence is saved under `results/ancilla/gate7/`.

Reproduce from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q
PYTHONPATH=src python -m quantum_project.ancilla.run --resume
PYTHONPATH=src python scripts/ancilla/audit_handoff.py
```

Use `--output-dir /tmp/ancilla-independent-reproduction` for a fresh numerical
run. Dependency/platform changes intentionally produce different cache keys.
`results/ancilla/runs/latest.json` selects the validated reports and figures.
See [RUNNING_EXPERIMENTS.md](RUNNING_EXPERIMENTS.md) for recovery and resource
controls, and [GATE5_REPORT.md](GATE5_REPORT.md) for the action-error derivation.

## Mathematical and physical assessment

Hamiltonians retain the detuning identity offset and the raising phase
exp(+i phi), using frequencies H/hbar in rad/s and time in seconds, or the
explicit equivalent h=H/(hbar Omega_ref), tau=Omega_ref*t convention. Complex
amplitude comparisons retain physical phase and apply no phase alignment.

The exact hybrid representation has dimension (N-m+1)*2**m. Its validity
requires permutation-invariant bulk coefficients at every time and a symmetric
initial bulk state. Compression and invariance checks against full space are
separate tests: B†HB=h alone would not prove that discarded sectors decouple.
The large-N sweep meets these structural assumptions. N=1024 uses dimension
2048 with one distinguished physical atom, not a syndrome ancilla. Its finest
time error estimate is 6.586e-5 in complex-state L2; norm drift is 2.842e-14.
These results cover the documented short sweep, not arbitrary long-time dynamics.

The ideal syndrome unitary P tensor I + Q tensor X is unitary because P and Q
are complementary orthogonal projectors. Conditional target fidelity, probability
of all successful syndromes and their product are reported separately.
Postselection changes asymmetric physical dynamics. High conditional fidelity
is not evidence of deterministic error correction or simulation of the original
Hamiltonian. The spin penalty G[Jmax(Jmax+1)-J²] vanishes in the symmetric sector;
the next sector is separated by NG. Its all-to-all physical implementation is
an assumption, not a demonstrated hardware resource saving.

Taylor exponential action and time integration have independent accuracy
controls. The action truncation bound uses the unitary integral remainder and
a Hermitian spectral-norm bound; it excludes rounding error. The midpoint
Richardson estimate is asymptotic, not a rigorous bound. Norm conservation alone
cannot certify phase or trajectory accuracy, hence the independent references
and refinement checks.

For the static approximate model, Duhamel's formula bounds the state error by
T*||HB-Bh|| for normalized initial states in B. This diagnostic explicitly uses
small full-space matrices (N<=8). The strict N=5 generic-disorder example reaches
dimension 32 before meeting its budget; useful compression fails in that case.
The implementation does not supply a scalable adaptive solver for general
time-dependent disorder. The residual need not decrease at every promotion.

## Uncertainty, scaling and feasibility

The error table distinguishes deterministic numerical error, physical leakage,
postselection success and approximation error. No sampled noise model or
experimental uncertainty was measured, so statistical error bars are omitted.
There is no cross-machine runtime benchmark or supported speedup claim.

Explicit projectors, syndrome oracles and adaptive residuals are small-system
references. The large-N solver avoids full-space expansion but scales
exponentially in the number of distinguished atoms. State/storage guards are
heuristics, not certified peak-memory bounds; subprocess timeouts supplement
checks between numerical kernels. The runner uses Unix advisory locking.

Compilation of the nonlocal controlled projector remains unresolved. Neither
constant-depth local implementation nor hardware-level fault tolerance follows
from these simulations. Permutation-control costs assume ideal SWAP/CNOT
equivalents; realistic pulse imperfections need separate experiments.

Recommended follow-up experiments are longer structured sweeps with tighter
phase-sensitive convergence, explicit noisy readout/reset models, finite-pulse
control comparisons and a genuinely scalable time-dependent residual estimator.
These are research extensions rather than missing validation of the present
documented models. Human scientific review is the next action.
