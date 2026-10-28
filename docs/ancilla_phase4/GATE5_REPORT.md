# Gate 5 — Sparse sweeps and controlled asymmetric approximation

Validated 2026-10-28 against upstream `c6b1cac624ba995966df1249467c7e1367d83b05`
on `extension/ancilla-hybrid-phase4`. This checkpoint extends the research
branch, leaving the Chapter 5 model and poster scripts unchanged.

## Independent numerical controls

`controlled_midpoint` adds sparse Taylor exponential action with a separately
specified `action_atol`. Unlike SciPy's `expm_multiply`, this backend exposes
an explicit truncation budget. For Hermitian H, after shifting its real trace,
the integral remainder is bounded by

\[
\|v\|\frac{(|h|\|H\|_2)^{p+1}}{(p+1)!}.
\]

The induced matrix 1-norm bounds the Hermitian spectral norm. Subdivision
keeps the expansion argument at most one. Local budgets sum to the supplied
action budget. The bound excludes floating-point roundoff and midpoint
time-ordering error. No state renormalisation hides either effect. Resource
guards cap matrix-vector products, trajectory storage and callback runtime;
runtime checks cannot interrupt a running numerical kernel.

This is a Taylor method, not a Krylov method. Tests compare it with dense
exponentials at positive, negative and zero times and with SciPy action on a
time-dependent model. Time integration, action budget and output sampling are
refined independently.

## Exact structured large-N model

Use tau=Omega_ref*t and h=H/(hbar*Omega_ref), Omega_ref=2*pi*1 MHz.
The bulk has N-1 atoms and one distinguished physical atom, initially all
in their ground states. The complete time-dependent profile is defined in
`scripts/ancilla/validate_gate5.py`; it includes a linear detuning sweep,
changing drive and phase, and nonzero pair couplings. Bulk permutation
symmetry holds at every time. The distinguished atom is not a measurement
ancilla. These trajectories solve the structured model, not generic disorder.

Sweep duration is tau=0.4, or about 63.66 ns. Time steps are 0.01, 0.005 and
0.0025; action budgets are 1e-11 and 1e-13, with 21/41 output samples. Every
refinement starts from the same initial state. Phase-sensitive L2 errors retain
the identity energy offset and use no global-phase alignment.

| N | Reduced dimension | Observed time order | Estimated finest time error | Norm drift |
| --- | --- | --- | --- | --- |
| 64 | 128 | 2.000004 | 4.3444e-6 | 1.08e-14 |
| 256 | 512 | 2.000004 | 1.6653e-5 | 1.53e-14 |
| 1024 | 2048 | 2.000004 | 6.5859e-5 | 2.84e-14 |

The time errors are asymptotic Richardson estimates, not rigorous bounds.
The N=1024 independent action comparison is 1.95e-14 in trajectory L2;
action-budget refinement is 1.32e-15 and output-grid refinement 1.70e-15.
The finest result meets the declared 1e-4 complex-state threshold; this short
sweep is not evidence of accuracy at arbitrarily long times or other profiles.
No speedup claim is made. These are deterministic numerical errors, not
statistical uncertainties or experimental error bars.

## Generic asymmetry: residual monitoring and failure example

`adaptive_static_grouping` is an explicitly small-N full-space reference.
It forms h_m=B†HB, measures R=HB-Bh_m, and promotes bulk atoms to the explicit
register until T*||R||_2 meets a requested error budget. For static Hermitian H
and initial states in the image of B, Duhamel's formula gives

\[
\|e^{-iHT}Bv-B e^{-ih_mT}v\|\le T\|R\|_2\|v\|.
\]

The bound excludes floating-point errors. It can be conservative and need
not decrease at every promotion. The implementation rejects N>8 rather than
presenting its full-space residual calculation as a scalable algorithm.

The N=5 example uses drive 1, zero phases/pair coupling and local detunings
[-0.8, 0.7, -0.4, 0.3, 1.1], T=0.7, and a fixed all-ground physical initial
state at every grouping. With no distinguished atom the actual state error is
0.1810, below the 0.5974 bound. The 1e-4 budget is met only at m=4, dimension
32, where the representation covers full space. This is a reproducible
failure of useful compression for the selected strict budget. It does not
show efficient exact simulation of arbitrary asymmetric qubits.

## Reproduce and inspect

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m pytest -q
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python scripts/ancilla/validate_gate5.py --output-dir /tmp/ancilla-gate5-reproduction
```

Use a new or empty output directory. Committed artifacts are
`results/ancilla/gate5/{validation.json,residuals.csv,observables.csv,adaptive.csv}`.
JSON records parameter choices, source/artifact SHA256 hashes and runtime
versions. Compact raw observable trajectories are stored; complex state arrays
are reproducible through the script. No high-N full-space array is constructed.
The retained data permit comparison of excitation fractions and distinguished
populations without committing large arrays.

Gate 6 remains: publication figures, broader cross-protocol comparison, and a
stable resumable experiment CLI with parameter-hash validation. Gate 7 review
has not begun. No hardware compilation or generic-disorder scaling result is
claimed by this checkpoint.
