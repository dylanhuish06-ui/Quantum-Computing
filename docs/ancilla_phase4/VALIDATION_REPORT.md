# Phase 4 mathematical and numerical validation

## Scope and provenance

This report covers the first implementation gates on
`extension/ancilla-hybrid-phase4`. Code in `src/quantum_project/ancilla` and its
initial tests was recovered from committed research prototypes `c8008f2` and
`a118bbd`, reviewed, and independently checked here. No branch merge or copying
of uncommitted predecessor reports/results was performed. The baseline ancestor
is `1eb445c`; all 345 inherited tests passed before adding Phase 4 code.

The milestone and its generated figures are untouched. A distinguished physical
atom in a hybrid basis belongs to the simulated system; it is not the temporary
measurement ancilla. No arbitrary disordered N-qubit simulation is claimed to
have polynomial exact cost.

## Conventions and assumptions

The physical generator is Hbar=H/hbar in rad/s. The computational basis is
binary, leftmost site most significant, with 0=ground and 1=Rydberg. Raising
has coefficient Omega_i exp(+i phi_i)/2. Detuning is Delta_i(1/2-n_i), including
the constant offset. N=2 matrices are checked against the original milestone.

Reproducible studies use h=Hbar/Omega_ref and tau=Omega_ref*t with
Omega_ref=2*pi*1e6 rad/s. A reported dimensionless interval 1.6 is therefore
0.254648 microseconds. Dimensionless G corresponds to physical G*Omega_ref;
the first spin gap is N*G*Omega_ref rad/s. Integration `dt` is in tau in these
studies; the public solver accepts seconds when supplied physical Hbar. A test
checks both unit conventions against the same complex states.

The default disordered profile for gates 1–3 is

- Omega_i/Omega_ref = 1+0.2 sin(2 tau),
- Delta_i/Omega_ref = -0.4+0.3 tau+(0.65+0.1 sin(tau))*linspace(-1,1,N)_i,
- phi_i = 0.23+0.3 tau,
- V_ij/Omega_ref = 0.17*(1+0.1 cos(tau)), for i != j.

All runs start in the ground state unless stated otherwise. This is a smooth,
noncommuting time-dependent model with mean-zero detuning inhomogeneity. It is
an illustrative mathematical experiment, not a species-specific device model.

## Operators and independent propagation (Gate 1)

The full generator is assembled from bit-string occupations and single-site
flips. Separate tests build tensor products with NumPy, compare every entry,
and check exact dense exponentials. The original N=2 Hamiltonian supplies a
third sign/phase convention check. B contains normalized Dicke columns;
P=BB† and Q=I-P. Both compression B†HB and invariance HB=Bh are tested for
symmetric models; merely satisfying compression does not imply invariance.

DOP853 integrates the complex Schrödinger ODE; the separate midpoint method
applies exp(-i h(t+dt/2) dt) using sparse exponential action. Dense SciPy expm
is a third, small-system backend. `expm_multiply` uses a scaled truncated
Taylor action, **not a tunable Lanczos/Krylov solver**. It exposes no independent
Krylov tolerance. Midpoint refinement measures time-ordering error, while dense
versus sparse action compares exponential kernels at the same time steps.
DOP853 tolerances and maximum steps are independently refined. No unitary
trajectory is renormalized and no phase alignment hides offset errors.

DOP853 output times are interpolation points. For discontinuities, use separate
calls on each smooth interval; merely adding a discontinuity to the output grid
is insufficient. A state's supplied time origin is `times[0]`.

At N=2–6, Gate 1 measured maximum DOP853 refinement error 1.446e-13 (complex
state L2), dense/sparse midpoint difference 4.862e-15, constant-generator error
7.557e-15 and probability-norm drift 8.904e-14. The midpoint observed order is
approximately 2.00004. The finest midpoint error is 1.296e-5: norm preservation
alone does not certify time-ordering accuracy. Full raw residuals and trajectories
are in `results/ancilla/gate1`.

The total spin implementation constructs J²=sum_alpha J_alpha² with
J_alpha=(1/2)sum_i sigma_alpha,i. Tests check all allowed j(j+1) eigenvalues,
maximal spin on every Dicke column and the penalty spectrum. For mean-zero
static detuning errors epsilon_i, direct reference tests verify
||QD|D_N^k>||² = k(N-k)*mean(epsilon²)/(N-1). From the ground state and uniform
drive Omega, the leading leakage is N*Omega²*mean(epsilon²)*t⁴/16. This last
expression is a short-time series, not a finite-time error bound.

## Ideal ancilla and selective Zeno (Gate 2)

The system-first joint basis realizes U=P tensor I+Q tensor X. Tests prepare
ancilla |0>, compare both output columns with P psi and Q psi, verify Born
weights and unitarity, and compare the partial trace with P rho P+Q rho Q.
Small-N enumeration of all measurement records independently checks the
trace-preserving nonselective density solver.

`selective_evolution` retains the unnormalized all-zero outcome branch, avoiding
loss of survival weights through normalization. Checks occur at checks[1:];
there is no automatic measurement at the initial time. At zero success the
conditional state/fidelity is undefined (NaN), while successful target
probability is zero. Symmetric invariant dynamics have success one.

The target is the unitary evolution under B†h(t)B, expanded by B. It is protected
target dynamics; the bare disordered evolution generally differs. Reported
conditional target fidelity F_c, all-success probability p_s, and successful
target probability p_s F_c are distinct. The latter is the probability of both
the all-success record and a terminal target projection. It is not the target
probability of the full nonselective ensemble including failed records.

Runs at N=2–6 compare 8, 32 and 128 checks at fixed duration tau=1.6, with
separate ODE and midpoint refinement at each fixed schedule. A check assumes
one ideal global oracle, readout and preparation of |0>. On the all-success
branch the measured ancilla is already |0>; failed outcomes are discarded.
Reusing one ancilla requires N_checks-1 reuses between checks. The oracle is
nonlocal: local gate count, synthesis error and depth are **unresolved**, not
one physical two-qubit gate. The optional spin-phase-register calculation
counts ideal controlled evolutions; it also leaves their compilation unresolved.

## Coherent gap and finite-rate symmetrisation (Gate 3)

Add G[Jmax(Jmax+1)I-J²]. It annihilates the symmetric sector and gives a first
positive gap N G (for N>1). This coherent all-to-all penalty has no discarded
outcomes; its target overlap is unconditional. It changes the original physical
Hamiltonian. Compare the same disorder, pulse duration and initial state as
Gate 2, using G/Omega_ref=0, 0.5, 2, 8. A larger gap can suppress leakage in this
experiment, but monotonic improvement under arbitrary time dependence, including
resonant disorder, is not a theorem established by this scan.

Optional permutation control applies each permutation, free evolution, then
its inverse. Full S_N averaging symmetrizes the static Hamiltonian; a cyclic
subgroup need not symmetrize all pair interactions. Finite cycles introduce
ordering errors. A separate N=3 static experiment counts swaps and the ideal
three-CNOT-per-SWAP cost, without a hardware fidelity or optimized depth claim.
Finite pulses with background dynamics on are tested independently at N=2;
shrinking the pulse width increases the required control norm pi/t_p. No finite
bandwidth hardware implementation is inferred.

## Exact hybrid representation (Gate 4)

The first M=N-m atoms form the bulk, followed by m distinguished physical atoms.
The ordered basis is |D_M^k> tensor |a>, dimension (M+1)*2^m. The sparse generator
contains the exact bulk diagonal Delta_b(M/2-k)+V_b*k(k-1)/2; the neighboring
bulk ladder entries Omega_b exp(+i phi_b)*sqrt((M-k)(k+1))/2; the distinguished
register Hamiltonian; and cross energy k*sum_a V_ba*n_a.

Exactness requires common bulk complex drive and detuning, uniform bulk pairs,
and site-independent coupling from bulk to each distinguished atom at every
time, with initial bulk in its symmetric sector. `HybridParameters` enforces
this structure. `validate_bulk_symmetry` rejects violating full parameters with
zero tolerance by default. Its explicitly nonzero `atol` option is an
approximation, for which the caller must assess error. Sampling an arbitrary
callback cannot prove its symmetry at unsampled times.

Validation spans N=2–8, m=0,1,2,3 where m<=N. Random complex initial states cover
nontrivial distinguished configurations. Tests also cover ground and fully
excited states. The report compares compression, invariance, independent full
ODE evolution, all local populations and a pair correlation. Three midpoint
grids check order two, and Richardson extrapolation checks phase-sensitive
convergence. `split_bulk_qubit` changes representation exactly; it cannot recover
past leakage discarded by an approximate evolution and is not by itself an
adaptive simulation algorithm.

## Reproduction and safeguards

From this worktree, with the project dependencies installed:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mpl-phase4 python -m pytest -q
python scripts/ancilla/validate_gates.py --gate 1 --output-dir /tmp/phase4-new-gate1
python scripts/ancilla/validate_gates.py --gate 2 --output-dir /tmp/phase4-new-gate2
python scripts/ancilla/validate_gates.py --gate 3 --output-dir /tmp/phase4-new-gate3
python scripts/ancilla/validate_gates.py --gate 4 --output-dir /tmp/phase4-new-gate4
```

Use a new output directory for each independent reproduction; the script refuses
to overwrite artifacts. The JSON is written only after numerical exit assertions
pass, and records source/artifact hashes, Git provenance, environment, parameters
and raw residuals. CSV and compressed NPZ retain residuals and complex states.
The single-threaded Linux script has a configurable 240-second alarm per gate;
individual propagations have evaluation and cooperative runtime budgets. No
background jobs or high-N full expansions are used. Full references are capped
at N=12, dense syndrome at N=8, density trajectories at dimension 64 and dense
expm at dimension 256. Hybrid construction and stored trajectories have 256 MiB
workspace estimates, which are conservative heuristics rather than measured RSS.
All recorded Codespaces elapsed times are provisional, not speed-up evidence.

Gate 5 (converged large-N dynamics and controlled adaptive approximation), Gate 6
(publication figures and resumable user CLI) and Gate 7 (final hand-off) remain
outstanding until explicitly verified in STATUS.md. Sparse construction tests
at N=1024 are not large-N dynamical results.

## Verified gate exit measurements (2026-10-23)

All 429 tests passed in 37.62 s after gates 1–4. Gate scripts passed sequentially:
5 reference cases, 15 measurement cases, 20 gap cases (plus the permutation
scan), and 27 hybrid N/m cases. Durations are provisional Codespaces observations.

| Check | Maximum measured error |
| --- | ---: |
| Gate 1 DOP853 tolerance/step refinement (state L2) | 1.446e-13 |
| Gate 1 dense vs sparse midpoint (state L2) | 4.862e-15 |
| Gate 2 fixed-schedule ODE refinement (branch L2) | 2.516e-14 |
| Gate 3 ODE refinement (state L2) | 8.142e-12 |
| Gate 3 midpoint vs ODE (state L2, dt=0.0025) | 8.536e-7 |
| Gate 4 compression/invariance (maximum element) | 3.553e-15 |
| Gate 4 full vs hybrid (state L2) | 4.705e-12 |
| Gate 4 local population | 2.231e-12 |
| Gate 4 pair correlation | 1.466e-12 |
| Gate 4 Richardson-extrapolated state L2 | 5.506e-10 |

At N=6, 8/32/128 selective measurements give success probabilities
0.9347156/0.9831395/0.9957521 and conditional target fidelities
0.9998097/0.9999865/0.9999991. Bare target probability is 0.6514078; bare leakage
is 0.3448744. At G/Omega_ref=8, maximum leakage is 0.00017543 and final target
probability 0.9998250. This demonstrates protection for the specified model,
not faithful recovery of its bare disordered dynamics. At N=3 the static
permutation approximation's state error decreases from 0.04260 (one cycle,
14 swaps) to 0.005271 (eight cycles, 112 swaps), consistent with a first-order
finite-cycle error for the chosen sequence.

Each numerical JSON records the exact source hashes when run. A later docstring
clarification about DOP853 discontinuity boundaries changes the propagation file
hash after gates 1–2 without changing executable code; those validated runs were
preserved rather than unnecessarily repeated.
