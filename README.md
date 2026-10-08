# Quantum Computing — PHYS3561

Dylan Huish’s Durham University MPhys Computing Project, based on Chapter 5,
“Quantum dynamics and quantum computing”. Runs 1–6 implement and audit the
static and time-dependent solvers. Run 6B resolves the source interpretation
and reproduces both milestone figures with validated numerical sampling.

The source specifies Ω/(2π) = 1 MHz and an interaction range of 0–100 MHz,
with a blue V = 0 curve and red interacting curves. It does not enumerate the
individual red-curve values. Both figures use the repository's representative
V/(2π) = 0, 0.5, 1, 5, 100 MHz sampling; these are numerical choices, not an
individual list from the booklet. Figure 5.2 additionally specifies the cubic
sweep, Δmax/(2π) = 50 MHz and duration 1 μs, with other parameters as Figure 5.1.
Constant φ = 0 is a phase convention: P00 and PΨ+ from |00> are invariant under
a common constant laser phase in this symmetric model.

Read the [milestone report](docs/milestone_report.md) and
[interview checklist](MILESTONE_CHECKLIST.md) before user review.

## Physics conventions

Every vector and matrix uses the fixed ordered basis
`|00>, |0r>, |r0>, |rr>`. The canonical vectors are `KET_00`, `KET_0R`,
`KET_R0`, and `KET_RR`; `PSI_PLUS` is `(KET_0R + KET_R0) / sqrt(2)`.
These exported arrays are read-only; use `.copy()` when a mutable array is needed.

`rydberg_hamiltonian(omega, delta, phase, V)` implements the supplied Chapter 5
Eq. (5.4), returning **Hbar = H / ħ**, a complex 4 × 4 array:

```text
        [ 2Δ        Ω exp(-iφ)  Ω exp(-iφ)   0          ]
Hbar =  [ Ω exp(iφ) 0          0           Ω exp(-iφ) ] / 2
        [ Ω exp(iφ) 0          0           Ω exp(-iφ) ]
        [ 0        Ω exp(iφ)  Ω exp(iφ)   -2Δ + 2V   ]
```

`omega`, `delta`, and `V` are finite real angular frequencies in **rad/s**;
`phase` is a finite real angle in **radians**. Signed frequencies are accepted.
Time is in **seconds**. Evolution therefore uses `exp(-i Hbar t)`.
The matrix convention and signs above are preserved without energy shifts.

In `quantum_project.units`, `mhz_to_angular_frequency(f)` returns `2π × 10⁶ × f`,
`angular_frequency_to_mhz(omega)` returns `omega / (2π × 10⁶)`, and
`microseconds_to_seconds(t)` returns `10⁻⁶ × t`. MHz means ordinary frequency
(cycles per microsecond), not angular frequency; for example, 1 MHz corresponds
to `2π × 10⁶ rad/s`. These helpers take finite real scalars.

`normalise_state(psi)` explicitly normalises a finite, nonzero, length-four
state vector. `state_probability(psi, basis_state)` and `fidelity(psi, target_state)`
both evaluate `|<target|psi>|²` for normalised pure states. They reject
unnormalised inputs rather than silently rescaling them. Probabilities can be
computed for any normalised projection target, including a Bell state.
Validation helpers `is_hermitian`, `is_normalised`, and `are_orthonormal` use an
absolute tolerance of `1e-12` with zero relative tolerance. Malformed shapes or
nonfinite inputs raise `ValueError`; valid inputs failing a check return `False`.
`are_orthonormal` takes one state per row, for example `BASIS_STATES`.

## Time-independent evolution

`evolve_static(psi0, times, Hbar)` implements Chapter 5 Eq. (5.6) for a constant
Hamiltonian. It calls `numpy.linalg.eigh` once. NumPy returns eigenvectors as
**columns** of `R`, so `Hbar = R @ diag(eigenvalues) @ R.conj().T` and:

```text
c = R.conj().T @ psi0
psi(t) = R @ (exp(-1j * eigenvalues * t) * c)
```

The implementation computes phases for all times together. Since output states
are stored as rows, reconstruction is `(phases * c) @ R.T`; the initial
projection uses the conjugate transpose. Complex amplitudes and their phases
are retained, without renormalising the output.

The initial state must be finite, normalised and have shape `(4,)`. `Hbar` must
be finite, Hermitian and have shape `(4, 4)`. `times` is a one-dimensional array
of finite real seconds measured from the initial state at **t = 0**, not from
the first requested time. Negative, repeated and unsorted times are accepted;
output follows the supplied order with shape `(len(times), 4)`. An empty time
array returns shape `(0, 4)`.

Using `Hbar = H / ħ` makes each eigenvalue an angular frequency in rad/s.
Multiplication by seconds gives the dimensionless phase in `exp(-i lambda t)`;
no extra factor of ħ or 2π belongs in the solver.

The independent validation reference `evolve_expm_reference` in
`tests/test_static_solver.py` directly applies `scipy.linalg.expm(-1j * Hbar * t)`
at each time. It does not call the eigenmethod. Tests compare complex amplitudes
with absolute tolerance `1e-11` and zero relative tolerance across all eight
zero/nonzero combinations of detuning, interaction and laser phase, three
initial states, and 81 times over 0–4 μs. The Rabi frequency is 1 MHz, detuning
is 0 or −0.4 MHz, interaction is 0 or 5 MHz (all converted to rad/s), and phase
is 0 or 0.73 radians. Run `python -m pytest -rP` to show the maximum measured
amplitude discrepancy and norm deviation. Roundoff in eigenvalues can accumulate
phase error for very large frequency–time products; the tolerance is validated
for these representative cases, not arbitrary propagation durations.

## Setup and tests

Use Python 3.10 or newer. From the repository root on Linux/macOS:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`.
`pyproject.toml` is the canonical dependency definition. NumPy, SciPy and
Matplotlib are runtime dependencies; the test extra contains pytest. SciPy is
used for exponential-step propagation and independent validation references.
No static checker is configured yet.

Example (after installation):

```python
from quantum_project import KET_0R, PSI_PLUS, rydberg_hamiltonian, state_probability
from quantum_project.units import mhz_to_angular_frequency

hbar = rydberg_hamiltonian(
    omega=mhz_to_angular_frequency(1.0),
    delta=0.0,
    phase=0.3,
    V=mhz_to_angular_frequency(10.0),
)
probability = state_probability(PSI_PLUS, KET_0R)  # 0.5
```

Static evolution example:

```python
import numpy as np
from quantum_project import KET_00, evolve_static
from quantum_project.units import microseconds_to_seconds

times = np.linspace(0, microseconds_to_seconds(4), 81)
states = evolve_static(KET_00, times, hbar)  # shape (81, 4)
populations = np.abs(states) ** 2
```

## Running the project

After setup, activate the environment and start the interactive launcher:

```bash
source .venv/bin/activate
python run_project.py
```

Or run a task directly:

```bash
python run_project.py tests
python run_project.py fig51
python run_project.py fig52
python run_project.py milestone
```

The complete milestone runs tests first, stops on any failure, generates both
figures, and lists the output files. Figure 5.2 uses the validated 131,072-interval
production grid and may take some time; the menu asks for confirmation, while a
direct command starts immediately. Use `python run_project.py results` to list
existing outputs or `python run_project.py --help` for usage. The `zeno` and
`symmetry` commands are placeholders for extensions that are not yet implemented.

## Time-dependent propagation

`evolve_trotter(psi0, time_grid, hamiltonian_function)` implements Eq. (5.7) with
`scipy.linalg.expm(-1j * Hbar(t_n) * (t_{n+1} - t_n))`. It evaluates the callback
separately at the **beginning of each interval**, then propagates chronologically.
The callback accepts seconds and returns H/ħ in rad/s; Ω, Δ, φ and V may all vary.

The grid must be finite, real, nonempty and strictly increasing; nonuniform
steps are supported. Here `psi0` is defined at `time_grid[0]`, and the output
has shape `(len(time_grid), 4)`, including the initial state. A one-point grid
returns just that state. The routine validates each Hamiltonian and preserves
complex amplitudes without renormalisation.

For constant H, the exponential product is exact apart from floating-point
error. For varying H(t), freezing the Hamiltonian at each left endpoint is a
piecewise-constant approximation with first-order global accuracy for smooth
time dependence. Smaller timesteps are needed to establish convergence, even
when the state norm is conserved. No midpoint method is substituted.

See [the method and convergence studies](docs/time_dependent_validation.md).
Run `python -m pytest tests/test_trotter.py -rP` to show the constant-H comparison
and detuning-ramp refinement results. For example:

```python
from quantum_project import evolve_trotter

time_grid = np.linspace(0, 1e-6, 321)

def hbar_at(t):
    return rydberg_hamiltonian(
        omega=mhz_to_angular_frequency(1.0),
        delta=mhz_to_angular_frequency(t / 1e-6 - 0.5),
        phase=0.37,
        V=mhz_to_angular_frequency(3.0),
    )

states = evolve_trotter(KET_00, time_grid, hbar_at)
```

## Figure 5.1 and analytical validation

```bash
python scripts/figure_5_1.py
```

The script uses the validated static solver, starting in `|00>` with
`Omega/(2π) = 1 MHz`, zero detuning and phase, and `V/(2π) = 0, 0.5, 1, 5, 100 MHz`.
It samples 0–1 microsecond at 4,001 points and saves two vertically arranged
panels for `P_00` and `P_Psi+` in `results/figures/figure_5_1.pdf` and `.png`.
Plotting is separate from the solver. The source range is sampled at the documented representative strengths;
exact colour density, typography and individual sampling are not source requirements.

See [the analytical derivations and validation procedure](docs/static_validation.md)
for the exact noninteracting and ideal blockade solutions. Additional small and
large interactions test convergence to those limits. Residual figures are saved
as `results/figures/static_residuals.pdf` and `.png`. Complex states, all four
basis populations, Bell probabilities and analytical residuals are retained in
`results/data/figure_5_1.npz`; parameters, maximum/RMS residuals, physical checks
and sampling-refinement results are in `results/data/static_validation.json`.
The script prints the residual table and fails if conservation, symmetry or
probability bounds are violated beyond `1e-12`.

Implementation lives under `src/quantum_project/`, with tests under `tests/`.
Deliberate result files are not ignored. No commit, push or tag is part of this run.

## Figure 5.2 sweep study

See [sweep profiles, reference method and accuracy criterion](docs/sweep_validation.md).
Regenerate the final milestone outputs with:

```bash
python scripts/figure_5_2.py --milestone --interaction-mhz 0 0.5 1 5 100 --phase 0
```

The command fixes Ω/(2π) = 1 MHz, Δmax/(2π) = 50 MHz and duration 1 μs.
It regenerates the full 100 MHz convergence benchmark and checks all other
sampled interactions at a common production timestep and half that timestep,
against an independently refined DOP853 reference. A stricter timestep is
selected automatically if required by any displayed case. The family check
requires P00 and PΨ+ absolute errors and step-halving differences ≤5e−5;
it does not establish the benchmark's largest passing timestep for every V.
The three panels are detuning in MHz, P00 and PΨ+.

Final figures, linear/cubic comparison, convergence results and raw data are
under `results/figures/` and `results/data/`. The historical Run 6 audit remains
under `results/audit_sweep/`. See the milestone report for measured results.

## Reproduction and validation

From the repository root on Linux/macOS:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install '.[test]'
python -m pytest -rP
python -m compileall -q src scripts tests
python -m pip check
python scripts/figure_5_1.py
python scripts/figure_5_2.py --milestone --interaction-mhz 0 0.5 1 5 100 --phase 0
git diff --check
git status --short --untracked-files=all
```

Both scripts save PDF/PNG figures and numerical data without clipping or
renormalising the propagated states. The sweep command can take several minutes.
For a restricted environment, set `MPLCONFIGDIR` to a writable temporary directory.
The automated phase tests cover static, cubic and linear evolution for all five
interactions and phases 0.37, π/2 and −1.2 radians, compared with zero phase.
No staging, commit, push or tag is part of this run.
