# Time-dependent propagation: Eq. (5.7)

`evolve_trotter(psi0, time_grid, hamiltonian_function)` implements the supplied
Chapter 5 update with a **left-endpoint Hamiltonian**:

\[
\psi_{n+1}=U_n\psi_n,\qquad
U_n=\exp[-i\bar H(t_n)(t_{n+1}-t_n)],\qquad \bar H=H/\hbar.
\]

The callback receives seconds and returns a finite Hermitian 4 × 4 matrix in
rad/s, in the basis \((|00\rangle,|0r\rangle,|r0\rangle,|rr\rangle)\).
The dimensionless exponent needs no extra ħ or 2π. The callback can construct
the Hamiltonian from any time-dependent combination of Ω, Δ, φ and V.

The grid is nonempty, one-dimensional, finite, real and strictly increasing.
Nonuniform steps and any finite starting time are allowed. `psi0` is the
normalised state **at `time_grid[0]`**, which differs from the static solver's
convention of an initial state at zero. Constant-H comparisons therefore use
`evolve_static(psi0, time_grid - time_grid[0], Hbar)`.

Output has shape `(len(time_grid), 4)` and contains the complex state at every
grid point. Its first row is an exact copy of the initial state. A one-point
grid requires no propagation and makes no callback calls. Repeated or decreasing
times are rejected. Every subsequent callback result is checked for shape,
finiteness and Hermiticity using the existing validation helpers.

## Exponentials and chronological order

`scipy.linalg.expm` computes a full matrix exponential for each interval, and
the callback is evaluated separately once at each interval's beginning.
The final grid time is not sampled, because it starts no interval. There is no
cached initial Hamiltonian and no midpoint substitution. No optional midpoint
method is implemented in this run.

Later propagators act after earlier ones:

\[
\psi_N=U_{N-1}U_{N-2}\cdots U_0\psi_0.
\]

In general, Hamiltonians at different times do not commute, so reversing this
order changes the answer. Tests apply a resonant drive followed by a detuning
interval and compare against independently derived product-state rotations and
diagonal phases. This also checks nonuniform spacing and absolute callback times.

Although the module is named `trotter`, the implementation does not split the
Hamiltonian into separately exponentiated terms. It uses the full Hamiltonian
and the time-slicing prescription of Eq. (5.7). For Hermitian Hbar, each exact
step is unitary. States are not renormalised during propagation; measured norm
errors remain visible.

## Constant Hamiltonian: validation before time dependence

For constant \(\bar H\), all interval propagators commute and

\[
\prod_n e^{-i\bar H\Delta t_n}
=e^{-i\bar H\sum_n\Delta t_n}
=e^{-i\bar H(t_N-t_0)}.
\]

Thus there is **no timestep truncation error for constant H**. Even coarse steps
should agree closely with the static eigensolver. Refinement need not make the
remaining floating-point error decrease monotonically: it also adds matrix
exponentials and matrix-vector products. A test demanding monotonic improvement
at the roundoff floor would be inappropriate.

The constant-H study was run and passed before any varying-H validation. It uses
Ω/(2π) = 1 MHz and the following (Δ/(2π) in MHz, V/(2π) in MHz, φ in radians):

- (0, 0, 0)
- (0.4, 5, 0.73)
- (−0.4, 100, π/2)

Each is tested with |00>, |Ψ+>, and the normalised vector
`[1+2j, -3j, 4, -5+1j]`, over 0–1 μs. Complex amplitudes are compared directly,
including phase, at every grid point. Reported errors are maxima over all these
cases and times.

| Steps | dt (ns) | Maximum amplitude difference | Maximum norm deviation |
| ---: | ---: | ---: | ---: |
| 8 | 125 | 5.762261e-14 | 3.841372e-14 |
| 32 | 31.25 | 1.021941e-13 | 4.640732e-14 |
| 128 | 7.8125 | 9.712520e-14 | 4.107825e-14 |
| 512 | 1.953125 | 5.058038e-14 | 1.887379e-14 |

Amplitude tolerance is 1e-11 and norm tolerance is 1e-12, both absolute with no
relative tolerance. Numerical results can vary slightly by platform/library.

## Varying Hamiltonian: approximation and convergence

For smooth Hbar(t), freezing Hbar at the left endpoint omits its variation over
the interval. The local error is O(dt²), leading to **first-order global error**
O(dt) over a fixed duration on uniform grids. For nonuniform grids the largest
step must resolve the Hamiltonian variation. Discontinuities should be placed
on grid boundaries. Norm conservation alone does not show that a varying-H
trajectory is accurate; a coarse product of unitary steps can still be wrong.

Two independent references check convergence:

1. A commuting, linear detuning with Ω = V = 0 and
   Δ(t) = Ωref t/T, where Ωref = 2π × 10⁶ rad/s and T = 1 μs.
   The integrated detuning is exactly Ωref t²/(2T). Multiplying the |00> and
   |rr> amplitudes by opposite phases gives the exact continuous solution.
   Refining through 20, 40, 80 and 160 steps approximately halves the error.
2. A noncommuting Rydberg detuning ramp with Ω/(2π) = 1 MHz, V/(2π) = 3 MHz,
   φ = 0.37 radians, Δ(t) = Ω(t/T − 1/2), T = 1 μs and initial |00>.
   An independent `scipy.integrate.solve_ivp` DOP853 reference solves
   `dpsi/dt = -1j * Hbar(t) @ psi` with rtol = 1e-12, atol = 1e-14 and
   maximum reference step T/100. It does not use the exponential-step solver.
   The error is the maximum complex-amplitude discrepancy over the grid.

Measured noncommuting detuning-ramp results:

| Steps | dt (ns) | Maximum amplitude difference from DOP853 |
| ---: | ---: | ---: |
| 20 | 50 | 5.870228e-2 |
| 40 | 25 | 2.951230e-2 |
| 80 | 12.5 | 1.475576e-2 |
| 160 | 6.25 | 7.379883e-3 |
| 320 | 3.125 | 3.689913e-3 |

The error approximately halves when dt halves, consistent with first-order
convergence. These finite errors demonstrate why the accuracy of the constant-H
case cannot be assumed for a changing detuning. This is a method-validation
study, not a choice of production timestep for Figure 5.2. The timestep for that
calculation must be checked against its own accuracy requirements.

All 206 tests pass, including the 169-test pre-change baseline. The largest
constant-H discrepancy is 1.021941e-13 and its largest state-norm deviation is
4.640732e-14. The static solver source was not modified.

To reproduce the studies and display their numerical output:

```bash
python -m pytest tests/test_trotter.py -rP
python -m pytest
```

The static solver and existing Figure 5.1 outputs are preserved. This run adds
propagation and validation only; it does not reproduce Figure 5.2 or complete
the full milestone.
