# Figure 5.2: sweep profiles and convergence procedure

## Source parameters and numerical choices

The user directly checked the Chapter 5 source for Run 6B. Figure 5.1 specifies
Ω/(2π) = 1 MHz, blue V = 0, and red curves spanning V/(2π) = 0–100 MHz.
It does not enumerate individual red-curve interactions. Figure 5.2 specifies
the cubic detuning below, Δmax/(2π) = 50 MHz and T = 1 μs, and carries over
the other Figure 5.1 parameters. The repository samples that range at
V/(2π) = 0, 0.5, 1, 5, 100 MHz for both figures. These are documented numerical
choices; no exact source colour density or individual curve list is claimed.

We choose constant φ = 0. With
D = diag(1, exp(iφ), exp(iφ), exp(2iφ)), Hbar(φ) = D Hbar(0) D†.
Since D is constant and D†|00> = |00>, ψφ(t) = D ψ0(t), including for varying
Δ(t) and for every exponential step. The |00> amplitude is unchanged and
the Bell projection gains only exp(iφ), so P00 and PΨ+ are invariant.
This argument requires a common constant phase and the specified initial state;
it does not apply without qualification to time-dependent phase or other states.
Automated tests compare phases 0.37, π/2 and −1.2 with zero for static, cubic
and linear evolution at every displayed V, including complex-state covariance.

The historical Run 6 benchmark is retained in `results/audit_sweep/`.
The final milestone outputs use the normal `results/` directories.

## Exact profiles and propagation

`src/quantum_project/sweeps.py` implements

\[
\Delta_{\rm cubic}(t)=8\Delta_{\max}\left(\frac{t-T/2}{T}\right)^3,
\qquad
\Delta_{\rm linear}(t)=2\Delta_{\max}\left(\frac tT-\frac12\right).
\]

Both equal −Δmax at zero, zero at T/2, and +Δmax at T. Functions accept scalar
or array times in seconds and return rad/s. They reject times outside [0,T],
nonfinite values and invalid durations. The requested Δ/(2π) plot divides
angular frequency by 2π × 10⁶ to display MHz.

The script uses the unchanged `evolve_trotter`: full `scipy.linalg.expm` at
each left-endpoint Hamiltonian, with no midpoint replacement. It computes
the four basis populations and the coherent projection onto |Ψ+>.
`F_Bell(t)` is exactly the same observable as `P_Psi_plus(t)`; their trajectory
errors are reported separately under the requested names but are identical.

## Reference and time-grid comparison

For each profile and interaction, the independent fine reference integrates
the complex Schrödinger equation with SciPy DOP853. First it uses rtol = 1e−11,
atol = 1e−13, maximum step T/2000, then rtol = 1e−13, atol = 1e−15, maximum
step T/4000. The refined reference is used for error measurements. Their
probabilities must agree within 1e−9 on a 40,001-point audit grid.

Each trial trajectory uses the requested exponential-step method. The reference
is evaluated at precisely the trial grid's times using DOP853's own continuous
interpolant. Trial probabilities are never interpolated, shifted or phase-aligned.
Reported suprema are maxima on the trial's sampling grid, not rigorous bounds
over continuous time. The fine grids around the acceptance boundary resolve
the physical oscillations; sampled peak-time uncertainty is approximately one
timestep and is reported separately.

For the rigorous 100 MHz benchmark, the trial step counts are 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024,
512 and 256 over 1 μs. This tests progressively larger timesteps, starting at
approximately 0.00763 ns and ending at 3.90625 ns. These candidates can be
refined if the finest fails. The largest passing timestep is **the largest
tested value**, with every finer tested value also passing, not a universal
maximum or an assertion about untested timesteps. Coarse accidental error
cancellation is not sufficient.

For each case the script reports final Bell fidelity, final fidelity error,
maximum absolute P00 trajectory error, maximum absolute PΨ+ trajectory error,
maximum Bell-fidelity difference and error in the maximum Bell fidelity.

## Operational four-significant-figure criterion

“Four significant figures” is ambiguous for an entire trajectory that crosses
zero: relative error at zero is undefined, and rounding identical values near
a decimal rounding boundary is stricter than a simple error threshold.
This project therefore defines the criterion explicitly:

- For a positive reference scale S, half a unit in the fourth significant
  place is ε(S) = 0.5 × 10^(floor(log10(S)) − 3). For example, S = 0.9 gives
  ε = 5e−5 while S = 0.09 gives ε = 5e−6; these are not four decimal places.
- Each trajectory's maximum absolute error must be no larger than ε(S), with
  S its reference peak probability. This is four-significant-place accuracy
  on the observable's scale, not pointwise relative accuracy near zeros.
- Final and maximum Bell fidelities must satisfy their own scalar ε(S)
  bounds **and** round identically to their reference values with `.4g`.
- All of the conditions must pass. Peak time is reported with its grid
  resolution but is not claimed accurate to four significant figures by this
  probability criterion.

The dense reference's peak is found by scanning 40,001 points and then locally
optimising around the highest interior sample, also checking both endpoints.
The trial maximum remains the sampled maximum; comparisons include the
resulting missed-peak error. The peak search is numerical, not visual.

Probabilities are never clipped or renormalised. Every trial and the reference
are checked for probability conservation, exchange symmetry and bounds with
absolute tolerance 1e−10. This allows accumulated roundoff over more than
100,000 unitary steps while remaining far below the probability error targets.

## Physical interpretation

The cubic profile has derivative
24 Δmax (t/T − 1/2)² / T, which vanishes at the central resonance.
The linear derivative is the constant 2 Δmax / T. The cubic pulse therefore
spends more time near resonance at the same total duration and endpoints,
which can improve excitation transfer for some parameters. This alone does
not prove adiabatic following or guarantee higher Bell fidelity.

Interaction also matters: for the specified Hamiltonian, the symmetric
singly-excited state and |rr> are resonant near Δ = V. A sweep reaching that
region can promote double excitation instead of retaining Bell population.
For V sufficiently beyond the sweep range, blockade suppresses this channel.
The finite initial detuning also means |00> is only approximately an energy
eigenstate when the drive is already on. Final fidelity, transient maximum
fidelity and the time of that maximum must therefore be reported separately.

The profiles can require different timesteps because their derivatives and
resulting trajectories differ. A higher fidelity is a physical outcome;
smaller discretisation error is a separate numerical issue established by
the convergence study.

## Reproduction and outputs

To reproduce the final milestone, run:

```text
python scripts/figure_5_2.py --milestone --interaction-mhz 0 0.5 1 5 100 --phase 0
```

The command saves:

- `results/figures/figure_5_2.pdf` and `.png`: cubic detuning, P00 and PΨ+
  in three vertically arranged panels, with time in microseconds.
- `results/figures/sweep_comparison.pdf` and `.png`: both profiles and populations.
- `results/figures/sweep_convergence.pdf` and `.png`: error/tolerance ratio and
  final fidelity error versus timestep. The final/peak rounding checks in the
  report also contribute to acceptance, though not to the plotted error ratio.
- `results/data/figure_5_2.npz`: selected trial times, complex states, detuning,
  all populations and matching reference values. Keys use profile and
  zero-based interaction index; `interaction_mhz` provides the mapping.
- `results/data/sweep_convergence.csv`: trial-step metrics and pass/fail flags.
- `results/data/sweep_validation.json`: parameters, library versions, reference
  checks, fidelity summaries and complete convergence records.

If no timestep passes, the data are saved but final figures are withheld and
the script raises an error instructing refinement. `--output-dir` can direct
trial runs elsewhere without replacing the deliberate project results.
`--case-label` records parameter provenance in the figures and JSON report.

## Production timestep verification across the family

Milestone mode first regenerates the full 100 MHz benchmark for both profiles.
Starting from the stricter of its selected timesteps, it checks each displayed V
at that step and half that step. At V < 100 MHz, both P00 and PΨ+ must agree
with a separately refined DOP853 reference to absolute error ≤5e−5 on each
trial grid, and the two trial trajectories must agree to ≤5e−5 at common times.
The 100 MHz case retains the full four-significant-place and scalar-rounding
criterion. If any case fails, the common production timestep is halved and
all cases are checked again. Thus the final step is verified for the displayed
family; the 100 MHz acceptance boundary is not asserted to be identical at
other V values. Family checks do not certify scalar four-significant-digit
rounding or a largest admissible timestep for the smaller interactions.

JSON separates `benchmark_convergence` from the two-grid production checks in
`convergence`. The CSV and convergence plot contain the full 100 MHz series
and the two-grid checks for smaller V; their acceptance scopes are distinct.
The plot's error/tolerance ratio retains the stricter benchmark metric for
comparison, even where family acceptance uses the 5e−5 plotting criterion.
Each saved production case uses the same timestep. Final and peak values,
reference errors and raw conservation diagnostics are recorded in JSON.
