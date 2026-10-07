# PHYS3561 Chapter 5 milestone readiness audit

**Run 6B: source interpretation resolved; final numerical validation recorded below.**
The user checked Chapter 5 directly: Figure 5.1 specifies Ω/(2π) = 1 MHz,
blue V = 0, and a red-curve interaction range 0–100 MHz, without enumerating
individual red-curve values. Figure 5.2 specifies the cubic sweep, 50 MHz
maximum detuning and 1 μs duration; other parameters are as Figure 5.1.
Both figures use the documented representative samples 0, 0.5, 1, 5, 100 MHz.
These are repository numerical choices, not individual values listed in the
booklet. No exact colour density or typography is claimed. There is no missing
source interaction list to recover, and no outstanding source-parameter blocker.

We choose constant φ = 0. With D = diag(1, exp(iφ), exp(iφ), exp(2iφ)),
Hbar(φ) = D Hbar(0) D† and ψφ(t) = D ψ0(t) because D†|00> = |00>.
The Bell projection gains a common phase exp(iφ), leaving its modulus squared
and P00 unchanged. This applies to the symmetric Hamiltonian and constant
common phase; varying phase or a different initial state needs separate analysis.
Automated tests check matrix equivalence, complex-state covariance and raw
probabilities for φ = 0.37, π/2, −1.2, all five V values, and static/cubic/linear
propagation. The original Hamiltonian and solvers are unchanged.

## 1. Physical model

Two identical coherently driven atoms each have states |0> and |r>. Double
Rydberg excitation adds an interaction energy ħV. The closed, four-dimensional
pure-state model includes no spontaneous decay, dephasing, motion or noise.
All milestone trajectories start in |00>. The Bell target is the symmetric
state |Ψ+> = (|0r> + |r0>)/√2, and fidelity is |<Ψ+|ψ>|².

## 2. Basis and Hamiltonian

The fixed basis is |00>, |0r>, |r0>, |rr>. The single implementation in
`src/quantum_project/hamiltonian.py` returns


\bar H=\frac12\begin{pmatrix}
2\Delta&\Omega e^{-i\phi}&\Omega e^{-i\phi}&0\\
\Omega e^{i\phi}&0&0&\Omega e^{-i\phi}\\
\Omega e^{i\phi}&0&0&\Omega e^{-i\phi}\\
0&\Omega e^{i\phi}&\Omega e^{i\phi}&-2\Delta+2V
\end{pmatrix}.
\]

The audit confirms Hermiticity, negative upper-coupling phase signs, positive
conjugate lower couplings, diagonal entries Δ, 0, 0, −Δ+V, and the interaction
only on |rr>. Tests compare every matrix entry, including φ = π/2. Canonical
states are orthonormal and the Bell state is normalised. No duplicate production
Hamiltonian or change to the supplied physics was found or introduced.

## 3. Units

The stored matrix is Hbar = H/ħ, so Ω, Δ and V are angular frequencies in rad/s;
phase is in radians and time is in seconds. Ordinary frequency f in MHz becomes
2π × 10⁶ f rad/s. Conversion back divides by the same factor. One microsecond
is 10⁻⁶ seconds. Thus exp(−i Hbar t) has a dimensionless exponent. No hidden ħ
or additional 2π appears in either propagator. Plotted MHz are Δ/(2π) or V/(2π),
with the necessary 10⁶ conversion.

## 4. Static numerical method

`evolve_static` calls `numpy.linalg.eigh` once per Hamiltonian. Its eigenvectors
are columns of R: Hbar = R diag(λ) R†. Initial coefficients are R†ψ0, acquire
phases exp(−iλt), and are reconstructed with R. Since returned states are rows,
the implementation uses `(phases * coefficients) @ R.T` for reconstruction.
It does not use a general eigensolver or diagonalise inside a time loop.

The initial state is defined at t = 0. Arbitrary requested times, including
negative or repeated times, are supported. Output is complex with shape
`(number_of_times, 4)`; inputs must be finite, dimensionally valid and normalised.

## 5. Static validation

Independent SciPy matrix exponentials agree with the static eigensolver to
3.389664e−14 maximum complex-amplitude difference in the representative Run 2
study. For Δ = φ = 0, analytical limits follow from independent single-atom
rotations at V = 0 and the symmetric two-state blockade system as |V/Ω| → ∞:

| Limit | P00 | PΨ+ |
| --- | --- | --- |
| V = 0 | cos⁴(Ωt/2) | ½ sin²(Ωt) |
| Ideal blockade | cos²(Ωt/√2) | sin²(Ωt/√2) |

For the 0–1 μs calculation, the measured probability residuals are:

| Case | Observable | Maximum absolute | RMS |
| --- | --- | ---: | ---: |
| V = 0 | P00 | 9.9920e−16 | 3.8206e−16 |
| V = 0 | PΨ+ | 5.5511e−16 | 1.8878e−16 |
| V/(2π) = 100 MHz vs ideal blockade | P00 | 8.1487e−5 | 3.8635e−5 |
| V/(2π) = 100 MHz vs ideal blockade | PΨ+ | 1.1551e−4 | 5.7081e−5 |

Near-zero and both signs of large interaction are included in the residual
study. Finite-V deviations are approximation errors, not propagation errors;
they decrease towards the limits. See `static_validation.md` for derivations,
all cases and doubled-grid checks.

All four populations sum to one within 1.9984e−15 across the static study.
The largest |P0r−Pr0| is 5.5511e−16. The raw probability range is
[0, 1.0000000000000013], consistent with roundoff. Both analysis scripts use
unclipped probabilities. The scalar fidelity helper bounds rounding excursions
after checking normalised inputs, but is not used to conceal audit residuals.

## 6. Figure 5.1

`scripts/figure_5_1.py` uses Ω/(2π) = 1 MHz, Δ = φ = 0, initial |00>, and
V/(2π) = 0, 0.5, 1, 5, 100 MHz. The 4,001-point grid spans 0–1 μs.
The upper panel shows P00; the lower shows PΨ+. A separate residual figure
compares exact noninteraction and ideal blockade. PDF, PNG, complex states,
basis probabilities and numerical residual reports are saved under `results/`.

The trends agree with the derived physics: independent atoms have Bell
probability at most 1/2, while strong blockade gives collective √2Ω oscillations
and nearly unit Bell probability. The individual curve samples are numerical choices within the source range,
not an enumerated source list; styling distinguishes blue noninteraction from
the red interacting family.

## 7. Time-dependent numerical method

`evolve_trotter` implements Eq. (5.7) directly:

\[
\psi_{n+1}=\exp[-i\bar H(t_n)(t_{n+1}-t_n)]\psi_n.
\]

`scipy.linalg.expm` acts on the full matrix. The callback is evaluated separately
once per interval at its left endpoint; later propagators multiply after earlier
ones. No term splitting or midpoint replacement is used. The state is defined
at the first grid time; the nonempty grid must be finite and strictly increasing.
Nonuniform steps are allowed. Every Hamiltonian is checked for Hermiticity.

Constant-H products are exact mathematically, so their errors are roundoff:
maximum discrepancy from the static solver is 1.021941e−13, with maximum norm
deviation 4.640732e−14. These errors need not decrease monotonically as the
number of steps increases. A noncommuting chronological-order test uses
independent analytical drive and detuning solutions. A separate detuning-ramp
test halves its error when dt halves, confirming first-order global convergence.
Neither solver silently renormalises or clips its evolving state.

## 8. Figure 5.2

The exact requested cubic detuning is implemented:

\[
\Delta(t)=8\Delta_{\max}[(t-T/2)/T]^3,
\quad\Delta_{\max}=2\pi\times50\times10^6\,\mathrm{rad/s},\quad T=1\,\mu s.
\]

It starts at −Δmax, crosses zero at T/2 and ends at +Δmax. Tests cover endpoints,
midpoint, symmetry and deterministic repeatability. `scripts/figure_5_2.py`
produces detuning, P00 and PΨ+ in three vertical panels and saves raw data.

The final milestone uses Ω/(2π) = 1 MHz, Δmax/(2π) = 50 MHz, T = 1 μs,
φ = 0 and representative V/(2π) = 0, 0.5, 1, 5, 100 MHz. The three panels
are Δ/(2π) in MHz, P00 and PΨ+, with consistent blue/noninteracting and
red/interacting colours in Figures 5.1 and 5.2. Final outputs and both-profile
raw states are in `results/figures/` and `results/data/`.
The historical single-interaction Run 6 audit in `results/audit_sweep/` remains
available; its original provenance labels describe that historical run.

## 9. Timestep convergence

The complete sweep script compares the original exponential-step solver against
a fine independent DOP853 reference. The reference is refined in tolerance and
maximum step and checked on 40,001 points. Its dense solution is sampled at
the exact trial grid times; trial populations are never interpolated or shifted.
The full 100 MHz benchmark trial step counts range from 131072 down to 256 over 1 μs.

For positive scale S, the four-significant-place absolute tolerance is
0.5 × 10^(floor(log10(S)) − 3). Trajectory errors use each observable's peak
scale, which avoids undefined relative errors at zeros. Final and peak Bell
fidelities use their own scales and must additionally round identically with
`.4g`. This does not guarantee pointwise relative accuracy near zero or
four-significant-figure accuracy of peak time. Maximum Bell-fidelity difference
equals maximum PΨ+ difference because they are the same observable.

The result is the largest **tested** timestep passing every criterion, with
all finer tested steps also passing. No claim is made about untested timesteps.
Run 6B regenerates the full 100 MHz benchmark in the normal `results/`
directories. For the other V values, a cheaper two-grid check compares the
production step and half that step against refined DOP853, requiring absolute
P00/PΨ+ errors ≤5e−5 and a common-grid step-halving difference ≤5e−5.
The common production step is reduced and all cases rechecked if necessary.
The benchmark's largest passing threshold and scalar-rounding guarantees are
not asserted for every V. See `sweep_validation.md` for both criteria.

For the **100 MHz, φ = 0 audit benchmark**, both profiles pass at
dt = 1.52587890625e−11 s = 0.0152587890625 ns (65,536 intervals).
The next tested step, 0.030517578125 ns, fails. The finest tested step,
0.00762939453125 ns, also passes. At the selected step:

| Sweep | Final Bell error | Maximum P00 error | Maximum PΨ+/Bell error | Peak Bell error |
| --- | ---: | ---: | ---: | ---: |
| Cubic | 5.256819e−7 | 3.763053e−5 | 3.762879e−5 | 3.229188e−8 |
| Linear | 5.196807e−7 | 4.923259e−5 | 4.922970e−5 | 2.405956e−5 |

Reference refinement changes probabilities by at most 8.779644e−13 (cubic)
and 8.864021e−13 (linear). Worst probability-sum error over all benchmark
timesteps is 1.549872e−12, below the 1e−10 long-propagation audit tolerance.
These data support the chosen criterion for this benchmark only.

## 10. Linear versus cubic sweep

The linear profile is Δ(t) = 2Δmax(t/T−1/2), with identical endpoints and duration.
The script compares final and maximum Bell fidelity, time of maximum, and both
probability trajectories for each interaction. The cubic derivative vanishes
at central resonance, while the linear derivative is constant. Spending more
time near resonance can improve transfer, but is not proof of adiabaticity or
a guarantee of greater final Bell fidelity. If Δ reaches V, double-excitation
transfer can matter. The 100 MHz benchmark provides a numerical example, not a universal profile ranking.

In the benchmark, cubic final fidelity rounds to 0.9985 and linear to 0.1921;
their maximum fidelities round to 0.9999 and 0.2479. The large physical difference
greatly exceeds their measured numerical errors. It is consistent with the
cubic sweep spending longer near the central avoided crossing, but no claim
of perfect adiabatic following or optimal control is made.

## 11. Numerical limitations

- Individual interaction samples are numerical choices within the source range;
  exact source colour density or typography is not asserted.
- Exponential steps conserve norm while still having finite discretisation
  error for varying H. Left-endpoint propagation is first order globally.
- Constant-H roundoff can grow with more steps; long frequency–time products
  can lose phase accuracy. No output renormalisation masks these effects.
- Sweep convergence maxima are sampled-grid maxima, not rigorous continuous
  suprema. Reference peak times use a dense scan and local optimisation;
  candidate peak times are limited by their grid spacing.
- Four significant places on a trajectory's peak scale differ from pointwise
  relative accuracy near zero. The criterion is explicitly documented.
- The model assumes ideal unitary closed-system dynamics.
- Dependencies have minimum versions rather than a lockfile; reports record
  the versions used. Numerical results, PDF metadata and bytewise output hashes
  can vary between environments. No static linter is configured.

The software audit found one production Hamiltonian and no dead experimental
scripts in the repository. Duplicate raw population calculation was replaced
by the existing shared helper; key sweep reference tolerances were named.
Remaining prints are requested analysis tables/progress, not debug output.
The core Hamiltonian and both propagators were not changed.

## 12. Reproduction and validation

From the repository root with the installed project environment:

```bash
source .venv/bin/activate
python -m pytest -rP
python -m compileall -q src scripts tests
python -m pip check
python scripts/figure_5_1.py
python scripts/figure_5_2.py --milestone --interaction-mhz 0 0.5 1 5 100 --phase 0
git diff --check
git status --short --untracked-files=all
```

Both scripts use raw probabilities; neither propagated states nor saved
probabilities are clipped or renormalised. Every exponential step validates
Hermiticity. Sweep trials enforce norm conservation, probability sum, exchange
symmetry and probability bounds within 1e−10, allowing accumulated roundoff
in long runs. Static probability checks use 1e−12. Independent references are
refined in tolerance and maximum step and checked on 40,001 points.
The artifact audit checks every numeric array for finiteness and recomputes
saved populations from complex amplitudes. Automated tests verify phase
invariance. Nothing is staged, committed, pushed or tagged.

Run 6 used a fresh wheel installation; Run 6B uses the existing `.venv` and
source checkout. The exact environment versions and numerical measurements
are stored with the results. The latest validation summary and full working-tree
inventory are recorded in `results/data/run6b_validation.json`.

## 13. Run 6B measured results

Completed on 2026-10-07 after resuming the interrupted final validation.
Figure 5.1 was regenerated from the current script. A fresh sweep run completed
both full 100 MHz benchmarks, but a server restart interrupted its family pass
before it replaced any saved sweep data. The complete earlier Run 6B sweep
NPZ/JSON/CSV outputs were recovered and audited; all sweep figures were rebuilt
from those arrays using the current plotting code. Figure 5.1 and Figure 5.2
use blue for V = 0 and red shades for the interacting samples.
All five PNG plots were visually inspected for labels, legends and layout.

All **254 tests passed**, including 45 phase-invariance comparisons across
static/cubic/linear evolution, five interactions and three nonzero phases.
The largest measured phase-induced P00/PΨ+ difference was
1.099121e-14. Byte-compilation, dependency consistency and whitespace
checks passed. Static generation and sweep-figure rebuilding emitted no warnings.
An initial pip-cache permission warning was resolved by using a writable
temporary cache; the repeated dependency check was clean. The artifact audit's
Matplotlib import also used an automatic temporary cache after a default-path
permission warning; all numerical checks passed.

The full 100 MHz benchmark retains the largest tested passing timestep
0.0152587890625 ns for both profiles. The common production grid uses
**131,072 intervals**, dt = **0.00762939453125 ns**,
as recorded in the recovered report. The number of refinement attempts was not
recorded in that report. Each of the ten profile/interaction cases passes at
that step and half that step. Across those checks, the largest
P00/PΨ+ error against refined DOP853 is 2.927660e-05; the largest
recorded family step-halving difference is 1.463893e-05, both below 5e−5.
Reference refinement changes probabilities by at most 8.868462e-13.
The full scalar-rounding criterion remains specific to the 100 MHz benchmark.

Production-grid Bell results (digits shown for reproducibility, not a claim
of equal significant-digit accuracy for every interaction):

| Sweep | V/(2π) (MHz) | Final Bell fidelity | Maximum Bell fidelity | Sampled peak time (μs) | Max P00/PΨ+ reference error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cubic | 0 | 0.303538227 | 0.500000000 | 0.578636169 | 1.531659e-05 |
| Cubic | 0.5 | 0.296336192 | 0.560295828 | 0.600357056 | 2.284155e-05 |
| Cubic | 1 | 0.360388694 | 0.642315616 | 0.627571106 | 2.927660e-05 |
| Cubic | 5 | 0.732787147 | 0.948544354 | 0.691406250 | 2.004609e-05 |
| Cubic | 100 | 0.998543580 | 0.999922981 | 0.834007263 | 1.880695e-05 |
| Linear | 0 | 0.182989934 | 0.228954401 | 0.586021423 | 2.343690e-05 |
| Linear | 0.5 | 0.181922073 | 0.227822883 | 0.585327148 | 2.356092e-05 |
| Linear | 1 | 0.182630416 | 0.226951277 | 0.584648132 | 2.363067e-05 |
| Linear | 5 | 0.169757061 | 0.229649026 | 0.582221985 | 2.403142e-05 |
| Linear | 100 | 0.192078292 | 0.247908423 | 0.585502625 | 2.461378e-05 |

The final artifact audit loaded every numeric array in both current NPZ files,
verified finiteness, checked all production grids, recomputed the basis and
coherent Bell probabilities directly from complex amplitudes, and reconciled
selected errors and fidelities with JSON and all 36 CSV rows.
The largest probability-reconstruction difference was
6.661338e-16; the largest saved-state
probability-sum error was 1.549871e-12,
and exchange-symmetry error was 7.821521e-14.
Static and sweep conservation checks satisfy their respective 1e−12 and 1e−10
tolerances without clipping or renormalisation.

`results/data/run6b_validation.json` records the validation script, commands,
test and generation output, library versions, array inventory, source and
artifact SHA-256 hashes, and full working-tree status. The historical Run 6
outputs remain in `results/audit_sweep/`. No files were staged, committed,
pushed or tagged. The interview checklist remains personal preparation.
