# Static two-atom validation

This calculation follows the Hamiltonian and figure requirements supplied for
Chapter 5. Run 6B's direct source clarification confirms Ω/(2π) = 1 MHz, a blue V = 0
curve and red interacting curves spanning 0–100 MHz. The source does not list
individual red-curve values. The repository's 0, 0.5, 1, 5, 100 MHz sampling
is a numerical choice, without claiming exact colour density or typography.

## Conventions and symmetric subspace

Use the ordered basis \((|00\rangle,|0r\rangle,|r0\rangle,|rr\rangle)\),
\(\bar H=H/\hbar\), and angular frequencies in rad/s. Throughout this run,
\(\Delta=\phi=0\), \(\Omega=2\pi\times10^6\,\mathrm{rad/s}\), and
\(|\psi(0)\rangle=|00\rangle\). Times supplied to the solver are seconds.

Identical driving preserves atom-exchange symmetry. Define
\(|S\rangle=(|0r\rangle+|r0\rangle)/\sqrt2=|\Psi^+\rangle\) and
\(|A\rangle=(|0r\rangle-|r0\rangle)/\sqrt2\). The antisymmetric state is
uncoupled and initially empty. In the basis \((|00\rangle,|S\rangle,|rr\rangle)\),

\[
\bar H_{\rm sym}=\begin{pmatrix}0&g&0\\g&0&g\\0&g&V\end{pmatrix},
\qquad g=\frac{\Omega}{\sqrt2}.
\]

The factors of \(\sqrt2\) follow by adding the two matrix elements
\(\Omega/2\) coherently. In the full numerical state, the Bell probability is
computed as \(|(a_{0r}+a_{r0})/\sqrt2|^2\), not as an incoherent sum of
populations. For this symmetric trajectory only, it equals
\(P_{0r}+P_{r0}\).

## Noninteracting limit

For \(V=0\), the atoms evolve independently:
\(\bar H=(\Omega/2)(\sigma_x\otimes I+I\otimes\sigma_x)\).
The two terms commute. A single atom starting in \(|0\rangle\) evolves to
\(c|0\rangle-is|r\rangle\), where
\(c=\cos(\Omega t/2)\), \(s=\sin(\Omega t/2)\). Taking the tensor product gives

\[
|\psi(t)\rangle=c^2|00\rangle-ics(|0r\rangle+|r0\rangle)-s^2|rr\rangle.
\]

Consequently,

\[
P_{00}=c^4,\quad P_{0r}=P_{r0}=c^2s^2,\quad P_{rr}=s^4,
\qquad P_{\Psi^+}=2c^2s^2=\tfrac12\sin^2(\Omega t).
\]

These are exact at zero interaction and are implemented directly in
`noninteracting_state`, independently of the eigensolver. They predict complete
double excitation at 0.5 μs, return to \(|00\rangle\) at 1 μs, and a maximum
Bell probability of 1/2. Small finite interactions are compared to the same
limit to distinguish approximation error from the zero-interaction roundoff.

For these real couplings, populations are even in \(V\): with
\(D=\mathrm{diag}(1,-1,1)\) in the symmetric subspace,
\(\bar H(-V)=-D\bar H(V)D\), so
\(\psi_{-V}(t)=D\psi_V(t)^*\) for this initial state. Thus probability corrections
near \(V=0\) start quadratically, at fixed finite time.

## Strong-blockade limit

Writing the symmetric amplitudes as \(a,b,d\), their equations are

\[
i\dot a=gb,\qquad i\dot b=g(a+d),\qquad i\dot d=gb+Vd.
\]

For \(|V|\gg|\Omega|\), the doubly excited state is far off resonance.
Eliminating it to leading order gives \(d\simeq-gb/V\), and the effective
Hamiltonian in \((|00\rangle,|S\rangle)\) is

\[
\bar H_{\rm eff}\simeq
\begin{pmatrix}0&g\\g&-g^2/V\end{pmatrix}.
\]

In the strict \(|V|\to\infty\) limit the light shift vanishes. Exponentiating
the remaining \(g\sigma_x\) gives

\[
|\psi_\infty(t)\rangle=\cos(gt)|00\rangle-i\sin(gt)|S\rangle,
\]

\[
P_{00}=\cos^2(\Omega t/\sqrt2),\quad
P_{\Psi^+}=\sin^2(\Omega t/\sqrt2),\quad
P_{0r}=P_{r0}=\tfrac12\sin^2(\Omega t/\sqrt2),\quad P_{rr}=0.
\]

These are implemented in `blockade_state`. The collective Rabi frequency is
\(\sqrt2\Omega\); the first ideal Bell maximum is at
\(t=\pi/(\sqrt2\Omega)=0.353553\ldots\) μs. At finite interaction, double
excitation is suppressed rather than identically zero. Its probability is of
order \((\Omega/V)^2\). Light shifts also change the phase and oscillation
frequency. For fixed \(\Omega t\) of order one the plotted population errors
are of order \((\Omega/V)^2\); this is not a uniform long-time error bound.
The leading limit is valid for either sign of \(V\). Amplitude errors need
not have the same scaling as probability errors.

## Reproducible calculation and residuals

Run from the repository root after `python -m pip install -e '.[test]'`:

```bash
python scripts/figure_5_1.py
python -m pytest
```

The two-panel figure uses \(V/(2\pi)=0,0.5,1,5,100\) MHz and 4,001 uniformly
spaced samples over 0–1 μs inclusive. The spacing is 0.25 ns, giving 40 samples
per fastest interaction period (10 ns at 100 MHz). The static solver evaluates
the spectral solution at each time; this spacing controls curve sampling,
not integration accuracy. Both panels use the same colors and line styles.

Residual cases use \(V/(2\pi)=0,0.001,0.01\) MHz against the noninteracting
limit and \(10,30,100,-100\) MHz against ideal blockade. For every one of the
four basis probabilities and the Bell probability, define

\[
r_k=P_k^{\rm numerical}-P_k^{\rm analytical},\qquad
r_{\max}=\max_k|r_k|,\qquad
r_{\rm RMS}=\sqrt{\frac1N\sum_{k=1}^{N}r_k^2}.
\]

These are discrete, equally weighted statistics including both endpoints.
Each residual case is also evaluated on 8,001 points (0.125 ns spacing) to
check that the reported maxima and RMS values are resolved. The residual plots
show the exact \(V=0\) and finite 100 MHz blockade comparisons separately,
since their error scales differ greatly.

Every simulated trajectory, including the refined grids, is checked for:

- \(|P_{00}+P_{0r}+P_{r0}+P_{rr}-1|\le10^{-12}\).
- \(|P_{0r}-P_{r0}|\le10^{-12}\).
- All basis and Bell probabilities in \([-10^{-12},1+10^{-12}]\).

The script raises an error if a check fails. Raw probabilities are neither
clipped nor renormalised, so violations cannot be hidden.

Outputs:

- `results/figures/figure_5_1.pdf` and `.png`: main two-panel figure.
- `results/figures/static_residuals.pdf` and `.png`: residual plots.
- `results/data/figure_5_1.npz`: seconds, physical parameters, complex states,
  all four basis populations, Bell probabilities, analytical states and
  signed residuals. Arrays are indexed by case, time, then basis/observable;
  label arrays and interaction arrays specify their order. Load with
  `np.load(path, allow_pickle=False)`.
- `results/data/static_validation.json`: parameters, library versions, physical
  checks and maximum/RMS residuals on both sampling grids. This is the
  machine-readable numerical report; the script also prints a compact table.

## Measured results for this run

Residuals on the 4,001-point grid (probabilities are dimensionless):

| Limit | V/(2π), MHz | P00 maximum | P00 RMS | PΨ+ maximum | PΨ+ RMS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Noninteracting | 0 | 9.9920e-16 | 3.8206e-16 | 5.5511e-16 | 1.8878e-16 |
| Noninteracting | 0.001 | 1.0876e-6 | 5.0965e-7 | 1.1037e-6 | 4.9930e-7 |
| Noninteracting | 0.01 | 1.0875e-4 | 5.0962e-5 | 1.1036e-4 | 4.9925e-5 |
| Blockade | 10 | 8.0598e-3 | 3.8313e-3 | 1.1892e-2 | 5.6715e-3 |
| Blockade | 30 | 9.0458e-4 | 4.2891e-4 | 1.2991e-3 | 6.3381e-4 |
| Blockade | 100 | 8.1487e-5 | 3.8635e-5 | 1.1551e-4 | 5.7081e-5 |
| Blockade | −100 | 8.1487e-5 | 3.8635e-5 | 1.1551e-4 | 5.7081e-5 |

The exact zero-interaction comparison agrees to floating-point precision.
Reducing a small interaction tenfold reduces the population error roughly
100-fold. Increasing the blockade interaction from 10 to 100 MHz likewise
reduces errors roughly 100-fold. The positive and negative strong interactions
give equal population dynamics within roundoff. These finite-interaction
residuals are physical approximation errors, not evidence of propagation drift.

For the two plotted observables at nonzero interaction, doubling the sample
count changes residual maxima by less than 5e-7 relative and RMS residuals by
less than 6.3e-5 relative. Relative changes at V = 0 are not informative because
those residuals are at the roundoff floor. The saved JSON contains both sets of
metrics for all five observables.

Across the figure, validation and refined-grid trajectories, the maximum
probability-normalisation error is 1.9984e-15 and the maximum exchange-symmetry
error is 5.5511e-16. The raw probability range is
[0, 1.0000000000000013], inside the 1e-12 tolerance. All 169 tests passed,
including the 144-test baseline and 25 new analytical/physical checks.

There is no discrepancy from the derived independent-atom or blockade
behaviour. Exact visual styling and individual red-curve sampling are not claimed. The original Hamiltonian and static solver are unchanged.
No time-dependent propagation or full-milestone completion is claimed.
