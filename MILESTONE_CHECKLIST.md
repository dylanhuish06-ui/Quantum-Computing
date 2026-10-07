# Milestone interview checklist

**Run 6B:** source interpretation resolved. Both figures use representative
V/(2π) = 0, 0.5, 1, 5, 100 MHz samples of the source's 0–100 MHz range;
the booklet does not list the individual red-curve values. Constant φ = 0
is a convention; P00 and PΨ+ are phase invariant from |00>.
The boxes below are personal interview preparation, not unfinished software tasks.

## Equations I can explain

- [ ] Eq. (5.4): basis order, conjugate phase signs, factors of 1/2,
  diagonals Δ and −Δ+V, interaction only on |rr>.
- [ ] Hbar = H/ħ, angular frequency versus Hz, MHz and microsecond conversions.
- [ ] Eq. (5.6): R exp(−i diag(λ)t) R†ψ0; eigenvectors are columns.
- [ ] Eq. (5.7): exp(−i Hbar(t_n)dt), evaluated at the left endpoint.
- [ ] |Ψ+> = (|0r>+|r0>)/√2; fidelity is a coherent projection.
- [ ] Explain H(φ) = D H(0) D† and why the milestone probabilities are invariant.
- [ ] Distinguish the full 100 MHz convergence criterion from family plot-stability checks.
- [ ] Derive independent-atom and blockade probability curves.
- [ ] Cubic and linear detuning endpoints, slopes and resonance crossing.

## Methods and checks I can explain

- [ ] `eigh` for a Hermitian Hamiltonian; one static diagonalisation.
- [ ] Full-matrix `expm` and chronological order for noncommuting steps.
- [ ] Static t = 0 convention versus the timestep solver's first-grid-time state.
- [ ] Independent `expm`, analytical limits and DOP853 reference checks.
- [ ] Norm, all-four-population sum, bounds and P0r = Pr0 symmetry checks.
- [ ] No hidden state renormalisation; raw probabilities used in audits.
- [ ] Four-significant-place criterion, scalar rounding and reference refinement.

## Plots I can reproduce

- [ ] Figure 5.1 and limiting-case residuals: `python scripts/figure_5_1.py`.
- [ ] Figure 5.2 with source parameters and documented representative V sampling.
- [ ] Linear/cubic comparison and timestep-convergence figure (same sweep command).
- [ ] Find data, full metrics and exact commands in `docs/milestone_report.md`.

## Likely numerical-error questions

- [ ] Why can norm conservation coexist with an inaccurate trajectory?
- [ ] Why does left-endpoint error approximately halve when dt halves?
- [ ] Why need constant-H error not decrease with smaller dt?
- [ ] Which residual is finite-blockade physics and which is numerical roundoff?
- [ ] How was the fine reference independently checked?
- [ ] Why do four significant figures differ from four decimal places?
- [ ] What happens to relative error near zero probability?
- [ ] What limits peak-time accuracy and the claim “largest timestep”?
- [ ] Why can reaching Δ ≈ V change Bell fidelity?

Before review: run all tests, regenerate the documented milestone cases, inspect warnings
and `git status --short`. Do not stage, commit, push or tag before explicit approval.
