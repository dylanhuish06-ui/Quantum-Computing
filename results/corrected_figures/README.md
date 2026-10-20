# Corrected figures — 51 interaction strengths

The original figures created in ChatGPT on 20 October 2026 were never
committed to GitHub or saved into this Codespace. This directory gives them a
**stable repository location after local generation**. It is separate from
`results/figures/` and does not overwrite the validated five-curve milestone.

## Generate in a Codespace

From the repository root:

```bash
python -m pip install -e '.[test]'
python scripts/generate_corrected_figures.py
```

Generated files (not precommitted):

- `results/corrected_figures/rydberg_figure_5_1_51_interactions.png`
- `results/corrected_figures/rydberg_figure_5_1_51_interactions.svg`
- `results/corrected_figures/rydberg_figure_5_2_51_interactions.png`
- `results/corrected_figures/rydberg_figure_5_2_51_interactions.svg`
- `results/corrected_figures/diagnostics.json`

Use `git add results/corrected_figures && git commit` and push after generation
to persist the PNG and SVG exports in the remote repository. Until then, they
exist only in your current Codespace or downloaded chat files.

## Scientific settings

- Two atoms, basis `|00>, |0r>, |r0>, |rr>`; original Chapter 5 Hamiltonian
  in angular-frequency units.
- `Omega/(2pi) = 1 MHz`, `V/(2pi) = 0, 1, ..., 50 MHz` (51 curves).
- Figure 5.1: static detuning `Delta = 0`, analytical eigenpropagation.
- Figure 5.2: cubic sweep `Delta/(2pi) = 50 MHz * (2t/T - 1)^3`,
  `T = 1 microsecond`; DOP853 with tight tolerances.
- `P_Psi+ = |(c_0r + c_r0)/sqrt(2)|^2`, the **coherent Bell-state
  projection**, not the sum of the independent populations.
- No probability clipping or post hoc normalization. Generated outputs are
  checked for state-norm drift and exchange symmetry.

The higher-detuning Bell population remains nonzero: reference values from
the 20 October re-creation at `T=1 microsecond` are about 0.3035 (`V=0`),
0.8135 (`V=10 MHz`), 0.8959 (`V=20 MHz`) and 0.9821 (`V=50 MHz`).
The generated JSON provides all checked final values.

These 51-curve plots are additional sampling, not a claim that the original
project brief specified every one of those interaction values.
