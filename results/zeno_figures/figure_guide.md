# Quantum Zeno figure guide

This directory is a figure library, not a poster or presentation. There are nine numbered figures and ten export sets (Figure 03 also has a single-column version). All numbered figures are exported as vector PDF, SVG with live text, and 400-DPI PNG. The contact sheet is only a selection aid and uses raster thumbnails.

## Reproduction

```bash
python -m pip install -r results/zeno_figures/requirements.txt -e '.[test]'
python scripts/generate_zeno_figures.py --all
python -m pytest tests/test_zeno_figures.py
```

Run from a checkout with the committed inputs. `--output-dir /tmp/zeno-reproduction` provides a separate output directory; `--validate-only` audits archived inputs without generating graphics. The manual GitHub Actions workflow regenerates to its temporary directory and uploads an artifact. No remote publishing, large-N simulation, or historical-data overwrite occurs. NumPy, SciPy, Matplotlib and Pillow (a Matplotlib dependency) are used; exact versions are recorded in `reproducibility_report.json`. Same-environment exports have fixed metadata and stable SVG IDs; different library/font versions can alter binary rendering. Numerical reproducibility is checked with stated absolute tolerances.

## Model and quantities

Basis: |00>, |0r>, |r0>, |rr>. In the repository convention, Hbar = H/hbar is

```
[[Delta, c, c, 0],
 [c*,    0, 0, c],
 [c*,    0, 0, c],
 [0,    c*, c*, V-Delta]],       c = (Omega/2) exp(-i phi).
```

Omega, Delta and V are angular frequencies in rad/s; U(t)=exp(-i Hbar t), with t in seconds. All plotted experiments use Omega/(2pi)=1 MHz, Delta=phi=0 and initial |00>. Unless explicitly varied, T=0.5 us. The exchange-antisymmetric state is dark; the remaining basis is |00>, |Psi+>=(|0r>+|r0>)/sqrt(2), |rr>. The coherent Bell population must not be mistaken for an incoherent sum for arbitrary states, although they coincide in this symmetric study.

For M>0, project at kT/M for k=1,...,M, including T; M=0 is unmeasured. The ground binary question uses P=diag(1,0,0,0), Q=I-P. Its non-selective channel is D(rho)=P rho P+Q rho Q. It preserves both diagonal blocks, including coherences inside Q. Switching the names/order of P and Q changes no non-selective dynamics. The double-excitation question instead uses R=|rr><rr| and I-R.

- **Non-selective final population:** p_P(T)=Tr(P rho_M), summed over all measurement histories, including departures followed by returns.
- **All-success probability:** S_M=|<00|U(T/M)|00>|^(2M), the joint probability of P at every event. This expression applies to the rank-one initial-state projector.
- **Post-selected state:** after a successful rank-one projection the normalized state is |00>, with conditional ground population exactly one. This is not S_M.
- **Subspace population:** for the double question the initial sector is rank-three I-R, with population 1-Prr. This is neither P00 nor the probability of all successful I-R outcomes.

From the actual Hamiltonian, Var(Hbar)=Omega^2/2 in |00>, independent of Delta, phi and V. Thus tau_Z=hbar/sqrt(Var(H))=sqrt(2)/|Omega|=225.079079 ns. For this static finite-dimensional Hermitian model, P(t)=1-(t/tau_Z)^2+O(t^4), stronger than the general requested O(t^3) remainder, because survival is even in t. No universal finite-time validity window follows from tau_Z alone; higher moments matter. Figure 04 calculates loss from the outgoing amplitudes to avoid cancellation in 1-P at very short times.

At fixed T and Hamiltonian, S_M=[1-a/M^2+O(M^-4)]^M with a=(T/tau_Z)^2, so log S_M=-a/M+O(M^-3). Both 1-S_M and 1-p_P have leading term a/M, but the non-selective population has return contributions at O(M^-2). The exponential exp(-a/M) is a large-M approximation, not an exact finite-M non-selective formula. Here a=pi^2/2. It is drawn only for M>=20 at V=0. The short-interval expansion must hold before applying it; a large interaction introduces additional spectral scales.

Interaction blockade detunes |rr> at large V, but permits collective motion between |00> and |Psi+>. Measuring the ground sector suppresses even this excitation. Measuring the double-excitation sector at high M also leaves internal motion: P00 tends to cos^2(Omega T/sqrt(2)) at Delta=0. These ideal numerical projections are not ancilla experiments or finite-strength continuous measurements. No statistical error bars are assigned to deterministic calculations.

For context on the distinction between freezing and evolution inside a multidimensional Zeno subspace, see P. Facchi and S. Pascazio, *Quantum Zeno dynamics: mathematical and physical aspects*, J. Phys. A 41 (2008) 493001, [author manuscript](https://arxiv.org/abs/0903.3297), [DOI](https://doi.org/10.1088/1751-8113/41/49/493001). The numerical results and asymptotic coefficients above were computed directly from this repository's model.

## Research inventory and preservation

The tracked Hamiltonian is `src/quantum_project/hamiltonian.py`; the milestone spectral state solver is `src/quantum_project/static_solver.py`. Existing Zeno implementations are `src/quantum_project/density_matrix.py` and `src/quantum_project/quantum_zeno.py`; studies and audit live in `scripts/quantum_zeno.py` and `scripts/audit_quantum_zeno.py`, with `tests/test_quantum_zeno.py`. These Zeno files and saved datasets were untracked at the start of this figure task. They were read and left unchanged.

`inputs/zeno_study.npz` is a byte-identical snapshot of the pre-existing `results/data/zeno_study.npz`, not a modified historical dataset. `inputs/provenance.json` records its original path and SHA-256, alongside text snapshots of the original generator, density/measurement routines, tests and audit documents. This makes the figure commit self-contained without staging the user's existing research changes. The independent figure generator imports only the tracked Hamiltonian; all saved complex matrices are revalidated against spectral unitaries and 16x16 channel powers. Plotting reads validated archived values rather than rerunning long M=10000 loops. New heatmaps, interaction controls, coherence examples and convergence diagnostics are small four-dimensional calculations saved in `derived_data.npz`.

The separate `symmetry-scaling/` worktree has a Dicke reduction, full-versus-reduced agreement data and nonuniform-interaction leakage data (`symmetry_agreement.npz`, `symmetry_broken.npz`, `symmetry_validation.json`). Inspection found no projection-protected large-N trajectories or validated frequency experiment there. Optional Figure 10 is therefore omitted. Symmetry leakage alone is not evidence of measurement protection. No large-N computation was run.

## Validation and limits

The report records full-complex-matrix reference errors, trace/Hermiticity/positivity, M=0 and M=1, exact projector properties, outcome-reversal equivalence, all-success ordering, short-time predictions, and large-M leading loss. Saved-array checks use a maximum absolute tolerance of 3e-11, allowing accumulated spectral roundoff at M=10000; no densities are renormalized, clipped or projected onto a positive cone. New RK4 convergence is checked against a separate spectral reference. Production propagation is exact for a static Hamiltonian up to floating-point error, so it has no integration-resolution bias to extrapolate. The manifest ties each graphic to specific dataset keys and hashes; file checks verify PNG dimensions/DPI, PDF signatures, parseable SVG with text and no embedded raster elements. Human visual inspection is recorded separately in `quality_control.md` after rendering, not claimed by the automated script.

Palette: navy ground population/all-success (context-specific), teal Bell population/non-selective protection, amber double-excitation/asymptotic reference, muted red alternative interaction/projector, grey unmeasured controls. State identities are consistent across population panels; comparison figures distinguish protocols by colour and line/marker style. Legends explicitly identify each quantity. White backgrounds, STIX mathematical glyphs, vector linework, and a sequential non-rainbow map support printing. Wide figures are 8.4 inches across; the Figure 03 single variant is 4.2 inches. PDFs retain embedded fonts; SVG live text requires compatible DejaVu/STIX fonts in the downstream editor.

## Figure-by-figure interpretation

### 01 — Measurement protocol (wide)

**Physical meaning:** Coherent intervals followed by the complete binary Lüders channel.

**Parameters:** Illustrative M=0,4,12; tau=T/M; general T.

**Main finding:** Increasing M shortens each coherent interval without altering the initial state.

**Axes:** Schematic time increases left to right; event spacings are uniform.

**Limitations:** Ideal mathematical projections, no apparatus, ancilla, finite pulse width or detection inefficiency.

**Suggested poster use:** Introduce the measurement model; pair with Figure 03.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** Vector schematic derived from the stated protocol; no numerical dataset.

**Exports:** [PDF](pdf/fig01_measurement_protocol.pdf), [SVG](svg/fig01_measurement_protocol.svg), [PNG](png/fig01_measurement_protocol.png).

### 02 — Population evolution (wide)

**Physical meaning:** Resolved populations in the symmetric three-state sector.

**Parameters:** Common defaults; V=0; M=0,2,10,100; 1001 time samples.

**Main finding:** Measurements of the ground-state sector inhibit excitation even at V=0, where blockade is absent.

**Axes:** Time in microseconds; populations dimensionless, common scales in all panels.

**Limitations:** The Bell population is a coherent projector expectation, not one atom’s excitation probability.

**Suggested poster use:** Main dynamics panel.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `archive:times_s`, `archive:trajectory_probabilities[0,0:4]`

**Exports:** [PDF](pdf/fig02_population_evolution.pdf), [SVG](svg/fig02_population_evolution.svg), [PNG](png/fig02_population_evolution.png).

### 03 — Zeno survival curve (wide)

**Physical meaning:** Final ground population and the joint probability of every ground outcome.

**Parameters:** Common defaults; V=0; counts 1 to 10000.

**Main finding:** Both losses approach (T/tau_Z)^2/M; return histories make non-selective population larger at finite M.

**Axes:** Logarithmic measurement count; right panel shows logarithmic loss to resolve near-unity probabilities.

**Limitations:** Exponential shown only at M>=20; asymptotic, not an exact formula for either finite-M curve. M=0 stated separately.

**Suggested poster use:** Flagship result; single-column version for narrow spaces.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `archive:measurement_counts`, `archive:final_probabilities[0,0]`, `archive:exact_all_yes[0]`

**Exports:** [PDF](pdf/fig03_zeno_survival_vs_measurements.pdf), [SVG](svg/fig03_zeno_survival_vs_measurements.svg), [PNG](png/fig03_zeno_survival_vs_measurements.png).

### 03 — Zeno survival curve (single)

**Physical meaning:** Final ground population and the joint probability of every ground outcome.

**Parameters:** Common defaults; V=0; counts 1 to 10000.

**Main finding:** Both losses approach (T/tau_Z)^2/M; return histories make non-selective population larger at finite M.

**Axes:** Logarithmic measurement count; right panel shows logarithmic loss to resolve near-unity probabilities.

**Limitations:** Exponential shown only at M>=20; asymptotic, not an exact formula for either finite-M curve. M=0 stated separately.

**Suggested poster use:** Flagship result; single-column version for narrow spaces.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `archive:measurement_counts`, `archive:final_probabilities[0,0]`, `archive:exact_all_yes[0]`

**Exports:** [PDF](pdf/fig03_zeno_survival_vs_measurements_single.pdf), [SVG](svg/fig03_zeno_survival_vs_measurements_single.svg), [PNG](png/fig03_zeno_survival_vs_measurements_single.png).

### 04 — Short-time quadratic law (wide)

**Physical meaning:** Survival loss of the initial pure state during a single unmeasured interval.

**Parameters:** Omega/(2pi)=1 MHz; Delta=phi=0; V/Omega=0,100; t=0.01 to 100 ns; initial |00>.

**Main finding:** The same leading coefficient holds at V/Omega=0 and 100; finite-time corrections need not be the same.

**Axes:** Nanoseconds on a logarithmic x axis; log-log loss at left, loss divided by the leading quadratic term at right.

**Limitations:** Static finite-dimensional Hermitian H makes P(t) even: next correction is O(t^4). tau_Z alone does not fix the useful expansion window.

**Suggested poster use:** Explain the mechanism beneath the main survival curve.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `derived:short_times_s`, `derived:short_loss`, `archive:variance_rad2_s2`

**Exports:** [PDF](pdf/fig04_short_time_quadratic_law.pdf), [SVG](svg/fig04_short_time_quadratic_law.svg), [PNG](png/fig04_short_time_quadratic_law.png).

### 05 — Zeno regime heatmap (wide)

**Physical meaning:** Non-selective ground-state population over total time and measurement count.

**Parameters:** Omega/(2pi)=1 MHz; Delta=phi=V=0; T=0.01 to 1 us; integer M=1 to 1000; initial |00>.

**Main finding:** Sustained protection across the time window requires more measurements as T increases.

**Axes:** Microseconds horizontally, logarithmic M vertically, colour gives final ground population. Amber contour marks 0.95.

**Limitations:** M is a count, not a rate: rate=M/T. Low-M revivals can also yield large final populations; see the matched M=0 panel.

**Suggested poster use:** Show the range of useful parameters.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `derived:heat_times_s`, `derived:heat_counts`, `derived:heat_P00`, `derived:heat_unmeasured_P00`

**Exports:** [PDF](pdf/fig05_zeno_regime_heatmap.pdf), [SVG](svg/fig05_zeno_regime_heatmap.svg), [PNG](png/fig05_zeno_regime_heatmap.png).

### 06 — Interaction dependence (wide)

**Physical meaning:** Matched measured and unmeasured endpoints as interaction strength changes.

**Parameters:** Common defaults; V/Omega=0 and 0.01 to 100; M=0,10,1000.

**Main finding:** Blockade suppresses |rr> while allowing collective single excitation; ground-sector measurements inhibit excitation even at V=0.

**Axes:** V/Omega dimensionless with a linear segment from 0 to 0.01 then logarithmic; populations at fixed T.

**Limitations:** Endpoint oscillations with V are coherent dynamics, not error bars or a monotonic blockade threshold.

**Suggested poster use:** Distinguish the two suppression mechanisms.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `derived:interaction_ratios`, `derived:interaction_counts`, `derived:interaction_populations`

**Exports:** [PDF](pdf/fig06_interaction_dependence.pdf), [SVG](svg/fig06_interaction_dependence.svg), [PNG](png/fig06_interaction_dependence.png).

### 07 — Density-matrix coherence (wide)

**Physical meaning:** Magnitudes of a physical evolved density matrix before/after the binary Lüders measurement, and magnitudes of removed entries.

**Parameters:** Omega/(2pi)=1 MHz; Delta=phi=V=0; pure |00> evolved for 0.15 us, then one binary measurement.

**Main finding:** Ground-to-excited coherences vanish, while all populations and within-Q coherences survive.

**Axes:** Row/column labels are basis indices |00>, |0r>, |r0>, |rr>; all three matrices use one colour scale.

**Limitations:** Magnitude plots do not show phases; complete complex matrices are saved. Not a measurement resolving all four basis states.

**Suggested poster use:** Explain the measurement superoperator.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `derived:coherence_before`, `derived:coherence_after`

**Exports:** [PDF](pdf/fig07_density_matrix_coherence.pdf), [SVG](svg/fig07_density_matrix_coherence.svg), [PNG](png/fig07_density_matrix_coherence.png).

### 08 — Numerical convergence (wide)

**Physical meaning:** Actual full-density Frobenius errors relative to an independently diagonalised reference.

**Parameters:** Common defaults; V/Omega=1; M=8; n=1,2,4,8,16,32,64; spectral reference.

**Main finding:** RK4 approaches fourth order; subdividing the exact exponential only changes accumulated floating-point roundoff.

**Axes:** Substeps per coherent interval; dimensionless Frobenius error. Left log-log, right log x and linear error divided by 1e-14.

**Limitations:** RK4 is an additional diagnostic, not the solver used to generate the production figures. Curves are numerical errors, not statistical uncertainties.

**Suggested poster use:** Validation inset or methods panel.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `derived:convergence_steps`, `derived:convergence_rk4_error`, `derived:convergence_expm_error`

**Exports:** [PDF](pdf/fig08_numerical_convergence.pdf), [SVG](svg/fig08_numerical_convergence.svg), [PNG](png/fig08_numerical_convergence.png).

### 09 — Measurement-subspace comparison (wide)

**Physical meaning:** Two genuine binary projector choices, with a common ground observable and the population of the sector containing the initial state.

**Parameters:** Common defaults; V=0; ground vs double-excitation binary projectors; M=1 to 10000.

**Main finding:** Protecting a multidimensional sector permits motion inside it. Avoiding |rr> does not freeze the initial ground state.

**Axes:** Log M; left ground population for both channels, right rank-one P00 or rank-three I-Prr sector population.

**Limitations:** The original any_excitation question only reverses the ground question’s outcomes; it is not an independent dynamics. Rank-three occupancy is not all-success probability.

**Suggested poster use:** Explain Zeno dynamics beyond rank-one freezing.

**Source script:** `scripts/generate_zeno_figures.py`. **Dataset keys:** `archive:measurement_counts`, `archive:final_probabilities[0,0]`, `archive:final_probabilities[0,2]`

**Exports:** [PDF](pdf/fig09_measurement_subspace_comparison.pdf), [SVG](svg/fig09_measurement_subspace_comparison.svg), [PNG](png/fig09_measurement_subspace_comparison.png).

Dataset prefix `archive:` refers to `inputs/zeno_study.npz`; `derived:` refers to `derived_data.npz`. Python-style array indices identify the exact subset. Figure numbering intentionally stops at 09 because the optional large-N protection study is unavailable.
