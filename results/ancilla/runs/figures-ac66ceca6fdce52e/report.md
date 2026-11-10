# Phase 4 comparison and reproducibility

All experiments were recomputed in this session environment. Raw CSV/NPZ data,
parameter recipes, software versions, local Git provenance and content hashes
are retained in the experiment directories linked by `latest.json`.

| Class | Validation / result | Interpretation |
| --- | --- | --- |
| Exact symmetry | Maximum full-reference L2 difference 2.576e-12 | Bulk symmetry and symmetric initial bulk required |
| Exact hybrid structured asymmetry | Maximum small-N reference difference 4.705e-12 | Distinguished physical atoms; invariant bulk |
| Postselected protected dynamics | Maximum ODE refinement difference 2.516e-14 | Physical target is projected dynamics; success probability is reported separately |
| Approximate adaptive dynamics | N=5 tight budget reaches dimension 32 | Static full-space diagnostic; generic compression can fail |

The large-N model uses N=64,256,1024, one distinguished physical atom and
tau=Omega_ref*t with Omega_ref=2*pi*1 MHz. At N=1024 the finest time-step error
estimate is 6.586e-05; norm drift is
2.842e-14. The time error is an asymptotic Richardson
estimate, not a rigorous error bar. Exponential action, time step and output
sampling are refined independently. Initial states and complete coefficient
profiles are specified by the content-hashed validation recipes.

`errors.csv` separates numerical residuals, physical leakage, postselection
probabilities and approximation error. These quantities must not be ranked as
if they measured the same task. They are deterministic; no statistical error
bars or experimental noise estimates are inferred. `comparison.json` retains
the full machine-readable table and provenance.

`accuracy.pdf/png` compares exact reductions with independent full-space
evolution and shows large-N second-order convergence. Exactly zero small-N
residuals are displayed at a labelled 1e-16 plotting floor; raw values are
retained unchanged in the error table. `protection.pdf/png`
compares ideal selective checks and coherent spin-gap protection at the same
N=6 disorder, duration and initial state. Conditional fidelity, all-check
success and successful-target probability are different observables.
`large_n.pdf/png` shows bulk and distinguished-site populations; almost
overlapping curves are a measured feature of this short sweep.
`adaptive.pdf/png` shows the conservative static Duhamel residual bound and
actual model error; neither bound nor error must decrease monotonically with
each promoted atom. Reaching dimension 32 eliminates compression at N=5.

The ideal syndrome oracle is nonlocal; its gate compilation depth and hardware
cost remain unresolved. The coherent total-spin penalty is not hardware-free.
No local hardware implementation, generic exact polynomial simulator or
speedup benchmark is claimed. CPU timings are provisional runtime records.

Repeat `PYTHONPATH=src python -m quantum_project.ancilla.run --resume` from the
repository. Only matching parameter/source/environment hashes and intact
artifact checksums are skipped. Invalid outputs are preserved with an invalid
suffix and rebuilt; interrupted pending folders remain available for inspection.
This is an experiment runner, not an autonomous code-generation command.
