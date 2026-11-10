# Phase 4 — Resumable Ancilla-Assisted Many-Qubit Programme

This is the canonical, session-independent master specification. Any Codex session must read this file **and** [STATUS.md](STATUS.md) before writing code. Resume using the one-line instruction in the user-facing response; do not assume prior chat context survives.

## Mission and scientific boundaries

Develop, validate, document, and benchmark a research-quality extension of the Chapter 5 Rydberg model that explores (1) symmetry-sector diagnostics and ancilla measurement, (2) selective Zeno stabilisation, (3) total-spin energy-gap and optional dynamical protection, (4) exact *hybrid Dicke + distinguished-qubit* time-dependent evolution for structured inhomogeneity, and (5) *controlled, explicitly approximate* adaptive treatments of more general asymmetry. Use the existing symmetric/full solver as the starting point; do **not** modify milestone physics, figures or poster scripts without an explicit separate request.

**Science guardrails:**
- A temporary ancilla used for syndrome extraction is **not** interchangeable with a distinguished physical qubit retained explicitly in a hybrid model. Clearly separate these approaches and resource costs.
- Postselection onto the symmetric sector changes generic symmetry-broken dynamics. Distinguish protected target dynamics from exact evolution under the original asymmetric Hamiltonian. Never promise exact polynomial simulation of arbitrary interacting qubits.
- An ideal global controlled-projector oracle is a mathematical proof of concept, not a constant-depth local gate implementation; calculate or explicitly leave unresolved its compilation cost.
- A hybrid Dicke reduction is exact **only if** permutations of the bulk atoms commute with the full time-dependent Hamiltonian and the initial bulk state lies in that symmetric subspace. Uniform bulk drives, detunings, bulk–bulk pair interactions, and a site-independent bulk coupling to each distinguished qubit suffice. Generic pair disorder or local bulk detuning disorder invalidates exactness.
- Historical prototype values mentioned in chat (N=4, N=6, N=1024) are **not validated repository results**. Recreate their numerical methods independently, check convergence and store actual reproducible outputs. Do not copy previously quoted numbers as measured findings.
- Use \(\bar H=H/\hbar\) with frequencies in rad/s and time in seconds, consistent with the milestone. Model in a named, documented dimensionless convention only if conversions are explicit and tested. Preserve \(\Delta(N/2-k)\) constant energy offsets when comparing complex amplitudes.
- No long-running autonomous/background jobs: each invocation works within the active session. Resume must rely on repository files/commits, not a chat summary.

## Mathematical model

\[
\bar H(t)=\sum_i\left[
\Delta_i(t)(\tfrac12-n_i)+
\tfrac{\Omega_i(t)}2(e^{i\phi_i(t)}\sigma_i^++e^{-i\phi_i(t)}\sigma_i^-)
\right]+\sum_{i<j}V_{ij}(t)n_i n_j .
\]

For \(N\) atoms, define \(|D_N^k\rangle\), the isometry \(B\) whose columns are Dicke states, \(P=BB^\dagger\), \(Q=I-P\), and \(J_{\max}=N/2\).

- Explicit full-space \(P\), \(Q\) are allowed for small-N validation, not for large N.
- Quantify \(\|Q\psi(t)\|^2\), \(\|Q\bar H(t)P\|\), conditional/unconditional target probabilities, and full-state error where applicable.
- Ideal ancilla syndrome unitary: \(U_{\mathrm{syn}}=P\otimes I_a+Q\otimes X_a\). Check unitarity and Born probabilities. Ancilla \(|0\rangle\) postselection returns \(P|\psi\rangle\) with its correct norm/probability; ancilla \(|1\rangle\) reports leakage.
- Compare selective Zeno measurements at varying interval and fixed numerical integration accuracy. Report both conditional fidelity and probability of obtaining all successful syndromes; count syndrome measurements and any reset requirements.
- Gap protection: \(\bar H_G(t)=\bar H(t)+G[J_{\max}(J_{\max}+1)I-J^2]\). The first total-spin-sector gap is \(NG\), but this all-to-all penalty is not hardware-free. Compare at fixed disorder, timings and target dynamics.
- Hybrid exact basis: \(\{|D_{N-m}^{k}\rangle\otimes|a_1\cdots a_m\rangle\}\), dimension \((N-m+1)2^m\). Implement sparse direct construction without creating a full \(2^N\) array except for small-N validation. Test \(\bar H_{\mathrm{hybrid}}=B^\dagger\bar H_{\mathrm{full}}B\) **and** invariance \(\bar H_{\mathrm{full}}B=B\bar H_{\mathrm{hybrid}}\) for admissible models. Validate observable trajectories and complex phases.
- For mean-zero static detuning imperfections \(D=\sum_i\epsilon_i n_i\), verify for Dicke states the small-N identity \(\|QD|D_N^k\rangle\|^2=[k(N-k)/(N-1)]\sigma_\epsilon^2\) for \(N>1\). Distinguish analytical upper bounds, short-time perturbation series and measured errors.
- Generic asymmetry: study low-rank/spatial-mode or adaptive-basis approximation only with residual estimates, convergence, failure examples and a full-space small-N reference. No claim of exactness without proof.

## Repository hygiene and continuation contract

- Work only on `extension/ancilla-hybrid-phase4` (or a same-purpose worktree rooted in it). It was forked from `extension/symmetry-scaling`. Never alter `main`, merge unrelated `extension/quantum-zeno`, or overwrite existing large-N results.
- At the **start of every invocation** inspect current branch, dirty state, recent commits, `docs/ancilla_phase4/MASTER_PLAN.md` and `docs/ancilla_phase4/STATUS.md`. Compare documented stage state to the actual files/tests, and reconcile discrepancies. If the branch is unavailable locally, fetch and check it out safely. If unrelated changes are present, preserve them and work separately.
- Stages are a directed dependency graph, not a demand to fit everything into one session. Choose the **earliest incomplete or failing gate** and work on it to a tangible tested milestone. Advance to the next gate only when the preceding exit tests pass; you may complete multiple stages in one invocation when feasible.
- Before ending **each invocation**, even partially through a stage: save deliverables, update STATUS.md with last verified commit, exact tests run, what passed/failed, parameter choices, artifact paths, known limitations and the *one immediate next action*. Commit coherent changes with descriptive messages and push the branch when authorised and available. Do not claim successful push if it failed. Keep generated large arrays out of Git unless size-appropriate; manifest hashes and recipes for reproducible local files.
- If tools/runtime/usage limits interrupt, prioritise a valid checkpoint; if a run cannot be verified, mark it `in_progress` or `blocked`, never `passed`. Repeat the identical one-line instruction in a new Codex session to continue.
- Scientific CPU timing must be labelled provisional in Codespaces and CI. No misleading 'wall-clock optimal' benchmarks. Respect resource guards and maximum runtime/memory budgets; do not launch high-N full-state expansions.

## Sequential gates

### Gate 0 — Baseline and audit
Read docs/symmetry_scaling_extension.md and underlying solver, units, tests, available scripts and performance safeguards. Record current base commit/branch, Python dependencies and baseline test counts in STATUS.md. Establish `src/quantum_project/ancilla/`, `tests/ancilla/`, `scripts/ancilla/`, `docs/ancilla_phase4/` and `results/ancilla/` only as needed. Preserve public interfaces. Exit: baseline tests pass and traceable plan.

### Gate 1 — Formal operators and independent time-dependent reference
Implement well-typed, tested small-N projector/isometry and total-spin-sector utilities; general time-dependent Hamiltonian with arbitrary local detunings/drives and pair couplings, matching existing signs and units. Build independent DOP853 and midpoint-exponential references. Tests: Hermiticity, norm conservation, known N=2 milestone Hamiltonian, basis ordering, `B†B=I`, `P²=P`, commutator/invariance, exact constant-H limit, solver grid refinement, complex phases. Exit: independent cross-solver agreement with measured tolerances.

### Gate 2 — Ancilla-syndrome and Zeno
Implement explicit ideal controlled-sector oracle for small N; verify unitary, branch norms, Born weights and conditional states, physically meaningful reset/report logic. Implement repeated measurements with deterministic conditional simulation and success product, independent of ODE/time-discretisation step sizes. Compare unprotected evolution, ideal projected dynamics and selective Zeno at N=2–6 with asymmetric time-dependent examples. Include exact success-vs-fidelity curves and tests for the limits `Pψ=ψ` and `Pψ=0`. Exit: reference-verified numerical results, honest oracle overhead caveat.

### Gate 3 — Gap/dynamical protection
Add spin-gap Hamiltonian and optional separate finite-rate dynamical symmetrisation. Derive the penalty spectrum and check its vanishing action on Dicke states. Separate coherent protection, projection and no-protection metrics. For dynamical averaging, check finite-cycle/Trotter errors and pulse costs; if unavailable, document and skip rather than fake results. Exit: unit tests and reproducible comparisons.

### Gate 4 — Exact hybrid Dicke + distinguished-register solver
Build sparse reduced generators with arbitrary `m` distinguished qubits (bounded in tests), time-dependent common bulk terms and site-independent bulk-to-each-defect interaction. Distinguish total N from bulk M=N-m; retain exact diagonal offsets. Validate at N=2..8 for m=0,1,2 (and m=3 if feasible), complex phases and nonsymmetric distinguished configurations. Compare to independent full \(2^N\) evolution in the physical sector, including invariance, populations, correlations, phase-sensitive amplitudes and extrapolated convergence. Test rejecting bulk-asymmetric inputs or clearly labelling them approximate. Exit: verified exactness under explicit assumptions.

### Gate 5 — Sparse large-N study and adaptive approximation
Time evolve selected admissible N=64,256,1024 (or smaller if safeguards bind) using sparse exponentiation / Krylov on the hybrid state, including a time-dependent sweep, excitation fractions and distinguished-site populations. Refine time step, Krylov tolerance, and output grid independently; check norm, finite states, probability bounds, reproducibility and raw report completeness. Implement one controlled approximate asymmetric correction model with residual monitoring and small-N full-space reference, rather than asserting arbitrary-system efficiency. Exit: converged trajectories and transparent failure regime.

### Gate 6 — Figures, comparison and reproducibility
Produce publication-ready PDF/PNG figures (legible axis units, colorblind-safe differentiated strokes, physically grounded legends, no overly dense panels) and CSV/NPZ + machine-readable JSON of parameters, environment, tolerances, git SHA, raw residuals and convergence. Provide an errors/statistics table and report that separates:
1. exact symmetry,
2. exact hybrid structured asymmetry,
3. postselected protected dynamics,
4. approximate adaptive dynamics.
Benchmark on same computer only, no unsupported speed-up claims. Clean documentation and one stable user CLI such as `python run_project.py ancilla --resume` **if** the existing launcher design supports it; otherwise add a documented `python -m quantum_project.ancilla.run --resume` entrypoint. A repeated CLI invocation must skip validated artifacts via parameter hashes and resume only incomplete experiments. Exit: all tests plus end-to-end smoke reproducibility checks.

### Gate 7 — Final review / hand-off
Re-run full tests, lint/type checks when configured, smoke runs, mathematical/physics audits, artifact checksums and documentation links. Finish a complete explanatory report with assumptions, error bars only where meaningful, scaling limits, quantum-gate feasibility and future experiments. Mark STATUS.md `ready_for_review` only if all mandatory gates pass. Do not merge to main without explicit user approval.

## Required checkpoint schema

STATUS.md should always have: branch and current head commit; last working tree check; per-gate state (`not_started`, `in_progress`, `passed`, `blocked`); evidence of passing tests and exact commands; last verified numerical residuals with units/meaning; reproducible result and figure paths; known problems and unsupported claims; immediate next action. Use factual UTC dates and exact SHA only when actually obtained.

## Execution output per invocation

Finish with a concise report: gate(s) addressed, files changed, tests run and their outcomes, commit/push status and link, scientific findings versus preliminary hypotheses, blockers, **one recommended next action**. Never substitute a motivational description for actual tool results.
