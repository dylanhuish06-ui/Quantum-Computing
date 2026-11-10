# Reusable Codex command

Paste exactly this instruction at the beginning of **each** Codex session, without editing:

```text
RESUME PHASE 4 — ANCILLA/HYBRID MANY-QUBIT SIMULATION. In dylanhuish06-ui/Quantum-Computing, check out extension/ancilla-hybrid-phase4 safely; read docs/ancilla_phase4/MASTER_PLAN.md and docs/ancilla_phase4/STATUS.md, inspect the actual repository state and recent commits, then execute the earliest unfinished or failing validation gate. Continue through later gates only when preceding tests pass and session resources allow. Implement physics and numerical tests, run independently validated experiments with convergence and resource safeguards, and save code, reproducible data, documentation, figures and a verified tests/errors report. Never overwrite unrelated work or change main, never conflate postselected protection with exact asymmetric dynamics, and never claim untested numerical results or background execution. Before this session ends, update STATUS.md with completed/in-progress gates, verified results, exact checks, commit SHA, blockers and one immediate next action; commit and push completed changes to the research branch. If interruption prevents completion, checkpoint partial work honestly. Finish with a concise summary and stop. On the next invocation of this SAME instruction, resume from STATUS.md rather than repeating completed work.
```

This command is intentionally *idempotent with respect to completed and verified gates*. It is not a command that runs unattended between invocations. A later invocation resumes by inspecting committed work, so push/checkpoint success is required for cross-machine continuation.
