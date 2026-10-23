# Phase 4 — Persistent Progress Ledger

**Canonical plan:** [MASTER_PLAN.md](MASTER_PLAN.md)

**Branch:** `extension/ancilla-hybrid-phase4`

**Current status:** Packaged for continuation; no scientific implementation has been undertaken on this branch by this packaging action.

| Gate | State | Evidence |
| --- | --- | --- |
| 0 — Baseline and audit | not_started | No run on branch yet |
| 1 — Operators/reference | not_started | — |
| 2 — Ancilla/Zeno | not_started | — |
| 3 — Gap/dynamical | not_started | — |
| 4 — Exact hybrid | not_started | — |
| 5 — Large-N/adaptive | not_started | — |
| 6 — Figures/CLI | not_started | — |
| 7 — Review | not_started | — |

## Baseline

Inherited from `extension/symmetry-scaling`. Baseline git SHA and test results **must be captured by Codex at first execution**; they have not been measured by the packaging task. Values quoted in earlier chats are not assumed to be reproducible.

## Last verified commit

Not yet recorded; record the current HEAD and the relevant tested commit after Gate 0.

## Verified tests

None run as part of the packaging action.

## Artifacts

Master specification and ledger only.

## Known limitations / open questions

- Implementations and plots referenced in earlier chat are not verified as committed repository assets.
- Exact hybrid reductions apply only to structured symmetry; ancilla projection is a protective operation, not generic exact propagation.
- Circuit/gate costs for a global symmetry oracle are not yet established.
- Large-N performance must be benchmarked under controlled memory and hardware conditions.

## Immediate next action

Run the standard repeatable Codex instruction. Audit branch/files/tests, reconcile this ledger, and complete Gate 0 before implementing Gate 1.

## Session log

- Packaging: created plan and ledger on isolated branch; no physics code altered.
