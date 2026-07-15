# Analysis Manifest

Registry of analyses and ideation sessions for the somalign package.

## Ideation sessions

- **2026-07-16-methods-improvement** — `/mycelium:ideas` session on improving the
  somalign package. 5 methods-focused personas (Statistical Physicist, Information
  Theorist, Topologist/Geometer, Causal Inference Researcher, Representation Learning
  Specialist), 10 ideas. Index: `analysis/ideas/2026-07-16-methods-improvement/00_index.md`.
  Recurring signal: principled epsilon/rho selection (3 independent ideas).

## Implementation plans

- **2026-07-16-ideas-1-9** — Code-level implementation plans for ideas #1–#9 from the
  ideation session. Master: `analysis/plans/2026-07-16-ideas-1-9/00_IMPLEMENTATION_PLAN.md`
  (shared infra S1–S3, phased build order, file-touch matrix, commit sequence). Nine
  standalone `plan_0N_*.md` files. Key finding: #1 and #3 unify into one
  `somalign_epsilon_sweep()`; three ideas share a cheap OT-only sweep primitive and
  three share anchor-displacement storage.
