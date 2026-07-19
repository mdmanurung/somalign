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

## Validation analyses

- **2026-07-19-other-t-improvement** — Batch1-to-Batch2 labelled BMV pilot
  validation for seven OtherT/abundance improvement approaches. Folder:
  `analysis/other-t-improvement/`. Key finding: low epsilon and T/NK marker
  weights improve deployable OtherT recall/F1, but no tested deployable setting
  improves both OtherT precision and recall over the `epsilon = 0.1` baseline.
  The SOMalign example notebook provides separate true repeat metadata
  (`repeated_samples.csv`): hard nearest-node metacluster CLR weighted-r is
  0.9322 and soft kNN is 0.9566 across 20 pilot-vs-BMV repeat samples.

- **2026-07-19-groundup-performance-review** — Multi-team ground-up review of
  SOMalign performance. Folder: `analysis/groundup-performance-review/`.
  Four read-only specialist passes covered mathematical/statistical formulation,
  OtherT/T/NK biology, validation/publication readiness, and runtime/API
  engineering. Key consensus: the next performance gains should prioritize
  calibrated decision logic, strict local T/NK handling for OtherT, soft
  posterior/frequency unification, automated release gates, and projection-memory
  hardening rather than replacing the OT solver.
