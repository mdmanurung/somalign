# Ideation session: somalign package improvements

**Date:** 2026-07-16
**Prompt:** "What else can we do to improve the somalign package?"
**Personas:** Statistical Physicist, Information Theorist, Topologist/Geometer,
Causal Inference Researcher, Representation Learning Specialist (2 ideas each)
**Grounding:** [CONTEXT.md](CONTEXT.md) + live read of `R/ot.R`, `R/fit.R`,
`R/anchored.R`, `R/utils.R`

10 ideas, grouped by effort. Every idea is grounded in specific somalign
functions/data structures. Full write-ups in the per-persona files.

## Low effort (quick wins)

| # | Persona | Idea | Hooks into | One-liner |
|---|---------|------|-----------|-----------|
| 2 | Stat. Physicist | Simulated-annealing Sinkhorn (epsilon cooling solver) | `.somalign_solve_internal_log`, `_solve_ot`, two_pass | Geometric epsilon cooling schedule as `solver="annealing"`; escapes rugged label-guided landscapes cold-start Sinkhorn gets trapped in |
| 3 | Info. Theorist | Mutual information as alignment-sharpness diagnostic + epsilon selector | `fit$transport_plan`, `diagnostics$nodes`, sensitivity_grid | `I(query node; ref node)` from the plan; per-node conditional entropy flags ambiguous nodes; info-rate-vs-cost elbow picks epsilon (`somalign_select_epsilon()`) |
| 8 | Causal Inference | Anchor exclusion-restriction test (Sargan–Hansen analog) | `fit$anchors`, batch subspace V, D | Permutation test that anchor displacements carry no signal orthogonal to V; validates the IV/negative-control assumption, flags biology×batch confounding |

## Medium effort

| # | Persona | Idea | Hooks into | One-liner |
|---|---------|------|-----------|-----------|
| 1 | Stat. Physicist | Epsilon phase-transition diagnostic + critical-epsilon estimator | `.somalign_solve_internal_log` (`lse_g`), `diagnostics$solver` | Plan row-entropy as order parameter; susceptibility peak locates critical epsilon; dual free energy (`log Z`) as new scalar diagnostic |
| 4 | Info. Theorist | Rate-distortion `outside_reference` via per-node surprisal | `somalign_results`, new `reference$node_var` | Calibrated chi-squared surprisal replaces ad-hoc quantile threshold; adds `outside_reference_pvalue` + `_top_marker` (pinpoints CD11c-style artifacts) |
| 5 | Topologist | Laplacian-regularized (smooth, curl-free) node-shift field | `shift_transform` hook in `.somalign_finish_fit`, `query$som$grid$pts` | Graph-Laplacian smoothing of the correction field over SOM lattice; `laplacian_lambda=0` default keeps backward compat |
| 6 | Topologist | Persistent-homology audit of node-shift topology | `somalign_diagnostics`, codebook point clouds | Persistence diagrams before/after correction detect population merging/erasure; `$topology` slot + `topology_warning` |
| 7 | Causal Inference | Sensitivity to unmeasured batch confounding (Rosenbaum-style) | `somalign_fit_anchored` subspace V, D, node_shifts | Bootstrap D → distribution over V → per-node correction CIs + "tipping angle" that reverses the correction (`somalign_subspace_sensitivity()`) |
| 9 | Repr. Learning | Learned diagonal Mahalanobis metric for batch-aware OT cost | `.somalign_pairwise_distance`, anchor D | Per-marker weights from D make batch directions cheap / biology expensive to transport; column whitening, no Sinkhorn changes |

## High effort (research-grade)

| # | Persona | Idea | Hooks into | One-liner |
|---|---------|------|-----------|-----------|
| 10 | Repr. Learning | Anchor-free batch subspace via contrastive multi-batch invariance | new `somalign_fit_multibatch()`, `.somalign_subspace_svd`, `shift_fn` | Estimate V from between-batch variance of pass-1 shifts across K batches → subspace correction *without anchors* |
| 6b | Topologist | (idea 6 at full H0+H1 with optional `TDA` dep) | — | Full persistent homology incl. loops; heavier, adds a suggested dependency |

## Cross-cutting themes

- **Principled epsilon/rho selection** shows up three times (ideas 1, 2, 3) from
  independent lenses — physics (critical temperature), info theory (rate-distortion
  elbow), and annealing. Strong signal this is the highest-value gap.
- **Better `outside_reference` / anomaly attribution** (idea 4) directly targets the
  real-world CD11c artifact from the BMV alignment.
- **Correction-field quality** (ideas 5, 6, 7) — smoothness, topology preservation,
  and uncertainty — all harden the barycentric correction the F2 fix just made a
  true Brenier map.
- **Anchor validity + anchor-free** (ideas 8, 10) address the biggest dependency of
  the signal-preserving subspace mode.
- **Metric learning** (idea 9) is the deepest lever: it changes what "close" means
  in every downstream step.
