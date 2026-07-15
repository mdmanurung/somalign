# Master implementation plan — somalign ideas #1–#9

**Date:** 2026-07-16
**Scope:** concrete, code-level implementation of the 9 methods-focused ideas from
`analysis/ideas/2026-07-16-methods-improvement/`. Each idea has a full standalone
plan (`plan_0N_*.md`, ~6,100 lines total); this file is the **integration layer**:
shared infrastructure, build order, file-touch matrix, cross-cutting risks, and an
atomic commit sequence.

**Guiding invariants (apply to every idea):**
- **Backward compatible.** Every new feature defaults to *exact* current behavior
  (new args default to `NULL`/`0`/`"euclidean"`; new diagnostic/result fields are
  additive; old reference/fit objects degrade gracefully with a message, not an error).
- **BiocCheck.** Exported function bodies ≤ 50 lines → delegate to internal `.somalign_*`.
- **F2 awareness.** `.somalign_pairwise_distance()` now returns **squared** Euclidean.
  Any idea needing true distances (#6) must `sqrt()`; the OT cost path (#9) stays squared.
  Projection/threshold distances use `.somalign_nearest_code` (separate, Euclidean) —
  never repurpose it (#4, #9).
- Tests are testthat edition-3; `devtools::document()` regenerates NAMESPACE/Rd.

---

## Idea overview

| # | Feature | New public API | Effort | Standalone plan |
|---|---------|----------------|--------|-----------------|
| 1 | Epsilon phase-transition diagnostic + critical-ε | `somalign_epsilon_sweep()` (+`print`/`plot`); `diagnostics$solver$log_Z` | Med | plan_01 |
| 2 | Annealing Sinkhorn (ε cooling solver) | `solver="annealing"` + `anneal_*` args | Low–Med | plan_02 |
| 3 | Mutual-information diagnostic + ε selector | `somalign_select_epsilon()`; `diagnostics$ot$mutual_information`, `$nodes$transport_entropy` | Low | plan_03 |
| 4 | Surprisal `outside_reference` (rate-distortion) | `reference$node_var`; `outside_reference_surprisal/_pvalue/_top_marker` cols | Med | plan_04 |
| 5 | Laplacian-smoothed node-shift field | `laplacian_lambda=0` arg | Med | plan_05 |
| 6 | Persistent-homology audit | `somalign_topology_audit()`; `diagnostics(topology=TRUE)$topology` | Med–High | plan_06 |
| 7 | Batch-subspace confounding sensitivity | `somalign_subspace_sensitivity()` | Med | plan_07 |
| 8 | Anchor exclusion-restriction test | `somalign_exclusion_test()` | Low | plan_08 |
| 9 | Learned Mahalanobis OT cost | `feature_weights=NULL` arg | Med | plan_09 |

---

## Phase 0 — Shared infrastructure (build FIRST; small, unblocks the rest)

Three primitives are reused by multiple ideas. Building them once, cleanly, avoids
three divergent copies. Each is a tiny, backward-compatible change.

### S1. Cheap OT-only sweep primitive — `.somalign_ot_sweep_one()` (`R/fit.R`)
Runs `.somalign_align_transport` (cost build + `.somalign_solve_ot`) for one
`(epsilon, rho_query, rho_ref)` **without** `.somalign_finish_fit`/per-cell
projection, returning the M×K plan plus cheap plan-level summaries (row-entropy,
expected cost, and — via S-hooks — MI and `log_Z`).
- Signature (canonical, from plan_03):
  `.somalign_ot_sweep_one(query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol, diagonal_boost = 0, label_mask = NULL)` → `list(plan, row_entropy, expected_cost, ...)`.
- **Consumers:** #1 (`somalign_epsilon_sweep`), #3 (`somalign_select_epsilon`,
  `sensitivity_grid` MI column). Both iterate this over an ε grid.
- **~30 lines.** Extract the pre-projection half of the current fit path; the sweep
  functions just `lapply` over the grid and rbind.

### S2. Dual warm-start in the log-domain solver (`R/ot.R`)
Add optional `f_init`/`g_init` to `.somalign_solve_internal_log()` and **return the
converged `f`/`g`** in its result list (currently discarded).
- **Consumers:** #2 (annealing hands potentials stage→stage); also lets S1 warm-start
  down an ε grid (10–50× fewer Sinkhorn iters on fine grids).
- **~5 lines**, no behavior change when `f_init=NULL` (current cold start).

### S3. Persist anchor displacements in the fit (`R/anchored.R`)
In `.somalign_anchored_dispatch()`, compute `D_scaled <- anchor_old_scaled -
anchor_new_scaled` once and store `fit$anchors$displacements <- D_scaled` (and thread
`variance_threshold` into `batch_subspace`).
- **Consumers:** #7 (bootstrap D), #8 (orthogonal-residual test), #9 (anchor-derived
  feature weights can read D directly).
- **~2 lines**, additive to `fit$anchors`.

> Phase 0 is ~40 lines total and each piece is independently committable.

---

## Consolidation recommendation (do before coding #1 and #3)

Ideas **#1 and #3 are the same sweep** over `.somalign_ot_sweep_one()` with different
readouts. **Build one `somalign_epsilon_sweep()`** (plan_01) whose returned `$table`
has *all* per-ε columns:

`epsilon | row_entropy | susceptibility | log_Z (free energy) | mutual_information | expected_cost | mean_conditional_entropy`

Then **`somalign_select_epsilon()` (plan_03) becomes a thin selector** on that table
(`method = c("critical","elbow","entropy_fraction")` → critical-ε = susceptibility
peak (#1), elbow = max-curvature of the MI/cost curve (#3)). This unifies the
physics and information-theory lenses into one object + one selector, instead of two
overlapping sweep implementations. **Net saving: one full sweep implementation and
one test file.** The standalone plans still apply for the individual diagnostics
(`log_Z`, per-node entropy) that also attach to a normal `somalign_fit`.

---

## Build order (phases; each phase is shippable)

**Phase 0 — infra:** S1, S2, S3 (above). ~40 lines, 3 commits.

**Phase 1 — quick wins (Low effort, high value):**
1. **#3 + #1 unified sweep** — `somalign_epsilon_sweep()` + `somalign_select_epsilon()`;
   attach `mutual_information`, `transport_entropy`, `log_Z` to normal-fit diagnostics.
   *(depends on S1, S2)* — this is the **highest-value gap** (3 personas converged on
   principled ε selection).
2. **#8 exclusion test** — `somalign_exclusion_test()`. *(depends on S3)* — smallest new
   feature; strong scientific-validity payoff for the anchor/subspace mode.
3. **#2 annealing solver** — `solver="annealing"`. *(depends on S2)* — recompute
   `tau_a/tau_b` per stage (ε-dependent!); helps the label-guided / small-ε regime.

**Phase 2 — Medium, independent:**
4. **#9 Mahalanobis cost** — `feature_weights` (explicit or `"anchor"` from S3's D).
   Guard: only the OT cost path; projection stays Euclidean; composes with cost_bonus.
5. **#4 surprisal outside_reference** — `reference$node_var` (compute before
   `som$data` is discarded in the `_from_som` path; watch 39.8M-cell memory) + new
   result columns; `_top_marker` pinpoints CD11c-style artifacts.
6. **#5 Laplacian shifts** — `laplacian_lambda`; reuse the `shift_transform` hook
   (`fit.R:126,143–145`); **compose** with the subspace transform (smooth-then-project),
   preserving the `correction_allowed` attribute.
7. **#7 subspace sensitivity** — `somalign_subspace_sensitivity()`. *(depends on S3)* —
   bootstrap D, principal angles, per-node CIs, tipping angle.

**Phase 3 — heaviest:**
8. **#6 persistent-homology audit** — `somalign_topology_audit()`; base-R H0 union-find
   (remember `sqrt()` — F2!), optional `TDA` in Suggests; `$topology` slot is opt-in.

---

## File-touch matrix

| File | S1 | S2 | S3 | #1 | #2 | #3 | #4 | #5 | #6 | #7 | #8 | #9 |
|------|----|----|----|----|----|----|----|----|----|----|----|----|
| `R/ot.R` | | ✎ | | ✎log_Z | ✎ | | | | | | | |
| `R/fit.R` | ✎new | | | ✎ | ✎plumb | ✎ | | ✎arg+hook | | | | ✎weight cost |
| `R/anchored.R` | | | ✎ | | ✎plumb | | | (compose) | | (read) | | ✎weights |
| `R/reference.R` | | | | | | | ✎node_var | | | | | |
| `R/results.R` | | | | | | | ✎cols | | | | | |
| `R/utils.R` | | | | | | ✎MI helpers | ✎surprisal | ✎laplacian | (sqrt dist) | ✎subspace boot | ✎ortho resid | ✎weights |
| `R/diagnostics.R` | | | | ✎sweep/print | | ✎MI fields | | | ✎topology | | | |
| `R/plot.R` | | | | ✎plot | | | | | | | | |
| new `R/*.R` | | | | | | | | | topology.R | sensitivity.R | exclusion.R | |
| `DESCRIPTION` | | | | | | | | | +Suggests TDA | | | |
| `tests/testthat/` | | | | + | + | + | + | + | + | + | + | + |

✎ = edit existing; "+" = new test file per idea.

---

## Cross-cutting risks & how the plans handle them

- **ε-dependence of `tau_a/tau_b`** (#2): the UOT proximal exponents change with ε;
  the annealing driver must recompute them each stage. Explicit in plan_02.
- **F2 squared distances** (#6, #9): #6 must `sqrt(.somalign_pairwise_distance(...))`
  for genuine PH; #9 deliberately keeps squared (whitened columns → weighted squared
  cost). Both plans flag it.
- **Memory at BMV scale** (#4): `node_var` for 39.8M cells must be computed from the
  SOM partition before `som$data` is dropped, reusing chunked projection. plan_04
  mitigates.
- **`shift_transform` composition** (#5 vs subspace #7/anchored): only one hook exists;
  Laplacian smoothing and subspace projection must **compose** (and both must restore
  the `correction_allowed` attribute). plan_05 specifies smooth-then-project.
- **Non-grid query SOMs** (#5): `query$som_query$grid$pts` may be absent for some
  `somalign_query_from_som` paths → `laplacian_lambda>0` errors cleanly / no-ops.
- **Subspace prerequisites** (#7/#8): require `correction="subspace"/"both"` fits;
  error clearly on plain `cost_bonus` fits.

---

## Suggested atomic commit sequence

Phase 0: `feat(ot): return + accept warm-start dual potentials (S2)` ·
`feat(fit): add cheap OT-only epsilon-sweep primitive (S1)` ·
`feat(anchored): persist scaled anchor displacements in fit$anchors (S3)`.

Phase 1: `feat(diagnostics): epsilon sweep with free-energy, MI, entropy + selector (#1,#3)` ·
`test: epsilon sweep + selector` ·
`feat(anchored): exclusion-restriction test for anchors (#8)` ·
`feat(ot): simulated-annealing epsilon-cooling solver (#2)`.

Phase 2: `feat(fit): learned Mahalanobis OT cost via feature_weights (#9)` ·
`feat(reference): per-node variance + surprisal outside_reference columns (#4)` ·
`feat(fit): Laplacian-smoothed node-shift field (#5)` ·
`feat(anchored): bootstrap subspace confounding sensitivity (#7)`.

Phase 3: `feat(diagnostics): persistent-homology topology audit (#6)` ·
`docs: add TDA to Suggests`.

Each `feat` commit pairs with its test; run `devtools::document()` +
`devtools::test()` per phase; `R CMD check`/BiocCheck at the end of each phase.

---

## Verification (per phase)
1. `devtools::document()` — clean NAMESPACE/Rd regen.
2. `devtools::test()` — full suite green incl. new tests; every default-off path proven
   byte-identical to pre-change fits (regression assertions in #2/#5/#9 tests).
3. `R CMD check` + BiocCheck — 0 errors; exported bodies ≤ 50 lines; no undeclared deps
   (TDA only in Suggests, guarded by `requireNamespace`).
4. Smoke: re-run a small anchored fit through the new diagnostics; confirm additive
   fields appear and old code paths are unchanged.

---

## Effort roll-up
- **Phase 0:** ~0.5 day (40 lines, high leverage).
- **Phase 1 (#1/#3, #8, #2):** ~3–4 days; unlocks principled ε selection + anchor validity.
- **Phase 2 (#9, #4, #5, #7):** ~5–7 days.
- **Phase 3 (#6):** ~1.5 days.
- **Total:** ~10–13 developer-days for all nine, if built in this order so the three
  shared primitives are written once.
