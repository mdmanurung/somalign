# Plan 09 — Learned Diagonal Mahalanobis Metric for Batch-Aware OT Cost

**Status:** Draft  
**Date:** 2026-07-16  
**Author:** mdmanurung  
**Idea source:** `analysis/ideas/2026-07-16-methods-improvement/05_representation-learning-specialist.md`, Idea 9

---

## 1. Summary

Replace the fixed squared-Euclidean OT cost in `.somalign_align_transport`
(fit.R line 210) with a **diagonal Mahalanobis cost** whose per-marker weights
are estimated from the anchor displacement matrix D. Markers that vary most
across the batch (high `var(D[,f])`) are batch-driven and get **low weight**
(cheap to transport); markers that are stable across the batch get **high
weight** (expensive to transport, so biology is preserved). The implementation
is a **pre-whitening** of both codebooks before passing them to the existing
`.somalign_pairwise_distance` — one helper function to estimate weights, one
helper to apply them, and a new argument on two exported functions. Sinkhorn,
barycentric correction, and all downstream code are unchanged. Projection and
threshold distances (`somalign_results`, `somalign_query`) call
`.somalign_nearest_code` / `.somalign_nearest_code_chunked` on the original
(unweighted) codebooks and are completely unaffected.

---

## 2. Public API

### 2a. New arguments

Add to **`somalign_fit`** (R/fit.R, line 76) and **`somalign_fit_anchored`**
(R/anchored.R, line 153):

```r
feature_weights = NULL
```

`feature_weights` is either:

- `NULL` (default) — pure squared-Euclidean cost, identical to current
  behaviour. Fully backward compatible.
- A named numeric vector of length `p` (one value per feature, names must
  match `colnames(query$codebook)`) — explicit weights supplied by the user.
  Useful when the user has prior knowledge of which markers carry the batch
  signal (e.g. from published gain-correction factors).
- The string `"anchor"` — auto-estimate from the anchor displacement matrix
  D using `.somalign_anchor_feature_weights()`. Only valid in
  `somalign_fit_anchored`; errors in `somalign_fit`.

Default is `NULL` in both functions, so no existing call breaks.

### 2b. Validation

Add a new validator in `R/utils.R`:

```r
.somalign_check_feature_weights <- function(fw, features) {
  if (is.null(fw) || identical(fw, "anchor")) return(invisible(fw))
  if (!is.numeric(fw) || !all(is.finite(fw)) || any(fw < 0))
    stop("`feature_weights` must be NULL, \"anchor\", or a non-negative numeric vector.",
         call. = FALSE)
  if (length(fw) != length(features))
    stop("`feature_weights` must have one entry per feature (", length(features), " expected).",
         call. = FALSE)
  if (!is.null(names(fw))) {
    missing <- setdiff(features, names(fw))
    if (length(missing) > 0)
      stop("`feature_weights` is missing names: ", paste(missing, collapse = ", "),
           call. = FALSE)
    fw <- fw[features]   # reorder to canonical column order
  } else {
    names(fw) <- features
  }
  invisible(fw)
}
```

---

## 3. Internal Helpers

### 3a. `.somalign_anchor_feature_weights(D, floor = 1e-2)`

**File:** `R/utils.R`  
**Signature:**

```r
.somalign_anchor_feature_weights <- function(D, floor = 1e-2) {
  # D: n_anchors x p matrix of anchor displacements
  #    = anchor_old_scaled - anchor_new_scaled (already in reference-scaled space)
  # floor: ridge term; prevents infinite weight for zero-variance markers.
  #        In reference-scaled space codebook entries are O(1), so 1e-2 is
  #        a sensible default (caps max weight at 100 x the min weight).
  v <- apply(D, 2, stats::var)         # length p; NA when n_anchors == 1
  v[is.na(v) | v < 0] <- 0            # safety: var() can return NA for n=1
  w <- 1 / (v + floor)                # high-variance markers -> low weight
  w / mean(w)                          # normalise: mean weight == 1
                                       # ensures cost_scale is comparable to
                                       # the unweighted case
}
```

The normalisation `w / mean(w)` keeps the median of the weighted cost matrix
on the same scale as the unweighted case, so `epsilon` values calibrated
without weights remain approximately valid (the median-normalisation in
`.somalign_prepare_cost` absorbs any residual scaling shift anyway).

### 3b. `.somalign_weighted_codebook(codebook, weights)`

**File:** `R/utils.R`  
**Signature:**

```r
.somalign_weighted_codebook <- function(codebook, weights) {
  # weights: named numeric vector length p, w_f > 0.
  # Returns codebook with column f scaled by sqrt(w_f).
  # Squared Euclidean distance in the scaled space equals
  #   sum_f w_f (q_if - r_jf)^2
  # which is the diagonal Mahalanobis cost.
  sweep(codebook, 2, sqrt(weights), "*")
}
```

This is the only change to how codebooks are passed into
`.somalign_pairwise_distance`. The existing function is completely unchanged.

### 3c. Integration in `.somalign_align_transport`

**File:** `R/fit.R`, current line 209–210:

```r
# CURRENT (line 210):
cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)

# NEW (replace lines 209-210 with):
if (!is.null(feature_weights)) {
  qcb <- .somalign_weighted_codebook(query$codebook,     feature_weights)
  rcb <- .somalign_weighted_codebook(reference$codebook, feature_weights)
  cost <- .somalign_pairwise_distance(qcb, rcb)
} else {
  cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)
}
```

Add `feature_weights = NULL` to `.somalign_align_transport`'s signature (line
205). Thread it from `somalign_fit` and `somalign_fit_anchored` through
`.somalign_align_transport`.

---

## 4. Data-Structure Changes

Store the resolved weights in `fit$diagnostics$cost_metric` so users can
inspect them and reproduce the run:

```r
# In .somalign_build_diagnostics (fit.R, line 264),
# add to the returned list:
cost_metric = list(
  feature_weights = feature_weights  # NULL or named numeric vector
)
```

Thread `feature_weights` into `.somalign_build_diagnostics` via
`.somalign_finish_fit` (new argument `feature_weights = NULL`).

For the anchored path, also store in `fit$anchors$feature_weights` the
resolved weight vector (even when estimated, not user-supplied) so the user
can retrieve `fit$anchors$feature_weights` for interpretation — which markers
were batch-drivers.

---

## 5. Algorithm Detail

### 5a. Weight formula

Given D (n_anchors × p), anchor displacements in reference-scaled space:

```
var_f = var(D[, f])                 # sample variance per marker
w_f   = 1 / (var_f + delta)        # delta = floor = 1e-2 (default)
w     = w / mean(w)                 # normalise mean to 1
```

### 5b. Whitened cost

With weights `w` (length p):

```
q_tilde[i, f] = sqrt(w_f) * q[i, f]   # query codebook row i, feature f
r_tilde[j, f] = sqrt(w_f) * r[j, f]   # reference codebook row j

cost[i, j] = sum_f (q_tilde[i,f] - r_tilde[j,f])^2
           = sum_f w_f * (q[i,f] - r[j,f])^2
```

This is computed by `.somalign_pairwise_distance(qcb, rcb)` after whitening.
No change to `.somalign_pairwise_distance` itself (utils.R line 381).

### 5c. Interaction with anchor cost_bonus

The anchor `cost_bonus` matrix (built by `.somalign_anchor_cost_bonus`,
anchored.R line 290) is subtracted from the **normalised** cost in
`.somalign_prepare_cost` (fit.R line 195):

```r
cost_normalized <- pmax(cost_normalized - cost_bonus, 0)
```

The bonus is a count-based matrix `rho_anchor * A / n_anchors` that is
independent of the feature-space metric — it is purely a node-pair
correspondence count. The Mahalanobis weighting enters before normalisation
(at line 210), and the bonus is subtracted after normalisation (line 195).
The two modifications are therefore **independent and composable**: the
feature weights reshape the geometry of the cost surface, while the cost bonus
biases routing toward anchor-supported pairs. No ordering conflict.

When `correction = "both"`, the `"subspace"` projection is applied to node
shifts **after** transport — it is also independent of the cost metric.

### 5d. Exact call trace (anchored path with `feature_weights = "anchor"`)

1. `somalign_fit_anchored(... feature_weights = "anchor")` (anchored.R)
2. `.somalign_validate_anchors` → returns `anchors_scaled` with
   `anchor_old_scaled` and `anchor_new_scaled`
3. `D <- anchors_scaled$anchor_old_scaled - anchors_scaled$anchor_new_scaled`
4. `fw <- .somalign_anchor_feature_weights(D, floor = 1e-2)`
5. `somalign_fit_anchored` passes `fw` to `.somalign_anchored_dispatch`
6. `.somalign_anchored_dispatch` passes `fw` to `.somalign_align_transport`
7. Inside `.somalign_align_transport`:
   - whitens both codebooks with `fw`
   - calls `.somalign_pairwise_distance` on whitened codebooks → `cost`
   - calls `.somalign_prepare_cost(cost, diagonal_boost, cost_bonus, label_mask)`
8. Sinkhorn, barycentric correction, projections proceed unchanged

---

## 6. Integration Points

### 6a. Exact location in fit.R

| Line (approx) | Change |
|---|---|
| 205 | Add `feature_weights = NULL` to `.somalign_align_transport` signature |
| 209–210 | Conditional whitening before `.somalign_pairwise_distance` |
| 264 | Add `feature_weights` to `.somalign_build_diagnostics` signature; store in `cost_metric` |
| 120 | Add `feature_weights = NULL` to `.somalign_finish_fit` signature |
| 109–112 | Pass `feature_weights` through to `.somalign_align_transport` and `.somalign_build_diagnostics` |

### 6b. Exact location in anchored.R

| Line (approx) | Change |
|---|---|
| 153 | Add `feature_weights = NULL` to `somalign_fit_anchored` |
| 183–188 | Resolve `feature_weights = "anchor"` → call `.somalign_anchor_feature_weights(D)` |
| 191 | Add `feature_weights` to `.somalign_anchored_dispatch` signature |
| 226–229 | Pass `feature_weights` to `.somalign_align_transport` |
| 231–244 | Store `feature_weights` in `anchors` list |

### 6c. somalign_fit (fit.R line 76)

```r
somalign_fit <- function(..., feature_weights = NULL) {
  # validate
  .somalign_check_feature_weights(feature_weights, colnames(query$codebook))
  if (identical(feature_weights, "anchor"))
    stop("`feature_weights = \"anchor\"` requires anchor data; use somalign_fit_anchored().",
         call. = FALSE)
  # pass through unchanged except feature_weights added to .somalign_align_transport call
}
```

### 6d. Projection distances stay Euclidean

`.somalign_nearest_code` (utils.R line 348) and
`.somalign_nearest_code_chunked` (utils.R line 359) are called from:

- `.somalign_project_samples` (fit.R line 469) — for `direct` and `corrected`
  projections used in `somalign_results`
- `.somalign_anchor_cost_bonus` (anchored.R line 300–305) — for mapping anchor
  cells to nodes

None of these receive or use `feature_weights`. They always operate on
unweighted codebooks in reference-scaled space. The Mahalanobis weights touch
**only** the `cost` matrix inside `.somalign_align_transport`.

---

## 7. Edge Cases

### 7a. Zero-variance marker (weight cap)

If `var(D[, f]) == 0` for all anchors (every anchor measured the same value
for marker `f`), the weight is `1 / (0 + delta) = 1 / delta`. With the default
`delta = 1e-2` this is `100` (in normalised units). That is large but finite and
bounded. No special handling is needed beyond the `floor` argument. Users who
want a lower cap can pass a larger `floor`.

If the entire column `D[, f]` is constant (e.g. a single anchor, so
`stats::var` returns `NA`), the guard `v[is.na(v)] <- 0` sets `var_f = 0`,
then `w_f = 1/delta`. Document this in the `floor` argument.

### 7b. Single anchor (n_anchors == 1)

`stats::var` returns `NA` for a length-1 vector. The guard above sets all
`v = 0`, so all weights equal `1/delta` before normalisation → all weights
become 1 after normalisation → the weighted cost is identical to unweighted.
A warning should be emitted:

```r
if (nrow(D) == 1L)
  warning("Only 1 anchor pair: feature weights cannot be estimated from variance; ",
          "using equal weights (equivalent to Euclidean cost). ",
          "Consider supplying explicit `feature_weights`.", call. = FALSE)
```

### 7c. Missing anchors when `feature_weights = "anchor"` on `somalign_fit`

Already handled in §6c: a clear error is thrown before any computation.

### 7d. `feature_weights` length mismatch

Caught by `.somalign_check_feature_weights` (§2b) before any computation.

### 7e. All-zero weight vector

If the user supplies a vector of all zeros, whitening collapses all columns to
zero and `.somalign_pairwise_distance` returns an all-zero cost matrix, which
causes `cost_scale = 1` (fallback in `.somalign_prepare_cost` line 186–188)
and a uniform Sinkhorn kernel — degenerate but not a crash. Add a check:

```r
if (all(fw == 0))
  stop("`feature_weights` must not be all zeros.", call. = FALSE)
```

in `.somalign_check_feature_weights`.

### 7f. Interaction with median cost normalisation

`.somalign_prepare_cost` (fit.R line 185) divides the cost by its median
positive entry:

```r
cost_scale <- stats::median(cost[cost > 0])
```

Because weights are normalised to `mean(w) == 1`, the median of the weighted
cost is approximately the same as the unweighted median (up to the change in
the cost distribution shape). This means the `epsilon` calibration from
unweighted runs remains approximately valid — users do not need to retune
`epsilon` when switching on feature weights. Exact equivalence is not
guaranteed (the distribution shifts), but the scale is controlled.

### 7g. Interaction with `diagonal_boost`

`diagonal_boost` operates on `cost_normalized` (line 190–193 of
`.somalign_prepare_cost`), after the median normalisation. It is independent
of `feature_weights`, which enters before normalisation. No conflict.

---

## 8. Tests

Add in `tests/testthat/test-mahalanobis-metric.R`:

### Test 1 — default is unchanged

```r
test_that("feature_weights = NULL produces identical fit to current default", {
  set.seed(42)
  mat <- rbind(matrix(rnorm(30 * 3, -2), ncol = 3),
               matrix(rnorm(30 * 3,  2), ncol = 3))
  colnames(mat) <- c("F1", "F2", "F3")
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.3, ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit_default  <- somalign_fit(qry, ref)
  fit_null     <- somalign_fit(qry, ref, feature_weights = NULL)
  expect_equal(fit_default$cost, fit_null$cost)
  expect_equal(fit_default$transport_plan, fit_null$transport_plan)
  expect_equal(fit_default$node_shifts,    fit_null$node_shifts)
})
```

### Test 2 — down-weighting a noisy marker changes the OT plan

```r
test_that("down-weighting a noisy marker shifts the OT plan", {
  set.seed(1)
  p <- 3L
  mat <- matrix(rnorm(40 * p), ncol = p, dimnames = list(NULL, paste0("F", 1:p)))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.5, ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  # Equal weights — should equal default
  fw_equal <- c(F1 = 1, F2 = 1, F3 = 1)
  fit_equal  <- somalign_fit(qry, ref, feature_weights = fw_equal)
  fit_default <- somalign_fit(qry, ref)
  expect_equal(fit_equal$cost, fit_default$cost)
  # Zero-out F1 weight — F1 should not contribute to cost
  fw_no_f1 <- c(F1 = 0, F2 = 1, F3 = 1)  # valid: not ALL zero
  # F1 whitened to 0 -> cost = only F2, F3 distances
  fit_no_f1 <- somalign_fit(qry, ref, feature_weights = fw_no_f1)
  # Manually compute expected cost using only F2, F3
  expected_cost <- .somalign_pairwise_distance(
    sweep(qry$codebook,     2, sqrt(fw_no_f1), "*"),
    sweep(ref$codebook,     2, sqrt(fw_no_f1), "*")
  )
  expect_equal(fit_no_f1$cost, expected_cost)
  # Transport plans differ between weighted and equal
  expect_false(isTRUE(all.equal(fit_no_f1$transport_plan, fit_default$transport_plan)))
})
```

### Test 3 — anchor-estimated weights: high-variance marker gets low weight

```r
test_that("anchor_feature_weights down-weights high-variance displacement marker", {
  set.seed(99)
  p <- 4L
  nm <- paste0("F", 1:p)
  mat <- matrix(rnorm(60 * p), ncol = p, dimnames = list(NULL, nm))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  # Batch effect: F1 has a large noisy component, F2-F4 are stable
  shift <- rbind(c(5, 0, 0, 0))  # only F1 shifts
  anc_idx <- 1:15
  anchor_old <- mat[anc_idx, , drop = FALSE]
  anchor_new <- mat[anc_idx, , drop = FALSE] + matrix(rnorm(15 * p) * c(2, 0.01, 0.01, 0.01),
                                                       ncol = p, byrow = FALSE)
  colnames(anchor_old) <- colnames(anchor_new) <- nm
  qry <- somalign_query(mat + shift[rep(1, nrow(mat)), ], ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit_anchored <- somalign_fit_anchored(qry, ref,
                                         anchor_old = anchor_old,
                                         anchor_new = anchor_new,
                                         feature_weights = "anchor")
  fw <- fit_anchored$anchors$feature_weights
  expect_true(fw["F1"] < fw["F2"])  # noisy F1 < stable F2
  expect_true(fw["F1"] < fw["F3"])
})
```

### Test 4 — projection distances unchanged vs default

```r
test_that("projection distances are identical regardless of feature_weights", {
  set.seed(7)
  mat <- matrix(rnorm(30 * 3), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.2, ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fw <- c(A = 0.1, B = 5, C = 1)
  fit_w <- somalign_fit(qry, ref, feature_weights = fw)
  fit_d <- somalign_fit(qry, ref)
  res_w <- somalign_results(fit_w)
  res_d <- somalign_results(fit_d)
  # Direct projection distances use .somalign_nearest_code, unaffected by weights
  expect_equal(res_w$outside_reference_distance, res_d$outside_reference_distance)
  expect_equal(res_w$old_som_unit, res_d$old_som_unit)
})
```

### Test 5 — explicit `feature_weights` vector is stored and retrieved

```r
test_that("explicit feature_weights are stored in diagnostics$cost_metric", {
  set.seed(3)
  mat <- matrix(rnorm(20 * 2), ncol = 2, dimnames = list(NULL, c("X", "Y")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fw <- c(X = 2, Y = 0.5)
  fit <- somalign_fit(qry, ref, feature_weights = fw)
  expect_equal(fit$diagnostics$cost_metric$feature_weights, fw[c("X", "Y")])
})
```

---

## 9. Documentation and NAMESPACE

### 9a. Roxygen param tags

Add to `somalign_fit` and `somalign_fit_anchored`:

```r
#' @param feature_weights Either `NULL` (default, squared-Euclidean cost),
#'   a named non-negative numeric vector with one entry per feature (explicit
#'   diagonal Mahalanobis weights), or the string `"anchor"` (only in
#'   `somalign_fit_anchored`: auto-estimates weights from the anchor
#'   displacement matrix D so batch-driven markers are down-weighted and
#'   biology-carrying markers are up-weighted). Weights are squared internally
#'   (`sqrt(w_f)` per-column scaling of both codebooks) before the squared
#'   Euclidean distance is computed, yielding cost
#'   \eqn{\sum_f w_f (q_{if} - r_{jf})^2}. Set to `NULL` to preserve the
#'   current behaviour. When estimated from anchors, the resolved weight vector
#'   is stored in `fit$anchors$feature_weights` and
#'   `fit$diagnostics$cost_metric$feature_weights`.
```

### 9b. NAMESPACE

No new exported functions. `.somalign_anchor_feature_weights` and
`.somalign_weighted_codebook` are internal helpers (dot-prefixed). No
`NAMESPACE` entry required.

### 9c. NEWS.md entry

```
# somalign 0.99.2 (development)

## New features

* `somalign_fit()` and `somalign_fit_anchored()` gain a `feature_weights`
  argument for a diagonal Mahalanobis OT cost. In the anchored path,
  `feature_weights = "anchor"` auto-estimates per-marker weights from the
  anchor displacement matrix so batch-driven markers are made cheap to
  transport and biology-carrying markers are penalised. Projection and
  threshold distances are unaffected.
```

---

## 10. BiocCheck Compliance (<=50-line exported bodies)

`somalign_fit` currently has ~33 lines of body (lines 89–118 inclusive). Adding
`feature_weights` validation (~3 lines) and its pass-through (~1 line) keeps the
exported function body well under 50 lines.

`somalign_fit_anchored` currently has ~37 lines (lines 170–189). Same
assessment — comfortably under 50.

All heavy logic is in internal helpers (`.somalign_anchor_feature_weights`,
`.somalign_weighted_codebook`, `.somalign_check_feature_weights`), which are
not exported and not subject to the 50-line check.

---

## 11. Effort, Risks, and Dependencies

### 11a. Effort

**Medium-small.** Estimated 2–3 hours of net coding:

- `R/utils.R`: ~25 new lines (`_anchor_feature_weights`, `_weighted_codebook`,
  `_check_feature_weights`)
- `R/fit.R`: ~15 modified/added lines (signature + conditional whitening +
  diagnostics threading)
- `R/anchored.R`: ~20 modified/added lines (signature + resolve + threading +
  anchors list)
- `tests/testthat/test-mahalanobis-metric.R`: ~90 lines (5 tests above)
- Roxygen + NEWS: ~15 lines

Total: ~165 lines changed/added.

### 11b. Risks

**Risk 1 — `epsilon` calibration shift.** The median-normalisation in
`.somalign_prepare_cost` absorbs most of the scale change introduced by
weights, but the shape of the cost distribution changes. Users who have
carefully tuned `epsilon` on the unweighted problem may need to re-tune after
enabling feature weights. Mitigated by: (a) the mean-normalisation of weights
to `mean(w) == 1`, which anchors the scale; (b) documenting the interaction
clearly; (c) `diagnostics$solver$cost_scale` already exposes the normalisation
factor for inspection.

**Risk 2 — Small anchor n.** With fewer anchors than features (p > n_anchors,
common in CyTOF: 30–50 markers, 10–50 anchor cells), per-marker variance
estimates are noisy. The `floor` parameter provides regularisation. For very
small anchor sets (n_anchors < 5), the Test 3 assertion may be unreliable.
Document a recommended minimum of n_anchors ≈ 20. For a future version,
cross-validated shrinkage (James-Stein toward the mean weight) could be added.

**Risk 3 — Interaction with F2 fix.** The F2 fix (switched cost from plain
Euclidean to squared Euclidean) restored the Brenier OT-map property of the
barycentric correction. The Mahalanobis extension `sum_f w_f (q-r)_f^2` is
still a squared (Bregman-like) cost and preserves the Brenier property as long
as weights are strictly positive and fixed. The `floor` term guarantees strict
positivity. No regression on F2.

**Risk 4 — Interaction with `cost_bonus`.** As described in §5c, the two
modifications are independent and commutative (weights reshape the geometry,
bonus modifies routing after normalisation). No conflict found in the call
trace. The risk is that a user sets both `feature_weights = "anchor"` and a
large `rho_anchor`, which concentrates the plan on anchor-supported routes AND
reshapes the geometry — the combined effect may be hard to interpret. Document
that these are complementary levers; the anchor feature weights refine the
metric globally while the cost bonus gives explicit routing preference.

**Risk 5 — `zero` weight entries.** The check `all(fw == 0)` guards the
degenerate all-zero case. A single zero weight (`F1 = 0`) collapses one
dimension from the cost but is mathematically valid and may be intentional
(e.g. a dead channel). Allow it with no warning; document that zeroing a
marker removes it from the OT cost entirely (that marker is free to shift
without cost, which may cause over-correction on that marker).

### 11c. Dependencies

- No new R package dependencies.
- Uses `stats::var` (base), `sweep` (base), `apply` (base).
- Builds on `.somalign_pairwise_distance` (utils.R line 381), which is stable
  post-F2.
- Orthogonal to the `subspace` correction path and `cost_bonus` path.
- No change to `somalign_results`, `somalign_diagnostics`, `somalign_query`,
  or `somalign_reference`.
- Compatible with `somalign_fit_two_pass`: if desired, `feature_weights` could
  be threaded to two-pass in a follow-up PR (not in scope here).

### 11d. Deepest lever

This is the deepest single-lever improvement available without architectural
change: it directly modifies what "distance" means in the OT objective,
concentrating transport mass on biologically coherent directions. Combined with
the F2 fix (which restored the proper squared-cost OT-map property) and the
anchor `cost_bonus` (which biases routing), the three mechanisms address
geometry (Mahalanobis), marginal balance (unbalanced KL), and routing
(cost_bonus) independently — a principled decomposition of the correction
problem.
