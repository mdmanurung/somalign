# Plan 04 — Rate-Distortion `outside_reference` via Per-Node Surprisal

**Date:** 2026-07-16  
**Feature branch:** `feat/surprisal-outside-reference` (proposed)  
**Effort:** Medium (~3 days implementation + 1 day tests)  
**Depends on:** No new R dependencies; `stats::pchisq` is already imported.

---

## 1. Summary

Replace the ad-hoc per-node Euclidean-distance quantile threshold with a
calibrated chi-squared surprisal score. For each query cell, sum the squared
z-scores across all markers under a diagonal-Gaussian model of its assigned
reference node:

```
s_i = sum_f  (x_{i,f} - mu_{k,f})^2 / sigma^2_{k,f}
```

where `k` is the assigned reference node, `mu_{k,f}` = `reference$codebook[k,f]`
(node centroid in reference-scaled space), and `sigma^2_{k,f}` =
`reference$node_var[k,f]` (within-node per-marker variance, stored at reference
construction time). Under the null (query cell drawn from the same diagonal
Gaussian as reference node `k`), `s_i ~ chi-squared(df = n_features)`. The
p-value `pchisq(s_i, df = n_features, lower.tail = FALSE)` is calibrated: at
significance level `alpha`, the expected false-positive rate is exactly `alpha`
(subject to the Gaussian approximation).

This adds three new columns to `somalign_results()` output:
`outside_reference_surprisal`, `outside_reference_pvalue`, and
`outside_reference_top_marker`. The existing `outside_reference_distance` and
`final_status` columns are **unchanged**. No existing argument defaults are
broken.

---

## 2. Public API

### 2a. Reference constructors — new `compute_node_var` argument

Add `compute_node_var = TRUE` to all three public reference constructors.
Default `TRUE` so new references automatically store `node_var`; setting
`FALSE` reproduces the pre-feature behaviour (useful when memory is tight or
reference data are unavailable, e.g. `somalign_reference_from_nodes`).

```r
# reference.R, line 75
somalign_reference <- function(som_ref,
                               data,
                               labels = NULL,
                               features = NULL,
                               center = NULL,
                               scale = NULL,
                               codebook_space = NULL,
                               quantile_probs = c(0.5, 0.9, 0.95, 0.99),
                               compute_node_var = TRUE)   # <-- NEW

# reference.R, line 21
somalign_train_reference <- function(data,
                                     labels = NULL,
                                     features = NULL,
                                     grid = NULL,
                                     rlen = 100,
                                     alpha = c(0.05, 0.01),
                                     compute_node_var = TRUE,  # <-- NEW
                                     ...)

# reference.R, line 327
somalign_reference_from_som <- function(som,
                                        center,
                                        scale,
                                        codebook_space,
                                        labels = c("codebook", "none"),
                                        quantile_probs = c(0.5, 0.9, 0.95, 0.99),
                                        distance_chunk_size = 1e6L,
                                        compute_node_var = TRUE)  # <-- NEW
```

`somalign_reference_from_nodes` does NOT receive `compute_node_var`. It
accepts an optional `node_var` matrix directly (see Section 2b), enabling
serialise-then-reload workflows.

### 2b. `somalign_reference_from_nodes` — optional `node_var` argument

```r
# reference.R, line 162
somalign_reference_from_nodes <- function(codebook,
                                          features,
                                          center,
                                          scale,
                                          node_masses = NULL,
                                          label_prob = NULL,
                                          distance_quantiles = NULL,
                                          global_distance_quantiles = NULL,
                                          node_var = NULL)  # <-- NEW
```

When `node_var = NULL` (default), the field is stored as `NULL` and the
surprisal path emits a one-time `message()` and returns `NA` columns (see
Section 7).

### 2c. `somalign_results()` — new columns (additive)

`somalign_results()` receives no new arguments. Three columns are always
appended to the returned data frame:

| Column | Type | Description |
|---|---|---|
| `outside_reference_surprisal` | `numeric` | Chi-squared statistic `s_i`; `NA` when `node_var` absent |
| `outside_reference_pvalue` | `numeric` | `pchisq(s_i, df, lower.tail=FALSE)`; `NA` when absent |
| `outside_reference_top_marker` | `character` | Name of marker with largest `z^2`; `NA` when absent |

An optional `outside_pvalue_threshold` argument (default `NULL`) adds a fourth
boolean column `outside_reference_pvalue_flag` when provided:

```r
somalign_results <- function(fit, data = NULL,
                             outside_pvalue_threshold = NULL)
```

This keeps the 0-argument call path identical to today.

---

## 3. Internal Helpers

### 3a. `.somalign_node_var()` — compute per-node variance at reference time

**File:** `R/utils.R` (alongside `.somalign_distance_quantiles`, line 414)

```r
# Compute per-node per-marker within-node variance in reference-scaled space.
# scaled_data : N x p matrix (reference cells, already scaled)
# units       : integer vector length N (node assignment for each cell)
# n_nodes     : integer — number of codebook rows
# var_floor   : scalar floor applied via pmax() (default 1e-8)
# Returns     : n_nodes x p matrix; rows for empty nodes get global variance.
.somalign_node_var <- function(scaled_data, units, n_nodes,
                               var_floor = 1e-8) {
  p <- ncol(scaled_data)
  feat_names <- colnames(scaled_data)

  # Global fallback: variance of each feature across all cells
  global_var <- pmax(apply(scaled_data, 2L, stats::var), var_floor)

  split_idx <- split(seq_len(nrow(scaled_data)),
                     factor(units, levels = seq_len(n_nodes)))

  node_rows <- lapply(split_idx, function(idx) {
    if (length(idx) < 2L) {
      global_var                           # fallback for 0- or 1-cell nodes
    } else {
      pmax(apply(scaled_data[idx, , drop = FALSE], 2L, stats::var),
           var_floor)
    }
  })

  out <- matrix(unlist(node_rows, use.names = FALSE),
                nrow = n_nodes, ncol = p, byrow = TRUE)
  colnames(out) <- feat_names
  out
}
```

### 3b. `.somalign_node_surprisal()` — per-cell surprisal at results time

**File:** `R/utils.R` (or a new `R/surprisal.R`; keep in `utils.R` to avoid
a new file for <30 lines)

```r
# Compute chi-squared surprisal for each query cell.
# scaled_data : N x p matrix (query cells in reference-scaled space)
# units       : integer vector length N (reference node assignment per cell)
# reference   : somalign_reference object
# Returns     : list(surprisal, pvalue, top_marker)
#               All three are length-N vectors; NA when node_var unavailable.
.somalign_node_surprisal <- function(scaled_data, units, reference) {
  node_var <- reference$node_var
  if (is.null(node_var)) {
    message(
      ".somalign_node_surprisal: `reference$node_var` is absent ",
      "(reference built before somalign >= 0.99.2 or with ",
      "`compute_node_var = FALSE`). Surprisal columns will be NA."
    )
    n <- nrow(scaled_data)
    return(list(
      surprisal  = rep(NA_real_, n),
      pvalue     = rep(NA_real_, n),
      top_marker = rep(NA_character_, n)
    ))
  }

  codebook <- reference$codebook
  p  <- ncol(scaled_data)
  df <- p

  # Per-cell per-marker squared z-score: (x - mu_k)^2 / sigma^2_k
  # All vectorised via matrix indexing — same approach as .somalign_som_cell_distances
  resid  <- scaled_data - codebook[units, , drop = FALSE]
  var_k  <- node_var[units, , drop = FALSE]         # N x p
  z2     <- resid^2 / var_k                         # N x p

  surprisal  <- rowSums(z2)                          # N — chi-squared statistic
  pvalue     <- stats::pchisq(surprisal, df = df, lower.tail = FALSE)
  top_marker <- colnames(scaled_data)[max.col(z2, ties.method = "first")]

  list(surprisal = surprisal, pvalue = pvalue, top_marker = top_marker)
}
```

---

## 4. Data-Structure Changes

### 4a. New field: `reference$node_var`

A numeric matrix of dimension `n_nodes x n_features` stored inside the
`somalign_reference` list. Values are within-node per-marker variances computed
in reference-scaled space (i.e., after `center`/`scale` are applied to the raw
data). Values are floored at `1e-8` via `pmax()`.

**When `node_var` is `NULL`:** old reference objects and objects constructed
via `somalign_reference_from_nodes(node_var = NULL)` have no field or `NULL`.
All surprisal code paths check `is.null(reference$node_var)` and degrade
gracefully (see Section 7).

### 4b. Where `node_var` is computed

**Path 1 — `somalign_reference()` (reference.R line 100)**

After `projected <- .somalign_nearest_code(scaled_data, codebook)` (line 100)
and before the `structure(list(...))` call (line 106):

```r
# reference.R — insert after line 104
node_var <- if (isTRUE(compute_node_var)) {
  .somalign_node_var(scaled_data, projected$unit, n_nodes)
} else {
  NULL
}
```

Add `node_var = node_var` to the `structure(list(...))` on line 106.

**Path 2 — `somalign_train_reference()` (reference.R line 36)**

`somalign_train_reference` delegates to `somalign_reference()` at line 36.
Pass `compute_node_var = compute_node_var` through to that call:

```r
somalign_reference(
  som_ref        = som_ref,
  data           = data,
  labels         = labels,
  features       = colnames(data),
  center         = scaling$center,
  scale          = scaling$scale,
  codebook_space = "reference_scaled",
  compute_node_var = compute_node_var   # <-- pass through
)
```

**Path 3 — `somalign_reference_from_som()` (reference.R line 405)**

Inside `somalign_reference_from_som`, after the chunked distance computation
at line 405, the same `X` (reference-scaled training data) and `unit` vectors
are available. Compute `node_var` before delegating to
`somalign_reference_from_nodes`:

```r
# reference.R — insert after line 408 (quantiles computation)
node_var <- if (isTRUE(compute_node_var)) {
  .somalign_node_var(X, unit, n_nodes)
} else {
  NULL
}
```

Pass `node_var = node_var` to the `somalign_reference_from_nodes(...)` call
at line 411. This requires the matching `node_var` argument to be accepted
there (Section 2b).

**Path 4 — `somalign_reference_from_nodes()`**

Validate and store the supplied `node_var`:

```r
# reference.R — after existing distance_quantiles validation
node_var <- .somalign_prepare_node_var(node_var, n_nodes, features)
```

Where:

```r
# utils.R
.somalign_prepare_node_var <- function(node_var, n_nodes, features) {
  if (is.null(node_var)) return(NULL)
  node_var <- as.matrix(node_var)
  storage.mode(node_var) <- "double"
  if (nrow(node_var) != n_nodes)
    stop("`node_var` must have one row per reference node.", call. = FALSE)
  if (!is.null(features) && !is.null(colnames(node_var))) {
    node_var <- node_var[, features, drop = FALSE]
  }
  if (any(!is.finite(node_var)) || any(node_var < 0))
    stop("`node_var` must contain non-negative finite values.", call. = FALSE)
  node_var
}
```

Add `node_var = node_var` to the returned `structure(list(...))` in all four
construction paths.

---

## 5. Algorithm

Full vectorised R — no loops over cells.

```r
# Minimal self-contained illustration (not the actual helper, which handles NULL)
.somalign_node_surprisal_core <- function(scaled_data, units, codebook, node_var) {
  p   <- ncol(scaled_data)
  # Residual of each cell from its assigned node centroid
  resid <- scaled_data - codebook[units, , drop = FALSE]   # N x p
  # Per-cell per-marker squared z-score
  var_k <- node_var[units, , drop = FALSE]                 # N x p
  z2    <- resid^2 / var_k                                 # N x p
  # Chi-squared statistic (sum of squared z-scores, df = p)
  surprisal  <- rowSums(z2)
  pvalue     <- pchisq(surprisal, df = p, lower.tail = FALSE)
  # Top-contributing marker per cell
  top_marker <- colnames(scaled_data)[max.col(z2, ties.method = "first")]
  list(surprisal = surprisal, pvalue = pvalue, top_marker = top_marker)
}
```

Memory: `resid`, `var_k`, `z2` are each `N x p` doubles. For the BMV dataset
(N = 39.8M, p ~ 30), that is ~9 GB each — matching the existing peak in
`.somalign_project_samples`. If memory is a concern, the same
`chunk_size`-loop pattern used in `.somalign_som_cell_distances` (utils.R
line 590) can be applied here: process `chunk_size` rows at a time and write
only the three length-N output vectors.

A chunked variant `.somalign_node_surprisal_chunked()` should be added
alongside `.somalign_node_surprisal()`, mirroring the
`.somalign_nearest_code_chunked` / `.somalign_nearest_code` pair:

```r
.somalign_node_surprisal_chunked <- function(scaled_data, units, reference,
                                             chunk_size = 10000L) {
  n     <- nrow(scaled_data)
  node_var <- reference$node_var
  if (is.null(node_var))
    return(.somalign_node_surprisal(scaled_data, units, reference))
  codebook <- reference$codebook
  p        <- ncol(scaled_data)
  surprisal  <- numeric(n)
  pvalue     <- numeric(n)
  top_marker <- character(n)
  for (s in seq(1L, n, by = chunk_size)) {
    idx <- s:min(s + chunk_size - 1L, n)
    res <- .somalign_node_surprisal_core(
      scaled_data[idx, , drop = FALSE], units[idx], codebook, node_var
    )
    surprisal[idx]  <- res$surprisal
    pvalue[idx]     <- res$pvalue
    top_marker[idx] <- res$top_marker
  }
  list(surprisal = surprisal, pvalue = pvalue, top_marker = top_marker)
}
```

---

## 6. Integration Points

### 6a. Where new columns are added in `somalign_results()`

`somalign_results()` (results.R line 33) currently assembles `out` from
`fit$projection$direct`. The scaled data and reference-unit assignments needed
for surprisal are already present at fit time:

- `fit$query$scaled_data` — N x p query cells in reference-scaled space
  (see fit.R line 256, same matrix passed to `.somalign_project_samples`)
- `direct$unit` — reference node assignment per cell (result of
  `.somalign_nearest_code`)

Add the surprisal computation at the top of `somalign_results()`, immediately
after the existing `direct` / `corrected` assignment block (after line 40):

```r
# results.R — insert after line 40
surpr <- .somalign_node_surprisal_chunked(
  fit$query$scaled_data,
  direct$unit,
  fit$reference,
  chunk_size = getOption("somalign.chunk_size", 10000L)
)
```

Then extend the `data.frame(...)` call (line 52–77) with:

```r
outside_reference_surprisal = surpr$surprisal,
outside_reference_pvalue    = surpr$pvalue,
outside_reference_top_marker = surpr$top_marker
```

These three columns follow `outside_reference_distance` (line 58) to keep OOD
columns together.

When `outside_pvalue_threshold` is non-NULL, append after building `out`:

```r
if (!is.null(outside_pvalue_threshold)) {
  .somalign_check_prob_scalar(outside_pvalue_threshold,
                              "outside_pvalue_threshold")
  out$outside_reference_pvalue_flag <-
    !is.na(out$outside_reference_pvalue) &
    out$outside_reference_pvalue < outside_pvalue_threshold
}
```

### 6b. Coexistence with `outside_reference_distance`

Both columns are always present in the output data frame. They answer
complementary questions:

- `outside_reference_distance` (logical): distance-quantile heuristic, unit
  is Euclidean in reference-scaled space; the current default flag.
- `outside_reference_surprisal` + `outside_reference_pvalue`: calibrated
  chi-squared score that weights per-marker deviations by within-node
  variance; interpretable p-value.

Neither column is deprecated. Documentation notes that users may prefer the
surprisal flag when per-marker anomalies (e.g. CD11c batch artifacts) are
suspected, and the distance flag when a simple radius-based threshold is
sufficient.

---

## 7. Edge Cases

### 7a. Nodes with zero or near-zero within-node variance

Applied inside `.somalign_node_var()`:

```r
pmax(apply(scaled_data[idx, , drop = FALSE], 2L, stats::var), var_floor)
```

`var_floor = 1e-8` mirrors the `pmax(..., 0)` patterns already used in
`.somalign_sinkhorn_kernel` and `.somalign_pairwise_distance`. The consequence
of flooring is that a marker whose reference node has near-zero variance will
produce large `z^2` for any query cell that expresses it at all — which is the
correct and desired behaviour (a binary marker uniformly 0 in reference should
flag query cells that are nonzero). The idea-document's alternative of
excluding such markers (using effective df) is a future enhancement; document
it as a known limitation in the Rd.

### 7b. Nodes with 0 or 1 assigned cell (empty/singleton nodes)

Handled in `.somalign_node_var()`:

```r
if (length(idx) < 2L) global_var
```

For 0 cells, `stats::var()` returns `NA`; for 1 cell, `NaN`. Both cases fall
back to the global (across all reference cells) per-marker variance, which is
the most conservative available estimate. The `var_floor` is then applied on
top.

### 7c. Old reference objects without `node_var`

```r
if (is.null(reference$node_var)) {
  message(
    ".somalign_node_surprisal: `reference$node_var` is absent ...",
    "Surprisal columns will be NA."
  )
  return(list(surprisal = NA_real_[seq_len(n)], ...))
}
```

`is.null(reference$node_var)` evaluates `TRUE` for both references serialised
before this feature (no `node_var` key in the list at all) and references
built with `compute_node_var = FALSE`. Existing user code that saves and
reloads a reference object with `saveRDS`/`readRDS` will continue to work;
the new result columns will simply be `NA`.

### 7d. Feature mismatch (query has different marker set)

`fit$query$scaled_data` is already aligned to `reference$features` by the
time `somalign_results()` is called (enforced in `somalign_fit`). No
additional guard needed here.

### 7e. Very large N (BMV-scale: 39.8M cells)

Use `.somalign_node_surprisal_chunked()` with `chunk_size` matching the
existing projection chunk size. Peak additional memory over the current path:
`3 * N * p * 8` bytes for `resid`, `var_k`, `z2` within a chunk — identical
to the existing `.somalign_project_samples` footprint.

---

## 8. Tests

Add a new file: `tests/testthat/test-surprisal.R`.

### Test 1 — cell at node centroid has surprisal ≈ 0 and high p-value

```r
test_that("cell at node centroid scores near-zero surprisal", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_false(is.null(ref$node_var))

  # Place a cell exactly at node 1's centroid
  centroid <- matrix(ref$codebook[1, ], nrow = 1,
                     dimnames = list(NULL, c("F1", "F2")))
  surpr <- somalign:::.somalign_node_surprisal_core(
    centroid, 1L, ref$codebook, ref$node_var
  )
  expect_equal(surpr$surprisal, 0, tolerance = 1e-10)
  expect_gt(surpr$pvalue, 0.99)
})
```

### Test 2 — cell offset on one marker has that marker as top_marker

```r
test_that("single-marker offset identifies top_marker correctly", {
  set.seed(7)
  mat <- matrix(rnorm(300), nrow = 150, ncol = 3,
                dimnames = list(NULL, c("CD3", "CD11c", "CD19")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)

  centroid <- ref$codebook[1, , drop = FALSE]
  # Offset only CD11c by a large amount
  shifted  <- centroid
  shifted[, "CD11c"] <- centroid[, "CD11c"] + 10

  surpr <- somalign:::.somalign_node_surprisal_core(
    shifted, 1L, ref$codebook, ref$node_var
  )
  expect_equal(surpr$top_marker, "CD11c")
  expect_lt(surpr$pvalue, 0.01)
})
```

### Test 3 — old reference (no node_var) returns NA columns without error

```r
test_that("old reference without node_var gives NA surprisal columns", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5,
           compute_node_var = FALSE)
  expect_null(ref$node_var)

  qry <- somalign_query(mat, ref,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  expect_message(
    res <- somalign_results(fit),
    "node_var.*absent"
  )
  expect_true(all(is.na(res$outside_reference_surprisal)))
  expect_true(all(is.na(res$outside_reference_pvalue)))
  expect_true(all(is.na(res$outside_reference_top_marker)))
})
```

### Test 4 — `compute_node_var = TRUE` (default) stores a conforming matrix

```r
test_that("node_var has correct dimensions and is positive", {
  set.seed(3)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  nv <- ref$node_var
  expect_equal(nrow(nv), nrow(ref$codebook))
  expect_equal(ncol(nv), length(ref$features))
  expect_equal(colnames(nv), ref$features)
  expect_true(all(nv > 0))
})
```

### Test 5 — `somalign_results()` new columns present and numeric

```r
test_that("somalign_results includes surprisal columns", {
  set.seed(5)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  res <- somalign_results(fit)
  expect_true("outside_reference_surprisal" %in% names(res))
  expect_true("outside_reference_pvalue"    %in% names(res))
  expect_true("outside_reference_top_marker" %in% names(res))
  expect_true(all(is.numeric(res$outside_reference_surprisal)))
  expect_true(all(res$outside_reference_pvalue >= 0 &
                  res$outside_reference_pvalue <= 1))
})
```

### Test 6 — `outside_pvalue_threshold` argument adds boolean flag column

```r
test_that("outside_pvalue_threshold adds pvalue_flag column", {
  set.seed(9)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref,
           grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  res <- somalign_results(fit, outside_pvalue_threshold = 0.05)
  expect_true("outside_reference_pvalue_flag" %in% names(res))
  expect_type(res$outside_reference_pvalue_flag, "logical")
  expect_equal(
    res$outside_reference_pvalue_flag,
    res$outside_reference_pvalue < 0.05
  )
})
```

---

## 9. Docs and NAMESPACE

### Roxygen

- `somalign_reference()`: add `@param compute_node_var Logical; if `TRUE`
  (default) per-node per-marker variances are computed from the reference
  cells assigned to each node and stored as `reference$node_var`. Set
  `FALSE` to skip (reduces memory by one `n_nodes × p` matrix; disables
  surprisal-based outside-reference detection in `somalign_results()`)`.
- `somalign_train_reference()`: same `@param` entry, forwarded.
- `somalign_reference_from_som()`: same.
- `somalign_reference_from_nodes()`: add `@param node_var Optional
  `n_nodes × p` matrix of per-node per-marker variances (reference-scaled
  space). Supply when deserializing a previously computed reference.
  `NULL` disables surprisal columns in `somalign_results()`.`
- `somalign_results()`: extend `@return` to document the three new columns
  and add `@param outside_pvalue_threshold`.

### New plot function (optional, medium-effort; can ship in a follow-on PR)

```r
#' @export
somalign_plot_surprisal <- function(fit, results,
                                    alpha = 0.01,
                                    ncells_sample = 50000L)
```

Add to `NAMESPACE` via `@export` roxygen tag. The body: histogram of
`results$outside_reference_surprisal` with `dchisq(x, df = length(fit$reference$features))`
density overlay; vertical line at `qchisq(1 - alpha, df = ...)`. Colored
by `results$final_status`. BiocCheck: body will be under 50 lines.

### NAMESPACE

`devtools::document()` / `roxygen2::roxygenise()` regenerates `NAMESPACE`
automatically. No manual edits needed. The two new internal helpers
(`.somalign_node_var`, `.somalign_node_surprisal`, etc.) are unexported and
do not appear in `NAMESPACE`.

---

## 10. BiocCheck (<=50-line exported bodies)

The only exported function bodies that change are:

- `somalign_results()` (results.R): currently ~55 lines including the `out`
  data.frame construction. Adding the surprisal call and three columns will
  push the body to ~65 lines. Extract the `data.frame(...)` construction
  into a private helper `.somalign_results_df(...)` so the exported body
  stays under 50 lines — consistent with the existing pattern of delegating
  to private helpers.
- `somalign_reference()`, `somalign_train_reference()`,
  `somalign_reference_from_som()`, `somalign_reference_from_nodes()`:
  each gains only 2–3 lines; all remain well under 50 lines.
- `somalign_plot_surprisal()` (if added): will be <= 40 lines.

---

## 11. Effort, Risks, and Dependencies

### Effort

| Task | Est. |
|---|---|
| `.somalign_node_var()` + `.somalign_prepare_node_var()` in utils.R | 0.5 d |
| Plumb `compute_node_var` through three reference constructors | 0.5 d |
| `.somalign_node_surprisal()` + `_chunked()` + `_core()` | 0.5 d |
| Wire into `somalign_results()` | 0.5 d |
| Tests (6 cases above) | 1.0 d |
| Roxygen + optional plot | 0.5 d |
| **Total** | **~3.5 d** |

### Risks

1. **Compatibility with `somalign_reference_from_som`** — This path stores a
   copy of the SOM data as `X` (reference.R line 379, extracted via
   `.somalign_extract_som_data`). The `X` matrix is already in memory for
   the distance computation and is available for `node_var` without
   additional I/O. However, after `ref$som_ref$data <- NULL` (line 427), the
   original data is discarded. The `node_var` computation must happen
   *before* line 427; the plan above places it correctly (after line 408,
   before the `somalign_reference_from_nodes` call at line 411).

2. **Memory at BMV scale** — Computing `node_var` requires iterating over all
   reference cells grouped by node. For 39.8M cells × 30 markers, the
   `split()` approach allocates one sub-matrix per node. The chunked approach
   used in `.somalign_som_cell_distances` does *not* apply directly to
   variance computation (variance needs the full node's cells). An alternative
   is a two-pass online variance (Welford's method, O(N×p) passes) but plain
   R `apply/var` over split indices is typically fast enough. If memory is
   tight, add a `node_var_chunk_size` argument analogous to
   `distance_chunk_size`.

3. **Gaussian approximation** — SOM nodes are not Gaussian, and marker
   expression is lognormal. The chi-squared p-value is anti-conservative for
   heavy-tailed distributions. Document this in Rd and mention the Q-Q plot
   diagnostic. The surprisal score is still useful as a *relative ranking*
   regardless of distributional correctness.

4. **`somalign_reference_from_nodes` backward compatibility** — This path
   accepts a pre-built `node_var` matrix. Users who serialise a reference
   with `saveRDS` and pass `node_var` explicitly will need to match the
   expected `n_nodes × p` dimension. The `_prepare_node_var` validator
   gives clear error messages.

5. **Test coverage for `_from_som` path** — `test-reference-from-som.R` does
   not currently test `node_var`. Add one test there (not in the new file)
   that confirms `node_var` is present and conforming after
   `somalign_reference_from_som(...)`.

### Dependencies

- No new R package dependencies. `stats::pchisq` is in base `stats`, already
  imported.
- `compute_node_var = FALSE` provides a clean opt-out for users who do not
  want the extra memory overhead (relevant for large production references).
