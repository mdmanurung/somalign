# Plan 06 — Persistent-Homology Audit of Node-Shift Topology

**Feature:** Idea #6 from `03_topologist-geometer.md`
**Status:** Design / not started
**Date:** 2026-07-16

---

## 1. Summary

SOM-level batch correction can silently merge or erase biological populations
(topological distortion). No current diagnostic catches this. This plan adds a
lightweight topology audit that computes the H0 persistence diagram of three
codebook point clouds — query, corrected-query, and reference — and reports how
many robustly-separated clusters survive at a biologically motivated distance
threshold. The result lives in a new `$topology` slot inside the existing
`diagnostics` list returned by `somalign_diagnostics(fit)`. A
`topology_warning` flag is raised when the corrected codebook has fewer
components than the query, indicating population merging. An optional richer
H0+H1 path is activated when the `TDA` package is available.

The implementation is pure-R for H0, adds no hard dependencies, is additive
(existing `diagnostics` structure is unchanged unless `topology_audit = TRUE`
is passed), and keeps exported function bodies under 50 lines (BiocCheck
requirement).

---

## 2. Public API

### 2a. `somalign_topology_audit()` — new exported function

```r
#' Compute persistent-homology topology audit for a somalign fit
#'
#' @param fit A `somalign_fit` object.
#' @param threshold Numeric scalar in reference-scaled marker space.
#'   H0 components with persistence (death - birth) greater than this value
#'   are counted as "robustly separated populations". `NULL` (default) derives
#'   the threshold automatically from `reference$distance_quantiles` as
#'   described in Details.
#' @param use_tda Logical. When `TRUE`, use `TDA::ripsDiag()` (H0 + H1) if
#'   the `TDA` package is available. Falls back silently to base-R H0 when
#'   `TDA` is absent. Default `FALSE`.
#' @param nodes Subset of nodes to include. One of `"correction_allowed"`
#'   (default, recommended) or `"all"`.
#' @return A named list of class `somalign_topology`; see Details.
#' @export
somalign_topology_audit <- function(fit,
                                    threshold = NULL,
                                    use_tda   = FALSE,
                                    nodes     = c("correction_allowed", "all")) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  nodes <- match.arg(nodes)
  .somalign_check_flag(use_tda, "use_tda")
  if (!is.null(threshold))
    .somalign_check_pos_scalar(threshold, "threshold")
  .somalign_topology_audit_impl(fit, threshold, use_tda, nodes)
}
```

File: `R/diagnostics.R` (append after line 262).

### 2b. Opt-in addition to `somalign_diagnostics()`

`somalign_diagnostics(fit)` (line 17–22 of `diagnostics.R`) currently just
returns `fit$diagnostics`. That return value is assembled in
`.somalign_build_diagnostics()` (`fit.R`, line 264). The `$topology` slot is
**not** added automatically during `somalign_fit()` (too expensive to compute
by default, and it is a post-hoc diagnostic). Instead:

```r
somalign_diagnostics <- function(fit, topology = FALSE, ...) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  diag <- fit$diagnostics
  if (isTRUE(topology))
    diag$topology <- somalign_topology_audit(fit, ...)
  diag
}
```

This is backward-compatible: default `topology = FALSE` leaves existing output
identical. Callers who want the slot pass `topology = TRUE`.

### 2c. Return structure of `somalign_topology_audit()`

```r
list(
  # Threshold used for "robustly separated" counting
  threshold         = <numeric scalar>,
  threshold_source  = "auto" | "user",   # how threshold was chosen

  # Per-cloud component counts at the threshold
  n_components_query      = <integer>,
  n_components_corrected  = <integer>,
  n_components_reference  = <integer>,

  # Change relative to query
  topology_delta          = <integer>,   # corrected - query; negative = merging

  # Raw H0 persistence diagrams (data frames: birth, death, persistence)
  diagram_query      = <data.frame>,
  diagram_corrected  = <data.frame>,
  diagram_reference  = <data.frame>,

  # Optional: bottleneck distance between corrected and reference H0 diagrams
  # NULL when TDA is absent
  bottleneck_h0            = <numeric | NULL>,

  # Optional: full TDA ripsDiag output (H0 + H1), NULL when TDA absent
  tda_query      = <list | NULL>,
  tda_corrected  = <list | NULL>,
  tda_reference  = <list | NULL>,

  # Warning flag
  topology_warning  = <logical>   # TRUE when |topology_delta| > 0
)
```

---

## 3. Internal Helpers

All helpers live in **`R/utils.R`** (append after line 399 — after
`.somalign_normalize_masses`).

### 3a. Union-find (base R)

```r
# Internal: initialize union-find parent vector
.uf_make <- function(n) seq_len(n)

# Internal: find root with path compression (modifies parent in-place via <<-)
.uf_find <- function(parent, i) {
  while (parent[i] != i) {
    parent[i] <<- parent[parent[i]]   # path halving
    i <- parent[i]
  }
  i
}

# Internal: union two components; return updated parent and whether merge happened
.uf_union <- function(parent, i, j) {
  ri <- .uf_find(parent, i)
  rj <- .uf_find(parent, j)
  if (ri == rj) return(list(parent = parent, merged = FALSE))
  parent[rj] <- ri
  list(parent = parent, merged = TRUE)
}
```

Note: `.uf_find` uses `<<-` for path compression; call it only inside a
closure or wrapper that owns `parent` in its local environment — see
`.somalign_h0_persistence` below.

### 3b. H0 persistence via single-linkage

Signature: `.somalign_h0_persistence(D)` where `D` is an M×M symmetric
distance matrix (Euclidean, NOT squared — take `sqrt` of
`.somalign_pairwise_distance` output). Returns a data frame of
`(birth=0, death, persistence)` for all M-1 merges, sorted by persistence
descending.

```r
.somalign_h0_persistence <- function(D) {
  M <- nrow(D)
  if (M <= 1L)
    return(data.frame(birth = numeric(0), death = numeric(0),
                      persistence = numeric(0)))
  # Upper-triangle edges sorted by distance
  idx <- which(upper.tri(D), arr.ind = TRUE)
  dists <- D[idx]
  ord   <- order(dists)
  idx   <- idx[ord, , drop = FALSE]
  dists <- dists[ord]

  parent <- seq_len(M)          # union-find state (local env for <<-)
  deaths <- numeric(M - 1L)     # we get M-1 merges
  k <- 0L
  for (e in seq_len(nrow(idx))) {
    ri <- .uf_find(parent, idx[e, 1L])
    rj <- .uf_find(parent, idx[e, 2L])
    if (ri != rj) {
      parent[rj] <- ri
      k <- k + 1L
      deaths[k] <- dists[e]
      if (k == M - 1L) break
    }
  }
  data.frame(birth = 0, death = deaths[seq_len(k)],
             persistence = deaths[seq_len(k)])
}
```

### 3c. Component count from diagram

```r
# Count H0 components with persistence > threshold.
# The "infinite" component (last to survive) is counted separately as 1;
# finite components are those with persistence > threshold.
.somalign_h0_n_components <- function(diagram, threshold, n_nodes) {
  if (n_nodes <= 0L) return(0L)
  if (nrow(diagram) == 0L) return(1L)  # single node or all identical
  sum(diagram$persistence > threshold) + 1L
}
```

### 3d. Threshold auto-selection from reference geometry

```r
.somalign_topo_threshold <- function(reference) {
  dq <- reference$distance_quantiles
  if (is.null(dq)) return(0.5)
  # 95th percentile per-node threshold is stored as a named vector or matrix
  q95 <- if (is.matrix(dq)) dq["95%", ] else dq["95%"]
  q95 <- q95[is.finite(q95)]
  if (length(q95) == 0L) return(0.5)
  stats::median(q95)
}
```

### 3e. Optional TDA path

```r
.somalign_tda_diagram <- function(CB, max_dim = 1L) {
  if (!requireNamespace("TDA", quietly = TRUE)) return(NULL)
  TDA::ripsDiag(CB, maxdimension = max_dim, maxscale = Inf,
                dist = "euclidean", library = "GUDHI")
}

.somalign_tda_bottleneck <- function(d1, d2) {
  if (is.null(d1) || is.null(d2)) return(NULL)
  if (!requireNamespace("TDA", quietly = TRUE)) return(NULL)
  # Extract H0 birth-death pairs
  h0_extract <- function(d) {
    diag <- d$diagram
    diag[diag[, "dimension"] == 0, c("Birth", "Death"), drop = FALSE]
  }
  TDA::bottleneck(h0_extract(d1), h0_extract(d2))
}
```

### 3f. Main implementation dispatcher

Signature in `R/diagnostics.R` (or a new `R/topology.R` — see §3g):

```r
.somalign_topology_audit_impl <- function(fit, threshold, use_tda, nodes) {
  allowed <- attr(fit$node_shifts, "correction_allowed")
  sel <- if (nodes == "correction_allowed") which(allowed) else
           seq_len(nrow(fit$query$codebook))

  CB_query     <- fit$query$codebook[sel, , drop = FALSE]
  CB_corrected <- CB_query + fit$node_shifts[sel, , drop = FALSE]
  CB_ref       <- fit$reference$codebook

  thresh <- if (is.null(threshold))
    .somalign_topo_threshold(fit$reference) else threshold
  thresh_source <- if (is.null(threshold)) "auto" else "user"

  # Euclidean distance matrices (NOT squared — PH operates on metric spaces)
  D_query  <- sqrt(.somalign_pairwise_distance(CB_query, CB_query))
  D_corr   <- sqrt(.somalign_pairwise_distance(CB_corrected, CB_corrected))
  D_ref    <- sqrt(.somalign_pairwise_distance(CB_ref, CB_ref))

  pd_q <- .somalign_h0_persistence(D_query)
  pd_c <- .somalign_h0_persistence(D_corr)
  pd_r <- .somalign_h0_persistence(D_ref)

  nq <- .somalign_h0_n_components(pd_q, thresh, nrow(CB_query))
  nc <- .somalign_h0_n_components(pd_c, thresh, nrow(CB_corrected))
  nr <- .somalign_h0_n_components(pd_r, thresh, nrow(CB_ref))

  tda_q <- tda_c <- tda_r <- bn <- NULL
  if (isTRUE(use_tda) && requireNamespace("TDA", quietly = TRUE)) {
    tda_q <- .somalign_tda_diagram(CB_query)
    tda_c <- .somalign_tda_diagram(CB_corrected)
    tda_r <- .somalign_tda_diagram(CB_ref)
    bn    <- .somalign_tda_bottleneck(tda_c, tda_r)
  }

  warn <- (nc - nq) != 0L
  if (warn) {
    direction <- if (nc < nq) "merged/erased" else "split"
    warning(sprintf(
      "topology_warning: corrected codebook has %d H0 component(s) vs %d in query ",
      nc, nq),
      sprintf("(delta = %+d; populations may have been %s). ",
              nc - nq, direction),
      "Inspect fit$diagnostics$topology for details.",
      call. = FALSE)
  }

  structure(
    list(threshold = thresh, threshold_source = thresh_source,
         n_components_query = nq, n_components_corrected = nc,
         n_components_reference = nr, topology_delta = nc - nq,
         diagram_query = pd_q, diagram_corrected = pd_c,
         diagram_reference = pd_r, bottleneck_h0 = bn,
         tda_query = tda_q, tda_corrected = tda_c, tda_reference = tda_r,
         topology_warning = warn),
    class = "somalign_topology")
}
```

### 3g. File placement

| Symbol | File |
|---|---|
| `.uf_make`, `.uf_find`, `.uf_union` | `R/utils.R` (after line 399) |
| `.somalign_h0_persistence` | `R/utils.R` (after union-find block) |
| `.somalign_h0_n_components` | `R/utils.R` |
| `.somalign_topo_threshold` | `R/utils.R` |
| `.somalign_tda_diagram`, `.somalign_tda_bottleneck` | `R/utils.R` |
| `.somalign_topology_audit_impl` | `R/diagnostics.R` (append) |
| `somalign_topology_audit` | `R/diagnostics.R` (append) |
| Updated `somalign_diagnostics` | `R/diagnostics.R` (replace lines 17–22) |

---

## 4. Data-Structure Changes

### 4a. New `$topology` slot (additive to existing `diagnostics`)

The existing list returned by `somalign_diagnostics(fit)` has four slots:
`$solver`, `$ot`, `$nodes`, `$projection`. A fifth slot `$topology` is added
only when `topology = TRUE` is passed:

```
fit$diagnostics$topology
├── threshold             <dbl>    # scale in reference-space Euclidean units
├── threshold_source      <chr>    # "auto" | "user"
├── n_components_query    <int>    # H0 count for query codebook
├── n_components_corrected<int>    # H0 count for corrected codebook
├── n_components_reference<int>    # H0 count for reference codebook
├── topology_delta        <int>    # corrected - query
├── diagram_query         <df>     # cols: birth, death, persistence (M-1 rows)
├── diagram_corrected     <df>     # same structure
├── diagram_reference     <df>     # same structure
├── bottleneck_h0         <dbl|NULL>
├── tda_query             <list|NULL>  # TDA::ripsDiag output
├── tda_corrected         <list|NULL>
├── tda_reference         <list|NULL>
└── topology_warning      <lgl>
```

### 4b. `topology_warning` semantics

`TRUE` when `topology_delta != 0L` (i.e., the number of robustly-separated
clusters changed after correction). The R warning message is emitted
automatically inside `.somalign_topology_audit_impl` so it appears in the
console even when the result is not inspected. `topology_delta < 0` means
populations were merged/erased; `topology_delta > 0` means populations were
spuriously split.

---

## 5. Algorithm

### 5a. H0 persistence via single-linkage (base R)

H0 persistent homology of a finite metric space is equivalent to the
single-linkage dendrogram. The algorithm is:

1. Compute the M×M symmetric Euclidean distance matrix `D` (sqrt of
   `.somalign_pairwise_distance` output — see §6 for the squared-distance
   interaction).
2. Enumerate all M(M-1)/2 edges; sort ascending by distance.
3. Run union-find (Kruskal's algorithm): process edges in order; when edge
   (i,j) merges two components, record `death = dist(i,j)`. This generates
   exactly M-1 merge events.
4. The H0 persistence diagram has M points: M-1 with `(birth=0, death=d_k)` and
   one essential class `(birth=0, death=Inf)`.
5. `persistence = death - birth = death` for all finite classes.
6. Count classes with `persistence > threshold`.

For M ≤ 1024 (typical SOM sizes: 64–400 nodes) the O(M² log M) sort dominates
and runs in milliseconds. For M=400, the edge list has 79 800 entries; sorting
and union-find is ~1 ms in base R.

**R snippet for the full H0 computation:**

```r
.somalign_h0_persistence <- function(D) {
  M <- nrow(D)
  if (M <= 1L) return(data.frame(birth=numeric(0), death=numeric(0),
                                 persistence=numeric(0)))
  idx   <- which(upper.tri(D), arr.ind = TRUE)
  dists <- D[idx]
  ord   <- order(dists)
  idx   <- idx[ord, , drop = FALSE]; dists <- dists[ord]
  parent <- seq_len(M)
  deaths <- numeric(M - 1L); k <- 0L
  for (e in seq_len(nrow(idx))) {
    ri <- .uf_find(parent, idx[e,1L])
    rj <- .uf_find(parent, idx[e,2L])
    if (ri != rj) { parent[rj] <- ri; k <- k+1L; deaths[k] <- dists[e]
      if (k == M-1L) break }
  }
  data.frame(birth=0, death=deaths[seq_len(k)],
             persistence=deaths[seq_len(k)])
}
```

Note: `.uf_find` uses `<<-` for path compression and must be called with
`parent` in the enclosing frame; in the loop above, `parent` is a local
variable in the same frame, so `<<-` reaches it correctly. An alternative
is to pass `parent` by reference through an environment object if the
`<<-` behavior is considered fragile in future refactors.

**Component count curve:** The number of H0 components as a function of
threshold `t` is `1 + sum(deaths > t)`. This is a step function that starts
at M (all nodes isolated) and decreases to 1 (fully connected). Plotting this
curve against `t` yields a barcode-like stability picture.

### 5b. Threshold selection from reference geometry

`reference$distance_quantiles` stores, for each reference node, the 95th
percentile of per-cell distances to that node. The median of these quantiles
gives a natural "within-population spread" scale. Populations separated by
more than this scale are considered robustly distinct. The formula is:

```r
threshold <- stats::median(reference$distance_quantiles["95%", ])
```

If `distance_quantiles` is a named vector rather than a matrix, use
`reference$distance_quantiles["95%"]`. Check the structure at runtime:
`is.matrix(dq)` to branch correctly (see `.somalign_topo_threshold` in §3d).

### 5c. Optional TDA path (H0 + H1)

When `use_tda = TRUE` and `TDA` is available:

```r
if (requireNamespace("TDA", quietly = TRUE)) {
  tda_q <- TDA::ripsDiag(CB_query, maxdimension = 1L, maxscale = Inf,
                          dist = "euclidean", library = "GUDHI")
  # bottleneck distance between corrected and reference H0 diagrams:
  bn <- TDA::bottleneck(
    h0_pairs(tda_c$diagram),
    h0_pairs(tda_r$diagram)
  )
}
```

where `h0_pairs` extracts the sub-matrix of rows with `dimension == 0`.

Graceful degradation: when `TDA` is absent, `tda_*` slots are `NULL`,
`bottleneck_h0` is `NULL`, and no error is raised. The base-R H0 result is
complete and self-contained. A `message()` (not `warning()`) is emitted once
when `use_tda = TRUE` but `TDA` is absent:

```r
if (isTRUE(use_tda) && !requireNamespace("TDA", quietly = TRUE))
  message("TDA package not available; falling back to base-R H0 only.")
```

---

## 6. Integration Points

### 6a. Squared-distance interaction (F2 fix)

`.somalign_pairwise_distance()` (`utils.R`, line 381) returns **squared**
Euclidean distances. Its docstring confirms: "Squared Euclidean distance. Used
only to build the OT cost matrix." Persistent homology requires a metric
(triangle inequality); squared distances violate it. Therefore:

```r
D_query <- sqrt(.somalign_pairwise_distance(CB_query, CB_query))
```

Always apply `sqrt()` before passing to `.somalign_h0_persistence()`. The
threshold is then in Euclidean (not squared) units, matching
`reference$distance_quantiles` which stores plain Euclidean distances (from
`.somalign_nearest_code`, line ~360 of `utils.R`). This is the single most
important integration note.

### 6b. Corrected codebook reconstruction

The corrected codebook is not stored directly in `fit`; reconstruct it as:

```r
CB_corrected <- fit$query$codebook + fit$node_shifts
```

For the `nodes = "correction_allowed"` case, subset both:

```r
sel <- which(attr(fit$node_shifts, "correction_allowed"))
CB_corrected <- fit$query$codebook[sel, ] + fit$node_shifts[sel, ]
CB_query     <- fit$query$codebook[sel, ]
```

This matches the recommendation in CONTEXT.md (`03_topologist-geometer.md`,
line 97): "PH should be computed only over the `correction_allowed = TRUE`
nodes, since the others are unchanged — zeroing corrected nodes' shifts avoids
spurious topology changes from the forced-zero entries."

### 6c. Reference-scaled space

Both `query$codebook` and `reference$codebook` are already in reference-scaled
space (the space used by the OT cost). No re-scaling is needed.

### 6d. Integration with `somalign_sensitivity_grid()`

Add `topology_delta` to `.somalign_grid_row_summary()` (diagnostics.R, line
123) to track topology changes across the ε/ρ grid. Requires passing
`topology = TRUE` to `somalign_diagnostics()` in that function; gate behind an
opt-in parameter `topology_audit = FALSE` on `somalign_sensitivity_grid()` to
keep the default fast:

```r
# in .somalign_grid_row_summary(), add:
topo_delta = if (isTRUE(topology_audit)) {
  ta <- somalign_topology_audit(fit)
  ta$topology_delta
} else NA_integer_
```

---

## 7. Edge Cases

| Scenario | Handling |
|---|---|
| `nrow(CB_query) == 1` | `.somalign_h0_persistence` returns empty data frame; `n_components = 1`; no warning. |
| All nodes identical (all distances = 0) | All M-1 merges happen at `death=0`; all components have `persistence=0`; every threshold > 0 gives `n_components = 1`. |
| `correction_allowed` all FALSE | `sel` is empty; early return with `n_components_* = 0` for query/corrected, compute reference normally; set `topology_warning = FALSE` with a message. |
| Very small SOM (M=4) | Algorithm works; `M-1 = 3` merges; no special-casing needed. |
| Threshold too large (> max pairwise dist) | All `persistence < threshold`; `n_components = 1`; not a bug — document in help. |
| Threshold too small (→ 0) | `n_components = M`; document as degenerate case. |
| TDA present but GUDHI backend missing | `TDA::ripsDiag` will error; wrap in `tryCatch` and fall back to base-R H0 with a warning. |
| High-dimensional marker space (p > 20) | H0 in R^p is valid but may be dominated by curse of dimensionality; document that users may want to run on PCA-reduced codebook; add a `pca_dims = NULL` parameter (future work) — not implemented in this plan. |

**Never auto-modify the fit.** This is a pure diagnostic. The function reads
`fit` and returns a new object; it does not write to `fit$diagnostics` in
place. If the warning fires, it is advisory only.

---

## 8. Tests

File: `tests/testthat/test-topology.R`

### Test 1 — two-cluster codebook reports 2 H0 components

```r
test_that(".somalign_h0_persistence detects 2 clusters", {
  # Two tight clusters far apart in 2D
  CB <- matrix(c(0,0, 0.1,0, 0,0.1,
                 10,10, 10.1,10, 10,10.1), ncol=2, byrow=TRUE)
  D  <- sqrt(somalign:::.somalign_pairwise_distance(CB, CB))
  pd <- somalign:::.somalign_h0_persistence(D)
  # Threshold between within-cluster spread (~0.14) and between-cluster (~14)
  n  <- somalign:::.somalign_h0_n_components(pd, threshold = 1, n_nodes = 6L)
  expect_equal(n, 2L)
})
```

### Test 2 — merging correction triggers topology_warning

```r
test_that("topology_warning fires when correction merges clusters", {
  set.seed(42)
  # Reference: two clusters
  ref_cb <- matrix(c(0,0, 0,0.1, 0.1,0,
                     5,5, 5,5.1, 5.1,5), ncol=2, byrow=TRUE)
  # Query: same two clusters (aligned)
  qry_cb <- ref_cb + 0.05
  # node_shifts that collapse both clusters to centroid (2.5, 2.5)
  shifts <- matrix(2.5, nrow=6, ncol=2) - qry_cb
  attr(shifts, "correction_allowed") <- rep(TRUE, 6L)

  # Build minimal fit-like structure for audit
  mock_fit <- structure(
    list(
      query     = list(codebook = qry_cb),
      reference = list(codebook = ref_cb,
                       distance_quantiles = c("95%" = 0.3)),
      node_shifts = shifts
    ),
    class = "somalign_fit"
  )
  # Expect a warning about topology
  expect_warning(
    ta <- somalign_topology_audit(mock_fit),
    regexp = "topology_warning"
  )
  expect_true(ta$topology_warning)
  expect_lt(ta$n_components_corrected, ta$n_components_query)
})
```

### Test 3 — TDA-absent path returns H0 only without error

```r
test_that("topology audit works without TDA package", {
  skip_if(requireNamespace("TDA", quietly = TRUE),
          "TDA is installed; test checks TDA-absent path only")
  set.seed(1)
  mat <- matrix(rnorm(20), nrow=10, ncol=2,
                dimnames=list(NULL, c("F1","F2")))
  ref <- somalign_train_reference(mat,
           grid=kohonen::somgrid(2,2,"hexagonal"), rlen=5)
  qry <- somalign_query(mat, ref,
           grid=kohonen::somgrid(2,2,"hexagonal"), rlen=5)
  fit <- somalign_fit(qry, ref)
  ta  <- expect_no_error(somalign_topology_audit(fit, use_tda=TRUE))
  expect_null(ta$bottleneck_h0)
  expect_null(ta$tda_query)
  expect_type(ta$n_components_corrected, "integer")
})
```

### Test 4 — single-node codebook does not error

```r
test_that(".somalign_h0_persistence handles M=1 without error", {
  CB <- matrix(1:2, nrow=1)
  D  <- sqrt(somalign:::.somalign_pairwise_distance(CB, CB))
  pd <- somalign:::.somalign_h0_persistence(D)
  expect_equal(nrow(pd), 0L)
  n  <- somalign:::.somalign_h0_n_components(pd, threshold=0.1, n_nodes=1L)
  expect_equal(n, 1L)
})
```

### Test 5 — `somalign_diagnostics(topology=TRUE)` is additive

```r
test_that("somalign_diagnostics topology=TRUE adds $topology without breaking defaults", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow=10, ncol=2,
                dimnames=list(NULL, c("F1","F2")))
  ref <- somalign_train_reference(mat,
           grid=kohonen::somgrid(2,2,"hexagonal"), rlen=5)
  qry <- somalign_query(mat, ref,
           grid=kohonen::somgrid(2,2,"hexagonal"), rlen=5)
  fit <- somalign_fit(qry, ref)
  d0  <- somalign_diagnostics(fit)
  d1  <- somalign_diagnostics(fit, topology=TRUE)
  expect_null(d0$topology)
  expect_s3_class(d1$topology, "somalign_topology")
  # All other slots unchanged
  expect_equal(d0$solver, d1$solver)
  expect_equal(d0$ot,     d1$ot)
})
```

---

## 9. Docs / NAMESPACE

### 9a. DESCRIPTION — add TDA to Suggests

```
Suggests:
    BiocParallel,
    BiocStyle,
    DiagrammeR,
    knitr,
    pkgdown,
    rmarkdown,
    TDA,
    testthat (>= 3.0.0),
    withr
```

`TDA` must NOT appear in `Imports`; it is a large Rcpp-based package with
system dependencies (CGAL). Using `requireNamespace("TDA", quietly = TRUE)`
at runtime satisfies Bioconductor policy for optional backends.

### 9b. NAMESPACE

`somalign_topology_audit` is exported via `@export` in its Roxygen block.
`print.somalign_topology` should also be exported (see below).

### 9c. print method

```r
#' @export
print.somalign_topology <- function(x, ...) {
  cat(sprintf(
    "somalign topology audit\n  threshold: %.4f (%s)\n",
    x$threshold, x$threshold_source))
  cat(sprintf(
    "  H0 components  query: %d  corrected: %d  reference: %d\n",
    x$n_components_query, x$n_components_corrected, x$n_components_reference))
  cat(sprintf("  topology_delta: %+d  warning: %s\n",
    x$topology_delta, x$topology_warning))
  invisible(x)
}
```

---

## 10. BiocCheck Compliance

BiocCheck requires exported function bodies ≤ 50 lines. `somalign_topology_audit`
has 12 lines. `somalign_diagnostics` (updated) has ~10 lines. All heavy logic
lives in internal helpers prefixed `.somalign_*` which are not exported and
not subject to the 50-line limit. The `print.somalign_topology` method is 6
lines. All exported bodies are well under 50 lines.

---

## 11. Effort, Risks, and Dependencies

### Effort

| Sub-task | Estimate |
|---|---|
| Union-find + `.somalign_h0_persistence` | 1–2 h |
| `.somalign_topo_threshold` | 0.5 h |
| `.somalign_topology_audit_impl` | 1–2 h |
| `somalign_topology_audit` exported function + Roxygen | 1 h |
| Updated `somalign_diagnostics` signature | 0.5 h |
| TDA optional path + graceful degradation | 1 h |
| Tests (5 test cases) | 2–3 h |
| Documentation + DESCRIPTION + NAMESPACE | 0.5 h |
| **Total** | **7–10 h** |

Effort is "medium" — the heaviest of the Idea #6 scope but tractable solo.
The implementation is straightforward once the squared-distance interaction
is handled correctly.

### Key Risks

1. **Squared-distance interaction (F2 fix).** `.somalign_pairwise_distance()`
   returns squared distances. Passing squared distances to
   `.somalign_h0_persistence()` without `sqrt()` would give wrong persistence
   scales and a wrong threshold. This is the single most likely implementation
   bug. Mitigation: the square-root is taken explicitly in
   `.somalign_topology_audit_impl` and documented with an inline comment.

2. **Threshold sensitivity.** The `topology_delta` is sensitive to the
   threshold parameter. Two codebook point clouds with slightly different
   geometry can flip between 2 and 3 components depending on threshold choice.
   Mitigation: auto-threshold is derived from `reference$distance_quantiles`
   (a biologically grounded scale already used elsewhere); allow user override;
   document the step-function nature of the count.

3. **`reference$distance_quantiles` structure.** The code in
   `.somalign_topo_threshold` branches on `is.matrix(dq)`. If the structure
   changes in future `somalign_reference` versions, the branch may fail.
   Mitigation: add a fallback default (`return(0.5)`) and a `tryCatch` wrapper.

4. **High-dimensional marker space.** In p > 20 dimensions, H0 components may
   not correspond to visually obvious clusters because nearest-neighbour
   distances concentrate. This is a known limitation of TDA in high dimensions,
   not a bug. Document clearly in `?somalign_topology_audit`; recommend
   projecting to PCA top-5 components for visualization (deferred feature).

5. **`<<-` in `.uf_find`.** The path-compression step uses `<<-` to modify
   `parent` in the enclosing scope. This is idiomatic R but unusual. If
   `.uf_find` is ever called outside a frame that owns `parent`, it will modify
   the wrong variable. Mitigation: keep all three union-find helpers tightly
   coupled; consider replacing with an environment-based approach in a future
   refactor.

6. **TDA system dependencies.** `TDA` requires CGAL/Boost at build time. It
   may not be installable on all HPC environments (e.g., the current
   `para-lipg-hpc` cluster). The optional path + `requireNamespace` check
   handles this correctly; the package must not fail R CMD check when `TDA` is
   absent.

### New hard dependencies

None. `stats` (already in `Imports`) is the only standard-library function
used (`stats::median`). `TDA` is `Suggests` only.
