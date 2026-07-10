# somalign patch implementation plan

## Status note

Several LOW-severity patches have already been applied to the working tree at the
time this plan was written: E10 (`@examples` blocks), E11 (DESCRIPTION author),
E12 (DESCRIPTION Suggests), E13 (`parallel` argument), and E14 (NEWS.md). Their
sections are included below for completeness with status **ALREADY APPLIED**. All
remaining patches must be applied in the order given.

---

## 1. Summary table

| Issue ID | Severity | File | Description | Review Status |
|----------|----------|------|-------------|---------------|
| A1a | HIGH | R/ot.R | Initialise `delta <- Inf` before Sinkhorn loop | Approved (rationale note only) |
| A1b | HIGH | R/ot.R | Warn on non-convergence after Sinkhorn loop | CONDITIONAL — two blockers (see section) |
| A2 | HIGH | R/ot.R | Underflow-fraction check + warning for kernel collapse | Approved |
| D5b | MEDIUM | tests/testthat/test-ot-warnings.R | Tests for A1b non-convergence warning | New file required |
| D5c | MEDIUM | tests/testthat/test-ot-warnings.R | Tests for A2 underflow warning | New file required |
| C7a | MEDIUM | R/results.R | Vectorize `.somalign_reference_top_labels` | Approved |
| C7d | MEDIUM | R/fit.R | Vectorize `.somalign_transfer_labels` | CONDITIONAL — second_label divergence (see section) |
| B7b | MEDIUM | R/utils.R | Vectorize `.somalign_distance_quantiles` | Approved |
| B7c | MEDIUM | R/utils.R | Vectorize `.somalign_label_probabilities` | CONDITIONAL — dead code in patch (see section) |
| B3 (utils) | MEDIUM | R/utils.R | Add `.somalign_nearest_code_chunked` | Approved |
| B3 (fit private) | MEDIUM | R/fit.R | `.somalign_project_samples` calls chunked function | Approved |
| B3 (fit public) | MEDIUM | R/fit.R | Add `chunk_size` parameter to `somalign_fit()` | Approved |
| B3 (threading) | MEDIUM | R/fit.R | Thread `chunk_size` through both projection call sites | Approved |
| D5a | MEDIUM | tests/testthat/test-chunked-projection.R | Tests for B3/B7b/B7c | CONDITIONAL — test assertion fix required (see section) |
| C4 | MEDIUM | R/fit.R | `message()` when `match_mass_ratio > 1` | Approved |
| C6 | MEDIUM | R/fit.R | `@details` roxygen section explaining marginal deviation | CONDITIONAL — old_snippet mismatch (see section) |
| D8 | MEDIUM | tests/testthat/test-training-integration.R | Test `parallel = TRUE` branch of sensitivity grid | New test required |
| A9 | LOW | R/ot.R | Suppress benign POT `reg_type='entropy'` warning | CONDITIONAL — scope too broad (see section) |
| E13 | LOW | R/diagnostics.R | `parallel` argument + `requireNamespace` guard | ALREADY APPLIED (guard missing — see section) |
| E10 | LOW | Multiple R files | `@examples \dontrun{}` blocks | ALREADY APPLIED |
| E11 | LOW | DESCRIPTION | Replace placeholder author | ALREADY APPLIED |
| E12 | LOW | DESCRIPTION | Add `bench`, `microbenchmark`, `parallel` to Suggests | ALREADY APPLIED |
| E14 | LOW | NEWS.md | Replace stub NEWS.md with full changelog | ALREADY APPLIED |

---

## 2. Implementation order

Apply patches in this sequence to minimise risk. Run `devtools::test()` after
each numbered batch before proceeding.

**Batch 1 — HIGH fixes (standalone, no dependencies)**
1. A1a — initialise `delta`
2. A1b — non-convergence warning (apply with blocker fixes incorporated)
3. A2 — kernel underflow warning

**Batch 2 — Tests for HIGH fixes**
4. D5b + D5c — `tests/testthat/test-ot-warnings.R` (new file)

**Batch 3 — Vectorization refactors (pure, no API change)**
5. C7a — vectorize `.somalign_reference_top_labels`
6. C7d — vectorize `.somalign_transfer_labels` (apply with warning fix incorporated)
7. B7b — vectorize `.somalign_distance_quantiles`
8. B7c — vectorize `.somalign_label_probabilities` (apply with dead-code fix incorporated)

**Batch 4 — Chunked projection (API-additive)**
9. B3 utils — add `.somalign_nearest_code_chunked`
10. B3 fit private — update `.somalign_project_samples`
11. B3 fit public — add `chunk_size` to `somalign_fit()` signature
12. B3 threading — thread `chunk_size` through both call sites
13. Run `devtools::document()` to regenerate NAMESPACE and man pages

**Batch 5 — Tests for Batch 3 + 4**
14. D5a — `tests/testthat/test-chunked-projection.R` (apply with test assertion fix)

**Batch 6 — Diagnostics and documentation**
15. C4 — message on `match_mass_ratio > 1`
16. C6 — `@details` roxygen section (apply with old_snippet fix and prose fix)
17. Run `devtools::document()`

**Batch 7 — Test for parallel branch**
18. D8 — add parallel test to `test-training-integration.R`

**Batch 8 — LOW (cleanup)**
19. A9 — scoped POT warning filter (apply with scope fix)
20. E13 guard — add `requireNamespace` guard (partially applied; see section)

---

## 3. Patch details

---

### Patch A1a: Initialise `delta` before the Sinkhorn loop (HIGH)

**File:** `R/ot.R`, lines 84–87

**Concern (suggestion — not a blocker):** The original rationale claims R
for-loops create a new scope (they do not). The real risk is the empty-loop case
(`max_iter = 0`) where the loop body never runs and `delta` would be undefined.
The patch is correct regardless; only the rationale is imprecise.

```r
# OLD
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  iterations <- max_iter
  for (iter in seq_len(max_iter)) {
```

```r
# NEW  (guards the max_iter = 0 empty-loop edge case)
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  delta <- Inf
  iterations <- max_iter
  for (iter in seq_len(max_iter)) {
```

No accompanying test or doc change required.

---

### Patch A1b: Non-convergence warning after Sinkhorn loop (HIGH)

**File:** `R/ot.R`, line 106

⚠️ **BLOCKER 1 — sprintf format mismatch:** The draft uses `"%d"` for
`max_iter`, but `max_iter` defaults to `1000` (a bare numeric/double in R, not
an integer literal). `sprintf("%d", 1000)` throws:
`invalid format "%d"; use format %f, %e, %g or %s for numeric objects`.
Fix: use `"%g"` (or cast `as.integer(max_iter)` explicitly).

⚠️ **BLOCKER 2 — false-positive on last-iteration convergence:** The guard
`iterations == max_iter` is true even when the solver converges on exactly the
final iteration (`iter == max_iter`), because `iterations` is not updated by the
`break` in that case — wait, actually looking at the code: `iterations <- iter;
break` *is* executed, so `iterations` would be `max_iter`, not `max_iter`.
The safe guard is `iterations == max_iter && delta >= tol`. This requires A1a
(`delta <- Inf` initialisation) to already be in place, which it will be.

```r
# OLD
  plan <- sweep(sweep(k, 1, u, "*"), 2, v, "*")
```

```r
# NEW
  if (iterations == max_iter && delta >= tol) {
    warning(
      sprintf(
        "Sinkhorn solver did not converge after %g iterations (final delta = %.3e). ",
        max_iter, delta
      ),
      "Consider increasing max_iter, raising epsilon, or reducing rho_query / rho_ref.",
      call. = FALSE
    )
  }

  plan <- sweep(sweep(k, 1, u, "*"), 2, v, "*")
```

**Required test additions:** See patch D5b.

---

### Patch A2: Kernel underflow fraction check (HIGH)

**File:** `R/ot.R`, lines 78–80

No adversarial blockers. The old_snippet matches exactly.

```r
# OLD
  tiny <- .Machine$double.xmin
  k <- exp(-cost / epsilon)
  k <- pmax(k, tiny)
```

```r
# NEW
  tiny <- .Machine$double.xmin
  k_raw <- exp(-cost / epsilon)
  underflow_fraction <- sum(k_raw < tiny) / length(k_raw)
  if (underflow_fraction > 0.01) {
    safe_eps <- signif(-max(cost) / log(.Machine$double.xmin), 3)
    warning(
      sprintf(
        "%.1f%% of Sinkhorn kernel entries underflowed (epsilon = %g). ",
        100 * underflow_fraction, epsilon
      ),
      sprintf(
        "Raise epsilon or reduce cost scale. Safe lower bound for epsilon: %g",
        safe_eps
      ),
      call. = FALSE
    )
  }
  k <- pmax(k_raw, tiny)
```

**Required test additions:** See patch D5c.

---

### Patches D5b + D5c: Tests for A1b and A2 (MEDIUM)

**File:** `tests/testthat/test-ot-warnings.R` (new file)

Create this file. It has no dependency on other pending patches.

```r
# NEW FILE: tests/testthat/test-ot-warnings.R

test_that("internal solver warns on non-convergence", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  expect_warning(
    somalign_fit(
      query_obj, ref,
      solver = "internal",
      epsilon = 0.1,
      max_iter = 1L,   # force non-convergence: integer to avoid double-format issues
      tol = 0
    ),
    "did not converge"
  )
})

test_that("internal solver does NOT warn when it converges on last iteration", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # With default tol (1e-7) and sufficient iterations, no warning expected
  expect_no_warning(
    somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1)
  )
})

test_that("internal solver warns when kernel underflows", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  # epsilon = 1e-300 forces all kernel entries to underflow
  expect_warning(
    somalign_fit(
      query_obj, ref,
      solver = "internal",
      epsilon = 1e-300
    ),
    "underflowed"
  )
})
```

**Note:** `max_iter = 1L` (integer literal) is used rather than `1` (double) to
avoid any confusion with the now-fixed `%g` format, and to ensure
`iterations == max_iter` triggers even if the solver takes one step that does
not converge. Use `tol = 0` to guarantee the single-iteration case never passes
the convergence check.

---

### Patch C7a: Vectorize `.somalign_reference_top_labels` (MEDIUM)

**File:** `R/results.R`, lines 80–91

No adversarial concerns. Old_snippet matches exactly.

```r
# OLD
  label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  label_names <- colnames(label_prob)
  for (i in seq_len(n_nodes)) {
    row <- label_prob[i, ]
    if (sum(row) > 0) {
      idx <- which.max(row)
      label[i] <- label_names[idx]
      confidence[i] <- row[idx]
    }
  }
  list(label = label, confidence = confidence)
```

```r
# NEW
  label_names <- colnames(label_prob)
  row_sums <- rowSums(label_prob)
  has_mass <- row_sums > 0
  idx <- max.col(label_prob, ties.method = "first")
  label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  label[has_mass] <- label_names[idx[has_mass]]
  confidence[has_mass] <- label_prob[cbind(which(has_mass), idx[has_mass])]
  list(label = label, confidence = confidence)
```

No accompanying doc changes required (internal function).

---

### Patch C7d: Vectorize `.somalign_transfer_labels` (MEDIUM)

**File:** `R/fit.R`, lines 161–182

⚠️ **WARNING — second_label divergence in single-dominant-label nodes:** The
original for-loop assigns `second_label[i] <- label_names[ord[2]]` even when
`row[ord[2]] == 0`, so `second_label` gets a label name with
`second_confidence == 0`. The vectorised replacement uses `max.col()` on a row
where the top column has been zeroed; when all remaining columns are also zero,
`max.col` returns column 1, making `second_label` equal to the top label. The
review recommendation is to explicitly set `second_label` to `NA_character_`
when `second_confidence == 0`, which is cleaner than the loop's behaviour. This
semantic change is intentional and must be documented in the commit message.

```r
# OLD
  probs <- correspondence %*% label_prob
  label_names <- colnames(label_prob)
  top_label <- rep(NA_character_, n_nodes)
  second_label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  second_confidence <- rep(NA_real_, n_nodes)
  entropy <- rep(NA_real_, n_nodes)

  for (i in seq_len(n_nodes)) {
    row <- as.numeric(probs[i, ])
    if (sum(row) > 0) {
      row <- row / sum(row)
      ord <- order(row, decreasing = TRUE)
      top_label[i] <- label_names[ord[1]]
      confidence[i] <- row[ord[1]]
      if (length(ord) > 1) {
        second_label[i] <- label_names[ord[2]]
        second_confidence[i] <- row[ord[2]]
      }
      entropy[i] <- .somalign_entropy(row)
    }
  }
```

```r
# NEW  (semantic change: second_label is NA when second_confidence == 0)
  probs <- correspondence %*% label_prob
  label_names <- colnames(label_prob)
  row_sums <- rowSums(probs)
  has_mass <- row_sums > 0
  probs_norm <- probs
  probs_norm[has_mass, ] <- probs[has_mass, , drop = FALSE] / row_sums[has_mass]
  top_idx <- max.col(probs_norm, ties.method = "first")
  top_label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  top_label[has_mass] <- label_names[top_idx[has_mass]]
  confidence[has_mass] <- probs_norm[cbind(which(has_mass), top_idx[has_mass])]
  probs_second <- probs_norm
  probs_second[cbind(seq_len(n_nodes), top_idx)] <- 0
  second_idx <- max.col(probs_second, ties.method = "first")
  second_label <- rep(NA_character_, n_nodes)
  second_confidence <- rep(NA_real_, n_nodes)
  if (ncol(probs) == 1L) {
    second_label <- rep(NA_character_, n_nodes)
    second_confidence <- rep(NA_real_, n_nodes)
  } else {
    second_label[has_mass] <- label_names[second_idx[has_mass]]
    second_confidence[has_mass] <- probs_second[cbind(which(has_mass), second_idx[has_mass])]
    # Cleaner than loop: set second_label to NA when there is no real second choice
    second_label[has_mass & (is.na(second_confidence) | second_confidence == 0)] <- NA_character_
  }
  entropy <- vapply(
    seq_len(n_nodes),
    function(i) if (has_mass[i]) .somalign_entropy(probs_norm[i, ]) else NA_real_,
    numeric(1)
  )
```

No accompanying doc changes required (internal function). Tests for the
`ncol == 1` guard and `second_confidence` correctness are covered by the
existing `test-ot-labels-results.R` label transfer test. A more targeted test
is recommended but not a blocker.

---

### Patch B7b: Vectorize `.somalign_distance_quantiles` (MEDIUM)

**File:** `R/utils.R`, lines 241–256

No adversarial blockers. Old_snippet matches exactly.

```r
# OLD
.somalign_distance_quantiles <- function(distances, units, n_nodes, probs) {
  names <- .somalign_quantile_names(probs)
  global <- stats::quantile(distances, probs = probs, names = FALSE, type = 7)
  names(global) <- names
  out <- matrix(NA_real_, nrow = n_nodes, ncol = length(probs))
  colnames(out) <- names
  for (i in seq_len(n_nodes)) {
    node_distances <- distances[units == i]
    if (length(node_distances) == 0) {
      out[i, ] <- global
    } else {
      out[i, ] <- stats::quantile(node_distances, probs = probs, names = FALSE, type = 7)
    }
  }
  list(node = out, global = global)
}
```

```r
# NEW
.somalign_distance_quantiles <- function(distances, units, n_nodes, probs) {
  names <- .somalign_quantile_names(probs)
  global <- stats::quantile(distances, probs = probs, names = FALSE, type = 7)
  names(global) <- names
  split_distances <- split(distances, factor(units, levels = seq_len(n_nodes)))
  node_rows <- lapply(split_distances, function(d) {
    if (length(d) == 0L) {
      global
    } else {
      stats::quantile(d, probs = probs, names = FALSE, type = 7)
    }
  })
  out <- matrix(unlist(node_rows, use.names = FALSE), nrow = n_nodes, ncol = length(probs), byrow = TRUE)
  colnames(out) <- names
  list(node = out, global = global)
}
```

No accompanying doc changes required (internal function).

---

### Patch B7c: Vectorize `.somalign_label_probabilities` (MEDIUM)

**File:** `R/utils.R`, lines 276–299

⚠️ **WARNING — dead code in draft patch:** The draft adds an `if (!any(valid))`
early-return block after computing `levels`. This block can never be reached
because the `if (length(levels) == 0)` guard above it already returns early
whenever all labels are NA (which makes `levels` empty). The dead block must be
removed to avoid misleading future readers.

Apply the patch below, which omits the dead code:

```r
# OLD
.somalign_label_probabilities <- function(labels, units, n_nodes) {
  if (is.null(labels)) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  if (length(labels) != length(units)) {
    stop("`labels` must have one value per row of reference data.", call. = FALSE)
  }
  labels <- as.character(labels)
  levels <- sort(unique(labels[!is.na(labels)]))
  if (length(levels) == 0) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  out <- matrix(0, nrow = n_nodes, ncol = length(levels))
  colnames(out) <- levels
  for (node in seq_len(n_nodes)) {
    node_labels <- labels[units == node]
    node_labels <- node_labels[!is.na(node_labels)]
    if (length(node_labels) > 0) {
      counts <- table(factor(node_labels, levels = levels))
      out[node, ] <- as.numeric(counts) / sum(counts)
    }
  }
  out
}
```

```r
# NEW  (dead code removed; vectorized with column-major linear indexing)
.somalign_label_probabilities <- function(labels, units, n_nodes) {
  if (is.null(labels)) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  if (length(labels) != length(units)) {
    stop("`labels` must have one value per row of reference data.", call. = FALSE)
  }
  labels <- as.character(labels)
  levels <- sort(unique(labels[!is.na(labels)]))
  if (length(levels) == 0) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  valid <- !is.na(labels)
  node_idx <- units[valid]
  lbl_idx <- as.integer(factor(labels[valid], levels = levels))
  combined_idx <- (lbl_idx - 1L) * n_nodes + node_idx
  raw_counts <- tabulate(combined_idx, nbins = n_nodes * length(levels))
  out <- matrix(raw_counts, nrow = n_nodes, ncol = length(levels))
  colnames(out) <- levels
  row_totals <- rowSums(out)
  nonzero <- row_totals > 0
  out[nonzero, , drop = FALSE] <- out[nonzero, , drop = FALSE] / row_totals[nonzero]
  out
}
```

No accompanying doc changes required (internal function).

---

### Patch B3 (utils): Add `.somalign_nearest_code_chunked` (MEDIUM)

**File:** `R/utils.R`, after line 206 (after the existing
`.somalign_nearest_code` function)

No adversarial blockers.

```r
# OLD  (end of .somalign_nearest_code — insert NEW immediately after)
.somalign_nearest_code <- function(x, codebook) {
  x <- as.matrix(x)
  codebook <- as.matrix(codebook)
  d2 <- outer(rowSums(x * x), rowSums(codebook * codebook), "+") -
    2 * tcrossprod(x, codebook)
  d2 <- pmax(d2, 0)
  unit <- max.col(-d2, ties.method = "first")
  distance <- sqrt(d2[cbind(seq_len(nrow(d2)), unit)])
  list(unit = as.integer(unit), distance = as.numeric(distance))
}
```

```r
# NEW  (append after .somalign_nearest_code, before .somalign_pairwise_distance)
.somalign_nearest_code_chunked <- function(x, codebook, chunk_size = 10000L) {
  x <- as.matrix(x)
  n <- nrow(x)
  # Short-circuit: Inf, NULL, or chunk_size >= n all use the single-call path
  if (is.null(chunk_size) || is.infinite(chunk_size) || n == 0L || chunk_size >= n) {
    return(.somalign_nearest_code(x, codebook))
  }
  chunk_size <- as.integer(chunk_size)
  unit <- integer(n)
  distance <- numeric(n)
  for (s in seq(1L, n, by = chunk_size)) {
    idx <- s:min(s + chunk_size - 1L, n)
    res <- .somalign_nearest_code(x[idx, , drop = FALSE], codebook)
    unit[idx] <- res$unit
    distance[idx] <- res$distance
  }
  list(unit = unit, distance = distance)
}
```

---

### Patch B3 (fit private): `.somalign_project_samples` calls chunked function (MEDIUM)

**File:** `R/fit.R`, lines 222–231

```r
# OLD
.somalign_project_samples <- function(scaled_data, reference) {
  projected <- .somalign_nearest_code(scaled_data, reference$codebook)
  threshold <- .somalign_thresholds(reference, projected$unit)
  list(
    unit = projected$unit,
    distance = projected$distance,
    threshold = threshold,
    outside = projected$distance > threshold
  )
}
```

```r
# NEW
.somalign_project_samples <- function(scaled_data, reference, chunk_size = 10000L) {
  projected <- .somalign_nearest_code_chunked(scaled_data, reference$codebook, chunk_size = chunk_size)
  threshold <- .somalign_thresholds(reference, projected$unit)
  list(
    unit = projected$unit,
    distance = projected$distance,
    threshold = threshold,
    outside = projected$distance > threshold
  )
}
```

---

### Patch B3 (fit public): Add `chunk_size` parameter to `somalign_fit()` (MEDIUM)

**File:** `R/fit.R`, lines 30–40

Also add the `@param` roxygen tag. Insert the tag just before the existing
`@return` line (or the `@details` block once C6 is applied).

```r
# OLD  (function signature)
somalign_fit <- function(query,
                         reference,
                         epsilon = 0.05,
                         rho_query = 1,
                         rho_ref = 1,
                         solver = c("auto", "pot", "internal"),
                         min_match_fraction = 0.05,
                         confidence_threshold = 0.6,
                         correction_min_mass = 1e-8,
                         max_iter = 1000,
                         tol = 1e-7) {
```

```r
# NEW
somalign_fit <- function(query,
                         reference,
                         epsilon = 0.05,
                         rho_query = 1,
                         rho_ref = 1,
                         solver = c("auto", "pot", "internal"),
                         min_match_fraction = 0.05,
                         confidence_threshold = 0.6,
                         correction_min_mass = 1e-8,
                         max_iter = 1000,
                         tol = 1e-7,
                         chunk_size = 10000L) {
```

Add the following roxygen tag in the parameter block of `somalign_fit()` (before
`@return`):

```r
#' @param chunk_size Integer. Number of samples to project per chunk when
#'   computing nearest reference node. Use `Inf` or `NULL` for no chunking
#'   (allocates a full n_samples x n_nodes matrix). Default `10000L`.
```

---

### Patch B3 (threading): Thread `chunk_size` through both projection call sites (MEDIUM)

**File:** `R/fit.R`, lines 82–84

```r
# OLD
  direct <- .somalign_project_samples(query$scaled_data, reference)
  corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]
  corrected <- .somalign_project_samples(corrected_matrix, reference)
```

```r
# NEW
  direct <- .somalign_project_samples(query$scaled_data, reference, chunk_size = chunk_size)
  corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]
  corrected <- .somalign_project_samples(corrected_matrix, reference, chunk_size = chunk_size)
```

After applying all B3 sub-patches, run `devtools::document()` to regenerate man
pages.

---

### Patch D5a: Tests for B3/B7b/B7c (MEDIUM)

**File:** `tests/testthat/test-chunked-projection.R` (new file)

⚠️ **BLOCKER — test assertion for node 2 is arithmetically wrong in the draft:**
The draft asserts `res[2, "A"] == 1/3` and `res[2, "B"] == 2/3` for node 2.
With `labels = c("A","B",NA,"A","B","C","A",NA,"C","B")` and
`units = c(1L,1L,2L,2L,3L,3L,4L,4L,1L,2L)`, the observations assigned to
unit 2 are positions 3 (NA, excluded), 4 (label "A"), and 10 (label "B").
That gives A=1, B=1, so both probabilities are 0.5. The draft comment
incorrectly counted position 2 (which has unit=1) as belonging to node 2.

Apply the corrected test below:

```r
# NEW FILE: tests/testthat/test-chunked-projection.R

test_that(".somalign_nearest_code_chunked matches unchunked for various chunk_sizes", {
  set.seed(42L)
  x <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  codebook <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)

  full          <- somalign:::.somalign_nearest_code(x, codebook)
  chunked3      <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = 3L)
  chunked_inf   <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = Inf)
  chunked_null  <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = NULL)
  chunked_big   <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = 10000L)

  expect_identical(chunked3$unit,     full$unit)
  expect_equal(   chunked3$distance,  full$distance)
  expect_identical(chunked_inf$unit,  full$unit)
  expect_identical(chunked_null$unit, full$unit)
  expect_identical(chunked_big$unit,  full$unit)
})

test_that(".somalign_distance_quantiles vectorised output matches loop output on data with empty nodes", {
  set.seed(7L)
  n_nodes   <- 5L
  distances <- abs(rnorm(20))
  # units 1-4 populated, node 5 empty
  units <- sample(1:4, 20, replace = TRUE)
  probs <- c(0.5, 0.9, 0.95, 0.99)

  res <- somalign:::.somalign_distance_quantiles(distances, units, n_nodes, probs)
  expect_equal(nrow(res$node), n_nodes)
  expect_equal(ncol(res$node), length(probs))
  # Empty node (5) should equal global quantile
  expect_equal(res$node[5, ], res$global)
  expect_null(rownames(res$node))
})

test_that(".somalign_label_probabilities vectorised output matches loop output", {
  set.seed(3L)
  n_nodes <- 4L
  labels  <- c("A", "B", NA,  "A", "B", "C", "A", NA,  "C", "B")
  units   <- c(1L,  1L,  2L,  2L,  3L,  3L,  4L,  4L,  1L,  2L)

  res <- somalign:::.somalign_label_probabilities(labels, units, n_nodes)
  expect_equal(nrow(res), n_nodes)
  expect_equal(colnames(res), c("A", "B", "C"))
  # All rows should sum to 1 (no empty nodes here)
  expect_equal(rowSums(res), rep(1, n_nodes))
  # Node 2: positions with unit==2 are indices 3 (NA, excluded), 4 (A), 10 (B)
  # -> A=1, B=1, so each is 0.5
  expect_equal(res[2, "A"], 0.5)
  expect_equal(res[2, "B"], 0.5)
  expect_equal(res[2, "C"], 0)
})
```

---

### Patch C4: Message when `match_mass_ratio > 1` (MEDIUM)

**File:** `R/fit.R`, lines 61–62

No adversarial blockers. Old_snippet matches exactly.

```r
# OLD
  match_mass_ratio <- ifelse(query$node_masses > 0, row_mass / query$node_masses, 0)
  match_fraction <- pmin(match_mass_ratio, 1)
```

```r
# NEW
  match_mass_ratio <- ifelse(query$node_masses > 0, row_mass / query$node_masses, 0)
  match_fraction <- pmin(match_mass_ratio, 1)
  n_over <- sum(match_mass_ratio > 1)
  if (n_over > 0) {
    message(sprintf(
      "somalign_fit: %d query node(s) have match_mass_ratio > 1 (max %.2f); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.",
      n_over, max(match_mass_ratio)
    ))
  }
```

No accompanying test or doc changes required beyond the roxygen detail in C6.

---

### Patch C6: `@details` roxygen section for `somalign_fit()` (MEDIUM)

**File:** `R/fit.R`

⚠️ **BLOCKER — old_snippet does not match the file:** The draft assumes
`#' @return A \`somalign_fit\` object.\n#' @export\nsomalign_fit <- function(query,`
is a contiguous string, but in the current `R/fit.R` (lines 19–30) the `@return`
tag (line 19) is followed by an `@examples` block (lines 20–28) before `@export`
(line 29). The string replacement will fail.

⚠️ **WARNING — parameter name inconsistency in draft prose:** The draft uses
bare `` `rho` `` in two places where the actual parameter names are `rho_query`
and `rho_ref`.

Use this corrected old_snippet and new_snippet:

```r
# OLD  (anchor on @return + beginning of @examples block)
#' @return A `somalign_fit` object.
#' @examples
#' \dontrun{
```

```r
# NEW  (insert @details before @return, fix parameter names throughout)
#' @details
#' The transport plan row sums will not equal `query$node_masses` exactly — this
#' is by design. Unbalanced optimal transport allows mass destruction, so some
#' query mass may be absorbed rather than transported. Deviation grows with lower
#' `rho_query` / `rho_ref` values and higher `epsilon`. At the defaults
#' (`rho_query = 1`, `rho_ref = 1`, `epsilon = 0.05`), row-sum deviation can
#' reach approximately 13%. Use `diagnostics$ot$max_row_mass_error` to quantify
#' the deviation in a given fit; for near-balanced data, increase `rho_query`
#' (e.g. `rho_query = 10`) to enforce tighter marginal constraints.
#'
#' @return A `somalign_fit` object.
#' @examples
#' \dontrun{
```

After applying, run `devtools::document()` to regenerate the man page.

---

### Patch D8: Test for `parallel = TRUE` in `somalign_sensitivity_grid()` (MEDIUM)

**File:** `tests/testthat/test-training-integration.R`

Append after the last existing test in the file. The test calls
`somalign_sensitivity_grid` with `parallel = TRUE` and verifies the result has
the same structure as the sequential call.

```r
# NEW  (append to tests/testthat/test-training-integration.R)

test_that("somalign_sensitivity_grid parallel = TRUE returns same structure as sequential", {
  skip_if_not_installed("kohonen")

  set.seed(99L)
  old <- rbind(
    matrix(rnorm(10 * 4, mean = -1), ncol = 4),
    matrix(rnorm(10 * 4, mean =  1), ncol = 4)
  )
  colnames(old) <- paste0("f", seq_len(ncol(old)))

  reference <- somalign_train_reference(
    old,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )
  query_obj <- somalign_query(
    old,
    reference,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )

  seq_result  <- somalign_sensitivity_grid(
    query_obj, reference,
    epsilon   = c(0.05, 0.1),
    rho_query = 1,
    rho_ref   = 1,
    solver    = "internal",
    parallel  = FALSE
  )
  par_result  <- somalign_sensitivity_grid(
    query_obj, reference,
    epsilon   = c(0.05, 0.1),
    rho_query = 1,
    rho_ref   = 1,
    solver    = "internal",
    parallel  = TRUE
  )

  expect_equal(nrow(par_result), nrow(seq_result))
  expect_equal(names(par_result), names(seq_result))
  expect_equal(par_result$epsilon, seq_result$epsilon)
})
```

---

### Patch A9: Suppress benign POT `reg_type='entropy'` UserWarning (LOW)

**File:** `R/ot.R`, lines 54–61

⚠️ **WARNING — filter scope is too broad:** The draft uses
`warnings.filterwarnings('ignore', category=UserWarning, module='ot')`, which
suppresses ALL `UserWarning` from POT for the entire Python session, including
POT's own Sinkhorn non-convergence warnings. This undermines the intent of A1b.
Additionally, the filter is added to global Python state on every call, which is
messy and redundant.

Apply a scoped, message-specific filter instead:

```r
# OLD
  plan <- fn(
    a = a,
    b = b,
    M = cost,
    reg = epsilon,
    reg_m = c(rho_query, rho_ref),
    reg_type = "entropy"
  )
```

```r
# NEW  (scoped, message-specific warning filter; does not silence POT convergence warnings)
  reticulate::py_run_string(
    "import warnings; warnings.filterwarnings('ignore', message='.*variable c.*', category=UserWarning, module='ot')"
  )
  plan <- fn(
    a = a,
    b = b,
    M = cost,
    reg = epsilon,
    reg_m = c(rho_query, rho_ref),
    reg_type = "entropy"
  )
```

**Note:** `py_run_string` sets a persistent Python-session filter, but the
message pattern `'.*variable c.*'` narrows suppression to only the specific
harmless "overwrites variable c" warning. This is acceptable. If stronger
isolation is desired in a future refactor, consider wrapping the call in a
Python context manager using `reticulate::py_eval`.

---

### Patch E13 guard: Add `requireNamespace` guard for `parallel` (LOW)

**File:** `R/diagnostics.R`, `somalign_sensitivity_grid()` function

**Status:** ALREADY APPLIED in `somalign_sensitivity_grid()` body — the
`parallel = FALSE` default and the `if (isTRUE(parallel))` dispatch are in
place. However, per CRAN policy, `parallel` is listed in `Suggests` and must
have a runtime availability guard when `parallel = TRUE` is used. The current
code calls `parallel::mclapply()` unconditionally when the branch is reached
without checking whether the package is available.

Add a guard at the top of the `if (isTRUE(parallel))` branch:

```r
# OLD
  if (isTRUE(parallel)) {
    rows <- parallel::mclapply(
      seq_len(nrow(grid)),
      .run_one,
      mc.cores = getOption("mc.cores", 1L)
    )
```

```r
# NEW
  if (isTRUE(parallel)) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop(
        "Package 'parallel' is required when parallel = TRUE. ",
        "Install it or set parallel = FALSE.",
        call. = FALSE
      )
    }
    rows <- parallel::mclapply(
      seq_len(nrow(grid)),
      .run_one,
      mc.cores = getOption("mc.cores", 1L)
    )
```

**Alternative:** Move `parallel` from `Suggests` to `Imports` in DESCRIPTION
(simpler, since `parallel` ships with every R installation). If you prefer this
route, also add `@importFrom parallel mclapply` in the roxygen block for
`somalign_sensitivity_grid()` and run `devtools::document()`.

---

### Patches E10, E11, E12, E14: Already applied

These patches are already in the working tree as confirmed by reading the files:
- `@examples \dontrun{}` blocks are present in all public functions.
- DESCRIPTION has the placeholder author with ORCID field.
- DESCRIPTION Suggests includes `bench`, `microbenchmark`, `parallel`.
- NEWS.md has the full version-0.0.0.9000 changelog.

No action required for these patches.

---

## 4. How to apply

**Tooling setup:** Keep a terminal with `devtools::load_all()` ready. All
patches touch internal R code and roxygen comments; no C/Fortran compilation is
needed.

**Batch workflow:**

1. Apply Batch 1 patches (A1a, A1b, A2) as three separate, focused edits to
   `R/ot.R`. Commit as one unit: `"fix: warn on Sinkhorn non-convergence and
   kernel underflow (A1a, A1b, A2)"`.

2. Create `tests/testthat/test-ot-warnings.R` (D5b + D5c). Run
   `devtools::test(filter = "ot-warnings")`. Commit alongside the HIGH fixes
   or as an immediate follow-up.

3. Apply Batch 3 vectorisation patches (C7a, C7d, B7b, B7c) as edits to
   `R/results.R`, `R/fit.R`, and `R/utils.R`. These are pure refactors with no
   API change. Run `devtools::test()` to confirm no regressions. Commit as one
   unit: `"perf: vectorize label and distance helpers (C7a, C7d, B7b, B7c)"`.

4. Apply all four B3 sub-patches to `R/utils.R` and `R/fit.R`. Run
   `devtools::document()` to pick up the new `@param chunk_size` roxygen tag.
   Run `devtools::test()`. Commit: `"feat: chunked nearest-node projection to
   cap peak memory (B3)"`.

5. Create `tests/testthat/test-chunked-projection.R` (D5a, with corrected test
   assertions). Run `devtools::test(filter = "chunked")`. Commit alongside B3.

6. Apply C4 and C6 to `R/fit.R`. For C6, use the corrected old_snippet that
   anchors on `@return ... @examples \dontrun{` rather than `@return ... @export`.
   Run `devtools::document()`. Run `devtools::test()`. Commit: `"docs: add
   @details for marginal deviation; add diagnostic message for ratio>1 (C4, C6)"`.

7. Append the D8 parallel test to `test-training-integration.R`. Run
   `devtools::test(filter = "training")`. Commit with C4/C6 or separately.

8. Apply A9 (scoped POT filter) and the E13 `requireNamespace` guard as
   standalone LOW fixes. These can be a single commit: `"fix: scope POT warning
   filter; guard parallel::mclapply availability (A9, E13)"`.

**After all batches:** Run `devtools::check()` to confirm no new NOTES or
WARNINGs. Pay particular attention to: (a) the ORCID placeholder in DESCRIPTION
must be replaced with a real value before CRAN submission; (b) the `parallel`
Suggests-without-guard NOTE should be resolved by the E13 guard or by moving
`parallel` to Imports.
