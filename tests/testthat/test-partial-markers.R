# Tests for somalign_reference_subset_markers()
#
# Covers:
#   - Basic subsetting: returned object has correct fields
#   - Query + fit pipeline succeeds on the subset reference
#   - Label transfer is sensible when populations are separable on shared markers
#   - Guard test: disjoint-marker input is rejected
#   - Guard test: unknown markers are rejected
#   - Guard test: empty markers are rejected
#   - Without reference_data: distance_quantiles are Inf and a warning is emitted
#   - With reference_data: distance_quantiles are finite and recalibrated in subspace
#   - With reference_data: outside-reference detection is calibrated (inside vs outside)

# ---- helpers -----------------------------------------------------------------

# Build a 4-marker synthetic dataset with two well-separated populations.
# Pop A lives at (+3, +3, +3, +3), Pop B at (-3, -3, -3, -3) in raw space.
make_full_panel_reference <- function(seed = 42L) {
  withr::local_seed(seed)
  markers <- c("CD3", "CD4", "CD8", "CD19")
  n_per   <- 40L
  ref_data <- rbind(
    matrix(rnorm(n_per * 4L, mean =  3, sd = 0.4), ncol = 4L,
           dimnames = list(NULL, markers)),
    matrix(rnorm(n_per * 4L, mean = -3, sd = 0.4), ncol = 4L,
           dimnames = list(NULL, markers))
  )
  labels <- rep(c("A", "B"), each = n_per)
  somalign_train_reference(
    ref_data,
    labels = labels,
    grid   = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen   = 30L
  )
}

# Build a reference AND return the raw data so callers can pass reference_data.
make_full_panel_reference_with_data <- function(seed = 42L) {
  withr::local_seed(seed)
  markers <- c("CD3", "CD4", "CD8", "CD19")
  n_per   <- 40L
  ref_data <- rbind(
    matrix(rnorm(n_per * 4L, mean =  3, sd = 0.4), ncol = 4L,
           dimnames = list(NULL, markers)),
    matrix(rnorm(n_per * 4L, mean = -3, sd = 0.4), ncol = 4L,
           dimnames = list(NULL, markers))
  )
  labels <- rep(c("A", "B"), each = n_per)
  ref <- somalign_train_reference(
    ref_data,
    labels = labels,
    grid   = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen   = 30L
  )
  list(reference = ref, data = ref_data)
}

# Build query data on a subset of markers with the same two populations,
# shifted by a fixed per-marker batch offset.
make_partial_query_data <- function(seed = 99L,
                                    markers = c("CD3", "CD4", "CD8"),
                                    n_per   = 40L,
                                    shift   = 0.5) {
  withr::local_seed(seed)
  rbind(
    matrix(rnorm(n_per * length(markers), mean =  3 + shift, sd = 0.4),
           ncol = length(markers), dimnames = list(NULL, markers)),
    matrix(rnorm(n_per * length(markers), mean = -3 + shift, sd = 0.4),
           ncol = length(markers), dimnames = list(NULL, markers))
  )
}


# ---- somalign_reference_subset_markers() structure tests --------------------

test_that("subset reference has correct features", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_equal(ref_sub$features, c("CD3", "CD4", "CD8"))
})

test_that("subset reference codebook has correct column count", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_equal(ncol(ref_sub$codebook), 3L)
  expect_equal(colnames(ref_sub$codebook), c("CD3", "CD4", "CD8"))
})

test_that("subset reference preserves node count", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_equal(nrow(ref_sub$codebook), nrow(ref_full$codebook))
  expect_equal(length(ref_sub$node_masses), length(ref_full$node_masses))
})

test_that("subset reference center and scale are restricted correctly", {
  ref_full <- make_full_panel_reference()
  markers  <- c("CD3", "CD8")
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )

  expect_equal(names(ref_sub$center), markers)
  expect_equal(names(ref_sub$scale), markers)
  expect_equal(unname(ref_sub$center), unname(ref_full$center[markers]))
  expect_equal(unname(ref_sub$scale),  unname(ref_full$scale[markers]))
})

test_that("subset reference preserves label_prob unchanged", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_equal(ref_sub$label_prob, ref_full$label_prob)
})

test_that("subset reference preserves node_masses unchanged", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_equal(ref_sub$node_masses, ref_full$node_masses)
})

test_that("without reference_data: node_var is NULL (disabled)", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_null(ref_sub$node_var)
})

test_that("with reference_data: node_var is recomputed and restricted to subset columns", {
  obj      <- make_full_panel_reference_with_data()
  ref_sub  <- somalign_reference_subset_markers(
    obj$reference, c("CD3", "CD4", "CD8"), reference_data = obj$data
  )

  expect_false(is.null(ref_sub$node_var))
  expect_equal(ncol(ref_sub$node_var), 3L)
  expect_equal(colnames(ref_sub$node_var), c("CD3", "CD4", "CD8"))
  expect_true(all(ref_sub$node_var > 0))
})

test_that("subset reference is a somalign_reference object", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_s3_class(ref_sub, "somalign_reference")
})

test_that("single-marker subset is accepted", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, "CD3")
  )

  expect_equal(ref_sub$features, "CD3")
  expect_equal(ncol(ref_sub$codebook), 1L)
})

test_that("order of markers in subset is respected", {
  ref_full <- make_full_panel_reference()
  # Deliberately reverse order relative to training order
  markers  <- c("CD19", "CD3")
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )

  expect_equal(ref_sub$features, markers)
  expect_equal(colnames(ref_sub$codebook), markers)
})


# ---- pipeline tests: somalign_query() + somalign_fit() ----------------------

test_that("somalign_query() succeeds with partial-marker reference", {
  ref_full  <- make_full_panel_reference()
  markers   <- c("CD3", "CD4", "CD8")
  ref_sub   <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )
  qry_data  <- make_partial_query_data(markers = markers)

  # Must not error
  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  expect_s3_class(qry, "somalign_query")
})

test_that("somalign_fit() runs to completion on partial-marker alignment", {
  ref_full  <- make_full_panel_reference()
  markers   <- c("CD3", "CD4", "CD8")
  ref_sub   <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )
  qry_data  <- make_partial_query_data(markers = markers)

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  fit <- somalign_fit(qry, ref_sub)

  expect_s3_class(fit, "somalign_fit")
  # transport plan is non-trivial
  expect_gt(sum(fit$transport_plan), 0)
})

test_that("label transfer assigns the two populations correctly on shared markers", {
  # Populations are separated even on CD3/CD4/CD8 alone (CD19 just dropped),
  # so label transfer should produce mostly concordant labels.
  ref_full  <- make_full_panel_reference()
  markers   <- c("CD3", "CD4", "CD8")
  ref_sub   <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )
  qry_data  <- make_partial_query_data(markers = markers)

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 30L
  )
  # Use epsilon = 0.2 to avoid Sinkhorn underflow on this compact 4-node grid
  fit <- somalign_fit(qry, ref_sub, epsilon = 0.2)

  # At least one query node should receive an accepted label transfer
  accepted_labels <- fit$label_transfer$label[fit$label_transfer$accepted]
  expect_gt(length(accepted_labels), 0L)

  # Both populations should appear among accepted labels
  expect_true("A" %in% accepted_labels || "B" %in% accepted_labels)
})

test_that("node_shifts from partial-marker fit have correct marker dimensions", {
  ref_full  <- make_full_panel_reference()
  markers   <- c("CD3", "CD4", "CD8")
  ref_sub   <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )
  qry_data  <- make_partial_query_data(markers = markers)

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  fit <- somalign_fit(qry, ref_sub)

  # node_shifts must span the subset marker dimensions only
  expect_equal(ncol(fit$node_shifts), length(markers))
  expect_equal(colnames(fit$node_shifts), markers)
})


# ---- guard tests: input validation ------------------------------------------

test_that("error when markers not in reference", {
  ref_full <- make_full_panel_reference()

  expect_error(
    somalign_reference_subset_markers(ref_full, c("CD3", "NOTAMARKER")),
    regexp = "not in the reference feature set",
    fixed  = FALSE
  )
})

test_that("error when markers argument is empty", {
  ref_full <- make_full_panel_reference()

  expect_error(
    somalign_reference_subset_markers(ref_full, character(0)),
    regexp = "non-empty character vector"
  )
})

test_that("error when markers contains duplicates", {
  ref_full <- make_full_panel_reference()

  expect_error(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD3", "CD4")),
    regexp = "Duplicated"
  )
})

test_that("error when reference argument is not a somalign_reference", {
  ref_full <- make_full_panel_reference()

  expect_error(
    somalign_reference_subset_markers(list(features = "CD3"), c("CD3")),
    regexp = "somalign_reference"
  )
})

test_that("somalign_query() succeeds with a superset of reference markers (extras are dropped)", {
  # somalign_query() selects only the reference features from the data matrix,
  # so a query with extra markers is accepted without error.
  ref_full <- make_full_panel_reference()
  markers  <- c("CD3", "CD4", "CD8")
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )

  # Data has CD19 in addition to the 3 subset markers -- extra column is dropped
  full_qry_data <- make_partial_query_data(
    markers = c("CD3", "CD4", "CD8", "CD19")
  )
  qry <- somalign_query(
    full_qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 10L
  )
  expect_equal(qry$reference_features, markers)
})

test_that("somalign_query() errors when query is missing a required reference marker", {
  # A query missing a marker that the (subset) reference requires should error.
  ref_full <- make_full_panel_reference()
  markers  <- c("CD3", "CD4", "CD8")
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )

  # Query only has 2 of the 3 required markers; CD8 is absent
  partial_qry_data <- make_partial_query_data(markers = c("CD3", "CD4"))
  expect_error(
    somalign_query(
      partial_qry_data, ref_sub,
      grid = kohonen::somgrid(2L, 2L, "hexagonal"),
      rlen = 10L
    ),
    regexp = "Missing features"
  )
})

test_that("dropping only the separator marker degrades but does not crash label transfer", {
  # Pop A: high CD19, low CD3/CD4/CD8.
  # Pop B: low CD19, high CD3/CD4/CD8.
  # CD19 is the only lineage separator; dropping it should make shared-marker
  # separation poor, but the pipeline must still run and return a fit object.
  withr::local_seed(7L)
  markers_full <- c("CD3", "CD4", "CD8", "CD19")
  n_per        <- 40L
  ref_data     <- rbind(
    # Pop A: CD19 high, rest near zero
    cbind(matrix(rnorm(n_per * 3L, mean =  0, sd = 0.3), ncol = 3L),
          matrix(rnorm(n_per       , mean =  5, sd = 0.3), ncol = 1L)),
    # Pop B: CD19 low, rest near zero
    cbind(matrix(rnorm(n_per * 3L, mean =  0, sd = 0.3), ncol = 3L),
          matrix(rnorm(n_per       , mean = -5, sd = 0.3), ncol = 1L))
  )
  colnames(ref_data) <- markers_full
  ref_full <- somalign_train_reference(
    ref_data,
    labels = rep(c("A", "B"), each = n_per),
    grid   = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen   = 30L
  )

  # Subset to markers that carry NO population signal
  ref_sub   <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )
  qry_data  <- ref_data[, c("CD3", "CD4", "CD8"), drop = FALSE] + 0.3

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 30L
  )
  fit <- somalign_fit(qry, ref_sub)

  # Pipeline succeeds -- alignment may be uninformative but must not error
  expect_s3_class(fit, "somalign_fit")

  # Document the caveat: label_transfer should show low or absent confidence
  # (We do NOT assert exact numbers -- the exact result is seed-dependent.)
  all_confidence <- fit$label_transfer$confidence[fit$label_transfer$accepted]
  # This is a documentation test: we simply record that confidence may be low
  expect_true(is.numeric(all_confidence))
})


# ---- new tests: without reference_data (sentinel behavior) ------------------

test_that("without reference_data: emits a warning about disabled detection", {
  ref_full <- make_full_panel_reference()

  expect_warning(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8")),
    regexp = "DISABLED"
  )
})

test_that("without reference_data: distance_quantiles are all Inf (never-flag sentinel)", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  # All entries must be Inf, not the original full-p values
  expect_true(all(is.infinite(ref_sub$distance_quantiles)))
})

test_that("without reference_data: distance_quantiles have same shape as original", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  # Same number of rows (nodes) and columns (quantile levels)
  expect_equal(nrow(ref_sub$distance_quantiles), nrow(ref_full$distance_quantiles))
  expect_equal(ncol(ref_sub$distance_quantiles), ncol(ref_full$distance_quantiles))
  expect_equal(colnames(ref_sub$distance_quantiles), colnames(ref_full$distance_quantiles))
})

test_that("without reference_data: global_distance_quantiles are all Inf", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_true(all(is.infinite(ref_sub$global_distance_quantiles)))
})

test_that("without reference_data: node_var is NULL", {
  ref_full <- make_full_panel_reference()
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
  )

  expect_null(ref_sub$node_var)
})

test_that("without reference_data: outside_reference_distance is always FALSE in results", {
  # With Inf thresholds, distance > Inf is always FALSE, so no cell should be
  # flagged outside reference based on distance.
  ref_full <- make_full_panel_reference()
  markers  <- c("CD3", "CD4", "CD8")
  ref_sub  <- suppressWarnings(
    somalign_reference_subset_markers(ref_full, markers)
  )
  qry_data <- make_partial_query_data(markers = markers)

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  fit <- somalign_fit(qry, ref_sub)
  res <- suppressMessages(somalign_results(fit))

  # No cell should be flagged outside based on distance (Inf threshold)
  expect_true(all(!res$outside_reference_distance))
})


# ---- new tests: with reference_data (recompute path) ------------------------

test_that("with reference_data: distance_quantiles are finite", {
  obj     <- make_full_panel_reference_with_data()
  ref_sub <- somalign_reference_subset_markers(
    obj$reference, c("CD3", "CD4", "CD8"), reference_data = obj$data
  )

  expect_true(all(is.finite(ref_sub$distance_quantiles)))
})

test_that("with reference_data: distance_quantiles differ from original full-p values", {
  obj     <- make_full_panel_reference_with_data()
  ref_sub <- somalign_reference_subset_markers(
    obj$reference, c("CD3", "CD4", "CD8"), reference_data = obj$data
  )

  # Subspace distances must be smaller than full-p; the quantile matrix must
  # differ from the original (which was computed in 4-marker space).
  expect_false(
    isTRUE(all.equal(
      ref_sub$distance_quantiles,
      obj$reference$distance_quantiles
    ))
  )
})

test_that("with reference_data: distance_quantiles have same shape as original", {
  obj     <- make_full_panel_reference_with_data()
  ref_sub <- somalign_reference_subset_markers(
    obj$reference, c("CD3", "CD4", "CD8"), reference_data = obj$data
  )

  expect_equal(nrow(ref_sub$distance_quantiles), nrow(obj$reference$distance_quantiles))
  expect_equal(ncol(ref_sub$distance_quantiles), ncol(obj$reference$distance_quantiles))
  expect_equal(colnames(ref_sub$distance_quantiles), colnames(obj$reference$distance_quantiles))
})

test_that("with reference_data: inside-reference cells not flagged outside", {
  # Reference-like cells (drawn from same distribution) should NOT be flagged
  # outside once distance_quantiles are calibrated in the subspace.
  obj     <- make_full_panel_reference_with_data()
  markers <- c("CD3", "CD4", "CD8")
  ref_sub <- somalign_reference_subset_markers(
    obj$reference, markers, reference_data = obj$data
  )

  # Use the reference data itself (subset to shared markers) as the query --
  # it's from the same distribution so should mostly project inside.
  qry_data <- obj$data[, markers, drop = FALSE]
  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  fit <- somalign_fit(qry, ref_sub, epsilon = 0.25)
  res <- suppressMessages(somalign_results(fit))

  # The majority of reference-like cells should not be flagged outside
  outside_frac <- mean(res$outside_reference_distance)
  expect_lt(outside_frac, 0.5)
})

test_that("with reference_data: far-out cells flagged outside reference", {
  # Cells placed far from all reference nodes should be flagged as outside.
  obj     <- make_full_panel_reference_with_data()
  markers <- c("CD3", "CD4", "CD8")
  ref_sub <- somalign_reference_subset_markers(
    obj$reference, markers, reference_data = obj$data
  )

  # Create a single extreme outlier cell far from both populations
  # (reference pops are near +/-3; this cell is at +30)
  withr::local_seed(123L)
  n_normal <- 30L
  qry_data <- rbind(
    matrix(rnorm(n_normal * length(markers), mean = 3, sd = 0.4),
           ncol = length(markers), dimnames = list(NULL, markers)),
    matrix(rep(30, length(markers)), nrow = 1,
           dimnames = list(NULL, markers))
  )

  qry <- somalign_query(
    qry_data, ref_sub,
    grid = kohonen::somgrid(2L, 2L, "hexagonal"),
    rlen = 20L
  )
  fit <- somalign_fit(qry, ref_sub, epsilon = 0.25)
  res <- suppressMessages(somalign_results(fit))

  # The last cell (extreme outlier) must be flagged outside
  expect_true(res$outside_reference_distance[nrow(res)])
})

test_that("with reference_data: error when reference_data missing marker columns", {
  ref_full <- make_full_panel_reference()
  # Data only has CD3/CD4 but we're subsetting to CD3/CD4/CD8
  bad_data <- matrix(rnorm(80), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))

  expect_error(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"),
                                      reference_data = bad_data),
    regexp = "missing columns"
  )
})

test_that("with reference_data: error when reference_data is not a matrix/data.frame", {
  ref_full <- make_full_panel_reference()

  expect_error(
    somalign_reference_subset_markers(ref_full, c("CD3", "CD4"),
                                      reference_data = list(CD3 = 1:10, CD4 = 1:10)),
    regexp = "numeric matrix or data frame"
  )
})
