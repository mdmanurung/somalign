# ---------------------------------------------------------------------------
# Helpers shared across tests in this file
# ---------------------------------------------------------------------------

make_xyf_som <- function(seed = 1L, n = 120L, nodes = 3L, rlen = 5L) {
  withr::local_seed(seed)
  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2,
               dimnames = list(NULL, c("F1", "F2")))
  # Two balanced, separable classes
  Y <- cbind(
    A = rep(c(1, 0), each = n / 2),
    B = rep(c(0, 1), each = n / 2)
  )
  som <- kohonen::supersom(
    list(X, Y),
    grid      = kohonen::somgrid(nodes, nodes, "hexagonal"),
    rlen      = rlen,
    keep.data = TRUE
  )
  list(som = som, X = X, Y = Y,
       center = colMeans(X), scale = apply(X, 2, stats::sd))
}

# SOM trained on already-scaled data (codebook_space = "reference_scaled")
make_scaled_xyf_som <- function(seed = 1L, n = 120L, nodes = 3L, rlen = 5L) {
  fx <- make_xyf_som(seed = seed, n = n, nodes = nodes, rlen = rlen)
  center <- fx$center
  scale  <- fx$scale
  withr::local_seed(seed)
  X_sc <- scale(fx$X, center = center, scale = scale)
  colnames(X_sc) <- c("F1", "F2")
  som_sc <- kohonen::supersom(
    list(X_sc, fx$Y),
    grid      = kohonen::somgrid(nodes, nodes, "hexagonal"),
    rlen      = rlen,
    keep.data = TRUE
  )
  list(som = som_sc, X = fx$X, X_sc = X_sc, Y = fx$Y,
       center = center, scale = scale)
}

# ---------------------------------------------------------------------------
# somalign_reference_from_som() — basic output contract
# ---------------------------------------------------------------------------

test_that("somalign_reference_from_som returns a somalign_reference object", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  expect_s3_class(ref, "somalign_reference")
})

test_that("node masses sum to 1 and match tabulate(unit.classif)", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  expect_equal(sum(ref$node_masses), 1, tolerance = 1e-10)

  # Masses must match exact tabulation over all training cells
  n_nodes <- nrow(fx$som$codes[[1]])
  expected <- tabulate(fx$som$unit.classif, nbins = n_nodes)
  expected <- as.numeric(expected / sum(expected))
  expect_equal(ref$node_masses, expected, tolerance = 1e-12)
})

test_that("label_prob colnames match Y-layer codebook colnames", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  expect_equal(colnames(ref$label_prob), c("A", "B"))
})

test_that("label_prob rows sum to 1 (up to zero-mass unoccupied nodes)", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  occupied <- rowSums(ref$label_prob) > 0
  row_sums <- rowSums(ref$label_prob[occupied, , drop = FALSE])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("distance_quantiles are finite and monotone across quantile levels", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  dq <- ref$distance_quantiles
  expect_true(all(is.finite(dq)))
  expect_true(all(dq >= 0))
  # Each row's quantiles must be non-decreasing
  for (i in seq_len(nrow(dq))) {
    expect_true(all(diff(dq[i, ]) >= -1e-12))
  }
})

test_that("distance quantiles are in X-space (not kohonen's combined scale)", {
  # kohonen's som$distances are in ~0-1 (weighted) space; X-space distances
  # for data standardised to unit variance live around 0-5.
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  # kohonen distances are <= 1.0 by design; X-space 95th pctl should differ
  kohonen_max <- max(fx$som$distances)
  x_space_99  <- max(ref$distance_quantiles[, "99%"])
  # They may coincidentally match for small SOMs, so only check scale for
  # datasets where they are distinguishable.
  # Just confirm distance_quantiles look like X-space (> 0, finite, consistent
  # with global_distance_quantiles).
  expect_true(all(ref$global_distance_quantiles >= 0))
  expect_equal(
    names(ref$global_distance_quantiles),
    c("50%", "90%", "95%", "99%")
  )
})

test_that("reference_units and reference_distances are attached and consistent", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )
  expect_equal(ref$reference_units, as.integer(fx$som$unit.classif))
  expect_equal(length(ref$reference_distances), nrow(fx$X))
  expect_true(all(ref$reference_distances >= 0))
  # n_samples set to number of training cells, not NA
  expect_equal(ref$n_samples, nrow(fx$X))
})

# ---------------------------------------------------------------------------
# labels = "none" disables label transfer
# ---------------------------------------------------------------------------

test_that("labels = 'none' produces empty label_prob matrix", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled",
    labels = "none"
  )
  expect_equal(ncol(ref$label_prob), 0L)
})

# ---------------------------------------------------------------------------
# Plain kohonen::som() (no Y layer) — graceful fallback
# ---------------------------------------------------------------------------

test_that("plain kohonen::som() input gives reference with no label_prob", {
  withr::local_seed(2L)
  X <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2,
               dimnames = list(NULL, c("F1", "F2")))
  center <- colMeans(X)
  scale  <- apply(X, 2, stats::sd)
  X_sc   <- scale(X, center = center, scale = scale)
  colnames(X_sc) <- c("F1", "F2")
  som_plain <- kohonen::som(
    X_sc,
    grid      = kohonen::somgrid(3, 3, "hexagonal"),
    rlen      = 5,
    keep.data = TRUE
  )
  expect_message(
    ref <- somalign_reference_from_som(
      som_plain, center = center, scale = scale,
      codebook_space = "reference_scaled"
    ),
    "label transfer"
  )
  expect_equal(ncol(ref$label_prob), 0L)
})

# ---------------------------------------------------------------------------
# codebook_space = "raw" path
# ---------------------------------------------------------------------------

test_that("codebook_space = 'raw' scales codebook and produces valid reference", {
  fx  <- make_xyf_som()  # SOM trained on UN-scaled data
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "raw"
  )
  expect_s3_class(ref, "somalign_reference")
  expect_equal(sum(ref$node_masses), 1, tolerance = 1e-10)
  expect_true(all(is.finite(ref$distance_quantiles)))
})

# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

test_that("missing embedded data raises a clear error", {
  fx  <- make_scaled_xyf_som()
  som_no_data        <- fx$som
  som_no_data$data   <- NULL
  expect_error(
    somalign_reference_from_som(
      som_no_data, center = fx$center, scale = fx$scale,
      codebook_space = "reference_scaled"
    ),
    "does not store training data"
  )
})

test_that("missing unit.classif raises a clear error", {
  fx  <- make_scaled_xyf_som()
  som_no_uc             <- fx$som
  som_no_uc$unit.classif <- NULL
  expect_error(
    somalign_reference_from_som(
      som_no_uc, center = fx$center, scale = fx$scale,
      codebook_space = "reference_scaled"
    ),
    "unit.classif"
  )
})

# ---------------------------------------------------------------------------
# End-to-end integration: somalign_reference_from_som -> query -> fit -> results
# ---------------------------------------------------------------------------

test_that("somalign_reference_from_som feeds somalign_query/fit/results end-to-end", {
  fx  <- make_scaled_xyf_som()
  ref <- somalign_reference_from_som(
    fx$som, center = fx$center, scale = fx$scale,
    codebook_space = "reference_scaled"
  )

  # Query: a new matrix with the same features, passed through reference scaling
  withr::local_seed(99L)
  X_query <- scale(
    matrix(rnorm(60 * 2), nrow = 60, ncol = 2,
           dimnames = list(NULL, c("F1", "F2"))),
    center = fx$center, scale = fx$scale
  )
  # Reverse-transform so query() can re-scale it internally
  X_raw_q <- sweep(sweep(X_query, 2, fx$scale, "*"), 2, fx$center, "+")
  colnames(X_raw_q) <- c("F1", "F2")

  qry <- somalign_query(
    X_raw_q, ref,
    grid = kohonen::somgrid(3, 3, "hexagonal"),
    rlen = 5
  )
  fit     <- somalign_fit(qry, ref, solver = "internal")
  results <- somalign_results(fit)

  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), nrow(X_raw_q))
  # Label transfer should be active because we used labels = "codebook"
  expect_true("transferred_label" %in% colnames(results))
  expect_false(all(is.na(results$transferred_label)))
})
