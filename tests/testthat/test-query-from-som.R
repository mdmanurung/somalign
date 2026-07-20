# Helper: build a tiny xyf SOM + reference + raw query data
make_query_som_fixtures <- function(seed = 42L,
                                    n_ref = 60L, n_query = 40L,
                                    p = 4L, n_labels = 2L,
                                    n_nodes_side = 2L) {
  withr::local_seed(seed)
  features <- paste0("F", seq_len(p))
  labels   <- rep(c("A", "B"), each = n_ref / 2L)
  grid     <- kohonen::somgrid(n_nodes_side, n_nodes_side, "hexagonal")

  X_ref <- matrix(rnorm(n_ref * p), nrow = n_ref, dimnames = list(NULL, features))
  Y_ref <- matrix(0, nrow = n_ref, ncol = n_labels,
                  dimnames = list(NULL, c("A", "B")))
  for (i in seq_along(labels)) Y_ref[i, labels[i]] <- 1

  som_ref <- kohonen::supersom(list(X_ref, Y_ref), grid = grid, rlen = 5L,
                               keep.data = TRUE)
  reference <- somalign_reference_from_som(
    som_ref,
    center = colMeans(X_ref),
    scale  = pmax(apply(X_ref, 2L, stats::sd), 1e-8),
    codebook_space = "raw"
  )

  X_query <- matrix(rnorm(n_query * p), nrow = n_query,
                    dimnames = list(NULL, features))
  center_q <- colMeans(X_query)
  scale_q  <- pmax(apply(X_query, 2L, stats::sd), 1e-8)
  X_query_scaled <- scale(X_query, center = center_q, scale = scale_q)
  som_query <- kohonen::som(X_query_scaled, grid = grid, rlen = 5L,
                            keep.data = TRUE)

  list(reference = reference, som_query = som_query, X_query = X_query,
       features = features, n_nodes = n_nodes_side^2L)
}

# ---------------------------------------------------------------------------
# Basic output structure
# ---------------------------------------------------------------------------

test_that("somalign_query_from_som returns a somalign_query with correct shape", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook_space = "reference_scaled")
  expect_s3_class(qry, "somalign_query")
  expect_equal(nrow(qry$data), nrow(fx$X_query))
  expect_equal(nrow(qry$scaled_data), nrow(fx$X_query))
  expect_equal(nrow(qry$codebook), fx$n_nodes)
  expect_equal(ncol(qry$codebook), length(fx$features))
  expect_equal(length(qry$sample_unit), nrow(fx$X_query))
  expect_equal(length(qry$sample_distance), nrow(fx$X_query))
  expect_true(all(is.na(qry$sample_distance)))
})

# ---------------------------------------------------------------------------
# sample_unit equals som$unit.classif
# ---------------------------------------------------------------------------

test_that("sample_unit matches som$unit.classif exactly", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook_space = "reference_scaled")
  expect_identical(qry$sample_unit, as.integer(fx$som_query$unit.classif))
})

# ---------------------------------------------------------------------------
# node_masses sum to 1 and match tabulate counts
# ---------------------------------------------------------------------------

test_that("node_masses sum to 1", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook_space = "reference_scaled")
  expect_equal(sum(qry$node_masses), 1, tolerance = 1e-10)
})

test_that("node_masses proportional to tabulate counts", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook_space = "reference_scaled")
  expected <- tabulate(fx$som_query$unit.classif, nbins = fx$n_nodes) /
    nrow(fx$X_query)
  expect_equal(qry$node_masses, expected, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# codebook argument (user-supplied pre-transformed codebook)
# ---------------------------------------------------------------------------

test_that("user-supplied codebook is used and column-reordered if needed", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  # Simulate a pre-transformed codebook (e.g. codes_ref_scaled from notebook)
  raw_cb <- kohonen::getCodes(fx$som_query, 1L)
  colnames(raw_cb) <- fx$features
  # Scale it with reference center/scale to get a "reference_scaled" codebook
  ref_cb_scaled <- scale(raw_cb,
                         center = fx$reference$center,
                         scale  = fx$reference$scale)
  # supply with shuffled columns to test reorder
  shuffled_cb <- ref_cb_scaled[, rev(fx$features), drop = FALSE]
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook = shuffled_cb,
                                 codebook_space = "reference_scaled")
  expect_equal(colnames(qry$codebook), fx$features)
})

test_that("codebook missing a feature raises an error", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  bad_cb <- matrix(0, nrow = fx$n_nodes, ncol = 1,
                   dimnames = list(NULL, "NOT_A_FEATURE"))
  expect_error(
    somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                            codebook = bad_cb),
    "missing feature"
  )
})

# ---------------------------------------------------------------------------
# codebook_space = "raw" applies center/scale
# ---------------------------------------------------------------------------

test_that("codebook_space='raw' scales the codebook by reference center/scale", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  raw_cb <- kohonen::getCodes(fx$som_query, 1L)
  colnames(raw_cb) <- fx$features

  qry_raw <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                     codebook = raw_cb,
                                     codebook_space = "raw")
  expected_cb <- scale(raw_cb,
                       center = fx$reference$center,
                       scale  = fx$reference$scale)
  expect_equal(qry_raw$codebook, expected_cb, tolerance = 1e-10,
               ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# Length mismatch between unit.classif and nrow(data)
# ---------------------------------------------------------------------------

test_that("unit.classif / data length mismatch raises an error", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  wrong_data <- fx$X_query[seq_len(nrow(fx$X_query) - 1L), , drop = FALSE]
  expect_error(
    somalign_query_from_som(fx$som_query, wrong_data, fx$reference),
    "must match"
  )
})

# ---------------------------------------------------------------------------
# Missing unit.classif raises an error
# ---------------------------------------------------------------------------

test_that("absent unit.classif raises a clear error", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  bad_som <- fx$som_query
  bad_som$unit.classif <- NULL
  expect_error(
    somalign_query_from_som(bad_som, fx$X_query, fx$reference),
    "unit.classif"
  )
})

# ---------------------------------------------------------------------------
# End-to-end: feeds somalign_fit without error
# ---------------------------------------------------------------------------

test_that("somalign_query_from_som feeds somalign_fit successfully", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()
  qry <- somalign_query_from_som(fx$som_query, fx$X_query, fx$reference,
                                 codebook_space = "reference_scaled")
  fit <- suppressMessages(suppressWarnings(
    somalign_fit(qry, fx$reference, epsilon = 0.1)
  ))
  expect_s3_class(fit, "somalign_fit")
  expect_true(is.finite(fit$diagnostics$ot$transport_mass))
})

# ---------------------------------------------------------------------------
# Reference output matches somalign_query() given same codebook
# ---------------------------------------------------------------------------

test_that("masses and codebook match somalign_query() when codebook is identical", {
  skip_if_not_installed("kohonen")
  fx <- make_query_som_fixtures()

  # somalign_query() projects cells onto som$codes[[1]] (in reference space)
  # to compute unit assignments.  somalign_query_from_som() uses unit.classif
  # directly.  When the SOM was trained in reference-scaled space they agree.
  withr::local_seed(1L)
  p <- 4L
  features <- paste0("F", seq_len(p))
  grid <- kohonen::somgrid(2L, 2L, "hexagonal")
  X <- matrix(rnorm(80L * p), nrow = 80L, dimnames = list(NULL, features))
  center <- colMeans(X)
  scale  <- pmax(apply(X, 2L, stats::sd), 1e-8)
  X_scaled <- scale(X, center = center, scale = scale)

  ref_mat <- matrix(rnorm(40L * p), nrow = 40L, dimnames = list(NULL, features))
  ref_mat_scaled <- scale(ref_mat, center = center, scale = scale)
  ref <- somalign_reference_from_nodes(
    codebook   = ref_mat_scaled[seq_len(4L), , drop = FALSE],
    features   = features,
    center     = center,
    scale      = scale,
    node_masses = rep(0.25, 4L)
  )

  # SOM trained IN reference-scaled space  →  unit.classif valid for that codebook
  som_q <- kohonen::som(X_scaled, grid = grid, rlen = 5L, keep.data = TRUE)

  qry_from <- somalign_query_from_som(som_q, X, ref,
                                      codebook_space = "reference_scaled")
  qry_std  <- somalign_query(X, ref, som_query = som_q,
                             codebook_space = "reference_scaled")

  expect_identical(qry_from$sample_unit,  qry_std$sample_unit)
  expect_equal(qry_from$node_masses, qry_std$node_masses, tolerance = 1e-10)
  expect_equal(qry_from$codebook,    qry_std$codebook,    tolerance = 1e-10)
  # Both constructors must expose label_prob (identical structure); for an
  # unlabelled SOM both are empty matrices.
  expect_true("label_prob" %in% names(qry_from))
  expect_equal(qry_from$label_prob, qry_std$label_prob)
})
