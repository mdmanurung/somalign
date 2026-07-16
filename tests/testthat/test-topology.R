## ---------------------------------------------------------------------------
## Tests for somalign_topology_audit() (Idea #6: persistent-homology topology
## audit). H0 via single-linkage/union-find; the F2 fix means
## .somalign_pairwise_distance() returns SQUARED distance, so sqrt() is
## required before persistent homology (which needs a genuine metric).
## ---------------------------------------------------------------------------

test_that(".somalign_h0_persistence detects 2 well-separated clusters", {
  cb <- matrix(c(0, 0, 0.1, 0, 0, 0.1,
                 10, 10, 10.1, 10, 10, 10.1), ncol = 2, byrow = TRUE)
  d <- sqrt(somalign:::.somalign_pairwise_distance(cb, cb))
  pd <- somalign:::.somalign_h0_persistence(d)
  n <- somalign:::.somalign_h0_n_components(pd, threshold = 1, n_nodes = 6L)
  expect_equal(n, 2L)
})

test_that(".somalign_h0_persistence handles a single node without error", {
  cb <- matrix(1:2, nrow = 1)
  d <- sqrt(somalign:::.somalign_pairwise_distance(cb, cb))
  pd <- somalign:::.somalign_h0_persistence(d)
  expect_equal(nrow(pd), 0L)
  n <- somalign:::.somalign_h0_n_components(pd, threshold = 0.1, n_nodes = 1L)
  expect_equal(n, 1L)
})

test_that("topology_warning fires when correction merges two clusters", {
  ref_cb <- matrix(c(0, 0, 0, 0.1, 0.1, 0,
                     5, 5, 5, 5.1, 5.1, 5), ncol = 2, byrow = TRUE)
  colnames(ref_cb) <- c("F1", "F2")
  qry_cb <- ref_cb + 0.05
  # node_shifts collapse both clusters onto their joint centroid.
  shifts <- matrix(2.5, nrow = 6, ncol = 2) - qry_cb
  colnames(shifts) <- c("F1", "F2")
  attr(shifts, "correction_allowed") <- rep(TRUE, 6L)

  mock_fit <- structure(
    list(
      query = list(codebook = qry_cb),
      reference = list(codebook = ref_cb, distance_quantiles = matrix(0.3, 1, 1, dimnames = list(NULL, "95%"))),
      node_shifts = shifts
    ),
    class = "somalign_fit"
  )
  expect_warning(ta <- somalign_topology_audit(mock_fit), "topology_warning")
  expect_true(ta$topology_warning)
  expect_lt(ta$n_components_corrected, ta$n_components_query)
  expect_equal(ta$topology_delta, ta$n_components_corrected - ta$n_components_query)
})

test_that("no topology_warning when clusters remain separated after correction", {
  ref_cb <- matrix(c(0, 0, 0, 0.1, 0.1, 0,
                     5, 5, 5, 5.1, 5.1, 5), ncol = 2, byrow = TRUE)
  colnames(ref_cb) <- c("F1", "F2")
  qry_cb <- ref_cb + 0.05
  shifts <- matrix(0, nrow = 6, ncol = 2)  # no correction at all
  colnames(shifts) <- c("F1", "F2")
  attr(shifts, "correction_allowed") <- rep(TRUE, 6L)
  mock_fit <- structure(
    list(
      query = list(codebook = qry_cb),
      reference = list(codebook = ref_cb, distance_quantiles = matrix(0.3, 1, 1, dimnames = list(NULL, "95%"))),
      node_shifts = shifts
    ),
    class = "somalign_fit"
  )
  ta <- somalign_topology_audit(mock_fit)
  expect_false(ta$topology_warning)
  expect_equal(ta$topology_delta, 0L)
})

test_that("topology audit works without TDA and degrades gracefully", {
  skip_if(requireNamespace("TDA", quietly = TRUE),
          "TDA is installed; this test checks the TDA-absent path only")
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  expect_message(ta <- somalign_topology_audit(fit, use_tda = TRUE), "TDA package not available")
  expect_null(ta$bottleneck_h0)
  expect_null(ta$tda_query)
  expect_type(ta$n_components_corrected, "integer")
})

test_that("somalign_diagnostics(topology = TRUE) is additive", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  d0 <- somalign_diagnostics(fit)
  d1 <- somalign_diagnostics(fit, topology = TRUE)
  expect_null(d0$topology)
  expect_s3_class(d1$topology, "somalign_topology")
  expect_equal(d0$solver, d1$solver)
  expect_equal(d0$ot, d1$ot)
})

test_that("print.somalign_topology does not error", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  ta <- somalign_topology_audit(fit)
  expect_output(print(ta))
})

test_that("somalign_topology_audit errors on a non-fit object", {
  expect_error(somalign_topology_audit(list()), "somalign_fit")
})
