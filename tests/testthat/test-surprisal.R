## ---------------------------------------------------------------------------
## Tests for the chi-squared surprisal outside_reference diagnostic (Idea #4).
## reference$node_var is computed at reference-construction time; surprisal
## is a calibrated alternative to the distance-quantile outside_reference
## flag, weighting per-marker deviations by within-node variance.
## ---------------------------------------------------------------------------

test_that("node_var has correct dimensions and is positive", {
  set.seed(3)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  nv <- ref$node_var
  expect_equal(nrow(nv), nrow(ref$codebook))
  expect_equal(ncol(nv), length(ref$features))
  expect_equal(colnames(nv), ref$features)
  expect_true(all(nv > 0))
})

test_that("cell at node centroid scores near-zero surprisal and high p-value", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_false(is.null(ref$node_var))

  centroid <- matrix(ref$codebook[1, ], nrow = 1, dimnames = list(NULL, c("F1", "F2")))
  surpr <- somalign:::.somalign_node_surprisal_core(centroid, 1L, ref$codebook, ref$node_var)
  expect_equal(surpr$surprisal, 0, tolerance = 1e-10)
  expect_gt(surpr$pvalue, 0.99)
})

test_that("single-marker offset identifies that marker as top_marker", {
  set.seed(7)
  mat <- matrix(rnorm(450), nrow = 150, ncol = 3,
                dimnames = list(NULL, c("CD3", "CD11c", "CD19")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)

  centroid <- ref$codebook[1, , drop = FALSE]
  shifted <- centroid
  shifted[, "CD11c"] <- centroid[, "CD11c"] + 10
  surpr <- somalign:::.somalign_node_surprisal_core(shifted, 1L, ref$codebook, ref$node_var)
  expect_equal(surpr$top_marker, "CD11c")
  expect_lt(surpr$pvalue, 0.01)
})

test_that("old reference without node_var gives NA surprisal columns and messages", {
  set.seed(1)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5,
                                  compute_node_var = FALSE)
  expect_null(ref$node_var)

  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  expect_message(res <- somalign_results(fit), "node_var.*absent")
  expect_true(all(is.na(res$outside_reference_surprisal)))
  expect_true(all(is.na(res$outside_reference_pvalue)))
  expect_true(all(is.na(res$outside_reference_top_marker)))
})

test_that("somalign_results includes surprisal columns with valid ranges", {
  set.seed(5)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  res <- somalign_results(fit)
  expect_true(all(c("outside_reference_surprisal", "outside_reference_pvalue",
                    "outside_reference_top_marker") %in% names(res)))
  expect_true(is.numeric(res$outside_reference_surprisal))
  expect_true(all(res$outside_reference_pvalue >= 0 & res$outside_reference_pvalue <= 1))
})

test_that("outside_pvalue_threshold adds a matching boolean flag column", {
  set.seed(9)
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit(qry, ref)
  res <- somalign_results(fit, outside_pvalue_threshold = 0.05)
  expect_true("outside_reference_pvalue_flag" %in% names(res))
  expect_type(res$outside_reference_pvalue_flag, "logical")
  expect_equal(res$outside_reference_pvalue_flag, res$outside_reference_pvalue < 0.05)

  res_default <- somalign_results(fit)
  expect_false("outside_reference_pvalue_flag" %in% names(res_default))
})

test_that("somalign_reference_from_som computes a conforming node_var", {
  skip_if_not_installed("kohonen")
  set.seed(11)
  n <- 80
  X <- matrix(rnorm(n * 2), nrow = n, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  center <- colMeans(X)
  scale <- apply(X, 2, stats::sd)
  X_scaled <- scale(X, center = center, scale = scale)
  som_scaled <- kohonen::som(X_scaled, grid = kohonen::somgrid(2, 2, "hexagonal"),
                             rlen = 5, keep.data = TRUE)
  ref <- somalign_reference_from_som(som_scaled, center = center, scale = scale,
                                     codebook_space = "reference_scaled")
  expect_false(is.null(ref$node_var))
  expect_equal(dim(ref$node_var), c(nrow(ref$codebook), length(ref$features)))
  expect_true(all(ref$node_var > 0))

  ref_off <- somalign_reference_from_som(som_scaled, center = center, scale = scale,
                                         codebook_space = "reference_scaled",
                                         compute_node_var = FALSE)
  expect_null(ref_off$node_var)
})

test_that("somalign_reference_from_nodes accepts and validates an explicit node_var", {
  cb <- matrix(c(0.1, 0.2, -0.1, 0.3, 0.4, -0.2, 0.0, 0.1),
              nrow = 4, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  nv <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 4, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_reference_from_nodes(
    codebook = cb, features = c("F1", "F2"),
    center = c(F1 = 0, F2 = 0), scale = c(F1 = 1, F2 = 1),
    node_var = nv
  )
  expect_equal(ref$node_var, nv)

  expect_error(
    somalign_reference_from_nodes(
      codebook = cb, features = c("F1", "F2"),
      center = c(F1 = 0, F2 = 0), scale = c(F1 = 1, F2 = 1),
      node_var = nv[1:2, , drop = FALSE]
    ),
    "one row per reference node"
  )
})
