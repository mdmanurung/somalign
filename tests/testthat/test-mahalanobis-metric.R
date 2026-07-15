## ---------------------------------------------------------------------------
## Tests for feature_weights (Idea #9: learned diagonal Mahalanobis OT cost).
## Weights whiten both codebooks before the squared-Euclidean cost is built;
## projection/threshold distances (.somalign_nearest_code) are untouched.
## ---------------------------------------------------------------------------

test_that("feature_weights = NULL produces an identical fit to the current default", {
  set.seed(42)
  mat <- rbind(matrix(rnorm(30 * 3, -2), ncol = 3),
               matrix(rnorm(30 * 3, 2), ncol = 3))
  colnames(mat) <- c("F1", "F2", "F3")
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.3, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit_default <- somalign_fit(qry, ref)
  fit_null <- somalign_fit(qry, ref, feature_weights = NULL)
  expect_equal(fit_default$cost, fit_null$cost)
  expect_equal(fit_default$transport_plan, fit_null$transport_plan)
  expect_equal(fit_default$node_shifts, fit_null$node_shifts)
  expect_null(fit_default$diagnostics$cost_metric$feature_weights)
})

test_that("zeroing a marker's weight removes it from the cost and shifts the plan", {
  set.seed(1)
  p <- 3L
  mat <- matrix(rnorm(40 * p), ncol = p, dimnames = list(NULL, paste0("F", 1:p)))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)

  fw_equal <- c(F1 = 1, F2 = 1, F3 = 1)
  fit_equal <- somalign_fit(qry, ref, feature_weights = fw_equal)
  fit_default <- somalign_fit(qry, ref)
  expect_equal(fit_equal$cost, fit_default$cost)

  fw_no_f1 <- c(F1 = 0, F2 = 1, F3 = 1)
  fit_no_f1 <- somalign_fit(qry, ref, feature_weights = fw_no_f1)
  expected_cost <- somalign:::.somalign_pairwise_distance(
    sweep(qry$codebook, 2, sqrt(fw_no_f1), "*"),
    sweep(ref$codebook, 2, sqrt(fw_no_f1), "*")
  )
  expect_equal(fit_no_f1$cost, expected_cost)
  expect_false(isTRUE(all.equal(fit_no_f1$transport_plan, fit_default$transport_plan)))
  expect_equal(fit_no_f1$diagnostics$cost_metric$feature_weights, fw_no_f1)
})

test_that("feature_weights = \"anchor\" down-weights the high-variance displacement marker", {
  skip_if_not_installed("kohonen")
  set.seed(99)
  p <- 4L
  nm <- paste0("F", 1:p)
  mat <- matrix(rnorm(60 * p), ncol = p, dimnames = list(NULL, nm))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  anc_idx <- 1:15
  anchor_old <- mat[anc_idx, , drop = FALSE]
  # F1 has large noisy per-anchor displacement (batch-driven); F2-F4 are stable.
  # sweep() applies the per-column multiplier AFTER the matrix is shaped;
  # multiplying the flat rnorm() vector first would recycle the 4-element
  # multiplier across the flat 60-element vector, scrambling it across columns.
  noise <- matrix(rnorm(15 * p), nrow = 15, ncol = p)
  anchor_new <- anchor_old + sweep(noise, 2, c(2, 0.01, 0.01, 0.01), "*")
  colnames(anchor_old) <- colnames(anchor_new) <- nm
  qry <- somalign_query(mat + matrix(c(5, 0, 0, 0), nrow(mat), p, byrow = TRUE), ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit <- somalign_fit_anchored(qry, ref, anchor_old = anchor_old, anchor_new = anchor_new,
                               feature_weights = "anchor")
  fw <- fit$anchors$feature_weights
  expect_true(fw[["F1"]] < fw[["F2"]])
  expect_true(fw[["F1"]] < fw[["F3"]])
  expect_equal(fit$diagnostics$cost_metric$feature_weights, fw)
})

test_that("feature_weights = \"anchor\" errors clearly on plain somalign_fit", {
  set.seed(1)
  mat <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_error(
    somalign_fit(qry, ref, feature_weights = "anchor"),
    "somalign_fit_anchored"
  )
})

test_that("projection distances are identical regardless of feature_weights", {
  set.seed(7)
  mat <- matrix(rnorm(30 * 3), ncol = 3, dimnames = list(NULL, c("A", "B", "C")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat + 0.2, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fw <- c(A = 0.1, B = 5, C = 1)
  fit_w <- somalign_fit(qry, ref, feature_weights = fw)
  fit_d <- somalign_fit(qry, ref)
  res_w <- somalign_results(fit_w)
  res_d <- somalign_results(fit_d)
  expect_equal(res_w$outside_reference_distance, res_d$outside_reference_distance)
  expect_equal(res_w$old_som_unit, res_d$old_som_unit)
})

test_that("feature_weights validation catches bad inputs", {
  set.seed(3)
  mat <- matrix(rnorm(20 * 2), ncol = 2, dimnames = list(NULL, c("X", "Y")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)

  expect_error(somalign_fit(qry, ref, feature_weights = c(1, 2, 3)), "one entry per feature")
  expect_error(somalign_fit(qry, ref, feature_weights = c(X = -1, Y = 1)), "non-negative")
  expect_error(somalign_fit(qry, ref, feature_weights = c(X = 0, Y = 0)), "all zeros")
  expect_error(somalign_fit(qry, ref, feature_weights = c(Z = 1, Y = 1)), "missing names")

  fw <- c(X = 2, Y = 0.5)
  fit <- somalign_fit(qry, ref, feature_weights = fw)
  expect_equal(fit$diagnostics$cost_metric$feature_weights, fw[c("X", "Y")])
})
