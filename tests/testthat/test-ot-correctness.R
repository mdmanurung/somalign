## Numerical-correctness tests for the OT solver and diagnostic outputs.
##
## These tests verify:
##   (a) solver determinism
##   (b) near-identity transport when query codebook == reference codebook
##   (c) presence of new diagnostic fields (converged, final_delta, cost_scale)
##   (d) known-truth label-transfer accuracy on well-separated synthetic data
##   (e) all-zero mass warning from the input validator

test_that("somalign_fit is deterministic: same inputs give identical transport plans", {
  skip_if_not_installed("kohonen")
  set.seed(11)
  mat <- rbind(
    matrix(rnorm(15 * 4, mean = -1), ncol = 4),
    matrix(rnorm(15 * 4, mean =  1), ncol = 4)
  )
  colnames(mat) <- paste0("f", seq_len(ncol(mat)))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  qry <- somalign_query(
    mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  fit1 <- somalign_fit(qry, ref)
  fit2 <- somalign_fit(qry, ref)
  expect_identical(fit1$transport_plan, fit2$transport_plan)
})

test_that("identical codebooks concentrate plan on diagonal", {
  skip_if_not_installed("kohonen")
  set.seed(22)
  mat <- rbind(
    matrix(rnorm(20 * 3, mean = -2), ncol = 3),
    matrix(rnorm(20 * 3, mean =  2), ncol = 3)
  )
  colnames(mat) <- paste0("f", seq_len(ncol(mat)))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  # Supply the reference SOM as the query SOM so codebooks are identical.
  qry <- somalign_query(mat, ref, som_query = ref$som_ref,
                        codebook_space = "reference_scaled")
  fit <- somalign_fit(qry, ref, epsilon = 0.1)
  plan <- fit$transport_plan
  n    <- min(nrow(plan), ncol(plan))
  diag_mass  <- sum(diag(plan[seq_len(n), seq_len(n)]))
  total_mass <- sum(plan)
  # diagonal should carry at least half the total transported mass
  expect_gt(diag_mass / max(total_mass, .Machine$double.xmin), 0.5)
})

test_that("diagnostics$solver contains converged, final_delta, and cost_scale fields", {
  skip_if_not_installed("kohonen")
  set.seed(33)
  mat <- matrix(rnorm(20 * 2), ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  qry <- somalign_query(
    mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  fit <- somalign_fit(qry, ref)
  s <- fit$diagnostics$solver
  expect_true("converged"    %in% names(s))
  expect_true("final_delta"  %in% names(s))
  expect_true("cost_scale"   %in% names(s))
  expect_true(is.logical(s$converged))
  expect_length(s$converged, 1L)
  expect_true(is.finite(s$final_delta) || is.nan(s$final_delta))
  expect_gt(s$cost_scale, 0)
})

test_that("known-truth label transfer achieves >80% accuracy on well-separated clusters", {
  skip_if_not_installed("kohonen")
  set.seed(44)
  n <- 50L
  # Two widely separated clusters; signal >> noise → correct node assignment.
  ref_data <- rbind(
    matrix(rnorm(n * 4, mean = -5, sd = 0.5), ncol = 4),
    matrix(rnorm(n * 4, mean =  5, sd = 0.5), ncol = 4)
  )
  colnames(ref_data) <- paste0("f", seq_len(ncol(ref_data)))
  labels <- rep(c("neg", "pos"), each = n)
  ref <- somalign_train_reference(
    ref_data, labels = labels,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 20
  )
  # Query: same cluster structure, tiny shift.
  qry_data <- rbind(
    matrix(rnorm(n * 4, mean = -5.1, sd = 0.5), ncol = 4),
    matrix(rnorm(n * 4, mean =  5.1, sd = 0.5), ncol = 4)
  )
  colnames(qry_data) <- colnames(ref_data)
  qry <- somalign_query(
    qry_data, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 20
  )
  # Use tight marginal constraints (high rho) so the coupling stays near-balanced.
  fit <- somalign_fit(qry, ref, rho_query = 10, rho_ref = 10)
  results <- somalign_results(fit)

  true_labels   <- rep(c("neg", "pos"), each = n)
  transferred   <- results$label
  accepted_idx  <- !is.na(transferred)
  # At minimum some labels should have been accepted
  expect_gt(sum(accepted_idx), 0L)
  accuracy <- mean(transferred[accepted_idx] == true_labels[accepted_idx])
  expect_gt(accuracy, 0.8)
})

test_that("all-zero query mass triggers a warning from .somalign_validate_ot_inputs", {
  cost <- matrix(c(0, 1, 1, 0, 0.5, 0.5), nrow = 3, ncol = 2)
  a <- rep(0, 3)           # all-zero query masses
  b <- c(0.5, 0.5)
  expect_warning(
    somalign:::.somalign_validate_ot_inputs(
      cost, a, b, epsilon = 0.5, rho_query = 1, rho_ref = 1
    ),
    "zero"
  )
})

test_that("codebook_space = 'raw' in somalign_query rescales codebook to reference space", {
  skip_if_not_installed("kohonen")
  set.seed(55)
  mat <- rbind(
    matrix(rnorm(12 * 3, mean = -1), ncol = 3),
    matrix(rnorm(12 * 3, mean =  1), ncol = 3)
  )
  colnames(mat) <- paste0("f", seq_len(ncol(mat)))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  # Build a SOM on raw (unscaled) data, then let somalign_query rescale it.
  raw_som <- kohonen::som(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry_raw  <- somalign_query(mat, ref, som_query = raw_som,
                             codebook_space = "raw")
  qry_ref  <- somalign_query(mat, ref,
                             grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  # Both codebooks should now be in the reference-scaled space (finite, same dims).
  expect_equal(ncol(qry_raw$codebook), length(ref$features))
  expect_true(all(is.finite(qry_raw$codebook)))
  # The rescaled codebook should differ from the raw one.
  expect_false(identical(qry_raw$codebook, raw_som$codes[[1]]))
})
