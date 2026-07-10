test_that("training, alignment, and sensitivity grid work on synthetic data", {
  skip_if_not_installed("kohonen")

  set.seed(42)
  old <- rbind(
    matrix(rnorm(20 * 40, mean = -1), ncol = 40),
    matrix(rnorm(20 * 40, mean = 1), ncol = 40)
  )
  colnames(old) <- paste0("f", seq_len(ncol(old)))
  labels <- rep(c("left", "right"), each = 20)

  reference <- somalign_train_reference(
    old,
    labels = labels,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )

  query <- rbind(
    old[1:10, ] + 0.1,
    old[31:40, ] + 0.3,
    matrix(rnorm(5 * 40, mean = 5), ncol = 40)
  )
  colnames(query) <- colnames(old)

  query_obj <- somalign_query(
    query,
    reference,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )
  fit <- somalign_fit(query_obj, reference, solver = "internal", epsilon = 0.1)
  results <- somalign_results(fit)

  expect_s3_class(reference, "somalign_reference")
  expect_s3_class(query_obj, "somalign_query")
  expect_equal(nrow(results), nrow(query))
  expect_true(any(results$outside_reference_distance))

  grid <- somalign_sensitivity_grid(
    query_obj,
    reference,
    epsilon = c(0.05, 0.1),
    rho_query = c(0.5, 1),
    rho_ref = 1,
    solver = "internal"
  )
  expect_equal(nrow(grid), 4L)
  expect_true(all(c("epsilon", "rho_query", "rho_ref", "transport_mass", "mean_match_fraction") %in% names(grid)))

  expect_error(
    somalign_sensitivity_grid(
      query_obj,
      reference,
      epsilon = numeric(),
      rho_query = 1,
      rho_ref = 1,
      solver = "internal"
    ),
    "epsilon"
  )
})

test_that("a separately trained query SOM can be corrected and projected to old nodes", {
  skip_if_not_installed("kohonen")

  set.seed(7)
  old <- rbind(
    matrix(rnorm(12 * 8, mean = -1), ncol = 8),
    matrix(rnorm(12 * 8, mean = 1), ncol = 8)
  )
  colnames(old) <- paste0("f", seq_len(ncol(old)))
  reference <- somalign_train_reference(
    old,
    labels = rep(c("low", "high"), each = 12),
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )

  query_data <- rbind(old[1:6, ] + 0.1, old[19:24, ] + 0.25)
  colnames(query_data) <- colnames(old)
  query_scaled_for_som <- sweep(
    sweep(query_data[, reference$features], 2, reference$center, "-"),
    2,
    reference$scale,
    "/"
  )
  query_som <- kohonen::som(
    query_scaled_for_som,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 5
  )

  query <- somalign_query(query_data, reference, som_query = query_som)
  fit <- somalign_fit(query, reference, solver = "internal")
  results <- somalign_results(fit)

  expect_equal(nrow(results), nrow(query_data))
  expect_true(all(is.finite(results$correction_norm)))
  expect_true(all(results$corrected_som_unit %in% seq_len(nrow(reference$codebook))))
  expect_true(all(results$old_som_unit %in% seq_len(nrow(reference$codebook))))
})

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

test_that("POT comparison is optional when reticulate can import ot.unbalanced", {
  skip_if_not(identical(Sys.getenv("SOMALIGN_RUN_POT_TESTS"), "true"))
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("ot.unbalanced"))

  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  pot <- somalign_fit(query_obj, ref, solver = "pot", epsilon = 0.1)
  internal <- somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1)

  expect_equal(dim(pot$transport_plan), dim(internal$transport_plan))
  expect_true(all(is.finite(pot$transport_plan)))
})
