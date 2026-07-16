## Supervised plan tuning against label accuracy (item 3).

skip_if_not_installed("kohonen")

make_tune_data <- function(n = 250L, sep = 3, seed = 1L) {
  withr::local_seed(seed)
  x <- rbind(
    matrix(rnorm(n * 3, -sep, 0.5), ncol = 3),
    matrix(rnorm(n * 3,  sep, 0.5), ncol = 3)
  )
  colnames(x) <- paste0("f", seq_len(3))
  list(x = x, labels = rep(c("low", "high"), each = n))
}

test_that("somalign_tune sweeps a data-frame grid and picks a best combo", {
  d <- make_tune_data()
  tuned <- somalign_tune(
    d$x, d$labels, grid = kohonen::somgrid(3, 3, "hexagonal"),
    param_grid = data.frame(epsilon = c(0.05, 0.1, 0.2)),
    k = 3, rlen = 20, metric = "mcc"
  )
  expect_s3_class(tuned, "somalign_tune")
  expect_equal(nrow(tuned$grid), 3L)
  expect_true(all(c("epsilon", "accuracy", "macro_f1", "mcc", "coverage", "ece")
                  %in% names(tuned$grid)))
  # best MCC row is the argmax of the mcc column
  expect_equal(tuned$best$mcc, max(tuned$grid$mcc))
  expect_true(tuned$best_params$epsilon %in% c(0.05, 0.1, 0.2))
})

test_that("somalign_tune minimises ECE when metric = 'ece'", {
  d <- make_tune_data(n = 150L)
  tuned <- somalign_tune(
    d$x, d$labels, grid = kohonen::somgrid(2, 2, "hexagonal"),
    param_grid = data.frame(epsilon = c(0.05, 0.2)),
    k = 2, rlen = 15, metric = "ece"
  )
  expect_equal(tuned$best$ece, min(tuned$grid$ece))
})

test_that("somalign_tune accepts a list-of-lists grid with feature_weights", {
  d <- make_tune_data(n = 150L)
  pg <- list(
    list(epsilon = 0.1),
    list(epsilon = 0.1, feature_weights = c(2, 1, 1))
  )
  tuned <- somalign_tune(
    d$x, d$labels, grid = kohonen::somgrid(2, 2, "hexagonal"),
    param_grid = pg, k = 2, rlen = 15
  )
  expect_equal(nrow(tuned$grid), 2L)
  expect_setequal(tuned$grid$feature_weights, c("none", "custom"))
  op <- withVisible(print(tuned))
  expect_false(op$visible)
})

test_that("somalign_tune requires epsilon in every combination", {
  d <- make_tune_data(n = 40L)
  expect_error(
    somalign_tune(d$x, d$labels, grid = kohonen::somgrid(2, 2, "hexagonal"),
                  param_grid = data.frame(rho_query = c(1, 2))),
    "must specify `epsilon`"
  )
})

test_that(".somalign_ot_sweep_one applies feature_weights to the cost", {
  # Weighting a feature by ~0 should change the plan vs unweighted.
  fx <- make_anchored_fixture()
  base <- somalign:::.somalign_ot_sweep_one(
    fx$qry, fx$ref, epsilon = 0.1, rho_query = 1, rho_ref = 1,
    solver = "internal", max_iter = 500, tol = 1e-7
  )
  weighted <- somalign:::.somalign_ot_sweep_one(
    fx$qry, fx$ref, epsilon = 0.1, rho_query = 1, rho_ref = 1,
    solver = "internal", max_iter = 500, tol = 1e-7,
    feature_weights = c(1, 1, 0.001)
  )
  expect_gt(max(abs(base$plan - weighted$plan)), 1e-6)
})
