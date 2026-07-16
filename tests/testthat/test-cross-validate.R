## Held-out cross-validation for label transfer (item 2).

skip_if_not_installed("kohonen")

make_cv_data <- function(n = 300L, sep = 3, seed = 1L) {
  withr::local_seed(seed)
  x <- rbind(
    matrix(rnorm(n * 3, -sep, 0.5), ncol = 3),
    matrix(rnorm(n * 3,  sep, 0.5), ncol = 3)
  )
  colnames(x) <- paste0("f", seq_len(3))
  list(x = x, labels = rep(c("low", "high"), each = n))
}

test_that("cross-validation recovers well-separated labels near-perfectly", {
  d <- make_cv_data()
  cv <- somalign_cross_validate(
    d$x, d$labels, grid = kohonen::somgrid(3, 3, "hexagonal"),
    k = 3, rlen = 20
  )
  expect_s3_class(cv, "somalign_cross_validation")
  expect_gt(cv$metrics$accuracy, 0.95)
  expect_gt(cv$metrics$mcc, 0.9)
  expect_equal(nrow(cv$per_fold), 3L)
  expect_equal(nrow(cv$predictions), nrow(d$x))
})

test_that("cross-validation returns metrics and calibration objects", {
  d <- make_cv_data(n = 200L)
  cv <- somalign_cross_validate(
    d$x, d$labels, grid = kohonen::somgrid(2, 2, "hexagonal"),
    k = 2, rlen = 15
  )
  expect_s3_class(cv$metrics, "somalign_label_metrics")
  expect_s3_class(cv$calibration, "somalign_calibration")
  expect_true(all(c("accuracy", "macro_f1", "mcc", "coverage") %in%
                    names(cv$per_fold)))
  op <- withVisible(print(cv))
  expect_false(op$visible)
})

test_that("stratified folds place every class in every fold", {
  labels <- rep(c("A", "B", "C"), times = c(30, 30, 30))
  withr::local_seed(1)
  folds <- somalign:::.somalign_stratified_folds(labels, k = 5)
  tab <- table(folds, labels)
  expect_true(all(tab > 0))          # no empty (fold, class) cell
  expect_equal(length(unique(folds)), 5L)
})

test_that("cross-validation validates label length", {
  d <- make_cv_data(n = 20L)
  expect_error(
    somalign_cross_validate(d$x, d$labels[1:5],
                            grid = kohonen::somgrid(2, 2, "hexagonal")),
    "one entry per row"
  )
})
