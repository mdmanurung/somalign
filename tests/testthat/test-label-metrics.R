## Label-transfer metrics and calibration (item 2).

test_that("multiclass MCC reduces to standard MCC in the 2-class case", {
  TP <- 5; TN <- 6; FP <- 2; FN <- 3
  std <- (TP * TN - FP * FN) /
    sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  truth <- c(rep("pos", TP + FN), rep("neg", TN + FP))
  pred  <- c(rep("pos", TP), rep("neg", FN), rep("neg", TN), rep("pos", FP))
  m <- somalign_label_metrics(pred, truth)
  expect_equal(m$mcc, std, tolerance = 1e-12)
})

test_that("accuracy, macro-F1 and per-class stats match hand computation", {
  truth <- rep(c("A", "B", "C"), each = 10)
  pred  <- truth
  pred[c(1, 12, 25)] <- c("B", "A", "A")     # 3 errors of 30
  m <- somalign_label_metrics(pred, truth)
  expect_equal(m$accuracy, 27 / 30)
  expect_equal(nrow(m$per_class), 3L)
  expect_equal(m$per_class$support, c(10, 10, 10))
  expect_true(m$macro_f1 > 0.85 && m$macro_f1 < 1)
  expect_equal(sum(m$confusion), 30)
})

test_that("perfect predictions give accuracy = MCC = macro_f1 = 1", {
  truth <- rep(c("A", "B", "C"), each = 5)
  m <- somalign_label_metrics(truth, truth)
  expect_equal(m$accuracy, 1)
  expect_equal(m$mcc, 1)
  expect_equal(m$macro_f1, 1)
})

test_that("single-class predictions yield MCC = 0 (undefined denominator)", {
  truth <- rep("A", 10)
  pred  <- rep("A", 10)
  m <- somalign_label_metrics(pred, truth)
  expect_equal(m$accuracy, 1)
  expect_equal(m$mcc, 0)             # denominator is 0 -> defined as 0
})

test_that("accepted gates scoring and coverage", {
  truth <- rep(c("A", "B"), each = 5)
  pred  <- truth
  accepted <- c(rep(TRUE, 8), FALSE, FALSE)
  m <- somalign_label_metrics(pred, truth, accepted = accepted)
  expect_equal(m$n, 8L)
  expect_equal(m$coverage, 0.8)
  expect_equal(m$accuracy, 1)                 # scored subset all correct
  expect_equal(m$accuracy_all, 0.8)           # abstentions counted wrong
})

test_that("NA predictions are treated as abstentions", {
  truth <- rep(c("A", "B"), each = 5)
  pred  <- truth
  pred[c(1, 2)] <- NA
  m <- somalign_label_metrics(pred, truth)
  expect_equal(m$n, 8L)
  expect_equal(m$coverage, 0.8)
})

test_that("label metrics validate input lengths", {
  expect_error(somalign_label_metrics(c("A", "B"), c("A")), "same length")
})

test_that("calibration: perfectly calibrated scores have near-zero ECE", {
  withr::local_seed(42)
  score <- runif(5000)
  correct <- runif(5000) < score
  cal <- somalign_calibration(score, correct, n_bins = 10)
  expect_lt(cal$ece, 0.03)
  expect_true(cal$brier > 0 && cal$brier < 0.25)
})

test_that("calibration: overconfident scores have large ECE", {
  cal <- somalign_calibration(rep(0.99, 1000), rep(c(TRUE, FALSE), 500))
  expect_gt(cal$ece, 0.45)
  expect_gt(cal$mce, 0.45)
})

test_that("calibration errors on empty input and filters NA", {
  expect_error(somalign_calibration(numeric(0), logical(0)), "No non-missing")
  cal <- somalign_calibration(c(0.9, NA, 0.8), c(TRUE, TRUE, NA))
  expect_equal(cal$n, 1L)
  expect_equal(cal$n_total, 3L)
  expect_equal(cal$coverage, 1 / 3)
})

test_that("print methods return their objects invisibly", {
  m <- somalign_label_metrics(rep("A", 4), rep("A", 4))
  om <- withVisible(print(m))
  expect_false(om$visible)
  cal <- somalign_calibration(c(0.9, 0.8), c(TRUE, TRUE))
  oc <- withVisible(print(cal))
  expect_false(oc$visible)
})
