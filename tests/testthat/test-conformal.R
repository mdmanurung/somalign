# Synthetic multiclass probabilities where the true class is favoured but the
# scores are deliberately imperfect. Split conformal guarantees coverage
# regardless of how well-calibrated the scores are, so these tests exercise the
# guarantee, not the quality of the underlying classifier.
.make_conformal_data <- function(n, classes, seed) {
  withr::with_seed(seed, {
    y <- sample(classes, n, replace = TRUE)
    p <- matrix(stats::runif(n * length(classes), 0, 0.5), n, length(classes),
                dimnames = list(NULL, classes))
    idx <- cbind(seq_len(n), match(y, classes))
    p[idx] <- p[idx] + 0.5
    list(p = p / rowSums(p), y = y)
  })
}

.coverage <- function(cp, truth) {
  mean(vapply(seq_along(truth), function(i) cp$sets[i, truth[i]], logical(1)))
}

test_that("split-conformal sets achieve marginal coverage >= 1 - alpha", {
  classes <- c("Tcell", "Bcell", "Mono", "NK")
  cal  <- .make_conformal_data(1500L, classes, seed = 1L)
  test <- .make_conformal_data(3000L, classes, seed = 2L)
  for (a in c(0.1, 0.2)) {
    cp <- somalign_conformal_labels(test$p, cal$p, cal$y, alpha = a)
    cov <- .coverage(cp, test$y)
    expect_gte(cov, 1 - a - 0.03)   # small finite-sample slack below the guarantee
    expect_lte(cov, 1)
    expect_equal(nrow(cp$sets), nrow(test$p))
    expect_equal(cp$set_size, as.integer(rowSums(cp$sets)))
  }
})

test_that("class-conditional (Mondrian) conformal also covers marginally", {
  classes <- c("Tcell", "Bcell", "Mono")
  cal  <- .make_conformal_data(1800L, classes, seed = 3L)
  test <- .make_conformal_data(3000L, classes, seed = 4L)
  cp <- somalign_conformal_labels(test$p, cal$p, cal$y, alpha = 0.1,
                                  class_conditional = TRUE)
  expect_length(cp$threshold, length(classes))
  expect_gte(.coverage(cp, test$y), 1 - 0.1 - 0.03)
})

test_that("class-conditional conformal warns when a calibration class is too thin", {
  classes <- c("A", "B", "C")
  y <- c(rep("A", 100L), rep("B", 100L), rep("C", 3L))     # class C < ceiling(1/0.1)
  withr::with_seed(9L, {
    p <- matrix(stats::runif(length(y) * 3, 0, 0.4), length(y), 3,
                dimnames = list(NULL, classes))
    i <- cbind(seq_along(y), match(y, classes))
    p[i] <- p[i] + 0.6
    p <- p / rowSums(p)
  })
  expect_warning(
    somalign_conformal_labels(p[1:10, , drop = FALSE], p, y, alpha = 0.1,
                              class_conditional = TRUE),
    "fewer than")
})

test_that("smaller alpha yields larger sets (fewer abstentions, more coverage)", {
  classes <- c("A", "B", "C")
  cal  <- .make_conformal_data(1200L, classes, seed = 5L)
  test <- .make_conformal_data(1200L, classes, seed = 6L)
  cp_tight <- somalign_conformal_labels(test$p, cal$p, cal$y, alpha = 0.30)
  cp_loose <- somalign_conformal_labels(test$p, cal$p, cal$y, alpha = 0.05)
  expect_gte(mean(cp_loose$set_size), mean(cp_tight$set_size))
  expect_gte(.coverage(cp_loose, test$y), .coverage(cp_tight, test$y) - 1e-9)
})

test_that("input validation rejects malformed arguments", {
  p <- matrix(0.5, 4, 2, dimnames = list(NULL, c("A", "B")))
  expect_error(somalign_conformal_labels(p, p, c("A", "B", "A", "B"), alpha = 0),
               "alpha")
  expect_error(somalign_conformal_labels(p, p, c("A", "B", "A", "C")),
               "must be a class column")
  q <- matrix(0.5, 4, 2, dimnames = list(NULL, c("A", "X")))
  expect_error(somalign_conformal_labels(q, p, c("A", "B", "A", "B")),
               "same class columns")
  u <- matrix(0.5, 4, 2)  # no colnames
  expect_error(somalign_conformal_labels(u, u, c("A", "B", "A", "B")),
               "class names")
  na <- matrix(c(0.5, NA, 0.5, 0.5), 2, 2, dimnames = list(NULL, c("A", "B")))
  expect_error(somalign_conformal_labels(na, na, c("A", "B")), "finite")
})

test_that("print.somalign_conformal is invisible and reports the breakdown", {
  classes <- c("A", "B", "C")
  cal  <- .make_conformal_data(400L, classes, seed = 7L)
  test <- .make_conformal_data(100L, classes, seed = 8L)
  cp <- somalign_conformal_labels(test$p, cal$p, cal$y)
  out <- withVisible(print(cp))
  expect_false(out$visible)
  expect_s3_class(cp, "somalign_conformal")
  # every cell is abstain XOR singleton XOR ambiguous
  expect_equal(sum(cp$abstain) + sum(cp$set_size == 1L) + sum(cp$set_size >= 2L),
               nrow(cp$sets))
})
