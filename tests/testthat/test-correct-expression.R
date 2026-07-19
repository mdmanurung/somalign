## ---------------------------------------------------------------------------
## Tests for somalign_correct_expression(): subspace-restricted, cell-level
## smoothed batch correction of query marker expression.
## make_subspace_fixture() is defined in helper-fixtures.R.
## ---------------------------------------------------------------------------

subspace_fit <- function(seed = 42L) {
  fx <- make_subspace_fixture(seed)
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    correction = "subspace"
  )
  list(fx = fx, fit = fit)
}

# (a) The correction pulls the batch-direction coordinate toward the reference.
test_that("correction reduces the batch-direction component", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, units = "scaled", bandwidth = 0.5)

  before <- mean(abs(sf$fit$query$scaled_data %*% sf$fx$b))
  after  <- mean(abs(corr %*% sf$fx$b))
  expect_lt(after, before)
})

# (b) Variation orthogonal to the batch subspace is preserved exactly.
test_that("orthogonal biology is preserved to machine precision", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, units = "scaled", bandwidth = 0.5)

  shift <- unclass(corr) - sf$fit$query$scaled_data
  V <- sf$fit$anchors$batch_subspace$V
  ortho <- shift - shift %*% V %*% t(V)
  expect_lt(norm(ortho, "F"), 1e-10)

  # the biology offset along cc (for the sub_idx cells) is untouched
  sub <- sf$fx$sub_idx
  cc_before <- sf$fit$query$scaled_data[sub, ] %*% sf$fx$cc
  cc_after  <- corr[sub, ] %*% sf$fx$cc
  expect_equal(as.vector(cc_after), as.vector(cc_before), tolerance = 1e-10)
})

# (c) Smoothing yields a more continuous field than the piecewise-constant one.
test_that("smoothing does not increase shift-magnitude dispersion", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  smooth <- somalign_correct_expression(sf$fit, units = "scaled",
                                        smooth = TRUE, bandwidth = 0.5)
  hard   <- somalign_correct_expression(sf$fit, units = "scaled", smooth = FALSE)

  s_smooth <- unclass(smooth) - sf$fit$query$scaled_data
  s_hard   <- unclass(hard)   - sf$fit$query$scaled_data
  expect_lte(var(rowSums(s_smooth^2)), var(rowSums(s_hard^2)))
})

# (d) A fit without a batch subspace is rejected with an informative message.
test_that("plain somalign_fit errors (no batch subspace)", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit_plain <- somalign_fit(fx$qry, fx$ref)
  expect_error(
    somalign_correct_expression(fit_plain),
    regexp = "batch subspace"
  )
})

# (e) Raw output is the exact back-transform of the scaled output.
test_that("raw and scaled outputs are consistent", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr_raw    <- somalign_correct_expression(sf$fit, units = "raw",    bandwidth = 0.5)
  corr_scaled <- somalign_correct_expression(sf$fit, units = "scaled", bandwidth = 0.5)

  center <- sf$fit$reference$center
  scale  <- sf$fit$reference$scale
  rebuilt <- sweep(sweep(unclass(corr_scaled), 2, scale, "*"), 2, center, "+")
  expect_equal(as.vector(unclass(corr_raw)), as.vector(rebuilt), tolerance = 1e-12)

  expect_equal(attr(corr_raw, "units"), "raw")
  expect_equal(attr(corr_scaled, "units"), "scaled")
})

# (f) k is silently clamped to the number of SOM nodes.
test_that("k larger than the SOM node count is clamped", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  m <- nrow(sf$fit$query$codebook)
  expect_no_error(
    corr <- somalign_correct_expression(sf$fit, units = "scaled",
                                        k = 8L, bandwidth = 0.5)
  )
  expect_equal(attr(corr, "k"), min(8L, m))
})

# confidence_gate = FALSE exercises the ungated weighting branch.
test_that("confidence_gate = FALSE stays subspace-restricted and finite", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, units = "scaled",
                                      confidence_gate = FALSE, bandwidth = 0.5)
  expect_true(all(is.finite(corr)))
  shift <- unclass(corr) - sf$fit$query$scaled_data
  V <- sf$fit$anchors$batch_subspace$V
  expect_lt(norm(shift - shift %*% V %*% t(V), "F"), 1e-10)
})

# Default bandwidth (NULL) derives a positive scale-adaptive value and corrects.
test_that("default bandwidth is derived and reduces the batch component", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, units = "scaled")   # bandwidth = NULL
  h <- attr(corr, "bandwidth")
  expect_true(is.finite(h) && h > 0)

  before <- mean(abs(sf$fit$query$scaled_data %*% sf$fx$b))
  after  <- mean(abs(corr %*% sf$fx$b))
  expect_lt(after, before)
})

# Fused chunking must be invariant to chunk_size.
test_that("corrected expression is invariant to chunk_size", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  a <- somalign_correct_expression(sf$fit, units = "scaled", bandwidth = 0.5,
                                   chunk_size = 1000000L)
  b <- somalign_correct_expression(sf$fit, units = "scaled", bandwidth = 0.5,
                                   chunk_size = 9L)
  expect_equal(unclass(a), unclass(b), tolerance = 1e-12, ignore_attr = TRUE)
})

# A two-pass fit carries no anchor-estimated subspace: somalign_correct_expression
# applies the full documented two-pass correction (residual + global), NOT a
# subspace-confined one. The $two_pass$batch_subspace diagnostic must not be used
# to censor the correction (it is documented as descriptive-only).
test_that("two-pass fit correction is applied in full and reduces the batch shift", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit_tp <- somalign_fit_two_pass(fx$qry, fx$ref)
  corr <- somalign_correct_expression(fit_tp, units = "scaled", bandwidth = 0.5)

  expect_s3_class(corr, "somalign_corrected_expression")
  expect_true(all(is.finite(corr)))

  # The correction reduces the batch-direction component of the query.
  before <- mean(abs(fit_tp$query$scaled_data %*% fx$b))
  after  <- mean(abs(corr %*% fx$b))
  expect_lt(after, before)

  # Regression lock: with smooth = FALSE each cell takes its nearest node's FULL
  # shift (disallowed nodes zeroed). If the correction were re-projected onto the
  # two-pass subspace this identity would break.
  ns <- fit_tp$node_shifts
  allowed <- attr(ns, "correction_allowed")
  if (is.null(allowed)) allowed <- rep(TRUE, nrow(ns))
  su <- fit_tp$query$sample_unit
  expected <- ns[su, , drop = FALSE]
  expected[!allowed[su], ] <- 0
  full <- somalign_correct_expression(fit_tp, units = "scaled", smooth = FALSE)
  expect_equal(unclass(full) - fit_tp$query$scaled_data, expected,
               ignore_attr = TRUE, tolerance = 1e-10)
})

# Output shape, dimnames, and class.
test_that("output has cell x marker shape with expected names and class", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, bandwidth = 0.5)
  expect_s3_class(corr, "somalign_corrected_expression")
  expect_true(is.matrix(corr))
  expect_equal(nrow(corr), nrow(sf$fit$query$scaled_data))
  expect_equal(colnames(corr), colnames(sf$fit$query$scaled_data))
  expect_equal(rownames(corr), sf$fit$query$sample_id)
})

test_that("print.somalign_corrected_expression does not error", {
  skip_if_not_installed("kohonen")
  sf <- subspace_fit()
  corr <- somalign_correct_expression(sf$fit, bandwidth = 0.5)
  expect_output(print(corr), "somalign_corrected_expression")
})
