test_that("somalign_quantile_normalize preserves dimensions", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 0.999)
  expect_equal(dim(out), dim(fx$qry_data))
  expect_equal(colnames(out), fx$ref$features)
})

test_that("somalign_quantile_normalize probs=0.999 maps empirical quantile to ~1.0", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 0.999)
  for (j in seq_len(ncol(out))) {
    q_out <- unname(stats::quantile(out[, j], probs = 0.999, na.rm = TRUE))
    expect_equal(q_out, 1.0, tolerance = 1e-9,
                 label = paste0("column ", colnames(out)[j]))
  }
})

test_that("somalign_quantile_normalize warns on zero quantile and leaves column unchanged", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  # Create data with one all-zero column to trigger the zero-quantile guard
  data_zero <- fx$qry_data
  data_zero[, 1L] <- 0
  expect_warning(
    out <- somalign_quantile_normalize(data_zero, fx$ref, probs = 0.999),
    regexp = "zero or non-finite"
  )
  # The all-zero column is divided by the substituted 1, so it stays unchanged
  expect_equal(out[, 1L], data_zero[, 1L])
  # Non-zero columns are still normalised
  q2 <- unname(stats::quantile(out[, 2L], probs = 0.999, na.rm = TRUE))
  expect_equal(q2, 1.0, tolerance = 1e-9)
})

test_that("somalign_quantile_normalize rejects invalid probs values", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  expect_error(somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 0),
               regexp = "strictly between 0 and 1")
  expect_error(somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 1),
               regexp = "strictly between 0 and 1")
  expect_error(somalign_quantile_normalize(fx$qry_data, fx$ref, probs = -0.1),
               regexp = "strictly between 0 and 1")
  expect_error(somalign_quantile_normalize(fx$qry_data, fx$ref, probs = "a"),
               regexp = "strictly between 0 and 1")
})

test_that("somalign_quantile_normalize output passes to somalign_query without error", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  normed <- somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 0.999)
  qry <- somalign_query(
    normed, fx$ref,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  expect_s3_class(qry, "somalign_query")
  expect_equal(colnames(qry$scaled_data), fx$ref$features)
})

test_that("somalign_quantile_normalize probs=0.5 maps median of each column to ~1.0", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_quantile_normalize(fx$qry_data, fx$ref, probs = 0.5)
  for (j in seq_len(ncol(out))) {
    med_out <- unname(stats::quantile(out[, j], probs = 0.5, na.rm = TRUE))
    expect_equal(med_out, 1.0, tolerance = 1e-9,
                 label = paste0("column ", colnames(out)[j]))
  }
})
