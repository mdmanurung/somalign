test_that("mapping confidence is bounded in (0, 1] and higher for near cells", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)

  conf <- somalign_mapping_confidence(fit, k = 3L)
  expect_length(conf, nrow(fit$query$scaled_data))
  expect_true(all(conf > 0 & conf <= 1 + 1e-9))
  expect_equal(attr(conf, "k"), 3L)
  expect_true(attr(conf, "reference_scale") > 0)
})

test_that("cells on the reference score higher than a distant cell (monotonicity)", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  cb <- fit$reference$codebook

  # Query = the reference nodes themselves (should map with high confidence)
  # plus one cell placed far outside the reference (should map poorly).
  far <- cb[1, ] + 50
  fit$query$scaled_data <- rbind(cb, matrix(far, nrow = 1,
                                            dimnames = list(NULL, colnames(cb))))
  fit$query$sample_id <- as.character(seq_len(nrow(fit$query$scaled_data)))

  conf <- somalign_mapping_confidence(fit, k = 3L)
  n <- nrow(cb)
  expect_gt(mean(conf[seq_len(n)]), conf[n + 1L])   # near > far on average
  expect_lt(conf[n + 1L], 0.2)                       # far cell is clearly low
  # the truly-outside cell is the least confident of all
  expect_lte(conf[n + 1L], min(conf[seq_len(n)]))
})

test_that("k is clamped to the number of reference nodes", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  conf <- somalign_mapping_confidence(fit, k = 10000L)
  expect_equal(attr(conf, "k"), nrow(fit$reference$codebook))
  expect_true(all(is.finite(conf)))
})
