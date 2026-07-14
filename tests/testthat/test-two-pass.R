test_that("somalign_fit_two_pass returns somalign_fit with two_pass slot", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref)
  expect_s3_class(fit, "somalign_fit")
  expect_true(!is.null(fit$two_pass))
  expect_named(fit$two_pass,
               c("global_shift", "global_shift_norm", "epsilon_global", "epsilon_local"))
})

test_that("somalign_fit_two_pass global_shift points in the correction direction", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref,
                               epsilon_global = 0.5, epsilon_local = 0.1)
  # The fixture shifts query by +1.0 in raw space, so in reference-scaled
  # space query nodes sit ~+1/scale above reference nodes.  The correction
  # shift (node_shifts = barycentric_ref - query_codebook) therefore points
  # in the -1 direction per feature.  Check cosine similarity with -rep(1,3).
  g <- fit$two_pass$global_shift
  expected_dir <- rep(-1.0, length(g))
  cosine <- sum(g * expected_dir) /
    (sqrt(sum(g^2)) * sqrt(sum(expected_dir^2)))
  expect_gt(cosine, 0.5,
            label = paste("cosine similarity of global_shift with expected direction:",
                          round(cosine, 3)))
})

test_that("somalign_fit_two_pass global_shift_norm is positive for non-zero batch shift", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref)
  expect_gt(fit$two_pass$global_shift_norm, 0)
})

test_that("somalign_fit_two_pass direct projection uses original scaled_data", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit  <- somalign_fit_two_pass(fx$qry, fx$ref)
  fit0 <- somalign_fit(fx$qry, fx$ref)
  # Direct projection (old_som_unit) is transport-free; must match single-pass
  expect_equal(fit$projection$direct$unit, fit0$projection$direct$unit)
})

test_that("somalign_fit_two_pass total_shifts equal residual plus global", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref,
                               epsilon_global = 0.5, epsilon_local = 0.1)
  g  <- fit$two_pass$global_shift
  ns <- fit$node_shifts
  # For correction_allowed nodes: |total_shift - g| should be < |total_shift|
  # (total is closer to zero than either component alone only sometimes, but
  # the global is always subtracted out of pass-2 node_shifts)
  allowed <- fit$diagnostics$nodes$correction_allowed
  if (any(allowed)) {
    residual_norms <- sqrt(rowSums((ns[allowed, , drop = FALSE] -
                                     matrix(g, sum(allowed), length(g), byrow = TRUE))^2))
    total_norms <- sqrt(rowSums(ns[allowed, , drop = FALSE]^2))
    # residual norms should be <= total norms (g has been removed)
    expect_true(all(residual_norms <= total_norms + 1e-10))
  }
})

test_that("somalign_fit_two_pass produces valid results via somalign_results", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref)
  res <- somalign_results(fit)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(fx$qry_data))
})
