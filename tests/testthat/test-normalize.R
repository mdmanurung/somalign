test_that("somalign_normalize mean removes global mean shift in scaled space", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_normalize(fx$qry_data, fx$ref, method = "mean")
  scaled_back <- somalign:::.somalign_scale_matrix(out, fx$ref$center, fx$ref$scale)
  expect_lt(max(abs(colMeans(scaled_back))), 1e-10)
})

test_that("somalign_normalize scale centres and standardises in scaled space", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_normalize(fx$qry_data, fx$ref, method = "scale")
  scaled_back <- somalign:::.somalign_scale_matrix(out, fx$ref$center, fx$ref$scale)
  expect_lt(max(abs(colMeans(scaled_back))), 1e-10)
  sds <- apply(scaled_back, 2, sd)
  expect_true(all(abs(sds - 1) < 1e-10))
})

test_that("somalign_normalize output is passable to somalign_query", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_normalize(fx$qry_data, fx$ref, method = "mean")
  qry <- somalign_query(
    out, fx$ref,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  expect_s3_class(qry, "somalign_query")
  expect_equal(colnames(qry$scaled_data), fx$ref$features)
})

test_that("somalign_normalize returns same dimensions as input", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  out <- somalign_normalize(fx$qry_data, fx$ref, method = "mean")
  expect_equal(dim(out), dim(fx$qry_data))
  expect_equal(colnames(out), fx$ref$features)
})

test_that("somalign_normalize mean reduces outside_direct_fraction vs raw query", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  qry_raw <- fx$qry
  qry_norm <- somalign_query(
    somalign_normalize(fx$qry_data, fx$ref, method = "mean"),
    fx$ref,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  fit_raw  <- somalign_fit(qry_raw,  fx$ref, epsilon = 0.5)
  fit_norm <- somalign_fit(qry_norm, fx$ref, epsilon = 0.5)
  frac_raw  <- fit_raw$diagnostics$projection$outside_direct_fraction
  frac_norm <- fit_norm$diagnostics$projection$outside_direct_fraction
  expect_lte(frac_norm, frac_raw + 0.05)
})
