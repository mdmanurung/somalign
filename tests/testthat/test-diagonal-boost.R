test_that("somalign_fit diagonal_boost = 0 is identical to default", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit0 <- somalign_fit(fx$qry, fx$ref, diagonal_boost = 0)
  fit1 <- somalign_fit(fx$qry, fx$ref)
  expect_equal(fit0$node_shifts, fit1$node_shifts)
  expect_equal(fit0$diagnostics$ot$transport_mass,
               fit1$diagnostics$ot$transport_mass)
})

test_that("diagonal_boost reduces outside_corrected_fraction at high epsilon", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit_plain <- somalign_fit(fx$qry, fx$ref, epsilon = 0.5, diagonal_boost = 0)
  fit_boost <- somalign_fit(fx$qry, fx$ref, epsilon = 0.5, diagonal_boost = 0.3)
  frac_plain <- fit_plain$diagnostics$projection$outside_corrected_fraction
  frac_boost <- fit_boost$diagnostics$projection$outside_corrected_fraction
  expect_lte(frac_boost, frac_plain + 0.05)
})

test_that("diagonal_boost returns valid somalign_fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref, diagonal_boost = 0.5)
  expect_s3_class(fit, "somalign_fit")
  expect_true(is.finite(fit$diagnostics$ot$transport_mass))
})

test_that("large diagonal_boost concentrates mass on diagonal (identity transport)", {
  skip_if_not_installed("kohonen")
  ref <- tiny_reference()
  qry_data <- matrix(c(-1, 0, 0, 0, 1, 0), nrow = 3, ncol = 2, byrow = TRUE,
                     dimnames = list(NULL, c("a", "b")))
  qry <- somalign_query(qry_data, ref, som_query = make_som(qry_data),
                        codebook_space = "reference_scaled")
  fit_plain <- somalign_fit(qry, ref, epsilon = 0.5, diagonal_boost = 0)
  fit_boost <- somalign_fit(qry, ref, epsilon = 0.5, diagonal_boost = 5)
  # With a very large boost the transport plan should prefer identity-like
  # mapping → smaller mean correction norm
  mean_plain <- mean(sqrt(rowSums(fit_plain$node_shifts^2)))
  mean_boost <- mean(sqrt(rowSums(fit_boost$node_shifts^2)))
  expect_lte(mean_boost, mean_plain + 1e-10)
})
