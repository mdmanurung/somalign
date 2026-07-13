## Tests for somalign_fit_anchored: anchor-regularized OT.
## make_anchored_fixture() lives in helper-fixtures.R.

test_that("somalign_fit_anchored returns correct classes and anchors slot", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_anchored(fx$qry, fx$ref,
                                anchor_old = fx$anchor_old,
                                anchor_new = fx$anchor_new,
                                rho_anchor = 1.0)
  expect_s3_class(fit, "somalign_anchored_fit")
  expect_s3_class(fit, "somalign_fit")
  expect_true(!is.null(fit$anchors))
  expect_equal(fit$anchors$n_anchors, nrow(fx$anchor_old))
  expect_equal(fit$anchors$rho_anchor, 1.0)
  expect_true(is.integer(fit$anchors$nodes_covered))
  expect_gte(fit$anchors$coverage_fraction, 0)
  expect_lte(fit$anchors$coverage_fraction, 1)
})

test_that("rho_anchor = 0 gives same plan as somalign_fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture(seed = 2L)
  fit_plain    <- somalign_fit(fx$qry, fx$ref)
  fit_anchored <- somalign_fit_anchored(fx$qry, fx$ref,
                                         anchor_old = fx$anchor_old,
                                         anchor_new = fx$anchor_new,
                                         rho_anchor = 0.0)
  expect_equal(fit_anchored$transport_plan, fit_plain$transport_plan,
               tolerance = 1e-10)
})

test_that("anchor regularization increases diagonal-mass when shift is known", {
  skip_if_not_installed("kohonen")
  # Reference SOM supplied as query SOM so codebooks match exactly;
  # anchors confirm the identity mapping.  rho_anchor > 0 should concentrate
  # even more mass on the diagonal than the baseline.
  fx <- make_anchored_fixture(seed = 3L)
  ref <- fx$ref

  qry_id <- somalign_query(fx$ref_data, ref,
                            som_query = ref$som_ref,
                            codebook_space = "reference_scaled")

  fit_plain <- somalign_fit(qry_id, ref, epsilon = 0.3)
  # Anchors all confirm identity: anchor_old == anchor_new in reference space
  fit_anc   <- somalign_fit_anchored(qry_id, ref,
                                      anchor_old = fx$ref_data[fx$anc_idx, ],
                                      anchor_new = fx$ref_data[fx$anc_idx, ],
                                      rho_anchor = 2.0,
                                      epsilon    = 0.3)

  plan_plain <- fit_plain$transport_plan
  plan_anc   <- fit_anc$transport_plan
  n <- min(nrow(plan_plain), ncol(plan_plain))
  diag_plain <- sum(diag(plan_plain[seq_len(n), seq_len(n)]))
  diag_anc   <- sum(diag(plan_anc[seq_len(n), seq_len(n)]))
  expect_gte(diag_anc, diag_plain - 1e-10)
})

test_that("somalign_results and somalign_diagnostics work on anchored fit", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture(seed = 4L)
  fit <- somalign_fit_anchored(fx$qry, fx$ref,
                                anchor_old = fx$anchor_old,
                                anchor_new = fx$anchor_new,
                                rho_anchor = 1.0)
  res  <- somalign_results(fit)
  diag <- somalign_diagnostics(fit)
  expect_s3_class(res, "data.frame")
  expect_true("corrected_som_unit" %in% names(res))
  expect_true(!is.null(diag$solver$converged))
})

test_that("somalign_fit_anchored validates mismatched anchor dimensions", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture(seed = 5L)
  bad_new <- fx$anchor_new[seq_len(nrow(fx$anchor_new) - 1L), , drop = FALSE]
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                           anchor_old = fx$anchor_old,
                           anchor_new = bad_new,
                           rho_anchor = 1.0),
    "same number of rows"
  )
})

test_that("somalign_fit_anchored validates wrong feature columns", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture(seed = 6L)
  bad <- fx$anchor_old
  colnames(bad) <- paste0("X", seq_len(ncol(bad)))
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                           anchor_old = bad,
                           anchor_new = fx$anchor_new,
                           rho_anchor = 1.0),
    "Missing features"
  )
})

test_that("somalign_fit_anchored rejects invalid rho_anchor", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture(seed = 7L)
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                           anchor_old = fx$anchor_old,
                           anchor_new = fx$anchor_new,
                           rho_anchor = -1.0),
    "rho_anchor"
  )
})
