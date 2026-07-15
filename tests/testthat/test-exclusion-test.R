## ---------------------------------------------------------------------------
## Tests for somalign_exclusion_test() (Idea #8: anchor exclusion-restriction
## test). The null permutes each FEATURE of the orthogonal residual
## independently across anchors (destroying coherent cross-feature
## structure while preserving each feature's own variance) -- NOT a row
## permutation, which is a mathematical no-op (row permutation is an
## orthogonal transform and leaves every singular value of a matrix exactly
## invariant; see the comment above .somalign_permutation_null in
## R/anchored.R for the full argument).
## ---------------------------------------------------------------------------

make_exclusion_fixture <- function(seed = 42L) {
  withr::local_seed(seed)
  p <- 3L
  b <- c(1, 0, 0)     # batch direction
  cc <- c(0, 1, 0)    # orthogonal biology direction
  ref_data <- matrix(rnorm(40L * p, 0, 0.3), ncol = p,
                     dimnames = list(NULL, paste0("F", seq_len(p))))
  batch_mag <- 3.0
  anc_idx <- seq(11L, 30L)
  anc_old <- ref_data[anc_idx, , drop = FALSE]
  anc_new <- anc_old + matrix(batch_mag * b, length(anc_idx), p, byrow = TRUE)
  qry_data <- ref_data + matrix(batch_mag * b, nrow(ref_data), p, byrow = TRUE)

  ref <- somalign_train_reference(ref_data, grid = kohonen::somgrid(2L, 2L, "hexagonal"),
                                  rlen = 10L)
  qry <- somalign_query(qry_data, ref, grid = kohonen::somgrid(2L, 2L, "hexagonal"),
                        rlen = 10L)
  list(ref = ref, qry = qry, anc_old = anc_old, anc_new = anc_new,
       b = b, cc = cc, batch_mag = batch_mag)
}

test_that("exclusion test passes for pure-batch anchors (no coherent orthogonal structure)", {
  skip_if_not_installed("kohonen")
  fx <- make_exclusion_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  result <- somalign_exclusion_test(fit, n_perm = 499L, seed = 1L)
  expect_s3_class(result, "somalign_exclusion_test")
  expect_gte(result$p_value, 0.05)
  expect_true(result$verdict %in% c("pass", "warn"))
  expect_equal(result$n_anchors, nrow(fx$anc_old))
  expect_equal(result$rank_used, fit$anchors$batch_subspace$rank)
})

test_that("exclusion test fails when anchors carry a coherent orthogonal (biology) direction", {
  skip_if_not_installed("kohonen")
  fx <- make_exclusion_fixture()
  withr::local_seed(2L)
  delta <- rnorm(nrow(fx$anc_old), mean = 0, sd = 0.5)
  anc_new_contaminated <- fx$anc_new + outer(delta, fx$cc)
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = anc_new_contaminated,
    rho_anchor = 1, correction = "subspace"
  ))
  result <- somalign_exclusion_test(fit, n_perm = 499L, seed = 1L)
  expect_lt(result$p_value, 0.05)
  expect_equal(result$verdict, "fail")
  # F2 (the cc direction) should dominate the residual
  expect_equal(unname(which.max(result$feature_residual_norm)), 2L)
})

test_that("exclusion test errors on a cost_bonus (non-subspace) fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit_bonus <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
    rho_anchor = 1, correction = "cost_bonus"
  ))
  expect_error(somalign_exclusion_test(fit_bonus), "subspace.*both")
})

test_that("exclusion test errors on an unanchored fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit_plain <- somalign_fit(fx$qry, fx$ref)
  expect_error(somalign_exclusion_test(fit_plain), "somalign_anchored_fit")
})

test_that("exclusion test return structure is complete", {
  skip_if_not_installed("kohonen")
  fx <- make_exclusion_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  result <- somalign_exclusion_test(fit, n_perm = 99L, seed = 2L)
  expect_true(all(c("sv_observed", "sv_null", "p_value", "null_quantiles",
                    "relative_stat", "rank_used", "n_anchors", "n_features",
                    "verdict", "feature_residual_norm") %in% names(result)))
  expect_length(result$sv_null, 99L)
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$verdict %in% c("pass", "warn", "fail"))
})

test_that("print.somalign_exclusion_test does not error", {
  skip_if_not_installed("kohonen")
  fx <- make_exclusion_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  result <- somalign_exclusion_test(fit, n_perm = 49L, seed = 1L)
  expect_output(print(result))
})

test_that("row permutation of the residual would be powerless (regression guard)", {
  # Direct unit check on the documented mathematical fact that motivated the
  # column-permutation design: row-permuting (or sign-flipping) a matrix
  # leaves its singular values exactly invariant.
  set.seed(5)
  R <- matrix(rnorm(30), 10, 3)
  sv1 <- svd(R, nu = 0, nv = 0)$d
  sv2 <- svd(R[sample.int(10), , drop = FALSE], nu = 0, nv = 0)$d
  expect_equal(sv1, sv2, tolerance = 1e-10)
})
