## ---------------------------------------------------------------------------
## Tests for somalign_subspace_sensitivity() (Idea #7: bootstrap sensitivity
## of the anchor-estimated batch subspace). Bootstraps D WITH replacement
## (a valid nonparametric bootstrap, unlike the WITHOUT-replacement
## permutation used by somalign_exclusion_test()'s null).
## ---------------------------------------------------------------------------

make_sensitivity_fixture <- function(seed = 42L, natural_sd = 0.3, batch_mag = 3.0) {
  withr::local_seed(seed)
  p <- 3L
  b <- c(1, 0, 0)
  ref_data <- matrix(rnorm(40L * p, 0, natural_sd), ncol = p,
                     dimnames = list(NULL, paste0("F", seq_len(p))))
  anc_idx <- seq(11L, 30L)
  anc_old <- ref_data[anc_idx, , drop = FALSE]
  anc_new <- anc_old + matrix(batch_mag * b, length(anc_idx), p, byrow = TRUE)
  qry_data <- ref_data + matrix(batch_mag * b, nrow(ref_data), p, byrow = TRUE)
  ref <- somalign_train_reference(ref_data, grid = kohonen::somgrid(2L, 2L, "hexagonal"), rlen = 10L)
  qry <- somalign_query(qry_data, ref, grid = kohonen::somgrid(2L, 2L, "hexagonal"), rlen = 10L)
  list(ref = ref, qry = qry, anc_old = anc_old, anc_new = anc_new)
}

test_that("tight anchors produce a stable subspace and large tipping angles", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  sens <- somalign_subspace_sensitivity(fit, n_boot = 100L, seed = 1L)
  expect_s3_class(sens, "somalign_subspace_sensitivity")
  expect_true(all(sens$subspace_angles < 10, na.rm = TRUE))
  ta <- sens$tipping_angle_deg[!is.na(sens$tipping_angle_deg)]
  expect_true(length(ta) > 0)
  expect_true(all(ta > 30))
})

test_that("small anchor count triggers the low-power warning", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old[1:2, , drop = FALSE],
    anchor_new = fx$anc_new[1:2, , drop = FALSE],
    rho_anchor = 1, correction = "subspace"
  ))
  expect_warning(
    somalign_subspace_sensitivity(fit, n_boot = 20L, seed = 1L),
    "n_anchors"
  )
})

test_that("somalign_subspace_sensitivity errors on a cost_bonus fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
    rho_anchor = 1, correction = "cost_bonus"
  ))
  expect_error(somalign_subspace_sensitivity(fit), "subspace.*both")
})

test_that("somalign_subspace_sensitivity errors on an unanchored fit", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  expect_error(somalign_subspace_sensitivity(fit), "somalign_anchored_fit")
})

test_that("return structure has all documented fields with correct dimensions", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  sens <- somalign_subspace_sensitivity(fit, n_boot = 30L, seed = 2L)
  m <- nrow(fit$node_shifts)
  p <- ncol(fit$node_shifts)
  r <- fit$anchors$batch_subspace$rank

  expect_equal(dim(sens$node_correction_ci), c(m, 2L))
  expect_equal(colnames(sens$node_correction_ci), c("lower", "upper"))
  expect_equal(dim(sens$node_shift_ci), c(m, p, 2L))
  expect_equal(dim(sens$subspace_angles), c(30L, r))
  expect_equal(length(sens$tipping_angle_deg), m)
  expect_equal(length(sens$anchor_leverage), nrow(fx$anc_old))
  expect_equal(sens$n_boot, 30L)
  expect_equal(sens$subspace_rank, r)
  expect_equal(sens$n_anchors, nrow(fx$anc_old))
})

test_that("rank-1 tipping angle matches the analytic formula", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture(seed = 7L)
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  sens <- somalign_subspace_sensitivity(fit, n_boot = 50L, seed = 1L)
  V <- fit$anchors$batch_subspace$V
  expect_equal(ncol(V), 1L)
  s_raw <- somalign:::.somalign_recover_raw_shifts(fit)
  allowed <- attr(fit$node_shifts, "correction_allowed")
  i <- which(allowed & !is.na(sens$tipping_angle_deg))[1]
  s_hat <- s_raw[i, ] / sqrt(sum(s_raw[i, ]^2))
  expected_angle <- asin(min(abs(sum(s_hat * V[, 1])), 1)) * 180 / pi
  expect_equal(sens$tipping_angle_deg[i], expected_angle, tolerance = 1e-6)
})

test_that("print.somalign_subspace_sensitivity does not error", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  sens <- somalign_subspace_sensitivity(fit, n_boot = 20L, seed = 1L)
  expect_output(print(sens))
})

test_that("seed does not leak into the caller's global RNG state", {
  skip_if_not_installed("kohonen")
  fx <- make_sensitivity_fixture()
  fit <- suppressWarnings(somalign_fit_anchored(
    fx$qry, fx$ref, anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  ))
  set.seed(123)
  before <- runif(1)
  set.seed(123)
  somalign_subspace_sensitivity(fit, n_boot = 20L, seed = 999L)
  after <- runif(1)
  expect_equal(before, after)
})
