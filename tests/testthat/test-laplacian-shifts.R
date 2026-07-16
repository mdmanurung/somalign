## ---------------------------------------------------------------------------
## Tests for laplacian_lambda (Idea #5: graph-Laplacian smoothing of the
## node-shift field). Wired through the existing shift_transform hook in
## .somalign_finish_fit(); composes with the anchored subspace projector as
## smooth -> project.
## ---------------------------------------------------------------------------

test_that("laplacian_lambda = 0 produces an identical fit to the default", {
  skip_if_not_installed("kohonen")
  set.seed(42)
  mat <- matrix(rnorm(60), nrow = 30, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  fit0 <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 0)
  expect_equal(fit0$node_shifts, fit_lap$node_shifts)
  expect_equal(attr(fit0$node_shifts, "correction_allowed"),
               attr(fit_lap$node_shifts, "correction_allowed"))
})

test_that("correction_allowed attribute is preserved after Laplacian smoothing", {
  skip_if_not_installed("kohonen")
  set.seed(7)
  mat <- matrix(rnorm(80), nrow = 40, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  qry <- somalign_query(mat + 1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  fit0 <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 0.5)
  expect_identical(attr(fit0$node_shifts, "correction_allowed"),
                   attr(fit_lap$node_shifts, "correction_allowed"))
})

test_that("disallowed nodes keep an exact zero shift after smoothing", {
  skip_if_not_installed("kohonen")
  set.seed(7)
  mat <- matrix(rnorm(80), nrow = 40, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  qry <- somalign_query(mat + 1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  # min_match_fraction = 0.5 deterministically disallows at least one node
  # for this fixture/seed, so the assertion below is always exercised.
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 0.5, min_match_fraction = 0.5)
  allowed <- attr(fit_lap$node_shifts, "correction_allowed")
  expect_true(any(!allowed))
  expect_true(all(fit_lap$node_shifts[!allowed, ] == 0))
})

test_that("large laplacian_lambda collapses shifts toward their mass-weighted mean", {
  skip_if_not_installed("kohonen")
  set.seed(13)
  mat <- matrix(rnorm(80), nrow = 40, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  qry <- somalign_query(mat + 1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  fit0 <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 1e6)
  raw_var <- stats::var(fit0$node_shifts[, 1])
  smooth_var <- stats::var(fit_lap$node_shifts[, 1])
  expect_lt(smooth_var, raw_var * 0.1)
})

test_that("invalid laplacian_lambda raises an informative error", {
  skip_if_not_installed("kohonen")
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2, dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_error(somalign_fit(qry, ref, laplacian_lambda = -1), "non-negative")
})

test_that("laplacian_lambda errors clearly when the query SOM has no grid coordinates", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  # make_som() returns a SOM-like object with a codebook but no $grid -- the
  # scenario laplacian_lambda must reject with a clear error.
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))
  expect_error(somalign_fit(qry, ref, laplacian_lambda = 0.5), "grid")
})

test_that("laplacian_lambda composes with subspace mode: smooth then project", {
  skip_if_not_installed("kohonen")
  # Inlined fixture (mirrors make_subspace_fixture() in
  # test-anchored-subspace.R, which is not reliably visible from this file
  # when tests run individually): a known batch direction b plus orthogonal
  # biology cc, with pure-batch anchors.
  withr::local_seed(42L)
  p <- 3L
  b <- c(1, 0, 0)
  ref_data <- matrix(rnorm(40L * p, 0, 0.5), ncol = p,
                     dimnames = list(NULL, paste0("F", seq_len(p))))
  batch_mag <- 2.0
  qry_data <- ref_data + matrix(batch_mag * b, nrow(ref_data), p, byrow = TRUE)
  anc_idx <- seq(11L, 30L)
  anc_old <- ref_data[anc_idx, , drop = FALSE]
  anc_new <- anc_old + matrix(batch_mag * b, length(anc_idx), p, byrow = TRUE)
  ref <- somalign_train_reference(ref_data, grid = kohonen::somgrid(2L, 2L, "hexagonal"), rlen = 10L)
  qry <- somalign_query(qry_data, ref, grid = kohonen::somgrid(2L, 2L, "hexagonal"), rlen = 10L)

  fit <- somalign_fit_anchored(
    qry, ref, anchor_old = anc_old, anchor_new = anc_new,
    rho_anchor = 1, correction = "subspace", laplacian_lambda = 0.2
  )
  expect_true(!is.null(attr(fit$node_shifts, "correction_allowed")))
  V <- fit$anchors$batch_subspace$V
  shifts_orth <- fit$node_shifts - fit$node_shifts %*% V %*% t(V)
  expect_lt(sqrt(mean(shifts_orth^2)), 1e-8)
})
