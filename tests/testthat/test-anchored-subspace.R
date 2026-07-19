## ---------------------------------------------------------------------------
## Tests for CellANOVA-inspired signal-preserving subspace correction
## (somalign_fit_anchored correction = "subspace" / "both")
## ---------------------------------------------------------------------------

# The make_subspace_fixture() helper now lives in helper-fixtures.R so it is
# shared with test-correct-expression.R.

# ---------------------------------------------------------------------------
# Backward compatibility: default mode is "cost_bonus", same as original
# ---------------------------------------------------------------------------
test_that("correction = 'cost_bonus' (default) is backward compatible", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit_default <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
    rho_anchor = 1
  )
  fit_explicit <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
    rho_anchor = 1, correction = "cost_bonus"
  )
  expect_equal(fit_default$node_shifts, fit_explicit$node_shifts)
  expect_equal(fit_default$anchors$nodes_covered,
               fit_explicit$anchors$nodes_covered)
  expect_null(fit_default$anchors$batch_subspace)
  expect_equal(fit_default$anchors$correction, "cost_bonus")
})

# ---------------------------------------------------------------------------
# Structure: subspace mode returns expected anchors fields
# ---------------------------------------------------------------------------
test_that("correction = 'subspace' exposes batch_subspace with V, rank, variance_explained", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace", variance_threshold = 0.9
  )
  expect_equal(fit$anchors$correction, "subspace")
  bs <- fit$anchors$batch_subspace
  expect_false(is.null(bs))
  p <- length(fx$ref$features)
  expect_true(is.matrix(bs$V))
  expect_equal(nrow(bs$V), p)
  expect_gte(ncol(bs$V), 1L)
  expect_equal(ncol(bs$V), bs$rank)
  expect_gte(bs$variance_explained, 0.9)
})

# ---------------------------------------------------------------------------
# Algebraic correctness: node_shifts lie in V_batch
# ---------------------------------------------------------------------------
test_that("subspace node_shifts are within the batch subspace (zero orthogonal component)", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace"
  )
  V <- fit$anchors$batch_subspace$V
  S <- fit$node_shifts
  S_perp <- S - S %*% V %*% t(V)
  # Frobenius norm of orthogonal component should be negligible
  expect_lt(norm(S_perp, "F"), 1e-10)
})

# ---------------------------------------------------------------------------
# Signal preservation: subspace preserves biology, cost_bonus erases it
# ---------------------------------------------------------------------------
test_that("subspace correction preserves query-only biology orthogonal to batch", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  anc_old <- fx$anc_old; anc_new <- fx$anc_new

  fit_sub <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = anc_old, anchor_new = anc_new,
    rho_anchor = 1, correction = "subspace"
  )
  fit_bon <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = anc_old, anchor_new = anc_new,
    rho_anchor = 1, correction = "cost_bonus"
  )

  # For the subpop nodes: measure the cc-component of the node shift.
  # Subspace mode should have near-zero cc-component (restricted to b-direction).
  # cost_bonus applies a full-space shift that changes the cc-component.
  sub_nodes <- unique(fx$qry$sample_unit[fx$sub_idx])
  cc <- fx$cc

  cc_sub <- fit_sub$node_shifts[sub_nodes, , drop = FALSE] %*% cc
  cc_bon <- fit_bon$node_shifts[sub_nodes, , drop = FALSE] %*% cc

  expect_lt(max(abs(cc_sub)), max(abs(cc_bon)) + 0.5)
})

# ---------------------------------------------------------------------------
# "both" mode: cost bonus + subspace projection
# ---------------------------------------------------------------------------
test_that("correction = 'both' applies cost bonus AND batch subspace restriction", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "both"
  )
  expect_equal(fit$anchors$correction, "both")
  expect_false(is.null(fit$anchors$batch_subspace))
  expect_gt(fit$anchors$coverage_fraction, 0)
  # Node shifts should lie in V_batch
  V <- fit$anchors$batch_subspace$V
  S <- fit$node_shifts
  expect_lt(norm(S - S %*% V %*% t(V), "F"), 1e-10)
})

# ---------------------------------------------------------------------------
# Validation: bad arguments error cleanly
# ---------------------------------------------------------------------------
test_that("somalign_fit_anchored validates correction and variance_threshold", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                          anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
                          correction = "bad_mode"),
    "cost_bonus"
  )
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                          anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
                          variance_threshold = 0),
    "variance_threshold"
  )
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                          anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
                          variance_threshold = 1.1),
    "variance_threshold"
  )
  expect_error(
    somalign_fit_anchored(fx$qry, fx$ref,
                          anchor_old = fx$anchor_old, anchor_new = fx$anchor_new,
                          variance_threshold = NA_real_),
    "variance_threshold"
  )
})

# ---------------------------------------------------------------------------
# Two-pass: batch_subspace diagnostic present, correction math unchanged
# ---------------------------------------------------------------------------
test_that("somalign_fit_two_pass batch_subspace diagnostic is present and read-only", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit <- somalign_fit_two_pass(fx$qry, fx$ref, variance_threshold = 0.9)
  bs <- fit$two_pass$batch_subspace
  expect_false(is.null(bs))
  expect_true(is.matrix(bs$V))
  p <- length(fx$ref$features)
  expect_equal(nrow(bs$V), p)
  expect_gte(bs$variance_explained, 0.9)
})

test_that("two-pass batch_subspace diagnostic does not change correction", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  fit_09 <- somalign_fit_two_pass(fx$qry, fx$ref, variance_threshold = 0.9)
  fit_10 <- somalign_fit_two_pass(fx$qry, fx$ref, variance_threshold = 1.0)
  # Different variance_threshold → different batch_subspace, but same correction
  expect_equal(fit_09$node_shifts, fit_10$node_shifts)
  expect_equal(fit_09$projection, fit_10$projection)
})

test_that("two-pass validates variance_threshold", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  expect_error(
    somalign_fit_two_pass(fx$qry, fx$ref, variance_threshold = 0),
    "variance_threshold"
  )
  expect_error(
    somalign_fit_two_pass(fx$qry, fx$ref, variance_threshold = 2),
    "variance_threshold"
  )
})
