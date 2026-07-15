# Helper: build a minimal reference + a matching query codebook
make_check_fixtures <- function(seed = 1L, nodes = 9L, features = c("F1","F2","F3")) {
  withr::local_seed(seed)
  p <- length(features)
  codebook <- matrix(rnorm(nodes * p), nrow = nodes, ncol = p,
                     dimnames = list(NULL, features))
  masses <- rep(1 / nodes, nodes)
  ref <- somalign_reference_from_nodes(
    codebook   = codebook,
    features   = features,
    center     = setNames(rep(0, p), features),
    scale      = setNames(rep(1, p), features),
    node_masses = masses
  )
  list(ref = ref, codebook = codebook, masses = masses, features = features)
}

# ---------------------------------------------------------------------------
# Pass case
# ---------------------------------------------------------------------------

test_that("identical codebooks give verdict=pass", {
  fx  <- make_check_fixtures()
  chk <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  expect_equal(chk$verdict, "pass")
  expect_equal(chk$n_critical_features, 0L)
  expect_equal(chk$n_warning_features,  0L)
  expect_true(all(chk$per_feature$flag == "ok"))
})

# ---------------------------------------------------------------------------
# Critical: zero range overlap
# ---------------------------------------------------------------------------

test_that("query shifted far outside reference gives verdict=critical", {
  fx <- make_check_fixtures()
  # Shift query codebook completely outside reference range on all features
  q_cb <- fx$codebook + 100
  expect_error(
    somalign_check_codebook_alignment(q_cb, fx$ref, stop_if_critical = TRUE),
    "zero range overlap"
  )
})

test_that("stop_if_critical=FALSE demotes critical to warning", {
  fx   <- make_check_fixtures()
  q_cb <- fx$codebook + 100
  expect_warning(
    chk <- somalign_check_codebook_alignment(q_cb, fx$ref, stop_if_critical = FALSE),
    "zero range overlap"
  )
  expect_equal(chk$verdict, "critical")
  expect_gt(chk$n_critical_features, 0L)
})

# ---------------------------------------------------------------------------
# Warning: partial mismatch on one feature
# ---------------------------------------------------------------------------

test_that("query shifted partially outside reference gives verdict=warning", {
  fx   <- make_check_fixtures()
  q_cb <- fx$codebook
  # Shift one feature so its range barely overlaps
  ref_range <- range(fx$codebook[, 1L])
  q_cb[, 1L] <- q_cb[, 1L] + (ref_range[2L] - ref_range[1L]) * 0.9
  chk <- suppressMessages(
    somalign_check_codebook_alignment(q_cb, fx$ref, stop_if_critical = FALSE)
  )
  expect_equal(chk$verdict, "warning")
})

# ---------------------------------------------------------------------------
# Return structure
# ---------------------------------------------------------------------------

test_that("result has expected fields and class", {
  fx  <- make_check_fixtures()
  chk <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  expect_s3_class(chk, "somalign_codebook_check")
  expect_true(all(c("per_feature", "cost_summary", "n_critical_features",
                    "n_warning_features", "verdict") %in% names(chk)))
  expect_equal(nrow(chk$per_feature), length(fx$features))
  expect_equal(chk$per_feature$feature, fx$features)
  expect_true(all(c("median_cost", "p95_cost", "cost_scale",
                    "fraction_near_eps") %in% names(chk$cost_summary)))
})

test_that("per_feature overlap_fraction is in [0, 1]", {
  fx  <- make_check_fixtures()
  chk <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  expect_true(all(chk$per_feature$overlap_fraction >= 0))
  expect_true(all(chk$per_feature$overlap_fraction <= 1))
})

test_that("cost_summary values are non-negative and finite", {
  fx  <- make_check_fixtures()
  chk <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  expect_true(all(is.finite(chk$cost_summary)))
  expect_true(all(chk$cost_summary >= 0))
})

# ---------------------------------------------------------------------------
# query_masses argument
# ---------------------------------------------------------------------------

test_that("query_masses length mismatch raises an error", {
  fx <- make_check_fixtures()
  expect_error(
    somalign_check_codebook_alignment(fx$codebook, fx$ref,
                                      query_masses = c(1, 2)),
    "query_masses.*length"
  )
})

test_that("supplying uniform query_masses gives same verdict as NULL", {
  fx  <- make_check_fixtures()
  n   <- nrow(fx$codebook)
  chk_null   <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  chk_masses <- somalign_check_codebook_alignment(fx$codebook, fx$ref,
                                                  query_masses = rep(1, n))
  expect_equal(chk_null$verdict, chk_masses$verdict)
  expect_equal(chk_null$per_feature$flag, chk_masses$per_feature$flag)
})

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

test_that("print.somalign_codebook_check runs without error", {
  fx  <- make_check_fixtures()
  chk <- somalign_check_codebook_alignment(fx$codebook, fx$ref)
  expect_output(print(chk), "verdict")
  expect_output(print(chk), "Median pairwise")
})

# ---------------------------------------------------------------------------
# Missing features raise clear error
# ---------------------------------------------------------------------------

test_that("query_codebook missing reference features raises an error", {
  fx   <- make_check_fixtures()
  q_cb <- fx$codebook[, 1L, drop = FALSE]  # only one feature
  expect_error(
    somalign_check_codebook_alignment(q_cb, fx$ref),
    "missing"
  )
})
