test_that("somalign_plot_mass_balance returns a ggplot for a valid fit", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  p   <- somalign_plot_mass_balance(fit)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "Query node mass")
  expect_identical(p$labels$y, "Transported mass")
  # one point per node
  d <- ggplot2::layer_data(p, which(vapply(p$layers, function(l)
    inherits(l$geom, "GeomPoint"), logical(1L))))
  n_nodes <- nrow(fx$qry$codebook)
  expect_equal(nrow(d), n_nodes)
})

test_that("somalign_plot_mass_balance rejects non-fit input", {
  expect_error(somalign_plot_mass_balance(list()), "somalign_fit", fixed = TRUE)
})

test_that("somalign_plot_match_fraction returns ggplot with node bars", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  p   <- somalign_plot_match_fraction(fit)
  expect_s3_class(p, "ggplot")
  bar_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), logical(1L))
  expect_true(any(bar_layers))
  d <- ggplot2::layer_data(p, which(bar_layers)[1L])
  expect_equal(nrow(d), nrow(fx$qry$codebook))
})

test_that("somalign_plot_correction returns ggplot coloured by correction_allowed", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  p   <- somalign_plot_correction(fit)
  expect_s3_class(p, "ggplot")
  bar_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), logical(1L))
  expect_true(any(bar_layers))
})

test_that("somalign_plot_outside_fraction returns ggplot with two bars", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  p   <- somalign_plot_outside_fraction(fit)
  expect_s3_class(p, "ggplot")
  bar_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomCol"), logical(1L))
  d <- ggplot2::layer_data(p, which(bar_layers)[1L])
  expect_equal(nrow(d), 2L)
})

test_that("somalign_plot_label_confusion returns ggplot when labels exist", {
  skip_if_not_installed("kohonen")
  withr::local_seed(1L)
  p_data <- 3L
  ref_data <- rbind(
    matrix(rnorm(30 * p_data, mean = -3, sd = 0.5), ncol = p_data),
    matrix(rnorm(30 * p_data, mean =  3, sd = 0.5), ncol = p_data)
  )
  colnames(ref_data) <- paste0("F", seq_len(p_data))
  labels <- rep(c("low", "high"), each = 30L)
  ref <- somalign_train_reference(
    ref_data, labels = labels,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  qry_data <- ref_data + 0.5
  qry <- somalign_query(qry_data, ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10)
  fit <- somalign_fit(qry, ref)
  p   <- somalign_plot_label_confusion(fit)
  expect_s3_class(p, "ggplot")
  tile_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomTile"), logical(1L))
  expect_true(any(tile_layers))
})

test_that("somalign_plot_label_confusion rejects out-of-range min_confidence", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  expect_error(
    somalign_plot_label_confusion(fit, min_confidence = 2),
    "must be in"
  )
  expect_error(
    somalign_plot_label_confusion(fit, min_confidence = -0.1),
    "must be in"
  )
})

test_that("somalign_plot_label_confusion errors with no accepted labels", {
  skip_if_not_installed("kohonen")
  withr::local_seed(1L)
  p_data <- 3L
  ref_data <- rbind(
    matrix(rnorm(30 * p_data, mean = -3, sd = 0.5), ncol = p_data),
    matrix(rnorm(30 * p_data, mean =  3, sd = 0.5), ncol = p_data)
  )
  colnames(ref_data) <- paste0("F", seq_len(p_data))
  labels <- rep(c("low", "high"), each = 30L)
  ref <- somalign_train_reference(
    ref_data, labels = labels,
    grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )
  qry <- somalign_query(ref_data + 0.5, ref,
                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10)
  fit <- somalign_fit(qry, ref)
  # min_confidence = 1 is valid but no cell achieves confidence exactly 1
  expect_error(
    somalign_plot_label_confusion(fit, min_confidence = 1),
    "No accepted transferred labels found"
  )
})

test_that("somalign_worst_nodes returns a data frame ordered by match_fraction", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  wn  <- somalign_worst_nodes(fit, n = 2L)
  expect_s3_class(wn, "data.frame")
  expect_lte(nrow(wn), 2L)
  expect_true("match_fraction" %in% names(wn))
  expect_true("top_ref_label" %in% names(wn))
  if (nrow(wn) == 2L)
    expect_lte(wn$match_fraction[1L], wn$match_fraction[2L])
})

test_that("somalign_plot_codebook_range returns ggplot with segment layers", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  chk <- somalign_check_codebook_alignment(fx$qry$codebook, fx$ref,
                                           stop_if_critical = FALSE)
  p   <- somalign_plot_codebook_range(chk)
  expect_s3_class(p, "ggplot")
  seg_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1L))
  expect_true(any(seg_layers))
  # two rows per feature (Reference + Query)
  d <- ggplot2::layer_data(p, which(seg_layers)[1L])
  expect_equal(nrow(d), 2L * ncol(fx$qry$codebook))
})

test_that("somalign_plot_codebook_range rejects non-check input", {
  expect_error(somalign_plot_codebook_range(list()), "somalign_codebook_check", fixed = TRUE)
})

test_that("somalign_plot_marker_distributions returns ggplot with density layer", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  p  <- somalign_plot_marker_distributions(fx$qry)
  expect_s3_class(p, "ggplot")
  dens_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomDensity"), logical(1L))
  expect_true(any(dens_layers))
})

test_that("somalign_plot_marker_distributions adds rug layer with reference", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  p  <- somalign_plot_marker_distributions(fx$qry, reference = fx$ref)
  expect_s3_class(p, "ggplot")
  rug_layers <- vapply(p$layers, function(l) inherits(l$geom, "GeomRug"), logical(1L))
  expect_true(any(rug_layers))
})

test_that("somalign_plot_marker_distributions adds second density with reference_data", {
  skip_if_not_installed("kohonen")
  fx   <- make_anchored_fixture()
  ref_scaled <- sweep(sweep(fx$ref_data, 2, fx$ref$center, "-"), 2, fx$ref$scale, "/")
  p    <- somalign_plot_marker_distributions(fx$qry, reference_data = ref_scaled)
  dens <- vapply(p$layers, function(l) inherits(l$geom, "GeomDensity"), logical(1L))
  expect_gte(sum(dens), 2L)
})

test_that("somalign_plot_marker_distributions respects features subset", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  p  <- somalign_plot_marker_distributions(fx$qry, features = "F1")
  expect_s3_class(p, "ggplot")
  d  <- ggplot2::layer_data(p, 1L)
  # only one panel's worth of data
  expect_true(all(as.character(p$data$feature) == "F1"))
})

test_that("somalign_plot_marker_distributions rejects unknown features", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  expect_error(
    somalign_plot_marker_distributions(fx$qry, features = "ZZZZ"),
    "Unknown features"
  )
})

test_that("somalign_plot_marker_distributions rejects non-query input", {
  expect_error(somalign_plot_marker_distributions(list()), "somalign_query", fixed = TRUE)
})

test_that("somalign_plot_match_fraction rejects bad threshold", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  expect_error(somalign_plot_match_fraction(fit, threshold = 1.5), "must be in")
  expect_error(somalign_plot_match_fraction(fit, threshold = -0.1), "must be in")
  expect_error(somalign_plot_match_fraction(fit, threshold = "high"), "single finite number")
  expect_error(somalign_plot_match_fraction(fit, threshold = c(0.05, 0.1)), "single finite number")
})

test_that("somalign_worst_nodes rejects bad n", {
  skip_if_not_installed("kohonen")
  fx  <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  expect_error(somalign_worst_nodes(fit, n = 0L), "positive integer")
  expect_error(somalign_worst_nodes(fit, n = -1L), "positive integer")
  expect_error(somalign_worst_nodes(fit, n = "all"), "positive integer")
})

test_that("somalign_plot_marker_distributions rejects bad downsample/seed/features", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture()
  expect_error(
    somalign_plot_marker_distributions(fx$qry, downsample = 0),
    "positive number"
  )
  expect_error(
    somalign_plot_marker_distributions(fx$qry, downsample = -10),
    "positive number"
  )
  expect_error(
    somalign_plot_marker_distributions(fx$qry, seed = NA_real_),
    "numeric scalar"
  )
  expect_error(
    somalign_plot_marker_distributions(fx$qry, features = 123),
    "character vector"
  )
})

test_that("somalign_plot_marker_distributions rejects reference_data with missing columns", {
  skip_if_not_installed("kohonen")
  fx       <- make_anchored_fixture()
  bad_data <- matrix(1:10, nrow = 5, dimnames = list(NULL, c("X1", "X2")))
  expect_error(
    somalign_plot_marker_distributions(fx$qry, reference_data = bad_data),
    "missing columns"
  )
})

test_that("somalign_plot_codebook_range rejects malformed per_feature", {
  chk <- structure(
    list(per_feature = data.frame(feature = "F1", ref_min = 0)),
    class = "somalign_codebook_check"
  )
  expect_error(somalign_plot_codebook_range(chk), "per_feature missing")
})

test_that(".somalign_downsample_rows preserves RNG state", {
  mat <- matrix(rnorm(100), nrow = 50)
  # State immediately before and after the call should be identical
  set.seed(42L)
  expected_next <- runif(1)
  set.seed(42L)
  .somalign_downsample_rows(mat, n = 10L, seed = 7L)
  actual_next <- runif(1)
  expect_equal(actual_next, expected_next)
})

test_that(".somalign_downsample_rows returns all rows when nrow <= n", {
  mat <- matrix(1:20, nrow = 5)
  out <- .somalign_downsample_rows(mat, n = 10L)
  expect_equal(nrow(out), 5L)
})
