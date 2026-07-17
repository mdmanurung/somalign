## ---------------------------------------------------------------------------
## Tests for somalign_soft_labels() / somalign_soft_frequencies():
## probabilistic k-NN projection of query cells onto reference labels.
## ---------------------------------------------------------------------------

soft_fixture <- function(seed = 1L) {
  withr::local_seed(seed)
  x <- rbind(matrix(rnorm(150 * 3, -3, 0.5), ncol = 3),
             matrix(rnorm(150 * 3,  3, 0.5), ncol = 3))
  colnames(x) <- paste0("m", seq_len(3))
  lab <- rep(c("low", "high"), each = 150)
  grid <- kohonen::somgrid(3, 3, "hexagonal")
  ref <- somalign_train_reference(x, labels = lab, grid = grid, rlen = 15)
  qry <- somalign_query(x, ref, grid = grid, rlen = 15)
  list(fit = somalign_fit(qry, ref), x = x, lab = lab)
}

test_that("soft labels form a per-cell distribution over labels", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  soft <- somalign_soft_labels(fx$fit, bandwidth = 0.5)
  expect_s3_class(soft, "somalign_soft_labels")
  expect_equal(nrow(soft), nrow(fx$x))
  expect_setequal(colnames(soft), c("low", "high"))
  rs <- rowSums(soft)
  expect_true(all(abs(rs - 1) < 1e-8 | abs(rs) < 1e-8))   # sums to 1 (or 0 if unlabelled nbrs)
})

test_that("soft argmax recovers the true label for well-separated classes", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  soft <- somalign_soft_labels(fx$fit, bandwidth = 0.5)
  pred <- colnames(soft)[max.col(soft, ties.method = "first")]
  expect_gt(mean(pred == fx$lab), 0.95)
})

test_that("unlabelled reference errors unless node_groups supplied", {
  skip_if_not_installed("kohonen")
  withr::local_seed(2)
  x <- matrix(rnorm(120 * 3), ncol = 3, dimnames = list(NULL, paste0("m", 1:3)))
  grid <- kohonen::somgrid(3, 3, "hexagonal")
  ref <- somalign_train_reference(x, grid = grid, rlen = 10)   # no labels
  qry <- somalign_query(x, ref, grid = grid, rlen = 10)
  fit <- somalign_fit(qry, ref)
  expect_error(somalign_soft_labels(fit), regexp = "carries no labels")
  # a custom node grouping works even without reference labels
  ng <- rep(c("A", "B", "C"), length.out = nrow(fit$reference$codebook))
  soft <- somalign_soft_labels(fit, node_groups = ng, bandwidth = 0.5)
  expect_setequal(colnames(soft), c("A", "B", "C"))
  expect_equal(nrow(soft), nrow(x))
})

test_that("k is clamped to the number of reference nodes", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  m <- nrow(fx$fit$reference$codebook)
  soft <- somalign_soft_labels(fx$fit, k = 50L, bandwidth = 0.5)
  expect_equal(attr(soft, "k"), min(50L, m))
})

test_that("soft frequencies aggregate soft labels by group and normalise", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  grp <- rep(c("s1", "s2", "s3"), length.out = nrow(fx$x))
  freq <- somalign_soft_frequencies(fx$fit, grp, bandwidth = 0.5)
  expect_s3_class(freq, "somalign_soft_frequencies")
  expect_setequal(rownames(freq), c("s1", "s2", "s3"))
  expect_true(all(abs(rowSums(freq) - 1) < 1e-8))

  # matches a manual rowsum-normalise of the soft-label matrix
  soft <- somalign_soft_labels(fx$fit, bandwidth = 0.5)
  manual <- rowsum(unclass(soft), grp)
  manual <- manual / rowSums(manual)
  expect_equal(unclass(freq)[c("s1", "s2", "s3"), ], manual[c("s1", "s2", "s3"), ],
               tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("soft labels are invariant to chunk_size (fused chunking)", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  a <- somalign_soft_labels(fx$fit, bandwidth = 0.5, chunk_size = 1000000L)
  b <- somalign_soft_labels(fx$fit, bandwidth = 0.5, chunk_size = 7L)
  expect_equal(unclass(a), unclass(b), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("soft frequencies errors on a non-fit object", {
  expect_error(somalign_soft_frequencies(list(), group = 1),
               regexp = "must be a somalign_fit")
})

test_that("soft frequencies checks group length and supports raw counts", {
  skip_if_not_installed("kohonen")
  fx <- soft_fixture()
  expect_error(somalign_soft_frequencies(fx$fit, group = c("s1", "s2")),
               regexp = "one entry per query cell")
  grp <- rep(c("s1", "s2"), length.out = nrow(fx$x))
  counts <- somalign_soft_frequencies(fx$fit, grp, normalize = FALSE, bandwidth = 0.5)
  # un-normalised soft counts per group sum to that group's cell count
  expect_equal(as.numeric(rowSums(counts))[order(rownames(counts))],
               as.numeric(table(grp))[order(names(table(grp)))], tolerance = 1e-8)
})
