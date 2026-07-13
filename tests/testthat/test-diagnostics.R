## Tests for somalign_som_stability

test_that("somalign_som_stability returns a data frame with correct shape", {
  skip_if_not_installed("kohonen")
  withr::local_seed(42L)
  p <- 2L
  mat <- rbind(
    matrix(rnorm(20 * p, mean = -3), ncol = p),
    matrix(rnorm(20 * p, mean =  3), ncol = p)
  )
  colnames(mat) <- c("A", "B")
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 10
  )

  result <- somalign_som_stability(
    mat, ref,
    som_seeds = 1:3,
    grid = kohonen::somgrid(2, 2, "hexagonal"),
    rlen = 10
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3L)
  expect_true("som_seed" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("transport_mass" %in% names(result))
  expect_true(all(result$som_seed == 1:3))
})

test_that("somalign_som_stability restores caller RNG state", {
  skip_if_not_installed("kohonen")
  withr::local_seed(99L)
  p <- 2L
  mat <- matrix(rnorm(20 * p), ncol = p, dimnames = list(NULL, c("A", "B")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )

  set.seed(123L)
  before <- .Random.seed
  somalign_som_stability(mat, ref, som_seeds = 1:2,
                         grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_equal(.Random.seed, before)
})

test_that("somalign_som_stability rejects empty som_seeds", {
  skip_if_not_installed("kohonen")
  withr::local_seed(1L)
  mat <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("A", "B")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  expect_error(
    somalign_som_stability(mat, ref, som_seeds = integer(0)),
    "non-empty"
  )
})

test_that("somalign_som_stability does not leak .Random.seed when none existed before", {
  skip_if_not_installed("kohonen")
  withr::local_seed(1L)
  mat <- matrix(rnorm(20), ncol = 2, dimnames = list(NULL, c("A", "B")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  rm(".Random.seed", envir = .GlobalEnv)
  somalign_som_stability(mat, ref, som_seeds = 1:2,
                         grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
