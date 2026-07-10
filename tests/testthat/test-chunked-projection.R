test_that(".somalign_nearest_code_chunked matches unchunked for various chunk_sizes", {
  set.seed(42L)
  x <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  codebook <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)

  full          <- somalign:::.somalign_nearest_code(x, codebook)
  chunked3      <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = 3L)
  chunked_inf   <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = Inf)
  chunked_null  <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = NULL)
  chunked_big   <- somalign:::.somalign_nearest_code_chunked(x, codebook, chunk_size = 10000L)

  expect_identical(chunked3$unit,     full$unit)
  expect_equal(   chunked3$distance,  full$distance)
  expect_identical(chunked_inf$unit,  full$unit)
  expect_identical(chunked_null$unit, full$unit)
  expect_identical(chunked_big$unit,  full$unit)
})

test_that(".somalign_distance_quantiles vectorised output matches loop output on data with empty nodes", {
  set.seed(7L)
  n_nodes   <- 5L
  distances <- abs(rnorm(20))
  units <- sample(1:4, 20, replace = TRUE)
  probs <- c(0.5, 0.9, 0.95, 0.99)

  res <- somalign:::.somalign_distance_quantiles(distances, units, n_nodes, probs)
  expect_equal(nrow(res$node), n_nodes)
  expect_equal(ncol(res$node), length(probs))
  expect_equal(res$node[5, ], res$global)
  expect_null(rownames(res$node))
})

test_that(".somalign_label_probabilities vectorised output matches loop output", {
  set.seed(3L)
  n_nodes <- 4L
  labels  <- c("A", "B", NA,  "A", "B", "C", "A", NA,  "C", "B")
  units   <- c(1L,  1L,  2L,  2L,  3L,  3L,  4L,  4L,  1L,  2L)

  res <- somalign:::.somalign_label_probabilities(labels, units, n_nodes)
  expect_equal(nrow(res), n_nodes)
  expect_equal(colnames(res), c("A", "B", "C"))
  expect_equal(rowSums(res), rep(1, n_nodes))
  # Node 2: positions with unit==2 are indices 3 (NA, excluded), 4 (A), 10 (B)
  # -> A=1, B=1, so each is 0.5
  expect_equal(unname(res[2, "A"]), 0.5)
  expect_equal(unname(res[2, "B"]), 0.5)
  expect_equal(unname(res[2, "C"]), 0)
})
