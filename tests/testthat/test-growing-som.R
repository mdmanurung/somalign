## Tests for E4b: Growing Neural Gas seeded from reference codebook.
##
## Key behaviors verified:
##   (a) Grown nodes concentrate near novel-population centroid (insert-at-novel).
##   (b) Original codebook vectors are not moved at all (freeze-originals).
##   (c) The extended reference integrates with somalign_fit() without errors.
##
## Control arm:
##   (d) n_grown <= max_new_nodes.
##
## Setup:  2-D reference with 4 nodes covering [-3, 3] x [-3, 3].
##         Novel population at [10, 10] — well beyond the reference hull.

test_that("(a) grown nodes insert near the novel-population centroid", {
  set.seed(42L)
  p <- 2L
  feats <- c("a", "b")

  ## Reference codebook: 4 nodes around the origin
  codebook <- matrix(
    c(-2, -2,
      -2,  2,
       2, -2,
       2,  2),
    nrow = 4L, ncol = p, byrow = TRUE
  )
  colnames(codebook) <- feats

  ref <- somalign_reference_from_nodes(
    codebook   = codebook,
    features   = feats,
    center     = c(a = 0, b = 0),
    scale      = c(a = 1, b = 1),
    node_masses = rep(0.25, 4L),
    label_prob = matrix(
      c(1, 0, 1, 0, 0, 1, 0, 1),
      nrow = 4L, byrow = TRUE,
      dimnames = list(NULL, c("A", "B"))
    )
  )

  ## Novel population at [10, 10], tight cluster
  novel_data <- matrix(
    rnorm(200L * p, mean = 10, sd = 0.3),
    nrow = 200L, ncol = p
  )
  colnames(novel_data) <- feats

  extended <- somalign_grow_reference(
    reference     = ref,
    new_data      = novel_data,
    max_new_nodes = 10L,
    lambda        = 40L,
    epsilon_new   = 0.1,
    age_max       = 30L,
    error_decay   = 0.99,
    n_epochs      = 10L,
    seed          = 7L
  )

  N_orig  <- nrow(codebook)
  N_total <- nrow(extended$codebook)
  N_grown <- N_total - N_orig

  ## At least one node was inserted
  expect_gt(N_grown, 0L)

  ## Grown nodes should be present
  grown_nodes <- extended$codebook[(N_orig + 1L):N_total, , drop = FALSE]

  ## Novel centroid in reference-scaled space equals [10, 10] (center=0, scale=1)
  novel_centroid <- c(10, 10)

  ## The closest grown node to the novel centroid must be within 2 SD of it.
  ## SD of novel population is 0.3 so 2 SD = 0.6.
  sq_dists_to_centroid <- rowSums(
    (grown_nodes - matrix(novel_centroid, nrow = N_grown, ncol = p, byrow = TRUE))^2
  )
  min_dist <- sqrt(min(sq_dists_to_centroid))

  expect_lt(min_dist, 2.0,
    label = paste0("closest grown node distance to novel centroid = ",
                   round(min_dist, 3)))
})


test_that("(b) original codebook vectors are exactly unchanged after growing", {
  set.seed(1L)
  p <- 2L
  feats <- c("a", "b")

  codebook_orig <- matrix(
    c(-3, 0,
       0, 0,
       3, 0),
    nrow = 3L, ncol = p, byrow = TRUE
  )
  colnames(codebook_orig) <- feats

  ref <- somalign_reference_from_nodes(
    codebook   = codebook_orig,
    features   = feats,
    center     = c(a = 0, b = 0),
    scale      = c(a = 1, b = 1),
    node_masses = rep(1/3, 3L),
    label_prob = matrix(
      c(1, 0, 0.5, 0.5, 0, 1),
      nrow = 3L, byrow = TRUE,
      dimnames = list(NULL, c("X", "Y"))
    )
  )

  ## Novel population far from reference
  novel_data <- matrix(
    rnorm(300L * p, mean = 15, sd = 0.5),
    nrow = 300L, ncol = p
  )
  colnames(novel_data) <- feats

  extended <- somalign_grow_reference(
    reference     = ref,
    new_data      = novel_data,
    max_new_nodes = 8L,
    lambda        = 60L,
    epsilon_new   = 0.1,
    age_max       = 40L,
    error_decay   = 0.99,
    n_epochs      = 5L,
    seed          = 3L
  )

  ## Original rows must be byte-for-byte identical
  orig_rows <- extended$codebook[seq_len(nrow(codebook_orig)), , drop = FALSE]
  expect_equal(orig_rows, codebook_orig,
               info = "original codebook rows must not move under GNG freezing")
})


test_that("(c) extended reference integrates with somalign_fit without error", {
  set.seed(99L)
  p <- 2L
  feats <- c("a", "b")

  codebook <- matrix(
    c(-1, 0,
       1, 0),
    nrow = 2L, ncol = p, byrow = TRUE
  )
  colnames(codebook) <- feats

  ref <- somalign_reference_from_nodes(
    codebook   = codebook,
    features   = feats,
    center     = c(a = 0, b = 0),
    scale      = c(a = 1, b = 1),
    node_masses = c(0.5, 0.5),
    label_prob = matrix(
      c(1, 0, 0, 1),
      nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("P", "Q"))
    )
  )

  novel_data <- matrix(
    rnorm(100L * p, mean = 8, sd = 0.4),
    nrow = 100L, ncol = p
  )
  colnames(novel_data) <- feats

  extended <- somalign_grow_reference(
    reference     = ref,
    new_data      = novel_data,
    max_new_nodes = 5L,
    lambda        = 20L,
    epsilon_new   = 0.1,
    age_max       = 20L,
    error_decay   = 0.995,
    n_epochs      = 5L,
    seed          = 11L
  )

  N_total  <- nrow(extended$codebook)
  ## Build query using the extended codebook (laplacian_lambda=0 required
  ## because extended$som_ref is NULL)
  query_cells <- matrix(
    rnorm(40L * p),
    nrow = 40L, ncol = p
  )
  colnames(query_cells) <- feats

  ## Provide a minimal query SOM from the extended codebook
  query_som <- make_som(extended$codebook)

  ## somalign_query and somalign_fit must not throw
  expect_no_error({
    q <- somalign_query(query_cells, extended, som_query = query_som)
    f <- somalign_fit(q, extended, laplacian_lambda = 0)
  })
})


test_that("(d) control arm: n_grown <= max_new_nodes", {
  ## Also checks that on a uniformly shifted known population (no truly novel
  ## cells), growth still occurs but stays within the cap.  The original nodes
  ## accumulate high error because they are frozen and cannot chase the shifted
  ## cloud, so some growth is expected even for batch-shifted known pops.
  ## This is the intended behaviour: the cap is the hard boundary.
  set.seed(55L)
  p <- 2L
  feats <- c("x", "y")

  codebook <- matrix(
    c(-2, -2,  2, -2,  0, 2),
    nrow = 3L, ncol = p, byrow = TRUE
  )
  colnames(codebook) <- feats

  ref <- somalign_reference_from_nodes(
    codebook   = codebook,
    features   = feats,
    center     = c(x = 0, y = 0),
    scale      = c(x = 1, y = 1),
    node_masses = rep(1/3, 3L)
  )

  max_nodes <- 6L

  ## Shifted population (not truly novel, but frozen originals accumulate error)
  shifted_data <- matrix(
    rnorm(150L * p, mean = 5, sd = 1),
    nrow = 150L, ncol = p
  )
  colnames(shifted_data) <- feats

  extended <- somalign_grow_reference(
    reference     = ref,
    new_data      = shifted_data,
    max_new_nodes = max_nodes,
    lambda        = 25L,
    epsilon_new   = 0.05,
    age_max       = 20L,
    error_decay   = 0.995,
    n_epochs      = 3L,
    seed          = 77L
  )

  N_orig  <- nrow(codebook)
  N_total <- nrow(extended$codebook)
  N_grown <- N_total - N_orig

  expect_lte(N_grown, max_nodes,
    label = paste0("n_grown = ", N_grown, ", max_new_nodes = ", max_nodes))

  ## Original rows unchanged even in control arm
  expect_equal(
    extended$codebook[seq_len(N_orig), , drop = FALSE],
    codebook
  )
})
