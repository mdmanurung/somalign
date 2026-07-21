# Ground truth: a query codebook that is the reference codebook, reordered by a
# known permutation and rigidly rotated into a *different* marker space. Rotation
# preserves all pairwise distances, so C1 = C2[perm, perm] exactly and the true
# Gromov-Wasserstein correspondence is unique (generic asymmetric points have no
# distance-preserving symmetry). GW must recover the permutation.
.gw_case <- function(m = 6L, p_dim = 3L, q_dim = 4L, seed = 123L) {
  withr::with_seed(seed, {
    ref_cb <- matrix(stats::rnorm(m * p_dim), m, p_dim)
    perm <- sample.int(m)
    # rotate the permuted reference into a q_dim marker space (>= p_dim) via a
    # random orthonormal embedding, so query "markers" differ from reference.
    E <- qr.Q(qr(matrix(stats::rnorm(p_dim * q_dim), p_dim, q_dim)))  # p x q, E E^T = I
    qry_cb <- ref_cb[perm, , drop = FALSE] %*% E                      # distances preserved
    list(ref_cb = ref_cb, qry_cb = qry_cb, perm = perm, m = m)
  })
}

test_that("entropic GW recovers a known node correspondence (distance-preserving map)", {
  cs <- .gw_case()
  C1 <- sqrt(somalign:::.somalign_pairwise_distance(cs$qry_cb, cs$qry_cb))
  C2 <- sqrt(somalign:::.somalign_pairwise_distance(cs$ref_cb, cs$ref_cb))
  # sanity: the intra-set structures match under the permutation
  expect_equal(unname(C1), unname(C2[cs$perm, cs$perm]), tolerance = 1e-8)

  p <- rep(1 / cs$m, cs$m)
  gw <- somalign:::.somalign_gromov_wasserstein(C1, C2, p, p, epsilon = 0.01,
                                                max_iter = 100L)
  recovered <- max.col(gw$coupling, ties.method = "first")
  expect_equal(recovered, cs$perm)                 # query node i -> ref node perm[i]
  expect_true(all(rowSums(gw$coupling) > 0))
})

test_that("somalign_fit_gw returns a coupling and transfers labels through it", {
  cs <- .gw_case(seed = 7L)
  labels <- paste0("c", seq_len(cs$m))             # one label per reference node
  lp <- diag(cs$m); colnames(lp) <- labels          # node i is purely label c i

  ref <- structure(list(codebook = cs$ref_cb, node_masses = rep(1 / cs$m, cs$m),
                        label_prob = lp),
                   class = "somalign_reference")
  qry <- structure(list(codebook = cs$qry_cb, node_masses = rep(1 / cs$m, cs$m)),
                   class = "somalign_query")

  fit <- somalign_fit_gw(qry, ref, epsilon = 0.01)
  expect_s3_class(fit, "somalign_gw_fit")
  expect_equal(dim(fit$coupling), c(cs$m, cs$m))
  expect_equal(rowSums(fit$correspondence), rep(1, cs$m), tolerance = 1e-6)
  # query node i should receive reference node perm[i]'s label
  expect_equal(fit$transferred_label, labels[cs$perm])

  out <- withVisible(print(fit))
  expect_false(out$visible)
})
