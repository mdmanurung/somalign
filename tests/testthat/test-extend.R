test_that("somalign_extend_reference produces a valid reference with correct dimensions", {
  ref <- tiny_reference()  # 3 nodes, features c("a","b"), labels A and B
  n_orig <- nrow(ref$codebook)  # 3

  # New node clearly in novel region
  new_cb <- matrix(c(2, 0), nrow = 1, dimnames = list(NULL, c("a", "b")))

  extended <- somalign_extend_reference(ref, new_cb, new_labels = "C")

  # (a) fit succeeds — tested more directly below; check structure first
  expect_s3_class(extended, "somalign_reference")
  expect_equal(ref$features, extended$features)

  # All per-node arrays have n_orig + n_new rows
  n_new   <- 1L
  n_total <- n_orig + n_new
  expect_equal(nrow(extended$codebook),          n_total)
  expect_equal(length(extended$node_masses),     n_total)
  expect_equal(nrow(extended$label_prob),        n_total)
  expect_equal(nrow(extended$distance_quantiles), n_total)

  # Masses renormalise to 1
  expect_equal(sum(extended$node_masses), 1, tolerance = 1e-12)

  # Original relative proportions preserved
  orig_ratio <- ref$node_masses[1] / ref$node_masses[3]
  ext_ratio  <- extended$node_masses[1] / extended$node_masses[3]
  expect_equal(orig_ratio, ext_ratio, tolerance = 1e-9)

  # New class column introduced; old columns still present
  expect_true(all(c("A", "B", "C") %in% colnames(extended$label_prob)))

  # New node is one-hot for C
  expect_equal(unname(extended$label_prob[n_total, "C"]), 1.0, tolerance = 1e-9)

  # Old nodes have zero probability for the new class C
  expect_true(all(extended$label_prob[seq_len(n_orig), "C"] == 0))
})


test_that("OT routes mass to new node when query is near it", {
  ref <- tiny_reference()  # 3 nodes at (-1,0), (0,0), (1,0)

  # New node at (2, 0)
  new_cb <- matrix(c(2, 0), nrow = 1, dimnames = list(NULL, c("a", "b")))
  extended <- somalign_extend_reference(ref, new_cb, new_labels = "C")

  # Query cells: three near the original right node, three near the new node
  query_mat <- matrix(
    c(
      -1.0, 0, # near left node (orig #1)
       0.0, 0, # near middle node (orig #2)
       1.0, 0, # near right node (orig #3)
       2.0, 0, # at new node (#4)
       2.1, 0,
       1.9, 0
    ),
    ncol = 2, byrow = TRUE
  )
  colnames(query_mat) <- c("a", "b")

  query_som_cb <- rbind(
    ref$codebook,
    new_cb
  )
  query_obj <- somalign_query(
    query_mat,
    extended,
    som_query = make_som(query_som_cb)
  )

  # (a) fit succeeds
  fit <- somalign_fit(query_obj, extended, solver = "internal")
  expect_s3_class(fit, "somalign_fit")

  # (b) OT routes mass to new node (column 4)
  n_total    <- nrow(extended$codebook)
  col_masses <- fit$diagnostics$ot$col_mass
  expect_gt(col_masses[n_total], 0)

  # (c) transferred_label can equal the new label "C"
  results <- somalign_results(fit)
  expect_true("C" %in% results$transferred_label | "C" %in% results$old_som_label)

  # (d) somalign_results() / direct projection / acceptance gate run without error
  expect_s3_class(results, "data.frame")
  expect_true("transferred_label" %in% names(results))
  expect_true("final_status"      %in% names(results))
})


test_that("somalign_extend_reference with multiple new nodes and no new labels", {
  ref <- tiny_reference()
  n_orig <- nrow(ref$codebook)

  new_cb <- matrix(
    c(3, 0,
      4, 0),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )

  # No new_labels: new nodes should inherit uniform distribution over A/B
  extended <- somalign_extend_reference(ref, new_cb)

  n_new   <- 2L
  n_total <- n_orig + n_new

  expect_equal(nrow(extended$codebook),           n_total)
  expect_equal(length(extended$node_masses),      n_total)
  expect_equal(nrow(extended$label_prob),         n_total)
  expect_equal(nrow(extended$distance_quantiles), n_total)
  expect_equal(sum(extended$node_masses), 1, tolerance = 1e-12)

  # New nodes have equal probability for A and B (uniform fallback)
  new_rows <- extended$label_prob[seq(n_orig + 1L, n_total), , drop = FALSE]
  expect_equal(new_rows[, "A"], c(0.5, 0.5), tolerance = 1e-9)
  expect_equal(new_rows[, "B"], c(0.5, 0.5), tolerance = 1e-9)

  # node_var absent in tiny_reference -> still absent after extension
  expect_null(extended$node_var)
})


test_that("somalign_extend_reference with soft label matrix and caller-supplied masses", {
  ref <- tiny_reference()
  n_orig <- nrow(ref$codebook)

  new_cb <- matrix(c(5, 0), nrow = 1, dimnames = list(NULL, c("a", "b")))

  # Soft probability row: 80% A, 20% D (new class)
  soft_lp <- matrix(c(0.8, 0.2), nrow = 1,
                    dimnames = list(NULL, c("A", "D")))

  # Caller supplies explicit mass
  extended <- somalign_extend_reference(
    ref, new_cb,
    new_labels      = soft_lp,
    new_node_masses = 0.1
  )

  n_total <- n_orig + 1L
  expect_equal(nrow(extended$codebook), n_total)
  expect_equal(sum(extended$node_masses), 1, tolerance = 1e-12)

  # New class D is present
  expect_true("D" %in% colnames(extended$label_prob))

  # New node: row-normalised soft probs (0.8/1.0, 0, 0.2/1.0)
  last_row <- extended$label_prob[n_total, ]
  expect_equal(unname(last_row["A"]), 0.8, tolerance = 1e-9)
  expect_equal(unname(last_row["D"]), 0.2, tolerance = 1e-9)
  expect_equal(unname(last_row["B"]), 0.0, tolerance = 1e-9)
})


test_that("somalign_extend_reference with node_var propagates to extended object", {
  # Build a reference that has node_var
  codebook <- matrix(c(-1, 0, 0, 0, 1, 0), nrow = 3, byrow = TRUE)
  colnames(codebook) <- c("a", "b")
  node_var_mat <- matrix(
    c(0.1, 0.2,
      0.3, 0.4,
      0.5, 0.6),
    nrow = 3, byrow = TRUE,
    dimnames = list(NULL, c("a", "b"))
  )
  ref_with_var <- somalign_reference_from_nodes(
    codebook = codebook,
    features = c("a", "b"),
    center   = c(a = 0, b = 0),
    scale    = c(a = 1, b = 1),
    label_prob = matrix(c(1, 0, 0.5, 0.5, 0, 1), nrow = 3,
                        dimnames = list(NULL, c("A", "B"))),
    node_var = node_var_mat
  )

  new_cb <- matrix(c(2, 0), nrow = 1, dimnames = list(NULL, c("a", "b")))
  extended <- somalign_extend_reference(ref_with_var, new_cb, new_labels = "C")

  expect_false(is.null(extended$node_var))
  expect_equal(nrow(extended$node_var), 4L)

  # New node var is column-mean of original rows
  expected_var_a <- mean(node_var_mat[, "a"])  # (0.1+0.3+0.5)/3
  expected_var_b <- mean(node_var_mat[, "b"])  # (0.2+0.4+0.6)/3
  expect_equal(unname(extended$node_var[4L, "a"]), expected_var_a, tolerance = 1e-9)
  expect_equal(unname(extended$node_var[4L, "b"]), expected_var_b, tolerance = 1e-9)
})
