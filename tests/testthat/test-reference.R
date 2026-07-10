test_that("reference objects preserve codebooks, masses, labels, and thresholds", {
  old <- rbind(
    c(-1, 0),
    c(-0.8, 0.1),
    c(0.9, 0),
    c(1.1, -0.1)
  )
  colnames(old) <- c("a", "b")
  labels <- c("A", "A", "B", "B")
  ref <- somalign_reference(make_som(scale(old)), old, labels = labels, codebook_space = "reference_scaled")

  expect_s3_class(ref, "somalign_reference")
  expect_equal(ref$features, c("a", "b"))
  expect_equal(nrow(ref$codebook), 4)
  expect_equal(sum(ref$node_masses), 1)
  expect_true(all(c("A", "B") %in% colnames(ref$label_prob)))
  expect_true(all(c("50%", "90%", "95%", "99%") %in% colnames(ref$distance_quantiles)))
})

test_that("node-level references reproduce direct projection behavior", {
  ref <- tiny_reference()
  query <- matrix(c(-1.2, 0, 0.1, 0, 1.3, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features

  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref, solver = "internal")
  results <- somalign_results(fit)

  expect_equal(results$old_som_unit, c(1L, 2L, 3L))
  expect_equal(round(results$old_som_distance, 6), c(0.2, 0.1, 0.3))
})

test_that("node-level references without distance thresholds do not mark samples inside", {
  codebook <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  colnames(codebook) <- c("a", "b")
  ref <- somalign_reference_from_nodes(
    codebook = codebook,
    features = colnames(codebook),
    center = c(a = 0, b = 0),
    scale = c(a = 1, b = 1)
  )
  query <- matrix(c(10, 10), ncol = 2)
  colnames(query) <- colnames(codebook)
  query_obj <- somalign_query(query, ref, som_query = codebook)
  fit <- somalign_fit(query_obj, ref, solver = "internal")
  results <- somalign_results(fit)

  expect_true(is.na(results$outside_reference_distance))
  expect_equal(results$final_status, "unknown_reference_distance")
})
