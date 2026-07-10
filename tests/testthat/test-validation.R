test_that("reference validation rejects malformed feature matrices", {
  unnamed <- matrix(rnorm(12), ncol = 3)
  expect_error(
    somalign_reference(make_som(unnamed), unnamed),
    "column names"
  )

  duplicated <- unnamed
  colnames(duplicated) <- c("a", "a", "b")
  expect_error(
    somalign_reference(make_som(duplicated), duplicated),
    "Duplicated"
  )

  zero_scale <- duplicated[, c(1, 3)]
  colnames(zero_scale) <- c("a", "b")
  zero_scale[, "b"] <- 1
  expect_error(
    somalign_reference(make_som(zero_scale), zero_scale),
    "zero variance"
  )
})

test_that("query validation uses saved reference feature order and scaling", {
  old <- matrix(c(1, 2, 2, 4, 3, 6, 4, 8), ncol = 2)
  colnames(old) <- c("a", "b")
  reference <- somalign_reference(make_som(scale(old)), old, codebook_space = "reference_scaled")

  query <- old[, c("b", "a")]
  query_obj <- somalign_query(query, reference, som_query = make_som(scale(old)))

  expected <- sweep(
    sweep(old[, reference$features], 2, reference$center, "-"),
    2,
    reference$scale,
    "/"
  )
  expect_equal(query_obj$scaled_data, expected)

  expect_error(
    somalign_query(query[, "b", drop = FALSE], reference, som_query = make_som(scale(old))),
    "Missing"
  )

  query_bad <- query
  query_bad[1, 1] <- Inf
  expect_error(
    somalign_query(query_bad, reference, som_query = make_som(scale(old))),
    "finite"
  )
})

test_that("pretrained query codebooks must carry feature names", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  unnamed_codebook <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)

  expect_error(
    somalign_query(query, ref, som_query = unnamed_codebook),
    "column names"
  )
})

test_that("existing reference SOM codebook coordinate system must be explicit", {
  old <- matrix(c(-1, 0, -0.8, 0.1, 0.9, 0, 1.1, -0.1), ncol = 2, byrow = TRUE)
  colnames(old) <- c("a", "b")

  expect_error(
    somalign_reference(make_som(scale(old)), old),
    "codebook_space"
  )

  raw_ref <- somalign_reference(
    make_som(old),
    old,
    codebook_space = "raw"
  )
  scaled_ref <- somalign_reference(
    make_som(scale(old)),
    old,
    codebook_space = "reference_scaled"
  )

  expect_equal(raw_ref$codebook, scaled_ref$codebook, tolerance = 1e-12)
})
