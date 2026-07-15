test_that("internal solver warns on non-convergence", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  expect_warning(
    somalign_fit(
      query_obj, ref,
      solver = "internal",
      epsilon = 0.1,
      max_iter = 1L,
      tol = .Machine$double.eps
    ),
    "did not converge"
  )
})

test_that("internal solver does NOT warn when it converges on last iteration", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  expect_no_warning(
    somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1)
  )
})

test_that("internal solver warns when kernel underflows", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  expect_warning(
    somalign_fit(
      query_obj, ref,
      solver = "internal",
      epsilon = 1e-300
    ),
    "underflowed"
  )
})
