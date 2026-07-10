test_that("internal OT solver returns finite coupling and diagnostics", {
  ref <- tiny_reference()
  query <- matrix(c(-1.1, 0, -0.9, 0, 0.95, 0, 1.2, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  fit <- somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1, rho_query = 1, rho_ref = 1)

  expect_s3_class(fit, "somalign_fit")
  expect_equal(dim(fit$transport_plan), c(2L, 3L))
  expect_true(all(is.finite(fit$transport_plan)))
  expect_true(all(fit$transport_plan >= 0))
  expect_equal(dim(fit$correspondence), c(2L, 3L))
  expect_named(somalign_diagnostics(fit), c("solver", "ot", "nodes", "projection"))
})

test_that("solver selection is pure R and keeps auto as a compatibility alias", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  fit_default <- expect_warning(
    somalign_fit(query_obj, ref, epsilon = 0.1),
    NA
  )
  expect_equal(somalign_diagnostics(fit_default)$solver$requested, "internal")
  expect_equal(somalign_diagnostics(fit_default)$solver$used, "internal")

  fit_auto <- expect_warning(
    somalign_fit(query_obj, ref, solver = "auto", epsilon = 0.1),
    NA
  )
  expect_equal(somalign_diagnostics(fit_auto)$solver$requested, "auto")
  expect_equal(somalign_diagnostics(fit_auto)$solver$used, "internal")

  expect_error(
    somalign_fit(query_obj, ref, solver = "pot", epsilon = 0.1),
    "internal.*auto|auto.*internal"
  )
})

test_that("OT parameters are validated before solver dispatch", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  expect_error(somalign_fit(query_obj, ref, epsilon = 0), "epsilon")
  expect_error(somalign_fit(query_obj, ref, rho_query = NA_real_), "rho_query")
  expect_error(somalign_fit(query_obj, ref, rho_ref = -1), "rho_ref")
})

test_that("label transfer exposes confidence, second best label, and low-match gating", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(0, 0), c(1, 0)))
  )

  fit <- somalign_fit(
    query_obj,
    ref,
    solver = "internal",
    epsilon = 0.05,
    min_match_fraction = 0.99,
    confidence_threshold = 0.8
  )

  expect_true(all(c("label", "confidence", "second_label", "entropy", "accepted") %in% names(fit$label_transfer)))
  expect_true(any(!fit$label_transfer$accepted))
  expect_true(any(is.na(fit$label_transfer$label[!fit$label_transfer$accepted])))
})

test_that("results keep direct projection as primary and corrected projection auxiliary", {
  ref <- tiny_reference()
  query <- matrix(c(-1.1, 0, 1.5, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1)
  results <- somalign_results(fit)

  expect_true(all(c(
    "old_som_unit",
    "old_som_distance",
    "outside_reference_distance",
    "final_status",
    "corrected_som_unit",
    "corrected_som_distance",
    "corrected_outside_reference_distance",
    "correction_norm"
  ) %in% names(results)))
  expect_equal(
    results$final_status,
    ifelse(results$outside_reference_distance, "outside_reference", "inside_reference")
  )
  expect_true(all(is.finite(results$correction_norm)))
})

test_that("weakly matched nodes receive zero correction", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  fit <- somalign_fit(
    query_obj,
    ref,
    solver = "internal",
    epsilon = 0.1,
    min_match_fraction = 2
  )

  expect_true(all(fit$diagnostics$nodes$correction_allowed == FALSE))
  expect_true(all(fit$projection$correction_norm == 0))
})
