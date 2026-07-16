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
  expect_named(somalign_diagnostics(fit),
               c("solver", "ot", "nodes", "projection", "cost_metric"))
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

test_that("somalign_results exposes second-best label and margin for triage", {
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
  results <- somalign_results(fit)

  expect_true(all(c(
    "transferred_label_second",
    "transferred_label_second_confidence",
    "transferred_label_margin"
  ) %in% names(results)))

  query_unit <- fit$query$sample_unit
  expect_equal(results$transferred_label_second, fit$label_transfer$second_label[query_unit])
  expect_equal(results$transferred_label_second_confidence, fit$label_transfer$second_confidence[query_unit])

  has_second <- !is.na(results$transferred_label_second_confidence)
  expect_equal(
    results$transferred_label_margin[has_second],
    fit$label_transfer$confidence[query_unit][has_second] -
      results$transferred_label_second_confidence[has_second]
  )
  no_second <- !has_second & !is.na(fit$label_transfer$confidence[query_unit])
  if (any(no_second)) {
    expect_equal(
      results$transferred_label_margin[no_second],
      fit$label_transfer$confidence[query_unit][no_second]
    )
  }
})

test_that("somalign_results errors when data has wrong row count", {
  ref <- tiny_reference()
  query <- matrix(c(-1.1, 0, 1.5, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref, solver = "internal", epsilon = 0.1)

  wrong <- data.frame(extra = seq_len(nrow(query) + 1L))
  expect_error(
    somalign_results(fit, data = wrong),
    "one row per query sample"
  )

  ok <- data.frame(extra = seq_len(nrow(query)))
  res <- somalign_results(fit, data = ok)
  expect_true("extra" %in% names(res))
})

test_that("nodes below mass threshold receive zero correction", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  query_obj <- somalign_query(
    query,
    ref,
    som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )

  # correction_min_mass far exceeds any node's transported mass (~0.5),
  # so all nodes must have correction_allowed = FALSE.
  fit <- somalign_fit(
    query_obj,
    ref,
    solver = "internal",
    epsilon = 0.1,
    correction_min_mass = 1e6
  )

  expect_true(all(fit$diagnostics$nodes$correction_allowed == FALSE))
  expect_true(all(fit$projection$correction_norm == 0))
})

## --------------------------------------------------------------------------
## Labels-first rescope (item 1): summary, print headline, include_correction
## --------------------------------------------------------------------------

test_that(".somalign_label_summary reports accepted fraction and class mix", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1.1, 0, 1.05, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, ref$features)),
    ref, som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref)
  s <- somalign:::.somalign_label_summary(fit)
  expect_true(s$enabled)
  expect_gte(s$accepted_fraction, 0)
  expect_lte(s$accepted_fraction, 1)
  expect_equal(s$n_cells, nrow(query_obj$scaled_data))
  expect_true(s$n_classes >= 1)
})

test_that(".somalign_label_summary reports disabled when reference has no labels", {
  fx <- make_anchored_fixture()
  fit <- somalign_fit(fx$qry, fx$ref)
  s <- somalign:::.somalign_label_summary(fit)
  expect_false(s$enabled)
})

test_that("summary.somalign_fit returns the fit invisibly", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1.1, 0, 1.05, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, ref$features)),
    ref, som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref)
  out <- withVisible(summary(fit))
  expect_false(out$visible)
  expect_identical(out$value, fit)
})

test_that("somalign_results(include_correction = FALSE) drops correction columns", {
  ref <- tiny_reference()
  query_obj <- somalign_query(
    matrix(c(-1.1, 0, 1.05, 0), ncol = 2, byrow = TRUE,
           dimnames = list(NULL, ref$features)),
    ref, som_query = make_som(rbind(c(-1, 0), c(1, 0)))
  )
  fit <- somalign_fit(query_obj, ref)
  full <- somalign_results(fit)
  lean <- somalign_results(fit, include_correction = FALSE)
  correction_cols <- c("corrected_som_unit", "corrected_som_distance",
                       "corrected_som_distance_threshold",
                       "corrected_outside_reference_distance", "correction_norm")
  expect_true(all(correction_cols %in% names(full)))
  expect_false(any(correction_cols %in% names(lean)))
  # label columns preserved
  expect_true(all(c("transferred_label", "transferred_label_confidence",
                    "transferred_label_margin") %in% names(lean)))
})
