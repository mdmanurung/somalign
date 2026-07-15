## ---------------------------------------------------------------------------
## Tests for solver = "annealing" (Idea #2: simulated-annealing Sinkhorn,
## geometric epsilon-cooling schedule with warm-started dual potentials)
## ---------------------------------------------------------------------------

test_that("annealing solver matches log_domain on an easy problem", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit_log <- somalign_fit(qry, ref, solver = "log_domain", epsilon = 0.2)
  fit_ann <- somalign_fit(qry, ref, solver = "annealing", epsilon = 0.2,
                          anneal_start = 5, anneal_stages = 5L)
  expect_equal(fit_ann$transport_plan, fit_log$transport_plan, tolerance = 1e-4)
})

test_that("anneal_stages = 1 is identical to a log_domain cold start", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit_log <- somalign_fit(qry, ref, solver = "log_domain", epsilon = 0.1)
  fit_ann <- somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1, anneal_stages = 1L)
  expect_equal(fit_ann$transport_plan, fit_log$transport_plan, tolerance = 1e-10)
})

test_that("annealing diagnostics expose schedule and per-stage info", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit <- somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1,
                      anneal_start = 10, anneal_stages = 4L)
  diag <- somalign_diagnostics(fit)
  expect_equal(diag$solver$used, "annealing")
  expect_equal(length(diag$solver$anneal_schedule), 4L)
  expect_equal(diag$solver$anneal_schedule[4], 0.1)
  expect_true(diag$solver$anneal_schedule[1] > diag$solver$anneal_schedule[4])
  expect_equal(length(diag$solver$anneal_stage_info), 4L)
  expect_true(is.finite(diag$solver$log_Z))
})

test_that("anneal_start < 1 is rejected", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))
  expect_error(
    somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1, anneal_start = 0.5),
    "`anneal_start` must be >= 1"
  )
})

test_that("anneal_factor >= 1 is rejected", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))
  expect_error(
    somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1, anneal_factor = 1.5),
    "anneal_factor"
  )
})

test_that("annealing is available on somalign_fit_anchored and somalign_fit_two_pass", {
  skip_if_not_installed("kohonen")
  fx <- make_anchored_fixture(seed = 11L)
  fit_anc <- somalign_fit_anchored(fx$qry, fx$ref,
                                   anchor_old = fx$anchor_old,
                                   anchor_new = fx$anchor_new,
                                   rho_anchor = 1, solver = "annealing",
                                   anneal_stages = 3L)
  expect_equal(somalign_diagnostics(fit_anc)$solver$used, "annealing")

  mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref2 <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  qry2 <- somalign_query(mat + 0.5, ref2, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  fit_tp <- somalign_fit_two_pass(qry2, ref2, solver = "annealing", anneal_stages = 3L)
  expect_equal(somalign_diagnostics(fit_tp)$solver$used, "annealing")
})
