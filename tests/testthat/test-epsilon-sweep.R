## ---------------------------------------------------------------------------
## Tests for somalign_epsilon_sweep() / somalign_select_epsilon() and the
## per-fit log_Z / mutual_information / transport_entropy diagnostics
## (Ideas #1 statistical-physics phase-transition + #3 information-theoretic
## mutual-information selector, unified into one sweep + selector).
## ---------------------------------------------------------------------------

local_sweep_fixture <- function(seed = 42L) {
  withr::local_seed(seed)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 4,
                dimnames = list(NULL, c("F1", "F2", "F3", "F4")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(4, 4, "hexagonal"),
                                  rlen = 10)
  qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(4, 4, "hexagonal"),
                        rlen = 10)
  list(ref = ref, qry = qry)
}

test_that("somalign_epsilon_sweep returns the correct class and columns", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 8))
  expect_s3_class(sw, "somalign_epsilon_sweep")
  expect_true(all(c("epsilon", "log_epsilon", "Phi", "susceptibility", "log_Z",
                     "mutual_information", "conditional_entropy_mean",
                     "expected_cost", "transport_mass", "iterations",
                     "converged") %in% names(sw$table)))
})

test_that("Phi is non-decreasing in epsilon and stays in (0, 1]", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 10))
  phi <- sw$table$Phi[is.finite(sw$table$Phi)]
  expect_true(all(diff(phi) >= -1e-6))
  expect_true(all(phi > 0 & phi <= 1))
})

test_that("epsilon_c and epsilon_rec are within the grid and related by 0.3x", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 10))
  expect_gte(sw$epsilon_c, min(sw$table$epsilon))
  expect_lte(sw$epsilon_c, max(sw$table$epsilon))
  expect_equal(sw$epsilon_rec, 0.3 * sw$epsilon_c)
})

test_that("log_Z is NA for internal solver and finite for log_domain", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw_int <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 5,
                                                    solver = "internal"))
  expect_true(all(is.na(sw_int$table$log_Z)))
  sw_log <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 5,
                                                    solver = "log_domain"))
  expect_true(all(is.finite(sw_log$table$log_Z)))
})

test_that("single-epsilon sweep messages and returns NA epsilon_c", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  expect_message(
    sw1 <- somalign_epsilon_sweep(fx$qry, fx$ref, epsilon_grid = 0.1),
    "fewer than 3 grid points"
  )
  expect_true(is.na(sw1$epsilon_c))
})

test_that("print and plot methods work on a somalign_epsilon_sweep", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw <- suppressWarnings(somalign_epsilon_sweep(fx$qry, fx$ref, n_grid = 6))
  expect_output(print(sw))
  p <- plot(sw)
  expect_s3_class(p, "ggplot")
})

test_that("somalign_fit diagnostics include log_Z, mutual_information, transport_entropy", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  fit <- somalign_fit(fx$qry, fx$ref, solver = "log_domain")
  diag <- somalign_diagnostics(fit)
  expect_true(is.finite(diag$solver$log_Z))
  expect_true(is.numeric(diag$ot$mutual_information))
  expect_true("transport_entropy" %in% names(diag$nodes))
  expect_equal(length(diag$nodes$transport_entropy), nrow(fx$qry$codebook))

  fit_int <- somalign_fit(fx$qry, fx$ref, solver = "internal")
  expect_true(is.na(somalign_diagnostics(fit_int)$solver$log_Z))
})

test_that("somalign_select_epsilon returns a value within the grid for all three methods", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  eps_grid <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
  for (m in c("critical", "elbow", "entropy_fraction")) {
    sel <- suppressWarnings(somalign_select_epsilon(fx$qry, fx$ref, epsilon = eps_grid,
                                                    method = m))
    expect_s3_class(sel, "somalign_epsilon_selection")
    expect_equal(sel$method, m)
    expect_true(nrow(sel$curve) == length(eps_grid))
    if (is.finite(sel$selected_epsilon)) {
      expect_gte(sel$selected_epsilon, min(eps_grid))
      expect_lte(sel$selected_epsilon, max(eps_grid))
    }
  }
})

test_that("somalign_epsilon_sweep accepts solver = 'annealing' and matches log_domain", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sw_ann <- somalign_epsilon_sweep(fx$qry, fx$ref, epsilon_grid = c(0.05, 0.1, 0.2),
                                   solver = "annealing")
  sw_log <- somalign_epsilon_sweep(fx$qry, fx$ref, epsilon_grid = c(0.05, 0.1, 0.2),
                                   solver = "log_domain")
  expect_equal(sw_ann$table$Phi, sw_log$table$Phi, tolerance = 1e-6)
  expect_true(all(sw_ann$table$converged))
})

test_that("somalign_select_epsilon accepts solver = 'annealing'", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  sel <- somalign_select_epsilon(fx$qry, fx$ref, epsilon = c(0.05, 0.1, 0.2),
                                 solver = "annealing")
  expect_s3_class(sel, "somalign_epsilon_selection")
})

test_that("somalign_sensitivity_grid includes a mutual_information column", {
  skip_if_not_installed("kohonen")
  fx <- local_sweep_fixture()
  grid <- somalign_sensitivity_grid(fx$qry, fx$ref, epsilon = c(0.1, 0.5),
                                    rho_query = 1, rho_ref = 1)
  expect_true("mutual_information" %in% names(grid))
})
