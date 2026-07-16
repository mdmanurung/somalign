## Regression tests for cross-solver agreement (internal vs log_domain).
##
## Background: zero-mass nodes produce f_i = -Inf (structurally correct --
## exp(-Inf) = 0 in the plan row) but abs(-Inf - -Inf) = NaN, which previously
## poisoned the convergence delta and caused log_domain to emit a spurious
## non-finite-delta warning even though both solvers reached the same transport
## plan. The fix restricts the delta computation to finite ratios only, so that
## zero-mass nodes are excluded from the convergence check without masking
## genuine divergence.

test_that("log_domain: zero-mass query node does not trigger non-finite-delta warning", {
  # Regression: this previously emitted
  # "Sinkhorn solver produced a non-finite iterate delta after N iterations."
  cost <- matrix(c(0.0, 1.0,
                   1.0, 0.0,
                   0.5, 0.5), nrow = 3)
  a <- c(0.6, 0.0, 0.4)   # second node has exact zero mass
  b <- c(0.5, 0.5)
  expect_no_warning(
    somalign:::.somalign_solve_internal(
      cost, a, b, epsilon = 0.5, rho_query = 1, rho_ref = 1,
      max_iter = 1000, tol = 1e-9, log_domain = TRUE
    )
  )
})

test_that("log_domain: zero-mass query node reports converged = TRUE", {
  cost <- matrix(c(0.0, 1.0,
                   1.0, 0.0,
                   0.5, 0.5), nrow = 3)
  a <- c(0.6, 0.0, 0.4)
  b <- c(0.5, 0.5)
  result <- somalign:::.somalign_solve_internal(
    cost, a, b, epsilon = 0.5, rho_query = 1, rho_ref = 1,
    max_iter = 1000, tol = 1e-9, log_domain = TRUE
  )
  expect_true(result$converged,
              info = "zero-mass node produces f=-Inf which is load-bearing, not a divergence")
  expect_true(is.finite(result$final_delta))
})

test_that("log_domain: zero-mass reference node does not trigger non-finite-delta warning", {
  # Mirror test: g_j = -Inf when b[j] = 0 (reference side zero mass)
  cost <- matrix(c(0, 1, 1, 0), nrow = 2)
  a <- c(0.5, 0.5)
  b <- c(0.0, 1.0)   # first reference node has zero mass
  expect_no_warning(
    somalign:::.somalign_solve_internal(
      cost, a, b, epsilon = 0.5, rho_query = 1, rho_ref = 1,
      max_iter = 1000, tol = 1e-9, log_domain = TRUE
    )
  )
})

test_that("log_domain and internal solvers agree on mixed zero-mass scenario", {
  # Both zero-mass query nodes and zero-mass reference nodes present.
  # Plans must agree to machine-epsilon tolerance regardless of which solver
  # reports delta convergence.
  withr::local_seed(123L)
  n <- 30L
  cost <- matrix(abs(rnorm(n * n)), nrow = n, ncol = n)
  a <- abs(rnorm(n)); a[c(3, 11, 22)] <- 0; a <- a / sum(a)
  b <- abs(rnorm(n)); b[c(5, 18)]     <- 0; b <- b / sum(b)

  result_int <- somalign:::.somalign_solve_internal(
    cost, a, b, epsilon = 0.3, rho_query = 2, rho_ref = 2,
    max_iter = 2000, tol = 1e-8, log_domain = FALSE
  )
  result_log <- somalign:::.somalign_solve_internal(
    cost, a, b, epsilon = 0.3, rho_query = 2, rho_ref = 2,
    max_iter = 2000, tol = 1e-8, log_domain = TRUE
  )

  expect_true(result_int$converged, info = "internal solver must converge")
  expect_true(result_log$converged, info = "log_domain solver must converge with zero-mass nodes present")

  rel_frob <- norm(result_log$plan - result_int$plan, "F") /
              max(norm(result_int$plan, "F"), .Machine$double.eps)
  expect_lt(rel_frob, 1e-5,
            label = "relative Frobenius diff between log_domain and internal plans")
})

test_that("log_domain: all-zero-mass input still warns (genuinely degenerate)", {
  # When ALL query masses are zero the iteration enters a NaN reset cycle and
  # never converges. The solver must still emit a "non-finite iterate delta"
  # warning (now delta = Inf because finite_ratios is always empty in the limit
  # cycle) and must NOT set converged = TRUE.
  # .somalign_solve_internal does not call the validator, so only the solver's
  # own "non-finite" warning is present here (no need to suppress anything).
  cost <- matrix(c(0, 1, 1, 0), nrow = 2)
  a <- c(0, 0)   # genuinely degenerate -- all masses zero
  b <- c(0.5, 0.5)
  result <- NULL
  expect_warning(
    result <- somalign:::.somalign_solve_internal(
      cost, a, b, epsilon = 0.5, rho_query = 1, rho_ref = 1,
      max_iter = 20, tol = 1e-9, log_domain = TRUE
    ),
    "non-finite"
  )
  expect_false(result$converged,
               info = "all-zero-mass is degenerate: solver must not report convergence")
})
