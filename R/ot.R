.somalign_solve_ot <- function(cost,
                               a,
                               b,
                               epsilon,
                               rho_query,
                               rho_ref,
                               solver,
                               max_iter,
                               tol) {
  requested_solver <- solver
  notes <- character()
  .somalign_validate_ot_inputs(cost, a, b, epsilon, rho_query, rho_ref)

  if (solver == "auto") {
    solver <- "internal"
    notes <- c(notes, "`solver = \"auto\"` uses the internal generalized Sinkhorn solver.")
  }
  internal <- .somalign_solve_internal(cost, a, b, epsilon, rho_query, rho_ref, max_iter, tol)
  plan <- internal$plan
  iterations <- internal$iterations

  list(
    plan = plan,
    solver = solver,
    requested_solver = requested_solver,
    notes = notes,
    iterations = iterations
  )
}

.somalign_solve_internal <- function(cost,
                                     a,
                                     b,
                                     epsilon,
                                     rho_query,
                                     rho_ref,
                                     max_iter,
                                     tol) {
  tiny <- .Machine$double.xmin
  k_raw <- exp(-cost / epsilon)
  underflow_fraction <- sum(k_raw < tiny) / length(k_raw)
  if (underflow_fraction > 0.01) {
    safe_eps <- signif(-max(cost) / log(.Machine$double.xmin), 3)
    warning(
      sprintf(
        "%.1f%% of Sinkhorn kernel entries underflowed (epsilon = %g). ",
        100 * underflow_fraction, epsilon
      ),
      sprintf(
        "Raise epsilon or reduce cost scale. Safe lower bound for epsilon: %g",
        safe_eps
      ),
      call. = FALSE
    )
  }
  k <- pmax(k_raw, tiny)
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_ref / (rho_ref + epsilon)

  u <- rep(1, length(a))
  v <- rep(1, length(b))
  delta <- Inf
  iterations <- max_iter
  for (iter in seq_len(max_iter)) {
    u_old <- u
    v_old <- v
    kv <- as.numeric(k %*% v)
    u <- (a / pmax(kv, tiny)) ^ tau_a
    ktu <- as.numeric(crossprod(k, u))
    v <- (b / pmax(ktu, tiny)) ^ tau_b
    u[!is.finite(u)] <- 0
    v[!is.finite(v)] <- 0
    delta <- max(
      abs(u - u_old) / pmax(1, abs(u_old)),
      abs(v - v_old) / pmax(1, abs(v_old))
    )
    if (is.finite(delta) && delta < tol) {
      iterations <- iter
      break
    }
  }

  if (iterations == max_iter && delta >= tol) {
    warning(
      sprintf(
        "Sinkhorn solver did not converge after %g iterations (final delta = %.3e). ",
        max_iter, delta
      ),
      "Consider increasing max_iter, raising epsilon, or reducing rho_query / rho_ref.",
      call. = FALSE
    )
  }

  plan <- sweep(sweep(k, 1, u, "*"), 2, v, "*")
  plan[!is.finite(plan)] <- 0
  plan <- pmax(plan, 0)
  list(plan = plan, iterations = iterations)
}

.somalign_validate_ot_inputs <- function(cost, a, b, epsilon, rho_query, rho_ref) {
  if (!is.matrix(cost) || any(!is.finite(cost)) || any(cost < 0)) {
    stop("`cost` must be a non-negative finite matrix.", call. = FALSE)
  }
  if (!is.numeric(a) || length(a) != nrow(cost) || any(!is.finite(a)) || any(a < 0)) {
    stop("Query masses must be non-negative and finite.", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != ncol(cost) || any(!is.finite(b)) || any(b < 0)) {
    stop("Reference masses must be non-negative and finite.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1 || !is.finite(epsilon) || epsilon <= 0) {
    stop("`epsilon` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(rho_query) || length(rho_query) != 1 || !is.finite(rho_query) || rho_query <= 0) {
    stop("`rho_query` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(rho_ref) || length(rho_ref) != 1 || !is.finite(rho_ref) || rho_ref <= 0) {
    stop("`rho_ref` must be a positive finite scalar.", call. = FALSE)
  }
  invisible(TRUE)
}

.somalign_row_normalize <- function(plan) {
  totals <- rowSums(plan)
  out <- plan
  nonzero <- totals > 0
  out[nonzero, ] <- out[nonzero, , drop = FALSE] / totals[nonzero]
  out[!nonzero, ] <- 0
  out
}
