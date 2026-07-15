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
  log_domain <- identical(solver, "log_domain")
  internal <- .somalign_solve_internal(cost, a, b, epsilon, rho_query, rho_ref,
                                       max_iter, tol, log_domain = log_domain)
  plan <- internal$plan
  iterations <- internal$iterations

  list(
    plan = plan,
    solver = solver,
    requested_solver = requested_solver,
    notes = notes,
    iterations = iterations,
    converged = internal$converged,
    final_delta = internal$final_delta
  )
}

.somalign_solve_internal <- function(cost,
                                     a,
                                     b,
                                     epsilon,
                                     rho_query,
                                     rho_ref,
                                     max_iter,
                                     tol,
                                     log_domain = FALSE) {
  if (log_domain) {
    return(.somalign_solve_internal_log(cost, a, b, epsilon, rho_query, rho_ref,
                                        max_iter, tol))
  }
  tiny <- .Machine$double.xmin
  k <- .somalign_sinkhorn_kernel(cost, epsilon, tiny)
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

  final_delta <- delta
  .somalign_warn_convergence(final_delta, iterations, max_iter, tol)
  converged <- is.finite(final_delta) && final_delta < tol

  plan <- sweep(sweep(k, 1, u, "*"), 2, v, "*")
  plan[!is.finite(plan)] <- 0
  plan <- pmax(plan, 0)
  list(plan = plan, iterations = iterations, converged = converged, final_delta = final_delta)
}

.somalign_logsumexp <- function(x) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0L) return(-Inf)
  m <- max(finite_x)
  m + log(sum(exp(x - m)))
}

# Log-domain unbalanced Sinkhorn. Avoids kernel underflow by working with log
# potentials f (M-vector) and g (K-vector) throughout, recovering the plan as
# P_ij = exp((f_i + g_j - C_ij) / eps) only at the end.
.somalign_solve_internal_log <- function(cost,
                                         a,
                                         b,
                                         epsilon,
                                         rho_query,
                                         rho_ref,
                                         max_iter,
                                         tol) {
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_ref   / (rho_ref   + epsilon)

  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  f <- numeric(length(a))
  g <- numeric(length(b))
  cost_over_eps <- cost / epsilon

  delta <- Inf
  iterations <- max_iter

  for (iter in seq_len(max_iter)) {
    f_old <- f
    g_old <- g

    # f_i = tau_a * (eps * log_a_i - eps * logsumexp_j((g_j - C_ij) / eps))
    M_g <- sweep(-cost_over_eps, 2, g / epsilon, "+")
    lse_g <- apply(M_g, 1, .somalign_logsumexp)
    f <- tau_a * (epsilon * log_a - epsilon * lse_g)
    f[is.nan(f)] <- 0

    # g_j = tau_b * (eps * log_b_j - eps * logsumexp_i((f_i - C_ij) / eps))
    M_f <- sweep(-cost_over_eps, 1, f / epsilon, "+")
    lse_f <- apply(M_f, 2, .somalign_logsumexp)
    g <- tau_b * (epsilon * log_b - epsilon * lse_f)
    g[is.nan(g)] <- 0

    delta <- max(
      abs(f - f_old) / pmax(1, abs(f_old)),
      abs(g - g_old) / pmax(1, abs(g_old))
    )
    if (is.finite(delta) && delta < tol) {
      iterations <- iter
      break
    }
  }

  final_delta <- delta
  .somalign_warn_convergence(final_delta, iterations, max_iter, tol)
  converged <- is.finite(final_delta) && final_delta < tol

  log_plan <- sweep(sweep(-cost_over_eps, 1, f / epsilon, "+"), 2, g / epsilon, "+")
  plan <- exp(log_plan)
  plan[!is.finite(plan)] <- 0
  plan <- pmax(plan, 0)
  list(plan = plan, iterations = iterations, converged = converged, final_delta = final_delta)
}

.somalign_sinkhorn_kernel <- function(cost, epsilon, tiny) {
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
  # A row that underflows entirely is flooded to a constant after flooring,
  # which destroys the cost ordering for that query node. The aggregate check
  # above can miss a single such row (e.g. 1/10000 entries), so flag it here.
  zero_rows <- which(rowSums(k_raw) == 0)
  if (length(zero_rows) > 0L) {
    shown <- zero_rows[seq_len(min(10L, length(zero_rows)))]
    warning(
      sprintf(
        "%d query node(s) have an all-underflow Sinkhorn kernel row (e.g. rows %s); ",
        length(zero_rows), paste(shown, collapse = ", ")
      ),
      "their transport-cost ordering is lost after flooring. ",
      "Use solver = \"log_domain\" to avoid this.",
      call. = FALSE
    )
  }
  pmax(k_raw, tiny)
}

.somalign_warn_convergence <- function(final_delta, iterations, max_iter, tol) {
  if (!is.finite(final_delta)) {
    warning(
      sprintf(
        "Sinkhorn solver produced a non-finite iterate delta after %d iterations. ",
        iterations
      ),
      "The solve may be degenerate (e.g. all-zero masses or extreme cost/epsilon ratio). ",
      "Check diagnostics$solver$final_delta and diagnostics$ot.",
      call. = FALSE
    )
  } else if (iterations == max_iter && final_delta >= tol) {
    warning(
      sprintf(
        "Sinkhorn solver did not converge after %g iterations (final delta = %.3e). ",
        max_iter, final_delta
      ),
      "Consider increasing max_iter, raising epsilon, or reducing rho_query / rho_ref.",
      call. = FALSE
    )
  }
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
  if (sum(a) == 0) {
    warning(
      "All query node masses are zero; the transport plan will be all zeros.",
      call. = FALSE
    )
  }
  if (sum(b) == 0) {
    warning(
      "All reference node masses are zero; the transport plan will be all zeros.",
      call. = FALSE
    )
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
