.somalign_solve_ot <- function(cost,
                               a,
                               b,
                               epsilon,
                               rho_query,
                               rho_ref,
                               solver,
                               max_iter,
                               tol,
                               anneal_start = 10,
                               anneal_factor = NULL,
                               anneal_stages = 10L) {
  requested_solver <- solver
  notes <- character()
  .somalign_validate_ot_inputs(cost, a, b, epsilon, rho_query, rho_ref)

  if (solver == "auto") {
    solver <- "internal"
    notes <- c(notes, "`solver = \"auto\"` uses the internal generalized Sinkhorn solver.")
  }

  if (identical(solver, "annealing")) {
    internal <- .somalign_solve_annealing(
      cost, a, b, epsilon, rho_query, rho_ref, max_iter, tol,
      anneal_start, anneal_factor, anneal_stages
    )
    return(list(
      plan = internal$plan,
      solver = "annealing",
      requested_solver = requested_solver,
      notes = notes,
      iterations = internal$iterations,
      converged = internal$converged,
      final_delta = internal$final_delta,
      log_Z = internal$log_Z,
      anneal_schedule = internal$schedule,
      anneal_stage_info = internal$stage_info
    ))
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
    final_delta = internal$final_delta,
    log_Z = internal$log_Z,
    anneal_schedule = NULL,
    anneal_stage_info = NULL
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
  list(plan = plan, iterations = iterations, converged = converged,
       final_delta = final_delta, log_Z = NA_real_)
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
                                         tol,
                                         f_init = NULL,
                                         g_init = NULL) {
  # tau_a/tau_b are recomputed from `epsilon` on every call -- this must never
  # be hoisted outside a per-epsilon call, since the annealing driver
  # (.somalign_solve_annealing) calls this function once per cooling stage
  # with a different epsilon each time.
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_ref   / (rho_ref   + epsilon)

  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  f <- if (!is.null(f_init)) f_init else numeric(length(a))
  g <- if (!is.null(g_init)) g_init else numeric(length(b))
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

    # Convergence is measured only over finite potentials. Zero-mass nodes
    # produce f_i = -Inf (or g_j = -Inf) permanently, which is structurally
    # correct (exp(-Inf) = 0 in the plan) but causes abs(-Inf - -Inf) = NaN,
    # poisoning max() even though all meaningful potentials have converged.
    # Restricting to is.finite() ratios mirrors how .somalign_logsumexp already
    # handles -Inf entries, and does not mask genuine divergence: a truly
    # degenerate solve produces all-non-finite ratios, so finite_ratios is
    # empty and delta falls through to the Inf guard (no convergence declared).
    f_ratio <- abs(f - f_old) / pmax(1, abs(f_old))
    g_ratio <- abs(g - g_old) / pmax(1, abs(g_old))
    finite_ratios <- c(f_ratio[is.finite(f_ratio)], g_ratio[is.finite(g_ratio)])
    delta <- if (length(finite_ratios) == 0L) Inf else max(finite_ratios)
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

  # Entropic (Sinkhorn) dual objective, symmetric in both marginals:
  #   log_Z = <a, f> + <b, g> - epsilon * sum_ij P_ij.
  # This is exact for balanced transport and an approximation under strong
  # unbalancing (it omits the rho KL-relaxation terms of the unbalanced dual);
  # it replaces the earlier row-only `epsilon * sum_i logsumexp_j(...)`, which
  # dropped the g/column-marginal contribution entirely. Zero-mass marginals
  # carry -Inf potentials, so guard the 0 * -Inf = NaN product by summing only
  # over positive-mass entries (the a>0/b>0 subsetting does this; sum() over an
  # empty vector is 0); `plan` is already finite (non-finite -> 0).
  fa <- sum(a[a > 0] * f[a > 0])
  gb <- sum(b[b > 0] * g[b > 0])
  log_Z <- fa + gb - epsilon * sum(plan)

  list(plan = plan, iterations = iterations, converged = converged,
       final_delta = final_delta, f = f, g = g, log_Z = log_Z)
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

# Geometric epsilon-cooling schedule from `anneal_start * epsilon_target` down
# to `epsilon_target` over `anneal_stages` steps. The last entry is clamped to
# exactly `epsilon_target` to avoid floating-point drift.
.somalign_cooling_schedule <- function(epsilon_target, anneal_start,
                                       anneal_stages, anneal_factor) {
  if (anneal_stages == 1L) return(epsilon_target)
  eps_0 <- anneal_start * epsilon_target
  r <- if (is.null(anneal_factor)) {
    (epsilon_target / eps_0) ^ (1 / (anneal_stages - 1))
  } else {
    anneal_factor
  }
  schedule <- eps_0 * r ^ seq(0, anneal_stages - 1)
  schedule[anneal_stages] <- epsilon_target
  schedule
}

# Simulated-annealing Sinkhorn: runs the log-domain solver across the cooling
# schedule, warm-starting each stage from the previous stage's dual
# potentials. Interior stages use a looser iteration budget and tolerance;
# the final stage (at the target epsilon) runs to full max_iter/tol.
.somalign_solve_annealing <- function(cost, a, b, epsilon,
                                      rho_query, rho_ref,
                                      max_iter, tol,
                                      anneal_start, anneal_factor,
                                      anneal_stages) {
  schedule <- .somalign_cooling_schedule(epsilon, anneal_start,
                                         anneal_stages, anneal_factor)
  n_stages <- length(schedule)
  inner_iter <- max(ceiling(max_iter / n_stages), 20L)
  inner_tol <- tol * 10
  total_iter <- 0L
  stage_info <- vector("list", n_stages)

  f <- NULL
  g <- NULL
  result <- NULL
  for (k in seq_len(n_stages)) {
    eps_k <- schedule[k]
    is_final <- (k == n_stages)
    result <- .somalign_solve_internal_log(
      cost = cost, a = a, b = b, epsilon = eps_k,
      rho_query = rho_query, rho_ref = rho_ref,
      max_iter = if (is_final) max_iter else inner_iter,
      tol = if (is_final) tol else inner_tol,
      f_init = f, g_init = g
    )
    f <- result$f
    g <- result$g
    total_iter <- total_iter + result$iterations
    stage_info[[k]] <- list(epsilon = eps_k, iterations = result$iterations,
                            converged = result$converged,
                            final_delta = result$final_delta)
  }

  list(
    plan = result$plan,
    iterations = total_iter,
    converged = result$converged,
    final_delta = result$final_delta,
    log_Z = result$log_Z,
    stage_info = stage_info,
    schedule = schedule
  )
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

# Derive the per-query-node correspondence quantities from a raw transport
# plan and the query node masses: the row-normalised correspondence, the
# transported row/column mass, the match-mass ratio (transported vs available
# node mass), and the match fraction (ratio clamped to <= 1). Shared by
# .somalign_align_transport() (the interactive fit path, which additionally
# emits a match_mass_ratio > 1 message) and .somalign_sweep_topology_row()
# (the silent epsilon-sweep diagnostic) so the two never drift.
.somalign_plan_to_correspondence <- function(plan, node_masses) {
  row_mass <- rowSums(plan)
  match_mass_ratio <- ifelse(node_masses > 0, row_mass / node_masses, 0)
  list(
    correspondence = .somalign_row_normalize(plan),
    row_mass = row_mass,
    col_mass = colSums(plan),
    match_mass_ratio = match_mass_ratio,
    match_fraction = pmin(match_mass_ratio, 1)
  )
}

# Mutual information I(query node; reference node) from the raw (possibly
# unnormalized) transport plan, plus per-query-node conditional entropy
# H(ref | query = i) -- both in bits (log2) -- and the plan's expected cost
# under `cost` (same units as `cost`; NA if `cost` is not supplied). Marginals
# are the *empirical* row/col sums of the plan (what was actually
# transported), not the input node masses -- correct under UOT, where mass
# can be destroyed.
.somalign_plan_mutual_information <- function(plan, cost = NULL) {
  tiny <- .Machine$double.eps
  total <- sum(plan)
  if (!is.finite(total) || total <= tiny) {
    m <- nrow(plan)
    return(list(
      mutual_information = NA_real_,
      conditional_entropy = rep(NA_real_, m),
      expected_cost = NA_real_
    ))
  }
  p_joint <- plan / total
  a_hat <- rowSums(p_joint)
  b_hat <- colSums(p_joint)
  outer_ab <- outer(a_hat, b_hat)
  nz <- p_joint > tiny & outer_ab > tiny
  mi_bits <- sum(p_joint[nz] * log2(p_joint[nz] / outer_ab[nz]))
  mi_bits <- max(mi_bits, 0)

  row_totals <- rowSums(plan)
  cond_ent <- vapply(seq_len(nrow(plan)), function(i) {
    if (row_totals[i] <= tiny) return(NA_real_)
    p_cond <- plan[i, ] / row_totals[i]
    nz_j <- p_cond > tiny
    if (!any(nz_j)) return(0)
    -sum(p_cond[nz_j] * log2(p_cond[nz_j]))
  }, numeric(1))

  expected_cost <- if (!is.null(cost)) sum(p_joint * cost) else NA_real_
  list(
    mutual_information = mi_bits,
    conditional_entropy = cond_ent,
    expected_cost = expected_cost
  )
}

# Perplexity-based order parameter: mean effective fraction of reference
# nodes used per query node, from per-node conditional entropy in bits.
# 2^H is the "perplexity" (effective node count) of the row-normalized plan;
# ranges from 1 (one node used) to K (all nodes used equally). Phi in
# (1/K, 1]; NA when every row is NA (e.g. all mass destroyed).
.somalign_order_parameter <- function(conditional_entropy, K) {
  perplexity <- 2 ^ conditional_entropy
  finite <- is.finite(perplexity)
  if (!any(finite)) return(NA_real_)
  mean(perplexity[finite]) / K
}
