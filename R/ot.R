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
    if (.somalign_pot_available()) {
      solver <- "pot"
    } else {
      solver <- "internal"
      notes <- c(notes, "POT unavailable; used internal generalized Sinkhorn solver.")
    }
  }

  if (solver == "pot") {
    if (!.somalign_pot_available()) {
      stop("`solver = \"pot\"` requires Python POT with `ot.unbalanced` importable.", call. = FALSE)
    }
    plan <- .somalign_solve_pot(cost, a, b, epsilon, rho_query, rho_ref)
    iterations <- NA_integer_
  } else {
    internal <- .somalign_solve_internal(cost, a, b, epsilon, rho_query, rho_ref, max_iter, tol)
    plan <- internal$plan
    iterations <- internal$iterations
  }

  list(
    plan = plan,
    solver = solver,
    requested_solver = requested_solver,
    notes = notes,
    iterations = iterations
  )
}

.somalign_pot_available <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  isTRUE(tryCatch(reticulate::py_module_available("ot.unbalanced"), error = function(e) FALSE))
}

.somalign_solve_pot <- function(cost, a, b, epsilon, rho_query, rho_ref) {
  ot <- reticulate::import("ot", delay_load = FALSE)
  fn <- ot$unbalanced$sinkhorn_unbalanced
  plan <- fn(
    a = a,
    b = b,
    M = cost,
    reg = epsilon,
    reg_m = c(rho_query, rho_ref),
    reg_type = "entropy"
  )
  plan <- as.matrix(plan)
  storage.mode(plan) <- "double"
  if (any(!is.finite(plan)) || any(plan < -sqrt(.Machine$double.eps))) {
    stop("POT returned an invalid transport plan.", call. = FALSE)
  }
  pmax(plan, 0)
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
  k <- exp(-cost / epsilon)
  k <- pmax(k, tiny)
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_ref / (rho_ref + epsilon)

  u <- rep(1, length(a))
  v <- rep(1, length(b))
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
