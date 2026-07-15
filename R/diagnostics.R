#' Extract somalign diagnostics
#'
#' @param fit A `somalign_fit` object.
#'
#' @return A named list of solver, OT, node, and projection diagnostics.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_diagnostics(fit)
#' @export
somalign_diagnostics <- function(fit) {
  if (!inherits(fit, "somalign_fit")) {
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  }
  fit$diagnostics
}

#' Run an OT sensitivity grid
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Numeric vector of entropic regularisation values.
#' @param rho_query Numeric vector of query-side mass relaxation values.
#' @param rho_ref Numeric vector of reference-side mass relaxation values.
#' @param solver Solver passed to `somalign_fit()`. `"auto"` is accepted as a
#'   compatibility alias for the internal pure-R solver.
#' @param min_match_fraction Minimum match fraction threshold passed to each
#'   `somalign_fit()` call. Default `0.05`.
#' @param confidence_threshold Minimum label confidence for accepted label
#'   transfer. Default `0.6`.
#' @param correction_min_mass Minimum OT mass for a node shift to be applied.
#'   Default `1e-8`.
#' @param max_iter Maximum Sinkhorn iterations. Default `1000`.
#' @param tol Sinkhorn convergence tolerance. Default `1e-7`.
#' @param chunk_size Integer. Number of samples per projection chunk.
#'   `NULL` processes all samples at once. Default `10000L`.
#' @param diagonal_boost Non-negative scalar added to same-node OT costs to
#'   discourage self-transport. Default `0`.
#' @param parallel Logical. When `TRUE`, grid rows are evaluated in parallel
#'   using [BiocParallel::bplapply()] with the registered
#'   `BiocParallel` back-end (see [BiocParallel::register()]). Configure the
#'   back-end before calling this function, e.g.
#'   `BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))`. When
#'   `FALSE` (default) a sequential for-loop is used, which is fully
#'   reproducible across platforms.
#' @param anneal_start,anneal_stages,anneal_factor Annealing-schedule tuning
#'   parameters, used only when `solver = "annealing"`. See [somalign_fit()].
#'
#' @return A data frame with one row per parameter combination, including a
#'   `mutual_information` column (bits; `diagnostics$ot$mutual_information`
#'   for that grid point's fit).
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' somalign_sensitivity_grid(qry, ref,
#'                           epsilon = c(0.05, 0.1),
#'                           rho_query = c(0.5, 1),
#'                           rho_ref = 1)
#' @export
somalign_sensitivity_grid <- function(query,
                                      reference,
                                      epsilon,
                                      rho_query,
                                      rho_ref,
                                      solver = c("internal", "auto", "log_domain", "annealing"),
                                      min_match_fraction = 0.05,
                                      confidence_threshold = 0.6,
                                      correction_min_mass = 1e-8,
                                      max_iter = 1000,
                                      tol = 1e-7,
                                      chunk_size = 10000L,
                                      diagonal_boost = 0,
                                      parallel = FALSE,
                                      anneal_start = 10,
                                      anneal_stages = 10L,
                                      anneal_factor = NULL) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  .somalign_check_prob_scalar(min_match_fraction, "min_match_fraction")
  .somalign_check_prob_scalar(confidence_threshold, "confidence_threshold")
  .somalign_check_pos_scalar(correction_min_mass, "correction_min_mass")
  .somalign_check_pos_int(max_iter, "max_iter")
  .somalign_check_pos_scalar(tol, "tol")
  .somalign_check_pos_int(chunk_size, "chunk_size", allow_null = TRUE)
  .somalign_check_nonneg_scalar(diagonal_boost, "diagonal_boost")
  .somalign_check_flag(parallel, "parallel")
  solver <- match.arg(solver, c("internal", "auto", "log_domain", "annealing"))
  if (identical(solver, "annealing"))
    .somalign_check_anneal_params(anneal_start, anneal_factor, anneal_stages)
  epsilon <- .somalign_validate_grid_vector(epsilon, "epsilon")
  rho_query <- .somalign_validate_grid_vector(rho_query, "rho_query")
  rho_ref <- .somalign_validate_grid_vector(rho_ref, "rho_ref")
  grid <- expand.grid(epsilon = epsilon, rho_query = rho_query,
                      rho_ref = rho_ref, KEEP.OUT.ATTRS = FALSE)

  # Pre-compute the direct projection once: it is identical for every grid point
  # (same scaled_data + reference codebook), saving (K-1) x O(N x nodes) passes.
  direct_cache <- .somalign_project_samples(query$scaled_data, reference,
                                            chunk_size = chunk_size)

  .run_one <- function(i) {
    transport <- .somalign_align_transport(
      query, reference,
      grid$epsilon[i], grid$rho_query[i], grid$rho_ref[i],
      solver, max_iter, tol,
      diagonal_boost = diagonal_boost,
      anneal_start = anneal_start, anneal_factor = anneal_factor,
      anneal_stages = anneal_stages
    )
    fit <- .somalign_finish_fit(
      query, reference, transport,
      min_match_fraction, confidence_threshold, correction_min_mass,
      chunk_size, grid$epsilon[i], grid$rho_query[i], grid$rho_ref[i],
      direct_cache = direct_cache
    )
    .somalign_grid_row_summary(fit, grid$epsilon[i], grid$rho_query[i], grid$rho_ref[i])
  }

  rows <- .somalign_run_grid(nrow(grid), .run_one, parallel)
  do.call(rbind, rows)
}

.somalign_grid_row_summary <- function(fit, epsilon, rho_query, rho_ref) {
  diag <- somalign_diagnostics(fit)
  data.frame(
    epsilon = epsilon,
    rho_query = rho_query,
    rho_ref = rho_ref,
    solver = diag$solver$used,
    transport_mass = diag$ot$transport_mass,
    mean_match_fraction = mean(diag$ot$match_fraction[is.finite(diag$ot$match_fraction)]),
    max_row_mass_error = diag$ot$max_row_mass_error,
    max_col_mass_error = diag$ot$max_col_mass_error,
    mutual_information = diag$ot$mutual_information,
    accepted_label_fraction = mean(fit$label_transfer$accepted),
    outside_direct_fraction = diag$projection$outside_direct_fraction,
    outside_corrected_fraction = diag$projection$outside_corrected_fraction
  )
}

.somalign_run_grid <- function(n, run_one, parallel) {
  if (isTRUE(parallel)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop(
        "Package 'BiocParallel' is required when parallel = TRUE. ",
        "Install it or set parallel = FALSE.",
        call. = FALSE
      )
    }
    rows <- BiocParallel::bplapply(
      seq_len(n),
      run_one
    )
  } else {
    rows <- vector("list", n)
    for (i in seq_len(n)) {
      rows[[i]] <- run_one(i)
    }
  }
  rows
}

.somalign_validate_grid_vector <- function(x, what) {
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x <= 0)) {
    stop("`", what, "` must be a non-empty numeric vector of positive finite values.", call. = FALSE)
  }
  x
}

#' Assess alignment stability across query SOM random seeds
#'
#' Trains a new query SOM for each seed and runs the full alignment pipeline,
#' holding the reference SOM fixed. The returned summary quantifies how much OT
#' alignment statistics vary with query SOM training randomness — the largest
#' uncontrolled variance source in the `somalign` workflow.
#'
#' @param query_data Numeric query data.
#' @param reference A `somalign_reference` object.
#' @param som_seeds Integer vector of random seeds used to train query SOMs.
#' @param epsilon Entropic regularisation passed to [somalign_fit()].
#' @param rho_query,rho_ref Mass-relaxation parameters passed to [somalign_fit()].
#' @param grid Optional `kohonen::somgrid()` for query SOM training.
#' @param rlen Number of SOM training iterations.
#' @param alpha Learning-rate schedule.
#' @param parallel Logical. Use [BiocParallel::bplapply()] when `TRUE`.
#' @param ... Additional arguments passed to [somalign_query()].
#'
#' @return A data frame with one row per seed containing key alignment summary
#'   statistics: `som_seed`, `transport_mass`, `mean_match_fraction`,
#'   `max_row_mass_error`, `accepted_label_fraction`,
#'   `outside_direct_fraction`, `outside_corrected_fraction`,
#'   `mean_correction_norm`, `converged`.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' somalign_som_stability(mat, ref, som_seeds = 1:3,
#'                        grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
#' @export
somalign_som_stability <- function(query_data,
                                   reference,
                                   som_seeds = seq_len(5L),
                                   epsilon = 0.1,
                                   rho_query = 1,
                                   rho_ref = 1,
                                   grid = NULL,
                                   rlen = 100,
                                   alpha = c(0.05, 0.01),
                                   parallel = FALSE,
                                   ...) {
  .somalign_check_reference(reference)
  if (!is.numeric(som_seeds) || length(som_seeds) == 0L) {
    stop("`som_seeds` must be a non-empty numeric vector.", call. = FALSE)
  }
  som_seeds <- as.integer(som_seeds)
  .somalign_check_pos_scalar(epsilon, "epsilon")
  .somalign_check_pos_scalar(rho_query, "rho_query")
  .somalign_check_pos_scalar(rho_ref, "rho_ref")
  .somalign_check_pos_int(rlen, "rlen")
  if (!is.numeric(alpha) || length(alpha) != 2L || any(!is.finite(alpha)) || any(alpha <= 0))
    stop("`alpha` must be a numeric vector of two positive finite values.", call. = FALSE)
  .somalign_check_flag(parallel, "parallel")

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  .run_one <- function(i) {
    seed <- som_seeds[i]
    set.seed(seed)
    qry <- somalign_query(query_data, reference, grid = grid,
                          rlen = rlen, alpha = alpha, ...)
    fit <- somalign_fit(qry, reference, epsilon = epsilon,
                        rho_query = rho_query, rho_ref = rho_ref)
    .somalign_stability_row_summary(fit, seed)
  }

  rows <- .somalign_run_grid(length(som_seeds), .run_one, parallel)
  do.call(rbind, rows)
}

.somalign_stability_row_summary <- function(fit, seed) {
  diag <- somalign_diagnostics(fit)
  allowed <- diag$nodes$correction_allowed
  data.frame(
    som_seed              = seed,
    transport_mass        = diag$ot$transport_mass,
    mean_match_fraction   = mean(diag$ot$match_fraction[is.finite(diag$ot$match_fraction)]),
    max_row_mass_error    = diag$ot$max_row_mass_error,
    accepted_label_fraction   = mean(fit$label_transfer$accepted),
    outside_direct_fraction   = diag$projection$outside_direct_fraction,
    outside_corrected_fraction = diag$projection$outside_corrected_fraction,
    mean_correction_norm  = if (any(allowed)) mean(diag$nodes$correction_norm[allowed]) else NA_real_,
    converged             = diag$solver$converged
  )
}

# ---------------------------------------------------------------------------
# Epsilon phase-transition sweep + principled epsilon selection
# ---------------------------------------------------------------------------

.somalign_default_eps_grid <- function(n = 25L) {
  exp(seq(log(1e-3), log(5), length.out = n))
}

# Finite-difference susceptibility d(Phi)/d(log_epsilon) on a (Phi, log_eps)
# grid sorted by log_eps ascending. NA at the first and last positions and
# wherever fewer than 3 points are available.
.somalign_susceptibility <- function(Phi, log_eps) {
  n <- length(Phi)
  chi <- rep(NA_real_, n)
  if (n < 3L) return(chi)
  idx <- 2:(n - 1)
  chi[idx] <- (Phi[idx + 1] - Phi[idx - 1]) / (log_eps[idx + 1] - log_eps[idx - 1])
  chi
}

.somalign_locate_epsilon_c <- function(table) {
  chi <- table$susceptibility
  finite_idx <- which(is.finite(chi))
  if (length(finite_idx) == 0L) return(NA_real_)
  k_star <- finite_idx[which.max(chi[finite_idx])]
  table$epsilon[k_star]
}

#' Epsilon phase-transition sweep for principled epsilon selection
#'
#' Runs the Sinkhorn OT solve across a log-spaced epsilon grid without the
#' per-cell projection step, computing the transport-plan order parameter
#' (mean fractional reference-node usage, a perplexity-based measure), its
#' susceptibility (rate of change with log epsilon), the dual free energy,
#' and the mutual information between query and reference nodes at each
#' grid point.
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon_grid Numeric vector of epsilon values. If `NULL` (default),
#'   a log-spaced grid of `n_grid` values from 1e-3 to 5 is used.
#' @param n_grid Integer. Grid size when `epsilon_grid = NULL`. Default `25`.
#' @param rho_query,rho_ref Mass relaxation parameters passed to the OT solver.
#' @param solver Sinkhorn solver. Default `"log_domain"` (required for
#'   `log_Z`; other solvers fill `log_Z` with `NA`).
#' @param max_iter,tol Sinkhorn convergence parameters.
#' @param diagonal_boost Non-negative cost reduction on nearest-reference-node
#'   entries. Default `0`.
#' @param label_guided Logical; see [somalign_fit()].
#' @param parallel Logical; see [somalign_sensitivity_grid()].
#'
#' @return A list of class `"somalign_epsilon_sweep"` with components
#'   `table`, `epsilon_c`, `epsilon_rec`, `cost_scale`. The `table` data
#'   frame has one row per epsilon value with columns `epsilon`,
#'   `log_epsilon`, `Phi`, `susceptibility`, `log_Z`, `mutual_information`,
#'   `conditional_entropy_mean`, `expected_cost`, `transport_mass`,
#'   `iterations`, `converged`.
#'
#' @details
#' The order parameter \eqn{\Phi(\epsilon)} is the mean effective fraction of
#' reference nodes used per query node: \eqn{\Phi = \mathrm{mean}_i(2^{H_i}) /
#' K}, where \eqn{H_i} is the (bit) conditional entropy of the row-normalised
#' transport plan for query node \eqn{i} and \eqn{K} is the number of
#' reference nodes. As \eqn{\epsilon \to 0}, \eqn{\Phi \to 1/K}; as
#' \eqn{\epsilon \to \infty}, \eqn{\Phi \to 1}. The susceptibility
#' \eqn{\chi = d\Phi/d\log\epsilon} peaks at the critical epsilon
#' \eqn{\epsilon_c}, marking the crossover between a localised
#' (transport-cost-dominated) and delocalised (entropy-dominated) plan.
#' `epsilon_rec` (\eqn{0.3 \, \epsilon_c}) keeps the plan in the ordered phase
#' with a safety margin.
#'
#' The sweep avoids the per-cell projection step, so it runs approximately
#' one OT solve per epsilon value -- much faster than
#' [somalign_sensitivity_grid()] for the same epsilon range.
#'
#' @seealso [somalign_fit()], [somalign_select_epsilon()],
#'   [somalign_sensitivity_grid()], [somalign_diagnostics()]
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat + 0.5, ref,
#'                       grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
#' sw <- somalign_epsilon_sweep(qry, ref, n_grid = 8)
#' sw$epsilon_c
#' sw$epsilon_rec
#' @export
somalign_epsilon_sweep <- function(query, reference,
                                   epsilon_grid = NULL, n_grid = 25L,
                                   rho_query = 1, rho_ref = 1,
                                   solver = c("log_domain", "internal", "auto"),
                                   max_iter = 1000, tol = 1e-7,
                                   diagonal_boost = 0, label_guided = FALSE,
                                   parallel = FALSE) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  .somalign_check_pos_scalar(rho_query, "rho_query")
  .somalign_check_pos_scalar(rho_ref, "rho_ref")
  .somalign_check_nonneg_scalar(diagonal_boost, "diagonal_boost")
  .somalign_check_flag(parallel, "parallel")
  solver <- match.arg(solver)

  label_mask <- .somalign_sweep_label_mask(query, reference, label_guided)
  if (is.null(epsilon_grid)) epsilon_grid <- .somalign_default_eps_grid(n_grid)
  epsilon_grid <- sort(unique(.somalign_validate_grid_vector(epsilon_grid, "epsilon_grid")))
  if (length(epsilon_grid) < 3L) {
    message("somalign_epsilon_sweep: fewer than 3 grid points; epsilon_c cannot be estimated reliably.")
  }

  K <- nrow(reference$codebook)
  .run_one <- function(i) {
    .somalign_sweep_table_row(
      query, reference, epsilon_grid[i], rho_query, rho_ref,
      solver, max_iter, tol, diagonal_boost, label_mask, K
    )
  }
  rows <- .somalign_run_grid(length(epsilon_grid), .run_one, parallel)
  table <- do.call(rbind, rows)
  table$susceptibility <- .somalign_susceptibility(table$Phi, table$log_epsilon)

  epsilon_c <- .somalign_locate_epsilon_c(table)
  epsilon_rec <- if (is.finite(epsilon_c)) 0.3 * epsilon_c else NA_real_

  structure(
    list(table = table, epsilon_c = epsilon_c, epsilon_rec = epsilon_rec,
         cost_scale = table$cost_scale[1]),
    class = "somalign_epsilon_sweep"
  )
}

.somalign_sweep_label_mask <- function(query, reference, label_guided) {
  if (!isTRUE(label_guided)) return(NULL)
  if (is.null(query$label_prob) || is.null(reference$label_prob)) {
    stop("label_guided = TRUE requires both query$label_prob and reference$label_prob.",
         call. = FALSE)
  }
  .somalign_build_label_mask(query$label_prob, reference$label_prob)
}

.somalign_sweep_table_row <- function(query, reference, eps, rho_query, rho_ref,
                                      solver, max_iter, tol, diagonal_boost,
                                      label_mask, K) {
  r <- .somalign_ot_sweep_one(query, reference, eps, rho_query, rho_ref,
                              solver, max_iter, tol,
                              diagonal_boost = diagonal_boost,
                              label_mask = label_mask)
  Phi <- .somalign_order_parameter(r$conditional_entropy, K)
  data.frame(
    epsilon = r$epsilon,
    log_epsilon = log(r$epsilon),
    Phi = Phi,
    log_Z = if (is.null(r$log_Z)) NA_real_ else r$log_Z,
    mutual_information = r$mutual_information,
    conditional_entropy_mean = mean(r$conditional_entropy, na.rm = TRUE),
    expected_cost = r$expected_cost,
    transport_mass = r$transport_mass,
    cost_scale = r$cost_scale,
    iterations = r$iterations,
    converged = r$converged
  )
}

#' @method print somalign_epsilon_sweep
#' @export
print.somalign_epsilon_sweep <- function(x, ...) {
  cat(sprintf("somalign epsilon sweep [%d points]\n", nrow(x$table)))
  cat(sprintf("  epsilon range      : %.4g - %.4g\n",
              min(x$table$epsilon), max(x$table$epsilon)))
  cat(sprintf("  critical epsilon   : %s\n",
              if (is.finite(x$epsilon_c)) sprintf("%.4g", x$epsilon_c) else "NA"))
  cat(sprintf("  recommended epsilon (0.3x critical): %s\n",
              if (is.finite(x$epsilon_rec)) sprintf("%.4g", x$epsilon_rec) else "NA"))
  cat(sprintf("  cost_scale         : %.4g\n", x$cost_scale))
  invisible(x)
}

#' @method plot somalign_epsilon_sweep
#' @export
plot.somalign_epsilon_sweep <- function(x, ...) {
  .somalign_plot_eps_sweep(x)
}

#' Select epsilon from an epsilon-sweep curve
#'
#' Sweeps a grid of epsilon values via [somalign_epsilon_sweep()] and selects
#' a recommended value using one of three criteria: the susceptibility
#' (critical-epsilon) peak, the mutual-information-vs-cost elbow, or the
#' smallest epsilon retaining a target fraction of maximum mutual information.
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Numeric vector of candidate epsilon values, or `NULL` to use
#'   the default grid. Default `c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0)`.
#' @param rho_query,rho_ref Marginal relaxations passed to the OT solver.
#' @param method Character. `"critical"` (default) uses the susceptibility
#'   peak's recommended epsilon (`0.3 * epsilon_c`). `"elbow"` uses the
#'   max-second-difference rule on the mutual-information-vs-cost curve.
#'   `"entropy_fraction"` chooses the largest (most-regularized) epsilon that
#'   still retains `entropy_fraction * max(mutual_information)`.
#' @param entropy_fraction Numeric in (0, 1]. Target fraction of maximum
#'   mutual information when `method = "entropy_fraction"`. Default `0.90`.
#' @param solver,max_iter,tol Sinkhorn solver parameters.
#' @param diagonal_boost,label_guided,parallel See [somalign_epsilon_sweep()].
#'
#' @return A list of class `"somalign_epsilon_selection"` with
#'   `selected_epsilon`, `curve` (the sweep's `table`), and `method`.
#' @seealso [somalign_epsilon_sweep()]
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' somalign_select_epsilon(qry, ref, epsilon = c(0.05, 0.1, 0.2))
#' @export
somalign_select_epsilon <- function(query, reference,
                                    epsilon = c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0),
                                    rho_query = 1, rho_ref = 1,
                                    method = c("critical", "elbow", "entropy_fraction"),
                                    entropy_fraction = 0.90,
                                    solver = c("log_domain", "internal", "auto"),
                                    max_iter = 1000, tol = 1e-7,
                                    diagonal_boost = 0, label_guided = FALSE,
                                    parallel = FALSE) {
  method <- match.arg(method)
  .somalign_check_prob_scalar(entropy_fraction, "entropy_fraction")

  sweep <- somalign_epsilon_sweep(
    query, reference, epsilon_grid = epsilon,
    rho_query = rho_query, rho_ref = rho_ref, solver = solver,
    max_iter = max_iter, tol = tol, diagonal_boost = diagonal_boost,
    label_guided = label_guided, parallel = parallel
  )

  selected <- switch(
    method,
    critical = if (is.finite(sweep$epsilon_rec)) sweep$epsilon_rec else NA_real_,
    elbow = .somalign_select_elbow(sweep$table),
    entropy_fraction = .somalign_select_entropy_fraction(sweep$table, entropy_fraction)
  )

  structure(
    list(selected_epsilon = selected, curve = sweep$table, method = method,
         sweep = sweep),
    class = "somalign_epsilon_selection"
  )
}

# Kneedle-style: sort by expected_cost ascending; find the index of maximum
# second difference of mutual_information (greatest deceleration of MI gain).
.somalign_select_elbow <- function(curve) {
  ord <- order(curve$expected_cost)
  mi <- curve$mutual_information[ord]
  if (length(mi) < 3L || all(is.na(mi))) {
    return(curve$epsilon[ord[which.max(mi)]])
  }
  d2 <- diff(diff(mi))
  elbow_idx <- which.max(d2) + 1L
  curve$epsilon[ord[elbow_idx]]
}

.somalign_select_entropy_fraction <- function(curve, fraction) {
  valid <- is.finite(curve$mutual_information)
  if (!any(valid)) return(NA_real_)
  target <- fraction * max(curve$mutual_information[valid])
  reached <- valid & curve$mutual_information >= target
  if (!any(reached)) return(curve$epsilon[which.max(curve$mutual_information)])
  curve$epsilon[max(which(reached))]
}
