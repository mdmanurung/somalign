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
#' @param parallel Logical. When `TRUE`, grid rows are evaluated in parallel
#'   using [BiocParallel::bplapply()] with the registered
#'   `BiocParallel` back-end (see [BiocParallel::register()]). Configure the
#'   back-end before calling this function, e.g.
#'   `BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))`. When
#'   `FALSE` (default) a sequential for-loop is used, which is fully
#'   reproducible across platforms.
#' @param ... Additional arguments passed to `somalign_fit()`.
#'
#' @return A data frame with one row per parameter combination.
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
                                      solver = c("internal", "auto"),
                                      min_match_fraction = 0.05,
                                      confidence_threshold = 0.6,
                                      correction_min_mass = 1e-8,
                                      max_iter = 1000,
                                      tol = 1e-7,
                                      chunk_size = 10000L,
                                      diagonal_boost = 0,
                                      parallel = FALSE) {
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
  solver <- match.arg(solver)
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
      diagonal_boost = diagonal_boost
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
