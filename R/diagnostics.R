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
                                      parallel = FALSE,
                                      ...) {
  solver <- match.arg(solver)
  epsilon <- .somalign_validate_grid_vector(epsilon, "epsilon")
  rho_query <- .somalign_validate_grid_vector(rho_query, "rho_query")
  rho_ref <- .somalign_validate_grid_vector(rho_ref, "rho_ref")
  grid <- expand.grid(
    epsilon = epsilon,
    rho_query = rho_query,
    rho_ref = rho_ref,
    KEEP.OUT.ATTRS = FALSE
  )

  .run_one <- function(i) {
    fit <- somalign_fit(
      query = query,
      reference = reference,
      epsilon = grid$epsilon[i],
      rho_query = grid$rho_query[i],
      rho_ref = grid$rho_ref[i],
      solver = solver,
      ...
    )
    diag <- somalign_diagnostics(fit)
    data.frame(
      epsilon = grid$epsilon[i],
      rho_query = grid$rho_query[i],
      rho_ref = grid$rho_ref[i],
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

  if (isTRUE(parallel)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop(
        "Package 'BiocParallel' is required when parallel = TRUE. ",
        "Install it or set parallel = FALSE.",
        call. = FALSE
      )
    }
    rows <- BiocParallel::bplapply(
      seq_len(nrow(grid)),
      .run_one
    )
  } else {
    rows <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      rows[[i]] <- .run_one(i)
    }
  }
  do.call(rbind, rows)
}

.somalign_validate_grid_vector <- function(x, what) {
  if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) || any(x <= 0)) {
    stop("`", what, "` must be a non-empty numeric vector of positive finite values.", call. = FALSE)
  }
  x
}
