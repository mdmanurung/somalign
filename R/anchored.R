#' Align a query SOM to a reference SOM using anchor sample pairs
#'
#' A variant of [somalign_fit()] for the case where a set of samples has been
#' measured in **both** the old batch (reference space) and the new batch
#' (query space). These *anchor pairs* are used to build a per-node-pair
#' correspondence count matrix, which is subtracted from the normalized OT
#' cost before the Sinkhorn solve. This makes transport along anchor-supported
#' routes cheaper, biasing the OT plan toward pairings that are consistent
#' with the observed per-sample batch displacement — while still solving a
#' valid optimal transport problem over the full codebook.
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Entropic regularisation strength (see [somalign_fit()]).
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver Sinkhorn solver variant. See [somalign_fit()].
#' @param anchor_old Numeric matrix (n_anchors × p). Old-batch measurements
#'   of the anchor samples, in the same feature space as `reference`. Also
#'   accepts a data frame of numeric columns.
#' @param anchor_new Numeric matrix (n_anchors × p). New-batch measurements
#'   of the **same** anchor samples, in the same feature space as `reference`.
#'   Rows of `anchor_old` and `anchor_new` must correspond to the same
#'   biological units. Also accepts a data frame of numeric columns.
#' @param rho_anchor Non-negative scalar. Controls how strongly anchor pairs
#'   bias the OT cost. At `rho_anchor = 0` the result equals [somalign_fit()].
#'   Larger values reduce the effective cost for anchor-supported node pairs,
#'   concentrating the transport plan on those routes. Typical range: 0.5–3.
#' @param min_match_fraction Minimum transported fraction for label transfer.
#' @param confidence_threshold Minimum top-label probability for label transfer.
#' @param correction_min_mass Minimum transported mass for a node correction.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Sinkhorn convergence tolerance.
#' @param chunk_size Integer. Samples projected per chunk. Default `10000L`.
#'
#' @details
#' **Cost modification.** Let \eqn{C} be the M×K codebook distance matrix
#' normalised by its median positive entry (as in [somalign_fit()]). The
#' anchor pairs are projected onto the query codebook (old batch) and
#' reference codebook (new batch), yielding a count matrix \eqn{A} where
#' \eqn{A_{kl}} is the number of anchor pairs whose old measurement maps to
#' query node \eqn{k} and new measurement maps to reference node \eqn{l}.
#' The modified cost is
#' \deqn{\tilde{C}_{kl} = \max\!\bigl(C_{kl} - \rho_{\mathrm{anchor}} \cdot
#'   A_{kl} / n_{\mathrm{anchors}},\; 0\bigr).}
#' Pairs with many anchor observations get cost reduced toward zero (free
#' transport), while uncovered pairs retain their original cost. Non-negativity
#' is enforced by the \eqn{\max(\cdot, 0)} clamp.
#'
#' **Clamp behaviour at large `rho_anchor`.** When the anchor bonus exceeds
#' \eqn{C_{kl}}, the effective cost is clamped to zero. All such pairs then
#' have identical effective cost and the transport mass among them is determined
#' by entropic regularisation alone rather than by relative anchor counts.
#' The clamp is required to keep costs non-negative; at very large `rho_anchor`
#' the plan for anchor-covered pairs becomes more entropic, not more
#' concentrated. A practical upper bound is `rho_anchor * max(A) / n_anchors
#' <= 1`, i.e., even the most-supported pair reduces cost by at most one
#' median-distance unit.
#'
#' **Fallback for uncovered nodes.** Query nodes with no anchor samples retain
#' their original pairwise costs, so the transport plan for those nodes is
#' determined entirely by the OT objective — the same as [somalign_fit()].
#' Inspect `$anchors$coverage_fraction` to see what fraction of query nodes
#' had at least one anchor pair.
#'
#' **Return value.** The object has class `c("somalign_anchored_fit",
#' "somalign_fit")`, so all downstream functions that accept a
#' `somalign_fit` object ([somalign_results()], [somalign_diagnostics()])
#' work unchanged. An additional `$anchors` list element is attached:
#' \describe{
#'   \item{`n_anchors`}{Number of anchor pairs supplied.}
#'   \item{`rho_anchor`}{The value of `rho_anchor` used.}
#'   \item{`nodes_covered`}{Number of query nodes with ≥ 1 anchor pair.}
#'   \item{`coverage_fraction`}{`nodes_covered / nrow(query$codebook)`.}
#' }
#'
#' @return A `somalign_anchored_fit` object (also inherits `somalign_fit`).
#' @seealso [somalign_fit()] for the unanchored variant.
#' @examples
#' set.seed(1)
#' p   <- 3L
#' mat <- rbind(
#'   matrix(rnorm(20 * p, mean = -2), ncol = p),
#'   matrix(rnorm(20 * p, mean =  2), ncol = p)
#' )
#' colnames(mat) <- paste0("F", seq_len(p))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' shifted <- mat + 0.5
#' qry <- somalign_query(shifted, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' # Use 10 samples as anchors measured in both batches
#' anc_idx <- 1:10
#' fit <- somalign_fit_anchored(qry, ref,
#'                               anchor_old = mat[anc_idx, , drop = FALSE],
#'                               anchor_new = shifted[anc_idx, , drop = FALSE],
#'                               rho_anchor = 1)
#' fit$anchors
#' @export
somalign_fit_anchored <- function(query,
                                   reference,
                                   epsilon = 0.5,
                                   rho_query = 1,
                                   rho_ref = 1,
                                   solver = c("internal", "log_domain", "auto"),
                                   anchor_old,
                                   anchor_new,
                                   rho_anchor = 1.0,
                                   min_match_fraction = 0.05,
                                   confidence_threshold = 0.6,
                                   correction_min_mass = 1e-8,
                                   max_iter = 1000,
                                   tol = 1e-7,
                                   chunk_size = 10000L) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  solver <- match.arg(solver, c("internal", "log_domain", "auto"))

  if (!is.numeric(rho_anchor) || length(rho_anchor) != 1L ||
      !is.finite(rho_anchor) || rho_anchor < 0) {
    stop("`rho_anchor` must be a non-negative finite scalar.", call. = FALSE)
  }
  if (rho_anchor == 0) {
    message("`rho_anchor = 0`: anchor pairs have no effect. Use `somalign_fit()` for equivalent results.")
  }

  anchors_scaled <- .somalign_validate_anchors(anchor_old, anchor_new, reference)
  cost_bonus <- .somalign_anchor_cost_bonus(
    anchors_scaled$anchor_old_scaled,
    anchors_scaled$anchor_new_scaled,
    query$codebook,
    reference$codebook,
    rho_anchor,
    chunk_size
  )

  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol,
    cost_bonus = cost_bonus$bonus
  )

  label_transfer <- .somalign_transfer_labels(
    correspondence = transport$correspondence,
    label_prob = reference$label_prob,
    match_fraction = transport$match_fraction,
    min_match_fraction = min_match_fraction,
    confidence_threshold = confidence_threshold
  )
  node_shifts <- .somalign_node_shifts(
    query_codebook = query$codebook,
    reference_codebook = reference$codebook,
    correspondence = transport$correspondence,
    row_mass = transport$row_mass,
    match_fraction = transport$match_fraction,
    min_match_fraction = min_match_fraction,
    correction_min_mass = correction_min_mass
  )
  projection <- .somalign_project_pair(query, reference, node_shifts, chunk_size)
  diagnostics <- .somalign_build_diagnostics(
    transport, query, reference, node_shifts, projection, epsilon, rho_query, rho_ref
  )
  .somalign_fit_warnings(diagnostics)

  fit <- .somalign_new_fit(
    query, reference, transport, label_transfer, node_shifts, projection, diagnostics,
    anchors = list(
      n_anchors         = nrow(anchor_old),
      rho_anchor        = rho_anchor,
      nodes_covered     = cost_bonus$nodes_covered,
      coverage_fraction = cost_bonus$coverage_fraction
    )
  )
  class(fit) <- c("somalign_anchored_fit", "somalign_fit")
  fit
}

.somalign_validate_anchors <- function(anchor_old, anchor_new, reference) {
  anchor_old <- .somalign_as_matrix(anchor_old, what = "anchor_old")
  anchor_new <- .somalign_as_matrix(anchor_new, what = "anchor_new")

  if (nrow(anchor_old) == 0L)
    stop("`anchor_old` and `anchor_new` must have at least one row.", call. = FALSE)
  if (nrow(anchor_old) != nrow(anchor_new))
    stop("`anchor_old` and `anchor_new` must have the same number of rows.",
         call. = FALSE)

  features <- reference$features
  if (!is.null(features)) {
    anchor_old <- .somalign_select_features(anchor_old, features, what = "anchor_old")
    anchor_new <- .somalign_select_features(anchor_new, features, what = "anchor_new")
  } else {
    if (ncol(anchor_old) != length(reference$center))
      stop("`anchor_old` must have the same number of columns as reference features.",
           call. = FALSE)
    if (ncol(anchor_new) != length(reference$center))
      stop("`anchor_new` must have the same number of columns as reference features.",
           call. = FALSE)
  }

  .somalign_validate_finite(anchor_old, what = "anchor_old")
  .somalign_validate_finite(anchor_new, what = "anchor_new")

  anchor_old_scaled <- sweep(anchor_old, 2, reference$center, "-")
  anchor_old_scaled <- sweep(anchor_old_scaled, 2, reference$scale, "/")
  anchor_new_scaled <- sweep(anchor_new, 2, reference$center, "-")
  anchor_new_scaled <- sweep(anchor_new_scaled, 2, reference$scale, "/")

  list(anchor_old_scaled = anchor_old_scaled,
       anchor_new_scaled = anchor_new_scaled)
}

.somalign_anchor_cost_bonus <- function(anchor_old_scaled, anchor_new_scaled,
                                         query_codebook, reference_codebook,
                                         rho_anchor, chunk_size) {
  M <- nrow(query_codebook)
  K <- nrow(reference_codebook)
  n_anchors <- nrow(anchor_old_scaled)

  old_units <- .somalign_nearest_code_chunked(
    anchor_old_scaled, query_codebook, chunk_size = chunk_size
  )$unit
  new_units <- .somalign_nearest_code_chunked(
    anchor_new_scaled, reference_codebook, chunk_size = chunk_size
  )$unit

  # Build M × K count matrix (column-major linear indexing)
  lin_idx <- (new_units - 1L) * M + old_units
  A <- matrix(tabulate(lin_idx, nbins = M * K), nrow = M, ncol = K)

  nodes_covered    <- sum(rowSums(A) > 0L)
  coverage_fraction <- nodes_covered / M

  bonus <- rho_anchor * (A / n_anchors)

  list(
    bonus             = bonus,
    nodes_covered     = nodes_covered,
    coverage_fraction = coverage_fraction
  )
}
