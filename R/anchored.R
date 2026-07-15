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
#' @param anchor_old Numeric matrix (n_anchors × p). Old-batch measurements
#'   of the anchor samples. Must be **raw (un-normalized) values in the same
#'   units and preprocessing pipeline as the data used to train `reference`**.
#'   Do not pre-center or pre-scale; this function applies `reference$center`
#'   and `reference$scale` internally. Also accepts a data frame of numeric
#'   columns.
#' @param anchor_new Numeric matrix (n_anchors × p). New-batch measurements
#'   of the **same** anchor samples. Must be **raw (un-normalized) values in
#'   the same units and preprocessing pipeline as `anchor_old`**. Do not
#'   pre-center or pre-scale; this function applies `reference$center` and
#'   `reference$scale` internally. Rows of `anchor_old` and `anchor_new` must
#'   correspond to the same biological units. Also accepts a data frame of
#'   numeric columns.
#' @param rho_anchor Non-negative scalar. Controls how strongly anchor pairs
#'   bias the OT cost. At `rho_anchor = 0` the result equals [somalign_fit()].
#'   Larger values reduce the effective cost for anchor-supported node pairs,
#'   concentrating the transport plan on those routes. Typical range: 0.5--3.
#'   Has no effect when `correction = "subspace"`.
#' @param epsilon Entropic regularisation strength (see [somalign_fit()]).
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver Sinkhorn solver variant. See [somalign_fit()].
#' @param min_match_fraction Minimum transported fraction for label transfer.
#' @param confidence_threshold Minimum top-label probability for label transfer.
#' @param correction_min_mass Minimum transported mass for a node correction.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Sinkhorn convergence tolerance.
#' @param chunk_size Integer. Samples projected per chunk. Default `10000L`.
#' @param correction Character. Correction strategy — one of
#'   `"cost_bonus"` (default), `"subspace"`, or `"both"`. See Details.
#' @param variance_threshold Numeric in (0, 1]. Cumulative singular-value-squared
#'   fraction for selecting the rank of the batch subspace. Default `0.9`
#'   (CellANOVA convention). Only used when `correction` is `"subspace"` or
#'   `"both"`.
#'
#' @details
#' **Correction modes.** Three strategies are available via the `correction`
#' argument.
#'
#' - `"cost_bonus"` (default, current behaviour): the anchor count matrix
#'   biases the OT cost so anchor-supported node pairs are cheaper; the
#'   resulting node shifts are applied to the full feature space.
#'
#' - `"subspace"`: a batch subspace \eqn{V_{\text{batch}}} is estimated by
#'   SVD of the anchor displacement matrix
#'   \eqn{D = X_{\text{old}} - X_{\text{new}}} (n_anchors × p). Because each
#'   row of \eqn{D} is a *same-biological-unit* before–after measurement, the
#'   dominant singular vectors isolate the true batch direction. Node shifts
#'   from a *plain* OT solve (no cost bonus) are then projected onto
#'   \eqn{V_{\text{batch}}}: only the batch-direction component is applied.
#'   Biological variation orthogonal to \eqn{V_{\text{batch}}} is preserved.
#'   A synthetic validation shows the orthogonal component survives at
#'   ~99.7% (1.496 vs ideal 1.500) while `"cost_bonus"` erases it.
#'   The rank \eqn{r} is the smallest index where the cumulative squared
#'   singular values reach `variance_threshold` (default 0.9).
#'   \eqn{D} is **not centred** — the mean batch direction is the dominant
#'   structure we want to capture.
#'
#' - `"both"`: applies the cost bonus to the OT solve *and* restricts the
#'   resulting shifts to \eqn{V_{\text{batch}}}.
#'
#' `"subspace"` and `"both"` expose `fit$anchors$batch_subspace` (a list with
#' `V`, `rank`, `variance_explained`). `"cost_bonus"` sets this to `NULL`.
#'
#' **Cost modification.** Let \eqn{C} be the M×K codebook distance matrix
#' normalised by its median positive entry (as in [somalign_fit()]). The
#' anchor pairs are projected onto the query codebook (old batch) and
#' reference codebook (new batch), yielding a count matrix \eqn{A} where
#' \eqn{A_{kl}} is the number of anchor pairs whose old measurement maps to
#' query node \eqn{k} and new measurement maps to reference node \eqn{l}.
#' (The query SOM was trained on new-batch data, so projecting the old-batch
#' anchor onto it identifies which query node the anchor occupied before the
#' batch shift; projecting the new-batch anchor onto the reference SOM
#' identifies the corresponding reference node after the shift.)
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
#'   \item{`correction`}{The correction mode: `"cost_bonus"`, `"subspace"`, or `"both"`.}
#'   \item{`nodes_covered`}{Number of query nodes with ≥ 1 anchor pair.}
#'   \item{`coverage_fraction`}{`nodes_covered / nrow(query$codebook)`.}
#'   \item{`batch_subspace`}{For `"subspace"` and `"both"` modes: a list with
#'     `V` (p × rank matrix), `rank` (integer), and `variance_explained`
#'     (cumulative variance at the selected rank). `NULL` for `"cost_bonus"`.}
#' }
#'
#' @return A `somalign_anchored_fit` object (also inherits `somalign_fit`).
#' @note At small `epsilon` with high anchor coverage the anchor bonus zeros out
#'   many entries of the normalised cost matrix, which sharpens the Sinkhorn
#'   kernel and can drive the remaining entries toward numerical underflow. If
#'   the solver warns about kernel underflow, pass `solver = "log_domain"`, which
#'   works in log-potential space and avoids the issue.
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
                                   anchor_old,
                                   anchor_new,
                                   rho_anchor = 1.0,
                                   epsilon = 0.1,
                                   rho_query = 1,
                                   rho_ref = 1,
                                   solver = c("internal", "log_domain", "auto"),
                                   min_match_fraction = 0.05,
                                   confidence_threshold = 0.6,
                                   correction_min_mass = 1e-8,
                                   max_iter = 1000,
                                   tol = 1e-7,
                                   chunk_size = 10000L,
                                   correction = c("cost_bonus", "subspace", "both"),
                                   variance_threshold = 0.9) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  solver <- match.arg(solver, c("internal", "log_domain", "auto"))
  correction <- match.arg(correction, c("cost_bonus", "subspace", "both"))
  if (!is.numeric(rho_anchor) || length(rho_anchor) != 1L ||
      !is.finite(rho_anchor) || rho_anchor < 0) {
    stop("`rho_anchor` must be a non-negative finite scalar.", call. = FALSE)
  }
  .somalign_check_pos_scalar(epsilon, "epsilon")
  .somalign_check_fit_params(rho_query, rho_ref, min_match_fraction,
                             confidence_threshold, correction_min_mass,
                             max_iter, tol, chunk_size)
  .somalign_check_unit_scalar(variance_threshold, "variance_threshold")
  anchors_scaled <- .somalign_validate_anchors(anchor_old, anchor_new, reference)
  .somalign_anchored_dispatch(
    query, reference, anchors_scaled, rho_anchor, epsilon, rho_query, rho_ref,
    solver, min_match_fraction, confidence_threshold, correction_min_mass,
    chunk_size, max_iter, tol, correction, variance_threshold
  )
}

.somalign_anchored_dispatch <- function(query, reference, anchors_scaled,
                                         rho_anchor, epsilon, rho_query, rho_ref,
                                         solver, min_match_fraction,
                                         confidence_threshold, correction_min_mass,
                                         chunk_size, max_iter, tol,
                                         correction, variance_threshold) {
  use_bonus    <- correction %in% c("cost_bonus", "both") && rho_anchor > 0
  use_subspace <- correction %in% c("subspace", "both")
  if (rho_anchor == 0 && correction == "cost_bonus") {
    message("`rho_anchor = 0`: anchor pairs have no effect. Use `somalign_fit()` for equivalent results.")
  }
  cb <- if (use_bonus) {
    .somalign_anchor_cost_bonus(anchors_scaled$anchor_old_scaled,
                                anchors_scaled$anchor_new_scaled,
                                query$codebook, reference$codebook,
                                rho_anchor, chunk_size)
  } else { list(bonus = NULL, nodes_covered = 0L, coverage_fraction = 0) }
  batch_sub <- if (use_subspace) {
    .somalign_batch_subspace(anchors_scaled$anchor_old_scaled,
                             anchors_scaled$anchor_new_scaled,
                             variance_threshold)
  } else { NULL }
  shift_fn <- if (use_subspace) {
    V <- batch_sub$V
    function(s) s %*% V %*% t(V)
  } else { NULL }
  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol,
    cost_bonus = cb$bonus
  )
  fit <- .somalign_finish_fit(
    query, reference, transport, min_match_fraction, confidence_threshold,
    correction_min_mass, chunk_size, epsilon, rho_query, rho_ref,
    shift_transform = shift_fn,
    anchors = list(
      n_anchors         = nrow(anchors_scaled$anchor_old_scaled),
      rho_anchor        = rho_anchor,
      correction        = correction,
      nodes_covered     = cb$nodes_covered,
      coverage_fraction = cb$coverage_fraction,
      batch_subspace    = batch_sub
    )
  )
  class(fit) <- c("somalign_anchored_fit", "somalign_fit")
  fit
}

.somalign_batch_subspace <- function(anchor_old_scaled, anchor_new_scaled,
                                      variance_threshold) {
  .somalign_subspace_svd(anchor_old_scaled - anchor_new_scaled, variance_threshold)
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
    ref_names <- names(reference$center)
    if (!is.null(ref_names)) {
      anchor_old <- .somalign_select_features(anchor_old, ref_names, what = "anchor_old")
      anchor_new <- .somalign_select_features(anchor_new, ref_names, what = "anchor_new")
    } else {
      if (ncol(anchor_old) != length(reference$center))
        stop("`anchor_old` must have the same number of columns as reference features.",
             call. = FALSE)
      if (ncol(anchor_new) != length(reference$center))
        stop("`anchor_new` must have the same number of columns as reference features.",
             call. = FALSE)
    }
  }

  .somalign_validate_finite(anchor_old, what = "anchor_old")
  .somalign_validate_finite(anchor_new, what = "anchor_new")

  list(
    anchor_old_scaled = .somalign_scale_matrix(anchor_old, reference$center, reference$scale),
    anchor_new_scaled = .somalign_scale_matrix(anchor_new, reference$center, reference$scale)
  )
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
