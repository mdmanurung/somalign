#' Mutual-nearest-neighbour novelty flag for query SOM nodes
#'
#' Computes a per-query-node "unmatched" flag based on mutual nearest-neighbour
#' (MNN) reciprocity between the query and reference codebooks.  This is a
#' **supplementary signal** intended for experimental evaluation (Experiment E2)
#' of whether MNN non-reciprocity adds novelty-detection value over the
#' unbalanced-OT leftover-mass signal (`match_fraction` / `diagnostics$nodes`).
#' It does **not** replace or modify the OT engine; the full Sinkhorn coupling,
#' soft-label posterior, and barycentric correction are unaffected.
#'
#' @section How reciprocity is computed:
#' Let `Q` be the query codebook and `R` the reference codebook, both in
#' **reference-scaled space** (the same coordinate system used to build the OT
#' cost matrix).  For each query node `q`, its nearest reference node is
#' `r(q) = argmin_j ||Q[q,] - R[j,]||`.  For each reference node `r`, its
#' nearest query node is `q(r) = argmin_i ||R[r,] - Q[i,]||`.  Query node `q`
#' is flagged `mnn_unmatched = TRUE` if and only if `q(r(q)) != q` — the
#' back-projection does not return to the originating query node.
#'
#' @section Batch-shift confound (critical for E2):
#' The flag is computed on the **pre-correction** (raw) query codebook.  Under a
#' population-specific batch shift, query nodes that represent populations already
#' present in the reference can break reciprocity purely because their codebook
#' vectors are displaced relative to the reference, even though no true novelty
#' exists.  This false-fire rate rises with shift magnitude and constitutes a
#' systematic confound that Experiment E2 is specifically designed to measure.
#' The OT-derived `match_fraction` is not confounded in the same way because the
#' Sinkhorn plan transports mass globally across the cost matrix.
#'
#' @section Feature weights:
#' When the OT fit used `feature_weights` (diagonal Mahalanobis cost), the OT
#' cost lives in `sqrt(w)`-scaled space, whereas this function always operates on
#' unweighted reference-scaled coordinates.  If you want exact alignment with the
#' OT cost space, apply `sqrt(feature_weights)` column scaling to both codebooks
#' before calling this function manually via `.somalign_nearest_code()`.  For the
#' E2 evaluation setting (unweighted default fits) this discrepancy is absent.
#'
#' @param fit A `somalign_fit` object.
#'
#' @return A named logical vector of length equal to the number of query SOM
#'   nodes.  `TRUE` indicates that node is MNN-unmatched with the reference
#'   (no reciprocal k=1 nearest neighbour) and is therefore a candidate novel
#'   node.  Names are `"q1"`, `"q2"`, ... matching query node indices.
#'
#' @seealso [somalign_mnn_novelty_cells()] to broadcast the per-node flag to
#'   individual cells; `fit$diagnostics$nodes$match_fraction` for the
#'   OT-derived leftover-mass signal.
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(100), nrow = 50, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' flags <- somalign_mnn_novelty(fit)
#' @export
somalign_mnn_novelty <- function(fit) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a `somalign_fit` object.", call. = FALSE)
  qcb <- fit$query$codebook
  rcb <- fit$reference$codebook
  if (is.null(qcb))
    stop("`fit$query$codebook` is NULL; cannot compute MNN novelty.", call. = FALSE)
  if (is.null(rcb))
    stop("`fit$reference$codebook` is NULL; cannot compute MNN novelty.", call. = FALSE)

  # k=1 nearest reference node for each query node
  r_of_q <- .somalign_nearest_code(qcb, rcb)$unit   # length n_query

  # k=1 nearest query node for each reference node
  q_of_r <- .somalign_nearest_code(rcb, qcb)$unit   # length n_ref

  # Reciprocity check: q is unmatched iff back-projecting r(q) does not return to q
  n_query <- nrow(qcb)
  mnn_unmatched <- q_of_r[r_of_q] != seq_len(n_query)

  names(mnn_unmatched) <- paste0("q", seq_len(n_query))
  mnn_unmatched
}


#' Broadcast the MNN novelty flag from query nodes to individual cells
#'
#' Maps the per-query-node `mnn_unmatched` vector returned by
#' [somalign_mnn_novelty()] to individual cells using each cell's query-SOM
#' node assignment (`fit$query$sample_unit`).
#'
#' @param fit A `somalign_fit` object.
#' @param mnn_unmatched A named logical vector of length `n_query_nodes`,
#'   as returned by [somalign_mnn_novelty()].
#'
#' @return A logical vector of length equal to the number of query cells.
#'   `TRUE` for cells assigned to an MNN-unmatched query node.
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(100), nrow = 50, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' flags <- somalign_mnn_novelty(fit)
#' cell_flags <- somalign_mnn_novelty_cells(fit, flags)
#' @export
somalign_mnn_novelty_cells <- function(fit, mnn_unmatched) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a `somalign_fit` object.", call. = FALSE)
  su <- fit$query$sample_unit
  if (is.null(su))
    stop("`fit$query$sample_unit` is NULL; cannot broadcast to cells.", call. = FALSE)
  n_query_nodes <- nrow(fit$query$codebook)
  if (!is.logical(mnn_unmatched) || length(mnn_unmatched) != n_query_nodes)
    stop("`mnn_unmatched` must be a logical vector of length equal to the number of query nodes (",
         n_query_nodes, ").", call. = FALSE)
  as.logical(mnn_unmatched[su])
}
