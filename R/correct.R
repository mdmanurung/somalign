#' Batch-correct query marker expression for downstream analysis
#'
#' Returns a cell-level (cells by markers) batch-corrected marker expression
#' matrix for the query cells, intended for downstream visualisation and
#' differential expression. The correction is restricted to the anchor-estimated
#' batch subspace and smoothed across each cell's nearest self-organising map
#' (SOM) nodes, so variation orthogonal to the batch direction is preserved.
#'
#' Unlike the per-node correction used internally (which is piecewise constant
#' across a node and contracts populations toward one another), this function
#' interpolates a smooth per-cell shift from the shifts of the k nearest SOM
#' nodes and confines it to the batch subspace. Cells therefore receive a
#' continuous correction, and structure orthogonal to the batch subspace is left
#' intact.
#'
#' @section Scope and limitations:
#' This output is an auxiliary correction aid, not the primary product of
#' `somalign`, which is label transfer (see [somalign_results()]). The
#' correction is restricted to the batch subspace and preserves orthogonal
#' variation, but within that subspace it still *reduces*, rather than fully
#' removes, the distance between populations; it does not undo genuine
#' over-merging. For comparing cell-type composition or abundance across
#' batches, use the direct projection columns from [somalign_results()] together
#' with a compositional (centred log-ratio) transform, not corrected expression.
#' Run [somalign_topology_audit()] before relying on this output to confirm that
#' correction is warranted for your data.
#'
#' @section Subspace restriction:
#' The correction is confined to the span of the batch subspace \eqn{V}
#' estimated from anchor displacements during fitting. Each cell shift is a
#' weighted average of node shifts that already lie in that span, so the cell
#' shift lies in it too; no post-smoothing re-projection is applied. Variation
#' orthogonal to \eqn{V} is untouched.
#'
#' @section Future extension:
#' A contraction-free variant would estimate the batch-shift field directly from
#' anchor displacements by kernel regression, rather than from the barycentric
#' node shifts. That path needs the scaled anchor positions, which a fit does not
#' store (only the anchor displacement matrix is retained), so it is not
#' available here.
#'
#' @param fit A `somalign_fit` from [somalign_fit_anchored()] with
#'   `correction = "subspace"` or `"both"`, or from [somalign_fit_two_pass()].
#'   A plain [somalign_fit()] or a `correction = "cost_bonus"` anchored fit
#'   carries no batch subspace and is rejected.
#' @param units One of `"raw"` (original expression units, the default) or
#'   `"scaled"` (reference-scaled units).
#' @param smooth Logical. When `TRUE` (default), smooths the correction across
#'   the k nearest SOM nodes with a Gaussian kernel. When `FALSE`, each cell
#'   takes its nearest node's shift directly (piecewise constant); this
#'   reproduces the correction `somalign` uses internally and is provided as a
#'   diagnostic baseline, not recommended for downstream analysis.
#' @param k Integer. Number of nearest SOM nodes used for smoothing, clamped to
#'   the number of query SOM nodes. Default `8L`.
#' @param bandwidth Positive scalar or `NULL`. Gaussian kernel bandwidth in
#'   reference-scaled space. `NULL` (default) uses the median nearest-neighbour
#'   distance of the SOM codebook, which adapts to the lattice spacing.
#' @param confidence_gate Logical. When `TRUE` (default), each node's kernel
#'   weight is multiplied by its transported match fraction, down-weighting nodes
#'   the transport plan could not align. Nodes whose correction is disallowed
#'   contribute zero weight in either case.
#' @param chunk_size Positive integer. Cells are processed in blocks of this size
#'   to bound peak memory. Default `10000L`.
#'
#' @return A numeric matrix of class
#'   `c("somalign_corrected_expression", "matrix")`, with one row per query cell
#'   and one column per marker. Row names are the query sample identifiers and
#'   column names are the reference features. Attributes `units`, `bandwidth`,
#'   `smooth`, and `k` record the settings used.
#'
#' @seealso [somalign_results()], [somalign_topology_audit()],
#'   [somalign_fit_anchored()], [somalign_fit_two_pass()]
#'
#' @examples
#' if (requireNamespace("kohonen", quietly = TRUE)) {
#'   set.seed(1)
#'   ref_x <- matrix(rnorm(60 * 3, 0, 0.5), ncol = 3,
#'                   dimnames = list(NULL, paste0("m", 1:3)))
#'   grid <- kohonen::somgrid(2, 2, "hexagonal")
#'   ref <- somalign_train_reference(ref_x, grid = grid, rlen = 10)
#'   shift <- matrix(c(2, 0, 0), nrow(ref_x), 3, byrow = TRUE)
#'   qry <- somalign_query(ref_x + shift, ref, grid = grid, rlen = 10)
#'   anc <- ref_x[1:20, ]
#'   fit <- somalign_fit_anchored(qry, ref, anchor_old = anc,
#'                                anchor_new = anc + shift[1:20, ],
#'                                correction = "subspace")
#'   expr <- somalign_correct_expression(fit)
#'   dim(expr)
#' }
#' @export
somalign_correct_expression <- function(fit,
                                        units = c("raw", "scaled"),
                                        smooth = TRUE,
                                        k = 8L,
                                        bandwidth = NULL,
                                        confidence_gate = TRUE,
                                        chunk_size = 10000L) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  units <- match.arg(units)
  .somalign_check_flag(smooth, "smooth")
  .somalign_check_pos_int(k, "k")
  if (!is.null(bandwidth)) .somalign_check_pos_scalar(bandwidth, "bandwidth")
  .somalign_check_flag(confidence_gate, "confidence_gate")
  .somalign_check_pos_int(chunk_size, "chunk_size")

  bsub <- .somalign_get_batch_subspace(fit)
  if (is.null(bsub))
    stop("`fit` does not contain a batch subspace. Re-fit with ",
         "`somalign_fit_anchored(correction = \"subspace\")` or ",
         "`somalign_fit_two_pass()`, and run `somalign_topology_audit()` first ",
         "to check that correction is warranted.", call. = FALSE)
  if (all(bsub$V == 0))
    warning("The batch subspace in `fit` has zero variance; all correction ",
            "shifts are zero.", call. = FALSE)

  shifts <- .somalign_correction_shifts(
    fit, smooth = smooth, k = k, bandwidth = bandwidth,
    confidence_gate = confidence_gate, chunk_size = chunk_size)

  corrected <- fit$query$scaled_data + shifts
  if (identical(units, "raw"))
    corrected <- sweep(sweep(corrected, 2, fit$reference$scale, "*"),
                       2, fit$reference$center, "+")

  .somalign_new_corrected_expression(
    corrected, fit, units, attr(shifts, "bandwidth"), smooth, attr(shifts, "k"))
}

# Resolve the batch subspace regardless of which fitting path produced it.
.somalign_get_batch_subspace <- function(fit) {
  if (!is.null(fit$anchors$batch_subspace)) return(fit$anchors$batch_subspace)
  if (!is.null(fit$two_pass$batch_subspace)) return(fit$two_pass$batch_subspace)
  NULL
}

# Per-cell correction shifts (N x p), in reference-scaled space. Dispatches
# between the smoothed field and the piecewise-constant nearest-node baseline,
# and records the bandwidth and effective k as attributes.
.somalign_correction_shifts <- function(fit, smooth, k, bandwidth,
                                        confidence_gate, chunk_size) {
  node_shifts <- fit$node_shifts
  allowed <- attr(node_shifts, "correction_allowed")
  if (is.null(allowed)) allowed <- rep(TRUE, nrow(node_shifts))
  if (!smooth) {
    su <- fit$query$sample_unit
    shifts <- node_shifts[su, , drop = FALSE]
    shifts[!allowed[su], ] <- 0
    attr(shifts, "bandwidth") <- NA_real_
    attr(shifts, "k") <- 1L
    return(shifts)
  }
  h <- if (is.null(bandwidth))
    .somalign_default_bandwidth(fit$query$codebook) else bandwidth
  k_eff <- min(as.integer(k), nrow(fit$query$codebook))
  shifts <- .somalign_smooth_cell_shifts(
    fit$query$scaled_data, fit$query$codebook, node_shifts, allowed,
    fit$diagnostics$ot$match_fraction, k_eff, h, confidence_gate, chunk_size)
  attr(shifts, "bandwidth") <- h
  attr(shifts, "k") <- k_eff
  shifts
}

# Median nearest-neighbour distance of the SOM codebook. Sets the default
# kernel bandwidth so it adapts to the lattice spacing in scaled space.
.somalign_default_bandwidth <- function(codebook) {
  if (nrow(codebook) < 2L) return(1)
  d2 <- .somalign_pairwise_distance(codebook, codebook)
  diag(d2) <- Inf
  h <- sqrt(stats::median(apply(d2, 1, min)))
  if (!is.finite(h) || h <= 0) 1 else h
}

# k nearest SOM nodes per cell, chunked over cells. Returns the node indices
# (N x k, ascending distance) and their squared distances (N x k).
.somalign_knn_codes_chunked <- function(scaled_data, codebook, k, chunk_size) {
  x <- as.matrix(scaled_data)
  n <- nrow(x)
  k <- min(as.integer(k), nrow(codebook))
  indices <- matrix(0L, n, k)
  sq_dist <- matrix(0, n, k)
  cs <- if (is.null(chunk_size) || is.infinite(chunk_size)) n else as.integer(chunk_size)
  cs <- max(1L, cs)
  for (s in seq(1L, n, by = cs)) {
    idx <- s:min(s + cs - 1L, n)
    d2 <- .somalign_pairwise_distance(x[idx, , drop = FALSE], codebook)
    ord <- t(apply(d2, 1L, order))[, seq_len(k), drop = FALSE]
    rows <- seq_len(nrow(d2))
    indices[idx, ] <- ord
    sq_dist[idx, ] <- d2[cbind(rep(rows, k), as.vector(ord))]
  }
  list(indices = indices, sq_dist = sq_dist)
}

# Stabilised Gaussian kernel weights times per-neighbour gate weights. The
# per-row minimum distance is subtracted before exponentiating, so the nearest
# node always has raw weight 1 and rows never underflow to all zeros.
.somalign_kernel_weights <- function(sq_dist, bandwidth, gate) {
  d2min <- apply(sq_dist, 1L, min)
  w <- exp(-(sq_dist - d2min) / (2 * bandwidth^2))
  w * gate
}

# Smooth per-cell correction field: kernel-weighted, confidence-gated average
# of the k nearest node shifts, computed in blocks to bound memory. Nodes whose
# correction is disallowed (or, with gating, poorly matched) contribute zero.
.somalign_smooth_cell_shifts <- function(scaled_data, codebook, node_shifts,
                                         allowed, match_fraction, k, bandwidth,
                                         confidence_gate, chunk_size) {
  n <- nrow(scaled_data)
  p <- ncol(node_shifts)
  node_conf <- if (confidence_gate && !is.null(match_fraction)) {
    mf <- match_fraction
    mf[!is.finite(mf)] <- 0
    mf
  } else {
    rep(1, nrow(node_shifts))
  }
  node_conf[!allowed] <- 0
  knn <- .somalign_knn_codes_chunked(scaled_data, codebook, k, chunk_size)
  gate <- matrix(node_conf[knn$indices], n, ncol(knn$indices))
  w <- .somalign_kernel_weights(knn$sq_dist, bandwidth, gate)
  wsum <- rowSums(w)
  shifts <- matrix(0, n, p)
  cs <- if (is.null(chunk_size) || is.infinite(chunk_size)) n else max(1L, as.integer(chunk_size))
  for (s in seq(1L, n, by = cs)) {
    idx <- s:min(s + cs - 1L, n)
    acc <- matrix(0, length(idx), p)
    for (j in seq_len(ncol(knn$indices)))
      acc <- acc + w[idx, j] * node_shifts[knn$indices[idx, j], , drop = FALSE]
    keep <- wsum[idx] > 0
    acc[keep, ] <- acc[keep, , drop = FALSE] / wsum[idx][keep]
    acc[!keep, ] <- 0
    shifts[idx, ] <- acc
  }
  shifts
}

# Attach dimnames, settings attributes, and class to the corrected matrix.
.somalign_new_corrected_expression <- function(x, fit, units, bandwidth, smooth, k) {
  rownames(x) <- fit$query$sample_id
  colnames(x) <- colnames(fit$query$scaled_data)
  attr(x, "units") <- units
  attr(x, "bandwidth") <- bandwidth
  attr(x, "smooth") <- smooth
  attr(x, "k") <- k
  class(x) <- c("somalign_corrected_expression", "matrix")
  x
}
