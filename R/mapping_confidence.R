#' Per-cell reference-mapping confidence
#'
#' Scores how well each query cell sits inside the reference's covered region, as
#' a **relative proximity score** in `(0, 1]`. It is a Symphony-style
#' mapping-quality signal (Kang et al., 2021, \doi{10.1038/s41467-021-25957-x}): a
#' query cell well inside the reference structure scores near 1, a cell far from
#' any reference node scores near 0. It is a heuristic kernel weight, **not** a
#' calibrated probability, a coverage level, or a probability of correct
#' assignment; scores are relative within one reference and are not comparable
#' across references of different node density.
#'
#' This is distinct from, and complementary to, the boolean outside-reference flag
#' in [somalign_results()], which thresholds the distance to a single assigned
#' node against that node's distance quantile. Here the score is a smooth function
#' of the mean distance to the cell's `k` nearest reference nodes, normalised by
#' the reference map's own intrinsic node spacing, so it captures *local density*
#' of reference coverage rather than a per-node cutoff.
#'
#' Concretely, let \eqn{d_k} be a query cell's mean Euclidean distance (in
#' reference-scaled coordinates) to its `k` nearest reference codebook nodes, and
#' \eqn{s} the reference's intrinsic scale (the median over reference nodes of the
#' mean distance to their `k` nearest neighbouring nodes). The score is
#' \eqn{\exp(-(d_k / s)^2)}: it approaches 1 as a cell moves well inside the
#' reference (\eqn{d_k \ll s}), equals \eqn{e^{-1} \approx 0.37} at the reference's
#' typical node spacing (\eqn{d_k = s}), and decays toward 0 as the cell moves
#' further away.
#'
#' @param fit A `somalign_fit` (uses `fit$query$scaled_data` and
#'   `fit$reference$codebook`).
#' @param k Number of nearest reference nodes to average over. Default `10L`;
#'   clamped to the number of reference nodes.
#' @param chunk_size Rows of the query scaled-data matrix scored at a time.
#'   Default `10000L`.
#'
#' @return A numeric vector in `(0, 1]`, one per query cell, named by
#'   `fit$query$sample_id`, with attributes `k` and `reference_scale`.
#' @examples
#' set.seed(1)
#' x <- rbind(matrix(rnorm(40, -2), ncol = 4), matrix(rnorm(40, 2), ncol = 4))
#' colnames(x) <- paste0("m", 1:4)
#' # (illustrative; a real fit comes from somalign_fit())
#' @export
somalign_mapping_confidence <- function(fit, k = 10L, chunk_size = 10000L) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  X <- fit$query$scaled_data
  cb <- fit$reference$codebook
  if (is.null(X) || is.null(cb))
    stop("`fit` must carry query$scaled_data and reference$codebook.", call. = FALSE)
  .somalign_validate_finite(cb, "fit$reference$codebook")
  .somalign_validate_finite(X, "fit$query$scaled_data")
  n_nodes <- nrow(cb)
  k <- max(1L, min(as.integer(k), n_nodes))

  # Reference intrinsic scale: median over nodes of the mean distance to the k
  # nearest *other* nodes.
  dnn <- sqrt(.somalign_pairwise_distance(cb, cb))
  diag(dnn) <- Inf
  k_ref <- min(k, n_nodes - 1L)
  if (k_ref < 1L) {
    ref_scale <- 1
  } else {
    per_node <- apply(dnn, 1, function(d) mean(sort(d)[seq_len(k_ref)]))
    ref_scale <- stats::median(per_node)
  }
  if (!is.finite(ref_scale) || ref_scale <= 0) ref_scale <- 1

  # Per-cell mean distance to its k nearest reference nodes, chunked.
  n <- nrow(X)
  dk <- numeric(n)
  starts <- seq(1L, n, by = chunk_size)
  for (st in starts) {
    idx <- seq.int(st, min(st + chunk_size - 1L, n))
    d <- sqrt(.somalign_pairwise_distance(X[idx, , drop = FALSE], cb))  # |idx| x nodes
    dk[idx] <- apply(d, 1, function(row) mean(sort(row)[seq_len(k)]))
  }

  score <- exp(-(dk / ref_scale)^2)
  names(score) <- fit$query$sample_id
  attr(score, "k") <- k
  attr(score, "reference_scale") <- ref_scale
  score
}
