#' Soft (probabilistic) label projection for query cells
#'
#' Projects each query cell onto the reference by a Gaussian kernel over its k
#' nearest reference SOM nodes, and returns a per-cell probability distribution
#' over labels (or any node-level grouping). This is the soft analogue of the
#' hard nearest-node assignment behind `old_som_label` in [somalign_results()].
#'
#' Hard projection assigns each cell to a single nearest node and inherits that
#' node's label, discarding the cell's position within the node's Voronoi region.
#' At a label boundary this makes assignment a discontinuous 0/1 decision, so a
#' small batch shift can flip a cell's label and a boundary cell contributes a
#' full unit to one label with no hedging. Soft projection instead spreads each
#' cell over its nearest nodes' labels, which removes that boundary discontinuity
#' and reduces the quantisation variance in downstream per-sample frequency
#' estimates. It changes the *frequency estimate*, not the most-likely label:
#' `max.col()` of a soft-label matrix typically matches the hard label.
#'
#' @section Coverage contract:
#' Soft projection applies **no** acceptance or out-of-reference gating. Unlike
#' the hard path in [somalign_results()] (which flags cells via
#' `outside_reference_distance` / `transferred_label_accepted` and can suppress
#' low-confidence transfers), every cell contributes to the soft distribution
#' according to its distance to the reference nodes alone. Cells that fall
#' outside the reference are still projected onto their nearest nodes. If you
#' need to exclude out-of-reference cells, filter them with [somalign_results()]
#' before aggregating, or subset the query.
#'
#' @param fit A `somalign_fit` object.
#' @param node_groups Optional node-level grouping to project onto instead of the
#'   reference labels. Either a length-`n_nodes` vector (one group per reference
#'   node, e.g. a node-to-metacluster map; converted to indicators) or an
#'   `n_nodes` by `n_groups` matrix of node-group memberships. When `NULL`
#'   (default), `fit$reference$label_prob` is used, and the reference must carry
#'   labels.
#' @param k Integer. Number of nearest reference nodes used for the kernel,
#'   clamped to the number of reference nodes. Default `8L`.
#' @param bandwidth Positive scalar or `NULL`. Gaussian kernel bandwidth in
#'   reference-scaled space. `NULL` (default) uses the median nearest-neighbour
#'   distance of the reference codebook.
#' @param chunk_size Positive integer. Cells are processed in blocks of this size
#'   to bound peak memory. Default `10000L`.
#'
#' @return A numeric matrix of class `c("somalign_soft_labels", "matrix")`, one
#'   row per query cell and one column per label/group, with rows summing to 1
#'   (a row is all-zero when a cell's nearest nodes carry no label mass). Row
#'   names are the query sample identifiers; attributes `k` and `bandwidth`
#'   record the settings used.
#'
#' @seealso [somalign_soft_frequencies()], [somalign_results()]
#' @examples
#' if (requireNamespace("kohonen", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- rbind(matrix(rnorm(90 * 3, -3, 0.5), ncol = 3),
#'              matrix(rnorm(90 * 3,  3, 0.5), ncol = 3))
#'   colnames(x) <- paste0("m", seq_len(3))
#'   lab <- rep(c("low", "high"), each = 90)
#'   grid <- kohonen::somgrid(3, 3, "hexagonal")
#'   ref <- somalign_train_reference(x, labels = lab, grid = grid, rlen = 10)
#'   qry <- somalign_query(x, ref, grid = grid, rlen = 10)
#'   fit <- somalign_fit(qry, ref)
#'   soft <- somalign_soft_labels(fit)
#'   head(soft)
#' }
#' @export
somalign_soft_labels <- function(fit, node_groups = NULL, k = 8L,
                                 bandwidth = NULL, chunk_size = 10000L) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  .somalign_check_pos_int(k, "k")
  if (!is.null(bandwidth)) .somalign_check_pos_scalar(bandwidth, "bandwidth")
  .somalign_check_pos_int(chunk_size, "chunk_size")

  soft <- .somalign_soft_project(fit, node_groups, k, bandwidth, chunk_size,
                                 group = NULL)
  rownames(soft) <- fit$query$sample_id
  class(soft) <- c("somalign_soft_labels", "matrix")
  soft
}

#' Per-group soft label frequencies for query cells
#'
#' Aggregates the per-cell soft label distributions from [somalign_soft_labels()]
#' by a grouping (typically a biological sample), giving a group-by-label matrix
#' of soft frequencies. Soft aggregation reduces the boundary/quantisation
#' variance of hard per-sample cluster proportions, which improves the
#' reproducibility of cluster-abundance profiles across batches (for example the
#' centred-log-ratio abundance comparison in the label-transfer vignette).
#'
#' @param fit A `somalign_fit` object.
#' @param group Vector of length equal to the number of query cells, giving each
#'   cell's group (e.g. `sample_id` or `fcs_filename`).
#' @param node_groups Optional node-level grouping passed to
#'   [somalign_soft_labels()] (e.g. a node-to-metacluster map). Default `NULL`
#'   uses the reference labels.
#' @param k,bandwidth,chunk_size Passed to [somalign_soft_labels()].
#' @param normalize Logical. When `TRUE` (default) each group's row is divided by
#'   its total so rows are frequencies summing to 1; when `FALSE` the raw summed
#'   soft memberships (soft counts) are returned, suitable for count-based
#'   differential-abundance models. Note that a group's soft counts sum to the
#'   number of that group's cells whose neighbours carry label mass, not
#'   necessarily its total cell count: cells all of whose k nearest nodes are
#'   unlabelled contribute a zero row. With a fully labelled reference the two
#'   coincide.
#'
#' @return A numeric matrix of class
#'   `c("somalign_soft_frequencies", "matrix")`, one row per group and one column
#'   per label/group, with attributes `k`, `bandwidth`, and `normalized`.
#'
#' @seealso [somalign_soft_labels()], [somalign_results()]
#' @examples
#' if (requireNamespace("kohonen", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- rbind(matrix(rnorm(90 * 3, -3, 0.5), ncol = 3),
#'              matrix(rnorm(90 * 3,  3, 0.5), ncol = 3))
#'   colnames(x) <- paste0("m", seq_len(3))
#'   lab <- rep(c("low", "high"), each = 90)
#'   grid <- kohonen::somgrid(3, 3, "hexagonal")
#'   ref <- somalign_train_reference(x, labels = lab, grid = grid, rlen = 10)
#'   qry <- somalign_query(x, ref, grid = grid, rlen = 10)
#'   fit <- somalign_fit(qry, ref)
#'   sample_id <- rep(c("s1", "s2", "s3"), length.out = nrow(x))
#'   somalign_soft_frequencies(fit, sample_id)
#' }
#' @export
somalign_soft_frequencies <- function(fit, group, node_groups = NULL, k = 8L,
                                      bandwidth = NULL, normalize = TRUE,
                                      chunk_size = 10000L) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  .somalign_check_pos_int(k, "k")
  if (!is.null(bandwidth)) .somalign_check_pos_scalar(bandwidth, "bandwidth")
  .somalign_check_flag(normalize, "normalize")
  .somalign_check_pos_int(chunk_size, "chunk_size")
  if (length(group) != nrow(fit$query$scaled_data))
    stop("`group` must have one entry per query cell.", call. = FALSE)

  agg <- .somalign_soft_project(fit, node_groups, k, bandwidth, chunk_size,
                                group = group)
  # Arithmetic below strips custom attributes, so capture and re-attach them.
  k_eff <- attr(agg, "k")
  h <- attr(agg, "bandwidth")
  if (normalize) agg <- agg / pmax(rowSums(agg), .Machine$double.eps)
  attr(agg, "k") <- k_eff
  attr(agg, "bandwidth") <- h
  attr(agg, "normalized") <- normalize
  class(agg) <- c("somalign_soft_frequencies", "matrix")
  agg
}

# Resolve the node-level probability/indicator matrix to project onto: the
# reference label_prob by default, or a user node grouping (a per-node vector
# turned into indicator columns, or an already-formed node-by-group matrix).
.somalign_resolve_node_prob <- function(fit, node_groups) {
  n_nodes <- nrow(fit$reference$codebook)
  if (is.null(node_groups)) {
    lp <- fit$reference$label_prob
    if (is.null(lp) || ncol(lp) == 0L)
      stop("`fit` reference carries no labels; supply `node_groups` or build a ",
           "labelled reference.", call. = FALSE)
    return(lp)
  }
  if (is.matrix(node_groups)) {
    if (nrow(node_groups) != n_nodes)
      stop("`node_groups` matrix must have one row per reference node.", call. = FALSE)
    if (!is.numeric(node_groups) && !is.logical(node_groups))
      stop("`node_groups` matrix must be numeric or logical.", call. = FALSE)
    if (any(!is.finite(node_groups)))
      stop("`node_groups` matrix must not contain missing values.", call. = FALSE)
    if (any(node_groups < 0))
      stop("`node_groups` matrix must not contain negative memberships.", call. = FALSE)
    if (is.null(colnames(node_groups)))
      colnames(node_groups) <- paste0("group", seq_len(ncol(node_groups)))
    return(node_groups)
  }
  if (length(node_groups) != n_nodes)
    stop("`node_groups` must have one entry per reference node.", call. = FALSE)
  f <- factor(node_groups)
  m <- matrix(0, n_nodes, nlevels(f), dimnames = list(NULL, levels(f)))
  ok <- !is.na(f)
  m[cbind(which(ok), as.integer(f[ok]))] <- 1
  m
}

# Resolve the node-level target, row-normalise it, and project query cells onto
# it with the shared fused k-NN kernel smoother (.somalign_knn_smooth, R/utils.R).
# `group = NULL` returns a per-cell N x C matrix; `group` supplied returns
# per-group soft counts (groups x C) without ever materialising the N x C matrix.
# Nodes with no label mass are gated out; cells whose k neighbours are all
# unlabelled get a zero row. Carries `k` and `bandwidth` as attributes.
.somalign_soft_project <- function(fit, node_groups, k, bandwidth, chunk_size,
                                   group) {
  node_prob <- .somalign_resolve_node_prob(fit, node_groups)
  mass <- rowSums(node_prob)
  ok <- mass > 0
  node_prob[ok, ] <- node_prob[ok, , drop = FALSE] / mass[ok]
  node_prob[!ok, ] <- 0
  cb <- fit$reference$codebook
  k_eff <- min(as.integer(k), nrow(cb))
  h <- if (is.null(bandwidth)) .somalign_default_bandwidth(cb) else bandwidth
  out <- .somalign_knn_smooth(fit$query$scaled_data, cb, node_prob,
                              as.numeric(ok), h, k_eff, chunk_size, group = group)
  attr(out, "k") <- k_eff
  attr(out, "bandwidth") <- h
  out
}
