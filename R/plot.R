utils::globalVariables(c(
  "query_mass", "transported_mass", "match_fraction",
  "query_node", "correction_norm", "correction_allowed",
  "Projection", "outside_pct",
  "old", "transferred", "pct",
  "src", "xmin", "xend", "xmax", "yend", "flag",
  "feature", "value"
))

# ---- internal helpers -------------------------------------------------------

.somalign_downsample_rows <- function(x, n, seed = 1L) {
  if (nrow(x) <= n) return(x)
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  x[sample(nrow(x), n), , drop = FALSE]
}

.somalign_melt_long <- function(mat) {
  nr  <- nrow(mat)
  nms <- colnames(mat)
  if (is.null(nms)) nms <- as.character(seq_len(ncol(mat)))
  data.frame(
    feature = factor(rep(nms, each = nr), levels = nms),
    value   = as.vector(mat),
    stringsAsFactors = FALSE
  )
}

.somalign_node_df <- function(fit) {
  node_df      <- as.data.frame(somalign_diagnostics(fit)$nodes)
  ref_top      <- .somalign_reference_top_labels(fit$reference)
  top_ref      <- max.col(fit$correspondence, ties.method = "first")
  node_df$top_ref_label <- ref_top$label[top_ref]
  node_df
}

.somalign_confusion_df <- function(results) {
  d <- results[
    !is.na(results$old_som_label) &
    !is.na(results$transferred_label) &
    results$transferred_label_accepted, ]
  if (nrow(d) == 0L)
    return(data.frame(old = character(), transferred = character(), pct = numeric()))
  tab  <- as.data.frame(table(old = d$old_som_label, transferred = d$transferred_label))
  tots <- stats::aggregate(Freq ~ old, data = tab, FUN = sum)
  names(tots)[2L] <- "total"
  tab  <- merge(tab, tots, by = "old")
  tab$pct <- tab$Freq / tab$total
  tab
}

.somalign_check_scalar <- function(x, nm, lo = -Inf, hi = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
    stop(sprintf("`%s` must be a single finite number.", nm), call. = FALSE)
  if (x < lo || x > hi)
    stop(sprintf("`%s` must be in [%g, %g].", nm, lo, hi), call. = FALSE)
  invisible(x)
}

.somalign_dist_data <- function(query, reference, reference_data, features, downsample, seed) {
  feats <- if (is.null(features)) query$reference_features else {
    miss <- setdiff(features, query$reference_features)
    if (length(miss))
      stop("Unknown features: ", paste(miss, collapse = ", "), call. = FALSE)
    features
  }
  if (!is.null(reference_data)) {
    if (!is.matrix(reference_data) && !is.data.frame(reference_data))
      stop("`reference_data` must be a numeric matrix or data frame.", call. = FALSE)
    miss_r <- setdiff(feats, colnames(reference_data))
    if (length(miss_r))
      stop("`reference_data` is missing columns: ", paste(miss_r, collapse = ", "),
           call. = FALSE)
  }
  q_sub  <- .somalign_downsample_rows(query$scaled_data[, feats, drop = FALSE], downsample, seed)
  long_q <- .somalign_melt_long(q_sub)
  long_q$source <- "Query (cells)"
  if (!is.null(reference_data)) {
    r_sub    <- .somalign_downsample_rows(reference_data[, feats, drop = FALSE], downsample, seed)
    long_ref <- .somalign_melt_long(r_sub)
    long_ref$source <- "Reference (cells)"
    list(query = long_q, ref_cells = long_ref, ref_cb = NULL)
  } else {
    ref_cb  <- if (!is.null(reference)) reference$codebook[, feats, drop = FALSE] else NULL
    long_cb <- if (!is.null(ref_cb)) .somalign_melt_long(ref_cb) else NULL
    list(query = long_q, ref_cells = NULL, ref_cb = long_cb)
  }
}

# ---- after-projection: plot functions ---------------------------------------

#' Plot node mass balance
#'
#' Scatter plot of query node mass vs transported mass, coloured by match
#' fraction. Points lying on the diagonal received all their mass; points
#' below it had mass destroyed by the unbalanced OT solver.
#'
#' @param fit A `somalign_fit` object.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_plot_mass_balance(fit)
#' @export
somalign_plot_mass_balance <- function(fit) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  nd <- .somalign_node_df(fit)
  ggplot2::ggplot(nd, ggplot2::aes(x = query_mass, y = transported_mass,
                                    colour = match_fraction)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_viridis_c(
      "Match\nfraction", option = "plasma", limits = c(0, 1)) +
    ggplot2::labs(title = "Node mass balance",
                  x = "Query node mass", y = "Transported mass") +
    ggplot2::theme_minimal()
}

#' Plot per-node match fraction
#'
#' Sorted bar chart of the match fraction for each query SOM node.
#' Nodes below `threshold` received too little mass from the OT plan and
#' are the primary candidates for inspection.
#'
#' @param fit A `somalign_fit` object.
#' @param threshold Numeric scalar. Threshold line drawn on the plot.
#'   Default `0.05`.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_plot_match_fraction(fit)
#' @export
somalign_plot_match_fraction <- function(fit, threshold = 0.05) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  .somalign_check_scalar(threshold, "threshold", lo = 0, hi = 1)
  nd <- .somalign_node_df(fit)
  ggplot2::ggplot(nd, ggplot2::aes(
    x = reorder(factor(query_node), match_fraction),
    y = match_fraction, fill = match_fraction)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.2) +
    ggplot2::geom_hline(
      yintercept = threshold, colour = "#d73027", linetype = "dashed") +
    ggplot2::scale_fill_viridis_c(
      option = "plasma", limits = c(0, 1), guide = "none") +
    ggplot2::labs(title = "Match fraction per query node",
                  subtitle = sprintf("Dashed = %.2f threshold", threshold),
                  x = "Query node (sorted)", y = "Match fraction") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

#' Plot per-node correction norms
#'
#' Bar chart of the correction vector length per query SOM node, coloured by
#' whether the correction was applied (green) or suppressed due to low mass
#' or match fraction (red).
#'
#' @param fit A `somalign_fit` object.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_plot_correction(fit)
#' @export
somalign_plot_correction <- function(fit) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  nd <- .somalign_node_df(fit)
  ggplot2::ggplot(nd, ggplot2::aes(
    x = reorder(factor(query_node), -correction_norm),
    y = correction_norm, fill = correction_allowed)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#1a9641", "FALSE" = "#d73027"),
      labels = c("TRUE" = "Applied", "FALSE" = "Suppressed")) +
    ggplot2::labs(title = "Correction norm per query node",
                  subtitle = "Green = applied; red = suppressed (too little mass or match fraction)",
                  x = "Query node (sorted by correction)", y = "Correction norm",
                  fill = "Correction") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}

#' Plot fraction of cells outside reference thresholds
#'
#' Compares the percentage of cells that fall outside the reference distance
#' threshold before (Direct) and after (Corrected) applying the OT correction
#' vectors.
#'
#' @param fit A `somalign_fit` object.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_plot_outside_fraction(fit)
#' @export
somalign_plot_outside_fraction <- function(fit) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  diag <- somalign_diagnostics(fit)
  proj_df <- data.frame(
    Projection = factor(
      c("Direct", "Corrected"), levels = c("Direct", "Corrected")),
    outside_pct = c(
      100 * diag$projection$outside_direct_fraction,
      100 * diag$projection$outside_corrected_fraction)
  )
  ggplot2::ggplot(proj_df,
                  ggplot2::aes(x = Projection, y = outside_pct, fill = Projection)) +
    ggplot2::geom_col(width = 0.5, colour = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", outside_pct)),
      vjust = -0.4, size = 4) +
    ggplot2::scale_fill_manual(
      values = c("Direct" = "#d73027", "Corrected" = "#1a9641"),
      guide  = "none") +
    ggplot2::labs(title = "% cells outside reference thresholds",
                  subtitle = "Corrected should be lower (or equal) to Direct",
                  y = "% outside threshold", x = NULL) +
    ggplot2::ylim(0, max(proj_df$outside_pct) * 1.2 + 5) +
    ggplot2::theme_minimal()
}

#' Plot label transfer confusion heatmap
#'
#' Row-normalised heatmap of old-to-transferred label pairs for accepted
#' cells. High values on the diagonal indicate coherent transfer; strong
#' off-diagonal entries warrant further inspection.
#'
#' @param fit A `somalign_fit` object.
#' @param min_confidence Minimum `transferred_label_confidence` to include.
#'   `NULL` (default) imposes no additional filter beyond acceptance.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_plot_label_confusion(fit)
#' @export
somalign_plot_label_confusion <- function(fit, min_confidence = NULL) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  if (!is.null(min_confidence))
    .somalign_check_scalar(min_confidence, "min_confidence", lo = 0, hi = 1)
  results <- somalign_results(fit)
  if (!is.null(min_confidence))
    results <- results[
      is.na(results$transferred_label_confidence) |
      results$transferred_label_confidence >= min_confidence, ]
  conf <- .somalign_confusion_df(results)
  if (nrow(conf) == 0L)
    stop("No accepted transferred labels found; cannot build confusion plot.",
         call. = FALSE)
  ggplot2::ggplot(conf, ggplot2::aes(x = transferred, y = old, fill = pct)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%", 100 * pct)),
      size = 4, fontface = "bold") +
    ggplot2::scale_fill_viridis_c(
      "Fraction", option = "plasma", limits = c(0, 1)) +
    ggplot2::labs(title = "Label transfer confusion (accepted cells)",
                  subtitle = "Row-normalised; diagonal = self-consistent transfer",
                  x = "Transferred label", y = "Old SOM label") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

#' Return worst-projecting query SOM nodes
#'
#' Returns the `n` query nodes with the lowest match fraction — the nodes
#' whose mass the OT solver could not route to the reference — sorted
#' ascending. Includes the dominant reference label each node maps to.
#'
#' @param fit A `somalign_fit` object.
#' @param n Number of nodes to return. Default `10`.
#'
#' @return A data frame with columns `query_node`, `query_mass`,
#'   `transported_mass`, `match_fraction`, `correction_allowed`,
#'   `correction_norm`, and `top_ref_label`, one row per node.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_worst_nodes(fit, n = 4)
#' @export
somalign_worst_nodes <- function(fit, n = 10L) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  n <- as.integer(n)
  if (is.na(n) || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)
  nd <- .somalign_node_df(fit)
  nd <- nd[order(nd$match_fraction), ]
  nd[seq_len(min(n, nrow(nd))), ]
}

# ---- before-projection: plot functions --------------------------------------

#' Plot query vs reference SOM code ranges per marker
#'
#' Visualises the min-to-max range of each marker's SOM codes for both the
#' query and reference codebooks. Colours indicate the overlap flag computed
#' by [somalign_check_codebook_alignment()]: `ok` (green), `warning` (orange),
#' or `critical` (red).
#'
#' @param check A `somalign_codebook_check` object, as returned by
#'   [somalign_check_codebook_alignment()].
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' chk <- somalign_check_codebook_alignment(qry$codebook, ref,
#'                                          stop_if_critical = FALSE)
#' somalign_plot_codebook_range(chk)
#' @export
somalign_plot_codebook_range <- function(check) {
  if (!inherits(check, "somalign_codebook_check"))
    stop("`check` must be a somalign_codebook_check object.", call. = FALSE)
  pf <- check$per_feature
  required_cols <- c("feature", "ref_min", "ref_max", "query_min", "query_max", "flag")
  miss_cols <- setdiff(required_cols, names(pf))
  if (length(miss_cols))
    stop("Malformed `check` object; per_feature missing: ",
         paste(miss_cols, collapse = ", "), call. = FALSE)
  segs <- rbind(
    data.frame(feature = pf$feature, src = "Reference",
               xmin = pf$ref_min, xmax = pf$ref_max,
               flag  = "ok", stringsAsFactors = FALSE),
    data.frame(feature = pf$feature, src = "Query",
               xmin = pf$query_min, xmax = pf$query_max,
               flag  = pf$flag, stringsAsFactors = FALSE)
  )
  segs$feature <- factor(segs$feature, levels = pf$feature)
  segs$src     <- factor(segs$src, levels = c("Reference", "Query"))
  ggplot2::ggplot(segs,
    ggplot2::aes(x = xmin, xend = xmax, y = src, yend = src, colour = flag)) +
    ggplot2::geom_segment(linewidth = 3, lineend = "round") +
    ggplot2::geom_point(ggplot2::aes(x = xmin), size = 2) +
    ggplot2::geom_point(ggplot2::aes(x = xmax), size = 2) +
    ggplot2::scale_colour_manual(
      "Overlap",
      values = c(ok = "#1a9641", warning = "#fd8d3c", critical = "#d73027")) +
    ggplot2::facet_wrap(~feature, scales = "free_x") +
    ggplot2::labs(title = "Query vs reference SOM code ranges per marker",
                  x = "Value (reference-scaled)", y = NULL) +
    ggplot2::theme_minimal()
}

#' Plot per-marker cell distributions before projection
#'
#' Density plot of query cells (reference-scaled, downsampled) faceted by
#' marker. When a `somalign_reference` object is supplied via `reference`, the
#' reference SOM code values for that marker are overlaid as a rug of node
#' prototypes (red tick marks). When raw reference cell data are available,
#' pass them as `reference_data` for a true cell-vs-cell density comparison.
#'
#' @param query A `somalign_query` object.
#' @param reference Optional `somalign_reference` object. When supplied, its
#'   SOM codebook values are shown as a rug of node prototypes.
#' @param reference_data Optional numeric matrix of reference cells in
#'   reference-scaled space (cells x features). When supplied, a second
#'   density curve is shown instead of the codebook rug. Takes precedence over
#'   `reference`.
#' @param features Character vector of features to plot. `NULL` (default) uses
#'   all features in `query`.
#' @param downsample Maximum number of cells to subsample from `query` (and
#'   `reference_data` when supplied) for plotting speed. Default `2000`.
#' @param seed Integer seed for the random subsample. Default `1`.
#'
#' @return A `ggplot` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' somalign_plot_marker_distributions(qry, reference = ref)
#' @export
somalign_plot_marker_distributions <- function(query, reference = NULL,
                                               reference_data = NULL,
                                               features = NULL,
                                               downsample = 2000L,
                                               seed = 1L) {
  if (!inherits(query, "somalign_query"))
    stop("`query` must be a somalign_query object.", call. = FALSE)
  if (!is.null(reference) && !inherits(reference, "somalign_reference"))
    stop("`reference` must be a somalign_reference object or NULL.", call. = FALSE)
  if (!is.null(features) && !is.character(features))
    stop("`features` must be a character vector or NULL.", call. = FALSE)
  if (!is.numeric(downsample) || length(downsample) != 1L ||
      is.na(downsample) || downsample < 1)
    stop("`downsample` must be a single positive number.", call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed))
    stop("`seed` must be a single numeric scalar.", call. = FALSE)
  d <- .somalign_dist_data(query, reference, reference_data, features, downsample, seed)
  p <- ggplot2::ggplot(d$query, ggplot2::aes(x = value)) +
    ggplot2::geom_density(fill = "#4575b4", colour = "#4575b4", alpha = 0.5) +
    ggplot2::facet_wrap(~feature, scales = "free") +
    ggplot2::labs(title = "Per-marker distributions (reference-scaled space)",
                  x = "Scaled value", y = "Density") +
    ggplot2::theme_minimal()
  if (!is.null(d$ref_cells)) {
    p <- p + ggplot2::geom_density(
      data = d$ref_cells, fill = "#d73027", colour = "#d73027", alpha = 0.4) +
      ggplot2::labs(subtitle = "Blue = query cells | Red = reference cells")
  } else if (!is.null(d$ref_cb)) {
    p <- p + ggplot2::geom_rug(
      data = d$ref_cb, colour = "#d73027", alpha = 0.8, sides = "b") +
      ggplot2::labs(subtitle = "Blue density = query cells | Red rug = reference SOM nodes")
  }
  p
}
