#' Return per-sample somalign results
#'
#' Direct reference projection columns are canonical. Corrected projection
#' columns are auxiliary and should be used for visualisation, annotation, and
#' triage rather than feature-level differential testing.
#'
#' `transferred_label` is the query node's top label transfer choice, an
#' argmax over `correspondence %*% reference$label_prob`. When two labels
#' receive close probability mass, this argmax can be brittle: a cell can
#' land on a node whose true dominant population is clean, yet still get the
#' runner-up label because it narrowly lost the top spot. `transferred_label_second`
#' and `transferred_label_margin` (top confidence minus second confidence,
#' from the same query node) expose this directly, so close calls can be
#' triaged downstream instead of silently trusted at face value.
#'
#' @param fit A `somalign_fit` object.
#' @param data Optional data frame to append after result columns.
#' @param outside_pvalue_threshold Optional numeric in \[0, 1\]. When
#'   supplied, adds a boolean `outside_reference_pvalue_flag` column,
#'   `TRUE` when `outside_reference_pvalue < outside_pvalue_threshold`.
#'   `NULL` (default) omits the column.
#' @param include_correction Logical. When `FALSE`, drops the auxiliary
#'   correction columns (`corrected_som_unit`, `corrected_som_distance`,
#'   `corrected_som_distance_threshold`, `corrected_outside_reference_distance`,
#'   `correction_norm`) for a label-transfer-focused result. The corrected
#'   coordinates are a diagnostic, not a batch-corrected expression product;
#'   see the package's topology audit for why they can over-merge populations.
#'   Default `TRUE` (unchanged behaviour).
#'
#' @return A data frame with direct and corrected projection columns, plus
#'   `transferred_label_second`, `transferred_label_second_confidence`, and
#'   `transferred_label_margin` for triaging low-margin label transfers, and
#'   `outside_reference_surprisal`, `outside_reference_pvalue`, and
#'   `outside_reference_top_marker`: a calibrated chi-squared alternative to
#'   `outside_reference_distance` that weights per-marker deviations from a
#'   query cell's assigned reference node by that node's own per-marker
#'   variance (`reference$node_var`). `NA` when the reference was built
#'   without `node_var` (e.g. `compute_node_var = FALSE`, or a reference
#'   built before this feature). The chi-squared calibration assumes a
#'   diagonal-Gaussian node model; it is anti-conservative for heavy-tailed
#'   (e.g. lognormal) marker distributions, but remains a useful *relative*
#'   ranking of anomalous cells regardless.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_results(fit)
#' @export
somalign_results <- function(fit, data = NULL, outside_pvalue_threshold = NULL,
                             include_correction = TRUE) {
  if (!inherits(fit, "somalign_fit")) {
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  }
  if (!is.null(outside_pvalue_threshold))
    .somalign_check_prob_scalar(outside_pvalue_threshold, "outside_pvalue_threshold")
  .somalign_check_flag(include_correction, "include_correction")

  direct <- fit$projection$direct
  li <- .somalign_results_label_info(fit, direct)
  surpr <- .somalign_node_surprisal_chunked(fit$query$scaled_data, direct$unit, fit$reference)
  out <- .somalign_results_df(fit, direct, li, surpr)

  if (!include_correction) {
    correction_cols <- c("corrected_som_unit", "corrected_som_distance",
                         "corrected_som_distance_threshold",
                         "corrected_outside_reference_distance", "correction_norm")
    out <- out[, setdiff(names(out), correction_cols), drop = FALSE]
  }

  if (!is.null(outside_pvalue_threshold)) {
    out$outside_reference_pvalue_flag <-
      !is.na(out$outside_reference_pvalue) &
      out$outside_reference_pvalue < outside_pvalue_threshold
  }
  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (nrow(data) != nrow(out)) {
      stop("`data` must have one row per query sample.", call. = FALSE)
    }
    out <- cbind(out, data)
  }
  out
}

# Cell-level label-transfer summary: the headline product of a fit. Maps the
# node-level fit$label_transfer decisions to query cells via sample_unit and
# summarises acceptance, confidence, margin, and the accepted class mix. Used
# by print()/summary() and available as a lightweight accessor. Returns a list
# with `enabled = FALSE` when the reference carried no labels.
.somalign_label_summary <- function(fit) {
  lp <- fit$reference$label_prob
  if (is.null(lp) || ncol(lp) == 0L) return(list(enabled = FALSE))
  lt <- fit$label_transfer
  unit <- fit$query$sample_unit
  accepted <- lt$accepted[unit]
  accepted[is.na(accepted)] <- FALSE
  conf <- lt$confidence[unit]
  second <- lt$second_confidence[unit]
  margin <- conf - ifelse(is.na(second), 0, second)
  labels <- lt$label[unit]
  labels[!accepted] <- NA_character_
  n_cells <- length(unit)
  list(
    enabled = TRUE,
    n_cells = n_cells,
    accepted_fraction = mean(accepted),
    n_accepted = sum(accepted),
    n_classes = length(unique(stats::na.omit(labels))),
    class_distribution = sort(table(labels), decreasing = TRUE),
    median_confidence = stats::median(conf[accepted], na.rm = TRUE),
    median_margin = stats::median(margin[accepted], na.rm = TRUE),
    confidence_quartiles = stats::quantile(conf[accepted], c(0.25, 0.5, 0.75), na.rm = TRUE)
  )
}

.somalign_results_label_info <- function(fit, direct) {
  query_unit <- fit$query$sample_unit
  transferred <- fit$label_transfer
  accepted <- transferred$accepted[query_unit]
  second_confidence <- transferred$second_confidence[query_unit]
  ref_top <- .somalign_reference_top_labels(fit$reference)
  list(
    query_unit = query_unit,
    label = transferred$label[query_unit],
    confidence = ifelse(accepted, transferred$confidence[query_unit], NA_real_),
    accepted = accepted,
    second_label = transferred$second_label[query_unit],
    second_confidence = second_confidence,
    margin = transferred$confidence[query_unit] -
      ifelse(is.na(second_confidence), 0, second_confidence),
    old_label = ref_top$label[direct$unit],
    old_confidence = ref_top$confidence[direct$unit]
  )
}

.somalign_results_df <- function(fit, direct, li, surpr) {
  corrected <- fit$projection$corrected
  data.frame(
    sample_id = fit$query$sample_id,
    query_som_unit = li$query_unit,
    old_som_unit = direct$unit,
    old_som_distance = direct$distance,
    old_som_distance_threshold = direct$threshold,
    outside_reference_distance = direct$outside,
    outside_reference_surprisal = surpr$surprisal,
    outside_reference_pvalue = surpr$pvalue,
    outside_reference_top_marker = surpr$top_marker,
    final_status = ifelse(
      is.na(direct$outside),
      "unknown_reference_distance",
      ifelse(direct$outside, "outside_reference", "inside_reference")
    ),
    old_som_label = li$old_label,
    old_som_label_confidence = li$old_confidence,
    corrected_som_unit = corrected$unit,
    corrected_som_distance = corrected$distance,
    corrected_som_distance_threshold = corrected$threshold,
    corrected_outside_reference_distance = corrected$outside,
    correction_norm = fit$projection$correction_norm,
    transferred_label = li$label,
    transferred_label_confidence = li$confidence,
    transferred_label_accepted = li$accepted,
    transferred_label_second = li$second_label,
    transferred_label_second_confidence = li$second_confidence,
    transferred_label_margin = li$margin
  )
}

.somalign_reference_top_labels <- function(reference) {
  label_prob <- reference$label_prob
  n_nodes <- nrow(reference$codebook)
  if (is.null(label_prob) || ncol(label_prob) == 0) {
    return(list(
      label = rep(NA_character_, n_nodes),
      confidence = rep(NA_real_, n_nodes)
    ))
  }
  label_names <- colnames(label_prob)
  row_sums <- rowSums(label_prob)
  has_mass <- row_sums > 0
  idx <- max.col(label_prob, ties.method = "first")
  label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  label[has_mass] <- label_names[idx[has_mass]]
  confidence[has_mass] <- label_prob[cbind(which(has_mass), idx[has_mass])]
  list(label = label, confidence = confidence)
}
