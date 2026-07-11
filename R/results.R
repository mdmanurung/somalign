#' Return per-sample somalign results
#'
#' Direct reference projection columns are canonical. Corrected projection
#' columns are auxiliary and should be used for visualisation, annotation, and
#' triage rather than feature-level differential testing.
#'
#' @param fit A `somalign_fit` object.
#' @param data Optional data frame to append after result columns.
#'
#' @return A data frame with direct and corrected projection columns.
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
somalign_results <- function(fit, data = NULL) {
  if (!inherits(fit, "somalign_fit")) {
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  }
  direct <- fit$projection$direct
  corrected <- fit$projection$corrected
  query_unit <- fit$query$sample_unit
  transferred <- fit$label_transfer
  label <- transferred$label[query_unit]
  accepted <- transferred$accepted[query_unit]
  confidence <- ifelse(accepted, transferred$confidence[query_unit], NA_real_)
  ref_top <- .somalign_reference_top_labels(fit$reference)
  old_label <- ref_top$label[direct$unit]
  old_confidence <- ref_top$confidence[direct$unit]

  out <- data.frame(
    sample_id = fit$query$sample_id,
    query_som_unit = query_unit,
    old_som_unit = direct$unit,
    old_som_distance = direct$distance,
    old_som_distance_threshold = direct$threshold,
    outside_reference_distance = direct$outside,
    final_status = ifelse(
      is.na(direct$outside),
      "unknown_reference_distance",
      ifelse(direct$outside, "outside_reference", "inside_reference")
    ),
    old_som_label = old_label,
    old_som_label_confidence = old_confidence,
    corrected_som_unit = corrected$unit,
    corrected_som_distance = corrected$distance,
    corrected_som_distance_threshold = corrected$threshold,
    corrected_outside_reference_distance = corrected$outside,
    correction_norm = fit$projection$correction_norm,
    transferred_label = label,
    transferred_label_confidence = confidence,
    transferred_label_accepted = accepted
  )

  if (!is.null(data)) {
    data <- as.data.frame(data)
    if (nrow(data) != nrow(out)) {
      stop("`data` must have one row per query sample.", call. = FALSE)
    }
    out <- cbind(out, data)
  }
  out
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
