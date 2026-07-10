#' Align a query SOM to a reference SOM
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Entropic regularisation strength.
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver `"auto"`, `"pot"`, or `"internal"`.
#' @param min_match_fraction Minimum transported fraction required before a
#'   query node label transfer is accepted.
#' @param confidence_threshold Minimum top-label probability required before a
#'   query node label transfer is accepted.
#' @param correction_min_mass Minimum transported node mass required before a
#'   correction shift is applied. Corrections also require the node match
#'   fraction to pass `min_match_fraction`.
#' @param max_iter Maximum internal Sinkhorn iterations.
#' @param tol Internal Sinkhorn convergence tolerance.
#'
#' @return A `somalign_fit` object.
#' @export
somalign_fit <- function(query,
                         reference,
                         epsilon = 0.05,
                         rho_query = 1,
                         rho_ref = 1,
                         solver = c("auto", "pot", "internal"),
                         min_match_fraction = 0.05,
                         confidence_threshold = 0.6,
                         correction_min_mass = 1e-8,
                         max_iter = 1000,
                         tol = 1e-7) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  solver <- match.arg(solver)

  cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)
  ot <- .somalign_solve_ot(
    cost = cost,
    a = query$node_masses,
    b = reference$node_masses,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_ref = rho_ref,
    solver = solver,
    max_iter = max_iter,
    tol = tol
  )
  plan <- ot$plan
  correspondence <- .somalign_row_normalize(plan)
  row_mass <- rowSums(plan)
  col_mass <- colSums(plan)
  match_mass_ratio <- ifelse(query$node_masses > 0, row_mass / query$node_masses, 0)
  match_fraction <- pmin(match_mass_ratio, 1)

  label_transfer <- .somalign_transfer_labels(
    correspondence = correspondence,
    label_prob = reference$label_prob,
    match_fraction = match_fraction,
    min_match_fraction = min_match_fraction,
    confidence_threshold = confidence_threshold
  )

  node_shifts <- .somalign_node_shifts(
    query_codebook = query$codebook,
    reference_codebook = reference$codebook,
    correspondence = correspondence,
    row_mass = row_mass,
    match_fraction = match_fraction,
    min_match_fraction = min_match_fraction,
    correction_min_mass = correction_min_mass
  )

  direct <- .somalign_project_samples(query$scaled_data, reference)
  corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]
  corrected <- .somalign_project_samples(corrected_matrix, reference)
  correction_norm <- sqrt(rowSums(node_shifts[query$sample_unit, , drop = FALSE]^2))

  diagnostics <- list(
    solver = list(
      requested = ot$requested_solver,
      used = ot$solver,
      notes = ot$notes,
      iterations = ot$iterations,
      epsilon = epsilon,
      rho_query = rho_query,
      rho_ref = rho_ref
    ),
    ot = list(
      transport_mass = sum(plan),
      row_mass = row_mass,
      col_mass = col_mass,
      query_mass = query$node_masses,
      reference_mass = reference$node_masses,
      match_fraction = match_fraction,
      match_mass_ratio = match_mass_ratio,
      max_row_mass_error = max(abs(row_mass - query$node_masses)),
      max_col_mass_error = max(abs(col_mass - reference$node_masses))
    ),
    nodes = data.frame(
      query_node = seq_len(nrow(query$codebook)),
      query_mass = query$node_masses,
      transported_mass = row_mass,
      match_fraction = match_fraction,
      correction_allowed = attr(node_shifts, "correction_allowed"),
      correction_norm = sqrt(rowSums(node_shifts^2))
    ),
    projection = list(
      outside_direct_fraction = mean(direct$outside),
      outside_corrected_fraction = mean(corrected$outside)
    )
  )

  structure(
    list(
      query = query,
      reference = reference,
      cost = cost,
      transport_plan = plan,
      correspondence = correspondence,
      label_transfer = label_transfer,
      node_shifts = node_shifts,
      projection = list(
        direct = direct,
        corrected = corrected,
        correction_norm = correction_norm
      ),
      diagnostics = diagnostics
    ),
    class = "somalign_fit"
  )
}

.somalign_transfer_labels <- function(correspondence,
                                      label_prob,
                                      match_fraction,
                                      min_match_fraction,
                                      confidence_threshold) {
  n_nodes <- nrow(correspondence)
  if (is.null(label_prob) || ncol(label_prob) == 0) {
    return(data.frame(
      query_node = seq_len(n_nodes),
      label = rep(NA_character_, n_nodes),
      confidence = rep(NA_real_, n_nodes),
      second_label = rep(NA_character_, n_nodes),
      second_confidence = rep(NA_real_, n_nodes),
      entropy = rep(NA_real_, n_nodes),
      match_fraction = match_fraction,
      accepted = rep(FALSE, n_nodes)
    ))
  }

  probs <- correspondence %*% label_prob
  label_names <- colnames(label_prob)
  top_label <- rep(NA_character_, n_nodes)
  second_label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  second_confidence <- rep(NA_real_, n_nodes)
  entropy <- rep(NA_real_, n_nodes)

  for (i in seq_len(n_nodes)) {
    row <- as.numeric(probs[i, ])
    if (sum(row) > 0) {
      row <- row / sum(row)
      ord <- order(row, decreasing = TRUE)
      top_label[i] <- label_names[ord[1]]
      confidence[i] <- row[ord[1]]
      if (length(ord) > 1) {
        second_label[i] <- label_names[ord[2]]
        second_confidence[i] <- row[ord[2]]
      }
      entropy[i] <- .somalign_entropy(row)
    }
  }

  accepted <- is.finite(match_fraction) &
    match_fraction >= min_match_fraction &
    is.finite(confidence) &
    confidence >= confidence_threshold
  top_label[!accepted] <- NA_character_

  data.frame(
    query_node = seq_len(n_nodes),
    label = top_label,
    confidence = confidence,
    second_label = second_label,
    second_confidence = second_confidence,
    entropy = entropy,
    match_fraction = match_fraction,
    accepted = accepted
  )
}

.somalign_node_shifts <- function(query_codebook,
                                  reference_codebook,
                                  correspondence,
                                  row_mass,
                                  match_fraction,
                                  min_match_fraction,
                                  correction_min_mass) {
  barycentric_reference <- correspondence %*% reference_codebook
  shifts <- matrix(0, nrow = nrow(query_codebook), ncol = ncol(query_codebook))
  colnames(shifts) <- colnames(query_codebook)
  strong <- row_mass >= correction_min_mass &
    is.finite(match_fraction) &
    match_fraction >= min_match_fraction &
    rowSums(correspondence) > 0
  shifts[strong, ] <- barycentric_reference[strong, , drop = FALSE] -
    query_codebook[strong, , drop = FALSE]
  attr(shifts, "correction_allowed") <- strong
  shifts
}

.somalign_project_samples <- function(scaled_data, reference) {
  projected <- .somalign_nearest_code(scaled_data, reference$codebook)
  threshold <- .somalign_thresholds(reference, projected$unit)
  list(
    unit = projected$unit,
    distance = projected$distance,
    threshold = threshold,
    outside = projected$distance > threshold
  )
}
