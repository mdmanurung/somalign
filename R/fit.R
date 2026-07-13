#' Align a query SOM to a reference SOM
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Entropic regularisation strength. The cost matrix is
#'   normalised by its median positive entry before computing the Sinkhorn
#'   kernel, so `epsilon` is approximately scale- and dimension-invariant.
#'   Values around `0.5` give meaningful regularisation for typical z-scored
#'   SOM codebooks; very small values (< 0.1) make the transport increasingly
#'   discrete. The normalisation scale is stored in
#'   `diagnostics$solver$cost_scale`.
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver Sinkhorn solver variant. `"internal"` (default) and `"auto"`
#'   both use the primal-domain scaling iteration. `"log_domain"` uses a
#'   numerically stable log-potential variant that avoids kernel underflow for
#'   small `epsilon` or high-dimensional codebooks; it is slower per iteration
#'   but tolerates cost/epsilon ratios that cause `"internal"` to warn.
#' @param min_match_fraction Minimum transported fraction required before a
#'   query node label transfer is accepted.
#' @param confidence_threshold Minimum top-label probability required before a
#'   query node label transfer is accepted.
#' @param correction_min_mass Minimum transported node mass required before a
#'   correction shift is applied. Corrections also require the node match
#'   fraction to pass `min_match_fraction`.
#' @param max_iter Maximum internal Sinkhorn iterations.
#' @param tol Internal Sinkhorn convergence tolerance.
#' @param chunk_size Integer. Number of samples to project per chunk when
#'   computing nearest reference node. Use `Inf` or `NULL` for no chunking
#'   (allocates a full n_samples x n_nodes matrix). Default `10000L`.
#'
#' @details
#' The transport plan row sums will not equal `query$node_masses` exactly -- this
#' is by design. Unbalanced optimal transport allows mass destruction, so some
#' query mass may be absorbed rather than transported. Deviation grows with lower
#' `rho_query` / `rho_ref` values and higher `epsilon`. Use
#' `diagnostics$ot$max_row_mass_error` to quantify the deviation in a given fit;
#' for near-balanced data, increase `rho_query` (e.g. `rho_query = 10`) to
#' enforce tighter marginal constraints. A warning is emitted automatically when
#' more than 50% of query mass is destroyed.
#'
#' The cost matrix is normalised by its median positive entry before the
#' Sinkhorn kernel is computed. This makes `epsilon` scale- and
#' dimension-invariant: the same value produces the same degree of regularisation
#' regardless of the number of features or the spread of codebook coordinates.
#' The raw normalisation factor is stored as `diagnostics$solver$cost_scale`.
#'
#' @return A `somalign_fit` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' @export
somalign_fit <- function(query,
                         reference,
                         epsilon = 0.5,
                         rho_query = 1,
                         rho_ref = 1,
                         solver = c("internal", "log_domain", "auto"),
                         min_match_fraction = 0.05,
                         confidence_threshold = 0.6,
                         correction_min_mass = 1e-8,
                         max_iter = 1000,
                         tol = 1e-7,
                         chunk_size = 10000L) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  solver <- match.arg(solver, c("internal", "log_domain", "auto"))

  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol
  )
  .somalign_finish_fit(
    query, reference, transport,
    min_match_fraction, confidence_threshold, correction_min_mass,
    chunk_size, epsilon, rho_query, rho_ref
  )
}

.somalign_finish_fit <- function(query, reference, transport,
                                 min_match_fraction, confidence_threshold,
                                 correction_min_mass, chunk_size,
                                 epsilon, rho_query, rho_ref,
                                 anchors = NULL) {
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
  .somalign_new_fit(
    query, reference, transport, label_transfer, node_shifts, projection, diagnostics,
    anchors = anchors
  )
}

.somalign_align_transport <- function(query, reference, epsilon, rho_query,
                                      rho_ref, solver, max_iter, tol,
                                      cost_bonus = NULL) {
  cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)
  cost_scale <- stats::median(cost[cost > 0])
  if (!is.finite(cost_scale) || cost_scale == 0) {
    cost_scale <- 1
  }
  cost_normalized <- cost / cost_scale
  if (!is.null(cost_bonus)) {
    cost_normalized <- pmax(cost_normalized - cost_bonus, 0)
  }
  ot <- .somalign_solve_ot(
    cost = cost_normalized,
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
  n_over <- sum(match_mass_ratio > 1)
  if (n_over > 0) {
    message(sprintf(
      "somalign_fit: %d query node(s) have match_mass_ratio > 1 (max %.2f); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.",
      n_over, max(match_mass_ratio)
    ))
  }
  list(
    cost = cost,
    cost_scale = cost_scale,
    ot = ot,
    plan = plan,
    correspondence = correspondence,
    row_mass = row_mass,
    col_mass = col_mass,
    match_mass_ratio = match_mass_ratio,
    match_fraction = match_fraction
  )
}

.somalign_project_pair <- function(query, reference, node_shifts, chunk_size) {
  direct <- .somalign_project_samples(query$scaled_data, reference, chunk_size = chunk_size)
  corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]
  corrected <- .somalign_project_samples(corrected_matrix, reference, chunk_size = chunk_size)
  correction_norm <- sqrt(rowSums(node_shifts[query$sample_unit, , drop = FALSE]^2))
  list(direct = direct, corrected = corrected, correction_norm = correction_norm)
}

.somalign_build_diagnostics <- function(transport, query, reference, node_shifts,
                                        projection, epsilon, rho_query,
                                        rho_ref) {
  ot <- transport$ot
  plan <- transport$plan
  row_mass <- transport$row_mass
  col_mass <- transport$col_mass
  direct <- projection$direct
  corrected <- projection$corrected
  list(
    solver = list(
      requested = ot$requested_solver,
      used = ot$solver,
      notes = ot$notes,
      iterations = ot$iterations,
      converged = ot$converged,
      final_delta = ot$final_delta,
      epsilon = epsilon,
      rho_query = rho_query,
      rho_ref = rho_ref,
      cost_scale = transport$cost_scale,
      rel_marginal_row_error = max(abs(row_mass - query$node_masses)) /
        max(sum(query$node_masses), .Machine$double.eps),
      rel_marginal_col_error = max(abs(col_mass - reference$node_masses)) /
        max(sum(reference$node_masses), .Machine$double.eps)
    ),
    ot = list(
      transport_mass = sum(plan),
      row_mass = row_mass,
      col_mass = col_mass,
      query_mass = query$node_masses,
      reference_mass = reference$node_masses,
      match_fraction = transport$match_fraction,
      match_mass_ratio = transport$match_mass_ratio,
      max_row_mass_error = max(abs(row_mass - query$node_masses)),
      max_col_mass_error = max(abs(col_mass - reference$node_masses))
    ),
    nodes = data.frame(
      query_node = seq_len(nrow(query$codebook)),
      query_mass = query$node_masses,
      transported_mass = row_mass,
      match_fraction = transport$match_fraction,
      correction_allowed = attr(node_shifts, "correction_allowed"),
      correction_norm = sqrt(rowSums(node_shifts^2))
    ),
    projection = list(
      outside_direct_fraction = mean(direct$outside),
      outside_corrected_fraction = mean(corrected$outside)
    )
  )
}

.somalign_fit_warnings <- function(diagnostics) {
  query_total_mass <- sum(diagnostics$ot$query_mass)
  if (query_total_mass > 0 && diagnostics$ot$transport_mass < 0.5 * query_total_mass) {
    warning(
      sprintf(
        "High mass destruction: only %.1f%% of query mass was transported. ",
        100 * diagnostics$ot$transport_mass / query_total_mass
      ),
      "Consider raising rho_query and rho_ref, or increasing epsilon.",
      call. = FALSE
    )
  }
  outside_frac <- diagnostics$projection$outside_direct_fraction
  if (is.finite(outside_frac) && outside_frac > 0.5) {
    warning(
      sprintf(
        "%.1f%% of query samples project outside reference distance thresholds. ",
        100 * outside_frac
      ),
      "This may indicate a distributional mismatch or a coordinate-space misconfiguration.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.somalign_new_fit <- function(query, reference, transport, label_transfer,
                              node_shifts, projection, diagnostics,
                              anchors = NULL) {
  fit <- structure(
    list(
      query = query,
      reference = reference,
      cost = transport$cost,
      transport_plan = transport$plan,
      correspondence = transport$correspondence,
      label_transfer = label_transfer,
      node_shifts = node_shifts,
      projection = list(
        direct = projection$direct,
        corrected = projection$corrected,
        correction_norm = projection$correction_norm
      ),
      diagnostics = diagnostics
    ),
    class = "somalign_fit"
  )
  if (!is.null(anchors)) {
    fit$anchors <- anchors
  }
  fit
}

.somalign_transfer_labels <- function(correspondence,
                                      label_prob,
                                      match_fraction,
                                      min_match_fraction,
                                      confidence_threshold) {
  n_nodes <- nrow(correspondence)
  if (is.null(label_prob) || ncol(label_prob) == 0) {
    return(.somalign_empty_label_transfer(n_nodes, match_fraction))
  }

  probs <- correspondence %*% label_prob
  label_names <- colnames(label_prob)
  row_sums <- rowSums(probs)
  has_mass <- row_sums > 0
  probs_norm <- probs
  probs_norm[has_mass, ] <- probs[has_mass, , drop = FALSE] / row_sums[has_mass]
  top_idx <- max.col(probs_norm, ties.method = "first")
  top_label <- rep(NA_character_, n_nodes)
  confidence <- rep(NA_real_, n_nodes)
  top_label[has_mass] <- label_names[top_idx[has_mass]]
  confidence[has_mass] <- probs_norm[cbind(which(has_mass), top_idx[has_mass])]
  second <- .somalign_second_labels(probs_norm, top_idx, has_mass, label_names)
  entropy <- vapply(
    seq_len(n_nodes),
    function(i) if (has_mass[i]) .somalign_entropy(probs_norm[i, ]) else NA_real_,
    numeric(1)
  )

  accepted <- is.finite(match_fraction) &
    match_fraction >= min_match_fraction &
    is.finite(confidence) &
    confidence >= confidence_threshold
  top_label[!accepted] <- NA_character_

  data.frame(
    query_node = seq_len(n_nodes),
    label = top_label,
    confidence = confidence,
    second_label = second$second_label,
    second_confidence = second$second_confidence,
    entropy = entropy,
    match_fraction = match_fraction,
    accepted = accepted
  )
}

.somalign_empty_label_transfer <- function(n_nodes, match_fraction) {
  data.frame(
    query_node = seq_len(n_nodes),
    label = rep(NA_character_, n_nodes),
    confidence = rep(NA_real_, n_nodes),
    second_label = rep(NA_character_, n_nodes),
    second_confidence = rep(NA_real_, n_nodes),
    entropy = rep(NA_real_, n_nodes),
    match_fraction = match_fraction,
    accepted = rep(FALSE, n_nodes)
  )
}

.somalign_second_labels <- function(probs_norm, top_idx, has_mass, label_names) {
  n_nodes <- nrow(probs_norm)
  second_label <- rep(NA_character_, n_nodes)
  second_confidence <- rep(NA_real_, n_nodes)
  if (ncol(probs_norm) == 1L) {
    return(list(second_label = second_label, second_confidence = second_confidence))
  }

  probs_second <- probs_norm
  probs_second[cbind(seq_len(n_nodes), top_idx)] <- 0
  second_idx <- max.col(probs_second, ties.method = "first")
  second_label[has_mass] <- label_names[second_idx[has_mass]]
  second_confidence[has_mass] <- probs_second[cbind(which(has_mass), second_idx[has_mass])]
  second_label[has_mass & (is.na(second_confidence) | second_confidence == 0)] <- NA_character_
  list(second_label = second_label, second_confidence = second_confidence)
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

.somalign_project_samples <- function(scaled_data, reference, chunk_size = 10000L) {
  projected <- .somalign_nearest_code_chunked(scaled_data, reference$codebook, chunk_size = chunk_size)
  threshold <- .somalign_thresholds(reference, projected$unit)
  list(
    unit = projected$unit,
    distance = projected$distance,
    threshold = threshold,
    outside = projected$distance > threshold
  )
}
