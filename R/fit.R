#' Align a query SOM to a reference SOM
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon Entropic regularisation strength. The cost matrix is
#'   normalised by its median positive entry before computing the Sinkhorn
#'   kernel, so `epsilon` is approximately scale- and dimension-invariant.
#'   The default `0.1` gives a sharp transport plan that preserves cell-type
#'   specificity for typical z-scored SOM codebooks. Larger values (0.3–0.5)
#'   produce smoother, more diffuse plans that can help convergence on noisy
#'   or high-dimensional data but dilute label posteriors and increase
#'   barycentric shrinkage in the corrected projection. Very small values
#'   (< 0.05) make the transport increasingly discrete and may require
#'   `solver = "log_domain"` for numerical stability. The normalisation scale
#'   is stored in `diagnostics$solver$cost_scale`.
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver Sinkhorn solver variant. `"internal"` (default) and `"auto"`
#'   both use the primal-domain scaling iteration. `"log_domain"` uses a
#'   numerically stable log-potential variant that avoids kernel underflow for
#'   small `epsilon` or high-dimensional codebooks; it is slower per iteration
#'   but tolerates cost/epsilon ratios that cause `"internal"` to warn.
#'   `"annealing"` runs the log-domain solver across a geometric epsilon
#'   cooling schedule (starting at `anneal_start * epsilon`, cooling to
#'   `epsilon` over `anneal_stages` stages), warm-starting each stage from the
#'   previous stage's dual potentials. Recommended for `label_guided` fits or
#'   any fit with small `epsilon` (< 0.05) where cold-start Sinkhorn is slow
#'   or non-convergent; never underflows, since it never exponentiates the
#'   kernel.
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
#' @param diagonal_boost Non-negative scalar. Amount by which to reduce the
#'   normalised OT cost for each query node's nearest reference node. A
#'   positive value makes the transport plan prefer identity-like mappings,
#'   shrinking over-correction when the two codebooks are already close. Zero
#'   (default) leaves the cost unchanged. Values around 0.1--0.5 are a
#'   reasonable starting point; very large values concentrate all mass on the
#'   diagonal and the plan degrades toward simple nearest-neighbour assignment.
#' @param label_guided Logical. When `TRUE`, uses `query$label_prob` and
#'   `reference$label_prob` to add a large cost penalty for node pairs whose
#'   dominant labels disagree, constraining OT to transport mass predominantly
#'   between concordant cell-type nodes. Nodes where the maximum label
#'   probability is below 0.5 are treated as unlabeled and are never penalized.
#'   Errors if `label_guided = TRUE` but either `label_prob` is `NULL`.
#' @param anneal_start Positive scalar >= 1. When `solver = "annealing"`, the
#'   starting epsilon is `anneal_start * epsilon`. Default `10`. Ignored when
#'   `solver != "annealing"`.
#' @param anneal_stages Positive integer. Number of cooling stages in the
#'   annealing schedule, including the final stage at the target `epsilon`.
#'   Default `10L`. A value of `1` degenerates to a cold-start log-domain
#'   solve. Ignored when `solver != "annealing"`.
#' @param anneal_factor Positive scalar < 1, or `NULL` (default). When not
#'   `NULL`, overrides the auto-computed per-stage cooling ratio. Ignored when
#'   `solver != "annealing"`.
#' @param feature_weights Either `NULL` (default, squared-Euclidean cost) or a
#'   named non-negative numeric vector with one entry per feature (explicit
#'   diagonal Mahalanobis weights on the OT cost). Weights are applied as
#'   `sqrt(w_f)` per-column scaling of both codebooks before the squared
#'   Euclidean distance is computed, yielding cost
#'   \eqn{\sum_f w_f (q_{if} - r_{jf})^2}. The resolved vector is stored in
#'   `fit$diagnostics$cost_metric$feature_weights`. Projection and threshold
#'   distances (`somalign_results()`) are unaffected -- weighting applies only
#'   to the OT cost. See [somalign_fit_anchored()] for `"anchor"`, which
#'   auto-estimates weights from anchor displacements.
#' @param laplacian_lambda Non-negative scalar. Graph-Laplacian regularisation
#'   strength for the node-shift field. When greater than zero, the M x p raw
#'   node shifts are smoothed by solving \eqn{(W + \lambda L)\,s^* = W\,s},
#'   where \eqn{W = \mathrm{diag}(\text{node\_masses})} (with
#'   `correction_allowed == FALSE` nodes zeroed out) and \eqn{L} is the graph
#'   Laplacian of the query SOM's hexagonal or rectangular neighbor graph.
#'   This penalises squared differences between adjacent-node shifts,
#'   producing a spatially coherent correction field instead of one where
#'   neighboring nodes can receive wildly different shifts from finite-sample
#'   OT noise. Default `0` (no smoothing, exact current behaviour). A natural
#'   starting range is `0.1`--`1.0` (same cost/squared-distance scale as
#'   `epsilon`); larger values increasingly collapse the field toward its
#'   mass-weighted mean. Requires the query SOM to carry 2-D grid coordinates
#'   (`query$som_query$grid$pts`, present for any `kohonen::som()`- or
#'   `kohonen::supersom()`-trained SOM); errors otherwise.
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
                         epsilon = 0.1,
                         rho_query = 1,
                         rho_ref = 1,
                         solver = c("internal", "log_domain", "auto", "annealing"),
                         min_match_fraction = 0.05,
                         confidence_threshold = 0.6,
                         correction_min_mass = 1e-8,
                         max_iter = 1000,
                         tol = 1e-7,
                         chunk_size = 10000L,
                         diagonal_boost = 0,
                         label_guided = FALSE,
                         anneal_start = 10,
                         anneal_stages = 10L,
                         anneal_factor = NULL,
                         feature_weights = NULL,
                         laplacian_lambda = 0) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  .somalign_check_pos_scalar(epsilon, "epsilon")
  .somalign_check_nonneg_scalar(diagonal_boost, "diagonal_boost")
  .somalign_check_fit_params(rho_query, rho_ref, min_match_fraction,
                             confidence_threshold, correction_min_mass,
                             max_iter, tol, chunk_size, label_guided)
  solver <- match.arg(solver, c("internal", "log_domain", "auto", "annealing"))
  if (identical(solver, "annealing"))
    .somalign_check_anneal_params(anneal_start, anneal_factor, anneal_stages)
  feature_weights <- .somalign_resolve_plain_feature_weights(
    feature_weights, colnames(query$codebook)
  )
  .somalign_check_nonneg_scalar(laplacian_lambda, "laplacian_lambda")
  label_mask <- .somalign_resolve_label_mask(query, reference, label_guided)
  shift_transform <- .somalign_make_laplacian_transform(query, laplacian_lambda)

  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol,
    diagonal_boost = diagonal_boost, label_mask = label_mask,
    anneal_start = anneal_start, anneal_factor = anneal_factor,
    anneal_stages = anneal_stages, feature_weights = feature_weights
  )
  .somalign_finish_fit(
    query, reference, transport,
    min_match_fraction, confidence_threshold, correction_min_mass,
    chunk_size, epsilon, rho_query, rho_ref, feature_weights = feature_weights,
    shift_transform = shift_transform
  )
}

.somalign_resolve_label_mask <- function(query, reference, label_guided) {
  if (!isTRUE(label_guided)) return(NULL)
  if (is.null(query$label_prob))
    stop("label_guided = TRUE but query$label_prob is NULL.", call. = FALSE)
  if (is.null(reference$label_prob))
    stop("label_guided = TRUE but reference$label_prob is NULL.", call. = FALSE)
  .somalign_build_label_mask(query$label_prob, reference$label_prob)
}

.somalign_resolve_plain_feature_weights <- function(feature_weights, features) {
  .somalign_check_feature_weights(feature_weights, features)
  if (identical(feature_weights, "anchor"))
    stop("`feature_weights = \"anchor\"` requires anchor data; use somalign_fit_anchored().",
         call. = FALSE)
  feature_weights
}

.somalign_finish_fit <- function(query, reference, transport,
                                 min_match_fraction, confidence_threshold,
                                 correction_min_mass, chunk_size,
                                 epsilon, rho_query, rho_ref,
                                 anchors = NULL,
                                 direct_cache = NULL,
                                 shift_transform = NULL,
                                 feature_weights = NULL) {
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
  if (!is.null(shift_transform)) {
    allowed <- attr(node_shifts, "correction_allowed")
    node_shifts <- shift_transform(node_shifts)
    attr(node_shifts, "correction_allowed") <- allowed
  }
  projection <- .somalign_project_pair(query, reference, node_shifts, chunk_size,
                                       direct_cache = direct_cache)
  diagnostics <- .somalign_build_diagnostics(
    transport, query, reference, node_shifts, projection, epsilon, rho_query, rho_ref,
    feature_weights = feature_weights
  )
  .somalign_fit_warnings(diagnostics)
  .somalign_new_fit(
    query, reference, transport, label_transfer, node_shifts, projection, diagnostics,
    anchors = anchors
  )
}

.somalign_build_label_mask <- function(query_label_prob, ref_label_prob) {
  if (!identical(colnames(query_label_prob), colnames(ref_label_prob))) {
    stop(
      "label_guided = TRUE requires query$label_prob and reference$label_prob ",
      "to have identical column names (same cell-type taxonomy). ",
      "query has: ", paste(colnames(query_label_prob), collapse = ", "), ". ",
      "reference has: ", paste(colnames(ref_label_prob), collapse = ", "), ".",
      call. = FALSE
    )
  }
  row_sum_q <- rowSums(query_label_prob)
  row_sum_r <- rowSums(ref_label_prob)
  q_norm <- query_label_prob / ifelse(row_sum_q > 0, row_sum_q, 1)
  r_norm <- ref_label_prob   / ifelse(row_sum_r > 0, row_sum_r, 1)
  q_dom <- max.col(q_norm, ties.method = "first")
  r_dom <- max.col(r_norm, ties.method = "first")
  q_unlab <- apply(q_norm, 1, max) < 0.5
  r_unlab <- apply(r_norm, 1, max) < 0.5
  mask <- outer(q_dom, r_dom, "!=")
  mask[q_unlab, ] <- FALSE
  mask[, r_unlab] <- FALSE
  mask
}

.somalign_prepare_cost <- function(cost, diagonal_boost, cost_bonus, label_mask) {
  cost_scale <- stats::median(cost[cost > 0])
  if (!is.finite(cost_scale) || cost_scale == 0) {
    cost_scale <- 1
  }
  cost_normalized <- cost / cost_scale
  if (diagonal_boost > 0) {
    nn_col <- max.col(-cost_normalized, ties.method = "first")
    idx <- cbind(seq_len(nrow(cost_normalized)), nn_col)
    cost_normalized[idx] <- pmax(cost_normalized[idx] - diagonal_boost, 0)
  }
  if (!is.null(cost_bonus)) {
    cost_normalized <- pmax(cost_normalized - cost_bonus, 0)
  }
  if (!is.null(label_mask)) {
    penalty <- max(cost_normalized) * 1e4
    cost_normalized[label_mask] <- cost_normalized[label_mask] + penalty
  }
  list(cost_normalized = cost_normalized, cost_scale = cost_scale)
}

.somalign_align_transport <- function(query, reference, epsilon, rho_query,
                                      rho_ref, solver, max_iter, tol,
                                      cost_bonus = NULL,
                                      diagonal_boost = 0,
                                      label_mask = NULL,
                                      anneal_start = 10,
                                      anneal_factor = NULL,
                                      anneal_stages = 10L,
                                      feature_weights = NULL) {
  if (!is.null(feature_weights)) {
    qcb <- .somalign_weighted_codebook(query$codebook, feature_weights)
    rcb <- .somalign_weighted_codebook(reference$codebook, feature_weights)
    cost <- .somalign_pairwise_distance(qcb, rcb)
  } else {
    cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)
  }
  prepared <- .somalign_prepare_cost(cost, diagonal_boost, cost_bonus, label_mask)
  cost_normalized <- prepared$cost_normalized
  cost_scale <- prepared$cost_scale
  ot <- .somalign_solve_ot(
    cost = cost_normalized,
    a = query$node_masses,
    b = reference$node_masses,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_ref = rho_ref,
    solver = solver,
    max_iter = max_iter,
    tol = tol,
    anneal_start = anneal_start,
    anneal_factor = anneal_factor,
    anneal_stages = anneal_stages
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

# Cheap OT-only sweep primitive shared by somalign_epsilon_sweep() and
# somalign_select_epsilon(): builds the cost matrix and solves the OT problem
# for a single epsilon, WITHOUT per-cell projection, node shifts, or label
# transfer. Deliberately bypasses .somalign_align_transport() (which emits a
# match_mass_ratio message meant for a single interactive fit -- noisy across
# a 25+ point epsilon grid).
.somalign_ot_sweep_one <- function(query, reference, epsilon,
                                   rho_query, rho_ref,
                                   solver, max_iter, tol,
                                   diagonal_boost = 0,
                                   label_mask = NULL,
                                   cost_bonus = NULL,
                                   anneal_start = 10,
                                   anneal_factor = NULL,
                                   anneal_stages = 10L) {
  cost <- .somalign_pairwise_distance(query$codebook, reference$codebook)
  prepared <- .somalign_prepare_cost(cost, diagonal_boost, cost_bonus, label_mask)
  ot <- .somalign_solve_ot(
    cost = prepared$cost_normalized,
    a = query$node_masses,
    b = reference$node_masses,
    epsilon = epsilon,
    rho_query = rho_query,
    rho_ref = rho_ref,
    solver = solver,
    max_iter = max_iter,
    tol = tol,
    anneal_start = anneal_start,
    anneal_factor = anneal_factor,
    anneal_stages = anneal_stages
  )
  mi <- .somalign_plan_mutual_information(ot$plan, cost = prepared$cost_normalized)
  list(
    epsilon = epsilon,
    plan = ot$plan,
    cost_scale = prepared$cost_scale,
    mutual_information = mi$mutual_information,
    conditional_entropy = mi$conditional_entropy,
    expected_cost = mi$expected_cost,
    transport_mass = sum(ot$plan),
    log_Z = ot$log_Z,
    iterations = ot$iterations,
    converged = ot$converged
  )
}

.somalign_project_pair <- function(query, reference, node_shifts, chunk_size,
                                   direct_cache = NULL) {
  direct <- if (!is.null(direct_cache)) {
    direct_cache
  } else {
    .somalign_project_samples(query$scaled_data, reference, chunk_size = chunk_size)
  }
  corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]
  corrected <- .somalign_project_samples(corrected_matrix, reference, chunk_size = chunk_size)
  correction_norm <- sqrt(rowSums(node_shifts[query$sample_unit, , drop = FALSE]^2))
  list(direct = direct, corrected = corrected, correction_norm = correction_norm)
}

.somalign_build_diagnostics <- function(transport, query, reference, node_shifts,
                                        projection, epsilon, rho_query,
                                        rho_ref, feature_weights = NULL) {
  ot <- transport$ot
  plan <- transport$plan
  row_mass <- transport$row_mass
  col_mass <- transport$col_mass
  direct <- projection$direct
  corrected <- projection$corrected
  mi_result <- .somalign_plan_mutual_information(
    plan, cost = transport$cost / transport$cost_scale
  )
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
      log_Z = ot$log_Z,
      anneal_schedule = ot$anneal_schedule,
      anneal_stage_info = ot$anneal_stage_info,
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
      max_col_mass_error = max(abs(col_mass - reference$node_masses)),
      mutual_information = mi_result$mutual_information
    ),
    nodes = .somalign_build_nodes_diag(query, transport, node_shifts, mi_result),
    projection = list(
      outside_direct_fraction = mean(direct$outside),
      outside_corrected_fraction = mean(corrected$outside)
    ),
    cost_metric = list(feature_weights = feature_weights)
  )
}

.somalign_build_nodes_diag <- function(query, transport, node_shifts, mi_result) {
  data.frame(
    query_node = seq_len(nrow(query$codebook)),
    query_mass = query$node_masses,
    transported_mass = transport$row_mass,
    match_fraction = transport$match_fraction,
    correction_allowed = attr(node_shifts, "correction_allowed"),
    correction_norm = sqrt(rowSums(node_shifts^2)),
    transport_entropy = mi_result$conditional_entropy
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

# Graph Laplacian L = D - A of a SOM's node neighbor graph, from 2-D grid
# coordinates. Two nodes are neighbors when their coordinate distance is 1
# (kohonen's unit lattice spacing, for both hexagonal and rectangular
# topologies); 1.01^2 absorbs floating-point rounding.
.somalign_som_laplacian <- function(grid) {
  if (is.null(grid) || is.null(grid$pts) || !is.matrix(grid$pts) || ncol(grid$pts) < 2L) {
    stop(
      "`laplacian_lambda > 0` requires the query SOM to have a kohonen grid ",
      "with 2-D node coordinates (grid$pts). Use a SOM trained with ",
      "kohonen::som() or kohonen::supersom().",
      call. = FALSE
    )
  }
  pts <- grid$pts
  d2 <- outer(pts[, 1L], pts[, 1L], "-")^2 + outer(pts[, 2L], pts[, 2L], "-")^2
  A <- (d2 > 0) & (d2 <= 1.01^2)
  storage.mode(A) <- "double"
  diag(rowSums(A)) - A
}

# Laplacian-regularised (Tikhonov) smoothing of node shifts: solves
# (W + lambda * L) x = W * shifts per feature column in one Cholesky
# factorisation, W = diag(node_masses) with correction_allowed == FALSE nodes
# zeroed out (they contribute no data term and are pulled toward their
# allowed neighbors' average). Disallowed-node rows of the *output* are then
# zeroed explicitly: .somalign_project_pair() applies node_shifts to every
# cell regardless of correction_allowed, so a disallowed node must keep an
# exact zero shift (matching current behavior) rather than receive an
# interpolated neighbor-average correction.
.somalign_smooth_shifts <- function(shifts, L, lambda, node_masses, correction_allowed) {
  m <- nrow(shifts)
  w <- node_masses
  if (!is.null(correction_allowed)) w[!correction_allowed] <- 0
  a_sys <- diag(w, nrow = m) + lambda * L
  a_sys <- a_sys + diag(m) * (.Machine$double.eps * max(abs(diag(a_sys))))
  rhs <- w * shifts
  ch <- tryCatch(chol(a_sys), error = function(e) NULL)
  out <- if (!is.null(ch)) chol2inv(ch) %*% rhs else solve(a_sys, rhs)
  if (!is.null(correction_allowed)) out[!correction_allowed, ] <- 0
  dimnames(out) <- dimnames(shifts)
  out
}

# Builds the shift_transform closure for laplacian_lambda > 0, or NULL
# (no-op) when lambda == 0 -- keeps somalign_fit()'s exported body short.
.somalign_make_laplacian_transform <- function(query, laplacian_lambda) {
  if (laplacian_lambda == 0) return(NULL)
  L <- .somalign_som_laplacian(query$som_query$grid)
  masses <- query$node_masses
  function(s) {
    ca <- attr(s, "correction_allowed")
    .somalign_smooth_shifts(s, L, laplacian_lambda, masses, correction_allowed = ca)
  }
}

# Composes an (optional) Laplacian smoother with an (optional) subspace
# projector as smooth -> project: the Laplacian operates over the full
# marker-space neighbor structure, and only then is the result restricted to
# the batch subspace V. Reversing the order would smooth an already
# rank-reduced field, which does not respect the marker-space geometry the
# Laplacian is defined over. Returns NULL when both inputs are NULL.
.somalign_compose_shift_transforms <- function(shift_fn_lap, shift_fn_sub) {
  if (!is.null(shift_fn_lap) && !is.null(shift_fn_sub)) {
    function(s) shift_fn_sub(shift_fn_lap(s))
  } else if (!is.null(shift_fn_lap)) {
    shift_fn_lap
  } else {
    shift_fn_sub
  }
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

#' Two-pass alignment decomposing correction into global and local components
#'
#' A two-stage variant of [somalign_fit()] that separates the batch correction
#' into a global shift (estimated at high regularisation) and a local residual
#' (refined at lower regularisation). This decomposition is most useful when
#' the batch effect has both a large uniform component and smaller
#' population-specific residuals.
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon_global Entropic regularisation for pass 1 (global). Higher
#'   values give a smoother, more diffuse transport plan that captures the
#'   mean batch shift while averaging out population-specific noise.
#'   Default `0.3`; should be larger than `epsilon_local`.
#' @param epsilon_local Entropic regularisation for pass 2 (local). Should
#'   be smaller than `epsilon_global` to refine residual node-level
#'   corrections. Default `0.1`.
#' @param rho_query Query-side unbalanced mass relaxation (both passes).
#' @param rho_ref Reference-side unbalanced mass relaxation (both passes).
#' @param solver Sinkhorn solver variant. See [somalign_fit()].
#' @param min_match_fraction Minimum transported fraction for corrections
#'   and label transfer. See [somalign_fit()].
#' @param confidence_threshold Minimum label confidence for transfer
#'   acceptance. See [somalign_fit()].
#' @param correction_min_mass Minimum transported mass for correction.
#'   See [somalign_fit()].
#' @param max_iter Maximum Sinkhorn iterations per pass.
#' @param tol Sinkhorn convergence tolerance.
#' @param chunk_size Integer. Number of samples per projection chunk.
#'   See [somalign_fit()].
#' @param label_guided Logical. When `TRUE`, applies a large cost penalty to
#'   node pairs with discordant dominant labels in both OT passes. See
#'   [somalign_fit()] for details.
#' @param variance_threshold Numeric in (0, 1]. Cumulative singular-value-squared
#'   fraction used to select the rank of the batch-subspace *diagnostic* stored in
#'   `$two_pass$batch_subspace`. Default `0.9`. Has no effect on the correction.
#' @param anneal_start,anneal_stages,anneal_factor Annealing-schedule tuning
#'   parameters, used only when `solver = "annealing"`. See [somalign_fit()].
#'   Applied independently within each pass (the warm start is within a pass,
#'   not across passes).
#'
#' @return A `somalign_fit` object with an additional `$two_pass` list
#'   containing `global_shift` (per-feature vector), `global_shift_norm`
#'   (Euclidean magnitude), `epsilon_global`, `epsilon_local`, and
#'   `batch_subspace` (a list with `V`, `rank`, `variance_explained` derived
#'   from the pass-1 correction field — **descriptive only**, not used for
#'   correction; may conflate batch effects with biology).
#'
#' @details
#' Pass 1 runs OT at `epsilon_global` between the original query codebook and
#' the reference codebook. The mass-weighted mean node shift across
#' correction-allowed nodes becomes the global shift `g`. Pass 2 runs OT at
#' `epsilon_local` between the globally shifted query codebook and the
#' reference, capturing residual population-specific displacements. The final
#' per-node correction is the residual plus `g`, so the total correction for
#' each cell equals its pass-2 barycentric target minus its original codebook
#' centroid.
#'
#' Direct projection (`old_som_unit`, `old_som_label`, `final_status`) is
#' computed from the original unshifted `query$scaled_data` and is unaffected
#' by the transport, preserving the transport-free primary result that
#' [somalign_fit()] guarantees.
#'
#' When the batch shift is predominantly global, `fit$two_pass$global_shift`
#' approximates the per-feature batch offset. When the shift is negligible,
#' the two-pass result converges toward a plain `somalign_fit()` at
#' `epsilon_local`.
#'
#' @seealso [somalign_fit()], [somalign_normalize()]
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit2 <- somalign_fit_two_pass(qry, ref)
#' @export
somalign_fit_two_pass <- function(query,
                                  reference,
                                  epsilon_global = 0.3,
                                  epsilon_local = 0.1,
                                  rho_query = 1,
                                  rho_ref = 1,
                                  solver = c("internal", "log_domain", "auto", "annealing"),
                                  min_match_fraction = 0.05,
                                  confidence_threshold = 0.6,
                                  correction_min_mass = 1e-8,
                                  max_iter = 1000,
                                  tol = 1e-7,
                                  chunk_size = 10000L,
                                  label_guided = FALSE,
                                  variance_threshold = 0.9,
                                  anneal_start = 10,
                                  anneal_stages = 10L,
                                  anneal_factor = NULL) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  .somalign_check_pos_scalar(epsilon_global, "epsilon_global")
  .somalign_check_pos_scalar(epsilon_local, "epsilon_local")
  .somalign_check_fit_params(rho_query, rho_ref, min_match_fraction,
                             confidence_threshold, correction_min_mass,
                             max_iter, tol, chunk_size, label_guided)
  .somalign_check_unit_scalar(variance_threshold, "variance_threshold")
  solver <- match.arg(solver, c("internal", "log_domain", "auto", "annealing"))
  if (identical(solver, "annealing"))
    .somalign_check_anneal_params(anneal_start, anneal_factor, anneal_stages)

  label_mask <- .somalign_resolve_label_mask(query, reference, label_guided)

  t1 <- .somalign_align_transport(query, reference, epsilon_global, rho_query,
                                   rho_ref, solver, max_iter, tol,
                                   label_mask = label_mask,
                                   anneal_start = anneal_start,
                                   anneal_factor = anneal_factor,
                                   anneal_stages = anneal_stages)
  ns1 <- .somalign_node_shifts(
    query_codebook    = query$codebook,
    reference_codebook = reference$codebook,
    correspondence    = t1$correspondence,
    row_mass          = t1$row_mass,
    match_fraction    = t1$match_fraction,
    min_match_fraction  = min_match_fraction,
    correction_min_mass = correction_min_mass
  )

  allowed1 <- attr(ns1, "correction_allowed")
  g <- if (any(allowed1)) {
    m <- query$node_masses[allowed1]
    colSums(ns1[allowed1, , drop = FALSE] * m) / sum(m)
  } else {
    stats::setNames(rep(0, ncol(query$codebook)), colnames(query$codebook))
  }

  g_mat <- matrix(g, nrow = nrow(query$codebook), ncol = length(g), byrow = TRUE)
  query2 <- query
  query2$codebook <- query$codebook + g_mat

  t2 <- .somalign_align_transport(query2, reference, epsilon_local, rho_query,
                                   rho_ref, solver, max_iter, tol,
                                   label_mask = label_mask,
                                   anneal_start = anneal_start,
                                   anneal_factor = anneal_factor,
                                   anneal_stages = anneal_stages)
  ns2 <- .somalign_node_shifts(
    query_codebook    = query2$codebook,
    reference_codebook = reference$codebook,
    correspondence    = t2$correspondence,
    row_mass          = t2$row_mass,
    match_fraction    = t2$match_fraction,
    min_match_fraction  = min_match_fraction,
    correction_min_mass = correction_min_mass
  )

  total_shifts <- ns2 + g_mat
  attr(total_shifts, "correction_allowed") <- attr(ns2, "correction_allowed")

  label_transfer <- .somalign_transfer_labels(
    correspondence     = t2$correspondence,
    label_prob         = reference$label_prob,
    match_fraction     = t2$match_fraction,
    min_match_fraction = min_match_fraction,
    confidence_threshold = confidence_threshold
  )

  projection <- .somalign_project_pair(query, reference, total_shifts, chunk_size)
  diagnostics <- .somalign_build_diagnostics(
    t2, query2, reference, total_shifts, projection, epsilon_local, rho_query, rho_ref
  )
  .somalign_fit_warnings(diagnostics)

  fit <- .somalign_new_fit(query, reference, t2, label_transfer, total_shifts,
                           projection, diagnostics)
  batch_sub <- if (any(allowed1)) {
    .somalign_subspace_svd(ns1[allowed1, , drop = FALSE], variance_threshold,
                           weights = query$node_masses[allowed1])
  } else {
    NULL
  }
  fit$two_pass <- list(
    global_shift      = g,
    global_shift_norm = sqrt(sum(g^2)),
    epsilon_global    = epsilon_global,
    epsilon_local     = epsilon_local,
    batch_subspace    = batch_sub
  )
  fit
}
