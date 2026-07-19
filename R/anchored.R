#' Align a query SOM to a reference SOM using anchor sample pairs
#'
#' A variant of [somalign_fit()] for the case where a set of samples has been
#' measured in **both** the old batch (reference space) and the new batch
#' (query space). These *anchor pairs* are used to build a per-node-pair
#' correspondence count matrix, which is subtracted from the normalized OT
#' cost before the Sinkhorn solve. This makes transport along anchor-supported
#' routes cheaper, biasing the OT plan toward pairings that are consistent
#' with the observed per-sample batch displacement — while still solving a
#' valid optimal transport problem over the full codebook.
#'
#' @param query A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param anchor_old Numeric matrix (n_anchors × p). Old-batch measurements
#'   of the anchor samples. Must be **raw (un-normalized) values in the same
#'   units and preprocessing pipeline as the data used to train `reference`**.
#'   Do not pre-center or pre-scale; this function applies `reference$center`
#'   and `reference$scale` internally. Also accepts a data frame of numeric
#'   columns.
#' @param anchor_new Numeric matrix (n_anchors × p). New-batch measurements
#'   of the **same** anchor samples. Must be **raw (un-normalized) values in
#'   the same units and preprocessing pipeline as `anchor_old`**. Do not
#'   pre-center or pre-scale; this function applies `reference$center` and
#'   `reference$scale` internally. Rows of `anchor_old` and `anchor_new` must
#'   correspond to the same biological units. Also accepts a data frame of
#'   numeric columns.
#' @param rho_anchor Non-negative scalar. Controls how strongly anchor pairs
#'   bias the OT cost. At `rho_anchor = 0` the result equals [somalign_fit()].
#'   Larger values reduce the effective cost for anchor-supported node pairs,
#'   concentrating the transport plan on those routes. Typical range: 0.5--3.
#'   Has no effect when `correction = "subspace"`.
#' @param epsilon Entropic regularisation strength (see [somalign_fit()]).
#' @param rho_query Query-side unbalanced mass relaxation.
#' @param rho_ref Reference-side unbalanced mass relaxation.
#' @param solver Sinkhorn solver variant. See [somalign_fit()].
#' @param min_match_fraction Minimum transported fraction for label transfer.
#' @param confidence_threshold Minimum top-label probability for label transfer.
#' @param correction_min_mass Minimum transported mass for a node correction.
#' @param max_iter Maximum Sinkhorn iterations.
#' @param tol Sinkhorn convergence tolerance.
#' @param chunk_size Integer. Samples projected per chunk. Default `10000L`.
#' @param correction Character. Correction strategy — one of
#'   `"cost_bonus"` (default), `"subspace"`, or `"both"`. See Details.
#' @param variance_threshold Numeric in (0, 1]. Cumulative singular-value-squared
#'   fraction for selecting the rank of the batch subspace. Default `0.9`
#'   (CellANOVA convention). Only used when `correction` is `"subspace"` or
#'   `"both"`.
#' @param anneal_start,anneal_stages,anneal_factor Annealing-schedule tuning
#'   parameters, used only when `solver = "annealing"`. See [somalign_fit()].
#' @param feature_weights Either `NULL` (default, squared-Euclidean cost), a
#'   named non-negative numeric vector of explicit per-feature weights (see
#'   [somalign_fit()]), or the string `"anchor"` -- auto-estimates weights
#'   from the anchor displacement matrix `D` via
#'   \eqn{w_f = 1 / (\mathrm{var}(D_{\cdot f}) + \delta)}, mean-normalised.
#'   Markers that vary most across the batch (large `var(D[, f])`, i.e.
#'   batch-driven) get low weight and are cheap to transport; markers stable
#'   across the batch get high weight and are expensive to transport,
#'   preserving biology. The resolved vector is stored in
#'   `fit$anchors$feature_weights` and
#'   `fit$diagnostics$cost_metric$feature_weights`. Composes independently
#'   with `correction`: the weights reshape the cost geometry, while
#'   `rho_anchor`/`correction` bias routing -- both act on the same
#'   underlying transport problem without conflict.
#' @param laplacian_lambda Non-negative scalar. Graph-Laplacian smoothing of
#'   the node-shift field; see [somalign_fit()]. When `correction` is
#'   `"subspace"` or `"both"`, smoothing is applied *before* the subspace
#'   projection (smooth in full marker space, then restrict to the batch
#'   subspace `V`) so the Laplacian neighbor structure is respected. Default
#'   `0` (no smoothing).
#'
#' @details
#' **Correction modes.** Three strategies are available via the `correction`
#' argument.
#'
#' - `"cost_bonus"` (default, current behaviour): the anchor count matrix
#'   biases the OT cost so anchor-supported node pairs are cheaper; the
#'   resulting node shifts are applied to the full feature space.
#'
#' - `"subspace"`: a batch subspace \eqn{V_{\text{batch}}} is estimated by
#'   SVD of the anchor displacement matrix
#'   \eqn{D = X_{\text{old}} - X_{\text{new}}} (n_anchors × p). Because each
#'   row of \eqn{D} is a *same-biological-unit* before–after measurement, the
#'   dominant singular vectors isolate the true batch direction. Node shifts
#'   from a *plain* OT solve (no cost bonus) are then projected onto
#'   \eqn{V_{\text{batch}}}: only the batch-direction component is applied.
#'   Biological variation orthogonal to \eqn{V_{\text{batch}}} is preserved.
#'   A synthetic validation shows the orthogonal component survives at
#'   ~99.7% (1.496 vs ideal 1.500) while `"cost_bonus"` erases it.
#'   The rank \eqn{r} is the smallest index where the cumulative squared
#'   singular values reach `variance_threshold` (default 0.9).
#'   \eqn{D} is **not centred** — the mean batch direction is the dominant
#'   structure we want to capture.
#'
#' - `"both"`: applies the cost bonus to the OT solve *and* restricts the
#'   resulting shifts to \eqn{V_{\text{batch}}}.
#'
#' `"subspace"` and `"both"` expose `fit$anchors$batch_subspace` (a list with
#' `V`, `rank`, `variance_explained`). `"cost_bonus"` sets this to `NULL`.
#'
#' **Topology preservation and epsilon.** Empirically, the primary driver of
#' topology/structure damage from batch correction is `epsilon`, not
#' `rho_anchor`. Higher epsilon blurs the transport plan across a wider
#' neighbourhood, causing the corrected codebook to collapse biologically
#' distinct populations (H0 component merging). Subspace-restricted modes
#' (`"subspace"` or `"both"`) substantially reduce merging at any given
#' epsilon because shifts are confined to the batch-variation subspace, leaving
#' orthogonal biological variation intact. As a result, choosing epsilon
#' involves a genuine trade-off: higher epsilon is more numerically stable for
#' the Sinkhorn solver, but lower epsilon preserves more topology. Before
#' committing to an epsilon, run
#' `somalign_epsilon_sweep(..., topology = TRUE)` alongside
#' [somalign_select_epsilon()] and inspect both the phase-transition criterion
#' and the `biggest_merge_mass_frac` column -- the two criteria can disagree,
#' especially at small epsilon near numerical instability.
#'
#' **Cost modification.** Let \eqn{C} be the M×K codebook distance matrix
#' normalised by its median positive entry (as in [somalign_fit()]). Each
#' anchor pair is projected onto both codebooks to build a count matrix
#' \eqn{A} where \eqn{A_{kl}} is the number of anchor pairs mapping to query
#' node \eqn{k} and reference node \eqn{l}. The query SOM was trained on
#' new-batch data, so the *new-batch* anchor measurement is projected onto the
#' query codebook to identify query node \eqn{k}; the reference SOM was trained
#' on old-batch data, so the *old-batch* anchor measurement is projected onto
#' the reference codebook to identify reference node \eqn{l}.
#' The modified cost is
#' \deqn{\tilde{C}_{kl} = \max\!\bigl(C_{kl} - \rho_{\mathrm{anchor}} \cdot
#'   A_{kl} / n_{\mathrm{anchors}},\; 0\bigr).}
#' Pairs with many anchor observations get cost reduced toward zero (free
#' transport), while uncovered pairs retain their original cost. Non-negativity
#' is enforced by the \eqn{\max(\cdot, 0)} clamp.
#'
#' **Clamp behaviour at large `rho_anchor`.** When the anchor bonus exceeds
#' \eqn{C_{kl}}, the effective cost is clamped to zero. All such pairs then
#' have identical effective cost and the transport mass among them is determined
#' by entropic regularisation alone rather than by relative anchor counts.
#' The clamp is required to keep costs non-negative; at very large `rho_anchor`
#' the plan for anchor-covered pairs becomes more entropic, not more
#' concentrated. A practical upper bound is `rho_anchor * max(A) / n_anchors
#' <= 1`, i.e., even the most-supported pair reduces cost by at most one
#' median squared-distance unit.
#'
#' **Fallback for uncovered nodes.** Query nodes with no anchor samples retain
#' their original pairwise costs, so the transport plan for those nodes is
#' determined entirely by the OT objective — the same as [somalign_fit()].
#' Inspect `$anchors$coverage_fraction` to see what fraction of query nodes
#' had at least one anchor pair.
#'
#' **Return value.** The object has class `c("somalign_anchored_fit",
#' "somalign_fit")`, so all downstream functions that accept a
#' `somalign_fit` object ([somalign_results()], [somalign_diagnostics()])
#' work unchanged. An additional `$anchors` list element is attached:
#' \describe{
#'   \item{`n_anchors`}{Number of anchor pairs supplied.}
#'   \item{`rho_anchor`}{The value of `rho_anchor` used.}
#'   \item{`correction`}{The correction mode: `"cost_bonus"`, `"subspace"`, or `"both"`.}
#'   \item{`nodes_covered`}{Number of query nodes with ≥ 1 anchor pair.}
#'   \item{`coverage_fraction`}{`nodes_covered / nrow(query$codebook)`.}
#'   \item{`batch_subspace`}{For `"subspace"` and `"both"` modes: a list with
#'     `V` (p × rank matrix), `rank` (integer), and `variance_explained`
#'     (cumulative variance at the selected rank). `NULL` for `"cost_bonus"`.}
#'   \item{`displacements`}{The scaled anchor displacement matrix
#'     \eqn{D = X_{\text{old,scaled}} - X_{\text{new,scaled}}}
#'     (n_anchors × p), always stored regardless of `correction` mode. Used by
#'     [somalign_subspace_sensitivity()] and [somalign_exclusion_test()].}
#' }
#'
#' @return A `somalign_anchored_fit` object (also inherits `somalign_fit`).
#' @note At small `epsilon` with high anchor coverage the anchor bonus zeros out
#'   many entries of the normalised cost matrix, which sharpens the Sinkhorn
#'   kernel and can drive the remaining entries toward numerical underflow. If
#'   the solver warns about kernel underflow, pass `solver = "log_domain"` or
#'   `solver = "annealing"`, both of which work in log-potential space and
#'   avoid the issue.
#' @seealso [somalign_fit()] for the unanchored variant.
#' @examples
#' set.seed(1)
#' p   <- 3L
#' mat <- rbind(
#'   matrix(rnorm(20 * p, mean = -2), ncol = p),
#'   matrix(rnorm(20 * p, mean =  2), ncol = p)
#' )
#' colnames(mat) <- paste0("F", seq_len(p))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' shifted <- mat + 0.5
#' qry <- somalign_query(shifted, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' # Use 10 samples as anchors measured in both batches
#' anc_idx <- 1:10
#' fit <- somalign_fit_anchored(qry, ref,
#'                               anchor_old = mat[anc_idx, , drop = FALSE],
#'                               anchor_new = shifted[anc_idx, , drop = FALSE],
#'                               rho_anchor = 1)
#' fit$anchors
#' @export
somalign_fit_anchored <- function(query,
                                   reference,
                                   anchor_old,
                                   anchor_new,
                                   rho_anchor = 1.0,
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
                                   correction = c("cost_bonus", "subspace", "both"),
                                   variance_threshold = 0.9,
                                   anneal_start = 10,
                                   anneal_stages = 10L,
                                   anneal_factor = NULL,
                                   feature_weights = NULL,
                                   laplacian_lambda = 0) {
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  solver <- match.arg(solver, c("internal", "log_domain", "auto", "annealing"))
  correction <- match.arg(correction, c("cost_bonus", "subspace", "both"))
  if (!is.numeric(rho_anchor) || length(rho_anchor) != 1L ||
      !is.finite(rho_anchor) || rho_anchor < 0) {
    stop("`rho_anchor` must be a non-negative finite scalar.", call. = FALSE)
  }
  .somalign_check_pos_scalar(epsilon, "epsilon")
  .somalign_check_fit_params(rho_query, rho_ref, min_match_fraction,
                             confidence_threshold, correction_min_mass,
                             max_iter, tol, chunk_size)
  .somalign_check_unit_scalar(variance_threshold, "variance_threshold")
  if (identical(solver, "annealing"))
    .somalign_check_anneal_params(anneal_start, anneal_factor, anneal_stages)
  feature_weights <- .somalign_check_feature_weights(feature_weights, colnames(query$codebook))
  .somalign_check_nonneg_scalar(laplacian_lambda, "laplacian_lambda")
  anchors_scaled <- .somalign_validate_anchors(anchor_old, anchor_new, reference)
  .somalign_anchored_dispatch(
    query, reference, anchors_scaled, rho_anchor, epsilon, rho_query, rho_ref,
    solver, min_match_fraction, confidence_threshold, correction_min_mass,
    chunk_size, max_iter, tol, correction, variance_threshold,
    anneal_start, anneal_factor, anneal_stages, feature_weights, laplacian_lambda
  )
}

.somalign_anchored_dispatch <- function(query, reference, anchors_scaled,
                                         rho_anchor, epsilon, rho_query, rho_ref,
                                         solver, min_match_fraction,
                                         confidence_threshold, correction_min_mass,
                                         chunk_size, max_iter, tol,
                                         correction, variance_threshold,
                                         anneal_start = 10, anneal_factor = NULL,
                                         anneal_stages = 10L,
                                         feature_weights = NULL,
                                         laplacian_lambda = 0) {
  use_bonus    <- correction %in% c("cost_bonus", "both") && rho_anchor > 0
  use_subspace <- correction %in% c("subspace", "both")
  if (rho_anchor == 0 && correction %in% c("cost_bonus", "both")) {
    message("`rho_anchor = 0`: the anchor cost bonus is inactive. ",
            if (correction == "cost_bonus")
              "Use `somalign_fit()` for equivalent results."
            else
              "Only the subspace correction is applied.")
  }
  cb <- if (use_bonus) {
    .somalign_anchor_cost_bonus(anchors_scaled$anchor_old_scaled,
                                anchors_scaled$anchor_new_scaled,
                                query$codebook, reference$codebook,
                                rho_anchor, chunk_size)
  } else {
    cov <- .somalign_anchor_coverage(anchors_scaled$anchor_new_scaled,
                                     query$codebook, chunk_size)
    list(bonus = NULL, nodes_covered = cov$nodes_covered,
         coverage_fraction = cov$coverage_fraction)
  }
  # Scaled anchor displacements: always computed and stored (regardless of
  # correction mode) so downstream diagnostics (somalign_subspace_sensitivity,
  # somalign_exclusion_test) can reuse them without re-deriving from raw
  # anchors. The batch subspace itself is only estimated when actually needed
  # for correction.
  d_scaled <- anchors_scaled$anchor_old_scaled - anchors_scaled$anchor_new_scaled
  batch_sub <- if (use_subspace) {
    .somalign_subspace_svd(d_scaled, variance_threshold)
  } else { NULL }
  shift_fn <- if (use_subspace) {
    V <- batch_sub$V
    function(s) s %*% V %*% t(V)
  } else { NULL }
  shift_fn_lap <- .somalign_make_laplacian_transform(query, laplacian_lambda)
  shift_fn <- .somalign_compose_shift_transforms(shift_fn_lap, shift_fn)
  fw <- if (identical(feature_weights, "anchor")) {
    .somalign_anchor_feature_weights(d_scaled)
  } else {
    feature_weights
  }
  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol,
    cost_bonus = cb$bonus,
    anneal_start = anneal_start, anneal_factor = anneal_factor,
    anneal_stages = anneal_stages, feature_weights = fw
  )
  fit <- .somalign_finish_fit(
    query, reference, transport, min_match_fraction, confidence_threshold,
    correction_min_mass, chunk_size, epsilon, rho_query, rho_ref,
    shift_transform = shift_fn, feature_weights = fw,
    anchors = list(
      n_anchors         = nrow(anchors_scaled$anchor_old_scaled),
      rho_anchor        = rho_anchor,
      correction        = correction,
      nodes_covered     = cb$nodes_covered,
      coverage_fraction = cb$coverage_fraction,
      batch_subspace    = batch_sub,
      variance_threshold = variance_threshold,
      displacements     = d_scaled,
      feature_weights   = fw
    )
  )
  class(fit) <- c("somalign_anchored_fit", "somalign_fit")
  fit
}

.somalign_validate_anchors <- function(anchor_old, anchor_new, reference) {
  anchor_old <- .somalign_as_matrix(anchor_old, what = "anchor_old")
  anchor_new <- .somalign_as_matrix(anchor_new, what = "anchor_new")

  if (nrow(anchor_old) == 0L)
    stop("`anchor_old` and `anchor_new` must have at least one row.", call. = FALSE)
  if (nrow(anchor_old) != nrow(anchor_new))
    stop("`anchor_old` and `anchor_new` must have the same number of rows.",
         call. = FALSE)

  features <- reference$features
  if (!is.null(features)) {
    anchor_old <- .somalign_select_features(anchor_old, features, what = "anchor_old")
    anchor_new <- .somalign_select_features(anchor_new, features, what = "anchor_new")
  } else {
    ref_names <- names(reference$center)
    if (!is.null(ref_names)) {
      anchor_old <- .somalign_select_features(anchor_old, ref_names, what = "anchor_old")
      anchor_new <- .somalign_select_features(anchor_new, ref_names, what = "anchor_new")
    } else {
      if (ncol(anchor_old) != length(reference$center))
        stop("`anchor_old` must have the same number of columns as reference features.",
             call. = FALSE)
      if (ncol(anchor_new) != length(reference$center))
        stop("`anchor_new` must have the same number of columns as reference features.",
             call. = FALSE)
    }
  }

  .somalign_validate_finite(anchor_old, what = "anchor_old")
  .somalign_validate_finite(anchor_new, what = "anchor_new")

  list(
    anchor_old_scaled = .somalign_scale_matrix(anchor_old, reference$center, reference$scale),
    anchor_new_scaled = .somalign_scale_matrix(anchor_new, reference$center, reference$scale)
  )
}

.somalign_anchor_cost_bonus <- function(anchor_old_scaled, anchor_new_scaled,
                                         query_codebook, reference_codebook,
                                         rho_anchor, chunk_size) {
  M <- nrow(query_codebook)
  K <- nrow(reference_codebook)
  n_anchors <- nrow(anchor_old_scaled)

  # The query SOM was trained on new-batch data and the reference SOM on
  # old-batch data, so the new-batch anchor identifies the query node (row) and
  # the old-batch anchor identifies the reference node (column).
  query_units <- .somalign_nearest_code_chunked(
    anchor_new_scaled, query_codebook, chunk_size = chunk_size
  )$unit
  ref_units <- .somalign_nearest_code_chunked(
    anchor_old_scaled, reference_codebook, chunk_size = chunk_size
  )$unit

  # Build M × K count matrix (column-major linear indexing)
  lin_idx <- (ref_units - 1L) * M + query_units
  A <- matrix(tabulate(lin_idx, nbins = M * K), nrow = M, ncol = K)

  nodes_covered    <- sum(rowSums(A) > 0L)
  coverage_fraction <- nodes_covered / M

  bonus <- rho_anchor * (A / n_anchors)

  list(
    bonus             = bonus,
    nodes_covered     = nodes_covered,
    coverage_fraction = coverage_fraction
  )
}

# Query-node coverage of the anchor set, independent of the cost bonus. Used by
# the subspace-only path where no bonus matrix is built. The new-batch anchor
# identifies the query node (the query SOM was trained on new-batch data).
.somalign_anchor_coverage <- function(anchor_new_scaled, query_codebook,
                                       chunk_size) {
  M <- nrow(query_codebook)
  query_units <- .somalign_nearest_code_chunked(
    anchor_new_scaled, query_codebook, chunk_size = chunk_size
  )$unit
  nodes_covered <- length(unique(query_units))
  list(nodes_covered = nodes_covered, coverage_fraction = nodes_covered / M)
}

# ---------------------------------------------------------------------------
# Anchor exclusion-restriction test (Sargan-Hansen analog)
# ---------------------------------------------------------------------------

# R = D (I_p - V V^T): projection of each row of D onto the orthogonal
# complement of span(V). Avoids materialising the p x p identity matrix.
.somalign_orthogonal_residual <- function(D, V) {
  D - D %*% V %*% t(V)
}

# Column-independent permutation null for the exclusion test: shuffle each
# COLUMN (feature) of the orthogonal residual R_obs independently across
# anchor rows. This preserves each feature's own marginal distribution (and
# hence the residual's total variance / trace) but destroys any cross-feature
# correlation -- i.e., any consistent multi-feature *direction* shared across
# anchors. The leading singular value is large only when R_obs's rows point
# in a coherent, shared direction (a structured, biology-like signal); if
# R_obs is unstructured per-feature noise, permuting each column
# independently barely changes the singular-value spectrum.
#
# NOTE ON A REJECTED ALTERNATIVE: permuting *rows* of R (or of D) is a
# mathematical no-op for this purpose. For any orthogonal row transform P
# (a row permutation matrix or a diagonal +-1 sign-flip matrix is orthogonal),
# (PR)^T (PR) = R^T P^T P R = R^T R, so every singular value of R is EXACTLY
# invariant under row permutation or row sign-flipping. A null built that way
# would silently have zero power. Per-column permutation is not an orthogonal
# row transform (it does not act as a single matrix multiplication on R from
# the left) and genuinely alters R^T R by decorrelating features, which is
# exactly what a valid null for "is there a coherent cross-feature direction"
# requires.
.somalign_permutation_null <- function(R_obs, n_perm) {
  n <- nrow(R_obs)
  p <- ncol(R_obs)
  sv_null <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    r_perm <- R_obs
    for (j in seq_len(p)) r_perm[, j] <- R_obs[sample.int(n), j]
    sv_null[b] <- svd(r_perm, nu = 0L, nv = 0L)$d[1L]
  }
  sv_null
}

#' Anchor exclusion-restriction test
#'
#' Permutation test of whether the anchor displacement matrix carries
#' *coherent, structured* signal orthogonal to the estimated batch subspace --
#' an overidentification-style check of the instrumental-variable assumption
#' underlying [somalign_fit_anchored()]'s `correction = "subspace"` mode:
#' anchors should isolate the batch direction, not biology.
#'
#' @section Role: A **diagnostic validating the correction path's** assumption,
#'   not a label-transfer diagnostic. A `fail` verdict means the corrected
#'   coordinates may mix biology into the batch subspace; it does not impugn the
#'   transferred labels, which are computed from the transport plan alone.
#'
#' @param fit A `somalign_anchored_fit` object from
#'   `somalign_fit_anchored(..., correction = "subspace")` or `"both"`.
#' @param n_perm Positive integer. Number of permutation replicates. Default `999L`.
#' @param seed Integer or `NULL`. RNG seed for reproducibility (restored on
#'   exit; does not leak into the caller's session). Default `1L`; `NULL`
#'   disables seeding.
#'
#' @details
#' The statistic is the leading singular value of the orthogonal residual
#' \eqn{R = D (I - V V^\top)}, where \eqn{D} is the scaled anchor displacement
#' matrix (`fit$anchors$displacements`) and \eqn{V} is the estimated batch
#' subspace (`fit$anchors$batch_subspace$V`). The null distribution is
#' generated by permuting each feature (column) of \eqn{R} independently
#' across anchors, which preserves each feature's own variance but destroys
#' any coherent cross-feature direction. A small p-value means \eqn{R}'s
#' features move together more than chance recombination would produce --
#' i.e., the anchors carry a real, structured direction that `correction =
#' "subspace"` is *not* removing (a batch-direction violation, or biology
#' leaking through the anchors). A large p-value supports the assumption that
#' anchors capture only batch structure.
#'
#' Row permutation of \eqn{R} or \eqn{D} would *not* work here: row
#' permutation is an orthogonal transformation and leaves every singular
#' value of a matrix exactly invariant, so it cannot detect anything.
#' Column-wise permutation is required because it is not an orthogonal row
#' transform and genuinely decorrelates features.
#'
#' The test has no power when `variance_threshold = 1` (V spans the full
#' feature space, so R is identically zero), and low power when
#' `n_anchors < 3 * rank`, or when the true structure is fully absorbed into
#' \eqn{V} (inspect `fit$anchors$batch_subspace$rank` and
#' `variance_explained`).
#'
#' @return A list of class `somalign_exclusion_test` with `sv_observed`,
#'   `sv_null`, `p_value`, `null_quantiles`, `relative_stat`, `rank_used`,
#'   `n_anchors`, `n_features`, `verdict` (`"pass"` if `p > 0.1`, `"warn"` if
#'   in `[0.05, 0.1]`, `"fail"` if `< 0.05`), and `feature_residual_norm`
#'   (per-feature `||R[, j]||`, identifying which markers drive a violation).
#' @seealso [somalign_fit_anchored()]
#' @examples
#' set.seed(1)
#' p <- 3L
#' mat <- rbind(
#'   matrix(rnorm(20 * p, mean = -2), ncol = p),
#'   matrix(rnorm(20 * p, mean =  2), ncol = p)
#' )
#' colnames(mat) <- paste0("F", seq_len(p))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' shifted <- mat + 0.5
#' qry <- somalign_query(shifted, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' anc_idx <- 1:10
#' fit <- somalign_fit_anchored(qry, ref,
#'                               anchor_old = mat[anc_idx, , drop = FALSE],
#'                               anchor_new = shifted[anc_idx, , drop = FALSE],
#'                               rho_anchor = 1, correction = "subspace")
#' somalign_exclusion_test(fit, n_perm = 199L)
#' @export
somalign_exclusion_test <- function(fit, n_perm = 999L, seed = 1L) {
  .somalign_check_exclusion_test_args(fit, n_perm, seed)
  bs <- fit$anchors$batch_subspace
  D <- fit$anchors$displacements
  V <- bs$V
  n <- nrow(D)
  r <- bs$rank
  if (n < r + 2L) {
    warning(sprintf(
      "n_anchors (%d) < rank + 2 (%d): the exclusion test has essentially no power. ",
      n, r + 2L), "Recommend >= 3 * rank anchor pairs.", call. = FALSE)
  }

  r_obs <- .somalign_orthogonal_residual(D, V)
  sv_obs <- svd(r_obs, nu = 0L, nv = 0L)$d[1L]
  sv_null <- .somalign_seeded_permutation_null(r_obs, n_perm, seed)

  p_value <- mean(sv_null >= sv_obs)
  nq <- stats::quantile(sv_null, probs = c(0.025, 0.5, 0.975), names = TRUE)
  relative_stat <- if (nq[["50%"]] > 0) sv_obs / nq[["50%"]] else NA_real_
  feature_residual_norm <- sqrt(colSums(r_obs^2))
  names(feature_residual_norm) <- colnames(D)
  verdict <- if (p_value > 0.10) "pass" else if (p_value >= 0.05) "warn" else "fail"

  structure(
    list(sv_observed = sv_obs, sv_null = sv_null, p_value = p_value,
         null_quantiles = nq, relative_stat = relative_stat,
         rank_used = r, n_anchors = n, n_features = ncol(D),
         verdict = verdict, feature_residual_norm = feature_residual_norm),
    class = "somalign_exclusion_test"
  )
}

.somalign_check_exclusion_test_args <- function(fit, n_perm, seed) {
  if (!inherits(fit, "somalign_anchored_fit"))
    stop("`fit` must be a somalign_anchored_fit object.", call. = FALSE)
  if (!fit$anchors$correction %in% c("subspace", "both"))
    stop("`somalign_exclusion_test` requires a fit with correction = 'subspace' or 'both'.",
         call. = FALSE)
  if (is.null(fit$anchors$batch_subspace))
    stop("No batch subspace found in `fit$anchors$batch_subspace`.", call. = FALSE)
  if (is.null(fit$anchors$displacements))
    stop("`fit$anchors$displacements` is NULL. Refit with a version of somalign ",
         "that stores anchor displacement matrices.", call. = FALSE)
  .somalign_check_pos_int(n_perm, "n_perm")
  if (!is.null(seed)) .somalign_check_pos_int(seed, "seed")
  invisible(NULL)
}

# Runs the permutation null with a seed that is local to this call: the
# caller's global RNG state is saved before seeding and restored on exit, so
# somalign_exclusion_test() does not leak RNG state into the session.
.somalign_seeded_permutation_null <- function(r_obs, n_perm, seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  .somalign_permutation_null(r_obs, n_perm)
}
