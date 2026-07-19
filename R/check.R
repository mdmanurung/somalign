#' Check query-reference codebook alignment before cell-level computation
#'
#' Compares a query SOM codebook (already in reference-scaled space) against the
#' reference codebook across three diagnostic dimensions.  The check is
#' O(nodes\eqn{^2 \times} p) — milliseconds on a 900-node SOM — and is designed
#' to surface distribution mismatches before any per-cell work begins.
#'
#' @details
#' \strong{Range overlap (per feature):} Does the query codebook's value range
#' for each marker intersect the reference codebook's range?  Zero overlap
#' means every query node sits entirely outside the reference for that marker:
#' a critical failure.  Less than 50\% overlap is a warning.
#'
#' \strong{Mass-weighted centroid drift (per feature):} The mass-weighted mean
#' of the query codebook minus that of the reference, expressed in units of the
#' reference codebook standard deviation.  Drift > 3 SDs flags a global batch
#' shift that the OT plan may not be able to absorb.
#'
#' \strong{Transport coverage (cost matrix preview):} Fraction of
#' query-reference codebook pairs whose normalised squared distance falls within
#' \eqn{3\varepsilon}.  Pairs outside this band contribute negligible weight to
#' the Sinkhorn kernel.  If fewer than 1\% of pairs are within
#' \eqn{3\varepsilon}, the transport plan will be near-singular and most query
#' mass may be destroyed; re-check coordinate alignment or raise epsilon.
#'
#' @param query_codebook Numeric matrix of query SOM codebook vectors in
#'   reference-scaled coordinate space (nodes \eqn{\times} features).  Column
#'   names must include all features in \code{reference$features}.
#' @param reference A \code{somalign_reference} object.
#' @param query_masses Optional numeric vector of query node masses (length
#'   \code{nrow(query_codebook)}).  Used for the mass-weighted centroid check.
#'   When \code{NULL}, uniform weights are assumed.
#' @param epsilon The OT regularisation parameter that will be passed to
#'   \code{\link{somalign_fit}()}.  Used to contextualise the cost matrix
#'   coverage check.  Default \code{0.1}.
#' @param stop_if_critical If \code{TRUE} (default), throw an error when any
#'   feature has zero range overlap.  Set to \code{FALSE} to emit a warning
#'   instead and return the diagnostics.
#'
#' @return A \code{somalign_codebook_check} list (returned invisibly) with:
#' \describe{
#'   \item{\code{per_feature}}{Data frame: one row per feature with
#'     \code{ref_min}, \code{ref_max}, \code{query_min}, \code{query_max},
#'     \code{overlap_fraction}, \code{centroid_drift}, \code{centroid_drift_sd},
#'     and \code{flag} (\code{"ok"} / \code{"warning"} / \code{"critical"}).}
#'   \item{\code{cost_summary}}{Named numeric vector:
#'     \code{median_cost} (raw median pairwise squared distance),
#'     \code{p95_cost}, \code{cost_scale} (normalisation factor used by
#'     \code{somalign_fit}), and \code{fraction_near_eps} (fraction of pairs
#'     within \eqn{3\varepsilon} of the normalised cost).}
#'   \item{\code{n_critical_features}}{Number of features with zero overlap.}
#'   \item{\code{n_warning_features}}{Number of features flagged as warning.}
#'   \item{\code{verdict}}{\code{"pass"}, \code{"warning"}, or
#'     \code{"critical"}.}
#' }
#'
#' @seealso [somalign_fit()], [somalign_query()]
#'
#' @export
somalign_check_codebook_alignment <- function(query_codebook,
                                              reference,
                                              query_masses  = NULL,
                                              epsilon       = 0.1,
                                              stop_if_critical = TRUE) {
  .somalign_validate_check_args(query_codebook, reference, epsilon,
                                stop_if_critical)
  query_codebook <- as.matrix(query_codebook)
  storage.mode(query_codebook) <- "double"

  features <- reference$features
  missing_f <- setdiff(features, colnames(query_codebook))
  if (length(missing_f) > 0L) {
    stop(
      "`query_codebook` is missing ", length(missing_f), " feature(s): ",
      paste(missing_f, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!identical(colnames(query_codebook), features)) {
    query_codebook <- query_codebook[, features, drop = FALSE]
  }
  ref_cb <- reference$codebook[, features, drop = FALSE]

  n_q <- nrow(query_codebook)
  n_r <- nrow(ref_cb)

  # --- normalise weights ---
  r_wt <- reference$node_masses
  if (is.null(query_masses)) {
    q_wt <- rep(1 / n_q, n_q)
  } else {
    q_wt <- as.numeric(query_masses)
    if (length(q_wt) != n_q) {
      stop(
        "`query_masses` length (", length(q_wt),
        ") must equal nrow(query_codebook) (", n_q, ").",
        call. = FALSE
      )
    }
    s <- sum(q_wt)
    if (!is.finite(s) || s <= 0) stop("`query_masses` must have a positive finite sum.", call. = FALSE)
    q_wt <- q_wt / s
  }

  # =========================================================================
  # Check 1 — per-feature range overlap
  # =========================================================================
  ref_min <- apply(ref_cb, 2L, min)
  ref_max <- apply(ref_cb, 2L, max)
  q_min   <- apply(query_codebook, 2L, min)
  q_max   <- apply(query_codebook, 2L, max)

  ref_range    <- ref_max - ref_min
  overlap      <- pmax(0, pmin(q_max, ref_max) - pmax(q_min, ref_min))
  overlap_frac <- ifelse(ref_range > 0, overlap / ref_range, 1.0)

  # =========================================================================
  # Check 2 — mass-weighted centroid drift
  # =========================================================================
  ref_centroid <- as.vector(r_wt %*% ref_cb)
  q_centroid   <- as.vector(q_wt %*% query_codebook)
  drift        <- q_centroid - ref_centroid

  occupied <- r_wt > 0
  ref_sd   <- apply(ref_cb[occupied, , drop = FALSE], 2L, stats::sd)
  ref_sd   <- pmax(ref_sd, 1e-8)
  drift_sd <- drift / ref_sd

  # per-feature verdict
  flag <- ifelse(
    overlap_frac == 0, "critical",
    ifelse(overlap_frac < 0.5 | abs(drift_sd) > 3, "warning", "ok")
  )

  per_feature <- data.frame(
    feature           = features,
    ref_min           = ref_min,
    ref_max           = ref_max,
    query_min         = q_min,
    query_max         = q_max,
    overlap_fraction  = overlap_frac,
    centroid_drift    = drift,
    centroid_drift_sd = drift_sd,
    flag              = flag,
    stringsAsFactors  = FALSE,
    row.names         = NULL
  )

  # =========================================================================
  # Check 3 — cost matrix preview (O(nodes^2 * p))
  # =========================================================================
  # Pairwise squared Euclidean distances between query and reference codebook
  d2 <- outer(rowSums(query_codebook^2), rowSums(ref_cb^2), "+") -
    2 * tcrossprod(query_codebook, ref_cb)
  d2 <- pmax(d2, 0)  # clamp floating-point negatives

  pos_vals   <- d2[d2 > 0]
  cost_scale <- if (length(pos_vals) > 0L) stats::median(pos_vals) else 1.0
  d2_norm    <- d2 / cost_scale

  median_cost    <- stats::median(d2)
  p95_cost       <- as.numeric(stats::quantile(d2, 0.95))
  # fraction of pairs with non-negligible Sinkhorn kernel weight exp(-d2_norm/eps)
  fraction_near  <- mean(d2_norm < 3 * epsilon)

  cost_summary <- c(
    median_cost       = median_cost,
    p95_cost          = p95_cost,
    cost_scale        = cost_scale,
    fraction_near_eps = fraction_near
  )

  # =========================================================================
  # Overall verdict
  # =========================================================================
  n_crit <- sum(flag == "critical")
  n_warn <- sum(flag == "warning")

  cost_critical <- fraction_near < 0.01
  verdict <- if (n_crit > 0) "critical" else if (n_warn > 0 || cost_critical) "warning" else "pass"

  # =========================================================================
  # Emit diagnostics
  # =========================================================================
  if (n_crit > 0L) {
    crit_feats <- per_feature$feature[per_feature$flag == "critical"]
    msg <- paste0(
      "somalign_check_codebook_alignment: ",
      n_crit, " feature(s) have zero range overlap between the query and ",
      "reference codebooks: ", paste(crit_feats, collapse = ", "), ". ",
      "Verify that the query codebook has been transformed into the same ",
      "coordinate space as the reference (center, scale, and any ",
      "winsorization bounds must be consistent between cohorts)."
    )
    if (isTRUE(stop_if_critical)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  if (n_warn > 0L) {
    warn_feats <- per_feature$feature[per_feature$flag == "warning"]
    message(
      "somalign_check_codebook_alignment: ",
      n_warn, " feature(s) show partial mismatch ",
      "(< 50% range overlap or centroid drift > 3 reference SDs): ",
      paste(warn_feats, collapse = ", "), ". ",
      "Label transfer accuracy may be reduced for these markers."
    )
  }

  if (cost_critical) {
    message(
      "somalign_check_codebook_alignment: only ",
      round(100 * fraction_near, 1L),
      "% of query-reference codebook pairs fall within 3\u03b5 ",
      "(epsilon = ", epsilon, "). ",
      "The Sinkhorn kernel will be near-singular and most query mass may be ",
      "destroyed. Consider re-checking coordinate alignment or raising epsilon."
    )
  }

  invisible(structure(
    list(
      per_feature         = per_feature,
      cost_summary        = cost_summary,
      n_critical_features = n_crit,
      n_warning_features  = n_warn,
      verdict             = verdict
    ),
    class = "somalign_codebook_check"
  ))
}

#' @method print somalign_codebook_check
#' @export
print.somalign_codebook_check <- function(x, ...) {
  v <- x$verdict
  cat(sprintf("somalign codebook alignment check  [verdict: %s]\n", v))
  cat(sprintf("  Features checked       : %d\n", nrow(x$per_feature)))
  cat(sprintf("  Critical (0%% overlap)  : %d\n", x$n_critical_features))
  cat(sprintf("  Warning  (partial)     : %d\n", x$n_warning_features))
  cat("\nCost matrix (", nrow(x$per_feature), "-feature space):\n", sep = "")
  cat(sprintf("  Median pairwise dist\u00b2  : %.4f\n",
              x$cost_summary[["median_cost"]]))
  cat(sprintf("  95th-pctile dist\u00b2      : %.4f\n",
              x$cost_summary[["p95_cost"]]))
  cat(sprintf("  Cost normalisation \u00d7   : %.4f\n",
              x$cost_summary[["cost_scale"]]))
  cat(sprintf("  Pairs within 3\u03b5        : %.1f%%\n",
              100 * x$cost_summary[["fraction_near_eps"]]))
  flagged <- x$per_feature[x$per_feature$flag != "ok", , drop = FALSE]
  if (nrow(flagged) > 0L) {
    cat("\nFlagged features:\n")
    show <- flagged[, c("feature", "overlap_fraction", "centroid_drift_sd", "flag"),
                    drop = FALSE]
    show$overlap_fraction  <- round(show$overlap_fraction,  3L)
    show$centroid_drift_sd <- round(show$centroid_drift_sd, 2L)
    print(show, row.names = FALSE)
  }
  invisible(x)
}

# Internal input-validation helper for somalign_check_codebook_alignment().
# Not exported; do not add validation logic inside compute helpers.
.somalign_validate_check_args <- function(query_codebook, reference,
                                          epsilon, stop_if_critical) {
  if (!inherits(reference, "somalign_reference")) {
    stop("`reference` must be a somalign_reference object.", call. = FALSE)
  }
  if (is.null(colnames(query_codebook))) {
    stop("`query_codebook` must have column names (one per feature).",
         call. = FALSE)
  }
  if (!is.numeric(query_codebook) && !is.integer(query_codebook)) {
    stop("`query_codebook` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0) {
    stop("`epsilon` must be a single finite positive number.", call. = FALSE)
  }
  if (!is.logical(stop_if_critical) || length(stop_if_critical) != 1L ||
      is.na(stop_if_critical)) {
    stop("`stop_if_critical` must be a scalar logical (TRUE or FALSE).",
         call. = FALSE)
  }
  invisible(NULL)
}
