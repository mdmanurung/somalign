#' Extend an existing reference with new SOM nodes
#'
#' Appends new nodes (e.g.\ query-derived "novel" prototypes) to a trained
#' reference **without retraining**.  The resulting object is a valid
#' `somalign_reference` and can be passed directly to [somalign_fit()].
#'
#' @section Codebook space:
#' `new_codebook` must be in **reference-scaled space** — the same space as
#' `reference$codebook`.  If your new node coordinates come from a raw-data SOM,
#' scale them first with
#' `scale(new_codebook, center = reference$center, scale = reference$scale)`.
#'
#' @section Mass strategy and OT marginal:
#' Node masses form the reference marginal `b` for the unbalanced Sinkhorn
#' solver.  After combining old and new masses, the full vector is
#' renormalised to sum to 1 by the constructor (`.somalign_normalize_masses()`),
#' so you may supply raw counts or unnormalised weights.
#'
#' When `new_node_masses` is `NULL`, each new node is assigned a mass equal to
#' the **mean of the existing node masses** (i.e.\ `mean(reference$node_masses)`
#' before renormalisation).  This gives each new node a "typical-node" share:
#' large enough for the unbalanced OT solver to route meaningful mass into those
#' columns, yet not so large as to redistribute the existing balance.  Rationale:
#' the solver's `v` update is `(b/ktu)^tau_b`; a very small `b_j` suppresses
#' transport into column `j` even when the corresponding cost is low, so a
#' "tiny mass" default would cause new nodes to be ignored at query time.
#'
#' After renormalisation the original nodes' relative proportions are preserved;
#' each original node's absolute mass shrinks by the factor
#' `n_orig / (n_orig + n_new)` approximately.
#'
#' @section Label widening:
#' When `new_labels` introduces class names not present in `reference$label_prob`,
#' the existing label probability matrix is widened with zero columns for the new
#' classes.  The new node rows receive the probabilities supplied via
#' `new_labels`; if `new_labels` is a character vector, each new node receives a
#' one-hot row for its class.  If `new_labels` is a matrix it is used directly
#' as the probability block for the new nodes (must have the same column names
#' or be column-compatible via name matching).
#'
#' @section Distance threshold fallback:
#' When `new_distance_quantiles` is `NULL`, each new node inherits
#' `reference$global_distance_quantiles` as its per-node threshold row.  This is
#' the same conservative fallback that `.somalign_thresholds()` applies to NA
#' entries at fit time, so the extended reference behaves identically for
#' "unknown" new nodes.  If `reference$global_distance_quantiles` is itself all
#' NA, new rows are filled with `Inf` (no outside-detection for those nodes).
#'
#' @param reference A `somalign_reference` object to extend.
#' @param new_codebook A numeric matrix with one row per new node and one column
#'   per feature.  Column names must match (a superset of) `reference$features`
#'   and will be reordered accordingly.
#' @param new_labels Optional label specification for the new nodes.  Either:
#'   \describe{
#'     \item{character vector}{Length `nrow(new_codebook)`.  Each element is the
#'       class name for that node; the row is one-hot encoded.}
#'     \item{matrix}{`nrow(new_codebook)` rows and named columns giving soft
#'       label probabilities.  Rows need not sum to 1; they will be
#'       row-normalised.}
#'     \item{`NULL`}{New nodes receive a uniform probability row over all
#'       existing (and any new) classes, or an empty label matrix when the
#'       reference carried no labels.}
#'   }
#' @param new_node_masses Optional non-negative numeric vector of length
#'   `nrow(new_codebook)`.  `NULL` (default) assigns each new node a mass equal
#'   to the mean existing node mass (see section *Mass strategy* above).
#' @param new_distance_quantiles Optional numeric matrix of shape
#'   `nrow(new_codebook)` x `ncol(reference$distance_quantiles)` with the same
#'   column names.  `NULL` (default) uses the existing
#'   `global_distance_quantiles` row for every new node.
#' @param ... Reserved for future arguments.
#'
#' @return A `somalign_reference` object with `n_original + n_new` rows in
#'   `codebook`, `node_masses`, `label_prob`, `distance_quantiles`, and
#'   (if present) `node_var`.
#'
#' @examples
#' ref <- somalign_reference_from_nodes(
#'   codebook = matrix(c(-1, 0, 1, 0), nrow = 2, ncol = 2,
#'                     dimnames = list(NULL, c("F1", "F2"))),
#'   features = c("F1", "F2"),
#'   center   = c(F1 = 0, F2 = 0),
#'   scale    = c(F1 = 1, F2 = 1),
#'   label_prob = matrix(c(1, 0, 0, 1), nrow = 2,
#'                       dimnames = list(NULL, c("A", "B")))
#' )
#' new_cb <- matrix(c(2, 0), nrow = 1, dimnames = list(NULL, c("F1", "F2")))
#' extended <- somalign_extend_reference(ref, new_cb, new_labels = "C")
#' nrow(extended$codebook)   # 3
#' @export
somalign_extend_reference <- function(reference,
                                      new_codebook,
                                      new_labels = NULL,
                                      new_node_masses = NULL,
                                      new_distance_quantiles = NULL,
                                      ...) {
  # ---- input validation -------------------------------------------------------
  if (!inherits(reference, "somalign_reference")) {
    stop("`reference` must be a somalign_reference object.", call. = FALSE)
  }
  new_codebook <- as.matrix(new_codebook)
  storage.mode(new_codebook) <- "double"
  if (is.null(colnames(new_codebook))) {
    if (ncol(new_codebook) == length(reference$features)) {
      colnames(new_codebook) <- reference$features
    } else {
      stop(
        "`new_codebook` has no column names and its column count (",
        ncol(new_codebook), ") does not match reference features (",
        length(reference$features), ").",
        call. = FALSE
      )
    }
  }
  missing_feats <- setdiff(reference$features, colnames(new_codebook))
  if (length(missing_feats) > 0) {
    stop(
      "`new_codebook` is missing features present in the reference: ",
      paste(missing_feats, collapse = ", "), ".",
      call. = FALSE
    )
  }
  new_codebook <- new_codebook[, reference$features, drop = FALSE]
  if (any(!is.finite(new_codebook))) {
    stop("`new_codebook` must contain finite values.", call. = FALSE)
  }

  n_orig <- nrow(reference$codebook)
  n_new  <- nrow(new_codebook)
  if (n_new < 1L) {
    stop("`new_codebook` must have at least one row.", call. = FALSE)
  }

  # ---- extended codebook ------------------------------------------------------
  extended_codebook <- rbind(reference$codebook, new_codebook)

  # ---- extended node masses ---------------------------------------------------
  # Default: mean existing mass per new node (gives each new node a
  # "typical-node" share so OT can route meaningful mass into those columns).
  if (is.null(new_node_masses)) {
    new_node_masses <- rep(mean(reference$node_masses), n_new)
  } else {
    new_node_masses <- as.numeric(new_node_masses)
    if (length(new_node_masses) != n_new ||
        any(!is.finite(new_node_masses)) ||
        any(new_node_masses < 0)) {
      stop(
        "`new_node_masses` must be a non-negative finite numeric vector of length ",
        n_new, " (one per new node).",
        call. = FALSE
      )
    }
  }
  # Combine; constructor will renormalise to sum 1
  extended_masses <- c(reference$node_masses, new_node_masses)

  # ---- extended label_prob ----------------------------------------------------
  extended_label_prob <- .somalign_extend_label_prob(
    reference$label_prob,
    new_labels,
    n_orig,
    n_new
  )

  # ---- extended distance_quantiles --------------------------------------------
  # Note: the stored distance_quantiles may contain NA (when the reference was
  # built without per-node thresholds; NA means "fall back to global" at fit
  # time).  somalign_reference_from_nodes() rejects NA in a supplied matrix, so
  # we replace NA with Inf before combining.  Inf is the canonical "no finite
  # threshold" sentinel: distance > Inf is always FALSE, so those nodes are
  # never flagged as outside-reference.  When global_distance_quantiles is
  # itself all NA (same "no-threshold" scenario) this is semantically identical.
  old_dq_raw <- reference$distance_quantiles  # [n_orig x q]
  old_dq <- old_dq_raw
  if (any(is.na(old_dq))) {
    old_dq[is.na(old_dq)] <- Inf
  }
  if (is.null(new_distance_quantiles)) {
    # Fallback: inherit global_distance_quantiles; replace NA with Inf
    global_q <- reference$global_distance_quantiles
    if (is.null(global_q) || all(!is.finite(global_q))) {
      fallback_row <- rep(Inf, ncol(old_dq))
    } else {
      fallback_row <- ifelse(is.finite(global_q), global_q, Inf)
    }
    new_dq <- matrix(
      rep(fallback_row, n_new),
      nrow  = n_new,
      ncol  = ncol(old_dq),
      byrow = TRUE,
      dimnames = list(NULL, colnames(old_dq))
    )
  } else {
    new_dq <- as.matrix(new_distance_quantiles)
    storage.mode(new_dq) <- "double"
    if (nrow(new_dq) != n_new) {
      stop(
        "`new_distance_quantiles` must have ", n_new,
        " rows (one per new node).",
        call. = FALSE
      )
    }
    if (ncol(new_dq) != ncol(old_dq)) {
      stop(
        "`new_distance_quantiles` must have ", ncol(old_dq),
        " columns to match reference$distance_quantiles.",
        call. = FALSE
      )
    }
    if (is.null(colnames(new_dq))) {
      colnames(new_dq) <- colnames(old_dq)
    }
    if (any(!is.finite(new_dq) & !is.infinite(new_dq))) {
      stop("`new_distance_quantiles` must be finite or Inf.", call. = FALSE)
    }
  }
  extended_dq <- rbind(old_dq, new_dq)

  # ---- extended node_var ------------------------------------------------------
  extended_node_var <- .somalign_extend_node_var(
    reference$node_var,
    n_new,
    reference$features
  )

  # ---- re-wrap through the clean constructor ----------------------------------
  # Pass global_distance_quantiles explicitly so it is not recomputed from
  # the (now potentially higher) column maxima of the extended matrix.
  out <- somalign_reference_from_nodes(
    codebook                  = extended_codebook,
    features                  = reference$features,
    center                    = reference$center,
    scale                     = reference$scale,
    node_masses               = extended_masses,
    label_prob                = extended_label_prob,
    distance_quantiles        = extended_dq,
    global_distance_quantiles = reference$global_distance_quantiles,
    node_var                  = extended_node_var
  )

  out
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Build the extended [n_total x n_classes_total] label_prob matrix.
# Handles: NULL new_labels, character vector, or matrix of soft probabilities.
# Always widens existing columns when new classes are introduced.
.somalign_extend_label_prob <- function(old_lp, new_labels, n_orig, n_new) {
  n_total <- n_orig + n_new

  # No labels in reference and none supplied -> stay label-free
  if ((is.null(old_lp) || ncol(old_lp) == 0L) && is.null(new_labels)) {
    return(matrix(numeric(0), nrow = n_total, ncol = 0L))
  }

  old_classes <- if (!is.null(old_lp)) colnames(old_lp) else character(0L)

  # Derive new node label_prob block
  if (is.null(new_labels)) {
    # New nodes get uniform over existing classes (or empty)
    if (length(old_classes) == 0L) {
      new_block <- matrix(numeric(0), nrow = n_new, ncol = 0L)
    } else {
      new_block <- matrix(
        1.0 / length(old_classes),
        nrow  = n_new,
        ncol  = length(old_classes),
        dimnames = list(NULL, old_classes)
      )
    }
    new_classes <- old_classes
  } else if (is.character(new_labels)) {
    if (length(new_labels) != n_new) {
      stop(
        "Character `new_labels` must have length ", n_new,
        " (one label per new node).",
        call. = FALSE
      )
    }
    all_classes <- union(old_classes, unique(new_labels))
    new_block <- matrix(
      0.0,
      nrow  = n_new,
      ncol  = length(all_classes),
      dimnames = list(NULL, all_classes)
    )
    for (i in seq_len(n_new)) {
      new_block[i, new_labels[i]] <- 1.0
    }
    new_classes <- all_classes
  } else {
    # Treat as soft probability matrix
    new_block <- as.matrix(new_labels)
    storage.mode(new_block) <- "double"
    if (nrow(new_block) != n_new) {
      stop(
        "Matrix `new_labels` must have ", n_new,
        " rows (one per new node).",
        call. = FALSE
      )
    }
    if (is.null(colnames(new_block))) {
      if (length(old_classes) > 0L && ncol(new_block) == length(old_classes)) {
        colnames(new_block) <- old_classes
      } else {
        colnames(new_block) <- paste0("label_", seq_len(ncol(new_block)))
      }
    }
    if (any(!is.finite(new_block)) || any(new_block < 0)) {
      stop(
        "Matrix `new_labels` must contain non-negative finite values.",
        call. = FALSE
      )
    }
    new_classes <- union(old_classes, colnames(new_block))
  }

  all_classes <- union(old_classes, new_classes)
  n_classes   <- length(all_classes)

  if (n_classes == 0L) {
    return(matrix(numeric(0), nrow = n_total, ncol = 0L))
  }

  out <- matrix(
    0.0,
    nrow     = n_total,
    ncol     = n_classes,
    dimnames = list(NULL, all_classes)
  )

  # Fill old rows
  if (!is.null(old_lp) && ncol(old_lp) > 0L) {
    out[seq_len(n_orig), colnames(old_lp)] <- old_lp
  }

  # Fill new rows
  new_row_idx <- seq(n_orig + 1L, n_total)
  if (ncol(new_block) > 0L) {
    out[new_row_idx, colnames(new_block)] <- new_block
  }

  out
}

# Extend node_var: if reference has node_var, default new rows to column-wise
# means of existing rows (same conservative fallback as nodes with <2 cells).
.somalign_extend_node_var <- function(old_nv, n_new, features) {
  if (is.null(old_nv)) return(NULL)
  # Column-wise means of existing node variances as fallback for new nodes
  col_means <- colMeans(old_nv)
  new_rows <- matrix(
    rep(col_means, n_new),
    nrow  = n_new,
    ncol  = ncol(old_nv),
    byrow = TRUE,
    dimnames = list(NULL, colnames(old_nv))
  )
  rbind(old_nv, new_rows)
}
