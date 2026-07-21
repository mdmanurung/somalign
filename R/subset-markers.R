#' Subset a reference to a shared marker panel
#'
#' Returns a new \code{somalign_reference} whose codebook, scaling vectors, and
#' feature list are restricted to the markers named in \code{markers}.  When
#' \code{reference_data} is supplied, the outside-reference detection statistics
#' (\code{distance_quantiles}, \code{global_distance_quantiles}, and
#' \code{node_var}) are \strong{recomputed} in the shared-marker subspace so
#' that distance thresholds remain calibrated after dimensionality reduction.
#' When \code{reference_data} is \code{NULL} (the default), those statistics
#' are set to safe sentinel values that disable outside-reference flagging, and
#' a warning is emitted.
#'
#' Use this helper when a query dataset was measured on a \emph{subset} of the
#' markers that the reference SOM was trained on.  Pass the returned reference
#' to \code{\link{somalign_query}()} together with query data that only contains
#' \code{markers}; the OT cost matrix and node-shift correction will then be
#' computed on the shared marker subspace.
#'
#' @section Outside-reference calibration:
#' \strong{Silent-failure warning.}  The full-\emph{p} \code{distance_quantiles}
#' stored in the original reference are computed in a higher-dimensional space.
#' Euclidean distances in the \emph{k}-marker subspace (\eqn{k < p}) are
#' uniformly smaller than in the full \emph{p}-marker space, so carrying those
#' thresholds forward uncorrected would mean that the outside-reference distance
#' flag and surprisal score are silently never triggered.  This function resolves
#' the ambiguity by either recomputing the statistics (when \code{reference_data}
#' is supplied) or disabling detection with an explicit warning (when it is not).
#'
#' Supplying \code{reference_data} enables calibrated outside-reference detection
#' in the subspace.  The raw reference cells in \code{reference_data} are scaled
#' with the \emph{subset} center and scale, projected to the subset codebook via
#' nearest-code assignment, and per-node distance quantiles and variances are
#' recomputed exactly as in the original reference constructor.
#'
#' @section Correctness caveat:
#' Removing markers from the reference codebook changes the geometry of the OT
#' cost matrix: the squared-Euclidean distance between two SOM nodes is
#' computed only over the retained dimensions.  If the dropped markers are
#' informative for separating cell populations (e.g.\ lineage-defining markers),
#' the inter-node distances will be compressed and cell-type nodes that were
#' far apart in full marker-space may appear close in the shared subspace.
#' Label transfer accuracy will then be limited by how well the shared markers
#' distinguish the reference populations.  \strong{Do not} silently accept the
#' subset reference when you know the dropped markers are lineage-defining;
#' inspect label confidence and entropy in the returned fit to assess how much
#' information was lost.
#'
#' Fully disjoint panels (zero shared markers) cannot be handled by subspace
#' projection and require Gromov-Wasserstein optimal transport, which is out of
#' scope for this function.
#'
#' @param reference A \code{somalign_reference} object, as returned by
#'   \code{\link{somalign_train_reference}()},
#'   \code{\link{somalign_reference}()}, or
#'   \code{\link{somalign_reference_from_nodes}()}.
#' @param markers A non-empty character vector of marker names.  Must be a
#'   subset of \code{reference$features}; an error is raised if any name is
#'   absent from the reference or if the vectors are disjoint.
#' @param reference_data Optional numeric matrix (or data frame) of the raw
#'   reference cells used to build \code{reference}.  Must contain at least the
#'   columns named in \code{markers}.  When supplied, \code{distance_quantiles},
#'   \code{global_distance_quantiles}, and \code{node_var} are recomputed in the
#'   \code{markers} subspace so that outside-reference detection remains
#'   calibrated.  When \code{NULL} (the default), those fields are set to
#'   \code{Inf} sentinels that disable distance-based outside-reference flagging,
#'   and \code{node_var} is set to \code{NULL} to disable surprisal scoring; a
#'   \code{warning()} is emitted in that case.
#'
#' @return A \code{somalign_reference} object with \code{$features} equal to
#'   \code{markers} (in the order given), and \code{$codebook},
#'   \code{$center}, \code{$scale}, and (if present and \code{reference_data}
#'   supplied) \code{$node_var} all restricted or recomputed for those columns.
#'   \code{$distance_quantiles} and \code{$global_distance_quantiles} are either
#'   recomputed from \code{reference_data} (calibrated) or set to \code{Inf}
#'   (detection disabled).
#'
#' @seealso [somalign_reference()], [somalign_reference_from_nodes()],
#'   [somalign_query()], [somalign_fit()]
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 4,
#'               dimnames = list(NULL, c("CD3", "CD4", "CD8", "CD19")))
#' ref_full <- somalign_train_reference(
#'   mat,
#'   labels = rep(c("T", "B"), each = 10),
#'   grid   = kohonen::somgrid(2, 2, "hexagonal"),
#'   rlen   = 5
#' )
#' # subset with recomputed distance quantiles
#' ref_sub <- somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"),
#'                                              reference_data = mat)
#' # subset with detection disabled (no reference_data supplied)
#' ref_sub_nodet <- suppressWarnings(
#'   somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
#' )
#' @export
somalign_reference_subset_markers <- function(reference, markers,
                                              reference_data = NULL) {
  # --- input validation ------------------------------------------------------
  if (!inherits(reference, "somalign_reference")) {
    stop(
      "`reference` must be a `somalign_reference` object.",
      call. = FALSE
    )
  }
  if (!is.character(markers) || length(markers) == 0L) {
    stop("`markers` must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(markers)) {
    stop(
      "Duplicated markers: ",
      paste(markers[duplicated(markers)], collapse = ", "),
      call. = FALSE
    )
  }

  missing_markers <- setdiff(markers, reference$features)
  if (length(missing_markers) > 0L) {
    stop(
      "The following markers are not in the reference feature set: ",
      paste(missing_markers, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (length(markers) == 0L) {
    stop(
      "`markers` and `reference$features` are disjoint (no shared markers). ",
      "Fully disjoint panels require Gromov-Wasserstein transport, which is ",
      "not supported by this function.",
      call. = FALSE
    )
  }

  # --- subset marker-indexed arrays ------------------------------------------
  codebook_sub <- reference$codebook[, markers, drop = FALSE]
  center_sub   <- reference$center[markers]
  scale_sub    <- reference$scale[markers]

  # --- recompute or disable distance quantiles / node_var -------------------
  n_nodes <- nrow(codebook_sub)

  if (!is.null(reference_data)) {
    # Validate that reference_data is a matrix/data.frame with the right columns
    if (!is.matrix(reference_data) && !is.data.frame(reference_data)) {
      stop("`reference_data` must be a numeric matrix or data frame.", call. = FALSE)
    }
    missing_ref_cols <- setdiff(markers, colnames(reference_data))
    if (length(missing_ref_cols) > 0L) {
      stop(
        "`reference_data` is missing columns required by `markers`: ",
        paste(missing_ref_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }

    # Select and order columns to match markers
    ref_mat <- as.matrix(reference_data[, markers, drop = FALSE])
    storage.mode(ref_mat) <- "double"

    # Scale in the subset subspace
    scaled_ref <- .somalign_scale_matrix(ref_mat, center_sub, scale_sub)

    # Project to subset codebook via nearest-code assignment
    projected <- .somalign_nearest_code_chunked(scaled_ref, codebook_sub)

    # Derive quantile_probs from existing colnames (e.g. "50%", "90%", ...)
    # This guarantees the recomputed matrix has the same shape as the original.
    existing_colnames <- colnames(reference$distance_quantiles)
    if (!is.null(existing_colnames) && length(existing_colnames) > 0L) {
      quantile_probs <- as.numeric(sub("%", "", existing_colnames, fixed = TRUE)) / 100
      quantile_probs <- quantile_probs[is.finite(quantile_probs)]
      if (length(quantile_probs) == 0L) {
        quantile_probs <- c(0.5, 0.9, 0.95, 0.99)
      }
    } else {
      quantile_probs <- c(0.5, 0.9, 0.95, 0.99)
    }

    # Recompute distance quantiles and node_var in the subspace
    quantiles <- .somalign_distance_quantiles(
      projected$distance, projected$unit, n_nodes, quantile_probs
    )
    node_var_sub <- .somalign_node_var(scaled_ref, projected$unit, n_nodes)

    dq_sub     <- quantiles$node
    global_dq  <- quantiles$global

  } else {
    # No reference_data supplied: disable detection with Inf sentinels
    warning(
      "Outside-reference distance detection and surprisal scoring are DISABLED ",
      "for this subset reference because `reference_data` was not supplied. ",
      "Full-p distance_quantiles are NOT valid in the ", length(markers),
      "-marker subspace. Supply `reference_data` (the original reference cells) ",
      "to recompute calibrated thresholds in the shared-marker subspace.",
      call. = FALSE
    )

    # Preserve shape/colnames of original matrix, set all values to Inf
    dq_sub <- reference$distance_quantiles
    dq_sub[] <- Inf

    # global_distance_quantiles: also set to Inf
    global_dq <- reference$global_distance_quantiles
    global_dq[] <- Inf

    # node_var: NULL to disable surprisal scoring
    node_var_sub <- NULL
  }

  # --- rebuild via the canonical constructor ---------------------------------
  somalign_reference_from_nodes(
    codebook                  = codebook_sub,
    features                  = markers,
    center                    = center_sub,
    scale                     = scale_sub,
    node_masses               = reference$node_masses,
    label_prob                = reference$label_prob,
    distance_quantiles        = dq_sub,
    global_distance_quantiles = global_dq,
    node_var                  = node_var_sub
  )
}
