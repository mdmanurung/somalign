#' Subset a reference to a shared marker panel
#'
#' Returns a new \code{somalign_reference} whose codebook, scaling vectors, and
#' feature list are restricted to the markers named in \code{markers}.  All
#' per-node arrays (\code{node_masses}, \code{label_prob},
#' \code{distance_quantiles}, \code{global_distance_quantiles}) are carried
#' through unchanged.
#'
#' Use this helper when a query dataset was measured on a \emph{subset} of the
#' markers that the reference SOM was trained on.  Pass the returned reference
#' to \code{\link{somalign_query}()} together with query data that only contains
#' \code{markers}; the OT cost matrix and node-shift correction will then be
#' computed on the shared marker subspace.
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
#'
#' @return A \code{somalign_reference} object with \code{$features} equal to
#'   \code{markers} (in the order given), and \code{$codebook},
#'   \code{$center}, \code{$scale}, and (if present) \code{$node_var} all
#'   restricted to those columns.
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
#' # subset to the two T-cell markers; CD19 dropped
#' ref_sub <- somalign_reference_subset_markers(ref_full, c("CD3", "CD4", "CD8"))
#' @export
somalign_reference_subset_markers <- function(reference, markers) {
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
  node_var_sub <- if (!is.null(reference$node_var))
    reference$node_var[, markers, drop = FALSE]
  else
    NULL

  # --- rebuild via the canonical constructor ---------------------------------
  somalign_reference_from_nodes(
    codebook                  = codebook_sub,
    features                  = markers,
    center                    = center_sub,
    scale                     = scale_sub,
    node_masses               = reference$node_masses,
    label_prob                = reference$label_prob,
    distance_quantiles        = reference$distance_quantiles,
    global_distance_quantiles = reference$global_distance_quantiles,
    node_var                  = node_var_sub
  )
}
