#' Train a reference SOM and build a somalign reference
#'
#' @param data Numeric matrix or data frame containing old/reference samples.
#' @param labels Optional labels, one per row of `data`.
#' @param features Optional feature names to use. Defaults to all columns.
#' @param grid Optional `kohonen::somgrid()` object.
#' @param rlen Number of SOM training iterations passed to [kohonen::som()].
#' @param alpha Learning-rate schedule passed to [kohonen::som()].
#' @param ... Additional arguments passed to [kohonen::som()].
#'
#' @return A `somalign_reference` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' labels <- rep(c("A", "B"), each = 5)
#' ref <- somalign_train_reference(mat, labels = labels,
#'                                 grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' @export
somalign_train_reference <- function(data,
                                     labels = NULL,
                                     features = NULL,
                                     grid = NULL,
                                     rlen = 100,
                                     alpha = c(0.05, 0.01),
                                     ...) {
  .somalign_check_som_train_args(data, labels, features, grid, rlen, alpha)
  data <- .somalign_prepare_feature_matrix(data, features, what = "data")
  scaling <- .somalign_compute_scaling(data)
  scaled <- .somalign_scale_matrix(data, scaling$center, scaling$scale)
  if (is.null(grid)) {
    grid <- .somalign_default_grid(nrow(data))
  }
  som_ref <- kohonen::som(scaled, grid = grid, rlen = rlen, alpha = alpha, ...)
  somalign_reference(
    som_ref = som_ref,
    data = data,
    labels = labels,
    features = colnames(data),
    center = scaling$center,
    scale = scaling$scale,
    codebook_space = "reference_scaled"
  )
}

#' Build a reference object from an existing SOM and old data
#'
#' @param som_ref A `kohonen` SOM object, or a SOM-like object containing a
#'   codebook matrix.
#' @param data Numeric old/reference data used to compute scaling, masses,
#'   labels, and distance thresholds.
#' @param labels Optional labels, one per row of `data`.
#' @param features Optional feature names to use. Defaults to all columns.
#' @param center Optional saved feature centers. Computed from `data` when
#'   omitted.
#' @param scale Optional saved feature scales. Computed from `data` when
#'   omitted.
#' @param codebook_space Coordinate system of the existing `som_ref` codebook.
#'   Use `"reference_scaled"` when the SOM was trained on `data` transformed
#'   with `center` and `scale`; use `"raw"` when the SOM was trained on raw
#'   feature values and should be transformed into reference-scaled space.
#' @param quantile_probs Distance quantiles used for outside-reference flags.
#'
#' @return A `somalign_reference` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' g <- kohonen::somgrid(2, 2, "hexagonal")
#' som_obj <- kohonen::som(scale(mat), grid = g, rlen = 5)
#' ref <- somalign_reference(som_obj, mat, labels = rep(c("A", "B"), each = 5),
#'                           codebook_space = "reference_scaled")
#' @export
somalign_reference <- function(som_ref,
                               data,
                               labels = NULL,
                               features = NULL,
                               center = NULL,
                               scale = NULL,
                               codebook_space = NULL,
                               quantile_probs = c(0.5, 0.9, 0.95, 0.99)) {
  .somalign_check_reference_args(data, labels, features, quantile_probs)
  data <- .somalign_prepare_feature_matrix(data, features, what = "data")
  features <- colnames(data)

  resolved <- .somalign_resolve_center_scale(center, scale, data)
  center <- resolved$center
  scale <- resolved$scale
  center <- .somalign_named_numeric(center, features, "center")
  scale <- .somalign_named_numeric(scale, features, "scale")
  .somalign_validate_scale(scale)

  codebook <- .somalign_get_codebook(som_ref, features = features, what = "som_ref")
  codebook_space <- .somalign_validate_codebook_space(codebook_space)
  if (identical(codebook_space, "raw")) {
    codebook <- .somalign_scale_matrix(codebook, center, scale)
  }
  scaled_data <- .somalign_scale_matrix(data, center, scale)
  projected <- .somalign_nearest_code(scaled_data, codebook)
  n_nodes <- nrow(codebook)
  node_masses <- .somalign_node_masses(projected$unit, n_nodes)
  label_prob <- .somalign_label_probabilities(labels, projected$unit, n_nodes)
  quantiles <- .somalign_distance_quantiles(projected$distance, projected$unit, n_nodes, quantile_probs)

  structure(
    list(
      som_ref = som_ref,
      features = features,
      center = center,
      scale = scale,
      codebook = codebook,
      node_masses = node_masses,
      label_prob = label_prob,
      distance_quantiles = quantiles$node,
      global_distance_quantiles = quantiles$global,
      reference_units = projected$unit,
      reference_distances = projected$distance,
      n_samples = nrow(data)
    ),
    class = "somalign_reference"
  )
}

.somalign_resolve_center_scale <- function(center, scale, data) {
  if (is.null(center) || is.null(scale)) {
    scaling <- .somalign_compute_scaling(data)
    if (is.null(center)) {
      center <- scaling$center
    }
    if (is.null(scale)) {
      scale <- scaling$scale
    }
  }
  list(center = center, scale = scale)
}

#' Build a reference object from saved node-level artifacts
#'
#' @param codebook Reference node codebook matrix.
#' @param features Feature names and order.
#' @param center Saved reference feature centers.
#' @param scale Saved reference feature scales.
#' @param node_masses Optional reference node masses.
#' @param label_prob Optional node-by-label probability matrix.
#' @param distance_quantiles Optional node-by-quantile distance matrix.
#' @param global_distance_quantiles Optional global reference distance
#'   quantiles.
#'
#' @return A `somalign_reference` object.
#' @examples
#' cb <- matrix(c(0.1, 0.2, -0.1, 0.3, 0.4, -0.2, 0.0, 0.1),
#'              nrow = 4, ncol = 2,
#'              dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_reference_from_nodes(
#'   codebook = cb,
#'   features = c("F1", "F2"),
#'   center   = c(F1 = 0, F2 = 0),
#'   scale    = c(F1 = 1, F2 = 1)
#' )
#' @export
somalign_reference_from_nodes <- function(codebook,
                                          features,
                                          center,
                                          scale,
                                          node_masses = NULL,
                                          label_prob = NULL,
                                          distance_quantiles = NULL,
                                          global_distance_quantiles = NULL) {
  .somalign_check_data_arg(codebook, what = "codebook")
  .somalign_check_opt_char(features, what = "features")
  codebook <- .somalign_validate_node_codebook(codebook, features)
  center <- .somalign_named_numeric(center, features, "center")
  scale <- .somalign_named_numeric(scale, features, "scale")
  .somalign_validate_scale(scale)

  n_nodes <- nrow(codebook)
  node_masses <- .somalign_normalize_masses(node_masses, n_nodes, "node_masses")
  label_prob <- .somalign_normalize_label_prob(label_prob, n_nodes)
  distance_quantiles <- .somalign_prepare_distance_quantiles(distance_quantiles, n_nodes)

  .somalign_warn_from_nodes(label_prob, distance_quantiles)
  global_distance_quantiles <- .somalign_resolve_global_quantiles(
    distance_quantiles,
    global_distance_quantiles
  )

  structure(
    list(
      som_ref = NULL,
      features = features,
      center = center,
      scale = scale,
      codebook = codebook,
      node_masses = node_masses,
      label_prob = label_prob,
      distance_quantiles = distance_quantiles,
      global_distance_quantiles = global_distance_quantiles,
      reference_units = NULL,
      reference_distances = NULL,
      n_samples = NA_integer_
    ),
    class = "somalign_reference"
  )
}

.somalign_validate_node_codebook <- function(codebook, features) {
  codebook <- .somalign_as_matrix(codebook, what = "codebook")
  if (!is.character(features) || length(features) == 0) {
    stop("`features` must be a non-empty character vector.", call. = FALSE)
  }
  if (length(unique(features)) != length(features)) {
    stop("Duplicated requested features: ", paste(features[duplicated(features)], collapse = ", "), call. = FALSE)
  }
  if (is.null(colnames(codebook))) {
    colnames(codebook) <- features
  }
  codebook <- .somalign_select_features(codebook, features, what = "codebook")
  .somalign_validate_finite(codebook, what = "codebook")
  codebook
}

.somalign_warn_from_nodes <- function(label_prob, distance_quantiles) {
  if (ncol(label_prob) == 0) {
    message(
      "somalign_reference_from_nodes: no label probabilities supplied; ",
      "label transfer will be disabled for this reference."
    )
  }
  if (all(!is.finite(distance_quantiles))) {
    message(
      "somalign_reference_from_nodes: distance quantiles not supplied; ",
      "outside-reference detection will be disabled for this reference."
    )
  }
}

.somalign_resolve_global_quantiles <- function(distance_quantiles,
                                              global_distance_quantiles) {
  if (is.null(global_distance_quantiles)) {
    global_distance_quantiles <- apply(distance_quantiles, 2, max)
  }
  global_distance_quantiles <- as.numeric(global_distance_quantiles)
  if (is.null(names(global_distance_quantiles))) {
    names(global_distance_quantiles) <- colnames(distance_quantiles)
  }
  global_distance_quantiles
}

#' Build a reference object directly from a trained kohonen SOM
#'
#' Constructs a \code{somalign_reference} without reprojecting any cells by
#' reusing information already stored inside the trained kohonen object:
#'
#' \describe{
#'   \item{X codebook}{\code{codes[[1]]} — the reference node positions in
#'     feature space, already used by the existing API.}
#'   \item{Node masses}{\code{tabulate(som$unit.classif)} — exact counts over
#'     *all* training cells, zero cost.}
#'   \item{Label probabilities}{\code{codes[[2]]} — the supervised Y-layer
#'     codebook from an \code{xyf}/\code{supersom} object.  Each row is a
#'     per-node distribution over cell-type labels.  Absent for plain
#'     \code{som()} objects, in which case label transfer is disabled.}
#'   \item{Distance quantiles}{Recomputed in reference-scaled X-space from the
#'     embedded \code{som$data[[1]]} using each cell's known
#'     \code{unit.classif} assignment.  This is O(N \eqn{\times} p) with no
#'     argmax and no O(N \eqn{\times} nodes) memory peak.}
#' }
#'
#' \strong{Note on partition semantics.}  On an \code{xyf}/\code{supersom}
#' trained with equal layer weights, \code{unit.classif} reflects a joint X+Y
#' (label-weighted) assignment, not a pure-X nearest-node assignment.  Node
#' masses therefore match the SOM's own supervised partition rather than
#' somalign's X-only geometry.  Distance quantiles are still computed in X-only
#' space so the outside-reference threshold is on the same scale as
#' \code{somalign_fit()}'s query distances.
#'
#' @param som A trained \code{kohonen} SOM object (from \code{kohonen::som()},
#'   \code{kohonen::supersom()}, or \code{kohonen::xyf()}) with
#'   \code{keep.data = TRUE} (the kohonen default).
#' @param center Named numeric vector of reference feature centers, one per
#'   feature.  Required.
#' @param scale Named numeric vector of reference feature scales, one per
#'   feature (must be strictly positive).  Required.
#' @param codebook_space Coordinate system of the SOM codebook.
#'   \code{"reference_scaled"} when the SOM was trained on data already
#'   transformed by \code{center} and \code{scale};
#'   \code{"raw"} when the SOM was trained on raw feature values that should
#'   be transformed into reference-scaled space before use.
#' @param labels \code{"codebook"} (default) reads the per-node label
#'   distribution from the Y-layer \code{codes[[2]]} and enables label
#'   transfer.  \code{"none"} skips label extraction and disables label
#'   transfer regardless of whether a Y-layer is present.
#' @param quantile_probs Quantile levels for per-node distance thresholds.
#'   Passed to \code{\link{somalign_reference_from_nodes}}.
#' @param distance_chunk_size Number of cells to process per chunk when
#'   computing X-space cell-to-node distances.  Reduce if memory is tight;
#'   increase for faster throughput.  Default 1e6.
#'
#' @return A \code{somalign_reference} object, identical in structure to the
#'   output of \code{\link{somalign_reference}} but built without reprojecting
#'   any cells.
#'
#' @seealso [somalign_reference()], [somalign_reference_from_nodes()],
#'   [somalign_query()], [somalign_fit()]
#'
#' @examples
#' set.seed(1)
#' n <- 120
#' X <- matrix(rnorm(n * 2), nrow = n, ncol = 2,
#'             dimnames = list(NULL, c("F1", "F2")))
#' Y <- cbind(A = rep(c(1, 0), each = n / 2),
#'            B = rep(c(0, 1), each = n / 2))
#' som_obj <- kohonen::supersom(list(X, Y),
#'                              grid = kohonen::somgrid(3, 3, "hexagonal"),
#'                              rlen = 5, keep.data = TRUE)
#' center <- colMeans(X)
#' scale  <- apply(X, 2, sd)
#' X_scaled <- scale(X, center = center, scale = scale)
#' som_scaled <- kohonen::supersom(list(X_scaled, Y),
#'                                 grid = kohonen::somgrid(3, 3, "hexagonal"),
#'                                 rlen = 5, keep.data = TRUE)
#' ref <- somalign_reference_from_som(som_scaled,
#'                                    center = center, scale = scale,
#'                                    codebook_space = "reference_scaled")
#' @export
somalign_reference_from_som <- function(som,
                                        center,
                                        scale,
                                        codebook_space,
                                        labels = c("codebook", "none"),
                                        quantile_probs = c(0.5, 0.9, 0.95, 0.99),
                                        distance_chunk_size = 1e6L) {
  .somalign_check_reference_from_som_args(quantile_probs, distance_chunk_size)
  labels <- match.arg(labels)

  # --- 1. X codebook and features -------------------------------------------
  codebook <- .somalign_get_codebook(som, what = "som")
  features <- colnames(codebook)
  n_nodes <- nrow(codebook)

  # --- 2. center / scale / codebook_space -----------------------------------
  codebook_space <- .somalign_validate_codebook_space(codebook_space)
  center <- .somalign_named_numeric(center, features, "center")
  scale  <- .somalign_named_numeric(scale,  features, "scale")
  .somalign_validate_scale(scale)

  if (identical(codebook_space, "raw")) {
    codebook <- .somalign_scale_matrix(codebook, center, scale)
  }

  # --- 3. Node masses from unit.classif (all training cells) ----------------
  unit <- som$unit.classif
  if (is.null(unit) || length(unit) == 0L) {
    stop(
      "`som$unit.classif` is absent or empty.  ",
      "Ensure the SOM was trained with `keep.data = TRUE`.",
      call. = FALSE
    )
  }
  if (max(unit) > n_nodes) {
    stop(
      "`som$unit.classif` contains node indices larger than the number of ",
      "codebook rows (", n_nodes, ").  The SOM may have been modified after ",
      "training.",
      call. = FALSE
    )
  }
  node_masses <- tabulate(unit, nbins = n_nodes)

  # --- 4. Label probabilities from Y-layer codebook -------------------------
  label_prob <- if (identical(labels, "codebook")) {
    .somalign_extract_label_codes(som)
  } else {
    NULL
  }

  # --- 5. X-space distances (chunked, O(N * p), no argmax) -----------------
  X <- .somalign_extract_som_data(som)
  if (!identical(as.integer(ncol(X)), as.integer(length(features)))) {
    stop(
      "`som$data[[1]]` has ", ncol(X), " columns but the X codebook has ",
      length(features), " features.  The SOM object may be inconsistent.",
      call. = FALSE
    )
  }
  if (!identical(as.integer(nrow(X)), as.integer(length(unit)))) {
    stop(
      "`som$data[[1]]` has ", nrow(X), " rows but `som$unit.classif` has ",
      length(unit), " entries.  The SOM object may be inconsistent.",
      call. = FALSE
    )
  }

  # If raw codebook, data must also be scaled before distance computation
  if (identical(codebook_space, "raw")) {
    X <- .somalign_scale_matrix(X, center, scale)
  }

  # Select features in the same order as the codebook (skip copy if already aligned)
  if (!identical(colnames(X), features)) {
    X <- X[, features, drop = FALSE]
  }

  d <- .somalign_som_cell_distances(X, codebook, unit,
                                    chunk_size = distance_chunk_size)

  quantiles <- .somalign_distance_quantiles(d, unit, n_nodes, quantile_probs)

  # --- 6. Delegate validation/normalisation to somalign_reference_from_nodes ---
  ref <- somalign_reference_from_nodes(
    codebook               = codebook,
    features               = features,
    center                 = center,
    scale                  = scale,
    node_masses            = node_masses,
    label_prob             = label_prob,
    distance_quantiles     = quantiles$node,
    global_distance_quantiles = quantiles$global
  )

  # --- 7. Enrich with per-cell fields (mirrors somalign_reference output) ---
  # Strip $data from the stored SOM copy: everything we needed is already
  # extracted above, and keeping the full N×p matrix (up to ~10 GB for 44.6M
  # cells) would double the memory footprint of the returned reference object.
  som_ref <- som
  som_ref$data <- NULL
  ref$som_ref             <- som_ref
  ref$reference_units     <- as.integer(unit)
  ref$reference_distances <- d
  ref$n_samples           <- length(unit)

  ref
}
