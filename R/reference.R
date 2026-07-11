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
