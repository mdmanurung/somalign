#' Prepare query data and attach or train a query SOM
#'
#' Query data are always scaled with the saved reference center and scale.
#'
#' @param data Numeric query data.
#' @param reference A `somalign_reference` object.
#' @param som_query Optional query SOM or SOM-like object with a codebook.
#' @param codebook_space Coordinate system of the `som_query` codebook. Only
#'   used when `som_query` is supplied. `"reference_scaled"` (default) assumes
#'   the codebook was already trained on query data scaled with
#'   `reference$center` and `reference$scale`; `"raw"` re-scales the codebook
#'   into reference-scaled space before use.
#' @param grid Optional `kohonen::somgrid()` object when `som_query` is omitted.
#' @param rlen Number of SOM training iterations passed to [kohonen::som()].
#' @param alpha Learning-rate schedule passed to [kohonen::som()].
#' @param features Optional feature names. Defaults to the reference feature
#'   order.
#' @param ... Additional arguments passed to [kohonen::som()].
#'
#' @return A `somalign_query` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' @export
somalign_query <- function(data,
                           reference,
                           som_query = NULL,
                           codebook_space = c("reference_scaled", "raw"),
                           grid = NULL,
                           rlen = 100,
                           alpha = c(0.05, 0.01),
                           features = NULL,
                           ...) {
  .somalign_check_reference(reference)
  features <- .somalign_query_features(features, reference)
  data <- .somalign_prepare_feature_matrix(data, features, what = "query data")
  scaled_data <- .somalign_scale_matrix(data, reference$center, reference$scale)

  query_som <- .somalign_resolve_query_som(
    som_query, codebook_space, scaled_data, grid, rlen, alpha, ...
  )
  som_query <- query_som$som_query
  codebook <- .somalign_query_codebook(
    som_query,
    reference,
    query_som$user_supplied_som,
    query_som$codebook_space
  )
  sample_map <- .somalign_nearest_code(scaled_data, codebook)

  .somalign_new_query(data, scaled_data, som_query, codebook, sample_map, reference)
}

.somalign_query_features <- function(features, reference) {
  if (is.null(features)) {
    return(reference$features)
  }
  if (!identical(features, reference$features)) {
    missing <- setdiff(reference$features, features)
    if (length(missing) > 0 || length(features) != length(reference$features)) {
      stop("`features` must match the reference feature set.", call. = FALSE)
    }
    features <- reference$features
  }
  features
}

.somalign_resolve_query_som <- function(som_query, codebook_space, scaled_data,
                                        grid, rlen, alpha, ...) {
  user_supplied_som <- !is.null(som_query)
  if (is.null(som_query)) {
    if (is.null(grid)) {
      grid <- .somalign_default_grid(nrow(scaled_data))
    }
    som_query <- kohonen::som(scaled_data, grid = grid, rlen = rlen, alpha = alpha, ...)
  } else {
    codebook_space <- match.arg(codebook_space, c("reference_scaled", "raw"))
  }
  list(
    som_query = som_query,
    codebook_space = codebook_space,
    user_supplied_som = user_supplied_som
  )
}

.somalign_query_codebook <- function(som_query, reference, user_supplied_som,
                                     codebook_space) {
  codebook <- .somalign_get_codebook(som_query, features = reference$features, what = "som_query")
  if (user_supplied_som && identical(codebook_space, "raw")) {
    codebook <- .somalign_scale_matrix(codebook, reference$center, reference$scale)
  }
  codebook
}

.somalign_new_query <- function(data, scaled_data, som_query, codebook,
                                sample_map, reference) {
  node_masses <- .somalign_node_masses(sample_map$unit, nrow(codebook))
  sample_id <- rownames(data)
  if (is.null(sample_id)) {
    sample_id <- as.character(seq_len(nrow(data)))
  }

  structure(
    list(
      data = data,
      scaled_data = scaled_data,
      som_query = som_query,
      codebook = codebook,
      node_masses = node_masses,
      sample_unit = sample_map$unit,
      sample_distance = sample_map$distance,
      sample_id = sample_id,
      reference_features = reference$features
    ),
    class = "somalign_query"
  )
}

#' Build a query object from a pre-trained kohonen SOM
#'
#' Creates a \code{somalign_query} by reusing the per-cell node assignments
#' (\code{som$unit.classif}) already computed during SOM training, bypassing
#' the O(N \eqn{\times} nodes) per-cell argmax that
#' \code{\link{somalign_query}()} would otherwise perform.
#'
#' @param som A trained kohonen SOM (output of \code{kohonen::som()} or
#'   \code{kohonen::supersom()}) with \code{$unit.classif} populated
#'   (\code{keep.data = TRUE} during training is \emph{not} required here, but
#'   \code{unit.classif} must be present).
#' @param data Numeric matrix of raw query cell values (cells \eqn{\times}
#'   features).  Must include all features in \code{reference$features}, and
#'   \code{nrow(data)} must equal \code{length(som$unit.classif)}.
#' @param reference A \code{somalign_reference} object.
#' @param codebook Optional numeric matrix of SOM codebook vectors in the
#'   coordinate space given by \code{codebook_space} (nodes \eqn{\times}
#'   features).  When \code{NULL} (default) the X-layer codebook
#'   \code{som$codes[[1]]} is used.  Supply an explicitly transformed codebook
#'   (e.g.\ after winsorisation and rescaling into reference-scaled space) when
#'   the SOM's native codebook has been post-processed before alignment.
#' @param codebook_space Coordinate system of the codebook.
#'   \code{"reference_scaled"} (default) assumes the codebook is already in the
#'   reference-scaled coordinate system.  \code{"raw"} applies
#'   \code{reference$center} and \code{reference$scale} before use.
#' @param features Optional character vector of feature names.  Defaults to
#'   \code{reference$features}.
#'
#' @details
#' \code{somalign_query_from_som()} reuses \code{som$unit.classif} directly as
#' the per-cell query-node assignment, so the O(N \eqn{\times} nodes) nearest-
#' code search is skipped entirely.  This is exact when the supplied
#' \code{codebook} was used as-is during SOM training.  If the codebook has
#' been post-processed (e.g.\ winsorised and rescaled into reference space),
#' the reused assignments are an approximation: cells near node boundaries may
#' flip under the non-linear transform, but in practice this affects only a
#' small fraction of cells at the distributional tails and does not measurably
#' impact OT alignment quality.
#'
#' Unlike \code{\link{somalign_query}()}, \code{sample_distance} is not
#' recomputed and is set to \code{NA} in the returned object.  The field is not
#' used by \code{\link{somalign_fit}()}.
#'
#' The returned object is identical in structure to a \code{somalign_query}
#' produced by \code{somalign_query()}, and is fully compatible with
#' \code{\link{somalign_fit}()}.
#'
#' @return A \code{somalign_query} object.
#'
#' @seealso [somalign_query()], [somalign_reference_from_som()], [somalign_fit()]
#'
#' @export
somalign_query_from_som <- function(som,
                                    data,
                                    reference,
                                    codebook = NULL,
                                    codebook_space = c("reference_scaled", "raw"),
                                    features = NULL) {
  .somalign_check_reference(reference)
  codebook_space <- match.arg(codebook_space)
  features <- .somalign_query_features(features, reference)
  data <- .somalign_prepare_feature_matrix(data, features, what = "query data")
  scaled_data <- .somalign_scale_matrix(data, reference$center, reference$scale)

  # --- SOM unit assignments --------------------------------------------------
  unit <- som$unit.classif
  if (is.null(unit) || length(unit) == 0L) {
    stop(
      "`som$unit.classif` is absent or empty. ",
      "The SOM must have been trained with `keep.data = TRUE`, ",
      "or `unit.classif` must otherwise be populated.",
      call. = FALSE
    )
  }
  if (length(unit) != nrow(data)) {
    stop(
      "`som$unit.classif` has ", length(unit), " entries but `data` has ",
      nrow(data), " rows. They must match.",
      call. = FALSE
    )
  }

  # --- Codebook --------------------------------------------------------------
  if (is.null(codebook)) {
    codebook <- .somalign_get_codebook(som, features = features, what = "som_query")
  } else {
    codebook <- as.matrix(codebook)
    storage.mode(codebook) <- "double"
    missing_f <- setdiff(features, colnames(codebook))
    if (length(missing_f) > 0L) {
      stop(
        "`codebook` is missing feature(s): ",
        paste(missing_f, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!identical(colnames(codebook), features)) {
      codebook <- codebook[, features, drop = FALSE]
    }
  }
  if (identical(codebook_space, "raw")) {
    codebook <- .somalign_scale_matrix(codebook, reference$center, reference$scale)
  }

  n_nodes <- nrow(codebook)
  if (max(unit) > n_nodes) {
    stop(
      "`som$unit.classif` contains indices up to ", max(unit),
      " but the codebook has only ", n_nodes, " rows.",
      call. = FALSE
    )
  }

  node_masses <- .somalign_node_masses(as.integer(unit), n_nodes)

  sample_id <- rownames(data)
  if (is.null(sample_id)) sample_id <- as.character(seq_len(nrow(data)))

  structure(
    list(
      data             = data,
      scaled_data      = scaled_data,
      som_query        = som,
      codebook         = codebook,
      node_masses      = node_masses,
      sample_unit      = as.integer(unit),
      sample_distance  = rep(NA_real_, nrow(data)),
      sample_id        = sample_id,
      reference_features = reference$features
    ),
    class = "somalign_query"
  )
}

#' Global pre-correction of query data to match the reference distribution
#'
#' An optional pre-processing step that shifts (or shifts and scales) query
#' data in reference-scaled coordinate space so that the per-marker query
#' means align with the reference coordinate origin. Passing the returned
#' matrix to [somalign_query()] trains the query SOM on pre-centred data,
#' reducing the global component of the batch shift that [somalign_fit()] then
#' needs to resolve via optimal transport.
#'
#' @param data Numeric matrix of query data, same format as the `data` argument
#'   to [somalign_query()].
#' @param reference A `somalign_reference` object.
#' @param method Normalisation method. `"mean"` (default) subtracts the
#'   per-marker query mean in reference-scaled space, removing a uniform
#'   location shift. `"scale"` additionally divides by the per-marker query
#'   standard deviation, removing a uniform scale shift as well.
#' @param features Optional character vector of feature names. Defaults to
#'   `reference$features`.
#'
#' @return A numeric matrix with the same number of rows as `data` and columns
#'   in `reference$features` order, expressed in the original (unscaled) units
#'   of `data`. Pass this matrix directly to [somalign_query()].
#'
#' @details
#' `somalign_normalize()` applies the same per-marker shift (and optionally
#' rescaling) to every cell, so population-specific batch effects remain for
#' [somalign_fit()] to resolve.
#'
#' **When not to use this function.** Mean-normalisation assumes the apparent
#' per-marker shift reflects instrument drift or reagent-lot differences
#' affecting all populations uniformly. If the shift reflects genuine
#' compositional differences between batches (e.g.\ different cell-type
#' frequencies), subtracting the global mean distorts the biology.
#'
#' @seealso [somalign_query()], [somalign_fit()]
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' shifted <- mat + 0.5
#' corrected <- somalign_normalize(shifted, ref)
#' @export
somalign_normalize <- function(data, reference,
                               method = c("mean", "scale"),
                               features = NULL) {
  .somalign_check_reference(reference)
  method <- match.arg(method)
  if (is.null(features)) features <- reference$features
  data <- .somalign_prepare_feature_matrix(data, features, what = "query data")

  scaled <- .somalign_scale_matrix(data, reference$center, reference$scale)

  query_means <- colMeans(scaled)
  scaled <- sweep(scaled, 2, query_means, "-")

  if (identical(method, "scale")) {
    sds <- apply(scaled, 2, stats::sd)
    sds[sds < .Machine$double.eps] <- 1
    scaled <- sweep(scaled, 2, sds, "/")
  }

  sweep(sweep(scaled, 2, reference$scale, "*"), 2, reference$center, "+")
}
