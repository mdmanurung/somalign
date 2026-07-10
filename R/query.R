#' Prepare query data and attach or train a query SOM
#'
#' Query data are always scaled with the saved reference center and scale.
#'
#' @param data Numeric query data.
#' @param reference A `somalign_reference` object.
#' @param som_query Optional query SOM or SOM-like object with a codebook. The
#'   codebook must be in the reference-scaled feature space, i.e. trained on
#'   query data transformed with `reference$center` and `reference$scale`.
#' @param grid Optional `kohonen::somgrid()` object when `som_query` is omitted.
#' @param rlen Number of SOM training iterations passed to [kohonen::som()].
#' @param alpha Learning-rate schedule passed to [kohonen::som()].
#' @param features Optional feature names. Defaults to the reference feature
#'   order.
#' @param ... Additional arguments passed to [kohonen::som()].
#'
#' @return A `somalign_query` object.
#' @examples
#' \dontrun{
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"))
#' qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"))
#' }
#' @export
somalign_query <- function(data,
                           reference,
                           som_query = NULL,
                           grid = NULL,
                           rlen = 100,
                           alpha = c(0.05, 0.01),
                           features = NULL,
                           ...) {
  .somalign_check_reference(reference)
  if (is.null(features)) {
    features <- reference$features
  }
  if (!identical(features, reference$features)) {
    missing <- setdiff(reference$features, features)
    if (length(missing) > 0 || length(features) != length(reference$features)) {
      stop("`features` must match the reference feature set.", call. = FALSE)
    }
    features <- reference$features
  }

  data <- .somalign_prepare_feature_matrix(data, features, what = "query data")
  scaled_data <- .somalign_scale_matrix(data, reference$center, reference$scale)

  if (is.null(som_query)) {
    if (is.null(grid)) {
      grid <- .somalign_default_grid(nrow(data))
    }
    som_query <- kohonen::som(scaled_data, grid = grid, rlen = rlen, alpha = alpha, ...)
  }

  codebook <- .somalign_get_codebook(som_query, features = reference$features, what = "som_query")
  sample_map <- .somalign_nearest_code(scaled_data, codebook)
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
