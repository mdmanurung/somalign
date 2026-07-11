.somalign_check_reference <- function(reference) {
  if (!inherits(reference, "somalign_reference")) {
    stop("`reference` must be a somalign_reference object.", call. = FALSE)
  }
  invisible(reference)
}

.somalign_check_query <- function(query) {
  if (!inherits(query, "somalign_query")) {
    stop("`query` must be a somalign_query object.", call. = FALSE)
  }
  invisible(query)
}

.somalign_as_matrix <- function(x, what = "data") {
  if (is.data.frame(x)) {
    non_numeric <- !vapply(x, is.numeric, logical(1))
    if (any(non_numeric)) {
      stop("`", what, "` must contain only numeric feature columns.", call. = FALSE)
    }
    x <- as.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  x
}

.somalign_validate_feature_names <- function(x, what = "data") {
  names <- colnames(x)
  if (is.null(names) || any(is.na(names)) || any(!nzchar(names))) {
    stop("`", what, "` must have non-empty column names.", call. = FALSE)
  }
  duplicated_names <- unique(names[duplicated(names)])
  if (length(duplicated_names) > 0) {
    stop(
      "Duplicated feature names in `", what, "`: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(names)
}

.somalign_select_features <- function(x, features, what = "data") {
  if (is.null(features)) {
    features <- colnames(x)
  }
  if (!is.character(features) || length(features) == 0) {
    stop("`features` must be a non-empty character vector.", call. = FALSE)
  }
  duplicated_features <- unique(features[duplicated(features)])
  if (length(duplicated_features) > 0) {
    stop(
      "Duplicated requested features: ",
      paste(duplicated_features, collapse = ", "),
      call. = FALSE
    )
  }
  missing <- setdiff(features, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Missing features in `", what, "`: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  x[, features, drop = FALSE]
}

.somalign_validate_finite <- function(x, what = "data") {
  if (!all(is.finite(x))) {
    stop("`", what, "` must contain only finite values.", call. = FALSE)
  }
  invisible(x)
}

.somalign_named_numeric <- function(x, features, what) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != length(features)) {
    stop("`", what, "` must be numeric with one value per feature.", call. = FALSE)
  }
  if (is.null(names(x))) {
    names(x) <- features
  }
  missing <- setdiff(features, names(x))
  if (length(missing) > 0) {
    stop(
      "`", what, "` is missing feature names: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  x <- x[features]
  if (!all(is.finite(x))) {
    stop("`", what, "` must contain only finite values.", call. = FALSE)
  }
  x
}

.somalign_compute_scaling <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  .somalign_validate_scale(scale)
  list(center = center, scale = scale)
}

.somalign_validate_scale <- function(scale) {
  if (!all(is.finite(scale)) || any(scale <= 0)) {
    stop("Reference features must have finite non-zero variance.", call. = FALSE)
  }
  invisible(scale)
}

.somalign_validate_codebook_space <- function(codebook_space) {
  if (is.null(codebook_space)) {
    stop(
      "`codebook_space` must be specified for existing reference SOMs. ",
      "Use \"reference_scaled\" for SOMs trained on reference-scaled data or ",
      "\"raw\" for SOMs trained on raw feature values.",
      call. = FALSE
    )
  }
  match.arg(codebook_space, c("reference_scaled", "raw"))
}

.somalign_scale_matrix <- function(x, center, scale) {
  sweep(sweep(x, 2, center, "-"), 2, scale, "/")
}

.somalign_prepare_feature_matrix <- function(x, features = NULL, what = "data") {
  x <- .somalign_as_matrix(x, what = what)
  .somalign_validate_feature_names(x, what = what)
  x <- .somalign_select_features(x, features, what = what)
  .somalign_validate_finite(x, what = what)
  x
}

.somalign_get_codebook <- function(som, features = NULL, what = "som") {
  codes <- .somalign_extract_codes(som, what)

  if (!is.null(features)) {
    if (is.null(colnames(codes))) {
      stop(
        "`", what, "` codebook must have column names matching the reference features.",
        call. = FALSE
      )
    }
    missing <- setdiff(features, colnames(codes))
    if (length(missing) > 0) {
      stop(
        "SOM codebook is missing features: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    codes <- codes[, features, drop = FALSE]
  } else {
    .somalign_validate_feature_names(codes, what = paste0(what, " codebook"))
  }
  .somalign_validate_finite(codes, what = paste0(what, " codebook"))
  codes
}

.somalign_extract_codes <- function(som, what) {
  codes <- NULL

  if (is.matrix(som) || is.data.frame(som)) {
    codes <- as.matrix(som)
  }

  if (is.null(codes) && requireNamespace("kohonen", quietly = TRUE)) {
    codes <- tryCatch(kohonen::getCodes(som), error = function(e) NULL)
  }

  if (is.null(codes) && !is.null(som$codes)) {
    codes <- som$codes
  }

  if (is.list(codes)) {
    if ("data" %in% names(codes)) {
      codes <- codes[["data"]]
    } else {
      codes <- codes[[1]]
    }
  }

  if (is.null(codes)) {
    stop("Could not extract a SOM codebook from `", what, "`.", call. = FALSE)
  }

  codes <- .somalign_as_matrix(codes, what = paste0(what, " codebook"))
  codes
}

.somalign_default_grid <- function(n_samples) {
  side <- max(2L, min(8L, ceiling(sqrt(max(4L, n_samples)) / 2)))
  kohonen::somgrid(xdim = side, ydim = side, topo = "hexagonal")
}

.somalign_nearest_code <- function(x, codebook) {
  x <- as.matrix(x)
  codebook <- as.matrix(codebook)
  d2 <- outer(rowSums(x * x), rowSums(codebook * codebook), "+") -
    2 * tcrossprod(x, codebook)
  d2 <- pmax(d2, 0)
  unit <- max.col(-d2, ties.method = "first")
  distance <- sqrt(d2[cbind(seq_len(nrow(d2)), unit)])
  list(unit = as.integer(unit), distance = as.numeric(distance))
}

.somalign_nearest_code_chunked <- function(x, codebook, chunk_size = 10000L) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (is.null(chunk_size) || is.infinite(chunk_size) || n == 0L || chunk_size >= n) {
    return(.somalign_nearest_code(x, codebook))
  }
  chunk_size <- as.integer(chunk_size)
  unit <- integer(n)
  distance <- numeric(n)
  for (s in seq(1L, n, by = chunk_size)) {
    idx <- s:min(s + chunk_size - 1L, n)
    res <- .somalign_nearest_code(x[idx, , drop = FALSE], codebook)
    unit[idx] <- res$unit
    distance[idx] <- res$distance
  }
  list(unit = unit, distance = distance)
}

.somalign_pairwise_distance <- function(x, y) {
  d2 <- outer(rowSums(x * x), rowSums(y * y), "+") - 2 * tcrossprod(x, y)
  sqrt(pmax(d2, 0))
}

.somalign_normalize_masses <- function(x, n, what) {
  if (is.null(x)) {
    x <- rep(1 / n, n)
  }
  if (!is.numeric(x) || length(x) != n || any(!is.finite(x)) || any(x < 0)) {
    stop("`", what, "` must be a non-negative finite numeric vector.", call. = FALSE)
  }
  total <- sum(x)
  if (total <= 0) {
    rep(1 / n, n)
  } else {
    as.numeric(x / total)
  }
}

.somalign_node_masses <- function(units, n_nodes) {
  counts <- tabulate(units, nbins = n_nodes)
  if (sum(counts) == 0) {
    rep(0, n_nodes)
  } else {
    counts / sum(counts)
  }
}

.somalign_quantile_names <- function(probs) {
  paste0(format(100 * probs, trim = TRUE, scientific = FALSE), "%")
}

.somalign_distance_quantiles <- function(distances, units, n_nodes, probs) {
  names <- .somalign_quantile_names(probs)
  global <- stats::quantile(distances, probs = probs, names = FALSE, type = 7)
  names(global) <- names
  split_distances <- split(distances, factor(units, levels = seq_len(n_nodes)))
  node_rows <- lapply(split_distances, function(d) {
    if (length(d) == 0L) {
      global
    } else {
      stats::quantile(d, probs = probs, names = FALSE, type = 7)
    }
  })
  out <- matrix(unlist(node_rows, use.names = FALSE), nrow = n_nodes, ncol = length(probs), byrow = TRUE)
  colnames(out) <- names
  list(node = out, global = global)
}

.somalign_thresholds <- function(reference, units, column = "95%") {
  q <- reference$distance_quantiles
  if (is.null(q) || ncol(q) == 0) {
    return(rep(NA_real_, length(units)))
  }
  if (!column %in% colnames(q)) {
    warning(
      sprintf(
        "Column '%s' not found in distance_quantiles; using '%s' instead.",
        column, colnames(q)[ncol(q)]
      ),
      call. = FALSE
    )
    column <- colnames(q)[ncol(q)]
  }
  thresholds <- q[units, column]
  missing <- !is.finite(thresholds)
  if (any(missing)) {
    global <- reference$global_distance_quantiles
    fallback <- if (!is.null(global) && column %in% names(global)) global[[column]] else NA_real_
    thresholds[missing] <- fallback
  }
  as.numeric(thresholds)
}

.somalign_label_probabilities <- function(labels, units, n_nodes) {
  if (is.null(labels)) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  if (length(labels) != length(units)) {
    stop("`labels` must have one value per row of reference data.", call. = FALSE)
  }
  labels <- as.character(labels)
  levels <- sort(unique(labels[!is.na(labels)]))
  if (length(levels) == 0) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  valid <- !is.na(labels)
  node_idx <- units[valid]
  lbl_idx <- as.integer(factor(labels[valid], levels = levels))
  combined_idx <- (lbl_idx - 1L) * n_nodes + node_idx
  raw_counts <- tabulate(combined_idx, nbins = n_nodes * length(levels))
  out <- matrix(raw_counts, nrow = n_nodes, ncol = length(levels))
  colnames(out) <- levels
  row_totals <- rowSums(out)
  nonzero <- row_totals > 0
  out[nonzero, ] <- out[nonzero, , drop = FALSE] / row_totals[nonzero]
  out
}

.somalign_normalize_label_prob <- function(label_prob, n_nodes) {
  if (is.null(label_prob)) {
    return(matrix(numeric(0), nrow = n_nodes, ncol = 0))
  }
  label_prob <- as.matrix(label_prob)
  storage.mode(label_prob) <- "double"
  if (nrow(label_prob) != n_nodes) {
    stop("`label_prob` must have one row per reference node.", call. = FALSE)
  }
  if (is.null(colnames(label_prob))) {
    colnames(label_prob) <- paste0("label_", seq_len(ncol(label_prob)))
  }
  if (any(!is.finite(label_prob)) || any(label_prob < 0)) {
    stop("`label_prob` must contain non-negative finite values.", call. = FALSE)
  }
  row_totals <- rowSums(label_prob)
  nonzero <- row_totals > 0
  label_prob[nonzero, ] <- label_prob[nonzero, , drop = FALSE] / row_totals[nonzero]
  label_prob
}

.somalign_prepare_distance_quantiles <- function(distance_quantiles, n_nodes) {
  if (is.null(distance_quantiles)) {
    out <- matrix(NA_real_, nrow = n_nodes, ncol = 4)
    colnames(out) <- c("50%", "90%", "95%", "99%")
    return(out)
  }
  distance_quantiles <- as.matrix(distance_quantiles)
  storage.mode(distance_quantiles) <- "double"
  if (nrow(distance_quantiles) != n_nodes) {
    stop("`distance_quantiles` must have one row per reference node.", call. = FALSE)
  }
  if (is.null(colnames(distance_quantiles))) {
    colnames(distance_quantiles) <- paste0(seq_len(ncol(distance_quantiles)))
  }
  if (any(!is.finite(distance_quantiles) & !is.infinite(distance_quantiles))) {
    stop("`distance_quantiles` must be finite or Inf.", call. = FALSE)
  }
  distance_quantiles
}

.somalign_entropy <- function(prob) {
  prob <- prob[is.finite(prob) & prob > 0]
  if (length(prob) == 0) {
    return(NA_real_)
  }
  -sum(prob * log(prob))
}
