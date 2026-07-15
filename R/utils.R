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

# ---------------------------------------------------------------------------
# Low-level scalar / type validators
# ---------------------------------------------------------------------------

.somalign_check_pos_scalar <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
    stop("`", nm, "` must be a single positive finite number.", call. = FALSE)
  invisible(x)
}

.somalign_check_nonneg_scalar <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0)
    stop("`", nm, "` must be a single non-negative finite number.", call. = FALSE)
  invisible(x)
}

.somalign_check_prob_scalar <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1)
    stop("`", nm, "` must be a single number in [0, 1].", call. = FALSE)
  invisible(x)
}

.somalign_check_pos_int <- function(x, nm, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(x))
  if (allow_null && length(x) == 1L && is.numeric(x) && is.infinite(x) && x > 0)
    return(invisible(x))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != round(x) || x < 1L)
    stop("`", nm, "` must be a single positive integer.", call. = FALSE)
  invisible(x)
}

.somalign_check_flag <- function(x, nm) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    stop("`", nm, "` must be TRUE or FALSE.", call. = FALSE)
  invisible(x)
}

.somalign_check_data_arg <- function(x, what = "data") {
  if (!is.matrix(x) && !is.data.frame(x))
    stop("`", what, "` must be a numeric matrix or data frame.", call. = FALSE)
  invisible(x)
}

.somalign_check_opt_char <- function(x, what) {
  if (!is.null(x) && (!is.character(x) || length(x) == 0L))
    stop("`", what, "` must be a non-empty character vector or NULL.", call. = FALSE)
  invisible(x)
}

.somalign_check_numeric_vec <- function(x, what) {
  if (!is.numeric(x) || length(x) == 0L || !all(is.finite(x)))
    stop("`", what, "` must be a finite numeric vector.", call. = FALSE)
  invisible(x)
}

.somalign_check_opt_grid <- function(grid) {
  if (!is.null(grid) && !inherits(grid, "somgrid"))
    stop("`grid` must be a `kohonen::somgrid()` object or NULL.", call. = FALSE)
  invisible(grid)
}

.somalign_check_unit_scalar <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 1)
    stop("`", nm, "` must be a single number in (0, 1].", call. = FALSE)
  invisible(x)
}

# Compute a low-rank batch subspace via variance-threshold SVD.
# M: n_obs × p matrix (displacement vectors); weights: optional length-n_obs mass weights.
# Returns list(V = p × r, rank = r, variance_explained = cumvar[r]).
.somalign_subspace_svd <- function(M, variance_threshold, weights = NULL) {
  if (!is.null(weights) && length(weights) > 0 && sum(weights) > 0) {
    w <- sqrt(weights / sum(weights))
    M <- M * w
  }
  tot <- sum(M^2)
  if (!is.finite(tot) || tot == 0) {
    V <- matrix(0, nrow = ncol(M), ncol = 1L)
    return(list(V = V, rank = 1L, variance_explained = 1))
  }
  sv <- svd(M, nu = 0L)
  cumvar <- cumsum(sv$d^2) / sum(sv$d^2)
  r <- which(cumvar >= variance_threshold)
  r <- if (length(r) == 0L) ncol(M) else r[[1L]]
  r <- max(1L, r)
  list(V = sv$v[, seq_len(r), drop = FALSE], rank = r, variance_explained = cumvar[[r]])
}

# Only validated when solver = "annealing"; ignored (and unvalidated)
# otherwise, since the args are inert for every other solver.
.somalign_check_anneal_params <- function(anneal_start, anneal_factor, anneal_stages) {
  if (!is.numeric(anneal_start) || length(anneal_start) != 1L ||
      !is.finite(anneal_start) || anneal_start < 1)
    stop("`anneal_start` must be >= 1.", call. = FALSE)
  if (!is.null(anneal_factor) &&
      (!is.numeric(anneal_factor) || length(anneal_factor) != 1L ||
       !is.finite(anneal_factor) || anneal_factor <= 0 || anneal_factor >= 1))
    stop("`anneal_factor` must be NULL or a single number in (0, 1) ",
         "(a per-stage cooling ratio).", call. = FALSE)
  .somalign_check_pos_int(anneal_stages, "anneal_stages")
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Bundle validators (one call per exported function)
# ---------------------------------------------------------------------------

# Shared OT tuning parameters for somalign_fit, somalign_fit_two_pass, and
# somalign_fit_anchored. label_guided is optional (anchored fit omits it).
.somalign_check_fit_params <- function(rho_query, rho_ref,
                                       min_match_fraction, confidence_threshold,
                                       correction_min_mass, max_iter, tol,
                                       chunk_size, label_guided = NULL) {
  .somalign_check_pos_scalar(rho_query, "rho_query")
  .somalign_check_pos_scalar(rho_ref, "rho_ref")
  .somalign_check_prob_scalar(min_match_fraction, "min_match_fraction")
  .somalign_check_prob_scalar(confidence_threshold, "confidence_threshold")
  .somalign_check_nonneg_scalar(correction_min_mass, "correction_min_mass")
  .somalign_check_pos_int(max_iter, "max_iter")
  .somalign_check_pos_scalar(tol, "tol")
  .somalign_check_pos_int(chunk_size, "chunk_size", allow_null = TRUE)
  if (!is.null(label_guided)) .somalign_check_flag(label_guided, "label_guided")
  invisible(NULL)
}

.somalign_check_som_train_args <- function(data, labels, features,
                                           grid, rlen, alpha,
                                           data_what = "data") {
  .somalign_check_data_arg(data, what = data_what)
  .somalign_check_opt_char(labels,   what = "labels")
  .somalign_check_opt_char(features, what = "features")
  .somalign_check_opt_grid(grid)
  .somalign_check_pos_int(rlen, "rlen")
  .somalign_check_numeric_vec(alpha, "alpha")
  invisible(NULL)
}

.somalign_check_reference_args <- function(data, labels, features,
                                           quantile_probs) {
  .somalign_check_data_arg(data, what = "data")
  .somalign_check_opt_char(labels,   what = "labels")
  .somalign_check_opt_char(features, what = "features")
  .somalign_check_numeric_vec(quantile_probs, "quantile_probs")
  invisible(NULL)
}

.somalign_check_reference_from_som_args <- function(quantile_probs,
                                                    distance_chunk_size) {
  .somalign_check_numeric_vec(quantile_probs, "quantile_probs")
  .somalign_check_pos_int(distance_chunk_size, "distance_chunk_size")
  invisible(NULL)
}

.somalign_check_query_from_som_args <- function(data, features) {
  .somalign_check_data_arg(data, what = "data")
  .somalign_check_opt_char(features, what = "features")
  invisible(NULL)
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

# Squared Euclidean distance. Used only to build the OT cost matrix
# (fit.R): squared cost makes the barycentric correction a Brenier optimal
# transport map. Cell-to-node projection distances use .somalign_nearest_code,
# which returns plain Euclidean distances for the distance-quantile thresholds.
.somalign_pairwise_distance <- function(x, y) {
  d2 <- outer(rowSums(x * x), rowSums(y * y), "+") - 2 * tcrossprod(x, y)
  pmax(d2, 0)
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
  # Only genuinely absent thresholds (NA) fall back to the global quantile.
  # An explicit Inf means "no finite threshold" (never flag this node) and must
  # be preserved: distance > Inf is FALSE, so such cells are never marked
  # outside the reference. `is.na(Inf)` is FALSE, so Inf survives here.
  missing <- is.na(thresholds)
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

# ---------------------------------------------------------------------------
# Helpers for somalign_reference_from_som()
# ---------------------------------------------------------------------------

# Extract the Y-layer (label) codebook from a kohonen xyf/supersom object.
# Returns NULL with a message when no second code layer is present (plain som).
.somalign_extract_label_codes <- function(som) {
  if (!is.list(som)) {
    return(NULL)
  }
  codes <- som$codes
  if (is.null(codes) || !is.list(codes) || length(codes) < 2L) {
    message(
      "somalign_reference_from_som: SOM has no second code layer; ",
      "label transfer will be disabled."
    )
    return(NULL)
  }
  yc <- codes[[2L]]
  if (!is.matrix(yc) && !is.data.frame(yc)) {
    message(
      "somalign_reference_from_som: second code layer is not a matrix; ",
      "label transfer will be disabled."
    )
    return(NULL)
  }
  yc <- as.matrix(yc)
  storage.mode(yc) <- "double"
  yc
}

# Extract the X-layer training data from a kohonen object.
# For supersom/xyf: som$data is a list; takes [[1]].
# For plain som:    som$data is also a list (one element).
# Errors clearly when data were not retained (keep.data = FALSE).
.somalign_extract_som_data <- function(som) {
  d <- som[["data"]]
  if (is.null(d)) {
    stop(
      "The SOM does not store training data (`som$data` is NULL). ",
      "Retrain with `keep.data = TRUE` (the kohonen default) or supply ",
      "distance quantiles manually via `somalign_reference_from_nodes()`.",
      call. = FALSE
    )
  }
  if (is.list(d)) {
    d <- d[[1L]]
  }
  d <- as.matrix(d)
  storage.mode(d) <- "double"
  d
}

# Compute per-cell Euclidean distance to each cell's already-known assigned
# node.  O(N * p) — no 900-way argmax, no O(N * nodes) matrix.
# Processed in chunks to bound peak memory.
.somalign_som_cell_distances <- function(X, codebook, unit, chunk_size) {
  n <- nrow(X)
  d <- numeric(n)
  chunk_size <- max(1L, as.integer(chunk_size))
  for (s in seq(1L, n, by = chunk_size)) {
    idx <- s:min(s + chunk_size - 1L, n)
    d[idx] <- sqrt(rowSums(
      (X[idx, , drop = FALSE] - codebook[unit[idx], , drop = FALSE]) ^ 2
    ))
  }
  d
}
