#' Cross-batch novelty gate for somalign
#'
#' Detects candidate novel cell populations in a multi-batch query by combining
#' a continuous per-cell novelty score with a cross-batch reproducibility
#' requirement.  A population is only minted as a novel candidate if it recurs
#' across \code{>= min_batches} independent query batches at consistent
#' reference-scaled coordinates.  Single-batch minting is intentionally
#' disallowed: a per-batch artifact (e.g., batch-specific displaced cells)
#' lands at different coordinates in each batch and is therefore rejected,
#' whereas a genuine biological novelty produces a consistent signal across
#' batches.
#'
#' @section Novelty score:
#'   The per-cell novelty score is the mean Euclidean distance from each query
#'   cell to its \code{k} nearest reference codebook nodes, computed in
#'   reference-scaled space.  This continuous score (higher = more novel) is
#'   empirically validated to separate novel from in-distribution cells with
#'   high recall at low FPR.  A boolean "outside-reference" flag is
#'   intentionally \emph{not} used: per-node thresholds loosen at the reference
#'   periphery and mask genuine novelty.
#'
#' @section Cross-batch safety gate:
#'   After identifying high-novelty tail cells, the algorithm clusters them
#'   within each batch and then looks for clusters that appear at consistent
#'   reference-scaled coordinates across \code{>= min_batches} batches.  Only
#'   cross-batch reproducible clusters are minted as candidates.  This rejects
#'   single-batch artifacts and technical outliers.
#'
#' @section Prototype deduplication:
#'   A spread-out novel population (e.g., NKG2C+ NK cells spanning a loose
#'   manifold) can yield multiple near-duplicate minted prototypes — separate
#'   cross-group connected components that happen to sit close together in
#'   reference-scaled space.  The post-mint dedup step agglomeratively merges
#'   such near-duplicates using single-linkage hierarchical clustering at
#'   height \code{merge_tol_factor * node_spacing}, where \code{node_spacing}
#'   is the median nearest-neighbour distance among codebook nodes.
#'
#'   \code{merge_tol_factor} (default 3.0) should be larger than
#'   \code{tol_factor} (default 1.5): \code{tol_factor} controls tight
#'   cross-batch matching ("same population across batches"), while
#'   \code{merge_tol_factor} controls the looser dedup criterion ("these
#'   minted candidates represent the same population").  A warning is emitted
#'   if \code{merge_tol_factor <= tol_factor}.
#'
#'   For each merged group of prototypes: the new prototype coordinate is the
#'   size-weighted mean of member coordinates; \code{size} is the sum of member
#'   sizes; and \code{n_groups_support} counts the number of DISTINCT groups in
#'   the UNION of member group sets (not the sum of individual counts — a group
#'   that supports two merged members is counted once).
#'
#'   Set \code{merge_tol_factor = NULL} or \code{0} to disable deduplication
#'   and reproduce the pre-dedup behaviour.
#'
#' @param fit A \code{somalign_fit} object.  Uses
#'   \code{fit$query$scaled_data} (cells x features, reference-scaled space)
#'   and \code{fit$reference$codebook}.
#' @param group A vector of length \code{nrow(fit$query$scaled_data)} assigning
#'   each query cell to a batch or sample (the reproducibility unit).  The
#'   unique values define the independent groups used for cross-batch matching.
#' @param score Optional numeric vector of length \code{nrow(fit$query$scaled_data)}
#'   providing a precomputed per-cell novelty score (higher = more novel).  If
#'   \code{NULL} (default), the score is computed as the mean Euclidean
#'   distance to the \code{k} nearest reference codebook nodes in
#'   reference-scaled space.
#' @param k Integer.  Number of nearest reference codebook nodes used to
#'   compute the novelty score.  Default \code{8L}.  Automatically clamped to
#'   \code{nrow(reference$codebook)} if the codebook is small.
#' @param tail_quantile Numeric in \code{(0, 1)}.  Cells with novelty score
#'   \code{>= quantile(score, tail_quantile)} are flagged as the high-novelty
#'   tail and clustered.  Default \code{0.90}.
#' @param n_clusters Integer.  Maximum number of k-means clusters per group in
#'   the tail.  The actual number is adaptive: it is capped so each cluster
#'   can plausibly contain at least \code{min_cluster} cells, and capped at
#'   the number of distinct coordinate rows (to prevent k-means from receiving
#'   more centers than distinct points).  Default \code{12L}.
#' @param min_cluster Integer.  Minimum number of tail cells a cluster must
#'   contain to be retained as a per-group candidate centroid.  Default
#'   \code{50L}.
#' @param min_batches Integer.  Minimum number of distinct groups a matched
#'   cluster component must span before being minted as a novel candidate.
#'   Default \code{2L}.  Set to \code{1} only if you have a single batch and
#'   understand the artifact risk.
#' @param tol_factor Positive numeric.  Two per-group centroids are linked in
#'   the cross-batch matching graph if their Euclidean distance (in
#'   reference-scaled space) is \code{<= tol_factor * reference_node_spacing},
#'   where \code{reference_node_spacing} is the median nearest-neighbour
#'   distance among codebook nodes.  Default \code{1.5}.
#' @param merge_tol_factor Positive numeric or \code{NULL}/\code{0} to disable.
#'   After cross-batch minting, minted prototypes closer than
#'   \code{merge_tol_factor * reference_node_spacing} (Euclidean) are
#'   agglomeratively merged (single-linkage) into a single deduplicated
#'   candidate.  Should be larger than \code{tol_factor} (a warning is emitted
#'   if not).  Default \code{3.0}.
#' @param chunk_size Integer.  Number of cells processed per chunk when
#'   computing the novelty score.  Increase for speed; decrease if memory is
#'   constrained.  Default \code{100000L}.
#'
#' @return An object of class \code{somalign_novelty_candidates}, a list with:
#'   \describe{
#'     \item{\code{prototypes}}{Matrix \code{[n_candidates x features]} of
#'       minted candidate centroids in reference-scaled space (colnames =
#'       feature names).  Has zero rows when no candidates are minted.}
#'     \item{\code{n_groups_support}}{Integer vector (length n_candidates):
#'       the number of distinct groups each candidate was observed in (after
#'       dedup, this is the union of the group sets of merged prototypes).}
#'     \item{\code{size}}{Integer vector (length n_candidates): total number
#'       of tail cells contributing to each candidate across all groups.}
#'     \item{\code{score}}{Numeric vector (length n_cells): the per-cell
#'       novelty score.}
#'     \item{\code{tail}}{Logical vector (length n_cells): \code{TRUE} for
#'       cells in the high-novelty tail.}
#'     \item{\code{n_groups}}{Integer: number of distinct groups.}
#'     \item{\code{params}}{Named list of the parameters used, including
#'       \code{merge_tol_factor}.}
#'   }
#'
#' @section Composing with somalign_extend_reference:
#'   The \code{prototypes} matrix is ready to pass directly to
#'   \code{\link{somalign_extend_reference}} as \code{new_codebook}:
#'   \preformatted{
#'   cand <- somalign_novelty_candidates(fit, group)
#'   extended_ref <- somalign_extend_reference(
#'     fit$reference,
#'     cand$prototypes,
#'     new_labels = paste0("novel_", seq_len(nrow(cand$prototypes)))
#'   )
#'   }
#'
#' @seealso \code{\link{somalign_extend_reference}} for grafting prototypes
#'   into the reference.
#'
#' @export
somalign_novelty_candidates <- function(
    fit,
    group,
    score            = NULL,
    k                = 8L,
    tail_quantile    = 0.90,
    n_clusters       = 12L,
    min_cluster      = 50L,
    min_batches      = 2L,
    tol_factor       = 1.5,
    merge_tol_factor = 3.0,
    chunk_size       = 100000L
) {
  # ---- Input validation -------------------------------------------------------
  if (!inherits(fit, "somalign_fit")) {
    stop("'fit' must be a somalign_fit object", call. = FALSE)
  }
  scaled_data <- fit$query$scaled_data
  codebook    <- fit$reference$codebook
  n_cells     <- nrow(scaled_data)
  n_features  <- ncol(scaled_data)
  features    <- colnames(scaled_data)

  if (length(group) != n_cells) {
    stop(
      sprintf(
        "'group' length (%d) must equal nrow(fit$query$scaled_data) (%d)",
        length(group), n_cells
      ),
      call. = FALSE
    )
  }
  if (!is.null(score)) {
    if (length(score) != n_cells) {
      stop(
        sprintf(
          "precomputed 'score' length (%d) must equal nrow(fit$query$scaled_data) (%d)",
          length(score), n_cells
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(score) || !all(is.finite(score))) {
      stop("'score' must be a finite numeric vector", call. = FALSE)
    }
  }

  .somalign_check_pos_int(k,           "k")
  .somalign_check_prob_scalar(tail_quantile, "tail_quantile")
  if (tail_quantile <= 0 || tail_quantile >= 1) {
    stop("'tail_quantile' must be in (0, 1)", call. = FALSE)
  }
  .somalign_check_pos_int(n_clusters,  "n_clusters")
  .somalign_check_pos_int(min_cluster, "min_cluster")
  .somalign_check_pos_int(min_batches, "min_batches")
  .somalign_check_pos_scalar(tol_factor, "tol_factor")
  .somalign_check_pos_int(chunk_size,  "chunk_size")

  # Validate merge_tol_factor: NULL or 0 = disabled; otherwise positive numeric.
  dedup_enabled <- TRUE
  if (is.null(merge_tol_factor) ||
      (is.numeric(merge_tol_factor) && length(merge_tol_factor) == 1L &&
       is.finite(merge_tol_factor) && merge_tol_factor == 0)) {
    dedup_enabled <- FALSE
    merge_tol_factor <- NULL
  } else {
    .somalign_check_pos_scalar(merge_tol_factor, "merge_tol_factor")
    if (merge_tol_factor <= tol_factor) {
      warning(
        sprintf(
          paste0(
            "`merge_tol_factor` (%.3g) <= `tol_factor` (%.3g): ",
            "dedup window is tighter than the cross-batch matching window. ",
            "Consider merge_tol_factor > tol_factor (e.g., %.3g)."
          ),
          merge_tol_factor, tol_factor, tol_factor * 2
        ),
        call. = FALSE
      )
    }
  }

  # ---- 1. Compute novelty score -----------------------------------------------
  if (is.null(score)) {
    k_eff <- min(k, nrow(codebook))
    score <- .somalign_knn_mean_distance(scaled_data, codebook, k_eff, chunk_size)
  }

  # ---- 2. Identify tail cells -------------------------------------------------
  score_threshold <- stats::quantile(score, tail_quantile)
  tail_flag       <- score >= score_threshold   # per-cell logical

  # Helper to build the params list (avoids repetition across return paths)
  .make_params <- function() {
    list(
      k                = k,
      tail_quantile    = tail_quantile,
      n_clusters       = n_clusters,
      min_cluster      = min_cluster,
      min_batches      = min_batches,
      tol_factor       = tol_factor,
      merge_tol_factor = merge_tol_factor
    )
  }

  # ---- 3. Per-group k-means on tail cells ------------------------------------
  groups       <- unique(group)
  n_groups     <- length(groups)
  all_centroids <- list()   # will hold: centroid coords + group + size

  for (g in groups) {
    in_group  <- which(group == g)
    tail_in_g <- in_group[tail_flag[in_group]]
    n_tail_g  <- length(tail_in_g)

    if (n_tail_g < 1L) next

    coords_g    <- scaled_data[tail_in_g, , drop = FALSE]
    n_distinct_g <- nrow(unique(coords_g))

    # Adaptive k: cap so each cluster can plausibly have min_cluster cells,
    # and cap at number of distinct coordinate rows to avoid k-means crash.
    k_g <- max(1L, min(n_clusters, n_distinct_g, n_tail_g %/% min_cluster))

    if (k_g < 1L) next

    # k-means (set.seed guard: caller must use withr::local_seed in tests)
    km <- tryCatch(
      stats::kmeans(coords_g, centers = k_g, nstart = 3L, iter.max = 50L),
      error = function(e) NULL
    )
    if (is.null(km)) next

    centers  <- km$centers          # [k_g x n_features]
    cl_sizes <- tabulate(km$cluster, nbins = k_g)

    keep <- which(cl_sizes >= min_cluster)
    if (length(keep) == 0L) next

    for (ci in keep) {
      all_centroids <- c(all_centroids, list(list(
        coords = centers[ci, , drop = FALSE],
        group  = g,
        size   = cl_sizes[ci]
      )))
    }
  }

  # ---- 4. Cross-group matching ------------------------------------------------
  # Reference node spacing (median nearest-neighbour distance among codebook nodes)
  node_spacing <- .somalign_default_bandwidth(codebook)
  link_tol     <- tol_factor * node_spacing

  n_cent <- length(all_centroids)

  # Empty prototype matrix with correct column names — returned when 0 candidates
  empty_prototypes <- matrix(
    numeric(0L), nrow = 0L, ncol = n_features,
    dimnames = list(NULL, features)
  )

  if (n_cent == 0L) {
    result <- structure(
      list(
        prototypes       = empty_prototypes,
        n_groups_support = integer(0L),
        size             = integer(0L),
        score            = score,
        tail             = tail_flag,
        n_groups         = n_groups,
        params           = .make_params()
      ),
      class = "somalign_novelty_candidates"
    )
    return(result)
  }

  # Build centroid-coord matrix for pairwise distance
  cent_mat <- do.call(rbind, lapply(all_centroids, function(x) x$coords))
  rownames(cent_mat) <- NULL

  # Pairwise Euclidean SQUARED distances between centroids
  d2_cent <- .somalign_pairwise_distance(cent_mat, cent_mat)
  diag(d2_cent) <- Inf   # exclude self-links

  # Build adjacency: only link centroids from DIFFERENT groups
  cent_groups <- vapply(all_centroids, function(x) as.character(x$group), character(1L))

  # Union-find connected components
  parent <- seq_len(n_cent)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  union <- function(i, j) {
    ri <- find(i); rj <- find(j)
    if (ri != rj) parent[ri] <<- rj
  }

  for (i in seq_len(n_cent - 1L)) {
    for (j in seq(i + 1L, n_cent)) {
      if (cent_groups[i] == cent_groups[j]) next   # same group — skip
      if (d2_cent[i, j] <= link_tol^2) {
        union(i, j)
      }
    }
  }

  # Collect components — track group SETS (not just counts) for dedup union logic
  comp_ids <- vapply(seq_len(n_cent), find, integer(1L))
  unique_comps <- unique(comp_ids)

  prototypes_list      <- list()
  n_groups_support_vec <- integer(0L)
  size_vec             <- integer(0L)
  groups_sets_list     <- list()   # character vectors: the groups each candidate covers

  for (cid in unique_comps) {
    members <- which(comp_ids == cid)
    member_groups <- unique(cent_groups[members])
    if (length(member_groups) < min_batches) next

    # Mass-weighted mean centroid
    sizes_m  <- vapply(members, function(idx) all_centroids[[idx]]$size, integer(1L))
    coords_m <- do.call(rbind, lapply(members, function(idx) all_centroids[[idx]]$coords))
    rownames(coords_m) <- NULL

    weights   <- sizes_m / sum(sizes_m)
    prototype <- colSums(coords_m * weights)

    prototypes_list      <- c(prototypes_list, list(prototype))
    n_groups_support_vec <- c(n_groups_support_vec, length(member_groups))
    size_vec             <- c(size_vec, sum(sizes_m))
    groups_sets_list     <- c(groups_sets_list, list(member_groups))
  }

  n_candidates <- length(prototypes_list)

  if (n_candidates == 0L) {
    result <- structure(
      list(
        prototypes       = empty_prototypes,
        n_groups_support = integer(0L),
        size             = integer(0L),
        score            = score,
        tail             = tail_flag,
        n_groups         = n_groups,
        params           = .make_params()
      ),
      class = "somalign_novelty_candidates"
    )
    return(result)
  }

  prototypes_mat <- do.call(rbind, prototypes_list)
  colnames(prototypes_mat) <- features
  rownames(prototypes_mat) <- NULL

  # ---- 5. Post-mint prototype deduplication -----------------------------------
  # Agglomeratively merge minted prototypes that are near-duplicates (i.e., the
  # same biological population fragmented into several components by k-means or
  # cross-batch matching).  Uses single-linkage at merge_tol_factor*node_spacing.
  # n_groups_support for a merged candidate is |UNION of member group sets|, NOT
  # the sum of individual counts (a group supporting two merged members counts once).
  if (dedup_enabled && n_candidates >= 2L) {
    merge_tol <- merge_tol_factor * node_spacing

    # stats::dist() computes TRUE Euclidean distances (not squared), so the
    # threshold is merge_tol directly.
    proto_dist <- stats::dist(prototypes_mat, method = "euclidean")
    hc         <- stats::hclust(proto_dist, method = "single")
    merge_ids  <- stats::cutree(hc, h = merge_tol)

    if (length(unique(merge_ids)) < n_candidates) {
      # At least one merge happened — rebuild outputs
      merged_protos    <- list()
      merged_n_support <- integer(0L)
      merged_size      <- integer(0L)

      for (mid in unique(merge_ids)) {
        mem <- which(merge_ids == mid)

        # Size-weighted mean of member prototypes
        sizes_m    <- size_vec[mem]
        weights_m  <- sizes_m / sum(sizes_m)
        coord_m    <- prototypes_mat[mem, , drop = FALSE]
        new_proto  <- colSums(coord_m * weights_m)

        # Union of group sets across merged members
        union_groups <- unique(unlist(groups_sets_list[mem]))
        new_n_support <- length(union_groups)

        merged_protos    <- c(merged_protos, list(new_proto))
        merged_n_support <- c(merged_n_support, new_n_support)
        merged_size      <- c(merged_size, sum(sizes_m))
      }

      # Rebuild outputs from merged results
      prototypes_mat       <- do.call(rbind, merged_protos)
      colnames(prototypes_mat) <- features
      rownames(prototypes_mat) <- NULL
      n_groups_support_vec  <- merged_n_support
      size_vec              <- merged_size
    }
  }

  structure(
    list(
      prototypes       = prototypes_mat,
      n_groups_support = n_groups_support_vec,
      size             = size_vec,
      score            = score,
      tail             = tail_flag,
      n_groups         = n_groups,
      params           = .make_params()
    ),
    class = "somalign_novelty_candidates"
  )
}

# Internal: chunked k-NN mean Euclidean distance from x rows to codebook
.somalign_knn_mean_distance <- function(x, codebook, k, chunk_size) {
  n        <- nrow(x)
  k_eff    <- min(k, nrow(codebook))
  score    <- numeric(n)
  starts   <- seq(1L, n, by = chunk_size)

  for (s in starts) {
    e       <- min(s + chunk_size - 1L, n)
    idx     <- s:e
    x_chunk <- x[idx, , drop = FALSE]
    d2      <- .somalign_pairwise_distance(x_chunk, codebook)
    sel     <- .somalign_knn_select(d2, k_eff)
    score[idx] <- rowMeans(sqrt(sel$sq_dist))
  }
  score
}
