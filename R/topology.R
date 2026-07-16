utils::globalVariables(c("persistence"))

# H0 persistent homology via single-linkage (Kruskal's algorithm): sort all
# M(M-1)/2 edges ascending by distance, union-find merges, record each
# merge's distance as a "death" (birth is always 0 for H0). Returns a data
# frame of (birth, death, persistence) for the M-1 finite classes; the one
# essential (never-dying) class is not included as a row -- it is implicitly
# the "+1" in .somalign_h0_n_components(). `D` must be a genuine (non-squared)
# Euclidean distance matrix: persistent homology requires a metric space, and
# squared distances violate the triangle inequality.
#
# find_root() is a nested closure that only READS `parent` from the
# enclosing frame (ordinary lexical scoping, no `<<-` needed); the actual
# union (parent[rj] <- ri) happens directly in this function's own loop, not
# inside a helper -- this avoids relying on `<<-` for correctness at all.
.somalign_h0_persistence <- function(D) {
  m <- nrow(D)
  if (m <= 1L) {
    return(data.frame(birth = numeric(0), death = numeric(0), persistence = numeric(0)))
  }
  idx <- which(upper.tri(D), arr.ind = TRUE)
  dists <- D[idx]
  ord <- order(dists)
  idx <- idx[ord, , drop = FALSE]
  dists <- dists[ord]

  parent <- seq_len(m)
  find_root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  deaths <- numeric(m - 1L)
  k <- 0L
  for (e in seq_len(nrow(idx))) {
    ri <- find_root(idx[e, 1L])
    rj <- find_root(idx[e, 2L])
    if (ri != rj) {
      parent[rj] <- ri
      k <- k + 1L
      deaths[k] <- dists[e]
      if (k == m - 1L) break
    }
  }
  data.frame(birth = 0, death = deaths[seq_len(k)], persistence = deaths[seq_len(k)])
}

# Per-node connected-component membership at a distance threshold, via the
# same single-linkage / Kruskal union-find as .somalign_h0_persistence():
# nodes joined through edges of length <= threshold share a label. Exposed so
# callers needing membership (e.g. mass-weighted component sizes) reuse this
# union-find instead of running a second, independent clustering. The number
# of distinct labels equals .somalign_h0_n_components() at the same threshold
# (both count MST edges of length > threshold, plus one), so the two never
# disagree at the boundary. `D` must be a genuine (non-squared) metric.
.somalign_h0_components <- function(D, threshold) {
  m <- nrow(D)
  if (m <= 1L) return(rep(1L, m))
  idx <- which(upper.tri(D), arr.ind = TRUE)
  keep <- D[idx] <= threshold
  idx <- idx[keep, , drop = FALSE]

  parent <- seq_len(m)
  find_root <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  for (e in seq_len(nrow(idx))) {
    ri <- find_root(idx[e, 1L])
    rj <- find_root(idx[e, 2L])
    if (ri != rj) parent[rj] <- ri
  }
  roots <- vapply(seq_len(m), find_root, integer(1))
  match(roots, unique(roots))
}

# Number of H0 components with persistence > threshold, plus the one
# essential (always-surviving) component.
.somalign_h0_n_components <- function(diagram, threshold, n_nodes) {
  if (n_nodes <= 0L) return(0L)
  if (nrow(diagram) == 0L) return(1L)
  sum(diagram$persistence > threshold) + 1L
}

# Auto-selects a topology threshold from the reference's own within-node
# distance spread: the median across nodes of the 95th-percentile per-node
# distance. Falls back to 0.5 if distance_quantiles is absent or degenerate.
.somalign_topo_threshold <- function(reference) {
  dq <- reference$distance_quantiles
  if (is.null(dq)) return(0.5)
  q95 <- tryCatch(
    if (is.matrix(dq)) dq[, "95%"] else dq[["95%"]],
    error = function(e) NA_real_
  )
  q95 <- q95[is.finite(q95)]
  if (length(q95) == 0L) return(0.5)
  stats::median(q95)
}

.somalign_tda_diagram <- function(cb, max_dim = 1L) {
  if (!requireNamespace("TDA", quietly = TRUE)) return(NULL)
  TDA::ripsDiag(cb, maxdimension = max_dim, maxscale = Inf,
               dist = "euclidean", library = "GUDHI")
}

.somalign_tda_bottleneck <- function(d1, d2) {
  if (is.null(d1) || is.null(d2) || !requireNamespace("TDA", quietly = TRUE)) return(NULL)
  h0_pairs <- function(d) {
    diag <- d$diagram
    diag[diag[, "dimension"] == 0, c("Birth", "Death"), drop = FALSE]
  }
  tryCatch(TDA::bottleneck(h0_pairs(d1), h0_pairs(d2)), error = function(e) NULL)
}

.somalign_topology_audit_impl <- function(fit, threshold, use_tda, nodes) {
  allowed <- attr(fit$node_shifts, "correction_allowed")
  sel <- if (identical(nodes, "correction_allowed")) which(allowed) else seq_len(nrow(fit$query$codebook))

  cb_query <- fit$query$codebook[sel, , drop = FALSE]
  cb_corrected <- cb_query + fit$node_shifts[sel, , drop = FALSE]
  cb_ref <- fit$reference$codebook

  thresh <- if (is.null(threshold)) .somalign_topo_threshold(fit$reference) else threshold
  thresh_source <- if (is.null(threshold)) "auto" else "user"

  # sqrt() is required: .somalign_pairwise_distance() returns SQUARED
  # Euclidean distance (the OT cost convention, since F2), but persistent
  # homology requires a genuine metric (triangle inequality).
  d_query <- sqrt(.somalign_pairwise_distance(cb_query, cb_query))
  d_corr <- sqrt(.somalign_pairwise_distance(cb_corrected, cb_corrected))
  d_ref <- sqrt(.somalign_pairwise_distance(cb_ref, cb_ref))

  pd_q <- .somalign_h0_persistence(d_query)
  pd_c <- .somalign_h0_persistence(d_corr)
  pd_r <- .somalign_h0_persistence(d_ref)

  nq <- .somalign_h0_n_components(pd_q, thresh, nrow(cb_query))
  nc <- .somalign_h0_n_components(pd_c, thresh, nrow(cb_corrected))
  nr <- .somalign_h0_n_components(pd_r, thresh, nrow(cb_ref))

  tda <- .somalign_topology_tda_slots(use_tda, cb_query, cb_corrected, cb_ref)
  warn <- (nc - nq) != 0L
  if (warn) .somalign_topology_warn(nc, nq)

  structure(
    list(threshold = thresh, threshold_source = thresh_source,
         n_components_query = nq, n_components_corrected = nc,
         n_components_reference = nr, topology_delta = nc - nq,
         diagram_query = pd_q, diagram_corrected = pd_c, diagram_reference = pd_r,
         bottleneck_h0 = tda$bottleneck_h0,
         tda_query = tda$tda_query, tda_corrected = tda$tda_corrected,
         tda_reference = tda$tda_reference, topology_warning = warn),
    class = "somalign_topology"
  )
}

.somalign_topology_tda_slots <- function(use_tda, cb_query, cb_corrected, cb_ref) {
  if (!isTRUE(use_tda)) {
    return(list(tda_query = NULL, tda_corrected = NULL, tda_reference = NULL, bottleneck_h0 = NULL))
  }
  if (!requireNamespace("TDA", quietly = TRUE)) {
    message("TDA package not available; falling back to base-R H0 only.")
    return(list(tda_query = NULL, tda_corrected = NULL, tda_reference = NULL, bottleneck_h0 = NULL))
  }
  tda_q <- .somalign_tda_diagram(cb_query)
  tda_c <- .somalign_tda_diagram(cb_corrected)
  tda_r <- .somalign_tda_diagram(cb_ref)
  list(tda_query = tda_q, tda_corrected = tda_c, tda_reference = tda_r,
       bottleneck_h0 = .somalign_tda_bottleneck(tda_c, tda_r))
}

.somalign_topology_warn <- function(nc, nq) {
  direction <- if (nc < nq) "merged/erased" else "split"
  warning(sprintf(
    "topology_warning: corrected codebook has %d H0 component(s) vs %d in query ", nc, nq),
    sprintf("(delta = %+d; populations may have been %s). ", nc - nq, direction),
    "Inspect fit$diagnostics$topology for details.",
    call. = FALSE)
}

#' Compute a persistent-homology topology audit for a somalign fit
#'
#' Computes the H0 (connected-component) persistence diagram of the query,
#' corrected-query, and reference codebooks and reports how many robustly
#' separated populations survive at a given distance threshold. Batch
#' correction that merges or erases a population shows up as a drop in the
#' number of H0 components between the query and corrected codebooks.
#'
#' @param fit A `somalign_fit` object.
#' @param threshold Numeric scalar in reference-scaled Euclidean units. H0
#'   components with persistence (death - birth) greater than this value are
#'   counted as robustly separated populations. `NULL` (default) derives the
#'   threshold from the median of `reference$distance_quantiles`'s 95th
#'   percentile column -- the reference's own within-population distance
#'   spread.
#' @param use_tda Logical. When `TRUE`, additionally compute H0 + H1 via
#'   `TDA::ripsDiag()` if the `TDA` package is installed (falls back silently
#'   to the base-R H0-only result, with a one-time message, when `TDA` is
#'   absent). Default `FALSE`.
#' @param nodes One of `"correction_allowed"` (default; recommended -- nodes
#'   with no correction contribute unchanged, potentially spurious topology)
#'   or `"all"`. Note that `somalign_epsilon_sweep(..., topology = TRUE)`
#'   reports its topology columns using `"all"` (the correction-allowed set is
#'   epsilon-dependent, so it cannot be held fixed across a sweep); pass
#'   `nodes = "all"` here to reproduce those numbers for a given fit.
#'
#' @return A list of class `somalign_topology`; see the source fields
#'   `threshold`, `threshold_source` (`"auto"`/`"user"`),
#'   `n_components_query`, `n_components_corrected`, `n_components_reference`,
#'   `topology_delta` (corrected minus query; negative means merging),
#'   `diagram_query`/`diagram_corrected`/`diagram_reference` (data frames of
#'   birth/death/persistence), `bottleneck_h0` and `tda_*` (`NULL` unless
#'   `use_tda = TRUE` and `TDA` is installed), and `topology_warning`
#'   (`TRUE` when `topology_delta != 0`; a warning is also emitted).
#'
#' @details
#' This is a pure diagnostic: it reads `fit` and returns a new object without
#' modifying it. It is not run automatically inside [somalign_fit()] (too
#' expensive to compute by default); call it directly or via
#' `somalign_diagnostics(fit, topology = TRUE)`.
#'
#' High-dimensional marker spaces (p > 20) are valid but nearest-neighbour
#' distances concentrate, which can make H0 components less visually
#' intuitive -- a known property of topological data analysis in high
#' dimensions, not a bug.
#'
#' @seealso [somalign_diagnostics()], [somalign_fit()]
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' fit <- somalign_fit(qry, ref)
#' somalign_topology_audit(fit)
#' @export
somalign_topology_audit <- function(fit, threshold = NULL, use_tda = FALSE,
                                    nodes = c("correction_allowed", "all")) {
  if (!inherits(fit, "somalign_fit"))
    stop("`fit` must be a somalign_fit object.", call. = FALSE)
  nodes <- match.arg(nodes)
  .somalign_check_flag(use_tda, "use_tda")
  if (!is.null(threshold)) .somalign_check_pos_scalar(threshold, "threshold")
  .somalign_topology_audit_impl(fit, threshold, use_tda, nodes)
}

#' Print a somalign_topology object
#'
#' @param x A `somalign_topology` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @export
print.somalign_topology <- function(x, ...) {
  cat(
    "<somalign_topology>\n",
    sprintf("  threshold: %.4f (%s)\n", x$threshold, x$threshold_source),
    sprintf("  H0 components  query: %d  corrected: %d  reference: %d\n",
            x$n_components_query, x$n_components_corrected, x$n_components_reference),
    sprintf("  topology_delta: %+d   warning: %s\n", x$topology_delta, x$topology_warning),
    sep = ""
  )
  invisible(x)
}
