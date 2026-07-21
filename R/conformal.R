#' Conformal prediction sets for label transfer
#'
#' Turns per-cell class probabilities into **conformal prediction sets** with a
#' distribution-free, finite-sample coverage guarantee. Given a labelled
#' calibration set, `somalign_conformal_labels()` calibrates a single threshold
#' (split conformal) so that, for a fresh exchangeable query cell, the true label
#' is contained in the returned set with probability at least `1 - alpha`. A cell
#' whose set is empty is an **abstention**; a set of size one is a confident
#' assignment; a larger set flags genuine ambiguity between labels.
#'
#' This complements the heuristic acceptance gate of [somalign_results()]: rather
#' than thresholding a single confidence, it returns calibrated label *sets* with
#' an explicit coverage level. The probability matrices are typically the per-cell
#' soft label memberships from [somalign_soft_labels()], but any cells-by-classes
#' matrix of non-negative scores works.
#'
#' The nonconformity score is `1 - p(true class)`; the set for a query cell is
#' every class whose probability is at least `1 - q`, where `q` is the
#' `ceiling((n + 1)(1 - alpha)) / n` empirical quantile of the calibration scores
#' (Angelopoulos & Bates, 2023, \doi{10.1561/2200000101}). With
#' `class_conditional = TRUE` the threshold is computed separately per true class
#' (Mondrian conformal), giving class-conditional rather than only marginal
#' coverage, which is preferable under class imbalance.
#'
#' @section Warning -- exchangeability and batch effects:
#' The coverage guarantee holds only when the calibration cells and the query
#' cells are **exchangeable** (drawn from the same distribution). This is exactly
#' what a batch effect breaks: if you calibrate on a labelled reference and apply
#' the threshold to a new batch, the guarantee is no longer exact, and the true
#' miscoverage can exceed `alpha` in a way the procedure cannot detect. Calibrate
#' on data as close as possible to the query distribution (for example a held-out,
#' same-batch labelled split, or reference cells after the query has been aligned
#' into reference-scaled space), and treat the coverage level as approximate under
#' residual shift. When a labelled same-distribution calibration set is
#' unavailable, `somalign_cross_validate()` provides one by construction.
#'
#' @param prob_query Numeric matrix, query cells by classes, of per-cell class
#'   probabilities (columns named by class). Rows need not sum to one.
#' @param prob_calibration Numeric matrix, calibration cells by classes, with the
#'   same class columns as `prob_query`.
#' @param truth_calibration Character/factor vector of true classes for the
#'   calibration cells; every value must be a column of `prob_calibration`.
#' @param alpha Target miscoverage in (0, 1); coverage is at least `1 - alpha`.
#'   Default `0.1`.
#' @param class_conditional Logical; if `TRUE`, calibrate a per-class threshold
#'   (Mondrian conformal) for class-conditional coverage. Default `FALSE`.
#'
#' @return An object of class `somalign_conformal`: a list with `sets` (a
#'   query-by-class logical matrix; `TRUE` = class is in the cell's set),
#'   `set_size` (per-cell integer), `abstain` (per-cell logical, `set_size == 0`),
#'   `threshold` (the calibrated `q`, scalar or per-class), `alpha`, `classes`,
#'   and `class_conditional`.
#' @examples
#' set.seed(1)
#' classes <- c("A", "B", "C")
#' # calibrated toy probabilities: true class gets the most mass on average
#' gen <- function(n) {
#'   y <- sample(classes, n, replace = TRUE)
#'   p <- matrix(stats::runif(n * 3, 0, 0.4), n, 3, dimnames = list(NULL, classes))
#'   p[cbind(seq_len(n), match(y, classes))] <-
#'     p[cbind(seq_len(n), match(y, classes))] + 0.6
#'   list(p = p / rowSums(p), y = y)
#' }
#' cal <- gen(500); qry <- gen(200)
#' cp <- somalign_conformal_labels(qry$p, cal$p, cal$y, alpha = 0.1)
#' cp
#' @export
somalign_conformal_labels <- function(prob_query, prob_calibration,
                                      truth_calibration, alpha = 0.1,
                                      class_conditional = FALSE) {
  prob_query <- as.matrix(prob_query)
  prob_calibration <- as.matrix(prob_calibration)
  if (!is.numeric(prob_query) || !is.numeric(prob_calibration))
    stop("`prob_query` and `prob_calibration` must be numeric matrices.", call. = FALSE)
  if (is.null(colnames(prob_calibration)) || is.null(colnames(prob_query)))
    stop("Both probability matrices must have class names as column names.", call. = FALSE)
  if (!identical(colnames(prob_query), colnames(prob_calibration)))
    stop("`prob_query` and `prob_calibration` must have the same class columns in the same order.",
         call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  if (any(!is.finite(prob_query)) || any(!is.finite(prob_calibration)))
    stop("`prob_query` and `prob_calibration` must contain only finite values.",
         call. = FALSE)
  classes <- colnames(prob_calibration)
  truth_calibration <- as.character(truth_calibration)
  if (length(truth_calibration) != nrow(prob_calibration))
    stop("`truth_calibration` must have one value per calibration cell.", call. = FALSE)
  if (!all(truth_calibration %in% classes))
    stop("Every value in `truth_calibration` must be a class column of `prob_calibration`.",
         call. = FALSE)

  # Nonconformity score: 1 - probability of the true class.
  n <- nrow(prob_calibration)
  true_idx <- match(truth_calibration, classes)
  scores <- 1 - prob_calibration[cbind(seq_len(n), true_idx)]

  # Split-conformal quantile: the ceiling((m + 1)(1 - alpha)) / m order statistic.
  # When the calibration set is too small for the requested alpha, the level is
  # not attainable, so the threshold is +Inf (every class admitted) -- valid, if
  # uninformative.
  conformal_q <- function(s) {
    m <- length(s)
    if (m == 0L) return(Inf)
    k <- ceiling((m + 1) * (1 - alpha))
    if (k > m) return(Inf)
    sort(s)[k]
  }

  if (isTRUE(class_conditional)) {
    # A class with fewer than ceiling(1/alpha) calibration cells cannot attain the
    # level, so its threshold is unconstrained (Inf) and it is admitted to *every*
    # cell's set, inflating set sizes. Warn so this is not silent.
    class_n <- tabulate(match(truth_calibration, classes), nbins = length(classes))
    thin <- classes[class_n < ceiling(1 / alpha)]
    if (length(thin))
      warning(sprintf(
        paste0("class-conditional conformal: class(es) %s have fewer than ",
               "ceiling(1/alpha) = %d calibration cells; their threshold is ",
               "unconstrained and admits every cell."),
        paste(thin, collapse = ", "), ceiling(1 / alpha)), call. = FALSE)
    q <- vapply(classes, function(cl) conformal_q(scores[truth_calibration == cl]),
                numeric(1))
    # Admit class c when p(c) >= 1 - q_c.
    keep <- sweep(prob_query, 2, 1 - q, FUN = ">=")
  } else {
    q <- conformal_q(scores)
    keep <- prob_query >= (1 - q)
  }
  storage.mode(keep) <- "logical"
  dimnames(keep) <- list(rownames(prob_query), classes)

  structure(
    list(
      sets = keep,
      set_size = as.integer(rowSums(keep)),
      abstain = rowSums(keep) == 0L,
      threshold = q,
      alpha = alpha,
      classes = classes,
      class_conditional = isTRUE(class_conditional)
    ),
    class = "somalign_conformal"
  )
}

#' Print a somalign_conformal object
#'
#' @param x A `somalign_conformal` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @method print somalign_conformal
#' @export
print.somalign_conformal <- function(x, ...) {
  n <- nrow(x$sets)
  sz <- x$set_size
  cat("<somalign_conformal>\n")
  cat(sprintf("  cells: %d  |  classes: %d  |  target coverage: %.0f%%%s\n",
              n, length(x$classes), 100 * (1 - x$alpha),
              if (x$class_conditional) "  (class-conditional)" else ""))
  cat(sprintf("  abstentions (empty set): %d (%.1f%%)\n",
              sum(x$abstain), 100 * mean(x$abstain)))
  cat(sprintf("  confident singletons:    %d (%.1f%%)\n",
              sum(sz == 1L), 100 * mean(sz == 1L)))
  cat(sprintf("  ambiguous (>=2 labels):  %d (%.1f%%)\n",
              sum(sz >= 2L), 100 * mean(sz >= 2L)))
  cat(sprintf("  mean set size: %.2f\n", mean(sz)))
  invisible(x)
}
