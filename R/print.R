#' Print a somalign_reference object
#'
#' @param x A \code{somalign_reference} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_reference <- function(x, ...) {
  cat(
    "<somalign_reference>\n",
    "  features: ", length(x$features), "\n",
    "  reference nodes: ", nrow(x$codebook), "\n",
    "  labelled nodes: ", sum(rowSums(x$label_prob) > 0), "\n",
    sep = ""
  )
  invisible(x)
}

#' Print a somalign_query object
#'
#' @param x A \code{somalign_query} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_query <- function(x, ...) {
  cat(
    "<somalign_query>\n",
    "  samples: ", nrow(x$scaled_data), "\n",
    "  features: ", ncol(x$scaled_data), "\n",
    "  query nodes: ", nrow(x$codebook), "\n",
    sep = ""
  )
  invisible(x)
}

#' Print a somalign_fit object
#'
#' Leads with the label-transfer result -- the primary product of a fit --
#' followed by the transport/solver line. Use [summary.somalign_fit()] for the
#' full label breakdown.
#'
#' @param x A \code{somalign_fit} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_fit <- function(x, ...) {
  cat("<somalign_fit>\n")
  .somalign_print_label_headline(x)
  cat(
    "  solver: ", x$diagnostics$solver$used,
    "  |  query nodes: ", nrow(x$query$codebook),
    "  |  reference nodes: ", nrow(x$reference$codebook),
    "  |  transport mass: ", signif(x$diagnostics$ot$transport_mass, 4), "\n",
    sep = ""
  )
  invisible(x)
}

# Shared one-line label-transfer headline for the fit print methods.
.somalign_print_label_headline <- function(x) {
  s <- .somalign_label_summary(x)
  if (!isTRUE(s$enabled)) {
    cat("  label transfer: disabled (reference has no labels)\n")
    return(invisible(NULL))
  }
  cat(sprintf(
    "  label transfer: %.1f%% of cells accepted across %d class(es); median confidence %.2f, median margin %.2f\n",
    100 * s$accepted_fraction, s$n_classes, s$median_confidence, s$median_margin
  ))
  invisible(NULL)
}

#' Summarise a somalign_fit's label transfer
#'
#' @param object A \code{somalign_fit} object.
#' @param ... Ignored.
#'
#' @return \code{object}, invisibly. Prints the accepted-cell fraction, class
#'   distribution, and confidence/margin summary.
#' @export
summary.somalign_fit <- function(object, ...) {
  cat("<somalign_fit> label-transfer summary\n")
  s <- .somalign_label_summary(object)
  if (!isTRUE(s$enabled)) {
    cat("  label transfer disabled (reference has no labels).\n")
    return(invisible(object))
  }
  cat(sprintf("  cells: %d  |  accepted: %d (%.1f%%)  |  classes: %d\n",
              s$n_cells, s$n_accepted, 100 * s$accepted_fraction, s$n_classes))
  cat(sprintf("  confidence quartiles (accepted): %.2f / %.2f / %.2f\n",
              s$confidence_quartiles[[1]], s$confidence_quartiles[[2]],
              s$confidence_quartiles[[3]]))
  cat(sprintf("  median margin (accepted): %.2f\n", s$median_margin))
  cat("  accepted class distribution:\n")
  cd <- s$class_distribution
  for (i in seq_along(cd)) {
    cat(sprintf("    %-20s %d\n", names(cd)[i], cd[[i]]))
  }
  invisible(object)
}

#' Print a somalign_anchored_fit object
#'
#' @param x A \code{somalign_anchored_fit} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_anchored_fit <- function(x, ...) {
  cat("<somalign_anchored_fit>\n")
  .somalign_print_label_headline(x)
  cat(
    "  solver: ", x$diagnostics$solver$used,
    "  |  query nodes: ", nrow(x$query$codebook),
    "  |  reference nodes: ", nrow(x$reference$codebook),
    "  |  transport mass: ", signif(x$diagnostics$ot$transport_mass, 4), "\n",
    "  anchors: ", x$anchors$n_anchors,
    " (", round(100 * x$anchors$coverage_fraction, 1), "% node coverage) ",
    "-- correction is a diagnostic, not a corrected-expression product\n",
    sep = ""
  )
  invisible(x)
}

#' Print a somalign_exclusion_test object
#'
#' @param x A \code{somalign_exclusion_test} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_exclusion_test <- function(x, ...) {
  cat(
    "<somalign_exclusion_test>\n",
    sprintf("  sv_observed = %.4f   p = %.4f   verdict = %s\n",
            x$sv_observed, x$p_value, x$verdict),
    sprintf("  null sv (2.5%%/50%%/97.5%%): %.4f / %.4f / %.4f\n",
            x$null_quantiles[["2.5%"]], x$null_quantiles[["50%"]],
            x$null_quantiles[["97.5%"]]),
    sprintf("  rank = %d   n_anchors = %d   n_features = %d\n",
            x$rank_used, x$n_anchors, x$n_features),
    sep = ""
  )
  invisible(x)
}

#' @method print somalign_soft_labels
#' @export
print.somalign_soft_labels <- function(x, ...) {
  cat(sprintf(
    "<somalign_soft_labels> [%d cells x %d labels]  k = %s  bandwidth = %.4g\n",
    nrow(x), ncol(x), format(attr(x, "k")), attr(x, "bandwidth")))
  invisible(x)
}

#' @method print somalign_soft_frequencies
#' @export
print.somalign_soft_frequencies <- function(x, ...) {
  cat(sprintf(
    "<somalign_soft_frequencies> [%d groups x %d labels]  %s  k = %s  bandwidth = %.4g\n",
    nrow(x), ncol(x),
    if (isTRUE(attr(x, "normalized"))) "frequencies" else "soft counts",
    format(attr(x, "k")), attr(x, "bandwidth")))
  invisible(x)
}

#' @method print somalign_corrected_expression
#' @export
print.somalign_corrected_expression <- function(x, ...) {
  cat(sprintf(
    "<somalign_corrected_expression> [%d cells x %d markers]  units = %s  smooth = %s  k = %s  bandwidth = %s\n",
    nrow(x), ncol(x), attr(x, "units"), attr(x, "smooth"),
    format(attr(x, "k")),
    if (is.na(attr(x, "bandwidth"))) "NA" else sprintf("%.4g", attr(x, "bandwidth"))))
  invisible(x)
}

#' Print a somalign_novelty_candidates object
#'
#' @param x A \code{somalign_novelty_candidates} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_novelty_candidates <- function(x, ...) {
  n_cand      <- nrow(x$prototypes)
  tail_frac   <- mean(x$tail)
  cat(sprintf(
    "<somalign_novelty_candidates>\n  candidates minted: %d  |  groups: %d  |  tail fraction: %.1f%%\n",
    n_cand, x$n_groups, 100 * tail_frac
  ))
  if (n_cand > 0L) {
    for (i in seq_len(n_cand)) {
      cat(sprintf(
        "  [%d] groups_support = %d  |  tail_cells = %d\n",
        i, x$n_groups_support[i], x$size[i]
      ))
    }
  }
  invisible(x)
}
