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
#' @param x A \code{somalign_fit} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_fit <- function(x, ...) {
  cat(
    "<somalign_fit>\n",
    "  solver: ", x$diagnostics$solver$used, "\n",
    "  query nodes: ", nrow(x$query$codebook), "\n",
    "  reference nodes: ", nrow(x$reference$codebook), "\n",
    "  transport mass: ", signif(x$diagnostics$ot$transport_mass, 4), "\n",
    sep = ""
  )
  invisible(x)
}

#' Print a somalign_anchored_fit object
#'
#' @param x A \code{somalign_anchored_fit} object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.somalign_anchored_fit <- function(x, ...) {
  cat(
    "<somalign_anchored_fit>\n",
    "  solver: ", x$diagnostics$solver$used, "\n",
    "  query nodes: ", nrow(x$query$codebook), "\n",
    "  reference nodes: ", nrow(x$reference$codebook), "\n",
    "  transport mass: ", signif(x$diagnostics$ot$transport_mass, 4), "\n",
    "  anchors: ", x$anchors$n_anchors,
    " (", round(100 * x$anchors$coverage_fraction, 1), "% node coverage)\n",
    sep = ""
  )
  invisible(x)
}
