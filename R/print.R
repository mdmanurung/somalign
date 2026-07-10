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
