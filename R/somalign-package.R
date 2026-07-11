#' somalign: align query SOMs to fixed reference SOMs
#'
#' The package keeps direct projection into a fixed reference `kohonen` SOM as
#' the conservative primary result and returns auxiliary optimal-transport
#' corrected projections for visualisation and annotation.
#'
#' @keywords internal
#' @importFrom kohonen getCodes som somgrid
#' @importFrom stats median quantile sd
"_PACKAGE"
