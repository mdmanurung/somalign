#' Subspace sensitivity analysis for anchored batch correction
#'
#' Bootstraps the anchor displacement matrix `D` to quantify how stable the
#' estimated batch subspace `V` is, and propagates that uncertainty to
#' per-node correction confidence intervals and "tipping angles" -- the
#' smallest rotation of `V` that would erase a node's correction.
#'
#' @param fit A `somalign_anchored_fit` with `correction = "subspace"` or
#'   `"both"`.
#' @param n_boot Positive integer. Number of bootstrap replicates of `D`.
#'   Default `200L`.
#' @param variance_threshold Numeric in (0, 1], or `NULL`. Variance threshold
#'   for SVD rank selection in bootstrap replicates. `NULL` (default) reuses
#'   the threshold from the original fit (`fit$anchors$variance_threshold`).
#' @param conf_level Numeric in (0, 1). Confidence level for the node-shift
#'   confidence intervals. Default `0.95`.
#' @param seed Integer or `NULL`. RNG seed for reproducibility (restored on
#'   exit; does not leak into the caller's session). Default `1L`; `NULL`
#'   disables seeding.
#'
#' @return A list of class `somalign_subspace_sensitivity`:
#' \describe{
#'   \item{`node_correction_ci`}{M x 2 matrix (lower, upper) -- bootstrap CI
#'     on each node's corrected-shift norm.}
#'   \item{`node_shift_ci`}{M x p x 2 array -- per-feature bootstrap CI.}
#'   \item{`subspace_angles`}{n_boot x rank matrix of principal angles
#'     (degrees) between the fitted `V` and each bootstrap `V_b`.}
#'   \item{`tipping_angle_deg`}{Length-M vector. For rank-1 `V`, the analytic
#'     angle \eqn{\arcsin(|\hat{s}_i \cdot v|)} between a node's unit raw
#'     shift and `v`; for rank > 1, a conservative proxy from the minimum
#'     singular value of the (1 x rank) projection. `NA` for disallowed or
#'     zero-norm nodes. Small (< 10 degrees) signals a fragile correction;
#'     large (> 45 degrees) indicates robustness.}
#'   \item{`anchor_leverage`}{Length-n_anchors vector: the maximum principal
#'     angle (degrees) between `V` and the leave-one-out subspace when that
#'     anchor is dropped -- a Cook's-distance analog for anchor influence.}
#'   \item{`n_boot`, `conf_level`, `subspace_rank`, `n_anchors`,
#'     `variance_threshold`}{Metadata.}
#' }
#'
#' @details
#' Raw pre-projection node shifts are recovered from stored fit components
#' (`fit$correspondence`, both codebooks), not from `fit$node_shifts`, which
#' stores *post-projection* shifts (`S \%*\% V \%*\% t(V)`) for subspace fits --
#' using it directly as the raw shift would double-project and produce
#' invalid tipping angles.
#'
#' Bootstrap replicate rank is fixed to the point-estimate rank of the
#' original fit (not re-selected per replicate), so principal angles compare
#' subspaces of the same dimension; if a replicate's own rank is lower,
#' its `V_b` is zero-padded.
#'
#' @seealso [somalign_fit_anchored()], [somalign_exclusion_test()]
#' @examples
#' set.seed(1)
#' p <- 3L
#' mat <- rbind(
#'   matrix(rnorm(20 * p, mean = -2), ncol = p),
#'   matrix(rnorm(20 * p, mean =  2), ncol = p)
#' )
#' colnames(mat) <- paste0("F", seq_len(p))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' shifted <- mat + 0.5
#' qry <- somalign_query(shifted, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                       rlen = 5)
#' anc_idx <- 1:10
#' fit <- somalign_fit_anchored(qry, ref,
#'                               anchor_old = mat[anc_idx, , drop = FALSE],
#'                               anchor_new = shifted[anc_idx, , drop = FALSE],
#'                               rho_anchor = 1, correction = "subspace")
#' somalign_subspace_sensitivity(fit, n_boot = 50L)
#' @export
somalign_subspace_sensitivity <- function(fit, n_boot = 200L,
                                          variance_threshold = NULL,
                                          conf_level = 0.95,
                                          seed = 1L) {
  .somalign_check_sensitivity_args(fit, n_boot, conf_level, seed)
  D <- fit$anchors$displacements
  V <- fit$anchors$batch_subspace$V
  r <- fit$anchors$batch_subspace$rank
  n_a <- nrow(D)
  vt <- if (is.null(variance_threshold)) fit$anchors$variance_threshold else variance_threshold
  if (n_a < r * 5L) {
    warning(sprintf(
      "n_anchors (%d) < rank * 5 (%d); the bootstrap distribution may be degenerate. ",
      n_a, r * 5L), "Consider adding more anchor pairs.", call. = FALSE)
  }

  s_raw <- .somalign_recover_raw_shifts(fit)
  allowed <- attr(fit$node_shifts, "correction_allowed")
  boot_res <- .somalign_seeded_bootstrap_subspace(D, V, r, n_boot, vt, seed)
  ci_res <- .somalign_node_correction_ci(s_raw, boot_res$V_boots, conf_level)

  structure(
    list(
      node_correction_ci = ci_res$correction_ci,
      node_shift_ci = ci_res$shift_ci,
      subspace_angles = boot_res$angle_mat,
      tipping_angle_deg = .somalign_tipping_angle(s_raw, V, allowed),
      anchor_leverage = .somalign_anchor_leverage(D, V, vt),
      n_boot = n_boot,
      conf_level = conf_level,
      subspace_rank = r,
      n_anchors = n_a,
      variance_threshold = vt
    ),
    class = "somalign_subspace_sensitivity"
  )
}

.somalign_check_sensitivity_args <- function(fit, n_boot, conf_level, seed) {
  if (!inherits(fit, "somalign_anchored_fit"))
    stop("`fit` must be a somalign_anchored_fit object.", call. = FALSE)
  if (!fit$anchors$correction %in% c("subspace", "both"))
    stop("`somalign_subspace_sensitivity` requires a fit with correction = 'subspace' or 'both'.",
         call. = FALSE)
  if (is.null(fit$anchors$batch_subspace) || is.null(fit$anchors$displacements))
    stop("`fit$anchors` is missing the batch subspace or displacement matrix. Refit with a ",
         "version of somalign that stores them.", call. = FALSE)
  .somalign_check_pos_int(n_boot, "n_boot")
  .somalign_check_prob_scalar(conf_level, "conf_level")
  if (!is.null(seed)) .somalign_check_pos_int(seed, "seed")
  invisible(NULL)
}

# Recovers pre-projection ("raw") node shifts from stored fit components.
# fit$node_shifts stores POST-projection shifts (S %*% V %*% t(V)) for
# subspace fits, so it cannot be used directly here -- doing so would
# double-project and invalidate every downstream tipping-angle/bootstrap
# calculation.
.somalign_recover_raw_shifts <- function(fit) {
  bary <- fit$correspondence %*% fit$reference$codebook
  s_raw <- bary - fit$query$codebook
  allowed <- attr(fit$node_shifts, "correction_allowed")
  s_raw[!allowed, ] <- 0
  s_raw
}

# Principal angles (degrees) between the column spaces of V1, V2 (both p x r).
# A bootstrap replicate whose D collapsed to (near-)zero variance yields an
# all-zero V2 column (see .somalign_subspace_svd's zero-D guard); such
# columns are reported at the maximal 90-degree angle rather than producing
# a meaningless SVD of an all-zero matrix.
.somalign_principal_angles <- function(V1, V2) {
  r <- ncol(V1)
  col_norms <- sqrt(colSums(V2^2))
  bad <- col_norms < .Machine$double.eps * 10
  angles <- rep(90, r)
  if (any(!bad)) {
    sv <- svd(crossprod(V1[, !bad, drop = FALSE], V2[, !bad, drop = FALSE]),
             nu = 0L, nv = 0L)$d
    angles[!bad] <- acos(pmin(pmax(sv, -1), 1)) * (180 / pi)
  }
  angles
}

# Bootstraps D (resampling anchor rows WITH replacement -- a valid
# nonparametric bootstrap, unlike the WITHOUT-replacement permutation used
# by somalign_exclusion_test()'s null, which would leave D's Gram matrix,
# and hence V, exactly unchanged), re-estimates V_b at the fixed
# point-estimate rank r for comparability, and records the principal angles
# to the original V.
.somalign_bootstrap_subspace <- function(D, V, r, n_boot, variance_threshold) {
  n_a <- nrow(D)
  p <- ncol(D)
  v_boots <- array(0, dim = c(p, r, n_boot))
  angle_mat <- matrix(NA_real_, nrow = n_boot, ncol = r)
  for (b in seq_len(n_boot)) {
    d_b <- D[sample.int(n_a, replace = TRUE), , drop = FALSE]
    sub_b <- .somalign_subspace_svd(d_b, variance_threshold)
    r_use <- min(r, sub_b$rank)
    v_b <- matrix(0, nrow = p, ncol = r)
    v_b[, seq_len(r_use)] <- sub_b$V[, seq_len(r_use), drop = FALSE]
    v_boots[, , b] <- v_b
    angle_mat[b, ] <- .somalign_principal_angles(V, v_b)
  }
  list(V_boots = v_boots, angle_mat = angle_mat)
}

# Runs the bootstrap with a seed local to this call: the caller's global RNG
# state is saved before seeding and restored on exit.
.somalign_seeded_bootstrap_subspace <- function(D, V, r, n_boot, variance_threshold, seed) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (!is.null(old_seed))
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  .somalign_bootstrap_subspace(D, V, r, n_boot, variance_threshold)
}

# Bootstrap confidence intervals on per-node corrected-shift norm and
# per-feature shift, from the raw shifts projected through each bootstrap V_b.
.somalign_node_correction_ci <- function(s_raw, v_boots, conf_level) {
  m <- nrow(s_raw)
  p <- ncol(s_raw)
  n_boot <- dim(v_boots)[3]
  alpha <- 1 - conf_level
  probs <- c(alpha / 2, 1 - alpha / 2)

  correction_boot <- matrix(NA_real_, nrow = n_boot, ncol = m)
  shift_boot <- array(NA_real_, dim = c(n_boot, m, p))
  for (b in seq_len(n_boot)) {
    v_b <- v_boots[, , b]
    corr_b <- s_raw %*% v_b %*% t(v_b)
    correction_boot[b, ] <- sqrt(rowSums(corr_b^2))
    shift_boot[b, , ] <- corr_b
  }

  correction_ci <- t(apply(correction_boot, 2, stats::quantile, probs = probs, na.rm = TRUE))
  colnames(correction_ci) <- c("lower", "upper")
  shift_ci <- array(NA_real_, dim = c(m, p, 2), dimnames = list(NULL, colnames(s_raw), c("lower", "upper")))
  for (j in seq_len(p)) {
    shift_ci[, j, ] <- t(apply(shift_boot[, , j], 2, stats::quantile, probs = probs, na.rm = TRUE))
  }
  list(correction_ci = correction_ci, shift_ci = shift_ci)
}

# Per-node tipping angle (degrees): the smallest rotation of V that would
# erase (rank-1: exactly reverse) node i's projected correction. Rank-1 has
# an analytic form; rank > 1 uses the minimum singular value of the (1 x r)
# projection as a conservative proxy. NA for disallowed or zero-norm nodes.
.somalign_tipping_angle <- function(s_raw, V, allowed) {
  m <- nrow(s_raw)
  r <- ncol(V)
  out <- rep(NA_real_, m)
  for (i in seq_len(m)) {
    if (!allowed[i]) next
    s <- s_raw[i, ]
    nrm <- sqrt(sum(s^2))
    if (nrm < .Machine$double.eps * 10) next
    s_hat <- s / nrm
    if (r == 1L) {
      out[i] <- asin(pmin(abs(sum(s_hat * V[, 1])), 1)) * (180 / pi)
    } else {
      sv_min <- min(svd(matrix(s_hat, nrow = 1) %*% V, nu = 0L, nv = 0L)$d)
      out[i] <- asin(pmin(sv_min, 1)) * (180 / pi)
    }
  }
  out
}

# Cook's-distance analog: leave-one-out influence of each anchor pair on V,
# measured as the maximum principal angle between the full-data V and the
# leave-one-out V_{-i}.
.somalign_anchor_leverage <- function(D, V, variance_threshold) {
  n_a <- nrow(D)
  lev <- numeric(n_a)
  for (i in seq_len(n_a)) {
    sub_loo <- .somalign_subspace_svd(D[-i, , drop = FALSE], variance_threshold)
    lev[i] <- max(.somalign_principal_angles(V, sub_loo$V))
  }
  lev
}

#' Print a somalign_subspace_sensitivity object
#'
#' @param x A `somalign_subspace_sensitivity` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#' @export
print.somalign_subspace_sensitivity <- function(x, ...) {
  ta <- x$tipping_angle_deg[!is.na(x$tipping_angle_deg)]
  cat(
    "<somalign_subspace_sensitivity>\n",
    sprintf("  n_anchors = %d   rank = %d   n_boot = %d   conf_level = %.2f\n",
            x$n_anchors, x$subspace_rank, x$n_boot, x$conf_level),
    sprintf("  median principal angle: %.1f deg\n", stats::median(x$subspace_angles[, 1])),
    sprintf("  median tipping angle:    %.1f deg (n=%d allowed nodes)\n",
            if (length(ta)) stats::median(ta) else NA_real_, length(ta)),
    sep = ""
  )
  invisible(x)
}
