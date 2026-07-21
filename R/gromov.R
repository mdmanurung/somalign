# Cross-panel SOM alignment via entropic Gromov-Wasserstein optimal transport.
#
# Gromov-Wasserstein (GW) aligns two point sets by matching their *intra-set*
# distance structures, so it needs no shared feature (marker) space -- unlike
# somalign_fit(), which requires the query and reference codebooks to live in the
# same reference-scaled coordinates. This makes GW a route to aligning a query
# SOM measured on a different panel/instrument to a fixed reference SOM.
#
# PROTOTYPE / SCOPE: this implements *balanced* entropic GW (all query mass is
# transported). That assumes the query and reference share roughly the same
# populations. Panels that add or drop populations need the *unbalanced* or
# *partial* GW relaxations used by Pamona (Cao et al. 2021) and SCOTv2 (Demetci
# et al. 2022); those are noted as the next step, not implemented here.

# Balanced entropic OT (log-domain Sinkhorn) for a fixed cost and marginals.
# Returns the transport plan. `cost` is shifted by its minimum for numerical
# stability; that is an additive constant and does not change the plan.
.somalign_gw_sinkhorn <- function(cost, p, q, epsilon, max_iter = 200L, tol = 1e-9) {
  cost <- cost - min(cost)
  logp <- log(p)
  logq <- log(q)
  Ke <- -cost / epsilon                 # log-kernel
  f <- numeric(length(p))
  g <- numeric(length(q))
  for (it in seq_len(max_iter)) {
    f_new <- epsilon * (logp - apply(sweep(Ke, 2, g / epsilon, "+"), 1,
                                     .somalign_logsumexp))
    g_new <- epsilon * (logq - apply(sweep(Ke, 1, f_new / epsilon, "+"), 2,
                                     .somalign_logsumexp))
    if (max(abs(f_new - f), abs(g_new - g)) < tol) {
      f <- f_new; g <- g_new; break
    }
    f <- f_new; g <- g_new
  }
  log_plan <- sweep(sweep(Ke, 1, f / epsilon, "+"), 2, g / epsilon, "+")
  plan <- exp(log_plan)
  plan[!is.finite(plan)] <- 0
  plan
}

# Project a near-feasible plan onto the exact transport polytope {P1 = p,
# P^T1 = q} (Altschuler, Weed & Rigollet, 2017). Alternating Sinkhorn leaves one
# marginal with a small residual at finite iterations / small epsilon; this
# rounding restores both marginals exactly at O(nm) cost, so the returned
# coupling is a valid balanced transport plan regardless of solver tolerance.
.somalign_round_transport <- function(P, p, q) {
  tiny <- .Machine$double.xmin
  P <- P * pmin(1, p / pmax(rowSums(P), tiny))        # scale each row down to <= p
  P <- sweep(P, 2, pmin(1, q / pmax(colSums(P), tiny)), "*")  # then columns to <= q
  err_r <- p - rowSums(P)
  err_c <- q - colSums(P)
  s <- sum(abs(err_r))
  if (s > 0) P <- P + outer(err_r, err_c) / s
  P
}

# Core entropic Gromov-Wasserstein. C1 (n x n) and C2 (m x m) are intra-set
# distance matrices; p, q are marginals. The GW coupling iteration (Peyre et al.
# 2016) reduces -- for the coupling -- to a Sinkhorn step on the pseudo-cost
# -2 * C1 %*% T %*% C2 (the additive row/column constants of the full GW tensor
# do not affect the plan and are dropped).
.somalign_gromov_wasserstein <- function(C1, C2, p, q, epsilon = 0.05,
                                         max_iter = 50L, tol = 1e-6,
                                         sinkhorn_max_iter = 200L) {
  # Scale-normalise so the pseudo-cost is O(1) and the Sinkhorn solve stays
  # well-conditioned for a fixed epsilon. This is NOT coupling-invariant in the
  # entropic problem: dividing C1, C2 by their maxima rescales the pseudo-cost by
  # 1 / (max(C1) max(C2)), equivalent to scaling the effective regularisation, so
  # the entropic plan differs from what the nominal `epsilon` alone would give.
  mx1 <- max(C1); mx2 <- max(C2)
  if (mx1 > 0) C1 <- C1 / mx1
  if (mx2 > 0) C2 <- C2 / mx2
  Tplan <- outer(p, q)                  # independent coupling as the starting point
  converged <- FALSE
  iter <- max_iter
  for (it in seq_len(max_iter)) {
    pseudo_cost <- -2 * (C1 %*% Tplan %*% C2)
    T_new <- .somalign_gw_sinkhorn(pseudo_cost, p, q, epsilon,
                                   max_iter = sinkhorn_max_iter)
    delta <- max(abs(T_new - Tplan))
    Tplan <- T_new
    if (delta < tol) {
      converged <- TRUE; iter <- it; break
    }
  }
  # Restore exact marginals on the returned coupling.
  Tplan <- .somalign_round_transport(Tplan, p, q)
  list(coupling = Tplan, iterations = iter, converged = converged)
}

#' Cross-panel SOM alignment by Gromov-Wasserstein optimal transport (prototype)
#'
#' Aligns a query SOM to a fixed reference SOM by matching their intra-codebook
#' distance structures with entropic Gromov-Wasserstein (GW) optimal transport,
#' so **no shared marker space is required**. This is the route to reusing a
#' reference across different panels or instruments, which the coordinate-matched
#' [somalign_fit()] cannot do.
#'
#' @section Prototype status:
#' This implements *balanced* GW, which transports all query mass and therefore
#' assumes the query and reference share roughly the same populations. Panels that
#' add or drop populations require the unbalanced/partial GW relaxations of Pamona
#' (\doi{10.1093/bioinformatics/btab594}) and SCOTv2 (\doi{10.1089/cmb.2022.0270});
#' those are the intended next step. Treat the correspondence as experimental.
#'
#' Entropic GW is **non-convex**: the alternating-Sinkhorn iteration converges to a
#' local optimum that depends on the initialisation (here the independent coupling
#' `outer(p, q)`) and on `epsilon`. `converged = TRUE` reports only that the outer
#' loop reached a fixed point, not that the alignment is globally optimal or
#' correct. Before trusting a single run on real cross-panel data, compare several
#' `epsilon` values (or add epsilon-annealing / restarts). `transferred_label_confidence`
#' is the top entry of the transport-weighted label posterior, not a calibrated
#' probability, and is inflated toward uniform as `epsilon` grows.
#'
#' @param query,reference `somalign_query` / `somalign_reference` objects (only
#'   their `$codebook` and `$node_masses` are used; the codebooks need not share
#'   markers or dimension).
#' @param epsilon Entropic regularisation strength for the GW solve. Default `0.05`.
#' @param max_iter,tol Outer GW iteration budget and convergence tolerance.
#'
#' @return An object of class `somalign_gw_fit`: `coupling` (query-node by
#'   reference-node transport plan), `correspondence` (row-stochastic coupling),
#'   `converged`, `iterations`, and, when the reference carries labels,
#'   `transferred_label` per query node from `correspondence %*% reference$label_prob`.
#' @examples
#' # (illustrative) align a reference codebook to a rotated copy of itself
#' set.seed(1)
#' cb <- matrix(rnorm(18), 6, 3)
#' rot <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
#' ref <- structure(list(codebook = cb, node_masses = rep(1/6, 6)),
#'                  class = "somalign_reference")
#' qry <- structure(list(codebook = cb %*% rot, node_masses = rep(1/6, 6)),
#'                  class = "somalign_query")
#' fit <- somalign_fit_gw(qry, ref)
#' @export
somalign_fit_gw <- function(query, reference, epsilon = 0.05,
                            max_iter = 50L, tol = 1e-6) {
  qcb <- query$codebook
  rcb <- reference$codebook
  if (is.null(qcb) || is.null(rcb))
    stop("`query` and `reference` must carry a `$codebook`.", call. = FALSE)
  p <- query$node_masses
  q <- reference$node_masses
  if (is.null(p)) p <- rep(1 / nrow(qcb), nrow(qcb))
  if (is.null(q)) q <- rep(1 / nrow(rcb), nrow(rcb))
  p <- p / sum(p)
  q <- q / sum(q)

  # Intra-codebook Euclidean distance matrices (marker spaces may differ).
  C1 <- sqrt(.somalign_pairwise_distance(qcb, qcb))
  C2 <- sqrt(.somalign_pairwise_distance(rcb, rcb))

  gw <- .somalign_gromov_wasserstein(C1, C2, p, q, epsilon = epsilon,
                                     max_iter = max_iter, tol = tol)
  Tplan <- gw$coupling
  row_mass <- rowSums(Tplan)
  correspondence <- Tplan / pmax(row_mass, .Machine$double.eps)

  out <- list(
    coupling = Tplan,
    correspondence = correspondence,
    converged = gw$converged,
    iterations = gw$iterations,
    epsilon = epsilon
  )
  lp <- reference$label_prob
  if (!is.null(lp) && ncol(lp) > 0L) {
    post <- correspondence %*% lp
    out$transferred_label <- colnames(lp)[max.col(post, ties.method = "first")]
    out$transferred_label_confidence <- apply(post, 1, max)
  }
  structure(out, class = "somalign_gw_fit")
}

#' Print a somalign_gw_fit object
#'
#' @param x A `somalign_gw_fit` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @method print somalign_gw_fit
#' @export
print.somalign_gw_fit <- function(x, ...) {
  cat("<somalign_gw_fit> (Gromov-Wasserstein cross-panel alignment, prototype)\n")
  cat(sprintf("  query nodes: %d  ->  reference nodes: %d\n",
              nrow(x$coupling), ncol(x$coupling)))
  cat(sprintf("  epsilon = %g   converged = %s (%d iters)\n",
              x$epsilon, x$converged, x$iterations))
  if (!is.null(x$transferred_label))
    cat(sprintf("  label transfer: %d query nodes labelled\n",
                length(x$transferred_label)))
  invisible(x)
}
