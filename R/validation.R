# Label-transfer validation and calibration.
#
# These functions measure whether transferred labels are correct and whether
# the reported confidence is honest -- the core question for somalign's stated
# aim (label transfer), which the correction/topology diagnostics do not answer.

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

# Confusion matrix with rows = true class, cols = predicted class, over the
# sorted union of observed classes (NA dropped -- callers filter abstentions
# before calling).
.somalign_confusion_matrix <- function(predicted, truth) {
  classes <- sort(unique(c(as.character(truth), as.character(predicted))))
  classes <- classes[!is.na(classes)]
  table(
    true = factor(truth, levels = classes),
    pred = factor(predicted, levels = classes)
  )
}

# Accuracy, macro-F1, multiclass MCC (Gorodkin R_K), and per-class
# precision/recall/F1 from a confusion matrix C (rows = true, cols = pred).
.somalign_confusion_metrics <- function(conf) {
  C <- matrix(as.numeric(conf), nrow = nrow(conf))
  classes <- rownames(conf)
  s <- sum(C)
  correct <- sum(diag(C))
  t_k <- rowSums(C)
  p_k <- colSums(C)
  accuracy <- if (s > 0) correct / s else NA_real_

  # Multiclass MCC (Gorodkin 2004 R_K): covariance of the confusion matrix.
  cov_tp <- correct * s - sum(t_k * p_k)
  cov_pp <- s^2 - sum(p_k^2)
  cov_tt <- s^2 - sum(t_k^2)
  denom <- sqrt(cov_pp * cov_tt)
  mcc <- if (is.finite(denom) && denom > 0) cov_tp / denom else 0

  diag_c <- diag(C)
  precision <- ifelse(p_k > 0, diag_c / p_k, 0)
  recall <- ifelse(t_k > 0, diag_c / t_k, 0)
  f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  present <- t_k > 0
  macro_f1 <- if (any(present)) mean(f1[present]) else NA_real_

  per_class <- data.frame(
    class = classes, precision = precision, recall = recall, f1 = f1,
    support = t_k, row.names = NULL, stringsAsFactors = FALSE
  )
  list(accuracy = accuracy, macro_f1 = macro_f1, mcc = mcc, per_class = per_class)
}

#' Label-transfer accuracy metrics
#'
#' Computes overall accuracy, macro-averaged F1, multiclass Matthews
#' correlation coefficient (Gorodkin's \eqn{R_K}), per-class
#' precision/recall/F1, the confusion matrix, and coverage (the fraction of
#' cells with an accepted prediction), given predicted and ground-truth labels.
#'
#' Metrics are computed on the *accepted* predictions only; `coverage` reports
#' what fraction of cells that is. `accuracy_all` additionally scores
#' abstentions (rejected or `NA` predictions) as wrong, for a
#' coverage-penalised view.
#'
#' @param predicted Character vector of predicted labels (`NA` = abstain).
#' @param truth Character vector of ground-truth labels, same length.
#' @param accepted Optional logical vector, same length. When supplied, only
#'   `accepted` (and non-`NA`) predictions are scored; the rest are abstentions.
#'   When `NULL` (default), all non-`NA` predictions are scored.
#'
#' @return A list of class `somalign_label_metrics` with `accuracy`,
#'   `macro_f1`, `mcc`, `per_class` (data frame), `confusion` (table), `n`
#'   (scored predictions), `coverage`, and `accuracy_all`.
#' @examples
#' truth <- rep(c("A", "B", "C"), each = 10)
#' pred  <- truth; pred[c(1, 12, 25)] <- "B"
#' somalign_label_metrics(pred, truth)
#' @export
somalign_label_metrics <- function(predicted, truth, accepted = NULL) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  if (length(predicted) != length(truth))
    stop("`predicted` and `truth` must have the same length.", call. = FALSE)
  if (is.null(accepted)) {
    accepted <- !is.na(predicted)
  } else {
    if (length(accepted) != length(truth))
      stop("`accepted` must have the same length as `truth`.", call. = FALSE)
    accepted <- as.logical(accepted) & !is.na(predicted)
  }
  coverage <- if (length(accepted)) mean(accepted) else NA_real_
  accuracy_all <- if (length(truth)) mean(accepted & predicted == truth) else NA_real_
  conf <- .somalign_confusion_matrix(predicted[accepted], truth[accepted])
  m <- .somalign_confusion_metrics(conf)
  structure(
    list(accuracy = m$accuracy, macro_f1 = m$macro_f1, mcc = m$mcc,
         per_class = m$per_class, confusion = conf, n = sum(accepted),
         coverage = coverage, accuracy_all = accuracy_all),
    class = "somalign_label_metrics"
  )
}

#' @method print somalign_label_metrics
#' @export
print.somalign_label_metrics <- function(x, ...) {
  cat("<somalign_label_metrics>\n")
  cat(sprintf("  accuracy = %.4f  macro_f1 = %.4f  MCC = %.4f\n",
              x$accuracy, x$macro_f1, x$mcc))
  cat(sprintf("  scored = %d  coverage = %.1f%%  accuracy_all = %.4f\n",
              x$n, 100 * x$coverage, x$accuracy_all))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------

#' Confidence calibration of label transfer
#'
#' Bins predictions by a confidence-like score in \[0, 1\] (e.g. the transfer
#' confidence or margin) and compares mean score to empirical accuracy per bin,
#' yielding a reliability table plus the expected (ECE) and maximum (MCE)
#' calibration error and the top-label Brier score. A well-calibrated model has
#' mean score ~ accuracy in every bin (ECE ~ 0).
#'
#' @param score Numeric vector in \[0, 1\]: the confidence/margin per prediction.
#' @param correct Logical vector, same length: whether each prediction was right.
#' @param n_bins Positive integer. Number of equal-width bins over \[0, 1\].
#'   Default `10L`.
#'
#' @return A list of class `somalign_calibration` with `table` (per-bin
#'   `score_mean`, `accuracy`, `n`), `ece`, `mce`, `brier`, and `n`.
#' @examples
#' set.seed(1)
#' score <- runif(200)
#' correct <- runif(200) < score        # perfectly calibrated by construction
#' somalign_calibration(score, correct)
#' @export
somalign_calibration <- function(score, correct, n_bins = 10L) {
  ok <- !is.na(score) & !is.na(correct)
  score <- pmin(pmax(as.numeric(score)[ok], 0), 1)
  correct <- as.numeric(as.logical(correct)[ok])
  if (length(score) == 0L)
    stop("No non-missing (score, correct) pairs to calibrate.", call. = FALSE)
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin <- cut(score, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  rows <- lapply(seq_len(n_bins), function(b) {
    sel <- bin == b
    if (!any(sel)) return(NULL)
    data.frame(bin = b, score_mean = mean(score[sel]),
               accuracy = mean(correct[sel]), n = sum(sel))
  })
  tbl <- do.call(rbind, rows)
  w <- tbl$n / sum(tbl$n)
  gap <- abs(tbl$accuracy - tbl$score_mean)
  structure(
    list(table = tbl, ece = sum(w * gap), mce = max(gap),
         brier = mean((score - correct)^2), n = length(score)),
    class = "somalign_calibration"
  )
}

#' @method print somalign_calibration
#' @export
print.somalign_calibration <- function(x, ...) {
  cat("<somalign_calibration>\n")
  cat(sprintf("  ECE = %.4f  MCE = %.4f  Brier = %.4f  (n = %d)\n",
              x$ece, x$mce, x$brier, x$n))
  cat("  reliability (score_mean -> accuracy, n):\n")
  for (i in seq_len(nrow(x$table))) {
    cat(sprintf("    %.2f -> %.2f  (%d)\n",
                x$table$score_mean[i], x$table$accuracy[i], x$table$n[i]))
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Cross-validation
# ---------------------------------------------------------------------------

# Stratified k-fold assignment: each class is spread round-robin across folds
# so every fold sees every (sufficiently common) class. Assumes the RNG has
# already been seeded by the caller.
.somalign_stratified_folds <- function(labels, k) {
  folds <- integer(length(labels))
  for (cl in unique(labels)) {
    idx <- which(labels == cl)
    idx <- idx[sample.int(length(idx))]
    folds[idx] <- rep_len(seq_len(k), length(idx))
  }
  folds
}

# One CV fold: train a reference on the training split, project the held-out
# split, transfer labels, and return per-cell (truth, predicted, confidence,
# margin, accepted).
.somalign_cv_fold <- function(fold, folds, data, labels, grid, rlen,
                              epsilon, solver, fit_args) {
  test_idx <- which(folds == fold)
  train_idx <- which(folds != fold)
  ref <- somalign_train_reference(data[train_idx, , drop = FALSE],
                                  labels = labels[train_idx], grid = grid, rlen = rlen)
  qry <- somalign_query(data[test_idx, , drop = FALSE], ref, grid = grid, rlen = rlen)
  fit <- do.call(somalign_fit,
                 c(list(qry, ref, epsilon = epsilon, solver = solver), fit_args))
  res <- somalign_results(fit, include_correction = FALSE)
  data.frame(fold = fold, truth = labels[test_idx],
             predicted = res$transferred_label,
             confidence = res$transferred_label_confidence,
             margin = res$transferred_label_margin,
             accepted = res$transferred_label_accepted,
             stringsAsFactors = FALSE)
}

#' Cross-validate label transfer on held-out cells
#'
#' Stratified k-fold cross-validation that measures how accurately somalign
#' transfers labels to cells it has not seen. Each fold trains a reference SOM
#' on the training split, projects the held-out split, and transfers labels;
#' every held-out cell has a real ground-truth label, so this estimates true
#' generalisation without any external labelled query data. Pooled results are
#' scored with [somalign_label_metrics()] and [somalign_calibration()].
#'
#' @param data Numeric cell-by-feature matrix.
#' @param labels Character vector of per-cell labels, `nrow(data)` long.
#' @param grid A `kohonen::somgrid` for the reference and query SOMs.
#' @param k Positive integer. Number of folds. Default `5L`.
#' @param stratify Logical. Stratify folds by class. Default `TRUE`.
#' @param epsilon,solver Passed to [somalign_fit()].
#' @param rlen SOM training iterations for both SOMs. Default `20L`.
#' @param n_bins Calibration bins. Default `10L`.
#' @param seed Integer or `NULL`. RNG seed, restored on exit. Default `1L`.
#' @param ... Further arguments forwarded to [somalign_fit()].
#'
#' @return A list of class `somalign_cross_validation` with `metrics`
#'   ([somalign_label_metrics()]), `calibration` ([somalign_calibration()]),
#'   `per_fold` (data frame), `predictions` (pooled per-cell data frame), `k`.
#' @examples
#' \donttest{
#' if (requireNamespace("kohonen", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- rbind(matrix(rnorm(200, -2), ncol = 2), matrix(rnorm(200, 2), ncol = 2))
#'   colnames(x) <- c("f1", "f2")
#'   lab <- rep(c("low", "high"), each = 100)
#'   cv <- somalign_cross_validate(x, lab,
#'     grid = kohonen::somgrid(2, 2, "hexagonal"), k = 3)
#'   cv$metrics$accuracy
#' }
#' }
#' @export
somalign_cross_validate <- function(data, labels, grid, k = 5L, stratify = TRUE,
                                    epsilon = 0.1, solver = "internal",
                                    rlen = 20L, n_bins = 10L, seed = 1L, ...) {
  data <- as.matrix(data)
  labels <- as.character(labels)
  if (nrow(data) != length(labels))
    stop("`labels` must have one entry per row of `data`.", call. = FALSE)
  if (!is.null(seed)) withr::local_seed(seed)

  folds <- if (isTRUE(stratify)) .somalign_stratified_folds(labels, k)
           else rep_len(seq_len(k), nrow(data))[sample.int(nrow(data))]
  fit_args <- list(...)

  preds <- lapply(seq_len(k), .somalign_cv_fold, folds = folds, data = data,
                  labels = labels, grid = grid, rlen = rlen, epsilon = epsilon,
                  solver = solver, fit_args = fit_args)
  pred_df <- do.call(rbind, preds)

  metrics <- somalign_label_metrics(pred_df$predicted, pred_df$truth, pred_df$accepted)
  calib <- somalign_calibration(pred_df$confidence,
                                pred_df$predicted == pred_df$truth, n_bins = n_bins)
  per_fold <- do.call(rbind, lapply(split(pred_df, pred_df$fold), function(d) {
    m <- somalign_label_metrics(d$predicted, d$truth, d$accepted)
    data.frame(fold = d$fold[1], accuracy = m$accuracy, macro_f1 = m$macro_f1,
               mcc = m$mcc, coverage = m$coverage)
  }))
  structure(
    list(metrics = metrics, calibration = calib, per_fold = per_fold,
         predictions = pred_df, k = k),
    class = "somalign_cross_validation"
  )
}

#' @method print somalign_cross_validation
#' @export
print.somalign_cross_validation <- function(x, ...) {
  cat(sprintf("<somalign_cross_validation> (%d folds)\n", x$k))
  cat(sprintf("  pooled: accuracy = %.4f  macro_f1 = %.4f  MCC = %.4f  coverage = %.1f%%\n",
              x$metrics$accuracy, x$metrics$macro_f1, x$metrics$mcc,
              100 * x$metrics$coverage))
  cat(sprintf("  calibration: ECE = %.4f  Brier = %.4f\n",
              x$calibration$ece, x$calibration$brier))
  invisible(x)
}
