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

# ---------------------------------------------------------------------------
# Supervised plan tuning
# ---------------------------------------------------------------------------

.somalign_get <- function(x, name, default) if (is.null(x[[name]])) default else x[[name]]

# Normalise a param grid into a list of named-list combinations. Accepts a
# data frame (each row -> combo of its scalar columns) or a list of named
# lists (used as-is). Every combo must specify `epsilon`.
.somalign_normalize_param_grid <- function(param_grid) {
  if (is.data.frame(param_grid)) {
    combos <- lapply(seq_len(nrow(param_grid)), function(i) as.list(param_grid[i, , drop = FALSE]))
  } else if (is.list(param_grid)) {
    combos <- param_grid
  } else {
    stop("`param_grid` must be a data frame or a list of named lists.", call. = FALSE)
  }
  if (!all(vapply(combos, function(p) "epsilon" %in% names(p), logical(1))))
    stop("Every parameter combination must specify `epsilon`.", call. = FALSE)
  combos
}

# Evaluate one parameter combination on a pre-built (query, reference) fold:
# OT-only solve + label transfer (no per-cell projection or node shifts),
# returning per-test-cell predicted label, acceptance, and confidence.
.somalign_tune_labels <- function(qry, ref, params, min_match_fraction,
                                  confidence_threshold, solver, max_iter, tol) {
  r <- .somalign_ot_sweep_one(
    qry, ref, epsilon = params$epsilon,
    rho_query = .somalign_get(params, "rho_query", 1),
    rho_ref = .somalign_get(params, "rho_ref", 1),
    solver = solver, max_iter = max_iter, tol = tol,
    diagonal_boost = .somalign_get(params, "diagonal_boost", 0),
    feature_weights = params$feature_weights
  )
  cc <- .somalign_plan_to_correspondence(r$plan, qry$node_masses)
  lt <- .somalign_transfer_labels(cc$correspondence, ref$label_prob,
                                  cc$match_fraction, min_match_fraction,
                                  confidence_threshold)
  unit <- qry$sample_unit
  accepted <- lt$accepted[unit]
  list(predicted = ifelse(accepted, lt$label[unit], NA_character_),
       accepted = accepted,
       confidence = ifelse(accepted, lt$confidence[unit], NA_real_))
}

# Score one combo's pooled cross-fold predictions into a one-row summary.
.somalign_tune_row <- function(params, preds, n_bins) {
  pdf <- do.call(rbind, preds)
  m <- somalign_label_metrics(pdf$predicted, pdf$truth, pdf$accepted)
  correct <- pdf$predicted == pdf$truth
  ece <- tryCatch(somalign_calibration(pdf$confidence, correct, n_bins = n_bins)$ece,
                  error = function(e) NA_real_)
  data.frame(
    epsilon = params$epsilon,
    rho_query = .somalign_get(params, "rho_query", 1),
    rho_ref = .somalign_get(params, "rho_ref", 1),
    diagonal_boost = .somalign_get(params, "diagonal_boost", 0),
    feature_weights = if (is.null(params$feature_weights)) "none" else "custom",
    accuracy = m$accuracy, macro_f1 = m$macro_f1, mcc = m$mcc,
    coverage = m$coverage, ece = ece, stringsAsFactors = FALSE
  )
}

#' Tune transport-plan knobs against label-transfer accuracy
#'
#' The supervised counterpart to [somalign_select_epsilon()]: instead of an
#' unsupervised plan-geometry criterion, it selects transport-plan parameters
#' by cross-validated label-transfer performance. For each parameter
#' combination it runs stratified k-fold CV (reusing pre-trained SOMs per fold
#' for efficiency -- plan knobs change only the OT solve, not the SOMs) and
#' scores pooled held-out predictions with [somalign_label_metrics()].
#'
#' Tunable knobs are those that shape the transport plan without needing
#' anchors or query labels: `epsilon`, `rho_query`, `rho_ref`,
#' `diagonal_boost`, and `feature_weights` (a numeric per-feature vector).
#' `label_guided` and `rho_anchor` are out of scope here (they require a
#' labelled query SOM or anchor pairs, respectively).
#'
#' Note that `"mcc"` and `"accuracy"` are scored on *accepted* predictions only,
#' so they can be inflated by settings that abstain on hard cells (higher
#' `epsilon` raises accuracy while dropping `coverage`). Always read the
#' `coverage` and `macro_f1` columns alongside the objective: `macro_f1` falls
#' when rare classes are abstained away, making it a more coverage-robust
#' target for imbalanced data.
#'
#' @param data Numeric cell-by-feature matrix.
#' @param labels Character vector of per-cell labels.
#' @param grid A `kohonen::somgrid`.
#' @param param_grid A data frame (one row per combination of the scalar knobs
#'   `epsilon`, `rho_query`, `rho_ref`, `diagonal_boost`) or a list of named
#'   lists (which may additionally carry `feature_weights`). Each must specify
#'   `epsilon`.
#' @param k Folds. Default `5L`.
#' @param metric Objective to optimise: `"mcc"` (default), `"macro_f1"`,
#'   `"accuracy"` (all maximised) or `"ece"` (minimised).
#' @param stratify,rlen,seed As in [somalign_cross_validate()].
#' @param min_match_fraction,confidence_threshold Label-acceptance gates,
#'   matching [somalign_fit()] defaults.
#' @param solver,max_iter,tol OT solver settings.
#' @param n_bins Calibration bins for the `ece` column. Default `10L`.
#'
#' @return A list of class `somalign_tune` with `best` (the winning combo as a
#'   one-row data frame), `best_params` (named list), `grid` (all combinations
#'   with their CV metrics), and `metric`.
#' @examples
#' \donttest{
#' if (requireNamespace("kohonen", quietly = TRUE)) {
#'   set.seed(1)
#'   x <- rbind(matrix(rnorm(200, -2), ncol = 2), matrix(rnorm(200, 2), ncol = 2))
#'   colnames(x) <- c("f1", "f2")
#'   lab <- rep(c("low", "high"), each = 100)
#'   tuned <- somalign_tune(x, lab, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'     param_grid = data.frame(epsilon = c(0.05, 0.1, 0.2)), k = 3)
#'   tuned$best_params$epsilon
#' }
#' }
#' @export
somalign_tune <- function(data, labels, grid, param_grid, k = 5L,
                          metric = c("mcc", "macro_f1", "accuracy", "ece"),
                          stratify = TRUE, rlen = 20L,
                          min_match_fraction = 0.05, confidence_threshold = 0.6,
                          solver = "internal", max_iter = 1000, tol = 1e-7,
                          n_bins = 10L, seed = 1L) {
  metric <- match.arg(metric)
  data <- as.matrix(data)
  labels <- as.character(labels)
  if (nrow(data) != length(labels))
    stop("`labels` must have one entry per row of `data`.", call. = FALSE)
  combos <- .somalign_normalize_param_grid(param_grid)
  if (!is.null(seed)) withr::local_seed(seed)

  folds <- if (isTRUE(stratify)) .somalign_stratified_folds(labels, k)
           else rep_len(seq_len(k), nrow(data))[sample.int(nrow(data))]
  fold_objs <- lapply(seq_len(k), function(f) {
    test_idx <- which(folds == f); train_idx <- which(folds != f)
    ref <- somalign_train_reference(data[train_idx, , drop = FALSE],
                                    labels = labels[train_idx], grid = grid, rlen = rlen)
    qry <- somalign_query(data[test_idx, , drop = FALSE], ref, grid = grid, rlen = rlen)
    list(qry = qry, ref = ref, truth = labels[test_idx])
  })

  rows <- lapply(combos, function(params) {
    preds <- lapply(fold_objs, function(fo) {
      ev <- .somalign_tune_labels(fo$qry, fo$ref, params, min_match_fraction,
                                  confidence_threshold, solver, max_iter, tol)
      data.frame(truth = fo$truth, predicted = ev$predicted,
                 accepted = ev$accepted, confidence = ev$confidence,
                 stringsAsFactors = FALSE)
    })
    .somalign_tune_row(params, preds, n_bins)
  })
  grid_tbl <- do.call(rbind, rows)

  objective <- grid_tbl[[metric]]
  best_idx <- if (identical(metric, "ece")) which.min(objective) else which.max(objective)
  structure(
    list(best = grid_tbl[best_idx, , drop = FALSE],
         best_params = combos[[best_idx]], grid = grid_tbl, metric = metric),
    class = "somalign_tune"
  )
}

#' @method print somalign_tune
#' @export
print.somalign_tune <- function(x, ...) {
  cat(sprintf("<somalign_tune> objective = %s over %d combination(s)\n",
              x$metric, nrow(x$grid)))
  cat("  best:\n")
  b <- x$best
  cat(sprintf("    epsilon=%.4g rho_query=%.4g rho_ref=%.4g diagonal_boost=%.4g fw=%s\n",
              b$epsilon, b$rho_query, b$rho_ref, b$diagonal_boost, b$feature_weights))
  cat(sprintf("    accuracy=%.4f macro_f1=%.4f MCC=%.4f coverage=%.1f%% ECE=%.4f\n",
              b$accuracy, b$macro_f1, b$mcc, 100 * b$coverage, b$ece))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Anchor benefit
# ---------------------------------------------------------------------------

# Transfer labels for one cost-bonus setting via the OT-only path (no per-cell
# projection or node shifts -- labels don't need them) and score them against
# ground truth on the eval subset.
.somalign_anchor_eval <- function(query, reference, cost_bonus, eval_mask,
                                  query_labels, min_match_fraction, confidence_threshold,
                                  epsilon, rho_query, rho_ref, solver, max_iter, tol, n_bins) {
  r <- .somalign_ot_sweep_one(
    query, reference, epsilon = epsilon, rho_query = rho_query, rho_ref = rho_ref,
    solver = solver, max_iter = max_iter, tol = tol, cost_bonus = cost_bonus
  )
  cc <- .somalign_plan_to_correspondence(r$plan, query$node_masses)
  lt <- .somalign_transfer_labels(cc$correspondence, reference$label_prob,
                                  cc$match_fraction, min_match_fraction, confidence_threshold)
  unit <- query$sample_unit
  accepted <- lt$accepted[unit]
  predicted <- ifelse(accepted, lt$label[unit], NA_character_)[eval_mask]
  confidence <- ifelse(accepted, lt$confidence[unit], NA_real_)[eval_mask]
  accepted <- accepted[eval_mask]
  truth <- query_labels[eval_mask]
  m <- somalign_label_metrics(predicted, truth, accepted)
  ece <- tryCatch(somalign_calibration(confidence, predicted == truth, n_bins = n_bins)$ece,
                  error = function(e) NA_real_)
  data.frame(accuracy = m$accuracy, macro_f1 = m$macro_f1, mcc = m$mcc,
             coverage = m$coverage, ece = ece)
}

#' Quantify the label-transfer benefit of anchor (repeat) samples
#'
#' Measures how much anchor samples -- QC/repeat specimens run in both the
#' reference and query batches -- improve label transfer, by sweeping the anchor
#' cost-bonus strength `rho_anchor` and scoring the transferred labels against a
#' *known* query label. `rho_anchor = 0` is the exact no-anchor baseline
#' (plain [somalign_fit()]); positive values bias the transport plan toward
#' anchor-supported node pairs, exactly as `somalign_fit_anchored(correction =
#' "cost_bonus")`. Anchors influence labels only through this cost bonus, so
#' this single sweep captures their entire label-transfer effect (the
#' `"subspace"` correction leaves labels identical to `rho_anchor = 0`).
#'
#' This is a *validation* tool: it needs ground-truth query labels (e.g. an
#' independently gated held-out batch), which production label transfer does
#' not have. It reuses the fixed reference/query SOMs and a single anchor count
#' matrix across the sweep, recomputing only the OT solve, so it is cheap.
#'
#' @param query A `somalign_query` object.
#' @param reference A labelled `somalign_reference` object.
#' @param query_labels Character ground-truth labels, one per query cell.
#' @param anchor_old,anchor_new Paired anchor cell matrices (repeat samples in
#'   the reference and query batches respectively), same number of rows.
#' @param rho_grid Non-negative anchor strengths to sweep. Default
#'   `c(0, 0.5, 1, 2, 5)`; include `0` for the no-anchor baseline.
#' @param metric Objective for `best`/`lift`: `"mcc"` (default), `"macro_f1"`,
#'   or `"accuracy"` (all maximised).
#' @param eval_mask Optional logical vector, one per query cell, selecting the
#'   cells to score (e.g. exclude the anchor samples for a clean held-out
#'   measure). Default `NULL` scores all cells.
#' @param epsilon,rho_query,rho_ref,solver,max_iter,tol OT settings.
#' @param min_match_fraction,confidence_threshold Label-acceptance gates.
#' @param chunk_size Anchor-projection chunk size.
#' @param n_bins Calibration bins for the `ece` column.
#'
#' @return A list of class `somalign_anchor_benefit` with `grid` (one row per
#'   `rho_anchor`: accuracy, macro_f1, mcc, coverage, ece), `baseline` (the
#'   `rho_anchor = 0` row), `best`, `lift` (best minus baseline on `metric`),
#'   and `metric`.
#' @seealso [somalign_fit_anchored()], [somalign_cross_validate()]
#' @export
somalign_anchor_benefit <- function(query, reference, query_labels,
                                    anchor_old, anchor_new,
                                    rho_grid = c(0, 0.5, 1, 2, 5),
                                    metric = c("mcc", "macro_f1", "accuracy"),
                                    eval_mask = NULL, epsilon = 0.1,
                                    rho_query = 1, rho_ref = 1,
                                    min_match_fraction = 0.05,
                                    confidence_threshold = 0.6, solver = "internal",
                                    max_iter = 1000, tol = 1e-7,
                                    chunk_size = 10000L, n_bins = 10L) {
  metric <- match.arg(metric)
  if (is.null(reference$label_prob) || ncol(reference$label_prob) == 0L)
    stop("`reference` must carry labels for anchor-benefit evaluation.", call. = FALSE)
  n_q <- nrow(query$scaled_data)
  query_labels <- as.character(query_labels)
  if (length(query_labels) != n_q)
    stop("`query_labels` must have one entry per query cell.", call. = FALSE)
  if (is.null(eval_mask)) {
    eval_mask <- rep(TRUE, n_q)
  } else {
    eval_mask <- as.logical(eval_mask)
    if (length(eval_mask) != n_q)
      stop("`eval_mask` must have one entry per query cell.", call. = FALSE)
  }
  if (any(rho_grid < 0)) stop("`rho_grid` values must be non-negative.", call. = FALSE)

  anc <- .somalign_validate_anchors(anchor_old, anchor_new, reference)
  base_bonus <- .somalign_anchor_cost_bonus(
    anc$anchor_old_scaled, anc$anchor_new_scaled,
    query$codebook, reference$codebook, rho_anchor = 1, chunk_size = chunk_size
  )$bonus

  rows <- lapply(rho_grid, function(rho) {
    cb <- if (rho > 0) rho * base_bonus else NULL
    cbind(rho_anchor = rho,
          .somalign_anchor_eval(query, reference, cb, eval_mask, query_labels,
                                min_match_fraction, confidence_threshold, epsilon,
                                rho_query, rho_ref, solver, max_iter, tol, n_bins))
  })
  grid_tbl <- do.call(rbind, rows)
  baseline <- grid_tbl[grid_tbl$rho_anchor == 0, , drop = FALSE]
  best <- grid_tbl[which.max(grid_tbl[[metric]]), , drop = FALSE]
  lift <- if (nrow(baseline) > 0) best[[metric]] - baseline[[metric]][1] else NA_real_
  structure(list(grid = grid_tbl, baseline = baseline, best = best,
                 lift = lift, metric = metric),
            class = "somalign_anchor_benefit")
}

#' @method print somalign_anchor_benefit
#' @export
print.somalign_anchor_benefit <- function(x, ...) {
  cat(sprintf("<somalign_anchor_benefit> objective = %s over %d rho value(s)\n",
              x$metric, nrow(x$grid)))
  if (nrow(x$baseline) > 0)
    cat(sprintf("  baseline (rho=0): %s = %.4f  coverage = %.1f%%\n",
                x$metric, x$baseline[[x$metric]][1], 100 * x$baseline$coverage[1]))
  cat(sprintf("  best: rho=%.3g  %s = %.4f  (lift %+.4f)  coverage = %.1f%%  ECE = %.4f\n",
              x$best$rho_anchor, x$metric, x$best[[x$metric]], x$lift,
              100 * x$best$coverage, x$best$ece))
  invisible(x)
}
