#!/usr/bin/env Rscript
# =============================================================================
# BMV pilot: experimental validation of OtherT label-transfer improvements.
#
# This script keeps one labelled acquisition batch as the reference and treats a
# second labelled acquisition batch as a query with held-out truth. It evaluates
# seven improvement ideas discussed in the publication-readiness review:
# hard-vs-soft abundance frequencies, targeted OtherT tuning, epsilon changes,
# marker-aware transport costs, T/NK-local projection, confidence triage, and
# anchor-strength tuning.
# =============================================================================

.libPaths(c("/exports/para-lipg-hpc/mdmanurung/R/4.5", .libPaths()))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE) else NA_character_
SCRIPT_DIR <- if (is.na(SCRIPT)) getwd() else dirname(SCRIPT)
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(ROOT, quiet = TRUE)
  } else {
    library(somalign)
  }
  library(kohonen)
})

`%||%` <- function(x, y) if (length(x) == 0L || is.na(x) || !nzchar(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
n_sub <- as.integer(args[1] %||% "60000")
grid_side <- as.integer(args[2] %||% "15")
rlen <- as.integer(args[3] %||% "30")
seed <- as.integer(args[4] %||% "1")
set.seed(seed)

PILOT <- Sys.getenv(
  "SOMALIGN_PILOT_QS",
  "/exports/para-lipg-hpc/Xuran/data/pilot/aurora_exvivo/processed/aurora_exvivo_xyfclusters.qs"
)
OUT_DIR <- SCRIPT_DIR
OUT_RDS <- file.path(OUT_DIR, "other_t_improvement_results.rds")

FEATURES <- c("TCRgd_CD141", "CD8", "CD3", "CD19", "CD159c_NKG2C", "CD38",
              "CD57", "CD161", "CD1c", "HLA-DR", "CD11c", "CD25", "CD16",
              "CD185_CXCR5", "CD14", "CD303_TCRvd2", "CD4", "CD294_CRTH2",
              "CD21", "CD56", "FoxP3", "Tbet", "CD197_CCR7", "CD127",
              "CD45RA", "CD27", "CD7")
LABEL_COL <- Sys.getenv("SOMALIGN_LABEL_COL", "gate_lineage2")
BATCH_COL <- Sys.getenv("SOMALIGN_BATCH_COL", "batch_id")
REF_BATCH <- Sys.getenv("SOMALIGN_REF_BATCH", "Batch1")
QRY_BATCH <- Sys.getenv("SOMALIGN_QRY_BATCH", "Batch2")

T_BOUNDARY_LABELS <- c("CD4T", "CD8T", "OtherT", "NK")
T_FEATURES <- c("TCRgd_CD141", "CD8", "CD3", "CD159c_NKG2C", "CD38", "CD57",
                "CD161", "CD25", "CD16", "CD185_CXCR5", "CD4", "CD56",
                "FoxP3", "Tbet", "CD197_CCR7", "CD127", "CD45RA", "CD27",
                "CD7", "HLA-DR")
T_FEATURES <- intersect(T_FEATURES, FEATURES)

stamp <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

write_table <- function(x, name) {
  utils::write.csv(x, file.path(OUT_DIR, name), row.names = FALSE)
}

safe_num <- function(x) ifelse(length(x) && is.finite(x), x, NA_real_)

full_per_class_metrics <- function(predicted, truth, accepted = NULL, classes = NULL) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  if (is.null(accepted)) accepted <- !is.na(predicted)
  accepted <- as.logical(accepted) & !is.na(predicted)
  if (is.null(classes)) classes <- sort(unique(c(truth, predicted[accepted])))
  classes <- classes[!is.na(classes)]
  rows <- lapply(classes, function(cl) {
    tp <- sum(accepted & predicted == cl & truth == cl, na.rm = TRUE)
    fp <- sum(accepted & predicted == cl & truth != cl, na.rm = TRUE)
    fn <- sum(truth == cl & !(accepted & predicted == cl), na.rm = TRUE)
    support <- sum(truth == cl, na.rm = TRUE)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
    data.frame(
      class = cl,
      precision = precision,
      recall = recall,
      f1 = f1,
      support = support,
      predicted = tp + fp,
      abstained_truth = sum(truth == cl & !accepted, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

score_prediction <- function(approach, config, predicted, truth, accepted = NULL,
                             confidence = NULL, margin = NULL, subset_name = "all",
                             note = "", classes = NULL,
                             evidence_class = "deployable") {
  if (is.null(accepted)) accepted <- !is.na(predicted)
  accepted <- as.logical(accepted) & !is.na(predicted)
  pred_scored <- ifelse(accepted, predicted, NA_character_)
  m <- somalign_label_metrics(pred_scored, truth, accepted)
  full_pc <- full_per_class_metrics(predicted, truth, accepted, classes = classes)
  accepted_pc <- m$per_class
  names(accepted_pc)[names(accepted_pc) %in% c("precision", "recall", "f1", "support")] <-
    paste0("accepted_", names(accepted_pc)[names(accepted_pc) %in% c("precision", "recall", "f1", "support")])
  pc <- merge(full_pc, accepted_pc, by = "class", all.x = TRUE)
  pc$approach <- approach
  pc$config <- config
  pc$subset <- subset_name
  pc$evidence_class <- evidence_class
  pc$note <- note
  pc <- pc[, c("approach", "config", "subset", "class", "precision", "recall",
               "f1", "support", "predicted", "abstained_truth",
               "accepted_precision", "accepted_recall", "accepted_f1",
               "accepted_support", "evidence_class", "note")]

  pick <- function(cls, col) {
    i <- match(cls, pc$class)
    if (is.na(i)) NA_real_ else pc[[col]][i]
  }
  calib <- list(ece = NA_real_, brier = NA_real_)
  ok_cal <- accepted & !is.na(confidence)
  if (sum(ok_cal) > 0L) {
    calib <- tryCatch(
      somalign_calibration(confidence[ok_cal], predicted[ok_cal] == truth[ok_cal]),
      error = function(e) list(ece = NA_real_, brier = NA_real_)
    )
  }
  full_macro_f1 <- mean(pc$f1[pc$support > 0], na.rm = TRUE)
  full_macro_precision <- mean(pc$precision[pc$support > 0], na.rm = TRUE)
  full_macro_recall <- mean(pc$recall[pc$support > 0], na.rm = TRUE)
  summary <- data.frame(
    approach = approach,
    config = config,
    subset = subset_name,
    evidence_class = evidence_class,
    is_deployable = identical(evidence_class, "deployable"),
    n_total = length(truth),
    n_accepted = sum(accepted),
    coverage = mean(accepted),
    accepted_accuracy = m$accuracy,
    accepted_macro_f1 = m$macro_f1,
    accepted_mcc = m$mcc,
    accuracy_all = m$accuracy_all,
    full_macro_precision = full_macro_precision,
    full_macro_recall = full_macro_recall,
    full_macro_f1 = full_macro_f1,
    otherT_precision = pick("OtherT", "precision"),
    otherT_recall = pick("OtherT", "recall"),
    otherT_f1 = pick("OtherT", "f1"),
    nk_precision = pick("NK", "precision"),
    nk_recall = pick("NK", "recall"),
    nk_f1 = pick("NK", "f1"),
    cd4_precision = pick("CD4T", "precision"),
    cd4_recall = pick("CD4T", "recall"),
    cd4_f1 = pick("CD4T", "f1"),
    cd8_precision = pick("CD8T", "precision"),
    cd8_recall = pick("CD8T", "recall"),
    cd8_f1 = pick("CD8T", "f1"),
    median_confidence = stats::median(confidence[accepted], na.rm = TRUE),
    median_margin = stats::median(margin[accepted], na.rm = TRUE),
    ece = calib$ece,
    brier = calib$brier,
    note = note,
    stringsAsFactors = FALSE
  )
  list(summary = summary, per_class = pc)
}

labels_from_fit <- function(fit) {
  unit <- fit$query$sample_unit
  lt <- fit$label_transfer
  accepted <- lt$accepted[unit]
  accepted[is.na(accepted)] <- FALSE
  second_conf <- lt$second_confidence[unit]
  label <- lt$label[unit]
  confidence <- lt$confidence[unit]
  data.frame(
    query_som_unit = unit,
    transferred_label = ifelse(accepted, label, NA_character_),
    raw_transferred_label = label,
    transferred_label_accepted = accepted,
    transferred_label_confidence = ifelse(accepted, confidence, NA_real_),
    transferred_label_second = lt$second_label[unit],
    transferred_label_second_confidence = second_conf,
    transferred_label_margin = confidence - ifelse(is.na(second_conf), 0, second_conf),
    stringsAsFactors = FALSE
  )
}

label_prob_from_units <- function(labels, units, n_nodes, class_levels) {
  tab <- table(
    unit = factor(units, levels = seq_len(n_nodes)),
    label = factor(labels, levels = class_levels)
  )
  mat <- matrix(as.numeric(tab), nrow = nrow(tab), dimnames = dimnames(tab))
  rs <- rowSums(mat)
  mat[rs > 0, ] <- mat[rs > 0, , drop = FALSE] / rs[rs > 0]
  mat
}

make_t_marker_weights <- function(features) {
  w <- setNames(rep(1, length(features)), features)
  high <- intersect(c("CD3", "CD4", "CD8", "TCRgd_CD141", "CD56", "CD16",
                      "CD159c_NKG2C", "CD161", "FoxP3", "Tbet",
                      "CD197_CCR7", "CD45RA", "CD27", "CD7", "CD25"), names(w))
  low <- intersect(c("CD19", "CD21", "CD14", "CD1c", "CD303_TCRvd2", "CD11c"), names(w))
  w[high] <- 2
  w[low] <- 0.5
  w / mean(w)
}

make_anchors <- function(X_old, y_old, X_new, y_new, n_total) {
  classes <- intersect(sort(unique(y_old)), sort(unique(y_new)))
  per_cls <- max(1L, floor(n_total / length(classes)))
  ao <- list()
  an <- list()
  new_idx <- integer(0)
  for (cl in classes) {
    io <- which(y_old == cl)
    inew <- which(y_new == cl)
    k <- min(per_cls, length(io), length(inew))
    if (k < 1L) next
    so <- sample(io, k)
    sn <- sample(inew, k)
    ao[[cl]] <- X_old[so, , drop = FALSE]
    an[[cl]] <- X_new[sn, , drop = FALSE]
    new_idx <- c(new_idx, sn)
  }
  list(anchor_old = do.call(rbind, ao), anchor_new = do.call(rbind, an),
       query_anchor_idx = new_idx)
}

make_random_anchors <- function(X_old, X_new, n_total) {
  k <- min(n_total, nrow(X_old), nrow(X_new))
  old_idx <- sample.int(nrow(X_old), k)
  new_idx <- sample.int(nrow(X_new), k)
  list(
    anchor_old = X_old[old_idx, , drop = FALSE],
    anchor_new = X_new[new_idx, , drop = FALSE],
    query_anchor_idx = new_idx
  )
}

count_matrix <- function(group, value, levels) {
  group <- as.character(group)
  value <- as.character(value)
  levels <- as.character(levels)
  ok <- !is.na(group) & !is.na(value)
  group <- group[ok]
  value <- value[ok]
  tab <- table(
    group = factor(group, levels = sort(unique(group))),
    value = factor(value, levels = levels)
  )
  matrix(as.numeric(tab), nrow = nrow(tab), dimnames = dimnames(tab))
}

clr_rows <- function(mat, pseudocount = 0.5) {
  logm <- log(mat + pseudocount)
  sweep(logm, 1, rowMeans(logm), "-")
}

pairwise_clr_cor <- function(ref_counts, qry_counts, method, scope) {
  common <- intersect(rownames(ref_counts), rownames(qry_counts))
  if (length(common) == 0L) {
    return(data.frame())
  }
  ref_counts <- ref_counts[common, , drop = FALSE]
  qry_counts <- qry_counts[common, colnames(ref_counts), drop = FALSE]
  rclr <- clr_rows(ref_counts)
  qclr <- clr_rows(qry_counts)
  corr <- vapply(seq_along(common), function(i) {
    if (stats::sd(rclr[i, ]) == 0 || stats::sd(qclr[i, ]) == 0) return(NA_real_)
    stats::cor(rclr[i, ], qclr[i, ], use = "pairwise.complete.obs")
  }, numeric(1))
  data.frame(
    method = method,
    scope = scope,
    repeat_group = common,
    correlation = corr,
    weight = rowSums(ref_counts) + rowSums(qry_counts),
    stringsAsFactors = FALSE
  )
}

summarise_correlations <- function(pairs) {
  if (nrow(pairs) == 0L) return(data.frame())
  rows <- lapply(split(pairs, paste(pairs$method, pairs$scope, sep = "||")), function(d) {
    ok <- is.finite(d$correlation)
    data.frame(
      method = d$method[1],
      scope = d$scope[1],
      n_pairs = sum(ok),
      mean_correlation = mean(d$correlation[ok], na.rm = TRUE),
      median_correlation = stats::median(d$correlation[ok], na.rm = TRUE),
      weighted_mean_correlation = if (any(ok)) {
        stats::weighted.mean(d$correlation[ok], d$weight[ok], na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

choose_repeat_group <- function(meta, candidate_cols, ref_batch, qry_batch) {
  for (cc in candidate_cols) {
    g <- as.character(meta[[cc]])
    g[!nzchar(g) | is.na(g)] <- NA_character_
    ok <- !is.na(g)
    if (!any(ok)) next
    tab <- table(group = g[ok], batch = meta$batch[ok])
    if (!all(c(ref_batch, qry_batch) %in% colnames(tab))) next
    repeated <- rownames(tab)[tab[, ref_batch] > 0 & tab[, qry_batch] > 0]
    repeated <- repeated[rowSums(tab[repeated, , drop = FALSE]) >= 20]
    if (length(repeated) >= 2L) {
      return(list(group = g, col = cc, mode = "metadata", usable = repeated))
    }
  }
  NULL
}

make_pseudo_repeat_groups <- function(meta, ref_batch, qry_batch, n_groups = 8L) {
  g <- rep(NA_character_, nrow(meta))
  for (bb in c(ref_batch, qry_batch)) {
    for (cl in sort(unique(meta$label))) {
      ix <- which(meta$batch == bb & meta$label == cl)
      if (!length(ix)) next
      ix <- ix[sample.int(length(ix))]
      g[ix] <- sprintf("pseudo_%02d", rep_len(seq_len(n_groups), length(ix)))
    }
  }
  g
}

run_abundance_validation <- function(ref, fit, ref_y, ref_group, qry_group, baseline_res) {
  n_nodes <- nrow(ref$codebook)
  node_levels <- as.character(seq_len(n_nodes))
  label_levels <- colnames(ref$label_prob)
  drop_group <- "__somalign_drop_nonrepeat__"
  qry_group_soft <- as.character(qry_group)
  qry_group_soft[is.na(qry_group_soft)] <- drop_group

  ref_node <- count_matrix(ref_group, ref$reference_units, node_levels)
  qry_node_hard <- count_matrix(qry_group, baseline_res$old_som_unit, node_levels)
  qry_node_soft <- as.matrix(somalign_soft_frequencies(
    fit, group = qry_group_soft, node_groups = seq_len(n_nodes), normalize = FALSE
  ))
  qry_node_soft <- qry_node_soft[rownames(qry_node_soft) != drop_group, , drop = FALSE]
  qry_node_soft <- qry_node_soft[, node_levels, drop = FALSE]

  ref_label <- count_matrix(ref_group, ref_y, label_levels)
  qry_label_hard <- count_matrix(qry_group, baseline_res$old_som_label, label_levels)
  qry_label_soft <- as.matrix(somalign_soft_frequencies(
    fit, group = qry_group_soft, normalize = FALSE
  ))
  qry_label_soft <- qry_label_soft[rownames(qry_label_soft) != drop_group, , drop = FALSE]
  qry_label_soft <- qry_label_soft[, label_levels, drop = FALSE]

  pairs <- rbind(
    pairwise_clr_cor(ref_node, qry_node_hard, "hard_direct", "som_unit"),
    pairwise_clr_cor(ref_node, qry_node_soft, "soft_knn", "som_unit"),
    pairwise_clr_cor(ref_label, qry_label_hard, "hard_direct", "lineage"),
    pairwise_clr_cor(ref_label, qry_label_soft, "soft_knn", "lineage")
  )
  list(pairs = pairs, summary = summarise_correlations(pairs))
}

run_fit_config <- function(query, reference, truth, approach, config, args,
                           classes, subset_name = "all", note = "",
                           evidence_class = "deployable") {
  stamp(sprintf("fit: %s / %s", approach, config))
  fit <- do.call(somalign_fit, c(list(query = query, reference = reference), args))
  pred <- labels_from_fit(fit)
  score <- score_prediction(
    approach, config,
    predicted = pred$transferred_label,
    truth = truth,
    accepted = pred$transferred_label_accepted,
    confidence = pred$transferred_label_confidence,
    margin = pred$transferred_label_margin,
    subset_name = subset_name,
    note = note,
    classes = classes,
    evidence_class = evidence_class
  )
  list(fit = fit, pred = pred, score = score)
}

append_score <- function(scores, score) {
  scores$summary <- rbind(scores$summary, score$summary)
  scores$per_class <- rbind(scores$per_class, score$per_class)
  scores
}

stamp("loading labelled pilot data")
df <- qs2::qs_read(PILOT, nthreads = 8)
required <- c(FEATURES, LABEL_COL, BATCH_COL)
missing <- setdiff(required, names(df))
if (length(missing)) {
  stop("Input data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
}
labels_all <- as.character(df[[LABEL_COL]])
batch_all <- as.character(df[[BATCH_COL]])
keep <- !is.na(labels_all) & labels_all != "" &
  !is.na(batch_all) & batch_all %in% c(REF_BATCH, QRY_BATCH) &
  stats::complete.cases(df[, FEATURES])

idx0 <- which(keep)
strat <- interaction(batch_all[idx0], labels_all[idx0], drop = TRUE)
per <- max(1L, floor(n_sub / nlevels(strat)))
sub_idx <- unlist(lapply(split(idx0, strat), function(ix) {
  if (length(ix) <= per) ix else sample(ix, per)
}), use.names = FALSE)
sub_idx <- sample(sub_idx)

candidate_group_cols <- intersect(
  c("donor_id", "donor", "patient_id", "participant_id", "subject_id",
    "sample_id", "sample", "Sample", "sample_name", "fcs_filename",
    "filename", "file_name", "file", "run_id"),
  names(df)
)
meta <- data.frame(
  row_index = sub_idx,
  label = labels_all[sub_idx],
  batch = batch_all[sub_idx],
  stringsAsFactors = FALSE
)
for (cc in candidate_group_cols) {
  meta[[cc]] <- as.character(df[[cc]][sub_idx])
}

X <- as.matrix(df[sub_idx, FEATURES])
storage.mode(X) <- "double"
colnames(X) <- FEATURES
rm(df)
gc()

class_levels <- sort(unique(meta$label))
grid <- kohonen::somgrid(grid_side, grid_side, "hexagonal")
ref_i <- which(meta$batch == REF_BATCH)
qry_i <- which(meta$batch == QRY_BATCH)
ref_X <- X[ref_i, , drop = FALSE]
qry_X <- X[qry_i, , drop = FALSE]
ref_y <- meta$label[ref_i]
qry_y <- meta$label[qry_i]

stamp(sprintf(
  "subsampled %d cells: %s=%d, %s=%d",
  nrow(X), REF_BATCH, length(ref_i), QRY_BATCH, length(qry_i)
))
print(sort(table(meta$batch, meta$label)))

repeat_choice <- choose_repeat_group(meta, candidate_group_cols, REF_BATCH, QRY_BATCH)
if (is.null(repeat_choice)) {
  repeat_group <- make_pseudo_repeat_groups(meta, REF_BATCH, QRY_BATCH, n_groups = 8L)
  repeat_col <- "pseudo_balanced_group"
  repeat_mode <- "no_true_repeats_pseudo_balanced_control"
  usable_repeat_groups <- sort(unique(repeat_group[!is.na(repeat_group)]))
} else {
  repeat_group <- repeat_choice$group
  repeat_col <- repeat_choice$col
  repeat_mode <- repeat_choice$mode
  usable_repeat_groups <- repeat_choice$usable
}
meta$repeat_group <- repeat_group
ref_group <- meta$repeat_group[ref_i]
qry_group <- meta$repeat_group[qry_i]
keep_ref_group <- ref_group %in% usable_repeat_groups
keep_qry_group <- qry_group %in% usable_repeat_groups
ref_group[!keep_ref_group] <- NA_character_
qry_group[!keep_qry_group] <- NA_character_

metadata <- data.frame(
  pilot_file = PILOT,
  n_subsample_requested = n_sub,
  n_subsample_actual = nrow(X),
  ref_batch = REF_BATCH,
  query_batch = QRY_BATCH,
  n_reference = nrow(ref_X),
  n_query = nrow(qry_X),
  grid_side = grid_side,
  n_nodes = grid_side * grid_side,
  rlen = rlen,
  seed = seed,
  label_col = LABEL_COL,
  batch_col = BATCH_COL,
  repeat_mode = repeat_mode,
  repeat_group_col = repeat_col,
  n_repeat_groups = length(usable_repeat_groups),
  stringsAsFactors = FALSE
)
write_table(metadata, "run_metadata.csv")

stamp("training global reference and query SOMs")
ref <- somalign_train_reference(ref_X, labels = ref_y, grid = grid, rlen = rlen)
qry <- somalign_query(qry_X, ref, grid = grid, rlen = rlen)

t_weights <- make_t_marker_weights(FEATURES)
qry_oracle <- qry
qry_oracle$label_prob <- label_prob_from_units(
  qry_y, qry$sample_unit, nrow(qry$codebook), colnames(ref$label_prob)
)

scores <- list(summary = data.frame(), per_class = data.frame())
fits <- list()
predictions <- list()

global_configs <- list(
  list(
    approach = "A0_baseline",
    config = "global_epsilon_0.1",
    args = list(epsilon = 0.1),
    query = qry,
    evidence_class = "deployable",
    note = "Current shipped operating point."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.02",
    args = list(epsilon = 0.02),
    query = qry,
    evidence_class = "deployable",
    note = "Sharper OT plan, included in targeted OtherT grid."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.05",
    args = list(epsilon = 0.05),
    query = qry,
    evidence_class = "deployable",
    note = "Sharper OT plan, included in targeted OtherT grid."
  ),
  list(
    approach = "A3_epsilon_0.2",
    config = "global_epsilon_0.2",
    args = list(epsilon = 0.2),
    query = qry,
    evidence_class = "deployable",
    note = "Smoother higher-coverage-quality operating point from prior CV."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.5",
    args = list(epsilon = 0.5),
    query = qry,
    evidence_class = "deployable",
    note = "High-epsilon coverage/precision tradeoff control."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.05_rho_2_2",
    args = list(epsilon = 0.05, rho_query = 2, rho_ref = 2),
    query = qry,
    evidence_class = "deployable",
    note = "Higher UOT mass-retention pressure at epsilon 0.05."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.1_rho_2_2",
    args = list(epsilon = 0.1, rho_query = 2, rho_ref = 2),
    query = qry,
    evidence_class = "deployable",
    note = "Higher UOT mass-retention pressure at epsilon 0.1."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.1_diagonal_boost_0.1",
    args = list(epsilon = 0.1, diagonal_boost = 0.1),
    query = qry,
    evidence_class = "deployable",
    note = "Identity-neighbour cost preference, included in targeted grid."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.1_diagonal_boost_0.2",
    args = list(epsilon = 0.1, diagonal_boost = 0.2),
    query = qry,
    evidence_class = "deployable",
    note = "Identity-neighbour cost preference, included in targeted grid."
  ),
  list(
    approach = "A2_targeted_grid",
    config = "global_epsilon_0.1_diagonal_boost_0.5",
    args = list(epsilon = 0.1, diagonal_boost = 0.5),
    query = qry,
    evidence_class = "deployable",
    note = "Stronger identity-neighbour cost preference."
  ),
  list(
    approach = "A4_marker_weighted_cost",
    config = "t_marker_weights_epsilon_0.05",
    args = list(epsilon = 0.05, feature_weights = t_weights),
    query = qry,
    evidence_class = "deployable",
    note = "Explicit T/NK boundary marker weighting plus epsilon 0.05."
  ),
  list(
    approach = "A4_marker_weighted_cost",
    config = "t_marker_weights_epsilon_0.1",
    args = list(epsilon = 0.1, feature_weights = t_weights),
    query = qry,
    evidence_class = "deployable",
    note = "Explicit T/NK boundary marker weighting."
  ),
  list(
    approach = "A4_marker_weighted_cost",
    config = "t_marker_weights_epsilon_0.2",
    args = list(epsilon = 0.2, feature_weights = t_weights),
    query = qry,
    evidence_class = "deployable",
    note = "Explicit T/NK boundary marker weighting plus epsilon 0.2."
  ),
  list(
    approach = "extra_oracle_label_guided",
    config = "oracle_query_label_guided_epsilon_0.1",
    args = list(epsilon = 0.1, label_guided = TRUE),
    query = qry_oracle,
    evidence_class = "oracle_upper_bound",
    note = "Upper-bound only: query SOM labels are built from held-out truth."
  )
)

for (cfg in global_configs) {
  out <- run_fit_config(cfg$query, ref, qry_y, cfg$approach, cfg$config,
                        cfg$args, class_levels, note = cfg$note,
                        evidence_class = cfg$evidence_class)
  key <- cfg$config
  fits[[key]] <- out$fit
  predictions[[key]] <- out$pred
  scores <- append_score(scores, out$score)
}

baseline_fit <- fits[["global_epsilon_0.1"]]
baseline_pred <- predictions[["global_epsilon_0.1"]]
baseline_res <- somalign_results(baseline_fit, include_correction = FALSE)

stamp("validating deployable OtherT rescue combinations")
rescue_sources <- intersect(c("global_epsilon_0.02", "t_marker_weights_epsilon_0.05"),
                            names(predictions))
rescue_margin_thresholds <- c(0, 0.20, 0.50, 0.80)
rescue_boundary_labels <- c("CD4T", "CD8T", "NK", "OtherT")
for (src in rescue_sources) {
  alt <- predictions[[src]]
  boundary_candidate <- !baseline_pred$transferred_label_accepted |
    baseline_pred$raw_transferred_label %in% rescue_boundary_labels |
    baseline_pred$transferred_label_second %in% rescue_boundary_labels
  for (thr in rescue_margin_thresholds) {
    rescue <- boundary_candidate &
      alt$transferred_label_accepted &
      alt$raw_transferred_label == "OtherT" &
      alt$transferred_label_margin >= thr
    combo_label <- baseline_pred$raw_transferred_label
    combo_accepted <- baseline_pred$transferred_label_accepted
    combo_conf <- baseline_pred$transferred_label_confidence
    combo_margin <- baseline_pred$transferred_label_margin
    combo_label[rescue] <- alt$raw_transferred_label[rescue]
    combo_accepted[rescue] <- TRUE
    combo_conf[rescue] <- alt$transferred_label_confidence[rescue]
    combo_margin[rescue] <- alt$transferred_label_margin[rescue]
    sc <- score_prediction(
      "A2_plus_A4_otherT_rescue",
      sprintf("baseline_plus_%s_otherT_margin_ge_%0.2f", src, thr),
      combo_label,
      qry_y,
      accepted = combo_accepted,
      confidence = combo_conf,
      margin = combo_margin,
      note = sprintf(
        "Deployable combination: baseline labels plus %d OtherT rescues from %s.",
        sum(rescue), src
      ),
      classes = class_levels
    )
    scores <- append_score(scores, sc)
  }
}

stamp("validating hard-vs-soft repeat abundance correlations")
abundance <- run_abundance_validation(
  baseline_fit$reference, baseline_fit, ref_y,
  ref_group = ref_group,
  qry_group = qry_group,
  baseline_res = baseline_res
)
if (nrow(abundance$pairs)) {
  abundance$pairs$repeat_mode <- repeat_mode
  abundance$pairs$repeat_group_col <- repeat_col
  abundance$pairs$evidence_class <- if (identical(repeat_mode, "metadata")) {
    "true_repeat_metadata"
  } else {
    "pseudo_repeat_control"
  }
}
if (nrow(abundance$summary)) {
  abundance$summary$repeat_mode <- repeat_mode
  abundance$summary$repeat_group_col <- repeat_col
  abundance$summary$evidence_class <- if (identical(repeat_mode, "metadata")) {
    "true_repeat_metadata"
  } else {
    "pseudo_repeat_control"
  }
}
write_table(abundance$pairs, "abundance_correlation_pairs.csv")
write_table(abundance$summary, "abundance_correlation_summary.csv")

stamp("validating margin and transport-entropy triage")
node_entropy <- baseline_fit$diagnostics$nodes$transport_entropy[baseline_pred$query_som_unit]
triage_thresholds <- c(0.05, 0.10, 0.20, 0.30)
for (thr in triage_thresholds) {
  accepted <- baseline_pred$transferred_label_accepted &
    baseline_pred$transferred_label_margin >= thr
  sc <- score_prediction(
    "A6_margin_entropy_triage",
    sprintf("baseline_margin_ge_%0.2f", thr),
    baseline_pred$raw_transferred_label,
    qry_y,
    accepted = accepted,
    confidence = baseline_pred$transferred_label_confidence,
    margin = baseline_pred$transferred_label_margin,
    note = "Selective acceptance; abstentions count against full recall.",
    classes = class_levels
  )
  scores <- append_score(scores, sc)
}
finite_entropy <- is.finite(node_entropy) & baseline_pred$transferred_label_accepted
if (sum(finite_entropy) > 0L) {
  entropy_thr <- stats::quantile(node_entropy[finite_entropy], probs = c(0.5, 0.75), na.rm = TRUE)
  for (qname in names(entropy_thr)) {
    accepted <- baseline_pred$transferred_label_accepted &
      baseline_pred$transferred_label_margin >= 0.10 &
      node_entropy <= as.numeric(entropy_thr[[qname]])
    sc <- score_prediction(
      "A6_margin_entropy_triage",
      sprintf("baseline_margin_ge_0.10_entropy_le_%s", qname),
      baseline_pred$raw_transferred_label,
      qry_y,
      accepted = accepted,
      confidence = baseline_pred$transferred_label_confidence,
      margin = baseline_pred$transferred_label_margin,
      note = "Margin plus low transport-entropy selective acceptance.",
      classes = class_levels
    )
    scores <- append_score(scores, sc)
  }
}

stamp("validating T/NK-local projection and rerouting")
local_labels <- intersect(T_BOUNDARY_LABELS, class_levels)
ref_local_i <- ref_y %in% local_labels
qry_truth_local_i <- qry_y %in% local_labels
router_diagnostics <- data.frame()
if (sum(ref_local_i) > 100L && sum(qry_truth_local_i) > 100L) {
  local_side <- max(6L, min(grid_side, 12L))
  local_grid <- kohonen::somgrid(local_side, local_side, "hexagonal")
  local_ref <- somalign_train_reference(
    ref_X[ref_local_i, T_FEATURES, drop = FALSE],
    labels = ref_y[ref_local_i],
    grid = local_grid,
    rlen = rlen
  )
  local_qry_truth <- somalign_query(
    qry_X[qry_truth_local_i, T_FEATURES, drop = FALSE],
    local_ref,
    grid = local_grid,
    rlen = rlen
  )
  for (eps in c(0.1, 0.2)) {
    fit_local <- somalign_fit(local_qry_truth, local_ref, epsilon = eps)
    pred_local <- labels_from_fit(fit_local)
    sc <- score_prediction(
      "A5_tnk_local_projection",
      sprintf("truth_tnk_local_epsilon_%0.1f", eps),
      pred_local$transferred_label,
      qry_y[qry_truth_local_i],
      accepted = pred_local$transferred_label_accepted,
      confidence = pred_local$transferred_label_confidence,
      margin = pred_local$transferred_label_margin,
      subset_name = "truth_TNK",
      note = "Oracle subset: true T/NK cells are selected before local projection.",
      classes = local_labels,
      evidence_class = "oracle_subset"
    )
    scores <- append_score(scores, sc)
  }

  route <- baseline_pred$raw_transferred_label %in% local_labels |
    baseline_pred$transferred_label_second %in% local_labels
  if (sum(route) > 100L) {
    router_diagnostics <- data.frame(
      config = "global_top_or_second_TNK",
      n_query = length(route),
      n_routed = sum(route),
      route_fraction = mean(route),
      true_tnk_route_recall = mean(route[qry_truth_local_i]),
      false_route_rate_non_tnk = mean(route[!qry_truth_local_i]),
      route_precision_tnk = mean(qry_y[route] %in% local_labels),
      otherT_route_recall = mean(route[qry_y == "OtherT"]),
      stringsAsFactors = FALSE
    )
    local_qry_route <- somalign_query(
      qry_X[route, T_FEATURES, drop = FALSE],
      local_ref,
      grid = local_grid,
      rlen = rlen
    )
    fit_route <- somalign_fit(local_qry_route, local_ref, epsilon = 0.1)
    pred_route <- labels_from_fit(fit_route)
    combo_label <- baseline_pred$raw_transferred_label
    combo_accepted <- baseline_pred$transferred_label_accepted
    combo_conf <- baseline_pred$transferred_label_confidence
    combo_margin <- baseline_pred$transferred_label_margin
    combo_label[route] <- pred_route$raw_transferred_label
    combo_accepted[route] <- pred_route$transferred_label_accepted
    combo_conf[route] <- pred_route$transferred_label_confidence
    combo_margin[route] <- pred_route$transferred_label_margin
    sc_all <- score_prediction(
      "A5_tnk_local_projection",
      "global_with_tnk_local_router",
      combo_label,
      qry_y,
      accepted = combo_accepted,
      confidence = combo_conf,
      margin = combo_margin,
      note = "Production-like reroute when global top or second label is T/NK.",
      classes = class_levels
    )
    scores <- append_score(scores, sc_all)
    sc_t <- score_prediction(
      "A5_tnk_local_projection",
      "global_with_tnk_local_router",
      combo_label[qry_truth_local_i],
      qry_y[qry_truth_local_i],
      accepted = combo_accepted[qry_truth_local_i],
      confidence = combo_conf[qry_truth_local_i],
      margin = combo_margin[qry_truth_local_i],
      subset_name = "truth_TNK",
      note = "Deployable router, evaluated only on true T/NK cells.",
      classes = local_labels
    )
    scores <- append_score(scores, sc_t)
    sc_combo_margin <- score_prediction(
      "A5_plus_A6_local_margin",
      "global_tnk_router_margin_ge_0.05",
      combo_label,
      qry_y,
      accepted = combo_accepted & combo_margin >= 0.05,
      confidence = combo_conf,
      margin = combo_margin,
      note = "Local rerouting followed by light margin triage.",
      classes = class_levels
    )
    scores <- append_score(scores, sc_combo_margin)
  }
}
write_table(router_diagnostics, "router_diagnostics.csv")

stamp("validating anchor tuning and anchor-derived marker weights")
n_anchor <- min(4000L, max(100L, floor(0.2 * nrow(qry_X))))
anchors <- make_anchors(ref_X, ref_y, qry_X, qry_y, n_anchor)
anchor_eval_mask <- rep(TRUE, nrow(qry_X))
anchor_eval_mask[anchors$query_anchor_idx] <- FALSE
random_anchors <- make_random_anchors(ref_X, qry_X, n_anchor)
random_anchor_eval_mask <- rep(TRUE, nrow(qry_X))
random_anchor_eval_mask[random_anchors$query_anchor_idx] <- FALSE
rho_grid <- c(0, 1, 5, 20, 100, 500)
anchor_proxy_grid <- somalign_anchor_benefit(
  qry, ref, qry_y,
  anchor_old = anchors$anchor_old,
  anchor_new = anchors$anchor_new,
  rho_grid = rho_grid,
  metric = "macro_f1",
  eval_mask = anchor_eval_mask,
  epsilon = 0.1
)
anchor_random_grid <- somalign_anchor_benefit(
  qry, ref, qry_y,
  anchor_old = random_anchors$anchor_old,
  anchor_new = random_anchors$anchor_new,
  rho_grid = rho_grid,
  metric = "macro_f1",
  eval_mask = random_anchor_eval_mask,
  epsilon = 0.1
)
anchor_table <- rbind(
  transform(anchor_proxy_grid$grid,
            anchor_type = "proxy_same_lineage",
            evidence_class = "proxy_oracle_anchor"),
  transform(anchor_random_grid$grid,
            anchor_type = "random_mispaired_control",
            evidence_class = "negative_control")
)
anchor_table$metric <- "macro_f1"
anchor_table$n_anchor <- nrow(anchors$anchor_old)
write_table(anchor_table, "anchor_grid.csv")

sc_anchor_baseline <- score_prediction(
  "A7_anchor_tuning",
  "baseline_eval_nonanchor_subset",
  baseline_pred$raw_transferred_label[anchor_eval_mask],
  qry_y[anchor_eval_mask],
  accepted = baseline_pred$transferred_label_accepted[anchor_eval_mask],
  confidence = baseline_pred$transferred_label_confidence[anchor_eval_mask],
  margin = baseline_pred$transferred_label_margin[anchor_eval_mask],
  subset_name = "non_anchor_eval",
  note = "Baseline restricted to the held-out non-anchor subset.",
  classes = class_levels
)
scores <- append_score(scores, sc_anchor_baseline)

for (rho in unique(c(1, anchor_proxy_grid$best$rho_anchor))) {
  if (!is.finite(rho) || rho <= 0) next
  stamp(sprintf("anchored fit at rho_anchor=%s", rho))
  fit_anchor <- somalign_fit_anchored(
    qry, ref,
    anchor_old = anchors$anchor_old,
    anchor_new = anchors$anchor_new,
    rho_anchor = rho,
    correction = "cost_bonus",
    epsilon = 0.1
  )
  pred_anchor <- labels_from_fit(fit_anchor)
  sc <- score_prediction(
    "A7_anchor_tuning",
    sprintf("rho_anchor_%s", rho),
    pred_anchor$raw_transferred_label[anchor_eval_mask],
    qry_y[anchor_eval_mask],
    accepted = pred_anchor$transferred_label_accepted[anchor_eval_mask],
    confidence = pred_anchor$transferred_label_confidence[anchor_eval_mask],
    margin = pred_anchor$transferred_label_margin[anchor_eval_mask],
    subset_name = "non_anchor_eval",
    note = "Proxy same-lineage anchors; compare to non-anchor baseline row.",
    classes = class_levels,
    evidence_class = "proxy_oracle_anchor"
  )
  scores <- append_score(scores, sc)
}

fit_anchor_fw <- somalign_fit_anchored(
  qry, ref,
  anchor_old = anchors$anchor_old,
  anchor_new = anchors$anchor_new,
  rho_anchor = 0,
  correction = "cost_bonus",
  epsilon = 0.1,
  feature_weights = "anchor"
)
pred_anchor_fw <- labels_from_fit(fit_anchor_fw)
sc_anchor_fw <- score_prediction(
  "A4_marker_weighted_cost",
  "anchor_feature_weights_rho_0",
  pred_anchor_fw$raw_transferred_label[anchor_eval_mask],
  qry_y[anchor_eval_mask],
  accepted = pred_anchor_fw$transferred_label_accepted[anchor_eval_mask],
  confidence = pred_anchor_fw$transferred_label_confidence[anchor_eval_mask],
  margin = pred_anchor_fw$transferred_label_margin[anchor_eval_mask],
  subset_name = "non_anchor_eval",
  note = "Anchor-derived feature weights without anchor cost bonus.",
  classes = class_levels,
  evidence_class = "proxy_oracle_anchor"
)
scores <- append_score(scores, sc_anchor_fw)

targeted_all <- scores$summary[
  scores$summary$subset == "all" &
    scores$summary$approach %in% c("A0_baseline", "A2_targeted_grid",
                                   "A3_epsilon_0.2", "A4_marker_weighted_cost",
                                   "A2_plus_A4_otherT_rescue",
                                   "extra_oracle_label_guided"),
  ,
  drop = FALSE
]
targeted_all$otherT_weighted_score <- targeted_all$otherT_f1 * targeted_all$coverage
targeted_all <- targeted_all[order(-targeted_all$otherT_weighted_score,
                                   -targeted_all$otherT_f1,
                                   -targeted_all$full_macro_f1,
                                   -targeted_all$accuracy_all), ]
targeted_deployable <- targeted_all[targeted_all$is_deployable, , drop = FALSE]

write_table(scores$summary, "label_summary.csv")
write_table(scores$per_class, "per_class_metrics.csv")
write_table(targeted_deployable, "targeted_grid_ranked.csv")
write_table(targeted_all, "targeted_grid_ranked_all.csv")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))

saveRDS(
  list(
    metadata = metadata,
    label_summary = scores$summary,
    per_class_metrics = scores$per_class,
    targeted_grid_ranked = targeted_deployable,
    targeted_grid_ranked_all = targeted_all,
    abundance_pairs = abundance$pairs,
    abundance_summary = abundance$summary,
    router_diagnostics = router_diagnostics,
    anchor_table = anchor_table,
    anchor_proxy_grid = anchor_proxy_grid,
    anchor_random_grid = anchor_random_grid,
    class_counts = list(reference = sort(table(ref_y), decreasing = TRUE),
                        query = sort(table(qry_y), decreasing = TRUE)),
    feature_weights = list(t_marker_weights = t_weights,
                           anchor_feature_weights = fit_anchor_fw$anchors$feature_weights)
  ),
  OUT_RDS
)

stamp(sprintf("saved tables and RDS bundle -> %s", OUT_DIR))
stamp("top deployable targeted-grid rows")
print(utils::head(targeted_deployable[, c("approach", "config", "coverage", "accuracy_all",
                                          "full_macro_f1", "otherT_precision",
                                          "otherT_recall", "otherT_f1",
                                          "otherT_weighted_score")], 10))
if (nrow(abundance$summary)) {
  stamp("abundance correlation summary")
  print(abundance$summary)
}
stamp("complete")
