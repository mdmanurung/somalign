## somalign — Label transfer quality and full diagnostics
##
## Extends benchmark_real_data.R by adding biologically meaningful reference
## labels and evaluating somalign label transfer against independently gated
## batch2 populations.
##
## Strategy:
##   1. Gate batch1 cells into major immune populations using lineage markers.
##   2. Train reference SOM with those per-cell labels.
##   3. Run somalign_fit(); extract transferred_label, confidence, accepted.
##   4. Apply the same thresholds to batch2 ("ground truth" batch2 labels).
##   5. Compare somalign transferred_label to ground truth → accuracy / recall.
##   6. Print every diagnostic field from somalign_diagnostics() in full.
##
## Run with:
##   Rscript benchmarks/benchmark_label_transfer.R 2>&1 | tee benchmarks/label_transfer_output.log

suppressPackageStartupMessages({
  library(devtools)
  load_all(quiet = TRUE)
  library(kohonen)
  library(flowCore)
})

source("benchmarks/helpers.R")   # sep, fmt_pct, js_div, node_dist, ks_per_marker

## ── Macro-F1 (balanced metric, insensitive to class imbalance) ───────────────
macro_f1 <- function(true_lab, pred_lab) {
  labs <- union(unique(true_lab), unique(pred_lab))
  f1s <- vapply(labs, function(cl) {
    tp <- sum(true_lab == cl & pred_lab == cl, na.rm = TRUE)
    fp <- sum(true_lab != cl & pred_lab == cl, na.rm = TRUE)
    fn <- sum(true_lab == cl & pred_lab != cl, na.rm = TRUE)
    p  <- if (tp + fp == 0) 0 else tp / (tp + fp)
    r  <- if (tp + fn == 0) 0 else tp / (tp + fn)
    if (p + r == 0) 0 else 2 * p * r / (p + r)
  }, numeric(1))
  mean(f1s)
}

## ── Confusion / recall table ──────────────────────────────────────────────────
conf_table <- function(true_lab, pred_lab, title = "") {
  lvls  <- sort(unique(c(true_lab, pred_lab)))
  mat   <- table(True = factor(true_lab, lvls), Pred = factor(pred_lab, lvls))
  cat(title, "\n")
  print(mat)
  per_class_recall <- diag(mat) / rowSums(mat)
  cat("Per-class recall:\n")
  for (cl in names(per_class_recall))
    cat(sprintf("  %-12s %.1f%%\n", cl, 100 * per_class_recall[cl]))
  cat(sprintf("Overall accuracy: %.1f%%\n",
              100 * sum(diag(mat)) / sum(mat)))
  invisible(mat)
}

## ─────────────────────────────────────────────────────────────────────────────
## 1.  Nuñez 2023 — Full-Spectrum Flow
## ─────────────────────────────────────────────────────────────────────────────

sep("Nuñez 2023 — loading and gating")

fcs1 <- "data/Nunez2023/Nunez_PBMCs_batch1.fcs"
fcs2 <- "data/Nunez2023/Nunez_PBMCs_batch2.fcs"

if (!file.exists(fcs1)) {
  cat("SKIP: Nuñez FCS files not found.\n")
} else {

  ## ── 1a. Load and transform ──
  cofactor <- 2000
  load_nunez <- function(path) {
    ff  <- read.FCS(path, truncate_max_range = FALSE)
    pd  <- pData(parameters(ff))
    # map desc → channel name for the markers we need
    ch  <- setNames(pd$name, pd$desc)
    mat <- asinh(exprs(ff)[, grep("^FJComp-", pd$name, value = TRUE)] / cofactor)
    colnames(mat) <- pd$desc[match(grep("^FJComp-", pd$name, value = TRUE), pd$name)]
    colnames(mat)[is.na(colnames(mat)) | colnames(mat) == "-"] <- "unnamed"
    list(mat = mat, channel_map = ch)
  }

  b1 <- load_nunez(fcs1)$mat
  b2 <- load_nunez(fcs2)$mat
  cat(sprintf("Loaded: batch1=%d cells  batch2=%d cells  %d markers\n",
              nrow(b1), nrow(b2), ncol(b1)))

  ## ── 1b. Hierarchical gate → per-cell labels ──
  ## Thresholds are on arcsinh(x / 2000) scale.
  ## Values chosen at the natural valley of each bimodal marker.
  gate_nunez <- function(mat) {
    cd3  <- mat[, "CD3"]
    cd19 <- mat[, "CD19"]
    cd14 <- mat[, "CD14"]
    cd56 <- mat[, "CD56"]

    ## Thresholds on arcsinh(x / 2000) scale.
    ## Full-spectrum flow has high spectral spillover; thresholds set at the valley
    ## of the bimodal distributions (verified on batch1 quantiles).
    ## CD4/CD8 are too ambiguous in this panel for reliable sub-gating;
    ## only primary lineage gate is used.
    thr_cd3  <- 1.5   # clear valley between neg (~0) and T-cell bright (~3)
    thr_cd19 <- 2.0   # valley between background and B cells for CD3-neg cells
    thr_cd14 <- 3.0   # valley between background and monocytes for CD3-neg cells
    thr_cd56 <- 1.5   # valley between background and NK for CD3-neg cells

    t_pos  <- cd3 >= thr_cd3
    b_pos  <- !t_pos

    label <- rep("Other", nrow(mat))
    label[t_pos]                                      <- "T"
    label[b_pos & cd19 >= thr_cd19]                   <- "B"
    label[b_pos & cd14 >= thr_cd14]                   <- "Mono"
    label[b_pos & cd56 >= thr_cd56]                   <- "NK"
    label
  }

  labels_b1 <- gate_nunez(b1)
  labels_b2 <- gate_nunez(b2)   # independent "ground truth" for batch2

  cat("\nBatch1 gate counts:\n")
  print(sort(table(labels_b1), decreasing = TRUE))
  cat("\nBatch2 gate counts (ground truth):\n")
  print(sort(table(labels_b2), decreasing = TRUE))

  ## Drop marker columns not in the reference SOM (LD viability already excluded)
  markers_fl <- setdiff(colnames(b1), c("LD", "unnamed"))
  b1_m <- b1[, markers_fl]
  b2_m <- b2[, markers_fl]

  ## ── 1b2. Per-marker batch shift ──
  sep("Nuñez 2023 — per-marker batch shift (batch2 vs batch1)")
  n_feat_fl <- ncol(b1_m)
  mu1_fl <- colMeans(b1_m);  mu2_fl <- colMeans(b2_m)
  sd1_fl <- apply(b1_m, 2, sd);  sd2_fl <- apply(b2_m, 2, sd)
  pool_sd_fl <- sqrt((sd1_fl^2 + sd2_fl^2) / 2)
  z_fl <- (mu2_fl - mu1_fl) / pool_sd_fl
  cat(sprintf("  n_markers=%d  |z| mean=%.3f  median=%.3f  max=%.3f\n",
              n_feat_fl, mean(abs(z_fl)), median(abs(z_fl)), max(abs(z_fl))))
  cat(sprintf("  |z|>0.5: %d  |z|>1.0: %d  |z|>2.0: %d\n",
              sum(abs(z_fl) > 0.5), sum(abs(z_fl) > 1.0), sum(abs(z_fl) > 2.0)))
  ord_fl <- order(abs(z_fl), decreasing = TRUE)
  cat("  Top 5 shifted markers:\n")
  for (i in seq_len(min(5L, n_feat_fl)))
    cat(sprintf("    %-22s  z = %+.3f\n", names(z_fl)[ord_fl[i]], z_fl[ord_fl[i]]))
  cat("  (small |z| = weak batch effect; OT correction likely adds noise in this regime)\n")

  ## ── 1c. Train labeled reference SOM ──
  sep("Nuñez 2023 — training labeled reference")
  set.seed(42)
  ref_fl <- somalign_train_reference(
    b1_m,
    labels = labels_b1,
    grid   = somgrid(10, 10, "hexagonal"),
    rlen   = 100
  )
  cat("Reference SOM: 100 nodes,", length(unique(labels_b1)), "label classes\n")

  ## ── 1d. Query + fit ──
  sep("Nuñez 2023 — somalign fit")
  set.seed(123)
  query_fl <- somalign_query(b2_m, ref_fl,
                              grid = somgrid(10, 10, "hexagonal"), rlen = 100)
  fit_fl   <- somalign_fit(query_fl, ref_fl,
                            epsilon = 0.1, rho_query = 2, rho_ref = 2)
  res_fl   <- somalign_results(fit_fl)

  ## ── 1e. Full diagnostics ──
  sep("Nuñez 2023 — full somalign_diagnostics()")
  diag_fl <- somalign_diagnostics(fit_fl)

  cat("\n### solver ###\n")
  cat(sprintf("  used:                    %s\n",   diag_fl$solver$used))
  cat(sprintf("  converged:               %s\n",   diag_fl$solver$converged))
  cat(sprintf("  iterations:              %d\n",   diag_fl$solver$iterations))
  cat(sprintf("  final_delta:             %.3e\n", diag_fl$solver$final_delta))
  cat(sprintf("  rel_marginal_row_error:  %.4f\n", diag_fl$solver$rel_marginal_row_error))
  cat(sprintf("  rel_marginal_col_error:  %.4f\n", diag_fl$solver$rel_marginal_col_error))
  cat(sprintf("  cost_scale:              %.4f\n", diag_fl$solver$cost_scale))

  cat("\n### ot ###\n")
  cat(sprintf("  transport_mass:       %.4f\n", diag_fl$ot$transport_mass))
  cat(sprintf("  max_row_mass_error:   %.4f\n", diag_fl$ot$max_row_mass_error))
  mf <- diag_fl$ot$match_fraction
  cat(sprintf("  match_fraction  [min / q25 / med / q75 / max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(mf), quantile(mf,.25), median(mf), quantile(mf,.75), max(mf)))
  mmr <- diag_fl$ot$match_mass_ratio
  cat(sprintf("  match_mass_ratio [min / q25 / med / q75 / max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(mmr), quantile(mmr,.25), median(mmr), quantile(mmr,.75), max(mmr)))
  cat(sprintf("  nodes with match_fraction < 0.5:  %d / %d\n",
              sum(mf < 0.5), length(mf)))
  cat(sprintf("  nodes with match_fraction < 0.25: %d / %d\n",
              sum(mf < 0.25), length(mf)))

  cat("\n### nodes (per-node correction) ###\n")
  cn <- diag_fl$nodes$correction_norm
  cat(sprintf("  correction_norm [min / q25 / med / q75 / max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(cn, na.rm=T), quantile(cn,.25,na.rm=T), median(cn,na.rm=T),
              quantile(cn,.75,na.rm=T), max(cn,na.rm=T)))
  cat(sprintf("  per-feature norm (/ sqrt(%d) markers):  med=%.3f  max=%.3f\n",
              n_feat_fl, median(cn, na.rm=TRUE)/sqrt(n_feat_fl),
              max(cn, na.rm=TRUE)/sqrt(n_feat_fl)))
  cat("  (values ~0.5–1 = plausible per-marker batch offset; >2 = gross shift)\n")
  if (all(mf >= 1.0 - 1e-9)) {
    cat("\n  DIAGNOSTIC FLAG: match_fraction = 1.0 for ALL nodes.\n")
    cat("  Cause: all match_mass_ratio > 1 → OT over-transports to every query node.\n")
    cat("  Effect: match_fraction acceptance gate is INOPERATIVE at these parameters.\n")
    cat("  Remedy: raise rho_query/rho_ref (e.g. rho=10) to tighten marginals.\n")
  }

  cat("\n### projection ###\n")
  cat(sprintf("  outside_direct_fraction:    %.4f  (%s of cells)\n",
              diag_fl$projection$outside_direct_fraction,
              fmt_pct(diag_fl$projection$outside_direct_fraction)))
  cat(sprintf("  outside_corrected_fraction: %.4f  (%s of cells)\n",
              diag_fl$projection$outside_corrected_fraction,
              fmt_pct(diag_fl$projection$outside_corrected_fraction)))

  ## Reference node purity diagnostic
  cat("\n### reference node purity (max label_prob per node) ###\n")
  node_purity_fl <- apply(ref_fl$label_prob, 1, max)
  cat(sprintf("  [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(node_purity_fl), quantile(node_purity_fl,.25), median(node_purity_fl),
              quantile(node_purity_fl,.75), max(node_purity_fl)))
  cat(sprintf("  nodes with purity >= 0.6: %d / %d\n", sum(node_purity_fl >= 0.6), length(node_purity_fl)))
  cat(sprintf("  nodes with purity >= 0.4: %d / %d\n", sum(node_purity_fl >= 0.4), length(node_purity_fl)))

  ## Per-node label confidence (transferred, before cell mapping)
  tl_node_conf_fl <- fit_fl$label_transfer$confidence
  cat(sprintf("\n### per-node transferred label confidence ###\n"))
  cat(sprintf("  [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(tl_node_conf_fl, na.rm=T), quantile(tl_node_conf_fl,.25,na.rm=T),
              median(tl_node_conf_fl,na.rm=T), quantile(tl_node_conf_fl,.75,na.rm=T),
              max(tl_node_conf_fl, na.rm=T)))
  cat(sprintf("  nodes with confidence >= 0.6: %d / %d\n",
              sum(tl_node_conf_fl >= 0.6, na.rm=T), sum(is.finite(tl_node_conf_fl))))
  cat(sprintf("  nodes with confidence >= 0.4: %d / %d\n",
              sum(tl_node_conf_fl >= 0.4, na.rm=T), sum(is.finite(tl_node_conf_fl))))

  ## ── 1f. Cell-level result columns ──
  sep("Nuñez 2023 — cell-level result summary")

  cat("\nfinal_status breakdown:\n")
  print(sort(table(res_fl$final_status), decreasing = TRUE))

  cat(sprintf("\nTransferred label: accepted=%s  rejected=%s  (of %d cells)\n",
              fmt_pct(mean(res_fl$transferred_label_accepted, na.rm=TRUE)),
              fmt_pct(mean(!res_fl$transferred_label_accepted, na.rm=TRUE)),
              nrow(res_fl)))

  conf_df <- data.frame(
    accepted  = res_fl$transferred_label_accepted,
    tl_conf   = res_fl$transferred_label_confidence,
    mf_cell   = mf[res_fl$query_som_unit]   # match_fraction of cell's query node
  )
  cat(sprintf("Mean confidence — accepted:  %.3f\n",
              mean(conf_df$tl_conf[conf_df$accepted], na.rm=TRUE)))
  cat(sprintf("Mean confidence — rejected:  %.3f\n",
              mean(conf_df$tl_conf[!conf_df$accepted], na.rm=TRUE)))
  cat(sprintf("Mean node match_fraction — accepted: %.3f\n",
              mean(conf_df$mf_cell[conf_df$accepted], na.rm=TRUE)))
  cat(sprintf("Mean node match_fraction — rejected: %.3f\n",
              mean(conf_df$mf_cell[!conf_df$accepted], na.rm=TRUE)))

  cat("\ncorrection_norm distribution (per cell):\n")
  cn_cell <- res_fl$correction_norm[is.finite(res_fl$correction_norm)]
  cat(sprintf("  [min / q25 / med / q75 / max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(cn_cell), quantile(cn_cell,.25), median(cn_cell),
              quantile(cn_cell,.75), max(cn_cell)))

  cat("\nold_som_distance vs distance_threshold:\n")
  od <- res_fl$old_som_distance
  ot <- res_fl$old_som_distance_threshold
  cat(sprintf("  Cells within threshold: %s\n", fmt_pct(mean(od <= ot, na.rm=TRUE))))
  cat(sprintf("  Cells outside threshold: %s\n", fmt_pct(mean(od > ot, na.rm=TRUE))))

  ## ── 1g. Label transfer quality vs ground truth ──
  sep("Nuñez 2023 — label transfer accuracy vs independent batch2 gating")

  ## Restrict to accepted labels only
  accepted_idx <- which(res_fl$transferred_label_accepted)
  cat(sprintf("%d / %d cells have accepted transferred labels\n",
              length(accepted_idx), nrow(res_fl)))

  tl  <- res_fl$transferred_label[accepted_idx]
  gt  <- labels_b2[accepted_idx]
  old <- res_fl$old_som_label[accepted_idx]

  cat("\n-- transferred_label vs batch2 ground truth (accepted cells only) --\n")
  conf_table(gt, tl, "Confusion matrix (rows=GT, cols=transferred_label):")

  cat("\n-- old_som_label (direct NN) vs batch2 ground truth (accepted cells) --\n")
  conf_table(gt, old, "Confusion matrix (rows=GT, cols=old_som_label):")

  ## Agreement between old_som_label and transferred_label
  agree <- mean(tl == old)
  cat(sprintf("\nAgreement between old_som_label and transferred_label: %.1f%%\n",
              100 * agree))

  ## ── 1g1. Direct NN on ALL cells (no OT gate, parameter-free) ──
  cat(sprintf("\n-- old_som_label (direct NN) vs batch2 GT — ALL %d cells --\n", nrow(res_fl)))
  conf_table(labels_b2, res_fl$old_som_label, "")

  ## ── 1g_acc. Label transfer acceptance/rejection breakdown by true class ──
  cat("\n-- label transfer acceptance by true class --\n")
  cat("  (rejected = confidence guardrail triggered, NOT mislabeled)\n")
  cat(sprintf("  %-10s  %7s  %7s  %7s  %s\n",
              "class", "n_total", "n_acc", "acc_%", "precision_if_acc"))
  for (cl in sort(unique(labels_b2))) {
    mask_cl  <- labels_b2 == cl
    n_cl     <- sum(mask_cl)
    mask_acc <- mask_cl & res_fl$transferred_label_accepted
    n_acc    <- sum(mask_acc, na.rm = TRUE)
    if (n_acc > 0) {
      prec <- mean(res_fl$transferred_label[mask_acc] == cl, na.rm = TRUE)
      cat(sprintf("  %-10s  %7d  %7d  %6.1f%%  %.1f%%\n",
                  cl, n_cl, n_acc, 100 * n_acc / n_cl, 100 * prec))
    } else {
      cat(sprintf("  %-10s  %7d  %7d  %6.1f%%  [none accepted]\n",
                  cl, n_cl, n_acc, 100 * n_acc / n_cl))
    }
  }

  ## ── 1g2. Relaxed confidence threshold ──
  sep("Nuñez 2023 — label transfer [confidence_threshold=0.4]")
  cat("NOTE: relaxing the acceptance threshold from 0.6 to 0.4 shows what IS transferred\n")
  cat("  for borderline nodes (same OT plan, only the acceptance threshold differs).\n\n")
  fit_fl_low <- somalign_fit(query_fl, ref_fl,
                              epsilon = 0.1, rho_query = 2, rho_ref = 2,
                              confidence_threshold = 0.4)
  res_fl_low <- somalign_results(fit_fl_low)

  cat(sprintf("Transferred label accepted: %s  rejected: %s\n",
              fmt_pct(mean(res_fl_low$transferred_label_accepted, na.rm=TRUE)),
              fmt_pct(mean(!res_fl_low$transferred_label_accepted, na.rm=TRUE))))

  acc_fl_low <- which(res_fl_low$transferred_label_accepted)
  if (length(acc_fl_low) > 0) {
    tl_fl_low  <- res_fl_low$transferred_label[acc_fl_low]
    gt_fl_low  <- labels_b2[acc_fl_low]
    old_fl_low <- res_fl_low$old_som_label[acc_fl_low]
    cat(sprintf("Confidence (accepted) — mean: %.3f  median: %.3f\n",
                mean(res_fl_low$transferred_label_confidence[acc_fl_low], na.rm=TRUE),
                median(res_fl_low$transferred_label_confidence[acc_fl_low], na.rm=TRUE)))
    cat("\n-- confidence_threshold=0.4: transferred_label vs batch2 GT --\n")
    conf_table(gt_fl_low, tl_fl_low, "")
    cat("\n-- confidence_threshold=0.4: old_som_label vs batch2 GT --\n")
    conf_table(gt_fl_low, old_fl_low, "")
    cat(sprintf("Agreement old_som_label ↔ transferred_label (thr=0.4): %.1f%%\n",
                100 * mean(tl_fl_low == old_fl_low)))
  } else {
    cat("No cells accepted at threshold=0.4 either.\n")
  }

  ## ── 1i. Sensitivity sweep: epsilon × rho ──
  sep("Nuñez 2023 — sensitivity sweep (epsilon x rho)")
  cat("Scores on ACCEPTED cells only (same set for OT and NN — like-for-like).\n")
  cat(sprintf("  %-7s %-5s  %-6s  %-9s  %-9s  %-9s  %-9s  %-12s  %-10s\n",
              "epsilon", "rho", "acc_%", "OT_acc%", "NN_acc%",
              "OT_F1%", "NN_F1%", "trans_mass", "out_corr%"))
  for (eps_sw in c(0.05, 0.1, 0.5)) {
    for (rho_sw in c(2, 10)) {
      f_sw <- somalign_fit(query_fl, ref_fl,
                           epsilon = eps_sw, rho_query = rho_sw, rho_ref = rho_sw)
      r_sw <- somalign_results(f_sw)
      d_sw <- somalign_diagnostics(f_sw)
      acc_frac_sw  <- mean(r_sw$transferred_label_accepted, na.rm = TRUE)
      acc_cells_sw <- which(r_sw$transferred_label_accepted)
      if (length(acc_cells_sw) > 0) {
        gt_sw  <- labels_b2[acc_cells_sw]
        ot_lab <- r_sw$transferred_label[acc_cells_sw]
        nn_lab <- r_sw$old_som_label[acc_cells_sw]
        ot_acc_str <- sprintf("%.1f", 100 * mean(ot_lab == gt_sw, na.rm = TRUE))
        nn_acc_str <- sprintf("%.1f", 100 * mean(nn_lab == gt_sw, na.rm = TRUE))
        ot_f1_str  <- sprintf("%.1f", 100 * macro_f1(gt_sw, ot_lab))
        nn_f1_str  <- sprintf("%.1f", 100 * macro_f1(gt_sw, nn_lab))
      } else {
        ot_acc_str <- nn_acc_str <- ot_f1_str <- nn_f1_str <- "  —  "
      }
      cat(sprintf("  %-7.2f %-5g  %-6.1f  %-9s  %-9s  %-9s  %-9s  %-12.4f  %-10.1f\n",
                  eps_sw, rho_sw, 100 * acc_frac_sw,
                  ot_acc_str, nn_acc_str, ot_f1_str, nn_f1_str,
                  d_sw$ot$transport_mass,
                  100 * d_sw$projection$outside_corrected_fraction))
    }
  }
  cat("  'acc_%'   = cells receiving an accepted transferred label\n")
  cat("  'OT_acc%' = overall accuracy of somalign on accepted cells\n")
  cat("  'NN_acc%' = overall accuracy of direct NN on the SAME accepted cells (like-for-like)\n")
  cat("  'OT_F1%'  = macro-F1 of somalign on accepted cells (immune to class imbalance)\n")
  cat("  'NN_F1%'  = macro-F1 of direct NN on the same accepted cells\n")
  cat("  High variation in OT_F1% or acc_% = unstable label transfer.\n")

  ## ── 1h. Per-label match_fraction (how well each reference pop is covered) ──
  sep("Nuñez 2023 — per-label OT coverage")
  ref_node_label <- colnames(ref_fl$label_prob)[apply(ref_fl$label_prob, 1, which.max)]
  for (cl in sort(unique(ref_node_label))) {
    idx <- which(ref_node_label == cl)
    mf_cl <- mf[idx]
    cat(sprintf("  %-10s  nodes=%2d  mean_mf=%.3f  min_mf=%.3f\n",
                cl, length(idx), mean(mf_cl), min(mf_cl)))
  }
}

## ─────────────────────────────────────────────────────────────────────────────
## 2.  CytoNorm CyTOF
## ─────────────────────────────────────────────────────────────────────────────

sep("CytoNorm CyTOF — loading and gating")

cyto_ok <- requireNamespace("CytoNorm", quietly = TRUE)
if (!cyto_ok) {
  cat("SKIP: CytoNorm not installed.\n")
} else {

  cyto_dir    <- system.file("extdata", package = "CytoNorm")
  all_files   <- list.files(cyto_dir, pattern = "\\.fcs$", full.names = TRUE)
  ref_files   <- sort(grep("_1\\.fcs$", all_files, value = TRUE))
  query_files <- sort(grep("_2\\.fcs$", all_files, value = TRUE))

  ## Channel selection — same as benchmark_real_data.R
  ref_ff0  <- read.FCS(ref_files[1], truncate_max_range = FALSE)
  pd       <- pData(parameters(ref_ff0))
  is_di    <- grepl("Di$", pd$name)
  has_desc <- !is.na(pd$desc)
  is_bc    <- grepl("^BC[0-9]", pd$desc) & !is.na(pd$desc)
  is_dna   <- pd$name %in% c("Ir191Di", "Ir193Di")
  mk_cyto  <- pd$name[is_di & has_desc & !is_bc & !is_dna]
  desc_map <- setNames(pd$desc, pd$name)   # channel → antibody name

  cofactor_cyto <- 5
  read_mat <- function(files) {
    mats <- lapply(files, function(f) {
      ff  <- read.FCS(f, truncate_max_range = FALSE)
      asinh(exprs(ff)[, mk_cyto] / cofactor_cyto)
    })
    m <- do.call(rbind, mats)
    colnames(m) <- desc_map[mk_cyto]
    m
  }

  c1 <- read_mat(ref_files)
  c2 <- read_mat(query_files)
  cat(sprintf("Loaded: batch1=%d  batch2=%d  markers=%d  (cofactor=%d)\n",
              nrow(c1), nrow(c2), ncol(c1), cofactor_cyto))

  ## ── 2a2. Per-marker batch shift ──
  sep("CytoNorm CyTOF — per-marker batch shift")
  n_feat_cyto <- ncol(c1)
  mu1_c <- colMeans(c1);  mu2_c <- colMeans(c2)
  sd1_c <- apply(c1, 2, sd);  sd2_c <- apply(c2, 2, sd)
  pool_sd_c <- sqrt((sd1_c^2 + sd2_c^2) / 2)
  z_c <- (mu2_c - mu1_c) / pool_sd_c
  cat(sprintf("  n_markers=%d  |z| mean=%.3f  median=%.3f  max=%.3f\n",
              n_feat_cyto, mean(abs(z_c)), median(abs(z_c)), max(abs(z_c))))
  cat(sprintf("  |z|>0.5: %d  |z|>1.0: %d  |z|>2.0: %d\n",
              sum(abs(z_c) > 0.5), sum(abs(z_c) > 1.0), sum(abs(z_c) > 2.0)))
  ord_c <- order(abs(z_c), decreasing = TRUE)
  cat("  Top 5 shifted markers:\n")
  for (i in seq_len(min(5L, n_feat_cyto)))
    cat(sprintf("    %-22s  z = %+.3f\n", names(z_c)[ord_c[i]], z_c[ord_c[i]]))
  cat("  NOTE: batches are DIFFERENT patients, not same sample re-measured;\n")
  cat("  'batch shift' here includes genuine inter-patient biology.\n")

  ## ── 2b. Gate batch1 and batch2 using lineage markers ──
  gate_cytof <- function(mat) {
    # arcsinh(x/5) > 0.5 ≈ raw > 2.5; sensible for CyTOF lineage markers
    thr <- 0.5
    cd3  <- mat[, "CD3"]
    cd4  <- mat[, "CD4"]
    cd8  <- mat[, "CD8a"]
    cd19 <- mat[, "CD19"]
    cd14 <- mat[, "CD14"]
    cd56 <- mat[, "CD56"]
    cd16 <- mat[, "CD16"]
    hladr<- mat[, "HLADR"]

    label <- rep("Other", nrow(mat))
    label[cd3 >= thr & cd4 >= thr & cd8 < thr]  <- "CD4.T"
    label[cd3 >= thr & cd8 >= thr & cd4 < thr]  <- "CD8.T"
    label[cd3 >= thr & cd4 >= thr & cd8 >= thr] <- "DP.T"   # double-positive (rare)
    label[cd3  < thr & cd19 >= thr]              <- "B"
    label[cd3  < thr & cd14 >= thr]              <- "Mono"
    label[cd3  < thr & cd56 >= thr & cd16 >= thr] <- "NK"
    label[cd3  < thr & hladr >= thr & cd14 < thr & cd19 < thr] <- "DC"
    label
  }

  labels_c1 <- gate_cytof(c1)
  labels_c2 <- gate_cytof(c2)

  cat("\nBatch1 gate counts:\n")
  print(sort(table(labels_c1), decreasing = TRUE))
  cat("\nBatch2 gate counts (ground truth):\n")
  print(sort(table(labels_c2), decreasing = TRUE))

  ## ── 2c. Train labeled reference ──
  sep("CytoNorm CyTOF — training labeled reference")
  set.seed(42)
  ref_cyto <- somalign_train_reference(
    c1,
    labels = labels_c1,
    grid   = somgrid(6, 6, "hexagonal"),
    rlen   = 100
  )
  cat("Reference: 36 nodes,", length(unique(labels_c1)), "label classes\n")

  ## ── 2d. Standard fit ──
  set.seed(123)
  query_cyto <- somalign_query(c2, ref_cyto,
                                grid = somgrid(6, 6, "hexagonal"), rlen = 100)
  fit_plain  <- somalign_fit(query_cyto, ref_cyto,
                              epsilon = 0.1, rho_query = 2, rho_ref = 2)
  res_plain  <- somalign_results(fit_plain)

  ## Anchored fit — all 6000 cells (3 patients × 2 batches) as anchor pairs
  fit_anc   <- somalign_fit_anchored(
    query_cyto, ref_cyto,
    anchor_old = c1, anchor_new = c2,
    rho_anchor = 1.0,
    epsilon = 0.1, rho_query = 2, rho_ref = 2
  )
  res_anc <- somalign_results(fit_anc)

  ## ── 2e. Full diagnostics — standard fit ──
  sep("CytoNorm CyTOF — full somalign_diagnostics() [standard]")
  diag_plain <- somalign_diagnostics(fit_plain)

  cat("\n### solver ###\n")
  cat(sprintf("  converged:               %s\n",   diag_plain$solver$converged))
  cat(sprintf("  iterations:              %d\n",   diag_plain$solver$iterations))
  cat(sprintf("  final_delta:             %.3e\n", diag_plain$solver$final_delta))
  cat(sprintf("  rel_marginal_row_error:  %.4f\n", diag_plain$solver$rel_marginal_row_error))
  cat(sprintf("  rel_marginal_col_error:  %.4f\n", diag_plain$solver$rel_marginal_col_error))
  cat(sprintf("  cost_scale:              %.4f\n", diag_plain$solver$cost_scale))

  cat("\n### ot ###\n")
  cat(sprintf("  transport_mass:      %.4f\n", diag_plain$ot$transport_mass))
  cat(sprintf("  max_row_mass_error:  %.4f\n", diag_plain$ot$max_row_mass_error))
  mf_p <- diag_plain$ot$match_fraction
  cat(sprintf("  match_fraction  [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(mf_p), quantile(mf_p,.25), median(mf_p), quantile(mf_p,.75), max(mf_p)))
  cat(sprintf("  nodes with match_fraction < 0.5:  %d / %d\n", sum(mf_p < 0.5), length(mf_p)))

  cat("\n### nodes ###\n")
  cn_p <- diag_plain$nodes$correction_norm
  cat(sprintf("  correction_norm [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(cn_p,na.rm=T), quantile(cn_p,.25,na.rm=T), median(cn_p,na.rm=T),
              quantile(cn_p,.75,na.rm=T), max(cn_p,na.rm=T)))
  cat(sprintf("  per-feature norm (/ sqrt(%d) markers):  med=%.3f  max=%.3f\n",
              n_feat_cyto, median(cn_p, na.rm=TRUE)/sqrt(n_feat_cyto),
              max(cn_p, na.rm=TRUE)/sqrt(n_feat_cyto)))
  cat("  (values ~0.5–1 = plausible per-marker batch offset; >2 = gross shift)\n")
  if (all(mf_p >= 1.0 - 1e-9)) {
    cat("\n  DIAGNOSTIC FLAG: match_fraction = 1.0 for ALL nodes — gate INOPERATIVE.\n")
  }

  cat("\n### projection ###\n")
  cat(sprintf("  outside_direct_fraction:    %s\n",
              fmt_pct(diag_plain$projection$outside_direct_fraction)))
  cat(sprintf("  outside_corrected_fraction: %s\n",
              fmt_pct(diag_plain$projection$outside_corrected_fraction)))

  ## Reference node purity diagnostic
  cat("\n### reference node purity (max label_prob per node) ###\n")
  node_purity <- apply(ref_cyto$label_prob, 1, max)
  cat(sprintf("  [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(node_purity), quantile(node_purity,.25), median(node_purity),
              quantile(node_purity,.75), max(node_purity)))
  cat(sprintf("  nodes with purity >= 0.6: %d / %d\n", sum(node_purity >= 0.6), length(node_purity)))
  cat(sprintf("  nodes with purity >= 0.4: %d / %d\n", sum(node_purity >= 0.4), length(node_purity)))
  cat("  (NOTE: transferred_label_accepted requires transferred confidence >= 0.6;\n")
  cat("   lowering epsilon sharpens the OT plan and raises per-node confidence.)\n")

  ## Per-node label confidence (transferred, before cell mapping)
  tl_node_conf <- fit_plain$label_transfer$confidence
  cat(sprintf("\n### per-node transferred label confidence ###\n"))
  cat(sprintf("  [min/q25/med/q75/max]: %.3f / %.3f / %.3f / %.3f / %.3f\n",
              min(tl_node_conf, na.rm=T), quantile(tl_node_conf,.25,na.rm=T),
              median(tl_node_conf,na.rm=T), quantile(tl_node_conf,.75,na.rm=T),
              max(tl_node_conf, na.rm=T)))
  cat(sprintf("  nodes with confidence >= 0.6: %d / %d\n",
              sum(tl_node_conf >= 0.6, na.rm=T), sum(is.finite(tl_node_conf))))
  cat(sprintf("  nodes with confidence >= 0.4: %d / %d\n",
              sum(tl_node_conf >= 0.4, na.rm=T), sum(is.finite(tl_node_conf))))

  ## ── 2f. Label transfer — standard ──
  sep("CytoNorm CyTOF — label transfer [standard fit]")

  cat(sprintf("Transferred label accepted: %s  rejected: %s\n",
              fmt_pct(mean(res_plain$transferred_label_accepted, na.rm=TRUE)),
              fmt_pct(mean(!res_plain$transferred_label_accepted, na.rm=TRUE))))
  cat("\nfinal_status:\n")
  print(sort(table(res_plain$final_status), decreasing = TRUE))

  acc_p <- which(res_plain$transferred_label_accepted)
  tl_p  <- res_plain$transferred_label[acc_p]
  gt_p  <- labels_c2[acc_p]
  old_p <- res_plain$old_som_label[acc_p]

  cat(sprintf("\nConfidence — mean: %.3f  median: %.3f  min: %.3f\n",
              mean(res_plain$transferred_label_confidence[acc_p], na.rm=TRUE),
              median(res_plain$transferred_label_confidence[acc_p], na.rm=TRUE),
              min(res_plain$transferred_label_confidence[acc_p], na.rm=TRUE)))

  cat("\n-- Standard: transferred_label vs batch2 GT (accepted only) --\n")
  conf_table(gt_p, tl_p, "")

  cat("\n-- Standard: old_som_label (direct NN) vs batch2 GT --\n")
  conf_table(gt_p, old_p, "")

  cat(sprintf("Agreement old_som_label ↔ transferred_label: %.1f%%\n",
              100 * mean(tl_p == old_p)))

  ## ── 2f1. Direct NN on ALL cells ──
  cat(sprintf("\n-- old_som_label (direct NN) vs batch2 GT — ALL %d cells --\n", nrow(res_plain)))
  conf_table(labels_c2, res_plain$old_som_label, "")

  ## ── 2f_acc. Label transfer acceptance/rejection by true class ──
  cat("\n-- label transfer acceptance by true class (standard fit) --\n")
  cat("  (rejected = confidence guardrail triggered, NOT mislabeled)\n")
  cat(sprintf("  %-10s  %7s  %7s  %7s  %s\n",
              "class", "n_total", "n_acc", "acc_%", "precision_if_acc"))
  for (cl in sort(unique(labels_c2))) {
    mask_cl  <- labels_c2 == cl
    n_cl     <- sum(mask_cl)
    mask_acc <- mask_cl & res_plain$transferred_label_accepted
    n_acc    <- sum(mask_acc, na.rm = TRUE)
    if (n_acc > 0) {
      prec <- mean(res_plain$transferred_label[mask_acc] == cl, na.rm = TRUE)
      cat(sprintf("  %-10s  %7d  %7d  %6.1f%%  %.1f%%\n",
                  cl, n_cl, n_acc, 100 * n_acc / n_cl, 100 * prec))
    } else {
      cat(sprintf("  %-10s  %7d  %7d  %6.1f%%  [none accepted]\n",
                  cl, n_cl, n_acc, 100 * n_acc / n_cl))
    }
  }

  ## ── 2g. Label transfer — anchored ──
  sep("CytoNorm CyTOF — label transfer [anchored fit]")
  diag_anc <- somalign_diagnostics(fit_anc)
  cat(sprintf("Anchor coverage: %d / %d nodes (%.0f%%)\n",
              fit_anc$anchors$nodes_covered,
              nrow(ref_cyto$codebook),
              100 * fit_anc$anchors$coverage_fraction))
  cat(sprintf("Anchored transport_mass: %.4f  (plain: %.4f)\n",
              diag_anc$ot$transport_mass, diag_plain$ot$transport_mass))

  acc_a <- which(res_anc$transferred_label_accepted)
  tl_a  <- res_anc$transferred_label[acc_a]
  gt_a  <- labels_c2[acc_a]

  cat(sprintf("\nTransferred label accepted: %s  rejected: %s\n",
              fmt_pct(mean(res_anc$transferred_label_accepted, na.rm=TRUE)),
              fmt_pct(mean(!res_anc$transferred_label_accepted, na.rm=TRUE))))

  cat("\n-- Anchored: transferred_label vs batch2 GT (accepted only) --\n")
  conf_table(gt_a, tl_a, "")

  ## ── 2g2. Rerun with relaxed confidence threshold ──
  sep("CytoNorm CyTOF — label transfer [confidence_threshold=0.4]")
  cat("NOTE: relaxing the acceptance threshold from 0.6 to 0.4 shows what IS transferred\n")
  cat("  for nodes whose per-node confidence falls between the two thresholds.\n\n")
  fit_low <- somalign_fit(query_cyto, ref_cyto,
                           epsilon = 0.1, rho_query = 2, rho_ref = 2,
                           confidence_threshold = 0.4)
  res_low <- somalign_results(fit_low)

  cat(sprintf("Transferred label accepted: %s  rejected: %s\n",
              fmt_pct(mean(res_low$transferred_label_accepted, na.rm=TRUE)),
              fmt_pct(mean(!res_low$transferred_label_accepted, na.rm=TRUE))))

  acc_l <- which(res_low$transferred_label_accepted)
  tl_l  <- res_low$transferred_label[acc_l]
  gt_l  <- labels_c2[acc_l]
  old_l <- res_low$old_som_label[acc_l]

  if (length(acc_l) > 0) {
    cat(sprintf("Confidence (accepted) — mean: %.3f  median: %.3f\n",
                mean(res_low$transferred_label_confidence[acc_l], na.rm=TRUE),
                median(res_low$transferred_label_confidence[acc_l], na.rm=TRUE)))
    cat("\n-- confidence_threshold=0.4: transferred_label vs batch2 GT --\n")
    conf_table(gt_l, tl_l, "")
    cat("\n-- confidence_threshold=0.4: old_som_label vs batch2 GT --\n")
    conf_table(gt_l, old_l, "")
    cat(sprintf("Agreement old_som_label ↔ transferred_label (thr=0.4): %.1f%%\n",
                100 * mean(tl_l == old_l)))
  }

  ## ── 2j. Sensitivity sweep: epsilon × rho (CytoNorm) ──
  sep("CytoNorm CyTOF — sensitivity sweep (epsilon x rho)")
  cat("Scores on ACCEPTED cells only (same set for OT and NN — like-for-like).\n")
  cat(sprintf("  %-7s %-5s  %-6s  %-9s  %-9s  %-9s  %-9s  %-12s  %-10s\n",
              "epsilon", "rho", "acc_%", "OT_acc%", "NN_acc%",
              "OT_F1%", "NN_F1%", "trans_mass", "out_corr%"))
  for (eps_sw in c(0.05, 0.1, 0.5)) {
    for (rho_sw in c(2, 10)) {
      f_sw_c <- somalign_fit(query_cyto, ref_cyto,
                              epsilon = eps_sw, rho_query = rho_sw, rho_ref = rho_sw)
      r_sw_c <- somalign_results(f_sw_c)
      d_sw_c <- somalign_diagnostics(f_sw_c)
      acc_frac_c  <- mean(r_sw_c$transferred_label_accepted, na.rm = TRUE)
      acc_cells_c <- which(r_sw_c$transferred_label_accepted)
      if (length(acc_cells_c) > 0) {
        gt_c   <- labels_c2[acc_cells_c]
        ot_c   <- r_sw_c$transferred_label[acc_cells_c]
        nn_c   <- r_sw_c$old_som_label[acc_cells_c]
        ot_acc_str_c <- sprintf("%.1f", 100 * mean(ot_c == gt_c, na.rm = TRUE))
        nn_acc_str_c <- sprintf("%.1f", 100 * mean(nn_c == gt_c, na.rm = TRUE))
        ot_f1_str_c  <- sprintf("%.1f", 100 * macro_f1(gt_c, ot_c))
        nn_f1_str_c  <- sprintf("%.1f", 100 * macro_f1(gt_c, nn_c))
      } else {
        ot_acc_str_c <- nn_acc_str_c <- ot_f1_str_c <- nn_f1_str_c <- "  —  "
      }
      cat(sprintf("  %-7.2f %-5g  %-6.1f  %-9s  %-9s  %-9s  %-9s  %-12.4f  %-10.1f\n",
                  eps_sw, rho_sw, 100 * acc_frac_c,
                  ot_acc_str_c, nn_acc_str_c, ot_f1_str_c, nn_f1_str_c,
                  d_sw_c$ot$transport_mass,
                  100 * d_sw_c$projection$outside_corrected_fraction))
    }
  }
  cat("  NOTE: CytoNorm batches are different patients; low accuracy reflects\n")
  cat("  genuine inter-patient biology, not method failure.\n")
  cat("  'OT_F1%'/'NN_F1%' are macro-F1 (equal weight per class, reveals minority-class collapse).\n")
  cat("  High variation in OT_F1% = unstable label transfer.\n")

  ## ── 2h. Per-label OT coverage ──
  sep("CytoNorm CyTOF — per-label OT coverage")
  ref_node_label_c <- colnames(ref_cyto$label_prob)[
    apply(ref_cyto$label_prob, 1, which.max)]
  cat("Standard fit:\n")
  for (cl in sort(unique(ref_node_label_c))) {
    idx <- which(ref_node_label_c == cl)
    cat(sprintf("  %-10s  nodes=%2d  mean_mf=%.3f  min_mf=%.3f\n",
                cl, length(idx), mean(mf_p[idx]), min(mf_p[idx])))
  }

  ## ── 2i. OT transport plan — label-level mass flow ──
  sep("CytoNorm CyTOF — OT mass flow (query node → reference label)")
  # Query nodes have no labels; derive each node's dominant reference label
  # from the transport plan × reference label_prob.
  tp <- fit_plain$transport_plan   # n_query_nodes × n_ref_nodes
  query_label_mass   <- tp %*% ref_cyto$label_prob   # n_query × n_labels
  query_node_label_c <- colnames(ref_cyto$label_prob)[
    apply(query_label_mass, 1, which.max)]

  ref_labels_cyto <- sort(unique(ref_node_label_c))
  qll <- factor(query_node_label_c, ref_labels_cyto)
  rll <- factor(ref_node_label_c,   ref_labels_cyto)

  flow <- matrix(0, length(ref_labels_cyto), length(ref_labels_cyto),
                 dimnames = list(Query_dominant = ref_labels_cyto,
                                 Ref_dominant   = ref_labels_cyto))
  for (i in seq_len(nrow(tp))) {
    ql <- as.character(qll[i])
    for (j in seq_len(ncol(tp))) {
      rl <- as.character(rll[j])
      flow[ql, rl] <- flow[ql, rl] + tp[i, j]
    }
  }
  flow_norm <- sweep(flow, 1, rowSums(flow) + 1e-10, "/")
  cat("Fraction of query-node mass transported to each reference label:\n")
  print(round(flow_norm, 3))
}

cat("\n── Done ────────────────────────────────────────────────────────\n")
