#!/usr/bin/env Rscript
# =============================================================================
# BMV pilot: label-transfer validation and supervised tuning report
#
# Demonstrates the item-2/item-3 validation harness on real data. The BMV
# query batch has no ground-truth labels (that is what somalign transfers), so
# the honest way to measure label-transfer accuracy is held-out cross-
# validation WITHIN the labelled reference (pilot) batch: hold out gated pilot
# cells, train a reference SOM on the rest, transfer labels to the held-out
# cells, and compare to their biological gate. This measures how faithfully the
# transport plan recovers a known label across the SOM discretisation.
#
# Labels: `gate_lineage2` (6 immune lineages) -- the same populations the v3
# topology analysis found being merged by the (separate) correction path.
#
# Usage:  Rscript bmv_label_validation.R [n_subsample]
# Output: analysis/label-validation/bmv_label_validation_results.rds + console.
# =============================================================================

.libPaths("/exports/para-lipg-hpc/mdmanurung/R/4.5")
suppressPackageStartupMessages({
  library(somalign)
  library(kohonen)
})

args <- commandArgs(trailingOnly = TRUE)
n_sub <- if (length(args) >= 1) as.integer(args[1]) else 120000L
seed  <- 1L
set.seed(seed)

PILOT <- "/exports/para-lipg-hpc/Xuran/data/pilot/aurora_exvivo/processed/aurora_exvivo_xyfclusters.qs"
OUT   <- "/exports/para-lipg-hpc/mdmanurung/somalign/analysis/label-validation/bmv_label_validation_results.rds"

FEATURES <- c("TCRgd_CD141", "CD8", "CD3", "CD19", "CD159c_NKG2C", "CD38", "CD57",
              "CD161", "CD1c", "HLA-DR", "CD11c", "CD25", "CD16", "CD185_CXCR5",
              "CD14", "CD303_TCRvd2", "CD4", "CD294_CRTH2", "CD21", "CD56",
              "FoxP3", "Tbet", "CD197_CCR7", "CD127", "CD45RA", "CD27", "CD7")
LABEL_COL <- "gate_lineage2"

cat(sprintf("[%s] loading pilot data ...\n", format(Sys.time(), "%H:%M:%S")))
df <- qs2::qs_read(PILOT, nthreads = 8)
cat(sprintf("  %d cells x %d cols\n", nrow(df), ncol(df)))

labels_all <- as.character(df[[LABEL_COL]])
keep <- !is.na(labels_all) & labels_all != "" & stats::complete.cases(df[, FEATURES])

# Stratified subsample for tractable, repeated cross-validation.
idx_by_class <- split(which(keep), labels_all[keep])
per_class <- max(1L, floor(n_sub / length(idx_by_class)))
sub_idx <- unlist(lapply(idx_by_class, function(ix)
  if (length(ix) <= per_class) ix else sample(ix, per_class)), use.names = FALSE)
X <- as.matrix(df[sub_idx, FEATURES])
y <- labels_all[sub_idx]
rm(df); gc()
cat(sprintf("  subsampled %d cells across %d classes:\n", nrow(X), length(unique(y))))
print(sort(table(y), decreasing = TRUE))

grid <- kohonen::somgrid(15, 15, "hexagonal")   # 225 nodes, as in the real run

# --- 1. Cross-validated accuracy + calibration at the shipped epsilon = 0.1 ---
cat(sprintf("\n[%s] 5-fold cross-validation at epsilon = 0.1 ...\n",
            format(Sys.time(), "%H:%M:%S")))
cv <- somalign_cross_validate(X, y, grid = grid, k = 5, epsilon = 0.1,
                              rlen = 40, seed = seed)
print(cv)
cat("\nPer-class (pooled):\n"); print(cv$metrics$per_class)
cat("\nReliability:\n"); print(cv$calibration$table)

# --- 2. Supervised epsilon tuning against MCC ---
cat(sprintf("\n[%s] tuning epsilon against MCC ...\n", format(Sys.time(), "%H:%M:%S")))
tuned <- somalign_tune(
  X, y, grid = grid,
  param_grid = data.frame(epsilon = c(0.02, 0.05, 0.1, 0.2, 0.5)),
  k = 5, rlen = 40, metric = "mcc", seed = seed
)
print(tuned)
cat("\nFull tuning grid:\n"); print(tuned$grid)

shipped <- tuned$grid[tuned$grid$epsilon == 0.1, ]
best    <- tuned$best
cat(sprintf("\nShipped epsilon=0.1 MCC = %.4f  vs  tuned epsilon=%.3g MCC = %.4f  (delta %+.4f)\n",
            shipped$mcc, best$epsilon, best$mcc, best$mcc - shipped$mcc))

# Interpretation: MCC/accuracy are scored on accepted cells only, so large
# epsilon can "win" by abstaining on hard cells (watch the coverage column).
# macro_f1 is coverage-robust because abstaining on rare classes tanks their
# recall. The best operating point balances MCC, coverage, and calibration.
best_f1 <- tuned$grid[which.max(tuned$grid$macro_f1), ]
cat(sprintf("Best macro_f1 = %.4f at epsilon=%.3g (coverage %.1f%%, ECE %.4f) -- a more coverage-robust choice.\n",
            best_f1$macro_f1, best_f1$epsilon, 100 * best_f1$coverage, best_f1$ece))

saveRDS(list(cv = cv, tuned = tuned, n_subsample = nrow(X),
             classes = sort(table(y), decreasing = TRUE)), OUT)
cat(sprintf("\n[%s] saved -> %s\n", format(Sys.time(), "%H:%M:%S"), OUT))
