#!/usr/bin/env Rscript
# End-to-end notebook smoke test with 1M-cell subset.
# Run from repo root:
#   /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
#     scripts/test_notebook_subset.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(qs2)
  library(kohonen)
  library(ggplot2)
  library(patchwork)
  pkgload::load_all(".", quiet = TRUE)   # dev somalign
})

# =============================================================================
# Paths
# =============================================================================
PATH_PILOT_DATA   <- "/exports/para-lipg-hpc/Xuran/data/pilot/aurora_exvivo/processed/aurora_exvivo_xyfclusters.qs"
PATH_PILOT_SOM    <- "/exports/para-lipg-hpc/Xuran/results/pilot/aurora_exvivo/SOM/model/xyf_som_model.rds"
PATH_BMV_DATA     <- "/exports/para-lipg-hpc/Xuran/data/bmv/aurora_exvivo/processed/aurora_exvivo_cycombine_normalized.qs"
PATH_BMV_SOM      <- "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/cluster/model/xyf_som_model.rds"
PATH_PARAMS_PILOT <- "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/mapping/batch_correction/winsor_scaler_params_pilot.rds"
PATH_PARAMS_QUERY <- "/exports/para-lipg-hpc/Xuran/results/bmv/aurora_exvivo/mapping/batch_correction/winsor_scaler_params_new.rds"

# =============================================================================
# Marker definitions (extracted from pilot, then pilot data is freed)
# =============================================================================
cat("=== Extracting markers from pilot data ===\n")
dat_pilot <- qs_read(PATH_PILOT_DATA, nthreads = 10)
cat(sprintf("pilot: %d rows\n", nrow(dat_pilot)))

get_marker_columns <- function(df) {
  non_marker_cols <- c(
    "rownum","Original","batch_id","Time","begin","time_sec","time","Zombie NIR-A",
    "sec_from_start","cell_id","AF1-A","AF-A","AF Index","peacoqc-passed",
    "OmiqFilter","filename","LIVE DEAD NIR-A","time_cut","fcs_filename","condition",
    "SSC-H","SSC-A","FSC-H","FSC-A","SSC-B-H","SSC-B-A","Orig_Row_Number",
    "Viability","Autofluorescence","Autofluorescence Index","Original_ID",
    "GateName","chmi_qPCR","control","study_id","timepoint_chron","batch",
    "GateNameNum","som_node","id","original_id","gate","gate_lineage2",
    "som_node","som_meta30","xyf_node","xyf_cluster","xyf_cluster_merged",
    "annotation"
  )
  setdiff(colnames(df), non_marker_cols)
}
CELL_STATE_MARKERS <- c(
  "GLUT1","CD71","LDHA","G6PD","ATP5A","CS","CPT1A","ATGL","GLS",
  "CD279_PD-1","CD279_PD1","PD-1","CD152_CTLA4","CTLA-4",
  "IL-1b","IL-2","IL-6","IL-10","IFNg","TNF","Th2"
)
markers          <- get_marker_columns(dat_pilot)
cellTypeMarkers  <- setdiff(markers, CELL_STATE_MARKERS)
cat(sprintf("Cell-type markers: %d\n", length(cellTypeMarkers)))
rm(dat_pilot); gc()

# =============================================================================
# Winsorization params + query codebook in reference-scaled space
# =============================================================================
params_pilot <- readRDS(PATH_PARAMS_PILOT)
params_query <- readRDS(PATH_PARAMS_QUERY)

som_bmv <- readRDS(PATH_BMV_SOM)
cat(sprintf("som_bmv unit.classif: %d entries\n", length(som_bmv$unit.classif)))

codes_query <- som_bmv$codes[[1]] |>
  as.data.frame() |>
  rename(TCRgd_CD141 = CD141, CD303_TCRvd2 = CD303, FoxP3 = Foxp3)

features  <- names(params_pilot[[1]])
mad_q     <- params_query$mad[features]
med_q     <- params_query$median[features]
low_ref   <- params_pilot$lower[features]
upp_ref   <- params_pilot$upper[features]
med_ref   <- params_pilot$median[features]
mad_ref   <- params_pilot$mad[features]

codes_query      <- as.matrix(codes_query[, features])
codes_unscaled   <- sweep(sweep(codes_query, 2, mad_q, "*"), 2, med_q, "+")
codes_wins_ref   <- pmin(pmax(codes_unscaled,
                              rep(low_ref, each = nrow(codes_unscaled))),
                         rep(upp_ref, each = nrow(codes_unscaled)))
codes_ref_scaled <- sweep(sweep(codes_wins_ref, 2, med_ref, "-"), 2, mad_ref, "/")
colnames(codes_ref_scaled) <- features
cat(sprintf("codes_ref_scaled: %d nodes × %d features\n",
            nrow(codes_ref_scaled), ncol(codes_ref_scaled)))
rm(codes_query, codes_unscaled, codes_wins_ref); gc()

# =============================================================================
# Reference from pilot SOM (all 44.6M cells) — frees som_pilot afterwards
# =============================================================================
cat("\n=== Building reference ===\n")
som_pilot <- readRDS(PATH_PILOT_SOM)
cat(sprintf("som_pilot unit.classif: %d entries\n", length(som_pilot$unit.classif)))
reference_test <- somalign_reference_from_som(
  som_pilot,
  center = params_pilot$median[cellTypeMarkers],
  scale  = pmax(params_pilot$mad[cellTypeMarkers], 1e-8),
  codebook_space = "reference_scaled"
)
cat(sprintf("reference: %d nodes, label_prob cols: %s\n",
            nrow(reference_test$codebook),
            paste(colnames(reference_test$label_prob), collapse = ", ")))
rm(som_pilot); gc()

# =============================================================================
# Pre-flight: codebook alignment check (using full SOM masses)
# =============================================================================
cat("\n=== Pre-flight codebook check ===\n")
query_masses_est <- tabulate(som_bmv$unit.classif, nrow(codes_ref_scaled)) /
  length(som_bmv$unit.classif)
cb_check <- somalign_check_codebook_alignment(
  query_codebook = codes_ref_scaled,
  reference      = reference_test,
  query_masses   = query_masses_est,
  epsilon        = 0.1
)
print(cb_check)

# =============================================================================
# Load bmv data; filter to labeled cells; free the full dataset
# =============================================================================
cat("\n=== Loading and filtering bmv data ===\n")
dat_bmv <- qs_read(PATH_BMV_DATA, nthreads = 10)
cat(sprintf("bmv: %d rows\n", nrow(dat_bmv)))
dat_bmv <- dat_bmv |>
  rename(TCRgd_CD141 = CD141,
         CD303_TCRvd2 = CD303,
         FoxP3 = Foxp3,
         `CD279_PD-1` = CD279_PD1)
dat_bmv_labeled <- dat_bmv |> filter(!OmiqFilter == 0)
cat(sprintf("bmv labeled: %d rows | som_bmv unit.classif: %d entries\n",
            nrow(dat_bmv_labeled), length(som_bmv$unit.classif)))
rm(dat_bmv); gc()

# =============================================================================
# RECONSTRUCT training subset — the bmv SOM was trained on a 40% random sample
# (set.seed(42), group_by(fcs_filename, OmiqFilter), slice_sample(prop=0.4))
# so unit.classif[i] corresponds to row i of that reconstruction.
# =============================================================================
cat("\n=== Reconstructing bmv SOM training subset (40% of labeled cells) ===\n")
set.seed(42)
dat_bmv_train <- dat_bmv_labeled |>
  group_by(fcs_filename, OmiqFilter) |>
  slice_sample(prop = 0.4) |>
  ungroup()
cat(sprintf("Reconstructed training subset: %d rows (SOM unit.classif: %d entries)\n",
            nrow(dat_bmv_train), length(som_bmv$unit.classif)))
stopifnot(nrow(dat_bmv_train) == length(som_bmv$unit.classif))
rm(dat_bmv_labeled); gc()

# =============================================================================
# SUBSET: 1M cells from the training subset + matching unit.classif
# =============================================================================
cat("\n=== Subsetting to 1M cells ===\n")
N_SUBSET <- 1e6L
set.seed(123)
subset_idx   <- sort(sample.int(nrow(dat_bmv_train), N_SUBSET))
dat_bmv_sub  <- dat_bmv_train[subset_idx, ]
unit_sub     <- som_bmv$unit.classif[subset_idx]
rm(dat_bmv_train); gc()

som_bmv_sub <- som_bmv
som_bmv_sub$unit.classif <- unit_sub
cat(sprintf("Subset: %d rows, unit.classif: %d entries\n",
            nrow(dat_bmv_sub), length(som_bmv_sub$unit.classif)))

# =============================================================================
# Query from SOM (1M cell subset)
# =============================================================================
cat("\n=== Building query (somalign_query_from_som) ===\n")
query <- somalign_query_from_som(
  som_bmv_sub,
  as.matrix(dat_bmv_sub[, cellTypeMarkers]),
  reference_test,
  codebook       = codes_ref_scaled,
  codebook_space = "reference_scaled"
)
rm(dat_bmv_sub, som_bmv_sub, som_bmv, codes_ref_scaled); gc()
cat(sprintf("Query: %d cells × %d features, %d nodes\n",
            nrow(query$data), ncol(query$data), nrow(query$codebook)))

# =============================================================================
# Fit
# =============================================================================
cat("\n=== somalign_fit ===\n")
fit <- somalign_fit(query, reference_test)

# =============================================================================
# Diagnostics
# =============================================================================
cat("\n=== Diagnostics ===\n")
diag <- somalign_diagnostics(fit)
cat(sprintf("Solver:                   %s (converged: %s)\n",
            diag$solver$used, diag$solver$converged))
cat(sprintf("Transport mass:           %.4f\n",    diag$ot$transport_mass))
cat(sprintf("Mean match fraction:      %.3f\n",
            mean(diag$ot$match_fraction[is.finite(diag$ot$match_fraction)])))
cat(sprintf("Max row mass error:       %.2e\n",    diag$ot$max_row_mass_error))
cat(sprintf("Max col mass error:       %.2e\n",    diag$ot$max_col_mass_error))
cat(sprintf("Outside direct %%:         %.1f%%\n",
            100 * diag$projection$outside_direct_fraction))
cat(sprintf("Outside corrected %%:      %.1f%%\n",
            100 * diag$projection$outside_corrected_fraction))
cat(sprintf("Correction-allowed nodes: %d / %d\n",
            sum(diag$nodes$correction_allowed), length(diag$nodes$correction_allowed)))
if (any(diag$nodes$correction_allowed)) {
  cat(sprintf("Mean correction norm:     %.4f\n",
              mean(diag$nodes$correction_norm[diag$nodes$correction_allowed])))
}

match_frac <- diag$ot$match_fraction
cat(sprintf("\nMatch fraction summary (ref nodes):\n"))
print(summary(match_frac[is.finite(match_frac)]))
cat(sprintf("Nodes with match_fraction < 0.05: %d\n",
            sum(is.finite(match_frac) & match_frac < 0.05)))

# =============================================================================
# Results
# =============================================================================
cat("\n=== Results ===\n")
sample_metadata <- data.frame(batch = rep("bmv", nrow(query$data)))
results <- somalign_results(fit, data = sample_metadata)
cat(sprintf("Label transfer acceptance rate: %.1f%%\n",
            100 * mean(results$transferred_label_accepted, na.rm = TRUE)))
cat("Label distribution (accepted):\n")
print(table(results$transferred_label[results$transferred_label_accepted]))

cat("\n=== head(results) ===\n")
print(head(results[, c(
  "sample_id",
  "query_som_unit",
  "old_som_unit",
  "old_som_label",
  "outside_reference_distance",
  "final_status",
  "corrected_som_unit",
  "correction_norm",
  "transferred_label",
  "transferred_label_confidence",
  "transferred_label_accepted"
)]))

cat("\n=== Diagnostics & Visualization ===\n")

OUT_DIR <- "scripts/diagnostics"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
PDF_PATH <- file.path(OUT_DIR, "projection_diagnostics.pdf")

# Helper for themed plots
theme_diag <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(color = "grey40", size = 9),
          strip.text = element_text(face = "bold"))
}

node_df <- as.data.frame(diag$nodes)
lt_df   <- fit$label_transfer   # per-node label transfer info
row_err <- abs(diag$ot$row_mass - diag$ot$query_mass)
col_err <- abs(diag$ot$col_mass - diag$ot$reference_mass)

# Annotate nodes with their reference label (majority label of top reference node)
top_ref_node <- apply(fit$correspondence, 1, which.max)
ref_labels   <- reference_test$label_prob
if (!is.null(ref_labels)) {
  node_df$top_ref_label <- ref_labels[top_ref_node, ] |>
    apply(1, function(x) if (max(x) >= 0.5) colnames(ref_labels)[which.max(x)] else "unlabeled")
} else {
  node_df$top_ref_label <- factor(top_ref_node)
}

pdf(PDF_PATH, width = 14, height = 10, onefile = TRUE)

# =============================================================================
# PAGE 1: OT solver quality
# =============================================================================

# 1a. Node mass balance scatter
p1a <- ggplot(node_df, aes(x = query_mass, y = transported_mass, color = match_fraction)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_point(alpha = 0.65, size = 1.8) +
  scale_color_viridis_c("Match\nfraction", option = "plasma", limits = c(0, 1)) +
  labs(title = "Node mass balance",
       subtitle = "Dashed = perfect balance; color = match fraction",
       x = "Query node mass", y = "Transported mass") +
  theme_diag()

# 1b. Match fraction histogram with worst-node callouts
n_low <- sum(node_df$match_fraction < 0.05)
p1b <- ggplot(node_df, aes(x = match_fraction)) +
  geom_histogram(bins = 50, fill = "#2166ac", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, color = "#d73027", linetype = "dashed") +
  annotate("text", x = 0.06, y = Inf, label = "threshold = 0.05",
           color = "#d73027", size = 3, hjust = 0, vjust = 1.5) +
  labs(title = "Match fraction per query node",
       subtitle = sprintf("%d / %d nodes below 0.05  |  median = %.3f",
                          n_low, nrow(node_df), median(node_df$match_fraction)),
       x = "Match fraction", y = "# nodes") +
  theme_diag()

# 1c. Marginal error distributions (row vs col)
err_df <- rbind(
  data.frame(error = row_err, marginal = "Row (query)"),
  data.frame(error = col_err, marginal = "Col (reference)")
)
p1c <- ggplot(err_df, aes(x = error, fill = marginal)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.70,
                 color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c("Row (query)" = "#2166ac", "Col (reference)" = "#d73027")) +
  scale_x_log10(labels = scales::label_scientific()) +
  labs(title = "Marginal mass error",
       subtitle = sprintf("Max row: %.1e  |  Max col: %.1e  |  Converged: %s",
                          max(row_err), max(col_err), diag$solver$converged),
       x = "Absolute error (log scale)", y = "# nodes", fill = NULL) +
  theme_diag()

# 1d. Worst query nodes by match fraction
worst_nodes <- node_df |>
  arrange(match_fraction) |>
  head(20) |>
  mutate(node_label = sprintf("node %d\n(%s)", query_node, top_ref_label))

p1d <- ggplot(worst_nodes, aes(x = reorder(node_label, match_fraction), y = match_fraction,
                                fill = match_fraction)) +
  geom_col(color = "white") +
  geom_hline(yintercept = 0.05, color = "#d73027", linetype = "dashed") +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1), guide = "none") +
  coord_flip() +
  labs(title = "20 worst query nodes (match fraction)",
       x = NULL, y = "Match fraction") +
  theme_diag()

print((p1a | p1b) / (p1c | p1d) +
        plot_annotation(title = "PAGE 1 — Optimal Transport solver quality",
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))

# =============================================================================
# PAGE 2: Correction quality
# =============================================================================

# 2a. Correction norm histogram
p2a <- ggplot(node_df, aes(x = correction_norm, fill = correction_allowed)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.75,
                 color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c("TRUE" = "#1a9641", "FALSE" = "#d73027"),
                    labels = c("TRUE" = "Allowed", "FALSE" = "Suppressed")) +
  labs(title = "Node correction norm distribution",
       subtitle = sprintf("Mean (allowed): %.2f  |  Max: %.2f  |  %d / %d nodes corrected",
                          mean(node_df$correction_norm[node_df$correction_allowed]),
                          max(node_df$correction_norm),
                          sum(node_df$correction_allowed), nrow(node_df)),
       x = "Correction norm", y = "# nodes", fill = "Correction") +
  theme_diag()

# 2b. Match fraction vs correction norm
p2b <- ggplot(node_df, aes(x = match_fraction, y = correction_norm,
                             color = correction_allowed)) +
  geom_point(alpha = 0.55, size = 1.8) +
  scale_color_manual(values = c("TRUE" = "#1a9641", "FALSE" = "#d73027"),
                     labels = c("TRUE" = "Allowed", "FALSE" = "Suppressed")) +
  labs(title = "Match fraction vs correction magnitude",
       subtitle = "Nodes with low match fraction and large corrections are worth inspecting",
       x = "Match fraction", y = "Correction norm", color = "Correction") +
  theme_diag()

# 2c. Direct vs corrected outside-reference fraction
proj_cmp <- data.frame(
  Projection = factor(c("Direct", "Corrected"), levels = c("Direct", "Corrected")),
  outside_pct = c(100 * diag$projection$outside_direct_fraction,
                  100 * diag$projection$outside_corrected_fraction)
)
p2c <- ggplot(proj_cmp, aes(x = Projection, y = outside_pct, fill = Projection)) +
  geom_col(width = 0.5, color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", outside_pct)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("Direct" = "#d73027", "Corrected" = "#1a9641"), guide = "none") +
  labs(title = "% cells outside reference thresholds",
       subtitle = "Lower is better after correction",
       y = "% outside threshold", x = NULL) +
  ylim(0, max(proj_cmp$outside_pct) * 1.15) +
  theme_diag()

# 2d. Worst nodes by correction norm (large corrections that were allowed)
worst_corr <- node_df |>
  filter(correction_allowed) |>
  arrange(desc(correction_norm)) |>
  head(20) |>
  mutate(node_label = sprintf("node %d\n(%s)", query_node, top_ref_label))

p2d <- ggplot(worst_corr, aes(x = reorder(node_label, correction_norm),
                               y = correction_norm, fill = match_fraction)) +
  geom_col(color = "white") +
  scale_fill_viridis_c("Match\nfraction", option = "plasma", limits = c(0, 1)) +
  coord_flip() +
  labs(title = "20 largest corrections (allowed nodes)",
       x = NULL, y = "Correction norm") +
  theme_diag()

print((p2a | p2b) / (p2c | p2d) +
        plot_annotation(title = "PAGE 2 — Correction quality",
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))

# =============================================================================
# PAGE 3: Label transfer quality
# =============================================================================

# 3a. Old → transferred label confusion heatmap (accepted cells only)
confusion <- results |>
  filter(transferred_label_accepted) |>
  count(old_som_label, transferred_label) |>
  group_by(old_som_label) |>
  mutate(pct = n / sum(n)) |>
  ungroup()

p3a <- ggplot(confusion, aes(x = transferred_label, y = old_som_label, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * pct)), size = 3.5, fontface = "bold") +
  scale_fill_viridis_c("Fraction", option = "plasma", limits = c(0, 1)) +
  labs(title = "Label transfer confusion (accepted cells)",
       subtitle = "Row-normalized: each row sums to 100%",
       x = "Transferred label", y = "Old SOM label") +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 3b. Acceptance rate by old label
accept_by_label <- results |>
  group_by(old_som_label) |>
  summarise(n = n(),
            rate = mean(transferred_label_accepted, na.rm = TRUE),
            .groups = "drop")

p3b <- ggplot(accept_by_label, aes(x = reorder(old_som_label, rate),
                                    y = 100 * rate, fill = rate)) +
  geom_col(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%\nn=%s", 100 * rate, format(n, big.mark = ","))),
            hjust = -0.07, size = 3.2) +
  scale_fill_viridis_c(option = "mako", direction = -1, guide = "none") +
  coord_flip() +
  labs(title = "Acceptance rate by old label",
       x = NULL, y = "Acceptance rate (%)") +
  ylim(0, 115) +
  theme_diag()

# 3c. Confidence distribution by transferred label (violin)
conf_acc <- results |>
  filter(transferred_label_accepted, !is.na(transferred_label_confidence))

p3c <- ggplot(conf_acc, aes(x = transferred_label, y = transferred_label_confidence,
                             fill = transferred_label)) +
  geom_violin(trim = TRUE, scale = "width", color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.12, color = "white", outlier.shape = NA, coef = 0) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Label confidence distribution (accepted)",
       x = NULL, y = "Confidence") +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# 3d. Worst SOM nodes by acceptance rate (nodes where < 50% of cells got a label)
node_accept <- results |>
  group_by(query_som_unit) |>
  summarise(n_cells = n(),
            rate = mean(transferred_label_accepted, na.rm = TRUE),
            .groups = "drop") |>
  left_join(
    data.frame(query_som_unit = node_df$query_node,
               correction_norm = node_df$correction_norm,
               match_fraction = node_df$match_fraction,
               top_ref_label = node_df$top_ref_label),
    by = "query_som_unit"
  )

worst_accept <- node_accept |>
  filter(n_cells >= 50) |>
  arrange(rate) |>
  head(20) |>
  mutate(node_label = sprintf("node %d (%s)\nn=%d", query_som_unit, top_ref_label, n_cells))

p3d <- ggplot(worst_accept, aes(x = reorder(node_label, rate), y = 100 * rate,
                                 fill = match_fraction)) +
  geom_col(color = "white") +
  scale_fill_viridis_c("Match\nfraction", option = "plasma", limits = c(0, 1)) +
  coord_flip() +
  labs(title = "20 nodes with lowest label acceptance (≥50 cells)",
       x = NULL, y = "Acceptance rate (%)") +
  theme_diag()

print((p3a | p3b) / (p3c | p3d) +
        plot_annotation(title = "PAGE 3 — Label transfer quality",
                        theme = theme(plot.title = element_text(face = "bold", size = 14))))

dev.off()

# =============================================================================
# Worst-performer summary tables (printed to console)
# =============================================================================
cat("\n--- 10 worst query nodes by match fraction ---\n")
print(node_df |>
  arrange(match_fraction) |>
  head(10) |>
  select(query_node, query_mass, transported_mass, match_fraction,
         correction_allowed, correction_norm, top_ref_label))

cat("\n--- 10 nodes with largest corrections (allowed) ---\n")
print(node_df |>
  filter(correction_allowed) |>
  arrange(desc(correction_norm)) |>
  head(10) |>
  select(query_node, query_mass, match_fraction, correction_norm, top_ref_label))

cat("\n--- Cell counts by final_status ---\n")
print(table(results$final_status))

cat(sprintf("\nPlots saved → %s\n", PDF_PATH))

cat("\n=== DONE ===\n")
