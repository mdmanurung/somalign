#!/usr/bin/env Rscript
# =============================================================================
# Audit the true repeat-sample frequency-correlation evidence from the SOMalign
# example notebook tree.
#
# The labelled Batch1-vs-Batch2 pilot used by other_t_improvement_experiment.R
# does not contain the notebook's repeat-pair metadata. The notebook evidence is
# a separate pilot-A versus BMV/query-B design joined by repeated_samples.csv.
# This script records that source metadata and summarizes the saved hard-vs-soft
# CLR reproducibility artifacts without recomputing the 40M-cell projection.
# =============================================================================

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
SCRIPT <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
} else {
  NA_character_
}
SCRIPT_DIR <- if (is.na(SCRIPT)) getwd() else dirname(SCRIPT)

PROJECT_DIR <- Sys.getenv(
  "SOMALIGN_NOTEBOOK_PROJECT_DIR",
  "/exports/para-lipg-hpc/mdmanurung/bmv_pilot_cytof_integration"
)
SOMALIGN_DIR <- Sys.getenv(
  "SOMALIGN_NOTEBOOK_SOMALIGN_DIR",
  file.path(PROJECT_DIR, "data/processed/dataset1/somalign")
)
REPEAT_CSV <- Sys.getenv(
  "SOMALIGN_REPEAT_CSV",
  file.path(PROJECT_DIR, "repeated_samples.csv")
)

write_table <- function(x, name) {
  utils::write.csv(x, file.path(SCRIPT_DIR, name), row.names = FALSE)
}

read_required_rds <- function(name) {
  path <- file.path(SOMALIGN_DIR, name)
  if (!file.exists(path)) {
    stop("Required artifact is missing: ", path, call. = FALSE)
  }
  readRDS(path)
}

stamp <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

stamp("reading repeat-pair metadata")
rep <- utils::read.csv(REPEAT_CSV, stringsAsFactors = FALSE, check.names = FALSE)
names(rep)[names(rep) == ""] <- "row_id"
required <- c("fcs_filename", "batch", "sample_id")
missing <- setdiff(required, names(rep))
if (length(missing)) {
  stop("Repeat table is missing required columns: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

batch_counts <- table(rep$sample_id, rep$batch)
if (!all(c("pilot", "new") %in% colnames(batch_counts))) {
  stop("Repeat table must contain both 'pilot' and 'new' batches.", call. = FALSE)
}
shared <- rownames(batch_counts)[batch_counts[, "pilot"] > 0 & batch_counts[, "new"] > 0]
if (!length(shared)) {
  stop("No sample_id has both pilot and new repeat files.", call. = FALSE)
}

pilot_map <- with(rep[rep$batch == "pilot", ], setNames(fcs_filename, sample_id))
new_map <- split(rep$fcs_filename[rep$batch == "new"], rep$sample_id[rep$batch == "new"])

stamp("reading dedicated soft-grid artifact")
soft_grid <- read_required_rds("clr_soft_projection.rds")
artifact_samples <- soft_grid$samples
if (is.null(artifact_samples)) {
  artifact_samples <- shared
}
if (!setequal(artifact_samples, shared)) {
  stop("Repeat metadata sample IDs do not match clr_soft_projection.rds samples.",
       call. = FALSE)
}

pair_table <- data.frame(
  sample_id = artifact_samples,
  n_pilot_files = as.integer(batch_counts[artifact_samples, "pilot"]),
  n_new_files = as.integer(batch_counts[artifact_samples, "new"]),
  pilot_fcs_filename = unname(pilot_map[artifact_samples]),
  new_fcs_filename = vapply(new_map[artifact_samples], paste, character(1), collapse = " | "),
  stringsAsFactors = FALSE
)
write_table(pair_table, "notebook_repeat_sample_metadata.csv")

stamp("reading projection-quality report artifact")
pqr <- read_required_rds("projection_quality_report.rds")
r <- pqr$clr
if (is.null(r$meta_hard_sample) || is.null(r$meta_soft_sample) ||
    is.null(r$lin_hard_sample) || is.null(r$lin_soft_sample)) {
  stop("projection_quality_report.rds does not contain the expected clr vectors.",
       call. = FALSE)
}

summary_rows <- data.frame(
  source = "projection_quality_report.rds",
  evidence_class = "true_repeat_metadata",
  resolution = c("metacluster", "lineage"),
  hard_method = "hard nearest-node",
  soft_method = "somalign_soft_frequencies k=8",
  n_repeat_samples = length(artifact_samples),
  hard_median_clr_weighted_r = c(median(r$meta_hard_sample, na.rm = TRUE),
                                 median(r$lin_hard_sample, na.rm = TRUE)),
  soft_median_clr_weighted_r = c(median(r$meta_soft_sample, na.rm = TRUE),
                                 median(r$lin_soft_sample, na.rm = TRUE)),
  median_difference_delta = c(
    median(r$meta_soft_sample, na.rm = TRUE) - median(r$meta_hard_sample, na.rm = TRUE),
    median(r$lin_soft_sample, na.rm = TRUE) - median(r$lin_hard_sample, na.rm = TRUE)
  ),
  median_pairwise_delta = c(median(r$meta_soft_sample - r$meta_hard_sample, na.rm = TRUE),
                            median(r$lin_soft_sample - r$lin_hard_sample, na.rm = TRUE)),
  samples_soft_gt_hard = c(sum(r$meta_soft_sample > r$meta_hard_sample, na.rm = TRUE),
                           sum(r$lin_soft_sample > r$lin_hard_sample, na.rm = TRUE)),
  split_half_ceiling_median = c(median(r$ceiling_sample, na.rm = TRUE), NA_real_),
  stringsAsFactors = FALSE
)
write_table(summary_rows, "notebook_repeat_correlation_summary.csv")

pair_vectors <- data.frame(
  sample_id = artifact_samples,
  metacluster_hard = as.numeric(r$meta_hard_sample),
  metacluster_soft_k8 = as.numeric(r$meta_soft_sample),
  metacluster_delta = as.numeric(r$meta_soft_sample - r$meta_hard_sample),
  lineage_hard = as.numeric(r$lin_hard_sample),
  lineage_soft_k8 = as.numeric(r$lin_soft_sample),
  lineage_delta = as.numeric(r$lin_soft_sample - r$lin_hard_sample),
  evidence_class = "true_repeat_metadata",
  stringsAsFactors = FALSE
)
write_table(pair_vectors, "notebook_repeat_correlation_pairs.csv")
soft_grid_summary <- soft_grid$results
soft_grid_summary$source <- "clr_soft_projection.rds"
soft_grid_summary$evidence_class <- "true_repeat_metadata"
write_table(soft_grid_summary, "notebook_repeat_soft_grid_summary.csv")

source_artifacts <- data.frame(
  artifact = c(
    "repeated_samples.csv",
    "projection_quality_report.rds",
    "clr_soft_projection.rds",
    "projection_quality_report.md",
    "clr_soft_projection.R",
    "projection_quality_report.R"
  ),
  path = c(
    REPEAT_CSV,
    file.path(SOMALIGN_DIR, "projection_quality_report.rds"),
    file.path(SOMALIGN_DIR, "clr_soft_projection.rds"),
    file.path(PROJECT_DIR, "nbs/projection_quality_report.md"),
    file.path(PROJECT_DIR, "nbs/clr_soft_projection.R"),
    file.path(PROJECT_DIR, "nbs/projection_quality_report.R")
  ),
  exists = file.exists(c(
    REPEAT_CSV,
    file.path(SOMALIGN_DIR, "projection_quality_report.rds"),
    file.path(SOMALIGN_DIR, "clr_soft_projection.rds"),
    file.path(PROJECT_DIR, "nbs/projection_quality_report.md"),
    file.path(PROJECT_DIR, "nbs/clr_soft_projection.R"),
    file.path(PROJECT_DIR, "nbs/projection_quality_report.R")
  )),
  stringsAsFactors = FALSE
)
write_table(source_artifacts, "notebook_repeat_source_artifacts.csv")

cat("\n===== notebook true-repeat correlation audit =====\n")
print(summary_rows, row.names = FALSE, digits = 4)
cat(sprintf("\nRepeat table: %d paired sample IDs, %d rows (%d pilot, %d new)\n",
            length(artifact_samples), nrow(rep), sum(rep$batch == "pilot"), sum(rep$batch == "new")))
cat(sprintf("Wrote audit CSVs -> %s\n", SCRIPT_DIR))
stamp("complete")
