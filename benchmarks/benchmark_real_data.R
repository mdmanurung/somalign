## somalign — Real-data batch correction benchmark
##
## Evaluates somalign alignment quality on two public cytometry datasets:
##
##   1. Nuñez 2023 full-spectrum flow PBMC
##      Same donor measured on two consecutive days — pure technical batch
##      effect, ground truth = identical biology.
##      Source: exampledata.scverse.org (downloaded to data/Nunez2023/)
##
##   2. CytoNorm CyTOF example (Gates et al.)
##      3 patients × 2 staining batches (_1 = train, _2 = validation).
##      Same patients in both batches → natural anchor samples for
##      somalign_fit_anchored().
##      Source: system.file("extdata", package = "CytoNorm")
##
## Batch correction metric: corrected-JS, the JS divergence between the
##   reference node-occupancy distribution and the corrected_som_unit
##   distribution.  Lower = better.  The direct projection (old_som_unit)
##   always equals the baseline by construction — both assign cells to the
##   nearest reference codebook node — so it is reported only as a sanity check,
##   not as a somalign performance metric.
##
## Run with:
##   Rscript benchmarks/benchmark_real_data.R 2>&1 | tee benchmarks/real_data_output.log

suppressPackageStartupMessages({
  library(devtools)
  load_all(quiet = TRUE)
  library(kohonen)
  library(flowCore)
})

## ── Helpers ──────────────────────────────────────────────────────────────────

source("benchmarks/helpers.R")   # sep, js_div, node_dist, ks_per_marker

## ── 1. Nuñez 2023 — Full-Spectrum Flow ───────────────────────────────────────

sep("Nuñez 2023 — Full-Spectrum Flow PBMC")

fcs_dir  <- file.path("data", "Nunez2023")
fcs1     <- file.path(fcs_dir, "Nunez_PBMCs_batch1.fcs")
fcs2     <- file.path(fcs_dir, "Nunez_PBMCs_batch2.fcs")

if (!file.exists(fcs1) || !file.exists(fcs2)) {
  cat("SKIP: Nuñez FCS files not found at", fcs_dir, "\n")
  cat("      Download from exampledata.scverse.org/cytovi/\n")
} else {

  ## Load and select markers (FJComp-* channels, skip viability)
  batch1_ff <- read.FCS(fcs1, truncate_max_range = FALSE)
  batch2_ff <- read.FCS(fcs2, truncate_max_range = FALSE)

  all_params <- pData(parameters(batch1_ff))$name
  markers_fl <- grep("^FJComp-", all_params, value = TRUE)
  markers_fl <- setdiff(markers_fl, "FJComp-Zombie NIR-A")   # viability

  ## arcsinh transform: cofactor 2000 (CytoVI recommendation for this dataset)
  cofactor_fl <- 2000
  batch1 <- asinh(exprs(batch1_ff)[, markers_fl] / cofactor_fl)
  batch2 <- asinh(exprs(batch2_ff)[, markers_fl] / cofactor_fl)
  colnames(batch1) <- colnames(batch2) <- markers_fl
  cat(sprintf("Cells: batch1=%d  batch2=%d  markers=%d  (cofactor=%d)\n",
              nrow(batch1), nrow(batch2), length(markers_fl), cofactor_fl))

  ## Reference SOM — trained on batch1
  cat("Training reference SOM (10×10, rlen=100)...\n")
  set.seed(42)
  ref_fl <- somalign_train_reference(
    batch1,
    grid = somgrid(10, 10, "hexagonal"),
    rlen = 100
  )
  n_nodes_fl  <- nrow(ref_fl$codebook)
  ref_fl_dist <- node_dist(ref_fl$som$unit.classif, n_nodes_fl)

  ## Baseline: project batch2 directly onto reference SOM via kohonen::map
  ## (no somalign alignment — pure batch effect)
  batch2_scaled    <- sweep(sweep(batch2, 2, ref_fl$center, "-"), 2, ref_fl$scale, "/")
  baseline_units   <- kohonen::map(ref_fl$som, batch2_scaled)$unit.classif
  baseline_fl_dist <- node_dist(baseline_units, n_nodes_fl)

  ## somalign alignment
  cat("Training query SOM and fitting somalign...\n")
  set.seed(123)
  query_fl <- somalign_query(
    batch2, ref_fl,
    grid = somgrid(10, 10, "hexagonal"),
    rlen = 100
  )
  fit_fl  <- somalign_fit(query_fl, ref_fl, epsilon = 0.1, rho_query = 2, rho_ref = 2)
  res_fl  <- somalign_results(fit_fl)
  diag_fl <- somalign_diagnostics(fit_fl)

  direct_fl_dist    <- node_dist(res_fl$old_som_unit,       n_nodes_fl)
  corrected_fl_dist <- node_dist(res_fl$corrected_som_unit, n_nodes_fl)

  ## Per-marker KS statistic (batch effect size, independent of somalign)
  ks_fl <- suppressWarnings(ks_per_marker(batch1, batch2))
  names(ks_fl) <- markers_fl
  cat(sprintf("Per-marker KS (batch1 vs batch2): mean=%.3f  median=%.3f  max=%.3f\n",
              mean(ks_fl), median(ks_fl), max(ks_fl)))
  cat("Top 5 markers by KS statistic:\n")
  top5_fl <- sort(ks_fl, decreasing = TRUE)[seq_len(min(5L, length(ks_fl)))]
  for (nm in names(top5_fl)) cat(sprintf("  %-42s %.3f\n", nm, top5_fl[nm]))

  sep("Nuñez 2023 — Alignment quality")
  cat(sprintf("Solver converged: %s   final_delta: %.2e\n",
              diag_fl$solver$converged, diag_fl$solver$final_delta))
  cat(sprintf("Transport mass: %.3f   mean match_fraction: %.3f\n",
              diag_fl$ot$transport_mass,
              mean(diag_fl$ot$match_fraction[is.finite(diag_fl$ot$match_fraction)])))

  ## old_som_unit == baseline by construction (both assign to nearest ref node).
  ## Direct-JS is reported only as a sanity check; corrected-JS is the metric.
  cat("\nJS divergence vs reference node distribution:\n")
  cat(sprintf("  Baseline / direct (old_som_unit, by construction equal): %.4f\n",
              js_div(ref_fl_dist, baseline_fl_dist)))
  cat(sprintf("  somalign corrected (corrected_som_unit):                 %.4f\n",
              js_div(ref_fl_dist, corrected_fl_dist)))

  ## Gate B diagnostic: corrected-JS vs epsilon
  sep("Nuñez 2023 — Corrected-JS vs epsilon (Gate B)")
  cat("Metric: JS(reference, corrected_som_unit_distribution).  Lower = better.\n")
  cat("direct_JS == baseline by construction for all epsilon rows.\n")
  cat(sprintf("  %-8s  %-12s  %-12s  %-10s  %-10s\n",
              "epsilon", "corrected_JS", "direct_JS", "converged", "outside_corr%"))
  baseline_js <- js_div(ref_fl_dist, baseline_fl_dist)
  for (eps in c(0.05, 0.1, 0.2, 0.5)) {
    fit_eps  <- somalign_fit(query_fl, ref_fl, epsilon = eps, rho_query = 2, rho_ref = 2,
                             solver = "log_domain", max_iter = 2000)
    diag_eps <- somalign_diagnostics(fit_eps)
    cor_dist_eps <- node_dist(fit_eps$projection$corrected$unit, n_nodes_fl)
    cat(sprintf("  %-8.3f  %-12.4f  %-12.4f  %-10s  %-10.1f\n",
                eps,
                js_div(ref_fl_dist, cor_dist_eps),
                baseline_js,
                as.character(diag_eps$solver$converged),
                100 * diag_eps$projection$outside_corrected_fraction))
  }

  ## Sensitivity grid (epsilon × rho)
  sep("Nuñez 2023 — Sensitivity grid")
  sg_fl <- somalign_sensitivity_grid(
    query_fl, ref_fl,
    epsilon   = c(0.1, 0.5, 1.0),
    rho_query = c(0.5, 2.0),
    rho_ref   = 2.0
  )
  print(sg_fl[, c("epsilon", "rho_query", "transport_mass",
                   "mean_match_fraction", "outside_direct_fraction")])

  ## SOM seed stability (5 seeds)
  sep("Nuñez 2023 — SOM seed stability")
  cat("Re-training query SOM for 5 seeds...\n")
  stab_fl <- somalign_som_stability(
    batch2, ref_fl,
    som_seeds = 1:5,
    grid      = somgrid(10, 10, "hexagonal"),
    rlen      = 100,
    epsilon   = 0.1,
    rho_query = 2,
    rho_ref   = 2
  )
  print(stab_fl[, c("som_seed", "transport_mass", "mean_match_fraction",
                     "mean_correction_norm", "converged")])
  cat(sprintf("transport_mass SD: %.4f   mean_correction_norm SD: %.4f\n",
              sd(stab_fl$transport_mass),
              sd(stab_fl$mean_correction_norm, na.rm = TRUE)))
}

## ── 2. CytoNorm CyTOF example — Anchor comparison ───────────────────────────

sep("CytoNorm CyTOF — 3 patients × 2 staining batches")

cytonorm_ok <- requireNamespace("CytoNorm", quietly = TRUE)
if (!cytonorm_ok) {
  cat("SKIP: CytoNorm not installed.\n")
  cat("      remotes::install_github('saeyslab/CytoNorm')\n")
} else {
  cyto_dir   <- system.file("extdata", package = "CytoNorm")
  all_files  <- list.files(cyto_dir, pattern = "\\.fcs$", full.names = TRUE)

  ## _1 = Train/reference batch;  _2 = Validation/query batch
  ref_files   <- grep("_1\\.fcs$", all_files, value = TRUE)
  query_files <- grep("_2\\.fcs$", all_files, value = TRUE)
  # Match patient order so anchor_old[i] and anchor_new[i] are same patient
  ref_files   <- sort(ref_files)
  query_files <- sort(query_files)
  patient_ids <- sub(".*PTLG([0-9]+).*", "PTLG\\1", basename(ref_files))

  cat("Reference batch files (_1):", basename(ref_files), "\n")
  cat("Query batch files     (_2):", basename(query_files), "\n")

  ## Select CyTOF antibody channels: Di channels with non-NA desc,
  ## excluding barcodes (Pd*), DNA (Ir191Di/Ir193Di), beadDist, Time,
  ## Event_length, and Y89Di / I127Di / Ba138Di / Ce140Di / Pt195Di (no desc)
  ref_ff0 <- read.FCS(ref_files[1], truncate_max_range = FALSE)
  pd      <- pData(parameters(ref_ff0))
  is_di   <- grepl("Di$", pd$name)
  has_desc <- !is.na(pd$desc)
  is_bc   <- grepl("^BC[0-9]", pd$desc) & !is.na(pd$desc)
  is_dna  <- pd$name %in% c("Ir191Di", "Ir193Di")
  is_bead <- pd$name == "beadDist"
  markers_cyto <- pd$name[is_di & has_desc & !is_bc & !is_dna & !is_bead]
  cat(sprintf("Antibody channels selected: %d\n", length(markers_cyto)))

  ## Load and arcsinh transform (cofactor 5 for CyTOF)
  cofactor_cyto <- 5
  read_fcs_matrix <- function(files) {
    mats <- lapply(files, function(f) {
      ff  <- read.FCS(f, truncate_max_range = FALSE)
      asinh(exprs(ff)[, markers_cyto] / cofactor_cyto)
    })
    do.call(rbind, mats)
  }

  cat("Loading batch1 (_1) files...\n")
  batch1_cyto <- read_fcs_matrix(ref_files)
  cat("Loading batch2 (_2) files...\n")
  batch2_cyto <- read_fcs_matrix(query_files)
  colnames(batch1_cyto) <- colnames(batch2_cyto) <- markers_cyto
  cat(sprintf("Cells: batch1=%d  batch2=%d  markers=%d  (cofactor=%d)\n",
              nrow(batch1_cyto), nrow(batch2_cyto),
              length(markers_cyto), cofactor_cyto))

  ## Per-marker KS statistic
  ks_cyto <- suppressWarnings(ks_per_marker(batch1_cyto, batch2_cyto))
  names(ks_cyto) <- markers_cyto
  cat(sprintf("Per-marker KS (batch1 vs batch2): mean=%.3f  median=%.3f  max=%.3f\n",
              mean(ks_cyto), median(ks_cyto), max(ks_cyto)))

  ## Reference SOM — trained on batch1 (_1 files combined)
  cat("Training reference SOM (6×6, rlen=100)...\n")
  set.seed(42)
  ref_cyto <- somalign_train_reference(
    batch1_cyto,
    grid = somgrid(6, 6, "hexagonal"),
    rlen = 100
  )
  n_nodes_cyto  <- nrow(ref_cyto$codebook)
  ref_cyto_dist <- node_dist(ref_cyto$som$unit.classif, n_nodes_cyto)

  ## Baseline: direct projection of batch2 onto reference SOM
  batch2c_scaled    <- sweep(sweep(batch2_cyto, 2, ref_cyto$center, "-"), 2, ref_cyto$scale, "/")
  baseline_cyto_units <- kohonen::map(ref_cyto$som, batch2c_scaled)$unit.classif
  baseline_cyto_dist  <- node_dist(baseline_cyto_units, n_nodes_cyto)

  ## Query SOM
  cat("Training query SOM...\n")
  set.seed(123)
  query_cyto <- somalign_query(
    batch2_cyto, ref_cyto,
    grid = somgrid(6, 6, "hexagonal"),
    rlen = 100
  )

  ## Standard fit (no anchor information)
  cat("Fitting somalign_fit()...\n")
  fit_plain  <- somalign_fit(query_cyto, ref_cyto, epsilon = 0.1, rho_query = 2, rho_ref = 2)
  res_plain  <- somalign_results(fit_plain)
  diag_plain <- somalign_diagnostics(fit_plain)

  ## Anchor pairs: same patients in both batches
  ## Each patient has 1000 cells in each batch.  Stack patients in the same
  ## order so anchor_old[i,] and anchor_new[i,] are from the same patient.
  ## (Not 1:1 cell-level pairing, but same population in old vs new batch.)
  n_per_file  <- nrow(batch1_cyto) / length(ref_files)   # 1000
  anchor_old  <- batch1_cyto
  anchor_new  <- batch2_cyto

  cat("Fitting somalign_fit_anchored()...\n")
  fit_anc  <- somalign_fit_anchored(
    query_cyto, ref_cyto,
    anchor_old = anchor_old,
    anchor_new = anchor_new,
    rho_anchor = 1.0,
    epsilon    = 0.1,
    rho_query  = 2,
    rho_ref    = 2
  )
  res_anc  <- somalign_results(fit_anc)
  diag_anc <- somalign_diagnostics(fit_anc)

  ## Node occupancy distributions
  plain_direct_dist    <- node_dist(res_plain$old_som_unit,       n_nodes_cyto)
  plain_corrected_dist <- node_dist(res_plain$corrected_som_unit, n_nodes_cyto)
  anc_direct_dist      <- node_dist(res_anc$old_som_unit,         n_nodes_cyto)
  anc_corrected_dist   <- node_dist(res_anc$corrected_som_unit,   n_nodes_cyto)

  sep("CytoNorm CyTOF — Standard vs anchored fit")
  cat(sprintf("Standard   converged: %s   transport_mass: %.3f\n",
              diag_plain$solver$converged, diag_plain$ot$transport_mass))
  cat(sprintf("Anchored   converged: %s   transport_mass: %.3f\n",
              diag_anc$solver$converged,  diag_anc$ot$transport_mass))
  cat(sprintf("Anchor coverage: %d nodes / %.0f%%\n",
              fit_anc$anchors$nodes_covered,
              100 * fit_anc$anchors$coverage_fraction))

  ## direct (old_som_unit) == baseline by construction; corrected_som_unit is the metric.
  cat("\nJS divergence vs reference node distribution:\n")
  cat(sprintf("  Baseline / direct (old_som_unit, equal by construction): %.4f\n",
              js_div(ref_cyto_dist, baseline_cyto_dist)))
  cat(sprintf("  Standard  corrected (corrected_som_unit):                %.4f\n",
              js_div(ref_cyto_dist, plain_corrected_dist)))
  cat(sprintf("  Anchored  corrected (corrected_som_unit):                %.4f\n",
              js_div(ref_cyto_dist, anc_corrected_dist)))

  ## Correction vector comparison
  plain_norms <- sqrt(rowSums(fit_plain$node_shifts^2))
  anc_norms   <- sqrt(rowSums(fit_anc$node_shifts^2))
  cat(sprintf("\nMean correction norm — standard: %.4f   anchored: %.4f\n",
              mean(plain_norms), mean(anc_norms)))

  ## rho_anchor sweep
  sep("CytoNorm CyTOF — rho_anchor sweep")
  rhos <- c(0, 0.5, 1.0, 1.5, 2.0)
  rho_results <- lapply(rhos, function(rho) {
    f    <- somalign_fit_anchored(
      query_cyto, ref_cyto,
      anchor_old = anchor_old, anchor_new = anchor_new,
      rho_anchor = rho, epsilon = 0.1, rho_query = 2, rho_ref = 2
    )
    dcor <- node_dist(f$projection$corrected$unit, n_nodes_cyto)
    data.frame(
      rho_anchor     = rho,
      transport_mass = round(f$diagnostics$ot$transport_mass, 5),
      js_corrected   = round(js_div(ref_cyto_dist, dcor), 5),
      mean_corr_norm = round(mean(sqrt(rowSums(f$node_shifts^2))), 5)
    )
  })
  print(do.call(rbind, rho_results))

  ## SOM seed stability
  sep("CytoNorm CyTOF — SOM seed stability")
  stab_cyto <- somalign_som_stability(
    batch2_cyto, ref_cyto,
    som_seeds = 1:5,
    grid      = somgrid(6, 6, "hexagonal"),
    rlen      = 100,
    epsilon   = 0.1,
    rho_query = 2,
    rho_ref   = 2
  )
  print(stab_cyto[, c("som_seed", "transport_mass", "mean_match_fraction",
                       "mean_correction_norm", "converged")])
  cat(sprintf("transport_mass SD: %.4f   mean_correction_norm SD: %.4f\n",
              sd(stab_cyto$transport_mass),
              sd(stab_cyto$mean_correction_norm, na.rm = TRUE)))
}

cat("\n── Done ────────────────────────────────────────────────────────\n")
