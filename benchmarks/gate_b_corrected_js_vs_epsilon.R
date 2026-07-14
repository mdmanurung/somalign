## Gate B: corrected-JS vs epsilon curve (Nuñez 2023)
##
## Adds the one measurement missing from benchmark_real_data.R:
## JS-divergence of the corrected projection across epsilon values.
##
## If corrected-JS → baseline as epsilon → small  → shrinkage confirmed
## If corrected-JS stays bad at small epsilon      → re-projection bug
##
## Run with:
##   Rscript benchmarks/gate_b_corrected_js_vs_epsilon.R 2>&1 \
##     | tee benchmarks/gate_b_output.log

suppressPackageStartupMessages({
  library(somalign)
  library(kohonen)
  library(flowCore)
})

## ── Helpers ──────────────────────────────────────────────────────────────────

source("benchmarks/helpers.R")   # js_div, node_dist

## ── Load Nuñez 2023 data ─────────────────────────────────────────────────────

fcs_dir <- file.path("data", "Nunez2023")
fcs1    <- file.path(fcs_dir, "Nunez_PBMCs_batch1.fcs")
fcs2    <- file.path(fcs_dir, "Nunez_PBMCs_batch2.fcs")

if (!file.exists(fcs1) || !file.exists(fcs2)) {
  stop("Nuñez FCS files not found at ", fcs_dir,
       "\nDownload from exampledata.scverse.org/cytovi/")
}

batch1_ff <- read.FCS(fcs1, truncate_max_range = FALSE)
batch2_ff <- read.FCS(fcs2, truncate_max_range = FALSE)

all_params  <- pData(parameters(batch1_ff))$name
markers_fl  <- grep("^FJComp-", all_params, value = TRUE)
markers_fl  <- setdiff(markers_fl, "FJComp-Zombie NIR-A")
cofactor_fl <- 2000

batch1 <- asinh(exprs(batch1_ff)[, markers_fl] / cofactor_fl)
batch2 <- asinh(exprs(batch2_ff)[, markers_fl] / cofactor_fl)
colnames(batch1) <- colnames(batch2) <- markers_fl
cat(sprintf("Cells: batch1=%d  batch2=%d  markers=%d\n",
            nrow(batch1), nrow(batch2), length(markers_fl)))

## ── Train reference SOM ───────────────────────────────────────────────────────

cat("Training reference SOM (10×10, rlen=100)...\n")
set.seed(42)
ref_fl      <- somalign_train_reference(batch1, grid = somgrid(10, 10, "hexagonal"), rlen = 100)
n_nodes_fl  <- nrow(ref_fl$codebook)
ref_fl_dist <- node_dist(ref_fl$som$unit.classif, n_nodes_fl)

## Baseline (identical to direct by construction — kept as sanity reference)
batch2_scaled  <- sweep(sweep(batch2, 2, ref_fl$center, "-"), 2, ref_fl$scale, "/")
baseline_units <- kohonen::map(ref_fl$som, batch2_scaled)$unit.classif
baseline_js    <- js_div(ref_fl_dist, node_dist(baseline_units, n_nodes_fl))
cat(sprintf("Baseline JS (direct NN, no alignment): %.4f\n", baseline_js))

## ── Train query SOM once ─────────────────────────────────────────────────────

cat("Training query SOM (10×10, rlen=100)...\n")
set.seed(123)
query_fl <- somalign_query(batch2, ref_fl, grid = somgrid(10, 10, "hexagonal"), rlen = 100)

## ── Gate B sweep ─────────────────────────────────────────────────────────────

epsilons <- c(0.05, 0.1, 0.2, 0.5)
rho      <- 2.0   # same as benchmark_real_data.R

cat("\n── corrected-JS vs epsilon (rho_query=rho_ref=", rho, ") ──\n", sep = "")
cat(sprintf("%-8s  %-12s  %-12s  %-8s  %-8s  %-12s\n",
            "epsilon", "corrected_JS", "direct_JS", "conv", "iters", "outside_corr%"))
cat(strrep("-", 72), "\n")

for (eps in epsilons) {
  fit <- somalign_fit(query_fl, ref_fl,
                      epsilon   = eps,
                      rho_query = rho,
                      rho_ref   = rho,
                      solver    = "log_domain",   # stable at small epsilon
                      max_iter  = 2000)
  res  <- somalign_results(fit)
  diag <- somalign_diagnostics(fit)

  cor_dist <- node_dist(res$corrected_som_unit, n_nodes_fl)
  dir_dist <- node_dist(res$old_som_unit,       n_nodes_fl)

  cor_js   <- js_div(ref_fl_dist, cor_dist)
  dir_js   <- js_div(ref_fl_dist, dir_dist)

  cat(sprintf("%-8.3f  %-12.4f  %-12.4f  %-8s  %-8d  %-12.1f\n",
              eps,
              cor_js,
              dir_js,
              ifelse(diag$solver$converged, "yes", "no"),
              diag$solver$iterations,
              100 * diag$projection$outside_corrected_fraction))
}

cat("\nConclusion: if corrected_JS decreases steadily as epsilon decreases,\n")
cat("            shrinkage is confirmed (the fix is lowering default epsilon\n")
cat("            and/or using a debiased barycenter).\n")
cat("            If corrected_JS stays high at eps=0.05, investigate\n")
cat("            .somalign_node_shifts / .somalign_project_pair for a bug.\n")
cat("\n── Done ────────────────────────────────────────────────────────\n")
