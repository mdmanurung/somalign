## Benchmark: ground-truth batch shift recovery
##
## Injects a known additive shift into reference data to simulate a batch
## effect, then runs the full somalign pipeline and measures how well the
## corrected projection recovers the original (pre-shift) reference node
## assignments.  Also demonstrates somalign_som_stability() to quantify
## run-to-run variance from query SOM training randomness.
##
## Run interactively or via:
##   Rscript benchmarks/shift_recovery.R

suppressPackageStartupMessages({
  library(somalign)
  library(kohonen)
})

set.seed(42)

## ── 1. Reference data: two well-separated clusters ───────────────────────────

n_per_cluster <- 200
p             <- 6      # features

ref_data <- rbind(
  matrix(rnorm(n_per_cluster * p, mean = -3, sd = 0.8), ncol = p),
  matrix(rnorm(n_per_cluster * p, mean =  3, sd = 0.8), ncol = p)
)
colnames(ref_data) <- paste0("F", seq_len(p))
ref_labels <- rep(c("neg", "pos"), each = n_per_cluster)

ref_grid <- somgrid(5, 5, "hexagonal")
reference <- somalign_train_reference(
  ref_data, labels = ref_labels, grid = ref_grid, rlen = 200
)
cat("Reference trained:", nrow(reference$codebook), "nodes\n")

## ── 2. Query data: same structure + known batch shift ────────────────────────

true_shift <- rep(1.5, p)   # additive shift on all features

query_data <- ref_data + matrix(true_shift, nrow = nrow(ref_data), ncol = p, byrow = TRUE)

## Ground truth: which reference node should each query sample map to?
##   After removing the shift, the query sample lies at the original position.
##   Direct projection of the un-shifted query == "true" reference assignment.
unshifted_scaled <- sweep(ref_data, 2, reference$center, "-")
unshifted_scaled <- sweep(unshifted_scaled, 2, reference$scale, "/")
nearest_unshifted <- somalign:::.somalign_nearest_code(unshifted_scaled, reference$codebook)
true_unit <- nearest_unshifted$unit

## ── 3. somalign alignment ─────────────────────────────────────────────────────

set.seed(123)
query <- somalign_query(query_data, reference, grid = somgrid(5, 5, "hexagonal"), rlen = 200)
fit   <- somalign_fit(query, reference, epsilon = 0.1, rho_query = 2, rho_ref = 2)

diag  <- somalign_diagnostics(fit)
res   <- somalign_results(fit)

cat(sprintf("\nSinkhorn converged: %s  (final_delta = %.2e)\n",
            diag$solver$converged, diag$solver$final_delta))
cat(sprintf("Transport mass: %.3f   mean match_fraction: %.3f\n",
            diag$ot$transport_mass,
            mean(diag$ot$match_fraction[is.finite(diag$ot$match_fraction)])))
cat(sprintf("rel_marginal_row_error: %.4f   rel_marginal_col_error: %.4f\n",
            diag$solver$rel_marginal_row_error,
            diag$solver$rel_marginal_col_error))

## ── 4. Recovery: compare direct vs corrected node assignment accuracy ─────────

direct_unit    <- res$old_som_unit
corrected_unit <- res$corrected_som_unit

direct_acc    <- mean(direct_unit    == true_unit)
corrected_acc <- mean(corrected_unit == true_unit)

cat(sprintf("\nNode assignment accuracy vs ground truth:\n"))
cat(sprintf("  Direct (no correction):  %.1f%%\n", 100 * direct_acc))
cat(sprintf("  Corrected (OT):          %.1f%%\n", 100 * corrected_acc))
cat(sprintf("  Improvement:             %+.1f pp\n", 100 * (corrected_acc - direct_acc)))

## ── 5. Mean correction norm vs true shift magnitude ──────────────────────────

true_shift_norm <- sqrt(sum(true_shift^2))   # in raw space
# Correction is applied in reference-scaled space; rescale for comparison
shift_scaled <- true_shift / reference$scale
true_shift_norm_scaled <- sqrt(sum(shift_scaled^2))

mean_corr_norm <- mean(res$correction_norm[is.finite(res$correction_norm)])
cat(sprintf("\nTrue shift ||Δ|| (scaled space): %.3f\n", true_shift_norm_scaled))
cat(sprintf("Mean correction_norm:             %.3f\n", mean_corr_norm))
cat(sprintf("Ratio (correction / true shift):  %.3f\n",
            mean_corr_norm / max(true_shift_norm_scaled, 1e-10)))

## ── 6. Sensitivity grid ───────────────────────────────────────────────────────

cat("\n── Sensitivity grid (epsilon × rho) ──\n")
grid_result <- somalign_sensitivity_grid(
  query, reference,
  epsilon   = c(0.1, 0.5, 1.0),
  rho_query = c(0.5, 2.0),
  rho_ref   = c(0.5, 2.0)
)
print(grid_result[, c("epsilon", "rho_query", "rho_ref",
                       "transport_mass", "mean_match_fraction",
                       "outside_direct_fraction", "outside_corrected_fraction")])

## ── 7. SOM seed stability ─────────────────────────────────────────────────────

cat("\n── Query SOM seed stability (10 seeds) ──\n")
stability <- somalign_som_stability(
  query_data, reference,
  som_seeds = 1:10,
  epsilon   = 0.1,
  rho_query = 2,
  rho_ref   = 2,
  grid      = somgrid(5, 5, "hexagonal"),
  rlen      = 200
)
print(stability)
cat(sprintf(
  "\ntransport_mass  range: [%.4f, %.4f]  SD: %.4f\n",
  min(stability$transport_mass), max(stability$transport_mass),
  sd(stability$transport_mass)
))
cat(sprintf(
  "mean_correction_norm  range: [%.4f, %.4f]  SD: %.4f\n",
  min(stability$mean_correction_norm, na.rm = TRUE),
  max(stability$mean_correction_norm, na.rm = TRUE),
  sd(stability$mean_correction_norm, na.rm = TRUE)
))

## ── 8. Log-domain solver: verify consistency ──────────────────────────────────

set.seed(123)
query2   <- somalign_query(query_data, reference, grid = somgrid(5, 5, "hexagonal"), rlen = 200)
fit_log  <- somalign_fit(query2, reference, epsilon = 0.1, rho_query = 2, rho_ref = 2,
                          solver = "log_domain")
res_log  <- somalign_results(fit_log)

corr_acc_log <- mean(res_log$corrected_som_unit == true_unit)
cat(sprintf("\nLog-domain solver corrected accuracy: %.1f%%\n", 100 * corr_acc_log))
cat(sprintf("Plan difference from internal solver: max|P_log - P_int| = %.2e\n",
            max(abs(fit_log$transport_plan - fit$transport_plan))))
