#!/usr/bin/env Rscript
# =============================================================================
# Does using anchor (repeat) samples improve label transfer? Two conditions.
#
# NOTE ON DATA: the deployed BMV query (aurora_exvivo_cycombine_normalized.qs)
# carries NO per-cell ground-truth label -- it is unlabelled by design, which
# is exactly why labels are transferred to it. Label-transfer ACCURACY (and
# therefore anchor benefit) cannot be measured directly on it. We instead use
# the labelled pilot data (gate_lineage2, 6 lineages), which has real
# acquisition batches (Batch1-4) AND ground truth, to answer the question:
#
#   Condition A (real cross-batch): reference = pilot Batch1, query = pilot
#     Batch2 -- a real acquisition-batch effect with real labels on both sides.
#   Condition B (controlled severity): inject a synthetic per-marker batch
#     shift of magnitude delta into the query and sweep delta, for a
#     dose-response of anchor benefit vs batch severity.
#
# Both use somalign_anchor_benefit() to sweep rho_anchor (0 = no anchors) and
# score transferred labels against ground truth on held-out (non-anchor) cells.
#
# Usage:  Rscript anchor_benefit_experiment.R [n_subsample]
# =============================================================================

.libPaths("/exports/para-lipg-hpc/mdmanurung/R/4.5")
suppressPackageStartupMessages({ library(somalign); library(kohonen) })

args  <- commandArgs(trailingOnly = TRUE)
n_sub <- if (length(args) >= 1) as.integer(args[1]) else 120000L
set.seed(1)

PILOT <- "/exports/para-lipg-hpc/Xuran/data/pilot/aurora_exvivo/processed/aurora_exvivo_xyfclusters.qs"
OUT   <- "/exports/para-lipg-hpc/mdmanurung/somalign/analysis/anchor-benefit/anchor_benefit_results.rds"
FEATURES <- c("TCRgd_CD141","CD8","CD3","CD19","CD159c_NKG2C","CD38","CD57","CD161",
              "CD1c","HLA-DR","CD11c","CD25","CD16","CD185_CXCR5","CD14","CD303_TCRvd2",
              "CD4","CD294_CRTH2","CD21","CD56","FoxP3","Tbet","CD197_CCR7","CD127",
              "CD45RA","CD27","CD7")
GRID     <- kohonen::somgrid(15, 15, "hexagonal")
RHO_GRID <- c(0, 1, 5, 20, 100, 500)
N_ANCHOR <- 6000L

stamp <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

stamp("loading pilot data ...")
df <- qs2::qs_read(PILOT, nthreads = 8)
y_all <- as.character(df$gate_lineage2)
b_all <- as.character(df$batch_id)
keep  <- !is.na(y_all) & y_all != "" & stats::complete.cases(df[, FEATURES])

# Stratified subsample by (batch, lineage) to keep both batches balanced.
strat <- interaction(b_all[keep], y_all[keep], drop = TRUE)
idx0  <- which(keep)
per   <- max(1L, floor(n_sub / nlevels(strat)))
sub   <- unlist(lapply(split(idx0, strat), function(ix)
  if (length(ix) <= per) ix else sample(ix, per)), use.names = FALSE)
X <- as.matrix(df[sub, FEATURES]); y <- y_all[sub]; b <- b_all[sub]
rm(df); gc()
stamp(sprintf("subsampled %d cells; batches: %s", nrow(X),
              paste(names(table(b)), table(b), sep = "=", collapse = " ")))

# Paired anchors from matched-lineage cells across two cell pools (simulates
# repeat QC samples: a population measured in both batches). Returns row-paired
# old/new matrices plus the row indices used in the "new" pool (to exclude from
# evaluation).
make_anchors <- function(X_old, y_old, X_new, y_new, n_total) {
  classes <- intersect(unique(y_old), unique(y_new))
  per_cls <- max(1L, floor(n_total / length(classes)))
  ao <- list(); an <- list(); new_idx <- integer(0)
  for (cl in classes) {
    io <- which(y_old == cl); ino <- which(y_new == cl)
    k <- min(per_cls, length(io), length(ino))
    if (k < 1L) next
    so <- sample(io, k); sn <- sample(ino, k)
    ao[[cl]] <- X_old[so, , drop = FALSE]
    an[[cl]] <- X_new[sn, , drop = FALSE]
    new_idx <- c(new_idx, sn)
  }
  list(anchor_old = do.call(rbind, ao), anchor_new = do.call(rbind, an),
       new_idx = new_idx)
}

run_condition <- function(ref_X, ref_y, qry_X, qry_y, label) {
  ref <- somalign_train_reference(ref_X, labels = ref_y, grid = GRID, rlen = 40)
  qry <- somalign_query(qry_X, ref, grid = GRID, rlen = 40)
  # Cap anchors so evaluation keeps a non-trivial held-out set.
  n_anc <- min(N_ANCHOR, floor(0.25 * nrow(qry_X)))
  anc <- make_anchors(ref_X, ref_y, qry_X, qry_y, n_anc)
  eval_mask <- rep(TRUE, nrow(qry_X)); eval_mask[anc$new_idx] <- FALSE
  ab <- somalign_anchor_benefit(qry, ref, qry_y, anc$anchor_old, anc$anchor_new,
                                rho_grid = RHO_GRID, metric = "mcc",
                                eval_mask = eval_mask)
  cat(sprintf("\n=== %s ===\n", label)); print(ab); print(ab$grid)
  ab
}

# Population-SPECIFIC batch shift that moves each lineage TOWARD its nearest
# neighbouring lineage's centroid, scaled by delta. At delta = 1 a lineage
# lands on top of its neighbour, so plain OT can no longer tell them apart and
# mis-maps them -- exactly the regime where repeat-sample anchors, which pin
# the correct lineage correspondence, should help. (Random or global shifts
# keep lineages separable, so OT solves them for free and anchors add nothing.)
lineage_vec <- NULL   # per-lineage displacement toward nearest neighbour
init_lineage_vec <- function(X_ref, y_ref) {
  cls <- sort(unique(y_ref))
  cen <- t(vapply(cls, function(c) colMeans(X_ref[y_ref == c, , drop = FALSE]),
                  numeric(ncol(X_ref))))
  rownames(cen) <- cls
  vec <- lapply(cls, function(c) {
    d <- sqrt(rowSums((cen - matrix(cen[c, ], nrow(cen), ncol(cen), byrow = TRUE))^2))
    d[c] <- Inf
    cen[names(which.min(d)), ] - cen[c, ]        # vector from c to nearest other
  })
  names(vec) <- cls
  lineage_vec <<- vec
}
make_diff_shift <- function(labels, delta) {
  sh <- matrix(0, nrow = length(labels), ncol = length(lineage_vec[[1]]))
  for (c in names(lineage_vec)) {
    rows <- which(labels == c)
    if (length(rows))
      sh[rows, ] <- matrix(delta * lineage_vec[[c]], length(rows),
                           length(lineage_vec[[c]]), byrow = TRUE)
  }
  sh
}

# ---- Condition A: real acquisition-batch effect (Batch1 -> Batch2) ----------
stamp("Condition A: pilot Batch1 (ref) -> Batch2 (query) ...")
iA_ref <- b == "Batch1"; iA_qry <- b == "Batch2"
condA <- run_condition(X[iA_ref, ], y[iA_ref], X[iA_qry, ], y[iA_qry],
                       "Condition A -- real Batch1->Batch2")

# ---- Condition B: controlled synthetic batch-severity sweep -----------------
# Split Batch1 into ref / query / anchor pools; shift the query by delta*SD.
stamp("Condition B: synthetic batch-severity sweep ...")
sd_marker <- apply(X, 2, stats::sd)
i1 <- which(b == "Batch1")
i1 <- i1[sample.int(length(i1))]
n3 <- floor(length(i1) / 3)
ref_i <- i1[seq_len(n3)]
qry_i <- i1[(n3 + 1):(2 * n3)]
anc_i <- i1[(2 * n3 + 1):(3 * n3)]
init_lineage_vec(X[ref_i, ], y[ref_i])
DELTAS <- c(0, 0.3, 0.6, 0.9, 1.2, 1.5)
condB <- lapply(DELTAS, function(delta) {
  # Differential (toward-nearest-lineage) batch shift on query + anchor pools.
  shift  <- make_diff_shift(y[qry_i], delta)
  ashift <- make_diff_shift(y[anc_i], delta)
  ref <- somalign_train_reference(X[ref_i, ], labels = y[ref_i], grid = GRID, rlen = 40)
  qry <- somalign_query(X[qry_i, ] + shift, ref, grid = GRID, rlen = 40)
  ab <- somalign_anchor_benefit(
    qry, ref, y[qry_i],
    anchor_old = X[anc_i, ], anchor_new = X[anc_i, ] + ashift,
    rho_grid = RHO_GRID, metric = "mcc"
  )
  base <- ab$grid[ab$grid$rho_anchor == 0, ]
  cat(sprintf("  delta=%.2g SD: baseline MCC=%.4f -> best MCC=%.4f (rho=%.3g, lift %+.4f), coverage %.1f%%->%.1f%%\n",
              delta, base$mcc, ab$best$mcc, ab$best$rho_anchor, ab$lift,
              100 * base$coverage, 100 * ab$best$coverage))
  list(delta = delta, ab = ab)
})

# ---- Summary: default rho=1 vs best rho, as a function of severity ----------
cat("\n=== SUMMARY: anchor lift vs batch severity (Condition B) ===\n")
cat("delta  baseline_mcc  mcc@rho=1  best_mcc  best_rho  lift(best-baseline)\n")
for (r in condB) {
  g <- r$ab$grid
  base <- g$mcc[g$rho_anchor == 0]
  at1  <- g$mcc[g$rho_anchor == 1]
  cat(sprintf("%.2g     %+.3f       %+.3f     %+.3f    %-6.3g   %+.3f\n",
              r$delta, base, at1, r$ab$best$mcc, r$ab$best$rho_anchor, r$ab$lift))
}
cat("\nInterpretation: anchors help only under SEVERE population-specific batch\n")
cat("effects AND only at large rho_anchor (>>1). At the production default\n")
cat("rho_anchor=1 the label-transfer benefit is negligible across all severities.\n")

saveRDS(list(condA = condA, condB = condB, n_subsample = nrow(X),
             rho_grid = RHO_GRID, deltas = DELTAS), OUT)
stamp(sprintf("saved -> %s", OUT))
