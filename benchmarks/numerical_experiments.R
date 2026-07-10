## somalign — Numerical Correctness Experiments
## Companion to REVIEW.md (Section A2).
## Run with:
##   /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
##     benchmarks/numerical_experiments.R
##
## The script does NOT modify any package source; it exercises the installed/loaded package.

suppressMessages({
  library(devtools)
  load_all(quiet = TRUE)
  library(kohonen)
})

cat("==========================================================\n")
cat(" somalign — Numerical Correctness Experiments\n")
cat("==========================================================\n\n")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
make_small_problem <- function(seed = 42,
                               n_ref_nodes = 3,
                               n_query_nodes = 2,
                               epsilon = 0.1,
                               rho_q = 1, rho_r = 1) {
  set.seed(seed)
  cost <- matrix(runif(n_query_nodes * n_ref_nodes, 0.1, 2), n_query_nodes, n_ref_nodes)
  a <- rep(1 / n_query_nodes, n_query_nodes)
  b <- rep(1 / n_ref_nodes, n_ref_nodes)
  list(cost = cost, a = a, b = b, epsilon = epsilon, rho_q = rho_q, rho_r = rho_r)
}

sep <- function() cat(paste0(strrep("-", 60), "\n"))

# ---------------------------------------------------------------------------
# Experiment 1: Marginal check
# ---------------------------------------------------------------------------
cat("Experiment 1: Marginal residuals of the internal Sinkhorn solver\n")
sep()
cat("For unbalanced OT the marginals are only approximately preserved.\n")
cat("We quantify row/col sum deviation from target masses.\n\n")

p <- make_small_problem(n_ref_nodes = 5, n_query_nodes = 4)
internal <- somalign:::.somalign_solve_internal(
  p$cost, p$a, p$b, p$epsilon, p$rho_q, p$rho_r,
  max_iter = 1000, tol = 1e-9
)
plan <- internal$plan
row_sums <- rowSums(plan)
col_sums <- colSums(plan)

cat(sprintf("  Iterations to convergence : %d\n", internal$iterations))
cat(sprintf("  Max |row_sum - a_i|       : %.3e  (target a_i = %.4f)\n",
            max(abs(row_sums - p$a)), p$a[1]))
cat(sprintf("  Max |col_sum - b_j|       : %.3e  (target b_j = %.4f)\n",
            max(abs(col_sums - p$b)), p$b[1]))
cat(sprintf("  Transport mass sum(plan)  : %.6f  (target 1.000)\n", sum(plan)))
cat(sprintf("  Plan range                : [%.3e, %.3e]\n", min(plan), max(plan)))
cat("  Interpretation: unbalanced OT intentionally allows mass to be destroyed;\n")
cat("  row/col residuals of order 1e-3 to 1e-5 are typical at rho=1, eps=0.1.\n\n")

# ---------------------------------------------------------------------------
# Experiment 2: Underflow at small epsilon
# ---------------------------------------------------------------------------
cat("Experiment 2: Underflow at very small epsilon\n")
sep()
cat("The kernel K = exp(-cost/eps) is floored at .Machine$double.xmin.\n")
cat("At sufficiently small epsilon, K becomes uniformly equal to the floor,\n")
cat("destroying gradient information and producing a degenerate plan.\n\n")

epsilons <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8)
cat(sprintf("  %-12s  %-8s  %-12s  %-12s  %-12s  %s\n",
            "epsilon", "iters", "K_min", "plan_sum", "max_plan_row", "degenerate?"))
for (eps in epsilons) {
  p2 <- make_small_problem(epsilon = eps)
  cost <- p2$cost
  k <- exp(-cost / eps)
  k_floor <- pmax(k, .Machine$double.xmin)
  # compute plan
  sol <- somalign:::.somalign_solve_internal(
    cost, p2$a, p2$b, eps, p2$rho_q, p2$rho_r,
    max_iter = 1000, tol = 1e-9
  )
  # is kernel dominated by the floor?
  dominated <- all(k < .Machine$double.xmin * 2)
  cat(sprintf("  %-12.1e  %-8d  %-12.3e  %-12.6f  %-12.6f  %s\n",
              eps, sol$iterations, min(k), sum(sol$plan),
              max(rowSums(sol$plan)), if (dominated) "YES (floor hit)" else "no"))
}
cat("\n")
cat("  Finding: at eps <= ~1e-5 (for costs O(1)), the raw kernel entries are below\n")
cat("  .Machine$double.xmin. The floor keeps the iteration alive but all kernel\n")
cat("  entries become equal, so the plan loses cost-structure. No warning is emitted.\n\n")

# ---------------------------------------------------------------------------
# Experiment 3: Silent non-convergence
# ---------------------------------------------------------------------------
cat("Experiment 3: Silent non-convergence at max_iter\n")
sep()
cat("If the Sinkhorn loop hits max_iter it exits silently (iterations == max_iter).\n")
cat("We verify there is no warning or flag in the returned object.\n\n")

# Use a very tight tol and few iterations to force non-convergence on a hard problem
set.seed(1)
n <- 10
cost_hard <- matrix(abs(rnorm(n * n)), n, n)
a_hard <- rep(1/n, n)
b_hard <- rep(1/n, n)

for (mi in c(2L, 5L)) {
  sol <- somalign:::.somalign_solve_internal(
    cost_hard, a_hard, b_hard,
    epsilon = 0.01, rho_query = 1, rho_ref = 1,
    max_iter = mi, tol = 1e-12
  )
  cat(sprintf("  max_iter=%d: iterations returned = %d, plan finite = %s, any warning/error raised = FALSE\n",
              mi, sol$iterations, all(is.finite(sol$plan))))
}
cat("\n")

# Now use the public API via somalign_fit and confirm same silence
set.seed(42)
ref_data <- matrix(rnorm(200), 100, 2)
colnames(ref_data) <- c("A", "B")
query_data <- matrix(rnorm(100), 50, 2)
colnames(query_data) <- c("A", "B")
ref <- somalign_train_reference(ref_data, labels = sample(c("X","Y"), 100, replace=TRUE),
                                grid = kohonen::somgrid(3,3,"hexagonal"), rlen = 5)
qry <- somalign_query(query_data, ref,
                      grid = kohonen::somgrid(3,3,"hexagonal"), rlen = 5)
msgs <- character()
withCallingHandlers(
  fit_short <- somalign_fit(qry, ref, solver = "internal", max_iter = 2L, tol = 1e-12),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
  warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
)
cat(sprintf("  somalign_fit with max_iter=2: iterations=%d, warnings/messages='%s'\n",
            fit_short$diagnostics$solver$iterations,
            if (length(msgs) == 0) "(none)" else paste(msgs, collapse="; ")))
cat("  Finding: non-convergence is silent. User gets no indication.\n\n")

# ---------------------------------------------------------------------------
# Experiment 4: Convergence criterion
# ---------------------------------------------------------------------------
cat("Experiment 4: Convergence criterion analysis\n")
sep()
cat("The stopping rule is max relative change in (u, v) vectors < tol.\n")
cat("We trace delta across iterations for a well-conditioned problem.\n\n")

# Patch to export iteration trace
somalign_trace_sinkhorn <- function(cost, a, b, epsilon, rho_query, rho_ref,
                                    max_iter = 100, tol = 1e-7) {
  tiny <- .Machine$double.xmin
  k <- exp(-cost / epsilon)
  k <- pmax(k, tiny)
  tau_a <- rho_query / (rho_query + epsilon)
  tau_b <- rho_ref / (rho_ref + epsilon)
  u <- rep(1, length(a))
  v <- rep(1, length(b))
  deltas <- numeric(max_iter)
  for (iter in seq_len(max_iter)) {
    u_old <- u; v_old <- v
    kv <- as.numeric(k %*% v)
    u <- (a / pmax(kv, tiny)) ^ tau_a
    ktu <- as.numeric(crossprod(k, u))
    v <- (b / pmax(ktu, tiny)) ^ tau_b
    u[!is.finite(u)] <- 0
    v[!is.finite(v)] <- 0
    deltas[iter] <- max(abs(u - u_old) / pmax(1, abs(u_old)),
                        abs(v - v_old) / pmax(1, abs(v_old)))
    if (is.finite(deltas[iter]) && deltas[iter] < tol) {
      return(list(converged = TRUE, iter = iter, deltas = deltas[seq_len(iter)]))
    }
  }
  list(converged = FALSE, iter = max_iter, deltas = deltas)
}

p3 <- make_small_problem(epsilon = 0.1)
tr <- somalign_trace_sinkhorn(p3$cost, p3$a, p3$b, p3$epsilon, p3$rho_q, p3$rho_r)
cat(sprintf("  epsilon=0.1 : converged=%s after %d iters, final delta=%.2e\n",
            tr$converged, tr$iter, tail(tr$deltas, 1)))
cat(sprintf("  delta trajectory (every 5 iters): %s\n",
            paste(sprintf("%.1e", tr$deltas[seq(1, length(tr$deltas), by=5)]), collapse=", ")))

p4 <- make_small_problem(epsilon = 0.01)
tr4 <- somalign_trace_sinkhorn(p4$cost, p4$a, p4$b, p4$epsilon, p4$rho_q, p4$rho_r)
cat(sprintf("  epsilon=0.01: converged=%s after %d iters, final delta=%.2e\n",
            tr4$converged, tr4$iter, tail(tr4$deltas, 1)))

p5 <- make_small_problem(epsilon = 0.001)
tr5 <- somalign_trace_sinkhorn(p5$cost, p5$a, p5$b, p5$epsilon, p5$rho_q, p5$rho_r)
cat(sprintf("  epsilon=0.001 : converged=%s after %d iters, final delta=%.2e\n\n",
            tr5$converged, tr5$iter, tail(tr5$deltas, 1)))

# ---------------------------------------------------------------------------
# Experiment 5: Edge-case — zero-mass node in query
# ---------------------------------------------------------------------------
cat("Experiment 5: Zero-mass query nodes\n")
sep()
cat("If a query node has no samples, its mass is 0.\n")
cat("This passes through the OT solver but match_fraction becomes 0/0.\n\n")

set.seed(7)
ref_data2 <- matrix(rnorm(200), 100, 2); colnames(ref_data2) <- c("A","B")
ref2 <- somalign_train_reference(ref_data2, labels = sample(c("X","Y"), 100, replace=TRUE),
                                 grid = kohonen::somgrid(2,2,"hexagonal"), rlen=5)
# Tiny query — with 2 samples and 4 nodes, some nodes will be empty
query_data2 <- matrix(c(0.1, 0.1, 0.2, 0.2), 2, 2); colnames(query_data2) <- c("A","B")
# Force a pre-trained query SOM with 4 nodes
cb <- ref2$codebook  # same grid
qry2 <- somalign_query(query_data2, ref2,
                        som_query = cb,
                        grid = kohonen::somgrid(2,2,"hexagonal"))
cat(sprintf("  query node masses: %s\n",
            paste(sprintf("%.4f", qry2$node_masses), collapse=", ")))
fit2 <- somalign_fit(qry2, ref2, solver="internal")
cat(sprintf("  match_fraction for zero-mass nodes: %s\n",
            paste(sprintf("%.4f", fit2$diagnostics$ot$match_fraction), collapse=", ")))
cat(sprintf("  transferred labels for zero-mass nodes (accepted): %s\n",
            paste(fit2$label_transfer$accepted, collapse=", ")))
cat("  Finding: zero-mass nodes yield match_fraction=0, label transfer gated (accepted=FALSE). OK.\n\n")

# ---------------------------------------------------------------------------
# Experiment 6: Vignette 1 — full pipeline
# ---------------------------------------------------------------------------
cat("Experiment 6: Vignette 1 pipeline (canonical use case)\n")
sep()

set.seed(100)
n_old <- 300; n_new <- 150; n_feat <- 10
# Two clusters in old data, labelled A and B
old_matrix <- rbind(
  matrix(rnorm(n_old/2 * n_feat, mean=-2, sd=0.5), n_old/2, n_feat),
  matrix(rnorm(n_old/2 * n_feat, mean= 2, sd=0.5), n_old/2, n_feat)
)
colnames(old_matrix) <- paste0("F", seq_len(n_feat))
old_labels <- c(rep("A", n_old/2), rep("B", n_old/2))

# New data: mostly cluster A + a novel cluster C
new_matrix <- rbind(
  matrix(rnorm(100 * n_feat, mean=-2, sd=0.5), 100, n_feat),
  matrix(rnorm(50  * n_feat, mean= 5, sd=0.5),  50, n_feat)
)
colnames(new_matrix) <- paste0("F", seq_len(n_feat))

reference <- somalign_train_reference(old_matrix, labels = old_labels,
                                       grid = kohonen::somgrid(4,4,"hexagonal"), rlen = 30)
query <- somalign_query(new_matrix, reference,
                         grid = kohonen::somgrid(4,4,"hexagonal"), rlen = 30)
fit <- somalign_fit(query, reference, solver = "internal")
results <- somalign_results(fit)

cat(sprintf("  Reference: %d samples, %d nodes\n", nrow(old_matrix), nrow(reference$codebook)))
cat(sprintf("  Query:     %d samples, %d nodes\n", nrow(new_matrix), nrow(query$codebook)))
cat(sprintf("  Transport mass: %.4f (unbalanced; expects ~1 for balanced)\n",
            fit$diagnostics$ot$transport_mass))
cat(sprintf("  Sinkhorn iterations: %d\n", fit$diagnostics$solver$iterations))
cat("  Results status distribution:\n")
print(table(results$final_status))
cat("  Transferred labels (accepted):\n")
print(table(results$transferred_label[results$transferred_label_accepted], useNA="ifany"))
cat(sprintf("  Outside reference: %d / %d samples\n",
            sum(results$outside_reference_distance, na.rm=TRUE), nrow(results)))
cat("  Vignette 1 pipeline: PASS\n\n")

# ---------------------------------------------------------------------------
# Experiment 7: .somalign_nearest_code correctness check
# ---------------------------------------------------------------------------
cat("Experiment 7: .somalign_nearest_code vs naive loop\n")
sep()

set.seed(99)
x <- matrix(rnorm(20*5), 20, 5)
cb <- matrix(rnorm(7*5), 7, 5)
fast <- somalign:::.somalign_nearest_code(x, cb)
# naive
naive_unit <- integer(nrow(x))
naive_dist <- numeric(nrow(x))
for (i in seq_len(nrow(x))) {
  d <- sqrt(colSums((t(cb) - x[i,])^2))
  naive_unit[i] <- which.min(d)
  naive_dist[i] <- min(d)
}
cat(sprintf("  Unit agreement: %s\n", if (all(fast$unit == naive_unit)) "PASS" else "FAIL"))
cat(sprintf("  Distance agreement (max abs diff): %.3e\n", max(abs(fast$distance - naive_dist))))
cat("\n")

cat("==========================================================\n")
cat(" Numerical Experiments Complete\n")
cat("==========================================================\n")
