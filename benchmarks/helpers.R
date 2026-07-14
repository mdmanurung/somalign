## Shared helpers for the somalign benchmark scripts.
## Sourced from repo root, e.g. source("benchmarks/helpers.R").

sep <- function(title = "") {
  cat(sprintf("\n── %s %s\n", title, strrep("─", max(0L, 56L - nchar(title)))))
}

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

js_div <- function(p, q) {
  m  <- (p + q) / 2
  kl <- function(a, b) sum(ifelse(a > 0, a * log(a / b), 0))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

node_dist <- function(units, n_nodes) {
  v <- tabulate(units, nbins = n_nodes) / length(units)
  v + 1e-10   # smooth zeros for KL
}

ks_per_marker <- function(mat1, mat2) {
  vapply(seq_len(ncol(mat1)), function(j)
    ks.test(mat1[, j], mat2[, j])$statistic, numeric(1))
}
