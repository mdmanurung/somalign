# Plan 01 — Epsilon Phase-Transition Diagnostic and Critical-Epsilon Estimator

## 1. Summary

Add `somalign_epsilon_sweep()`, a cheap diagnostic that runs the Sinkhorn OT
solve across a log-spaced epsilon grid **without the per-cell projection step**,
computes the transport-plan row-entropy (order parameter Phi) and its
susceptibility d(Phi)/d(log epsilon) over the grid, and returns the susceptibility
peak as an estimated critical epsilon `epsilon_c`. Simultaneously, expose the dual
free energy (`log_Z`) as a new scalar field `diagnostics$solver$log_Z` in every
standard `somalign_fit()` call (requires a one-line addition to the log-domain
solver and a one-line addition to `.somalign_build_diagnostics`).

---

## 2. Public API

### 2a. New exported function

```r
somalign_epsilon_sweep(
  query,                      # somalign_query object
  reference,                  # somalign_reference object
  epsilon_grid = NULL,        # numeric vector; if NULL uses .somalign_default_eps_grid()
  n_grid       = 25L,         # used only when epsilon_grid = NULL
  rho_query    = 1,
  rho_ref      = 1,
  solver       = c("log_domain", "internal", "auto"),
  max_iter     = 1000,
  tol          = 1e-7,
  diagonal_boost = 0,
  label_guided   = FALSE,
  parallel       = FALSE
)
```

**Returns** a list of class `"somalign_epsilon_sweep"` with:

```
$table        # data.frame, one row per epsilon value:
              #   epsilon, log_epsilon, Phi, log_Z, iterations, converged
$epsilon_c    # scalar — epsilon at max susceptibility; NA if sweep has < 3 points
$epsilon_rec  # scalar — 0.3 * epsilon_c (recommended working epsilon)
$susceptibility  # numeric vector — d(Phi)/d(log_epsilon), NA at boundaries
$cost_scale   # scalar — median positive entry of the raw cost matrix (for reference)
```

`epsilon_grid` defaults to `NULL`, in which case `.somalign_default_eps_grid()`
generates a log-spaced vector from `1e-3` to `5` (25 points). The default solver
is `"log_domain"` because free-energy extraction requires the log potentials;
`"internal"` is accepted but skips the `log_Z` column (filled with `NA`).

No existing function signatures are changed. Default behavior of `somalign_fit()`
is preserved: `log_Z` appears in `diagnostics$solver` only when
`solver = "log_domain"`, otherwise `NA`.

### 2b. New S3 methods (minimal)

```r
print.somalign_epsilon_sweep(x, ...)   # prints table + epsilon_c recommendation
plot.somalign_epsilon_sweep(x, ...)    # delegates to .somalign_plot_eps_sweep()
```

---

## 3. Internal helpers

All new helpers go in the files indicated below.

| Helper | File | Purpose |
|--------|------|---------|
| `.somalign_default_eps_grid(n)` | `R/diagnostics.R` | `exp(seq(log(1e-3), log(5), length.out = n))` |
| `.somalign_plan_entropy(plan)` | `R/ot.R` | M-vector of per-row entropies; calls `.somalign_row_normalize()` then `.somalign_entropy()` |
| `.somalign_order_parameter(plan, K)` | `R/ot.R` | scalar Phi = mean(exp(H_i)) / K |
| `.somalign_susceptibility(Phi_vec, log_eps_vec)` | `R/diagnostics.R` | finite-difference d(Phi)/d(log_eps); length-preserving with NA padding at boundaries |
| `.somalign_log_Z_from_potentials(f, g, log_a, log_b, epsilon)` | `R/ot.R` | scalar dual free energy — formula in §5 |
| `.somalign_sweep_one_eps(cost_norm, a, b, eps, rho_q, rho_r, solver, max_iter, tol)` | `R/diagnostics.R` | runs `.somalign_solve_ot()` for a single epsilon; extracts Phi and log_Z; returns a one-row data.frame |
| `.somalign_plot_eps_sweep(sweep, ...)` | `R/plot.R` | ggplot2 two-panel: Phi vs log(eps) + susceptibility vs log(eps); vertical line at `epsilon_c` |

Signature details:

```r
# R/ot.R
.somalign_plan_entropy <- function(plan)
  # plan: M x K matrix
  # returns: numeric M-vector; NA for all-zero rows

.somalign_order_parameter <- function(plan, K)
  # K = ncol(plan) — the number of reference nodes
  # returns: scalar in (1/K, 1]

.somalign_log_Z_from_potentials <- function(f, g, log_a, log_b, epsilon,
                                             rho_query, rho_ref)
  # f: M-vector of converged log-potentials (as returned by log-domain solver)
  # g: K-vector of converged log-potentials
  # returns: scalar — dual free energy; NA if any input non-finite

# R/diagnostics.R
.somalign_default_eps_grid <- function(n = 25L)

.somalign_susceptibility <- function(Phi_vec, log_eps_vec)
  # finite differences on (Phi, log_eps); NA at first and last positions

.somalign_sweep_one_eps <- function(cost_norm, a, b, eps,
                                    rho_query, rho_ref, solver,
                                    max_iter, tol)
  # returns data.frame(epsilon, log_epsilon, Phi, log_Z, iterations, converged)
```

---

## 4. Data-structure changes

### 4a. `diagnostics$solver$log_Z` (every `somalign_fit` call)

`log_Z` is added to the `solver` sub-list built in `.somalign_build_diagnostics()`
(`R/fit.R`, lines 273–289). It holds:

- The scalar dual free energy when `solver = "log_domain"` (computed in
  `.somalign_solve_internal_log()` from converged potentials `f`, `g`).
- `NA_real_` for `solver = "internal"` or `"auto"`.

**Location in the return value of `somalign_diagnostics()`:**
```
fit$diagnostics$solver$log_Z   # scalar or NA_real_
```

**Computed from:** converged `f` (M-vector) and `g` (K-vector) returned by
`.somalign_solve_internal_log()`, together with `log_a = log(a)`,
`log_b = log(b)`, and `epsilon`, `rho_query`, `rho_ref`.

### 4b. `somalign_epsilon_sweep` return value

New S3 object — not embedded in `somalign_fit`; created and returned only
by `somalign_epsilon_sweep()`. No changes to existing `somalign_fit` S3
structure beyond the `$log_Z` scalar described above.

---

## 5. Algorithm

### 5a. Order parameter and susceptibility

For the converged transport plan `P` (M x K matrix, unnormalized):

1. **Row-normalize** to get the correspondence matrix `p`:
   ```r
   p <- .somalign_row_normalize(plan)  # already exists in R/ot.R, line 247
   ```

2. **Per-row entropy** (nats):
   ```r
   H_i <- .somalign_plan_entropy(plan)
   # H_i[i] = -sum_j p[i,j] * log(p[i,j])   (0*log0 = 0 by convention)
   # H_i[i] = NA   when row_i is all zero
   ```
   This calls the already-present `.somalign_entropy()` (`R/utils.R`, line 526)
   row-by-row.

3. **Effective reference nodes** per query node:
   ```r
   S_i <- exp(H_i)   # ranges from 1 (one node used) to K (all nodes equally)
   ```

4. **Order parameter** (fractional usage of reference codebook):
   ```r
   K <- ncol(plan)
   Phi <- mean(S_i[is.finite(S_i)]) / K   # in (1/K, 1]
   ```
   At epsilon → 0: Phi → 1/K. At epsilon → infinity: Phi → 1.

5. **Susceptibility** (finite difference over log-epsilon grid):
   ```r
   log_eps <- log(epsilon_grid)
   chi <- .somalign_susceptibility(Phi_vec, log_eps)
   # chi[k] = (Phi[k+1] - Phi[k-1]) / (log_eps[k+1] - log_eps[k-1])
   # chi[1] = chi[n] = NA
   ```

6. **Critical epsilon**:
   ```r
   k_star <- which.max(chi)
   epsilon_c <- epsilon_grid[k_star]
   epsilon_rec <- 0.3 * epsilon_c
   ```

### 5b. Dual free energy (log Z)

The log-domain Sinkhorn converges to dual potentials `f` (M-vector) and `g`
(K-vector). The dual objective of the UOT problem is:

```
log_Z = sum_i [tau_a * epsilon * log_a[i] + (1 - tau_a) * f[i]]
      + sum_j [tau_b * epsilon * log_b[j] + (1 - tau_b) * g[j]]
```

where `tau_a = rho_query / (rho_query + epsilon)` and
`tau_b = rho_ref / (rho_ref + epsilon)`.

Alternatively (simpler, equivalent at convergence) from the already-computed
`lse_g` vector inside `.somalign_solve_internal_log()` (line 123):

```
log_Z = epsilon * sum_i lse_g[i]
```

where `lse_g[i] = logsumexp_j((g_j - C_ij) / eps)` is already on the stack at
the end of each outer iteration. After convergence, this is the log-normalizer
of the Gibbs kernel summed over query nodes.

**Implementation: one-line addition to `.somalign_solve_internal_log()`**
(`R/ot.R`, after line 151):

```r
log_Z <- epsilon * sum(lse_g[is.finite(lse_g)])
```

Return it alongside `plan`:
```r
list(plan = plan, iterations = iterations, converged = converged,
     final_delta = final_delta, log_Z = log_Z)
```

Then propagate it up the call stack:
- `.somalign_solve_internal()` (line 83): when `log_domain = TRUE`, pass
  `log_Z` through; when `log_domain = FALSE`, add `log_Z = NA_real_`.
- `.somalign_solve_ot()` (line 24): include `log_Z` in the return list.
- `.somalign_build_diagnostics()` (`R/fit.R`, line 284): add
  `log_Z = transport$ot$log_Z` to the `solver` sub-list.

### 5c. Core sweep loop (R pseudocode)

```r
somalign_epsilon_sweep <- function(query, reference,
                                   epsilon_grid = NULL, n_grid = 25L,
                                   rho_query = 1, rho_ref = 1,
                                   solver = c("log_domain", "internal", "auto"),
                                   max_iter = 1000, tol = 1e-7,
                                   diagonal_boost = 0,
                                   label_guided = FALSE, parallel = FALSE) {
  solver <- match.arg(solver)
  ## --- validation (delegate to existing helpers) --------------------------
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  .somalign_check_pos_scalar(rho_query, "rho_query")
  .somalign_check_pos_scalar(rho_ref,   "rho_ref")

  ## --- prepare cost once ---------------------------------------------------
  label_mask <- .somalign_epsilon_sweep_label_mask(query, reference, label_guided)
  cost_raw   <- .somalign_pairwise_distance(query$codebook, reference$codebook)
  prepared   <- .somalign_prepare_cost(cost_raw, diagonal_boost,
                                       cost_bonus = NULL, label_mask = label_mask)
  cost_norm  <- prepared$cost_normalized
  cost_scale <- prepared$cost_scale

  ## --- epsilon grid --------------------------------------------------------
  if (is.null(epsilon_grid))
    epsilon_grid <- .somalign_default_eps_grid(n_grid)
  epsilon_grid <- sort(unique(.somalign_validate_grid_vector(epsilon_grid, "epsilon_grid")))

  ## --- sweep ---------------------------------------------------------------
  .run_one <- function(k)
    .somalign_sweep_one_eps(cost_norm, query$node_masses, reference$node_masses,
                            epsilon_grid[k], rho_query, rho_ref,
                            solver, max_iter, tol)

  rows  <- .somalign_run_grid(length(epsilon_grid), .run_one, parallel)
  table <- do.call(rbind, rows)

  ## --- susceptibility and epsilon_c ----------------------------------------
  table$susceptibility <- .somalign_susceptibility(table$Phi, log(table$epsilon))
  epsilon_c  <- .somalign_locate_epsilon_c(table)
  epsilon_rec <- if (is.finite(epsilon_c)) 0.3 * epsilon_c else NA_real_

  structure(
    list(table        = table,
         epsilon_c    = epsilon_c,
         epsilon_rec  = epsilon_rec,
         cost_scale   = cost_scale),
    class = "somalign_epsilon_sweep"
  )
}
```

`.somalign_sweep_one_eps()` body:

```r
.somalign_sweep_one_eps <- function(cost_norm, a, b, eps,
                                    rho_query, rho_ref, solver, max_iter, tol) {
  ot   <- .somalign_solve_ot(cost_norm, a, b, eps, rho_query, rho_ref,
                             solver, max_iter, tol)
  plan <- ot$plan
  Phi  <- .somalign_order_parameter(plan, ncol(plan))
  log_Z <- if (!is.null(ot$log_Z)) ot$log_Z else NA_real_
  data.frame(epsilon     = eps,
             log_epsilon = log(eps),
             Phi         = Phi,
             log_Z       = log_Z,
             iterations  = ot$iterations,
             converged   = ot$converged)
}
```

`.somalign_locate_epsilon_c()`:

```r
.somalign_locate_epsilon_c <- function(table) {
  chi <- table$susceptibility
  finite_idx <- which(is.finite(chi))
  if (length(finite_idx) == 0L) return(NA_real_)
  k_star <- finite_idx[which.max(chi[finite_idx])]
  table$epsilon[k_star]
}
```

---

## 6. Integration points

### Which existing functions call or expose this

| Existing function | Change |
|---|---|
| `.somalign_solve_internal_log()` (`R/ot.R`, line 96) | Add `log_Z` computation and return field |
| `.somalign_solve_internal()` (`R/ot.R`, line 35) | Pass `log_Z` through when `log_domain = TRUE`; set `NA_real_` otherwise |
| `.somalign_solve_ot()` (`R/ot.R`, line 1) | Include `log_Z` in returned list |
| `.somalign_build_diagnostics()` (`R/fit.R`, line 264) | Add `log_Z = transport$ot$log_Z` to `solver` sub-list |

### How default behavior is preserved

- `somalign_fit()` with default `solver = "internal"` will have
  `diagnostics$solver$log_Z = NA_real_`. No warnings are emitted.
- `somalign_epsilon_sweep()` does not call `.somalign_finish_fit()` or
  `.somalign_project_pair()`, so no per-cell projection occurs. It only calls
  `.somalign_prepare_cost()` (once), `.somalign_solve_ot()` (once per epsilon),
  and the new order-parameter helper.
- No existing function signatures are modified; all changes are additive.

---

## 7. Edge cases and numerical stability

| Scenario | Handling |
|---|---|
| All-zero query row in plan | `.somalign_plan_entropy()` returns `NA` for that row; `Phi` averages over finite rows only; if all rows are NA, `Phi = NA` and a warning is emitted |
| `log(0)` in entropy computation | `.somalign_entropy()` already filters `prob > 0` (line 527 of `R/utils.R`); no change needed |
| Single-epsilon input (`length(epsilon_grid) == 1`) | `susceptibility` is a single `NA`; `epsilon_c = NA`; a message advises using at least 3 points |
| Non-convergence at a grid point | Row records `converged = FALSE`, `Phi` is computed from the non-converged plan (may be unreliable); susceptibility peak search excludes NA susceptibility values but does not exclude non-converged rows — caller must inspect `table$converged` |
| `log_Z` overflow/underflow | `lse_g` uses `.somalign_logsumexp()` which is numerically stable; the final sum is over M finite values; if any `lse_g[i]` is non-finite it is excluded via `is.finite()` guard before summing |
| Small codebook (K = 4) | Phi is bounded in (0.25, 1.0); the transition may be a gentle monotone curve with no clear peak; `epsilon_c` is still returned as the argmax of susceptibility, but may simply be the first finite susceptibility value; the caller can inspect the plot to judge whether the peak is meaningful |
| Duplicate epsilon values in `epsilon_grid` | Deduplicated and sorted at the top of `somalign_epsilon_sweep()` with `sort(unique(...))` |

---

## 8. Tests

File: `tests/testthat/test-epsilon-sweep.R`

```r
# Setup shared across tests
local_setup <- function() {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 100, ncol = 4,
                dimnames = list(NULL, c("F1","F2","F3","F4")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(4, 4, "hexagonal"),
                                  rlen = 10)
  qry <- somalign_query(mat + 0.5, ref, grid = kohonen::somgrid(4, 4, "hexagonal"),
                        rlen = 10)
  list(ref = ref, qry = qry)
}
```

| Test name | What it asserts |
|---|---|
| `test_that("sweep returns somalign_epsilon_sweep class", ...)` | `inherits(sw, "somalign_epsilon_sweep")` |
| `test_that("table has correct columns", ...)` | `names(sw$table)` contains `epsilon`, `log_epsilon`, `Phi`, `log_Z`, `iterations`, `converged`, `susceptibility` |
| `test_that("Phi is monotone increasing in epsilon", ...)` | `all(diff(sw$table$Phi[is.finite(sw$table$Phi)]) >= -1e-6)` (Phi is non-decreasing as epsilon increases) |
| `test_that("Phi is in (1/K, 1]", ...)` | `all(sw$table$Phi[is.finite(sw$table$Phi)] > 0 & sw$table$Phi[is.finite(sw$table$Phi)] <= 1)` |
| `test_that("epsilon_c is within epsilon_grid range", ...)` | `sw$epsilon_c >= min(eps_grid) & sw$epsilon_c <= max(eps_grid)` |
| `test_that("epsilon_rec is 0.3 * epsilon_c", ...)` | `all.equal(sw$epsilon_rec, 0.3 * sw$epsilon_c)` |
| `test_that("log_Z is NA for internal solver", ...)` | `sw_int <- somalign_epsilon_sweep(qry, ref, solver = "internal"); all(is.na(sw_int$table$log_Z))` |
| `test_that("log_Z is finite for log_domain solver", ...)` | `sw_log <- somalign_epsilon_sweep(qry, ref, solver = "log_domain"); all(is.finite(sw_log$table$log_Z))` |
| `test_that("log_Z is monotone decreasing in epsilon", ...)` | `all(diff(sw_log$table$log_Z) <= 0)` — free energy decreases as regularization increases (plan spreads, energy rises but log_Z as defined here decreases) |
| `test_that("single-epsilon sweep returns NA epsilon_c with message", ...)` | `expect_message(sw1 <- somalign_epsilon_sweep(qry, ref, epsilon_grid = 0.1)); is.na(sw1$epsilon_c)` |
| `test_that("fit diagnostics$solver$log_Z is NA for internal solver", ...)` | `fit <- somalign_fit(qry, ref, solver = "internal"); is.na(fit$diagnostics$solver$log_Z)` |
| `test_that("fit diagnostics$solver$log_Z is finite for log_domain solver", ...)` | `fit <- somalign_fit(qry, ref, solver = "log_domain"); is.finite(fit$diagnostics$solver$log_Z)` |
| `test_that("all-zero mass query node does not break sweep", ...)` | Create a query with one zero-mass node; sweep completes without error; that row's contribution to Phi is excluded |
| `test_that("print.somalign_epsilon_sweep does not error", ...)` | `expect_output(print(sw))` |
| `test_that("plot.somalign_epsilon_sweep returns ggplot", ...)` | `p <- plot(sw); inherits(p, "ggplot")` |

---

## 9. Docs / NAMESPACE

### Roxygen skeleton for `somalign_epsilon_sweep()`

```r
#' Epsilon phase-transition sweep for principled epsilon selection
#'
#' Runs the Sinkhorn OT solve across a log-spaced epsilon grid without the
#' per-cell projection step, computes the transport-plan order parameter
#' (mean fractional reference-node usage) and its susceptibility
#' (d(Phi)/d(log epsilon)), and returns the susceptibility peak as a
#' principled critical epsilon.
#'
#' @param query A \code{somalign_query} object.
#' @param reference A \code{somalign_reference} object.
#' @param epsilon_grid Numeric vector of epsilon values. If \code{NULL}
#'   (default), a log-spaced grid of \code{n_grid} values from 1e-3 to 5 is
#'   used.
#' @param n_grid Integer. Grid size when \code{epsilon_grid = NULL}. Default 25.
#' @param rho_query,rho_ref Mass relaxation parameters passed to the OT solver.
#' @param solver Sinkhorn solver. Default \code{"log_domain"} (required for
#'   \code{log_Z}; \code{"internal"} fills \code{log_Z} with \code{NA}).
#' @param max_iter,tol Sinkhorn convergence parameters.
#' @param diagonal_boost Non-negative cost reduction on nearest-reference-node
#'   entries. Default 0.
#' @param label_guided Logical; see \code{\link{somalign_fit}}.
#' @param parallel Logical; see \code{\link{somalign_sensitivity_grid}}.
#'
#' @return A list of class \code{"somalign_epsilon_sweep"} with components
#'   \code{table}, \code{epsilon_c}, \code{epsilon_rec},
#'   \code{cost_scale}. The \code{table} data frame contains columns
#'   \code{epsilon}, \code{log_epsilon}, \code{Phi}, \code{log_Z},
#'   \code{iterations}, \code{converged}, \code{susceptibility}.
#'
#' @details
#' The order parameter \eqn{\Phi(\epsilon)} is the mean effective fraction of
#' reference nodes used per query node: \eqn{\Phi = \mathrm{mean}_i(\exp H_i) / K},
#' where \eqn{H_i = -\sum_j p_{ij} \log p_{ij}} is the Shannon entropy of the
#' row-normalised transport plan and \eqn{K} is the number of reference nodes.
#' As \eqn{\epsilon \to 0}, \eqn{\Phi \to 1/K}; as \eqn{\epsilon \to \infty},
#' \eqn{\Phi \to 1}. The susceptibility \eqn{\chi = d\Phi/d\log\epsilon}
#' peaks at the critical epsilon \eqn{\epsilon_c}, which marks the crossover
#' between a localised (transport-cost-dominated) and a delocalised
#' (entropy-dominated) plan. Working at \eqn{0.3\,\epsilon_c} keeps the plan
#' in the ordered phase with a safety margin against numerical instability.
#'
#' The sweep avoids the per-cell projection step
#' (\code{.somalign_project_pair}), so it runs approximately one OT solve per
#' epsilon value — typically 10–50x faster than
#' \code{\link{somalign_sensitivity_grid}} for the same epsilon range.
#'
#' @seealso \code{\link{somalign_fit}}, \code{\link{somalign_sensitivity_grid}},
#'   \code{\link{somalign_diagnostics}}
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
#'               dimnames = list(NULL, c("F1", "F2")))
#' ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
#'                                 rlen = 5)
#' qry <- somalign_query(mat + 0.5, ref,
#'                       grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
#' sw <- somalign_epsilon_sweep(qry, ref, n_grid = 10)
#' sw$epsilon_c
#' sw$epsilon_rec
#' plot(sw)
#'
#' @export
```

**`devtools::document()`** regenerates `NAMESPACE` and `man/somalign_epsilon_sweep.Rd`.
Additionally add `@export` to `print.somalign_epsilon_sweep` and
`plot.somalign_epsilon_sweep`.

### NAMESPACE entries generated by roxygen2

```
export(somalign_epsilon_sweep)
S3method(print, somalign_epsilon_sweep)
S3method(plot, somalign_epsilon_sweep)
```

---

## 10. BiocCheck constraints

`somalign_epsilon_sweep()` must stay at or below 50 lines. The body as sketched in §5c
has ~30 lines. All substantive computation is delegated to:
- `.somalign_default_eps_grid()` (2 lines)
- `.somalign_epsilon_sweep_label_mask()` (3 lines — thin wrapper for the
  existing label-mask logic, identical to the block in `somalign_fit()` lines
  99–107 of `R/fit.R`)
- `.somalign_sweep_one_eps()` (~10 lines)
- `.somalign_susceptibility()` (~6 lines)
- `.somalign_locate_epsilon_c()` (~5 lines)
- `.somalign_run_grid()` (already exists in `R/diagnostics.R`, line 140)

`print.somalign_epsilon_sweep` and `plot.somalign_epsilon_sweep` are each
≤15 lines; each delegates to an internal helper.

---

## 11. Effort, risks, and dependencies

### Effort

Medium (estimated 3–4 hours implementation + 1 hour tests):

- `log_Z` propagation in `R/ot.R` + `R/fit.R`: ~25 lines changed/added.
- New helpers in `R/ot.R` (`.somalign_plan_entropy`, `.somalign_order_parameter`,
  `.somalign_log_Z_from_potentials`): ~25 lines.
- New helpers + public function in `R/diagnostics.R`: ~80 lines.
- `R/plot.R` (`.somalign_plot_eps_sweep`, S3 methods): ~35 lines.
- Tests (`test-epsilon-sweep.R`): ~80 lines.
- Roxygen docs: ~50 lines.

### Risks

1. **Washed-out transition on small codebooks.** For a 2×2 SOM (K = 4), Phi is
   already in (0.25, 1.0) and the range is small; `epsilon_c` may be at the
   grid boundary. Mitigated by returning the full `table` and the plot so the
   user can judge whether the peak is meaningful; no hard recommendation is
   forced when the susceptibility range is below a threshold (e.g.,
   `max(chi, na.rm = TRUE) < 0.01` → warning).
2. **`log_Z` sign convention.** The dual objective can be defined with various
   sign conventions in the UOT literature. The plan uses the formulation
   consistent with the log-domain solver's existing `lse_g` accumulation.
   Tests must verify that `log_Z` decreases (becomes more negative) as epsilon
   increases for a fixed problem — if the sign is backwards, flip by negation.
3. **Non-convergence at small epsilon biases susceptibility peak.** At very
   small epsilon the solver may fail to converge within `max_iter`, making Phi
   unreliable and the finite-difference noisy. Mitigation: the test suite checks
   that `converged = FALSE` rows are flagged; the plot can shade or cross-mark
   non-converged grid points; the user can raise `max_iter` or truncate the grid.

### Dependencies on other ideas

**Idea #2 (Simulated Annealing Sinkhorn)** also calls `.somalign_solve_internal_log()`
and needs the `f_init` / `g_init` warm-start interface. These two ideas share the
same internal solver but do not conflict: the warm-start arguments are optional
(`NULL` default) and the `log_Z` return value is independent of warm-starting.
They can be developed in parallel, but the `log_Z` addition to the return list
of `.somalign_solve_internal_log()` should land **before** Idea #2 so that the
annealing solver's final temperature also reports `log_Z` for free.

**Idea #3** (if it also performs an OT-only sweep without projection) should
reuse `.somalign_sweep_one_eps()` and `.somalign_run_grid()` rather than
re-implementing the loop pattern. The shared primitive is
`.somalign_sweep_one_eps(cost_norm, a, b, eps, ...)` — design it to be
general enough that Ideas #2/#3 can call it with different cost matrices
(e.g., anchor-modified costs) without modification.
