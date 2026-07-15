# Plan 03 — Mutual Information as Alignment-Sharpness Diagnostic + Epsilon Selector

## 1. Summary

This plan adds mutual information (MI) I(query node; reference node) as a
first-class alignment diagnostic derived from the M×K UOT transport plan.  Three
components are introduced:

1. **Per-fit scalar MI** (`diagnostics$ot$mutual_information`, bits) — alignment
   sharpness in one number, comparable across datasets and SOM resolutions.
2. **Per-node conditional entropy** (`diagnostics$nodes$transport_entropy`, bits)
   — flags nodes whose transport mass is spread across many reference nodes and
   whose label transfer is therefore unreliable, independent of match fraction.
3. **`somalign_select_epsilon()`** — sweeps an epsilon grid (OT only, no cell
   re-projection), computes MI at each point, draws the information-vs-cost
   curve, and returns the elbow epsilon via a second-difference rule.

All changes are additive.  No existing fields are removed or renamed.

---

## 2. Public API

### 2.1 `somalign_select_epsilon()`

**File**: `R/diagnostics.R`

```r
#' Select epsilon from a mutual-information curve
#'
#' Sweeps a grid of epsilon values, computing the OT plan and mutual information
#' at each point without running the full alignment pipeline.  Returns the
#' epsilon value at the elbow of the I(query; reference) vs. expected transport
#' cost curve.
#'
#' @param query    A `somalign_query` object.
#' @param reference A `somalign_reference` object.
#' @param epsilon  Numeric vector of candidate epsilon values (>0).
#'   Default `c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0)`.
#' @param rho_query,rho_ref  Marginal relaxations passed to the OT solver.
#' @param method   Character.  `"elbow"` (default) uses the max-second-difference
#'   rule on the sorted I-vs-cost curve.  `"entropy_fraction"` chooses the
#'   smallest epsilon where I >= `entropy_fraction * max(I)`.
#' @param entropy_fraction  Numeric in (0, 1].  Target fraction of maximum MI
#'   when `method = "entropy_fraction"`.  Default `0.90`.
#' @param solver,max_iter,tol  Sinkhorn solver parameters.
#' @param diagonal_boost  Passed to `.somalign_ot_sweep_one()`.
#' @param label_mask  Optional M×K logical penalty matrix (from
#'   `.somalign_build_label_mask()`).
#' @param parallel  Logical.  Use BiocParallel when `TRUE`.
#'
#' @return A list of class `somalign_epsilon_selection`:
#'   \describe{
#'     \item{`selected_epsilon`}{Numeric scalar: the recommended epsilon.}
#'     \item{`curve`}{Data frame with columns `epsilon`, `mutual_information`
#'       (bits), `expected_cost`, `conditional_entropy_mean` (mean per-node
#'       H(ref|query), bits).}
#'     \item{`method`}{Character: which selection method was used.}
#'   }
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(20), nrow=10, ncol=2,
#'               dimnames=list(NULL, c("F1","F2")))
#' ref <- somalign_train_reference(mat, grid=kohonen::somgrid(2,2,"hexagonal"),
#'                                 rlen=5)
#' qry <- somalign_query(mat, ref, grid=kohonen::somgrid(2,2,"hexagonal"),
#'                       rlen=5)
#' somalign_select_epsilon(qry, ref, epsilon = c(0.05, 0.1, 0.2))
#' @export
somalign_select_epsilon <- function(
  query,
  reference,
  epsilon      = c(0.02, 0.05, 0.1, 0.2, 0.5, 1.0),
  rho_query    = 1,
  rho_ref      = 1,
  method       = c("elbow", "entropy_fraction"),
  entropy_fraction = 0.90,
  solver       = c("internal", "log_domain", "auto"),
  max_iter     = 1000L,
  tol          = 1e-7,
  diagonal_boost = 0,
  label_mask   = NULL,
  parallel     = FALSE
)
```

**Return structure** (`somalign_epsilon_selection`):

```
list(
  selected_epsilon        = <scalar>,
  curve = data.frame(
    epsilon                 = <numeric vector, length E>,
    mutual_information      = <numeric, bits>,
    expected_cost           = <numeric, normalized cost units>,
    conditional_entropy_mean = <numeric, bits>
  ),
  method = <character>
)
```

Body delegates entirely to `.somalign_epsilon_sweep()` and
`.somalign_select_elbow()` / `.somalign_select_entropy_fraction()`, keeping the
exported body under 30 lines (well within the 50-line BiocCheck limit).

### 2.2 Changes to `somalign_diagnostics()` / `somalign_sensitivity_grid()`

`somalign_diagnostics()` is unchanged — it is a pass-through (`fit$diagnostics`)
at line 18 of `R/diagnostics.R`.  The new fields arrive via
`.somalign_build_diagnostics()` and `.somalign_build_nodes_diag()` in `R/fit.R`.

`somalign_sensitivity_grid()` gains one new column `mutual_information` in its
returned data frame.  This is additive and backward compatible.

---

## 3. Internal Helpers

All new helpers live in **`R/ot.R`** unless noted.

### 3.1 `.somalign_plan_mutual_information(plan)`

```r
# File: R/ot.R
# Compute I(query node; reference node) from the raw (possibly unnormalized)
# transport plan.  Returns a list with:
#   $mutual_information      – scalar (bits)
#   $conditional_entropy     – numeric vector length M (bits), one per query node
#   $expected_cost           – scalar (same units as `cost_normalized`, optional
#                              arg; NA if cost is not supplied)
.somalign_plan_mutual_information <- function(plan, cost = NULL) {
  tiny <- .Machine$double.eps
  total <- sum(plan)
  if (total <= tiny) {
    M <- nrow(plan)
    return(list(
      mutual_information  = NA_real_,
      conditional_entropy = rep(NA_real_, M),
      expected_cost       = NA_real_
    ))
  }
  P     <- plan / total                          # joint, sums to 1
  a_hat <- rowSums(P)                            # empirical query marginal
  b_hat <- colSums(P)                            # empirical reference marginal
  # MI: sum over nonzero P_ij of P_ij * log2(P_ij / (a_hat_i * b_hat_j))
  outer_ab <- outer(a_hat, b_hat)
  nz       <- P > tiny & outer_ab > tiny
  mi_bits  <- sum(P[nz] * log2(P[nz] / outer_ab[nz]))
  mi_bits  <- max(mi_bits, 0)                    # numerical noise floor
  # Per-row conditional entropy H(ref | query = i) = -sum_j P(j|i) log2 P(j|i)
  row_totals <- rowSums(plan)
  cond_ent <- vapply(seq_len(nrow(plan)), function(i) {
    if (row_totals[i] <= tiny) return(NA_real_)
    p_cond <- plan[i, ] / row_totals[i]          # P(ref | query=i)
    nz_j   <- p_cond > tiny
    if (!any(nz_j)) return(0)
    -sum(p_cond[nz_j] * log2(p_cond[nz_j]))
  }, numeric(1))
  exp_cost <- if (!is.null(cost)) sum(P * cost) else NA_real_
  list(
    mutual_information  = mi_bits,
    conditional_entropy = cond_ent,
    expected_cost       = exp_cost
  )
}
```

### 3.2 `.somalign_ot_sweep_one(query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol, diagonal_boost, label_mask)`

This is the **shared lightweight sweep primitive** also used by ideas #1 and #2.
It runs `.somalign_align_transport()` but skips cell reprojection, label
transfer, and node shifts — only the OT plan and summary statistics matter.

```r
# File: R/fit.R  (internal, not exported)
.somalign_ot_sweep_one <- function(query, reference, epsilon,
                                   rho_query, rho_ref,
                                   solver, max_iter, tol,
                                   diagonal_boost = 0,
                                   label_mask = NULL) {
  transport <- .somalign_align_transport(
    query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol,
    diagonal_boost = diagonal_boost, label_mask = label_mask
  )
  mi_result <- .somalign_plan_mutual_information(
    transport$plan, cost = transport$cost / transport$cost_scale
  )
  list(
    epsilon             = epsilon,
    plan                = transport$plan,
    cost_scale          = transport$cost_scale,
    mutual_information  = mi_result$mutual_information,
    conditional_entropy = mi_result$conditional_entropy,
    expected_cost       = mi_result$expected_cost,
    transport_mass      = sum(transport$plan),
    converged           = transport$ot$converged
  )
}
```

Ideas #1 and #2 can call `.somalign_ot_sweep_one()` directly for their own
epsilon sweeps, sharing all cost-computation and OT logic without duplication.

### 3.3 `.somalign_epsilon_sweep(query, reference, epsilon_grid, rho_query, rho_ref, solver, max_iter, tol, diagonal_boost, label_mask, parallel)`

```r
# File: R/diagnostics.R
.somalign_epsilon_sweep <- function(query, reference, epsilon_grid,
                                    rho_query, rho_ref, solver,
                                    max_iter, tol, diagonal_boost,
                                    label_mask, parallel) {
  .run_one <- function(i) {
    .somalign_ot_sweep_one(
      query, reference, epsilon_grid[i],
      rho_query, rho_ref, solver, max_iter, tol,
      diagonal_boost = diagonal_boost, label_mask = label_mask
    )
  }
  rows <- .somalign_run_grid(length(epsilon_grid), .run_one, parallel)
  data.frame(
    epsilon                  = epsilon_grid,
    mutual_information       = vapply(rows, `[[`, numeric(1), "mutual_information"),
    expected_cost            = vapply(rows, `[[`, numeric(1), "expected_cost"),
    conditional_entropy_mean = vapply(rows, function(r) {
      mean(r$conditional_entropy, na.rm = TRUE)
    }, numeric(1))
  )
}
```

### 3.4 `.somalign_select_elbow(curve)` and `.somalign_select_entropy_fraction(curve, fraction)`

```r
# File: R/diagnostics.R
# Kneedle-style: sort by expected_cost ascending; find index of maximum
# second difference of mutual_information.
.somalign_select_elbow <- function(curve) {
  ord  <- order(curve$expected_cost)
  mi   <- curve$mutual_information[ord]
  if (length(mi) < 3L || all(is.na(mi))) {
    return(curve$epsilon[ord[which.max(mi)]])
  }
  d2 <- diff(diff(mi))                # second differences (length E-2)
  # Elbow = point of maximum upward curvature (largest positive d2)
  elbow_idx <- which.max(d2) + 1L     # +1: d2[k] refers to mi[k+1]
  curve$epsilon[ord[elbow_idx]]
}

.somalign_select_entropy_fraction <- function(curve, fraction) {
  valid <- is.finite(curve$mutual_information)
  if (!any(valid)) return(NA_real_)
  target <- fraction * max(curve$mutual_information[valid])
  reached <- valid & curve$mutual_information >= target
  if (!any(reached)) return(curve$epsilon[which.max(curve$mutual_information)])
  # Return the largest (most-regularized) epsilon that still meets target
  curve$epsilon[max(which(reached))]
}
```

---

## 4. Data-Structure Changes

### 4.1 `diagnostics$ot` — new field

| Field | Type | Description |
|---|---|---|
| `mutual_information` | numeric scalar (bits) | I(query node; reference node) from the joint plan P̃ = plan / sum(plan). |

### 4.2 `diagnostics$nodes` — new column

| Column | Type | Description |
|---|---|---|
| `transport_entropy` | numeric, length M (bits) | H(ref node | query node = i) = per-row conditional entropy of the row-normalized plan.  NA for zero-mass nodes. |

A node with high `transport_entropy` exports mass to many reference nodes —
label transfer will be unreliable even if `match_fraction` is adequate.  Flag
threshold: `transport_entropy > log2(K / 4)` (i.e., the conditional distribution
is more diffuse than uniform over K/4 reference nodes).

### 4.3 Computation: exact formulas

Let P (M×K) be the raw plan from `.somalign_solve_ot()`.

```
P_total  = sum(P)                        # total transported mass
P_tilde  = P / P_total                   # joint distribution, sums to 1
a_hat_i  = rowSums(P_tilde)              # empirical query marginal
b_hat_j  = colSums(P_tilde)             # empirical reference marginal

MI = sum_{i,j: P_tilde_ij > 0}
       P_tilde_ij * log2(P_tilde_ij / (a_hat_i * b_hat_j))   [bits]

H_i = -sum_{j: P_ij > 0}
        (P_ij / sum_j P_ij) * log2(P_ij / sum_j P_ij)        [bits, per row]
```

Note: marginals are *empirical* (from the plan), not the input masses `a`, `b`.
This is correct for UOT where mass is destroyed — the empirical marginals
reflect what was actually transported.

### 4.4 `somalign_sensitivity_grid()` — new column

The returned data frame gains `mutual_information` (bits) per grid row, computed
inside `.somalign_grid_row_summary()`.

---

## 5. Algorithm and Pseudocode

### 5.1 MI at a single epsilon (`.somalign_plan_mutual_information`)

```r
# Already given in Section 3.1.
# Key: floor with .Machine$double.eps before any log; cap MI at 0 from below.
```

### 5.2 Information-rate-vs-cost curve

For each epsilon in the grid (ascending), `.somalign_ot_sweep_one()` returns:
- `mutual_information` = I(query; reference) in bits
- `expected_cost` = E_P[C_normalized] = sum(P_tilde * C_norm)

Plot: x = expected_cost, y = mutual_information.  As epsilon ↑ the plan
diffuses: cost ↑ (mass is transported farther on average) and MI ↓ (plan
approaches product marginal).  As epsilon ↓ the plan sharpens: cost ↓ (mass
concentrates on cheap pairs) and MI ↑ (but overfits SOM topology noise).

The optimal operating point is the elbow of this Pareto frontier.

### 5.3 Elbow detection (Kneedle-style, no external package)

```r
# Sort curve by expected_cost (ascending = decreasing epsilon usually)
# Compute second differences of MI on that sorted sequence
# Elbow = index of maximum second difference (greatest deceleration of MI gain)
# Return curve$epsilon at that index
```

For grids with < 3 valid points, return the epsilon at maximum MI.

### 5.4 `entropy_fraction` rule

```r
target <- entropy_fraction * max(curve$mutual_information, na.rm = TRUE)
# Return the largest epsilon that still achieves >= target MI
# (most regularized acceptable point — conservative, less grid-sensitive)
```

---

## 6. Integration Points

### 6.1 `fit.R` — `.somalign_build_diagnostics()` (line 264–307)

Add one call to `.somalign_plan_mutual_information()` and plumb results into
both `diagnostics$ot` and `diagnostics$nodes`:

```r
# After computing plan, row_mass, col_mass (lines 267-270):
mi_result <- .somalign_plan_mutual_information(plan)

# In the `ot = list(...)` block (lines 291-300), append:
mutual_information = mi_result$mutual_information,

# In `.somalign_build_nodes_diag()` (line 309-318), add column:
transport_entropy = mi_result$conditional_entropy
```

`.somalign_build_nodes_diag()` currently receives `query`, `transport`,
`node_shifts`; it must also receive `mi_result`.  Signature change:

```r
# Old (line 309):
.somalign_build_nodes_diag <- function(query, transport, node_shifts)

# New:
.somalign_build_nodes_diag <- function(query, transport, node_shifts, mi_result)
```

Call site at line 301:
```r
nodes = .somalign_build_nodes_diag(query, transport, node_shifts, mi_result),
```

### 6.2 `diagnostics.R` — `.somalign_grid_row_summary()` (lines 123–138)

```r
# Add after `outside_corrected_fraction` line:
mutual_information = .somalign_plan_mutual_information(
  fit$transport_plan
)$mutual_information,
```

The `fit$transport_plan` field is already set by `.somalign_new_fit()` (line 353
of `fit.R`), so no structural change is needed in the grid function.

### 6.3 `diagnostics.R` — `somalign_sensitivity_grid()` cost-free OT sweep

The current grid function calls `.somalign_finish_fit()` which runs cell
reprojection.  For `somalign_select_epsilon()` this is too expensive (no cell
data needed for the MI curve).  The shared `.somalign_ot_sweep_one()` primitive
(Section 3.2) avoids that cost.  The existing `somalign_sensitivity_grid()`
continues to call `.somalign_finish_fit()` (unchanged) and merely picks up the
new `mutual_information` field from the diagnostics.

---

## 7. Edge Cases

| Case | Handling |
|---|---|
| `sum(plan) == 0` (all mass destroyed) | Return `NA_real_` for MI and `NA_real_` for all `conditional_entropy` entries. |
| Row i of plan is all-zero (zero-mass query node) | `conditional_entropy[i] = NA_real_`. Do not divide by zero. |
| Single reference node (K=1) | `b_hat` is scalar 1; outer product is `a_hat`; `P_tilde / outer` = `P_tilde / a_hat` row-wise; MI = 0 (no information about reference node). `conditional_entropy` is 0 for every row (only one destination). |
| Single query node (M=1) | `a_hat` scalar 1; MI = H(b_hat) (entropy of reference marginal under plan). `conditional_entropy[1]` = same. |
| `log(0)` | Guard with `nz <- P_tilde > tiny` before any `log2` call; 0 * log(0) contributes 0 by convention. |
| Plan entries not summing to 1 | Normalize explicitly via `P / sum(P)` inside helper; do not rely on caller. |
| `epsilon_grid` length 1 or 2 in `somalign_select_epsilon` | Second-difference elbow undefined; fall back to returning `epsilon_grid[which.max(mi)]` with a `message()`. |
| Negative MI from floating-point noise | `mi_bits <- max(mi_bits, 0)` as final guard. |

---

## 8. Tests

File: `tests/testthat/test-mutual-information.R`

### 8.1 MI = 0 for uniform plan

```r
test_that(".somalign_plan_mutual_information returns 0 for product-marginal plan", {
  # Uniform plan: P_ij = a_i * b_j  =>  MI = 0
  a <- c(0.5, 0.3, 0.2)
  b <- c(0.4, 0.4, 0.2)
  plan <- outer(a, b)  # already normalized to sum 1
  res  <- .somalign_plan_mutual_information(plan)
  expect_equal(res$mutual_information, 0, tolerance = 1e-10)
})
```

### 8.2 MI is maximal for a permutation plan

```r
test_that(".somalign_plan_mutual_information is maximized for a permutation plan", {
  # Perfect bijection: all mass on diagonal
  n <- 4L
  plan_perm <- diag(rep(0.25, n))   # permutation plan, uniform marginals
  plan_unif <- matrix(1/n^2, n, n)   # product of uniform marginals
  res_perm <- .somalign_plan_mutual_information(plan_perm)
  res_unif <- .somalign_plan_mutual_information(plan_unif)
  expect_gt(res_perm$mutual_information, res_unif$mutual_information)
  # For uniform marginals, max MI = log2(n)
  expect_equal(res_perm$mutual_information, log2(n), tolerance = 1e-10)
})
```

### 8.3 Conditional entropy flags an ambiguous node

```r
test_that("transport_entropy is high for a diffusely-mapped node", {
  # Node 1: all mass to reference node 1 (sharp) => H = 0
  # Node 2: mass split evenly over 4 reference nodes (diffuse) => H = log2(4) = 2 bits
  plan <- matrix(0, 2, 4)
  plan[1, 1] <- 1.0
  plan[2, ] <- rep(0.25, 4)
  res <- .somalign_plan_mutual_information(plan)
  expect_equal(res$conditional_entropy[1], 0, tolerance = 1e-10)
  expect_equal(res$conditional_entropy[2], log2(4), tolerance = 1e-10)
})
```

### 8.4 Zero-plan returns NA

```r
test_that(".somalign_plan_mutual_information returns NA for all-zero plan", {
  plan <- matrix(0, 3, 3)
  res  <- .somalign_plan_mutual_information(plan)
  expect_true(is.na(res$mutual_information))
  expect_true(all(is.na(res$conditional_entropy)))
})
```

### 8.5 `somalign_select_epsilon` returns a value within the grid

```r
test_that("somalign_select_epsilon returns an epsilon within the supplied grid", {
  set.seed(42)
  mat <- matrix(rnorm(60), 30, 2, dimnames = list(NULL, c("F1","F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"),
                                  rlen = 10)
  qry <- somalign_query(mat + 0.3, ref,
                        grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10)
  eps_grid <- c(0.05, 0.1, 0.2, 0.5)
  sel <- somalign_select_epsilon(qry, ref, epsilon = eps_grid)
  expect_s3_class(sel, "somalign_epsilon_selection")
  expect_true(sel$selected_epsilon %in% eps_grid)
  expect_equal(nrow(sel$curve), length(eps_grid))
  expect_true(all(c("epsilon","mutual_information","expected_cost") %in%
                    names(sel$curve)))
})
```

### 8.6 `diagnostics$ot$mutual_information` is present after `somalign_fit()`

```r
test_that("somalign_fit diagnostics include mutual_information scalar", {
  set.seed(1)
  mat <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("F1","F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2,2,"hexagonal"),
                                  rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2,2,"hexagonal"),
                        rlen = 5)
  fit  <- somalign_fit(qry, ref)
  diag <- somalign_diagnostics(fit)
  expect_true("mutual_information" %in% names(diag$ot))
  expect_true(is.numeric(diag$ot$mutual_information))
  expect_true("transport_entropy" %in% names(diag$nodes))
  expect_equal(length(diag$nodes$transport_entropy), nrow(qry$codebook))
})
```

### 8.7 `somalign_sensitivity_grid` includes `mutual_information` column

```r
test_that("somalign_sensitivity_grid includes mutual_information column", {
  set.seed(1)
  mat <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("F1","F2")))
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2,2,"hexagonal"),
                                  rlen = 5)
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2,2,"hexagonal"),
                        rlen = 5)
  grid <- somalign_sensitivity_grid(qry, ref,
                                    epsilon = c(0.1, 0.5),
                                    rho_query = 1, rho_ref = 1)
  expect_true("mutual_information" %in% names(grid))
})
```

---

## 9. Docs / NAMESPACE

### 9.1 NAMESPACE

Add one export line (in alphabetical order, after `somalign_results`):

```
export(somalign_select_epsilon)
```

No new `importFrom` lines are needed: `log2`, `sum`, `outer`, `vapply`, `diff`,
`which.max`, `max` are all base R.

### 9.2 Roxygen / man pages

- `somalign_select_epsilon.Rd` — generated from the roxygen block in Section 2.1.
- `somalign_diagnostics.Rd` — add to `@return` a note that `$ot$mutual_information`
  (scalar, bits) and `$nodes$transport_entropy` (per-node, bits) are included.
- `somalign_sensitivity_grid.Rd` — add `mutual_information` to the `@return`
  data frame description.

---

## 10. BiocCheck (<=50-line exported bodies)

`somalign_select_epsilon()` body:

```r
somalign_select_epsilon <- function(...) {
  # ~10 lines: argument validation via existing .somalign_check_* helpers
  .somalign_check_query(query)
  .somalign_check_reference(reference)
  epsilon  <- .somalign_validate_grid_vector(epsilon, "epsilon")
  method   <- match.arg(method, c("elbow", "entropy_fraction"))
  solver   <- match.arg(solver, c("internal", "log_domain", "auto"))
  .somalign_check_pos_scalar(entropy_fraction, "entropy_fraction")
  .somalign_check_flag(parallel, "parallel")

  # ~5 lines: sweep
  curve <- .somalign_epsilon_sweep(
    query, reference, epsilon, rho_query, rho_ref,
    solver, max_iter, tol, diagonal_boost, label_mask, parallel
  )

  # ~8 lines: select
  selected <- if (identical(method, "elbow")) {
    .somalign_select_elbow(curve)
  } else {
    .somalign_select_entropy_fraction(curve, entropy_fraction)
  }

  # ~6 lines: return
  structure(
    list(selected_epsilon = selected, curve = curve, method = method),
    class = "somalign_epsilon_selection"
  )
}
# Total: ~29 lines. Well under 50.
```

`somalign_diagnostics()` body remains 5 lines (unchanged).

---

## 11. Effort, Risks, and Dependencies

### 11.1 Effort

| Component | Lines of new R code (approx.) | Time |
|---|---|---|
| `.somalign_plan_mutual_information()` | ~35 | 1 h |
| `.somalign_ot_sweep_one()` | ~20 | 0.5 h |
| `.somalign_epsilon_sweep()` | ~15 | 0.5 h |
| `.somalign_select_elbow()` + `_entropy_fraction()` | ~20 | 0.5 h |
| `somalign_select_epsilon()` exported wrapper | ~30 | 0.5 h |
| Integration in `fit.R` (`_build_diagnostics`, `_build_nodes_diag`) | ~10 | 0.5 h |
| Integration in `diagnostics.R` (`_grid_row_summary`) | ~5 | 0.25 h |
| Tests (7 test cases) | ~80 | 1.5 h |
| Docs + NAMESPACE | ~20 | 0.5 h |
| **Total** | **~235** | **~5.75 h** |

Overall rating: **Low–Medium** (consistent with the idea document's "Low" rating;
slightly higher due to the elbow detection and shared sweep primitive).

### 11.2 Risks

| Risk | Severity | Mitigation |
|---|---|---|
| Near-zero plan entries generate `-Inf` in `log2` | Medium | Guard with `nz <- P > .Machine$double.eps` before every log call. |
| Second-difference elbow is noisy on small grids (< 5 points) | Low | Fall back to `which.max(MI)` with a `message()`; document that ≥5 grid points are recommended. |
| UOT mass destruction makes `sum(plan)` << 1 | Low | Normalize by `sum(plan)`, not by 1, inside helper; explicitly documented. |
| Confusion between per-node `transport_entropy` and `label_transfer$entropy` | Medium | Name column `transport_entropy` (not `entropy`) and document the distinction in `?somalign_diagnostics`. |

### 11.3 Shared primitive with ideas #1 and #2

**`.somalign_ot_sweep_one()` signature (canonical, to be used by ideas #1, #2, #3)**:

```r
.somalign_ot_sweep_one(
  query,          # somalign_query
  reference,      # somalign_reference
  epsilon,        # scalar
  rho_query,      # scalar
  rho_ref,        # scalar
  solver,         # character
  max_iter,       # integer
  tol,            # numeric
  diagonal_boost = 0,
  label_mask     = NULL
) -> list(
  epsilon, plan, cost_scale, mutual_information,
  conditional_entropy, expected_cost, transport_mass, converged
)
```

Ideas #1 and #2 call this same function and extract whichever fields they need
(e.g. idea #1 may extract `plan` for subspace diagnostics; idea #2 may use
`transport_mass`).  This avoids three independent reimplementations of the
lightweight epsilon sweep.
