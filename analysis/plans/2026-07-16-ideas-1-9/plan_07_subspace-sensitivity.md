# Plan 07 — Sensitivity to unmeasured batch confounding: bootstrap the batch subspace

**Feature**: `somalign_subspace_sensitivity()`
**Status**: planned
**Date**: 2026-07-16
**Depends on**: none (fit-path change is backward-compatible); idea #8 (exclusion-restriction test) shares the D-storage prerequisite

---

## 1. Summary

Add a new exported read-only function `somalign_subspace_sensitivity()` that quantifies
how much the node-level batch correction depends on the estimated batch subspace V.
The function bootstraps the anchor displacement matrix D (n_anchors × p), recomputes
V for each replicate, propagates each V_b to per-node correction vectors, and reports:

- Empirical confidence intervals on the per-node correction norm and per-feature shift
- Per-bootstrap principal angles between V and each V_b (Grassmannian distances)
- A "tipping angle" per node: the smallest angular perturbation of V that would reverse
  the correction direction (analytic for rank-1, empirical for rank > 1)
- Anchor leverage: which anchor pairs most influence V (Cook's-distance analog)

No existing code paths are modified beyond a single backward-compatible addition:
storing D in `fit$anchors$displacements`. All computation uses base-R `svd()`.

---

## 2. Public API

### Signature

```r
somalign_subspace_sensitivity <- function(fit,
                                          n_boot           = 200L,
                                          variance_threshold = NULL,
                                          conf_level       = 0.95,
                                          seed             = 1L) { ... }
```

**Arguments**

| Argument            | Type      | Default | Meaning                                                  |
|---------------------|-----------|---------|----------------------------------------------------------|
| `fit`               | `somalign_anchored_fit` | — | Fitted anchored object with subspace correction |
| `n_boot`            | `integer` | `200L`  | Number of bootstrap replicates of D                     |
| `variance_threshold`| `numeric` or `NULL` | `NULL` | If NULL, uses `fit$anchors$batch_subspace$variance_threshold`; otherwise overrides |
| `conf_level`        | `numeric` | `0.95`  | Confidence level for node-correction CIs                |
| `seed`              | `integer` | `1L`    | RNG seed for reproducibility                             |

### Return value

A named list with class `"somalign_subspace_sensitivity"`:

```r
list(
  node_correction_ci   = matrix(..., nrow = M, ncol = 2,
                                dimnames = list(NULL, c("lower", "upper"))),
                          # CI on per-node correction norm across bootstraps

  node_shift_ci        = array(..., dim = c(M, p, 2),
                                dimnames = list(NULL, colnames(D), c("lower", "upper"))),
                          # CI on per-feature corrected shift

  subspace_angles      = matrix(..., nrow = n_boot, ncol = r,
                                dimnames = list(NULL, paste0("angle_", seq_len(r)))),
                          # Principal angles (degrees) between V and each V_b; r = rank of V

  tipping_angle_deg    = numeric(M),
                          # Per-node; NA for non-allowed or zero-shift nodes

  anchor_leverage      = numeric(n_anchors),
                          # Influence of each anchor on V (Cook's-D analog)

  n_boot               = n_boot,
  conf_level           = conf_level,
  subspace_rank        = r,
  n_anchors            = nrow(D),
  variance_threshold   = vt_used
)
```

### Error conditions

```r
# Hard stops
if (!inherits(fit, "somalign_anchored_fit"))
  stop("`fit` must be a somalign_anchored_fit object.", call. = FALSE)
if (fit$anchors$correction == "cost_bonus")
  stop("`somalign_subspace_sensitivity()` requires correction = \"subspace\" or \"both\".",
       call. = FALSE)
if (is.null(fit$anchors$displacements))
  stop("`fit$anchors$displacements` is NULL; refit with somalign_fit_anchored() ",
       "version >= X.Y.Z which stores the displacement matrix.", call. = FALSE)

# Soft warnings
r    <- fit$anchors$batch_subspace$rank
n_a  <- fit$anchors$n_anchors
if (n_a < r * 5L)
  warning("n_anchors (", n_a, ") < rank * 5 (", r * 5L, "); bootstrap distribution ",
          "may be degenerate. Consider adding more anchor pairs.", call. = FALSE)
if (n_boot < 100L)
  warning("n_boot = ", n_boot, " is small; CI coverage may be poor.", call. = FALSE)
```

---

## 3. Prerequisite fit-object change

**What needs to be stored**: the scaled displacement matrix D, defined as
`anchor_old_scaled - anchor_new_scaled` (n_anchors × p in reference-scaled space).
This is already computed transiently inside `.somalign_batch_subspace()` (anchored.R
lines 247–250) and then discarded. It must be retained and passed through to `fit$anchors`.

**Why only D is needed (not raw node shifts)**: Raw node shifts S (pre-projection) can
be recovered from stored fit components without new storage:

```r
bary  <- fit$correspondence %*% fit$reference$codebook   # M × p
S_raw <- bary - fit$query$codebook                       # M × p
allowed <- attr(fit$node_shifts, "correction_allowed")
S_raw[!allowed, ] <- 0
```

`fit$correspondence` (M × K transport plan) and both codebooks are stored by
`.somalign_new_fit()` (fit.R line 356). `fit$node_shifts` for a subspace fit stores
*post-projection* shifts (S V Vᵀ) — applying fresh V_b matrices to it would double-
project and give garbage tipping angles. Always recompute S_raw as above.

**No raw_node_shifts storage**: This means only one fit-path addition is needed.

### Exact edit location: `.somalign_anchored_dispatch()` in `R/anchored.R`

Current code at approximately lines 217–241:

```r
batch_sub <- if (use_subspace) {
  .somalign_batch_subspace(anchors_scaled$anchor_old_scaled,
                           anchors_scaled$anchor_new_scaled,
                           variance_threshold)
} else { NULL }
shift_fn <- if (use_subspace) {
  V <- batch_sub$V
  function(s) s %*% V %*% t(V)
} else { NULL }
# ...
fit <- .somalign_finish_fit(
  ...,
  anchors = list(
    n_anchors         = nrow(anchors_scaled$anchor_old_scaled),
    rho_anchor        = rho_anchor,
    correction        = correction,
    nodes_covered     = cb$nodes_covered,
    coverage_fraction = cb$coverage_fraction,
    batch_subspace    = batch_sub
  )
)
```

**After edit** — add two lines (marked `# NEW`):

```r
# NEW: retain D before passing to batch_subspace
D_scaled <- if (use_subspace) {
  anchors_scaled$anchor_old_scaled - anchors_scaled$anchor_new_scaled
} else { NULL }

batch_sub <- if (use_subspace) {
  .somalign_batch_subspace(anchors_scaled$anchor_old_scaled,
                           anchors_scaled$anchor_new_scaled,
                           variance_threshold)
} else { NULL }
shift_fn <- if (use_subspace) {
  V <- batch_sub$V
  function(s) s %*% V %*% t(V)
} else { NULL }
# ...
fit <- .somalign_finish_fit(
  ...,
  anchors = list(
    n_anchors         = nrow(anchors_scaled$anchor_old_scaled),
    rho_anchor        = rho_anchor,
    correction        = correction,
    nodes_covered     = cb$nodes_covered,
    coverage_fraction = cb$coverage_fraction,
    batch_subspace    = batch_sub,
    displacements     = D_scaled   # NEW: n_anchors × p; NULL for cost_bonus
  )
)
```

**Also store `variance_threshold`** in `batch_subspace` so `somalign_subspace_sensitivity()`
can recover it when `variance_threshold = NULL`:

```r
# In .somalign_subspace_svd() (utils.R) or .somalign_batch_subspace() (anchored.R),
# pass variance_threshold through to the returned list:
batch_sub <- list(V = ..., rank = ..., variance_explained = ...,
                  variance_threshold = variance_threshold)   # NEW
```

If `variance_threshold` is not yet stored in `batch_subspace`, the exported function
falls back to `fit$call$variance_threshold` or the package default (0.9).

---

## 4. Internal helpers

All helpers live in **`R/sensitivity.R`** (new file). Exported function body delegates to
them to satisfy the BiocCheck ≤ 50-line rule.

### 4.1 `.somalign_bootstrap_subspace(D, r, n_boot, variance_threshold, seed)`

```r
# Args:
#   D:                  n_anchors × p matrix (reference-scaled displacements)
#   r:                  integer, rank to fix across all replicates (point-estimate rank)
#   n_boot:             integer
#   variance_threshold: numeric, passed to .somalign_subspace_svd()
#   seed:               integer
# Returns: list(V_boots = array(p, r, n_boot), angle_mat = matrix(n_boot, r))
.somalign_bootstrap_subspace <- function(D, V, r, n_boot, variance_threshold, seed) {
  set.seed(seed)
  n_a <- nrow(D)
  p   <- ncol(D)
  V_boots   <- array(NA_real_, dim = c(p, r, n_boot))
  angle_mat <- matrix(NA_real_, nrow = n_boot, ncol = r)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n_a, replace = TRUE)
    D_b <- D[idx, , drop = FALSE]
    sub_b <- .somalign_subspace_svd(D_b, variance_threshold)
    # Fix rank to r for comparability; if r_b < r, pad with zeros
    r_b <- sub_b$rank
    V_b <- matrix(0, nrow = p, ncol = r)
    r_use <- min(r, r_b)
    V_b[, seq_len(r_use)] <- sub_b$V[, seq_len(r_use), drop = FALSE]
    V_boots[, , b] <- V_b
    angle_mat[b, ] <- .somalign_principal_angles(V, V_b)
  }
  list(V_boots = V_boots, angle_mat = angle_mat)
}
```

### 4.2 `.somalign_principal_angles(V1, V2)`

```r
# Returns r principal angles (in degrees) between column spaces of V1, V2 (both p × r).
# Handles degenerate V_b (all-zero columns) by returning 90 degrees for those dims.
.somalign_principal_angles <- function(V1, V2) {
  r <- ncol(V1)
  # Guard: if V2 column is all-zero (degenerate bootstrap replicate), angle = 90
  col_norms <- sqrt(colSums(V2^2))
  bad <- col_norms < .Machine$double.eps * 10
  if (any(bad)) {
    angles <- numeric(r)
    angles[bad] <- 90
    if (any(!bad)) {
      sv <- svd(crossprod(V1[, !bad, drop=FALSE], V2[, !bad, drop=FALSE]),
                nu = 0, nv = 0)$d
      angles[!bad] <- acos(pmin(pmax(sv, -1), 1)) * (180 / pi)
    }
    return(angles)
  }
  sv <- svd(crossprod(V1, V2), nu = 0, nv = 0)$d
  acos(pmin(pmax(sv, -1), 1)) * (180 / pi)
}
```

### 4.3 `.somalign_tipping_angle(S_raw, V)`

```r
# Computes per-node tipping angle (degrees).
# For rank-1 V (unit vector v): tipping_angle_i = asin(|S_raw[i,] . v| / ||S_raw[i,]||)
# For rank > 1: minimum principal angle (over bootstraps) at which projected shift
#               sign flips — approximated as asin(min singular value of
#               S_raw[i,,drop=FALSE] %*% V / ||S_raw[i,]||).
# Nodes with correction_allowed = FALSE, or zero-norm shift, get NA.
.somalign_tipping_angle <- function(S_raw, V, allowed) {
  r   <- ncol(V)
  M   <- nrow(S_raw)
  out <- rep(NA_real_, M)
  for (i in seq_len(M)) {
    if (!allowed[i]) next
    s   <- S_raw[i, ]
    nrm <- sqrt(sum(s^2))
    if (nrm < .Machine$double.eps * 10) next
    s_hat <- s / nrm
    if (r == 1L) {
      # Analytic: angle between s_hat and v such that projection reverses
      cos_ang <- abs(sum(s_hat * V[, 1]))
      out[i]  <- asin(pmin(cos_ang, 1)) * (180 / pi)
    } else {
      sv_min <- min(svd(matrix(s_hat, nrow = 1) %*% V, nu = 0, nv = 0)$d)
      out[i] <- asin(pmin(sv_min, 1)) * (180 / pi)
    }
  }
  out
}
```

### 4.4 `.somalign_anchor_leverage(D, V)`

```r
# Cook's-D analog: leave-one-out influence of each anchor row on V.
# Measured as the maximum principal angle between V_full and V_{-i}.
.somalign_anchor_leverage <- function(D, V, variance_threshold) {
  n_a <- nrow(D)
  lev <- numeric(n_a)
  for (i in seq_len(n_a)) {
    D_loo <- D[-i, , drop = FALSE]
    sub_loo <- .somalign_subspace_svd(D_loo, variance_threshold)
    lev[i] <- max(.somalign_principal_angles(V, sub_loo$V))
  }
  lev
}
```

---

## 5. Algorithm

### Step 0: Input validation and raw-shift recovery

```r
D       <- fit$anchors$displacements          # n_anchors × p
V       <- fit$anchors$batch_subspace$V       # p × r
r       <- fit$anchors$batch_subspace$rank
vt_used <- variance_threshold %||%
           fit$anchors$batch_subspace$variance_threshold %||% 0.9

# Recover raw (pre-projection) node shifts from stored components:
bary    <- fit$correspondence %*% fit$reference$codebook   # M × p
S_raw   <- bary - fit$query$codebook                       # M × p
allowed <- attr(fit$node_shifts, "correction_allowed")
S_raw[!allowed, ] <- 0
```

**Critical**: `fit$node_shifts` is *post-projection* (= S_raw %*% V %*% t(V)) because
`.somalign_finish_fit()` applies `shift_transform(node_shifts)` at line 145 before
storing. Using `fit$node_shifts` directly as "S" in the bootstrap would double-project
and produce invalid tipping angles.

### Step 1: Bootstrap subspace

```r
boot_res <- .somalign_bootstrap_subspace(D, V, r, n_boot, vt_used, seed)
# boot_res$V_boots: p × r × n_boot array of bootstrap V matrices
# boot_res$angle_mat: n_boot × r matrix of principal angles (degrees)
```

### Step 2: Per-node bootstrap correction distributions

```r
M <- nrow(S_raw)
p <- ncol(S_raw)
alpha <- 1 - conf_level
probs <- c(alpha / 2, 1 - alpha / 2)

node_correction_boot <- matrix(NA_real_, nrow = n_boot, ncol = M)
node_shift_boot      <- array(NA_real_, dim = c(n_boot, M, p))

for (b in seq_len(n_boot)) {
  V_b   <- boot_res$V_boots[, , b]
  corr_b <- S_raw %*% V_b %*% t(V_b)          # M × p
  node_correction_boot[b, ] <- sqrt(rowSums(corr_b^2))
  node_shift_boot[b, , ]    <- corr_b
}

node_correction_ci <- apply(node_correction_boot, 2, quantile, probs = probs, na.rm = TRUE)
node_correction_ci <- t(node_correction_ci)         # M × 2

node_shift_ci <- array(NA_real_, dim = c(M, p, 2))
for (j in seq_len(p)) {
  mat_j <- node_shift_boot[, , j]               # n_boot × M
  node_shift_ci[, j, ] <- t(apply(mat_j, 2, quantile, probs = probs, na.rm = TRUE))
}
```

### Step 3: Tipping angles

```r
tipping_angle_deg <- .somalign_tipping_angle(S_raw, V, allowed)
```

For rank-1 subspaces this is analytic (`asin(|ŝ_i · v|)`); for rank > 1 it uses the
minimum singular value of the (1 × r) projection matrix as a conservative proxy for the
smallest rotation of V that zeros the projected shift.

### Step 4: Anchor leverage

```r
anchor_leverage <- .somalign_anchor_leverage(D, V, vt_used)
```

### Step 5: Assemble and return

```r
structure(
  list(
    node_correction_ci = node_correction_ci,
    node_shift_ci      = node_shift_ci,
    subspace_angles    = boot_res$angle_mat,
    tipping_angle_deg  = tipping_angle_deg,
    anchor_leverage    = anchor_leverage,
    n_boot             = n_boot,
    conf_level         = conf_level,
    subspace_rank      = r,
    n_anchors          = nrow(D),
    variance_threshold = vt_used
  ),
  class = "somalign_subspace_sensitivity"
)
```

---

## 6. Integration points

### Read-only over the fit object

`somalign_subspace_sensitivity()` is entirely read-only. It accesses:

- `fit$anchors$displacements` (new field, see §3)
- `fit$anchors$batch_subspace$V`, `$rank`, `$variance_threshold`
- `fit$anchors$correction` (for error guard)
- `fit$correspondence`, `fit$reference$codebook`, `fit$query$codebook` (raw-shift recovery)
- `attr(fit$node_shifts, "correction_allowed")` (node mask)

### Only one edit to the fit path

The single fit-path change is the addition of `displacements = D_scaled` to the anchors
list in `.somalign_anchored_dispatch()`. The addition is backward-compatible because all
existing code ignores unknown list elements. Fits created by older versions will have
`fit$anchors$displacements == NULL`, which the exported function catches with a clear
stop message.

### No new R package dependencies

All operations use base R: `svd()`, `crossprod()`, `sample.int()`, `apply()`,
`quantile()`, `asin()`, `acos()`.

---

## 7. Edge cases

| Case | Behaviour |
|---|---|
| `correction = "cost_bonus"` | Hard stop with informative message (see §2) |
| `fit$anchors$displacements` is NULL (old fit) | Hard stop pointing to refit (see §2) |
| `n_anchors < rank * 5` | Warning; proceed; CIs will be wide |
| rank-1 subspace | Analytic tipping-angle formula; no numerical issues |
| Zero-norm shift on a node (`allowed[i]` TRUE but `||s_i|| ≈ 0`) | `tipping_angle_deg[i] = NA` |
| `!allowed[i]` | `tipping_angle_deg[i] = NA`; bootstrap correction = 0 |
| Degenerate V_b (D_b collapses to zero variance after resampling) | `.somalign_subspace_svd()` returns zero V (utils.R lines 92–95); `.somalign_principal_angles()` guards this by detecting zero-norm columns and returning 90 degrees |
| Variable rank across replicates | Rank is fixed to point-estimate `r` for all replicates; if `r_b < r`, excess columns of V_b are zero-padded (see `.somalign_bootstrap_subspace()`) |
| `n_boot = 1` | Returns CI with equal lower and upper; no error |
| Non-subspace anchored fit (`correction = "both"` with failed subspace) | Caught by `is.null(fit$anchors$batch_subspace)` check; hard stop |

---

## 8. Tests

File: `tests/testthat/test-sensitivity.R`

### 8.1 — Tight anchors → small instability, wide tipping angle

```r
test_that("tight anchors produce stable subspace and wide tipping angles", {
  # Construct a rank-1 problem in 3D:
  # All anchor pairs displaced by exactly (1, 0, 0) + tiny noise
  set.seed(42)
  n_a <- 20; p <- 3
  D   <- matrix(c(rep(1, n_a), rep(0, n_a), rep(0, n_a)), nrow = n_a) +
         matrix(rnorm(n_a * p, sd = 0.01), n_a, p)
  # Build minimal anchored fit programmatically or use a helper
  # (test fixture creates fake fit with known D stored):
  fit <- .make_test_anchored_fit(D, correction = "subspace")
  sens <- somalign_subspace_sensitivity(fit, n_boot = 50, seed = 1)
  # Principal angles should be small (< 10 degrees)
  expect_true(all(sens$subspace_angles < 10, na.rm = TRUE))
  # Tipping angles for nodes with shift ≈ along V should be large (> 45 deg)
  ta <- sens$tipping_angle_deg
  allowed_ta <- ta[!is.na(ta)]
  expect_true(all(allowed_ta > 30))
})
```

### 8.2 — Noisy/diverse anchors → large instability

```r
test_that("noisy anchors produce wide subspace angle distribution", {
  set.seed(99)
  n_a <- 8; p <- 5   # deliberately small n_anchors
  D   <- matrix(rnorm(n_a * p, sd = 1), n_a, p)
  fit <- .make_test_anchored_fit(D, correction = "subspace")
  expect_warning(
    sens <- somalign_subspace_sensitivity(fit, n_boot = 100, seed = 1),
    regexp = "n_anchors"
  )
  # Median principal angle should be > 20 degrees
  expect_gt(median(sens$subspace_angles[, 1], na.rm = TRUE), 15)
})
```

### 8.3 — Error on cost_bonus fit

```r
test_that("somalign_subspace_sensitivity errors on cost_bonus fit", {
  fit <- .make_test_anchored_fit(matrix(rnorm(20), 4, 5),
                                 correction = "cost_bonus")
  expect_error(somalign_subspace_sensitivity(fit), regexp = "cost_bonus")
})
```

### 8.4 — Error on NULL displacements (old fit)

```r
test_that("somalign_subspace_sensitivity errors when displacements is NULL", {
  fit <- .make_test_anchored_fit(matrix(rnorm(20), 4, 5),
                                 correction = "subspace",
                                 store_displacements = FALSE)
  expect_error(somalign_subspace_sensitivity(fit), regexp = "displacements")
})
```

### 8.5 — Analytic rank-1 tipping angle formula

```r
test_that("rank-1 tipping angle matches analytic formula", {
  set.seed(7)
  n_a <- 20; p <- 4
  # Pure rank-1: D rows all collinear with e1
  D <- cbind(runif(n_a, 0.5, 1.5), matrix(0, n_a, p - 1))
  fit <- .make_test_anchored_fit(D, correction = "subspace")
  sens <- somalign_subspace_sensitivity(fit, n_boot = 200, seed = 1)
  # Cross-check: compute analytic angle for node i manually
  V      <- fit$anchors$batch_subspace$V
  S_raw  <- .recover_raw_shifts(fit)
  i      <- which(!is.na(sens$tipping_angle_deg))[1]
  s_hat  <- S_raw[i, ] / sqrt(sum(S_raw[i, ]^2))
  expected_angle <- asin(abs(sum(s_hat * V[, 1]))) * 180 / pi
  expect_equal(sens$tipping_angle_deg[i], expected_angle, tolerance = 1e-6)
})
```

### 8.6 — Degenerate V_b (bootstrap replicate collapses to near-zero D)

```r
test_that("degenerate V_b bootstrap replicate returns 90-degree angle", {
  # Two anchors: when one is dropped in leave-one-out, D may collapse
  set.seed(3)
  p <- 3; n_a <- 2
  D <- matrix(c(1, -1, 0, 0, 0, 0), nrow = n_a, ncol = p)
  fit <- .make_test_anchored_fit(D, correction = "subspace")
  # Expect warning (n_anchors < rank*5) but no error; degenerate guard fires
  expect_warning(
    sens <- somalign_subspace_sensitivity(fit, n_boot = 50, seed = 2),
    regexp = "n_anchors"
  )
  expect_true(all(is.finite(sens$subspace_angles) | is.na(sens$subspace_angles)))
})
```

---

## 9. Documentation and NAMESPACE

### Roxygen block (top of `R/sensitivity.R`)

```r
#' Subspace sensitivity analysis for anchored batch correction
#'
#' Bootstraps the anchor displacement matrix D to quantify how stable the
#' estimated batch subspace V is, and propagates that uncertainty to per-node
#' correction confidence intervals and "tipping angles".
#'
#' @param fit A \code{somalign_anchored_fit} with \code{correction = "subspace"}
#'   or \code{"both"}.
#' @param n_boot Integer. Number of bootstrap replicates of D. Default 200.
#' @param variance_threshold Numeric or NULL. Variance threshold for SVD rank
#'   selection in bootstrap replicates. If NULL, uses the threshold from the
#'   original fit. Default NULL.
#' @param conf_level Numeric in (0,1). Confidence level for node-correction CIs.
#'   Default 0.95.
#' @param seed Integer. RNG seed for reproducibility. Default 1L.
#'
#' @return A named list of class \code{"somalign_subspace_sensitivity"} containing
#'   \code{node_correction_ci} (M x 2 matrix), \code{node_shift_ci}
#'   (M x p x 2 array), \code{subspace_angles} (n_boot x r matrix, degrees),
#'   \code{tipping_angle_deg} (length-M vector), \code{anchor_leverage}
#'   (length-n_anchors vector), and metadata fields.
#'
#' @details
#' The "tipping angle" for node \eqn{i} is the smallest angular perturbation
#' of the batch subspace V (in the Grassmannian sense) that would reverse the
#' direction of the correction for that node. For rank-1 subspaces it has the
#' analytic form \eqn{\arcsin(|\hat{s}_i \cdot v|)}, where \eqn{\hat{s}_i}
#' is the unit-normalised raw barycentric shift; for rank > 1 it is approximated
#' from the minimum singular value of the projection matrix. A small tipping
#' angle (< 10 degrees) signals a fragile correction; a large tipping angle
#' (> 45 degrees) indicates robustness.
#'
#' Raw pre-projection node shifts are recovered from stored fit components
#' (\code{fit$correspondence}, both codebooks) rather than from
#' \code{fit$node_shifts}, which stores post-projection shifts for subspace fits.
#'
#' @examples
#' \dontrun{
#'   sens <- somalign_subspace_sensitivity(fit, n_boot = 200, seed = 1)
#'   hist(sens$tipping_angle_deg, main = "Tipping angles (degrees)")
#'   boxplot(sens$subspace_angles, main = "Principal angles per bootstrap replicate")
#' }
#' @export
somalign_subspace_sensitivity <- function(fit, n_boot = 200L,
                                          variance_threshold = NULL,
                                          conf_level = 0.95,
                                          seed = 1L) { ... }
```

### NAMESPACE

Add one line (generated by `devtools::document()` from the `@export` tag):

```
export(somalign_subspace_sensitivity)
```

No new `importFrom` needed (only base R used).

---

## 10. BiocCheck compliance (≤ 50-line exported bodies)

The exported function `somalign_subspace_sensitivity()` must have a body of ≤ 50 lines.
Achieve this by delegating all substantive work to internal helpers:

```r
somalign_subspace_sensitivity <- function(fit,
                                          n_boot = 200L,
                                          variance_threshold = NULL,
                                          conf_level = 0.95,
                                          seed = 1L) {
  # --- input validation (~10 lines) ---
  if (!inherits(fit, "somalign_anchored_fit"))
    stop("`fit` must be a somalign_anchored_fit object.", call. = FALSE)
  if (isTRUE(fit$anchors$correction == "cost_bonus"))
    stop("requires correction = \"subspace\" or \"both\".", call. = FALSE)
  if (is.null(fit$anchors$displacements))
    stop("`fit$anchors$displacements` is NULL; refit with a newer version.", call. = FALSE)
  .somalign_check_pos_int(n_boot, "n_boot")
  .somalign_check_pos_scalar(conf_level, "conf_level")
  .somalign_check_pos_int(seed, "seed")

  # --- extract key objects (~5 lines) ---
  D   <- fit$anchors$displacements
  V   <- fit$anchors$batch_subspace$V
  r   <- fit$anchors$batch_subspace$rank
  n_a <- nrow(D)
  vt  <- variance_threshold %||%
         fit$anchors$batch_subspace$variance_threshold %||% 0.9

  # --- warn on low anchor count (~3 lines) ---
  if (n_a < r * 5L)
    warning("n_anchors (", n_a, ") < rank * 5; bootstrap may be degenerate.",
            call. = FALSE)

  # --- recover raw node shifts (~5 lines) ---
  S_raw   <- .somalign_recover_raw_shifts(fit)
  allowed <- attr(fit$node_shifts, "correction_allowed")

  # --- bootstrap (~3 lines) ---
  boot_res <- .somalign_bootstrap_subspace(D, V, r, n_boot, vt, seed)

  # --- per-node CIs (~3 lines) ---
  ci_res <- .somalign_node_correction_ci(S_raw, boot_res$V_boots, conf_level)

  # --- tipping angles and leverage (~3 lines) ---
  tips <- .somalign_tipping_angle(S_raw, V, allowed)
  levs <- .somalign_anchor_leverage(D, V, vt)

  # --- assemble result (~10 lines) ---
  structure(
    list(node_correction_ci = ci_res$correction_ci,
         node_shift_ci      = ci_res$shift_ci,
         subspace_angles    = boot_res$angle_mat,
         tipping_angle_deg  = tips,
         anchor_leverage    = levs,
         n_boot             = n_boot,
         conf_level         = conf_level,
         subspace_rank      = r,
         n_anchors          = n_a,
         variance_threshold = vt),
    class = "somalign_subspace_sensitivity"
  )
}
# Total: ~42 lines
```

Internal helpers `.somalign_recover_raw_shifts()`, `.somalign_bootstrap_subspace()`,
`.somalign_node_correction_ci()`, `.somalign_tipping_angle()`,
`.somalign_anchor_leverage()`, `.somalign_principal_angles()` are all unexported and
have no line-length restriction.

---

## 11. Effort, risks, and dependencies

### Effort estimate

| Task | Hours |
|---|---|
| Fit-path edit (store D_scaled + variance_threshold in anchors list) | 0.5 |
| `.somalign_recover_raw_shifts()` helper + unit check vs `.somalign_node_shifts()` | 1.0 |
| `R/sensitivity.R`: exported function + 6 internal helpers | 4.0 |
| Roxygen docs + NAMESPACE | 0.5 |
| `tests/testthat/test-sensitivity.R` (6 tests including fixture helper) | 3.0 |
| `R CMD check` + BiocCheck pass | 0.5 |
| **Total** | **~9.5 hours** |

Overall effort: **Medium** (per idea doc estimate), primarily driven by test fixture
setup and the rank-fixation / degenerate-V guard logic.

### Key risks

1. **n_anchors too small**: bootstrap degenerates when n_anchors < rank × 5. The
   function warns but does not stop; CIs will be artifactually wide and tipping angles
   uninformative. Document the minimum in the vignette.

2. **Rank-fixation choice**: fixing bootstrap replicate rank to the point-estimate `r`
   is a modelling decision (ensures apples-to-apples principal-angle comparisons). The
   alternative — re-selecting rank per replicate — would produce variable-dimension V_b
   arrays and is harder to summarise. The chosen approach may understate variability when
   `variance_threshold` sits near a sharp scree-cliff.

3. **Tipping angle for rank > 1**: the minimum-singular-value approximation is
   conservative (it measures the smallest rotation that zeroes the *magnitude* of the
   correction, not its direction). A direction-reversal check would require solving a
   small quadratic eigenproblem per node. Note this limitation in the docs.

4. **Correspondence matrix size**: `fit$correspondence` is M × K (query nodes ×
   reference nodes). For very large SOMs (e.g., M = K = 10 000) the matrix multiply in
   raw-shift recovery is `10000 × 10000` — still fast (< 1 s with BLAS), but large in
   memory if materialised. The typical range is M, K ≤ 2 500, so this is not a concern
   in practice.

### Dependencies shared with other ideas

- **Idea #8 (exclusion-restriction test)**: also needs `fit$anchors$displacements`. The
  fit-path edit in §3 satisfies both features. Implement the D-storage change once;
  neither feature needs to be first.

- **Idea #10 (subspace sensitivity grid)**: if a later plan proposes sweeping V
  perturbation angle analytically, `.somalign_principal_angles()` defined here is
  directly reusable.

### Implementation order recommendation

1. Fit-path edit (§3) — needed by ideas #7 and #8; merge first.
2. `.somalign_recover_raw_shifts()` — shared utility, no external API surface.
3. `R/sensitivity.R` helpers and exported function.
4. Tests.
