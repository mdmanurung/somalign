# Plan 05 — Laplacian-Regularized Node Shifts

**Feature idea:** Idea 5 from `03_topologist-geometer.md`
**Status:** Draft
**Date:** 2026-07-16

---

## 1. Summary

The per-node correction vectors produced by `.somalign_node_shifts()` (`R/fit.R`,
line 449) are computed node-by-node with no spatial coupling; adjacent SOM nodes
can receive wildly different shifts from finite-sample OT noise. This plan wires
a graph-Laplacian smoother into the existing `shift_transform` hook
(`.somalign_finish_fit()`, lines 143–147) to post-process node shifts by solving
a Tikhonov-regularized system `(I + lambda * L) s_smooth = s` over the SOM
neighbor graph. The smoother penalises squared differences between adjacent
shifts, yielding a displacement field that is spatially coherent while still
fitting the OT-derived targets where data are dense. The feature is exposed via
`laplacian_lambda = 0` (default), which reproduces current behaviour exactly.

---

## 2. Public API

### 2.1 `somalign_fit()`

Add one argument to `somalign_fit()` (`R/fit.R`, line 76 signature):

```r
somalign_fit <- function(query,
                         reference,
                         ...,                   # existing args unchanged
                         laplacian_lambda = 0)
```

Validation (alongside the existing `check_*` calls at line 90):

```r
if (!is.numeric(laplacian_lambda) || length(laplacian_lambda) != 1L ||
    !is.finite(laplacian_lambda) || laplacian_lambda < 0)
  stop("`laplacian_lambda` must be a non-negative finite scalar.", call. = FALSE)
```

When `laplacian_lambda > 0`, build the Laplacian and construct the
`shift_transform` closure before calling `.somalign_finish_fit()`:

```r
shift_transform <- if (laplacian_lambda > 0) {
  L <- .somalign_som_laplacian(query$som_query$grid)
  masses <- query$node_masses
  function(s) .somalign_smooth_shifts(s, L, laplacian_lambda, masses)
} else {
  NULL   # exact current behaviour
}

.somalign_finish_fit(
  query, reference, transport,
  min_match_fraction, confidence_threshold, correction_min_mass,
  chunk_size, epsilon, rho_query, rho_ref,
  shift_transform = shift_transform
)
```

### 2.2 `somalign_fit_anchored()` (optional, recommended)

Add the same `laplacian_lambda = 0` argument to
`somalign_fit_anchored()` (`R/anchored.R`, line 153) and thread it through
`.somalign_anchored_dispatch()`. See Section 6 for composition with the subspace
transform.

### 2.3 Default and scale

`laplacian_lambda = 0` is the default. No existing call site changes.

The parameter shares the same scale as `epsilon` (both live in cost/squared-distance
space post-normalization). A natural starting range is `0.1`–`1.0`; documentation
should mention this analogy.

### 2.4 Composition with `shift_transform` in `somalign_fit_anchored()`

The subspace mode in `somalign_fit_anchored()` already sets:

```r
shift_fn <- function(s) s %*% V %*% t(V)   # (anchored.R line 224)
```

Laplacian smoothing and subspace projection are **independent** linear
post-processors on `node_shifts`. Their composition order matters:

- **Smooth first, then project** (`project(smooth(s))`): smoothing is done in
  full marker space before restricting to the batch subspace V. This is the
  correct order because the smoother should operate over all features (reflecting
  the true geometry of marker space), and then the subspace projection discards
  biology-orthogonal components as intended. Applying the subspace projection
  first and then smoothing would smooth an already-projected (rank-reduced) field,
  which is less geometrically meaningful.

Implementation in `.somalign_anchored_dispatch()`: if `laplacian_lambda > 0`
and a subspace `shift_fn` is also active, compose them explicitly:

```r
shift_fn_sub  <- shift_fn                    # existing subspace projection
shift_fn_lap  <- local({
  L <- .somalign_som_laplacian(query$som_query$grid)
  masses <- query$node_masses
  lam <- laplacian_lambda
  function(s) .somalign_smooth_shifts(s, L, lam, masses)
})
shift_transform <- if (use_subspace && laplacian_lambda > 0) {
  function(s) shift_fn_sub(shift_fn_lap(s))   # smooth -> project
} else if (laplacian_lambda > 0) {
  shift_fn_lap
} else if (use_subspace) {
  shift_fn_sub
} else {
  NULL
}
```

Pass `shift_transform` to `.somalign_finish_fit()` in place of the
current `shift_transform = shift_fn`.

---

## 3. Internal Helpers

Both helpers live in **`R/fit.R`** (keeping all shift post-processing in one
file). They are internal (no `@export`).

### 3.1 `.somalign_som_laplacian(grid)`

```r
# Build the graph Laplacian of a kohonen SOM hexagonal grid.
#
# @param grid  A kohonen::somgrid object. Must have a $pts component
#              (M x 2 matrix of 2-D hexagonal coordinates).
# @return      An M x M numeric matrix L = D - A, where A is the
#              binary adjacency matrix of the unit-distance neighbor graph.
.somalign_som_laplacian <- function(grid) {
  pts <- grid$pts
  if (is.null(pts) || !is.matrix(pts) || ncol(pts) < 2L) {
    stop(
      "SOM grid does not contain 2-D node coordinates ($pts). ",
      "Laplacian smoothing requires a kohonen SOM with grid$pts.",
      call. = FALSE
    )
  }
  M <- nrow(pts)
  # Pairwise squared Euclidean distances between grid coordinates.
  # O(M^2); negligible for M <= 1024.
  dx <- outer(pts[, 1L], pts[, 1L], "-")
  dy <- outer(pts[, 2L], pts[, 2L], "-")
  d2 <- dx^2 + dy^2
  # For hexagonal grids, unit-distance neighbors satisfy d^2 <= 1 + tol.
  # Use a threshold of 1.01^2 to accommodate floating-point hexagonal coords.
  A <- (d2 > 0) & (d2 <= 1.01^2)
  storage.mode(A) <- "double"
  D <- diag(rowSums(A))
  D - A
}
```

**File:** `R/fit.R` (append below `.somalign_node_shifts`, after line 467).

**Complexity:** O(M^2), dominated by the outer products. For M = 256 this is
~65 k operations; for M = 1024 it is ~1 M operations — negligible vs. Sinkhorn.

### 3.2 `.somalign_smooth_shifts(shifts, L, lambda, node_masses)`

```r
# Laplacian-regularised smoothing of node shifts (Tikhonov solve).
#
# Solves the M x M system (W + lambda * L) x = W s  per feature column, where
# W = diag(node_masses).  A single Cholesky factorisation is reused across
# all p right-hand sides.
#
# The correction_allowed attribute is NOT touched here; it is preserved by the
# wrapper in .somalign_finish_fit (lines 143-147) which saves and restores the
# attribute around any shift_transform call.
#
# @param shifts       M x p numeric matrix (raw node shifts).
# @param L            M x M Laplacian matrix from .somalign_som_laplacian().
# @param lambda       Non-negative scalar regularisation strength.
# @param node_masses  Length-M non-negative numeric vector.
# @return             M x p smoothed shift matrix (same dimnames as shifts).
.somalign_smooth_shifts <- function(shifts, L, lambda, node_masses) {
  M <- nrow(shifts)
  W <- diag(node_masses, nrow = M)
  A_sys <- W + lambda * L         # M x M, symmetric positive semi-definite
  # Add a small ridge for robustness when isolated nodes give A_sys singular.
  A_sys <- A_sys + diag(M) * (.Machine$double.eps * max(abs(diag(A_sys))))
  rhs   <- W %*% shifts           # M x p
  # Cholesky: one factorisation, p RHS solves
  ch    <- tryCatch(
    chol(A_sys),
    error = function(e) NULL
  )
  out <- if (!is.null(ch)) {
    chol2inv(ch) %*% rhs
  } else {
    # Fallback: base solve() (slower but always works)
    solve(A_sys, rhs)
  }
  dimnames(out) <- dimnames(shifts)
  out
}
```

**File:** `R/fit.R` (append after `.somalign_som_laplacian`).

---

## 4. Data-Structure Changes

No new fields are stored in the `somalign_fit` object by default. The smoothed
`node_shifts` matrix replaces the raw matrix in-place; downstream code
(`.somalign_project_pair()`, `.somalign_build_nodes_diag()`) sees the smoothed
values transparently.

**Optional diagnostic** (not required for MVP): add a `laplacian_roughness`
entry to `diagnostics$nodes` reporting per-node shift roughness before and after
smoothing. This can be a follow-up addition in `.somalign_build_nodes_diag()`:

```r
# Before smoothing (in .somalign_finish_fit):
roughness_raw <- if (!is.null(L_for_diag))
  sum(diag(t(raw_shifts) %*% L_for_diag %*% raw_shifts))
else NA_real_
```

Defer this diagnostic to a second pass; it requires threading `L` through
`.somalign_finish_fit()`, which changes the signature unnecessarily for the MVP.

---

## 5. Algorithm

### Neighbor graph

For a `kohonen::somgrid` with `topo = "hexagonal"`, `grid$pts` is an M × 2
matrix of hex-grid coordinates. Two nodes are neighbors if their Euclidean
coordinate distance equals 1.0 (unit lattice spacing). The threshold `d^2 <=
1.01^2` absorbs floating-point rounding (actual values are exact multiples of
`0.5` and `sqrt(3)/2 ≈ 0.866`). For rectangular grids (`topo = "rectangular"`),
neighbors are at distance 1 in the same coordinate system.

The binary adjacency matrix `A` (M × M), degree matrix `D = diag(rowSums(A))`,
and Laplacian `L = D - A` are all dense but tiny: at M = 256 they are 512 KB
each (doubles), well within memory.

### Smoothing problem

Let `s_d ∈ R^M` be the raw OT-derived shifts for feature `d`. The regularised
shifts solve:

```
minimise_x  ||W^{1/2} (x - s_d)||^2 + lambda * x^T L x
```

where `W = diag(node_masses)`. Nodes with `correction_allowed = FALSE` have
`node_masses` set to 0 in the construction (see **Section 7**), so they
contribute no data term and are pulled toward the neighbor-average of their
adjacent allowed nodes — a geometrically sensible fallback.

Closed-form solution:

```
(W + lambda * L) x* = W s_d
```

R implementation (from Section 3.2):

```r
A_sys <- diag(node_masses) + lambda * L   # M x M
rhs   <- diag(node_masses) %*% shifts     # M x p (all features at once)
out   <- chol2inv(chol(A_sys)) %*% rhs    # M x p smoothed shifts
```

One Cholesky of the M × M system `A_sys` is computed; all `p` right-hand sides
are solved simultaneously via the `chol2inv` product. For M = 256, p = 40 this
is one Cholesky (~O(M^3/3) ≈ 5.5 M flops) plus one M × p matrix multiply — far
cheaper than the Sinkhorn iterations.

### `shift_transform` signature

The hook in `.somalign_finish_fit()` (lines 143–147) is:

```r
if (!is.null(shift_transform)) {
  allowed <- attr(node_shifts, "correction_allowed")
  node_shifts <- shift_transform(node_shifts)
  attr(node_shifts, "correction_allowed") <- allowed
}
```

The transform receives the M × p matrix and must return an M × p matrix of the
same dimensions. It must **not** set or clear `correction_allowed` — the hook
restores the original attribute after the call. The smoother in
`.somalign_smooth_shifts()` respects this: it does not touch `correction_allowed`.

Crucially, the current code saves `allowed` *before* calling `shift_transform`
and re-attaches it *after*, so the smoother is free to return a matrix without
the attribute.

---

## 6. Integration Points

### `.somalign_finish_fit()` hook (lines 143–147, `R/fit.R`)

No change to `.somalign_finish_fit()` itself. The Laplacian smoother is passed
in as `shift_transform`, the same way the subspace projector is passed by
`somalign_fit_anchored()`.

### `somalign_fit()` call flow

```
somalign_fit()
  -> build shift_transform closure (if laplacian_lambda > 0)
       L <- .somalign_som_laplacian(query$som_query$grid)
       shift_transform <- function(s) .somalign_smooth_shifts(s, L, lambda, masses)
  -> .somalign_align_transport()
  -> .somalign_finish_fit(..., shift_transform = shift_transform)
       -> .somalign_node_shifts()     # raw shifts
       -> shift_transform(node_shifts) # Laplacian smooth
       -> restore correction_allowed
       -> .somalign_project_pair()   # uses smoothed shifts
```

### `somalign_fit_anchored()` — composition with subspace mode

Composition is handled in `.somalign_anchored_dispatch()` (`R/anchored.R`,
lines 222–225). The current code builds `shift_fn` there and passes it to
`.somalign_finish_fit()`. The `laplacian_lambda` argument must be threaded from
`somalign_fit_anchored()` → `.somalign_anchored_dispatch()`.

**Order: smooth first, then project into subspace V.** Rationale: the Laplacian
operates in full p-dimensional marker space; the subspace projection then
discards biology-orthogonal components. Reversing the order would smooth an
already-projected field, which is geometrically inconsistent (the Laplacian
neighborhood structure is defined by marker-space proximity, not subspace
proximity).

```r
# In .somalign_anchored_dispatch():
shift_fn_lap <- if (laplacian_lambda > 0) {
  L <- .somalign_som_laplacian(query$som_query$grid)
  masses <- query$node_masses
  lam <- laplacian_lambda
  function(s) .somalign_smooth_shifts(s, L, lam, masses)
} else NULL

shift_transform <- if (!is.null(shift_fn) && !is.null(shift_fn_lap)) {
  function(s) shift_fn(shift_fn_lap(s))   # lap smooth -> subspace project
} else if (!is.null(shift_fn_lap)) {
  shift_fn_lap
} else {
  shift_fn    # may be NULL if neither subspace nor laplacian
}
```

---

## 7. Edge Cases

### Non-grid SOMs / missing grid coordinates

`query$som_query$grid$pts` may be `NULL` if the user supplied a custom SOM
object via `somalign_query_from_som()` with an object that lacks `$grid`. Guard:

```r
if (is.null(query$som_query$grid) || is.null(query$som_query$grid$pts)) {
  stop(
    "`laplacian_lambda > 0` requires the query SOM to have a kohonen ",
    "grid with 2-D node coordinates (grid$pts). ",
    "Use a SOM trained with kohonen::som() or kohonen::supersom().",
    call. = FALSE
  )
}
```

Emit this check in `somalign_fit()` immediately after the `laplacian_lambda`
validation, before `.somalign_align_transport()`.

### `laplacian_lambda = 0` short-circuit

When `laplacian_lambda == 0`, `shift_transform` is `NULL`; the hook at lines
143–147 is skipped entirely. No Laplacian is constructed. Zero overhead.

### Isolated nodes

A node with no grid neighbors (degree 0 in the neighbor graph) has a 0 row/col
in `L`. The system `(W + lambda * L) x = W s` reduces to `w_i x_i = w_i s_i`
for that node, so the smoothed shift equals the raw shift regardless of `lambda`.
This is the correct behavior: an isolated node receives no spatial regularization.
The small ridge added in `.somalign_smooth_shifts()` (`eps * max(diag(A_sys))`)
prevents singularity from zero-mass nodes with degree 0.

### Interaction with `correction_allowed`

Nodes with `correction_allowed = FALSE` have raw shifts already set to 0 by
`.somalign_node_shifts()` (line 462–464). In the smoother, their `node_masses`
entries are used as-is (they may be small but non-zero). This means a
`correction_allowed = FALSE` node that has low but positive mass will pull
slightly toward the neighbor average.

**Alternative (more conservative):** zero out `node_masses` for
`!correction_allowed` nodes inside the smoother, making them purely passive
(receive neighbor average, contribute nothing to the system RHS):

```r
.somalign_smooth_shifts <- function(shifts, L, lambda, node_masses) {
  # correction_allowed is NOT passed in; the hook restores it.
  # We zero mass for nodes whose raw shift is identically zero in all features,
  # as a proxy for !correction_allowed.
  # Better: pass correction_allowed explicitly.
  ...
}
```

**Recommended implementation:** pass `correction_allowed` explicitly as a fourth
argument:

```r
.somalign_smooth_shifts <- function(shifts, L, lambda, node_masses,
                                    correction_allowed = NULL) {
  w <- node_masses
  if (!is.null(correction_allowed)) w[!correction_allowed] <- 0
  W <- diag(w, nrow = length(w))
  ...
}
```

Then the closure in `somalign_fit()` is:

```r
shift_transform <- function(s) {
  ca <- attr(s, "correction_allowed")
  .somalign_smooth_shifts(s, L, laplacian_lambda, masses, correction_allowed = ca)
}
```

This is clean: the attribute is read before the hook strips it, the smoother
uses it to zero out disallowed-node mass terms, and the hook restores it after.
Disallowed nodes receive the spatial interpolation from their allowed neighbors
but do not anchor the fit.

---

## 8. Tests

Add `tests/testthat/test-laplacian-shifts.R`.

### Test 1 — `lambda = 0` identity

```r
test_that("laplacian_lambda = 0 produces identical result to default", {
  skip_if_not_installed("kohonen")
  set.seed(42)
  mat <- matrix(rnorm(60), nrow = 30, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  qry <- somalign_query(
    mat + 0.5, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  fit0    <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 0)
  expect_equal(fit0$node_shifts, fit_lap$node_shifts)
  expect_equal(
    attr(fit0$node_shifts, "correction_allowed"),
    attr(fit_lap$node_shifts, "correction_allowed")
  )
})
```

### Test 2 — `correction_allowed` attribute is preserved

```r
test_that("correction_allowed attribute is preserved after Laplacian smoothing", {
  skip_if_not_installed("kohonen")
  set.seed(7)
  mat <- matrix(rnorm(80), nrow = 40, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  qry <- somalign_query(
    mat + 1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  fit0    <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 0.5)
  # correction_allowed must be identical (not changed by smoothing)
  expect_identical(
    attr(fit0$node_shifts, "correction_allowed"),
    attr(fit_lap$node_shifts, "correction_allowed")
  )
})
```

### Test 3 — Large `lambda` drives shifts toward global mean

```r
test_that("large laplacian_lambda collapses shifts toward their mass-weighted mean", {
  skip_if_not_installed("kohonen")
  set.seed(13)
  mat <- matrix(rnorm(80), nrow = 40, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  qry <- somalign_query(
    mat + 1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"), rlen = 10
  )
  fit0    <- somalign_fit(qry, ref)
  fit_lap <- somalign_fit(qry, ref, laplacian_lambda = 1e6)
  # Smoothed shifts should have much lower variance across nodes
  raw_var   <- var(fit0$node_shifts[, 1])
  smooth_var <- var(fit_lap$node_shifts[, 1])
  expect_lt(smooth_var, raw_var * 0.1)   # at least 10x variance reduction
})
```

### Test 4 — `laplacian_lambda` validation

```r
test_that("invalid laplacian_lambda raises informative error", {
  skip_if_not_installed("kohonen")
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                dimnames = list(NULL, c("F1", "F2")))
  ref <- somalign_train_reference(
    mat, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5
  )
  qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
  expect_error(somalign_fit(qry, ref, laplacian_lambda = -1), "non-negative")
  expect_error(somalign_fit(qry, ref, laplacian_lambda = "x"), "non-negative")
})
```

### Test 5 — Composition with subspace mode in `somalign_fit_anchored()`

```r
test_that("laplacian_lambda + subspace mode produces correction_allowed-preserving result", {
  skip_if_not_installed("kohonen")
  fx <- make_subspace_fixture()
  fit <- somalign_fit_anchored(
    fx$qry, fx$ref,
    anchor_old = fx$anc_old, anchor_new = fx$anc_new,
    rho_anchor = 1, correction = "subspace",
    laplacian_lambda = 0.2
  )
  expect_true(!is.null(attr(fit$node_shifts, "correction_allowed")))
  # Shifts restricted to batch subspace V (verify norm in orthogonal direction is near 0)
  V <- fit$anchors$batch_subspace$V
  shifts_orth <- fit$node_shifts - fit$node_shifts %*% V %*% t(V)
  expect_lt(sqrt(mean(shifts_orth^2)), 1e-8)
})
```

---

## 9. Docs / NAMESPACE

### Roxygen

Add `@param laplacian_lambda` to the `somalign_fit()` roxygen block (before
`@details`):

```
#' @param laplacian_lambda Non-negative scalar. Graph-Laplacian regularisation
#'   strength for the node-shift field. When greater than zero, the M × p
#'   raw node shifts are smoothed by solving \eqn{(W + \lambda L)\,s^* = W\,s},
#'   where \eqn{W = \mathrm{diag}(\text{node\_masses})} and \eqn{L} is the
#'   graph Laplacian of the SOM hexagonal neighbor graph. This penalises squared
#'   differences between adjacent-node shifts, producing a spatially coherent
#'   correction field. Default \code{0} (no smoothing, exact current behaviour).
#'   A natural starting range is \code{0.1}--\code{1.0}; larger values
#'   increasingly collapse the field toward its mass-weighted mean.
```

Same parameter block for `somalign_fit_anchored()` if added.

### NAMESPACE

No new exports. Both internal helpers (`.somalign_som_laplacian`,
`.somalign_smooth_shifts`) are unexported. NAMESPACE is unchanged.

---

## 10. BiocCheck

`somalign_fit()` currently has a body well under 50 lines (lines 76–118 = 43
lines including blank lines and comments). Adding 5–8 lines for the
`laplacian_lambda` validation and closure construction keeps it under the 50-line
exported-body limit.

If needed, extract the closure-building step:

```r
.somalign_make_laplacian_transform <- function(query, laplacian_lambda) {
  if (laplacian_lambda == 0) return(NULL)
  L <- .somalign_som_laplacian(query$som_query$grid)
  masses <- query$node_masses
  function(s) {
    ca <- attr(s, "correction_allowed")
    .somalign_smooth_shifts(s, L, laplacian_lambda, masses, correction_allowed = ca)
  }
}
```

This keeps the exported function at ~2 extra lines.

---

## 11. Effort, Risks, Dependencies

### Effort: Medium (~3–4 hours implementation + 1–2 hours tests)

- ~50 lines of new R code (two helpers + closure logic in `somalign_fit`).
- ~80 lines of new tests.
- Minor roxygen additions.

### Risks

1. **SOM objects without `$grid$pts`**: users who construct a `somalign_query`
   via `somalign_query_from_som()` with a synthetic SOM lacking a proper kohonen
   grid will get a clear error at `laplacian_lambda > 0`. No silent wrong result.

2. **`lambda` selection**: too large collapses all corrections to a single global
   vector (topology-preserving but under-correcting). Document the analogy to
   `epsilon` and recommend `somalign_sensitivity_grid()` extension for
   lambda-sweeps as future work.

3. **Interaction with `correction_allowed`**: the recommended implementation
   (pass `correction_allowed` into `.somalign_smooth_shifts()` and zero mass for
   disallowed nodes) is more conservative than using raw `node_masses`, but
   means disallowed nodes receive interpolated shifts. If this is undesirable,
   zero out the smoothed shifts for disallowed nodes after the solve — but the
   current hook already re-enforces `correction_allowed` as an attribute, not as
   a zero-mask on the shift values. Review `.somalign_project_pair()` to confirm
   whether disallowed nodes' shift values are ever applied to cells: line 258
   applies `node_shifts[query$sample_unit, , drop = FALSE]` to all cells, so
   disallowed-node shifts **are** applied. This means zeroing disallowed nodes'
   smoothed shifts (to preserve the existing `correction_min_mass` semantics) is
   important. Add a post-solve zero-out step:

   ```r
   if (!is.null(correction_allowed)) out[!correction_allowed, ] <- 0
   ```

4. **Shares `shift_transform` hook with Idea 7** (if another idea also uses
   this hook). See composition logic in Section 6. Both transforms can be
   composed as nested function calls — no architectural conflict, but the
   composition order must be documented.

5. **Dense M × M Laplacian**: for M > 2000 the dense matrix becomes > 32 MB.
   For current typical use (M = 64–400) this is not a problem. Add a warning for
   M > 1500 suggesting use of the `Matrix` package sparse Laplacian as future
   optimization.

### Dependencies

No new R package dependencies. `chol()`, `chol2inv()`, `diag()`, and `solve()`
are all base R. The `Matrix` package (already a Bioconductor suggestion) can be
used for sparse Laplacian at larger M, but is not required for the MVP.
