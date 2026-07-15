# Implementation Plan: Simulated-Annealing Sinkhorn Solver (`solver = "annealing"`)

**Date:** 2026-07-16  
**Idea source:** `analysis/ideas/2026-07-16-methods-improvement/01_statistical-physicist.md`, Idea 2  
**Effort estimate:** Low-Medium (~2 days)

---

## 1. Summary

Add `solver = "annealing"` to `somalign_fit`, `somalign_fit_anchored`,
`somalign_fit_two_pass`, and `somalign_sensitivity_grid`. When selected, the
solver cools from a hot start (`anneal_start * epsilon`) down to the target
`epsilon` over `anneal_stages` geometric steps, warm-starting each stage from
the previous stage's converged log-domain dual potentials `f` and `g`. This
escapes local minima in rugged cost landscapes (especially `label_guided = TRUE`
fits with the `max(cost) * 1e4` penalty) and converges faster on hard
(small-epsilon) problems than isothermal `log_domain`. Backward compatibility is
preserved: existing `"internal"`, `"log_domain"`, and `"auto"` choices are
unchanged.

---

## 2. Public API

### 2.1 New `solver` choice

All four exported fit functions gain `"annealing"` as an additional valid value
for their `solver` argument. The `match.arg` call at the top of each function
body is the only edit needed; call sites below `.somalign_solve_ot` need no
changes because the dispatch lives entirely inside that function and the new
internal helper.

#### `somalign_fit` — `R/fit.R` line 81 / 97

```r
# Before (line 81):
solver = c("internal", "log_domain", "auto"),
# After:
solver = c("internal", "log_domain", "auto", "annealing"),

# Before (line 97):
solver <- match.arg(solver, c("internal", "log_domain", "auto"))
# After:
solver <- match.arg(solver, c("internal", "log_domain", "auto", "annealing"))
```

New tuning parameters added to `somalign_fit`'s signature after `tol`:

```r
anneal_start  = 10,    # eps_0 = anneal_start * epsilon  (multiplier >= 1)
anneal_factor = NULL,  # NULL → computed from anneal_start and anneal_stages
anneal_stages = 10L    # number of cooling stages (integer >= 1)
```

These are only used when `solver = "annealing"` and are silently ignored
otherwise, so existing code that does not pass them is unaffected.

#### `somalign_fit_anchored` — `R/anchored.R` line 161 / 172

Same `match.arg` edit; same three new tuning args inserted after `tol`.

#### `somalign_fit_two_pass` — `R/fit.R` line 561 / 578

Same `match.arg` edit; same three new tuning args inserted after `tol`.
Note: `somalign_fit_two_pass` calls `.somalign_align_transport` twice (lines
590–591 and 615–616). Both calls receive the `solver` argument unchanged — the
annealing driver will run independently for each pass with the pass-specific
epsilon (`epsilon_global`, `epsilon_local`). The warm-start is *within* each
pass, not across passes.

#### `somalign_sensitivity_grid` — `R/diagnostics.R` line 72 / 91

```r
# Before (line 72):
solver = c("internal", "auto"),
# After:
solver = c("internal", "auto", "log_domain", "annealing"),

# Before (line 91):
solver <- match.arg(solver)
# After:
solver <- match.arg(solver, c("internal", "auto", "log_domain", "annealing"))
```

Same three new tuning args; pass-through to `.somalign_align_transport` via the
`.run_one` closure (line 104).

### 2.2 Tuning arg semantics and defaults

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `anneal_start` | positive scalar | `10` | Hot epsilon = `anneal_start * epsilon`. 10× means starting at ~10 times the target temperature. |
| `anneal_stages` | positive integer | `10L` | Number of discrete cooling steps including the final step at `epsilon`. |
| `anneal_factor` | positive scalar or `NULL` | `NULL` | Explicit per-stage ratio `eps_{k+1}/eps_k`. If `NULL`, computed as `(1/anneal_start)^(1/(anneal_stages-1))` — exactly the geometric schedule. Provided for power users who want to override the auto-schedule. |

Validation: `anneal_start >= 1`, `anneal_stages >= 1L`. If `anneal_stages == 1`,
the schedule degenerates to a cold-start `log_domain` solve at `epsilon` (see
§7 edge cases).

### 2.3 Passing tuning args through the call stack

`.somalign_align_transport` gains three new `...`-style arguments (or explicit
formals). The cleanest approach is explicit formals because the function is
internal and all call sites are known:

```r
# R/fit.R — .somalign_align_transport signature (line 205):
.somalign_align_transport <- function(query, reference, epsilon, rho_query,
                                      rho_ref, solver, max_iter, tol,
                                      cost_bonus = NULL,
                                      diagonal_boost = 0,
                                      label_mask = NULL,
                                      anneal_start = 10,
                                      anneal_factor = NULL,
                                      anneal_stages = 10L) {
  ...
  ot <- .somalign_solve_ot(
    cost = cost_normalized, a = query$node_masses, b = reference$node_masses,
    epsilon = epsilon, rho_query = rho_query, rho_ref = rho_ref,
    solver = solver, max_iter = max_iter, tol = tol,
    anneal_start = anneal_start, anneal_factor = anneal_factor,
    anneal_stages = anneal_stages
  )
  ...
}
```

`.somalign_solve_ot` gains matching formals with the same defaults.

---

## 3. Internal Helpers

All new code lives in `R/ot.R`. Three changes:

### 3.1 Warm-start interface for `.somalign_solve_internal_log`

**File:** `R/ot.R`, current lines 96–152.

Add two optional arguments `f_init` and `g_init` (default `NULL`):

```r
.somalign_solve_internal_log <- function(cost,
                                         a,
                                         b,
                                         epsilon,
                                         rho_query,
                                         rho_ref,
                                         max_iter,
                                         tol,
                                         f_init = NULL,  # NEW
                                         g_init = NULL)  # NEW
{
  tau_a <- rho_query / (rho_query + epsilon)   # recomputed from epsilon — critical
  tau_b <- rho_ref   / (rho_ref   + epsilon)   # idem

  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  # Warm-start: use supplied potentials if provided, else zero-initialise
  f <- if (!is.null(f_init)) f_init else numeric(length(a))   # NEW
  g <- if (!is.null(g_init)) g_init else numeric(length(b))   # NEW

  cost_over_eps <- cost / epsilon
  ...  # remainder of function unchanged
}
```

The change is exactly two lines of initialisation (replacing lines 110–111 in
current `ot.R`) plus two new formal arguments. The convergence loop, delta
calculation, and plan reconstruction are untouched.

**Critical correctness point:** `tau_a` and `tau_b` are computed at lines
104–105 of the current function from the *epsilon passed in at that call*, not
stored from a previous call. This is already correct. When the annealing driver
calls `.somalign_solve_internal_log` at each stage with a different `epsilon`,
the function recomputes `tau_a = rho_query / (rho_query + eps_k)` and
`tau_b = rho_ref / (rho_ref + eps_k)` automatically. The warm-start potentials
`f_init`/`g_init` come from the previous (hotter) epsilon stage, so they are
not at equilibrium for the new epsilon — this is intentional and is the entire
point of warm starting. No special handling is needed.

### 3.2 Cooling-schedule helper `.somalign_cooling_schedule`

**File:** `R/ot.R` (new function, ~15 lines).

```r
.somalign_cooling_schedule <- function(epsilon_target, anneal_start,
                                       anneal_stages, anneal_factor) {
  if (anneal_stages == 1L) return(epsilon_target)
  eps_0 <- anneal_start * epsilon_target
  if (is.null(anneal_factor)) {
    # Geometric: eps_k = eps_0 * r^k, r chosen so eps_{n-1} = epsilon_target
    r <- (epsilon_target / eps_0) ^ (1 / (anneal_stages - 1))
  } else {
    r <- anneal_factor
  }
  # k = 0, 1, ..., anneal_stages-1
  schedule <- eps_0 * r ^ seq(0, anneal_stages - 1)
  # Clamp final value to exact epsilon_target to avoid floating-point drift
  schedule[anneal_stages] <- epsilon_target
  schedule
}
```

This returns a length-`anneal_stages` numeric vector, decreasing from
`anneal_start * epsilon` to `epsilon`.

### 3.3 Annealing driver `.somalign_solve_annealing`

**File:** `R/ot.R` (new function, ~45 lines).

```r
.somalign_solve_annealing <- function(cost, a, b, epsilon,
                                       rho_query, rho_ref,
                                       max_iter, tol,
                                       anneal_start, anneal_factor,
                                       anneal_stages) {
  schedule <- .somalign_cooling_schedule(epsilon, anneal_start,
                                         anneal_stages, anneal_factor)
  n_stages  <- length(schedule)
  # Per-stage iteration budget: interior stages use a looser budget;
  # the final stage runs to full tolerance.
  inner_iter   <- max(ceiling(max_iter / n_stages), 20L)
  inner_tol    <- tol * 10   # looser inner tolerance for interior stages
  total_iter   <- 0L
  stage_info   <- vector("list", n_stages)

  f <- NULL   # warm-start potentials; NULL triggers zero-init in first stage
  g <- NULL

  for (k in seq_len(n_stages)) {
    eps_k    <- schedule[k]
    is_final <- (k == n_stages)
    iter_k   <- if (is_final) max_iter else inner_iter
    tol_k    <- if (is_final) tol       else inner_tol

    result <- .somalign_solve_internal_log(
      cost     = cost,
      a        = a,
      b        = b,
      epsilon  = eps_k,
      rho_query = rho_query,
      rho_ref   = rho_ref,
      max_iter = iter_k,
      tol      = tol_k,
      f_init   = f,
      g_init   = g
    )
    f <- result$f   # carry potentials forward (see §3.1 — need to expose f, g)
    g <- result$g
    total_iter <- total_iter + result$iterations
    stage_info[[k]] <- list(epsilon = eps_k, iterations = result$iterations,
                            converged = result$converged,
                            final_delta = result$final_delta)
  }

  # The final result is at eps_k = epsilon; its plan and convergence are definitive
  list(
    plan          = result$plan,
    iterations    = total_iter,
    converged     = result$converged,
    final_delta   = result$final_delta,
    stage_info    = stage_info,   # per-stage diagnostics list
    schedule      = schedule      # full epsilon schedule used
  )
}
```

**Dependency:** the driver needs `f` and `g` from `.somalign_solve_internal_log`'s
return value. Currently that function only returns `plan`, `iterations`,
`converged`, `final_delta` (lines 147–151 of `ot.R`). We must add `f` and `g`
to the return list:

```r
# At the end of .somalign_solve_internal_log — replace current return:
list(plan = plan, iterations = iterations, converged = converged,
     final_delta = final_delta, f = f, g = g)   # f and g are NEW
```

This is a backward-compatible addition; existing callers ignore extra list
elements.

### 3.4 Dispatch in `.somalign_solve_ot`

**File:** `R/ot.R`, lines 1–33.

```r
.somalign_solve_ot <- function(cost, a, b, epsilon, rho_query, rho_ref,
                               solver, max_iter, tol,
                               anneal_start = 10,       # NEW
                               anneal_factor = NULL,    # NEW
                               anneal_stages = 10L) {   # NEW
  requested_solver <- solver
  notes <- character()
  .somalign_validate_ot_inputs(cost, a, b, epsilon, rho_query, rho_ref)

  if (solver == "auto") {
    solver <- "internal"
    notes <- c(notes, "`solver = \"auto\"` uses the internal generalized Sinkhorn solver.")
  }

  if (solver == "annealing") {                           # NEW BRANCH
    internal <- .somalign_solve_annealing(
      cost, a, b, epsilon, rho_query, rho_ref,
      max_iter, tol, anneal_start, anneal_factor, anneal_stages
    )
    return(list(
      plan             = internal$plan,
      solver           = "annealing",
      requested_solver = requested_solver,
      notes            = notes,
      iterations       = internal$iterations,
      converged        = internal$converged,
      final_delta      = internal$final_delta,
      stage_info       = internal$stage_info,   # NEW diagnostic field
      schedule         = internal$schedule      # NEW diagnostic field
    ))
  }

  log_domain <- identical(solver, "log_domain")
  internal <- .somalign_solve_internal(cost, a, b, epsilon, rho_query, rho_ref,
                                       max_iter, tol, log_domain = log_domain)
  list(
    plan = internal$plan,
    solver = solver,
    requested_solver = requested_solver,
    notes = notes,
    iterations = internal$iterations,
    converged = internal$converged,
    final_delta = internal$final_delta
  )
}
```

The existing `log_domain` / `internal` path is completely unchanged.

---

## 4. Data-Structure Changes

### 4.1 `diagnostics$solver` additions

`.somalign_build_diagnostics` (currently `R/fit.R` lines 264–307) reads from
`ot` (the return of `.somalign_solve_ot`). Two new optional fields are plumbed
through:

```r
# In .somalign_build_diagnostics, inside the solver = list(...) block:
solver = list(
  ...                                              # all existing fields
  anneal_schedule  = ot$schedule,    # NULL for non-annealing solvers
  anneal_stage_info = ot$stage_info  # NULL for non-annealing solvers
)
```

`ot$schedule` and `ot$stage_info` are `NULL` when the non-annealing path was
used, so `somalign_diagnostics` callers see no change unless they explicitly
inspect these fields.

### 4.2 Per-stage info structure

`stage_info` is a list of length `anneal_stages`, each element a named list:

```r
list(
  epsilon     = <numeric scalar>,   # temperature for this stage
  iterations  = <integer>,          # Sinkhorn iterations consumed
  converged   = <logical>,          # whether inner tol was met
  final_delta = <numeric>           # final relative change in potentials
)
```

This lets users diagnose which stage consumed most iterations or failed to
converge at its inner tolerance.

---

## 5. Algorithm Details

### 5.1 Geometric cooling schedule

Given target `epsilon`, multiplier `anneal_start` (default 10), and stages
`n = anneal_stages` (default 10):

```
eps_0 = anneal_start * epsilon     # e.g. 1.0 when epsilon = 0.1
r     = (epsilon / eps_0)^(1/(n-1))  # = (1/10)^(1/9) ≈ 0.774
eps_k = eps_0 * r^k,  k = 0,...,n-1
eps_{n-1} = epsilon (exact, clamped)
```

For `epsilon = 0.1`, `anneal_start = 10`, `anneal_stages = 10`:
schedule = {1.0, 0.774, 0.599, 0.464, 0.359, 0.278, 0.215, 0.167, 0.129, 0.1}

### 5.2 Warm-start handoff

At stage `k`, the log-domain solver converges to `f^(k)` and `g^(k)` at
temperature `eps_k`. These are passed directly as `f_init`/`g_init` to stage
`k+1`. Because Sinkhorn fixed-point equations are `C^0` in epsilon (the
potentials shift continuously as temperature changes), the warm start is
near-equilibrium, reducing the iterations needed at `eps_{k+1}` substantially
compared to zero initialisation.

The potentials at temperature `eps_k` satisfy:
```
f_i = tau_a(eps_k) * (eps_k * log a_i - eps_k * lse_j((g_j - C_ij)/eps_k))
```
When epsilon changes to `eps_{k+1}`, `tau_a`, `tau_b`, and `cost_over_eps`
all change. The warm potentials are *not* at equilibrium for the new epsilon
but they are far better than zeros because the optimal coupling structure
(which query node maps to which reference node) changes slowly along the
geometric schedule.

### 5.3 Convergence per stage vs. final

- **Interior stages** (k = 1, ..., n-2): run for at most `ceil(max_iter /
  n_stages)` iterations, convergence criterion `10 * tol`. The looser tolerance
  avoids wasting iterations equilibrating fully at an intermediate temperature;
  partial convergence is sufficient for a good warm start.
- **Final stage** (k = n-1, epsilon = target): runs for at most `max_iter`
  iterations with full `tol`. This stage issues the standard
  `.somalign_warn_convergence` warning if `max_iter` is exhausted.

Interior stage non-convergence is recorded in `stage_info[[k]]$converged` but
does **not** emit a warning — this is expected and normal.

### 5.4 Pseudocode of the driver loop

```r
schedule <- .somalign_cooling_schedule(epsilon, anneal_start, anneal_stages, anneal_factor)
f <- NULL; g <- NULL; total_iter <- 0L
for k in seq_along(schedule):
    eps_k    <- schedule[k]
    is_final <- (k == length(schedule))
    result   <- .somalign_solve_internal_log(
                    cost, a, b, eps_k, rho_query, rho_ref,
                    max_iter = if (is_final) max_iter else ceil(max_iter / n_stages),
                    tol      = if (is_final) tol      else tol * 10,
                    f_init   = f, g_init = g)
    f <- result$f; g <- result$g
    total_iter <- total_iter + result$iterations
plan <- result$plan    # final-stage plan at target epsilon
```

---

## 6. Integration Points

### 6.1 Call chain summary

```
somalign_fit (fit.R:76)
  └─ .somalign_align_transport (fit.R:205)   ← gains anneal_* args
       └─ .somalign_solve_ot (ot.R:1)        ← gains anneal_* args; new "annealing" branch
            └─ .somalign_solve_annealing (ot.R: new)
                 └─ .somalign_solve_internal_log (ot.R:96)  ← gains f_init/g_init
```

`somalign_fit_anchored` → `.somalign_anchored_dispatch` → `.somalign_align_transport`
(same chain; `anneal_*` args forwarded through `somalign_fit_anchored` →
`.somalign_anchored_dispatch` signature).

`somalign_fit_two_pass` → two calls to `.somalign_align_transport`; both
receive the same `solver` and `anneal_*` args.

`somalign_sensitivity_grid` → `.run_one` closure → `.somalign_align_transport`;
same chain.

### 6.2 `.somalign_anchored_dispatch` signature

**File:** `R/anchored.R` line 191. Add `anneal_start`, `anneal_factor`,
`anneal_stages` to its formals and the single `.somalign_align_transport` call
inside it.

### 6.3 Existing solver choices unchanged

The `solver == "annealing"` dispatch is an early-return branch that runs before
the existing `log_domain <- identical(solver, "log_domain")` line (current
`ot.R` line 18). The primal `"internal"` and log-domain paths are untouched.

---

## 7. Edge Cases and Numerical Stability

### 7.1 `tau_a` / `tau_b` recomputation — critical

`tau_a = rho_query / (rho_query + epsilon)` and `tau_b = rho_ref / (rho_ref +
epsilon)` appear at lines 104–105 of `.somalign_solve_internal_log` and are
recomputed from the `epsilon` argument on every call. Because the annealing
driver calls the function with a different `epsilon` at each stage, `tau_a` and
`tau_b` are automatically correct at each temperature. **No special handling is
needed** — but this is the most common mistake in naive implementations that
cache `tau_a`/`tau_b` outside the loop. The plan author must confirm that any
future refactor of `.somalign_solve_internal_log` that hoists `tau_a`/`tau_b`
to an outer scope would break the annealing driver.

### 7.2 Schedule overshoot

If `anneal_factor` is supplied by the user and is > 1 (a cooling-rate > 1 is a
warming schedule), the solver would diverge. Validate:

```r
# In .somalign_cooling_schedule:
if (!is.null(anneal_factor) && anneal_factor >= 1)
  stop("`anneal_factor` must be < 1 (a per-stage cooling ratio).", call. = FALSE)
```

Additionally, if `anneal_start < 1` the hot epsilon would be below the target,
meaning no annealing occurs. Validate:

```r
if (!is.numeric(anneal_start) || anneal_start < 1)
  stop("`anneal_start` must be >= 1.", call. = FALSE)
```

Both checks belong in the `somalign_fit` validation block alongside
`.somalign_check_fit_params`, only when `solver == "annealing"`.

### 7.3 Single-stage degenerate case (`anneal_stages = 1L`)

`.somalign_cooling_schedule` returns `c(epsilon_target)` — a length-1 vector.
The driver loop runs once with `is_final = TRUE`, full `max_iter` and `tol`,
and `f_init = NULL` (zero initialisation). This is identical to a plain
`log_domain` solve. The `stage_info` list has one element. No special-casing
needed; the code handles it naturally.

### 7.4 Interaction with the F3 all-underflow warning

The F3 warning is emitted by `.somalign_sinkhorn_kernel` in the *primal* solver
path (`ot.R:154–187`). The annealing solver always uses the log-domain path and
never calls `.somalign_sinkhorn_kernel`. Consequently, F3 cannot fire during an
annealing solve — which is one of the reasons to prefer `"annealing"` over
`"internal"` for small-epsilon problems. No code change needed; document this in
the `@param solver` roxygen.

### 7.5 `cost_over_eps` recomputation

Inside `.somalign_solve_internal_log`, `cost_over_eps <- cost / epsilon` is
computed at line 112 (current). Because the function is called anew at each
stage with a fresh `epsilon`, this recomputation is automatic and correct.

### 7.6 NaN potentials after warm start

When `f_init` values from the previous stage are very large (which can happen
if `rho_query` is small and the plan is nearly deterministic), the initial
`lse_g` computation at the new (cooler) temperature can produce `-Inf`. The
existing `f[is.nan(f)] <- 0` guard (line 125 of `ot.R`) handles this. No
additional guard is required.

### 7.7 `anneal_start = 1` (no cooling)

Produces a single-step schedule equal to `epsilon`. Equivalent to a standard
log-domain cold-start. Valid; no special treatment.

---

## 8. Tests

All tests belong in a new file `tests/testthat/test-annealing-solver.R`.
Fixtures use the existing `tiny_reference()` and `make_som()` helpers from
`tests/testthat/helper-fixtures.R`.

### 8.1 Easy problem: annealing ≈ log_domain

```r
test_that("annealing solver gives approximately the same plan as log_domain on easy problem", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit_log  <- somalign_fit(qry, ref, solver = "log_domain",  epsilon = 0.2)
  fit_ann  <- somalign_fit(qry, ref, solver = "annealing",   epsilon = 0.2,
                           anneal_start = 5, anneal_stages = 5L)

  # Transport plans should be close (not identical due to different convergence paths)
  expect_equal(fit_ann$transport_plan, fit_log$transport_plan, tolerance = 1e-4)
})
```

### 8.2 Diagnostics structure

```r
test_that("annealing diagnostics contain schedule and stage_info", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit <- somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1,
                      anneal_start = 10, anneal_stages = 4L)
  diag <- somalign_diagnostics(fit)

  expect_equal(length(diag$solver$anneal_schedule), 4L)
  expect_equal(diag$solver$anneal_schedule[4], 0.1)
  expect_true(diag$solver$anneal_schedule[1] > diag$solver$anneal_schedule[4])
  expect_equal(length(diag$solver$anneal_stage_info), 4L)
})
```

### 8.3 Hard problem: cold-start warns, annealing does not (label-guided)

```r
test_that("annealing converges where cold-start log_domain warns on label-guided hard problem", {
  skip_if_not_installed("kohonen")
  withr::local_seed(99)
  # Build a well-separated two-cluster problem with label probabilities.
  # The label-guided penalty (max_cost * 1e4) makes the cost landscape rugged.
  mat <- rbind(matrix(rnorm(30 * 4, mean = -3), ncol = 4),
               matrix(rnorm(30 * 4, mean =  3), ncol = 4))
  colnames(mat) <- paste0("f", 1:4)
  ref <- somalign_train_reference(mat, grid = kohonen::somgrid(3, 3, "hexagonal"),
                                  rlen = 20)
  # Assign deterministic label_prob to reference (two classes)
  n_ref <- nrow(ref$codebook)
  lp_ref <- matrix(0, nrow = n_ref, ncol = 2,
                   dimnames = list(NULL, c("A", "B")))
  lp_ref[seq_len(ceiling(n_ref / 2)), "A"] <- 1
  lp_ref[seq(ceiling(n_ref / 2) + 1, n_ref), "B"] <- 1
  ref$label_prob <- lp_ref

  qry <- somalign_query(mat + 0.1, ref, grid = kohonen::somgrid(3, 3, "hexagonal"),
                        rlen = 20)
  lp_qry <- lp_ref   # same structure for query
  qry$label_prob <- lp_qry

  # Cold-start log_domain at very small epsilon may warn or not converge well
  expect_warning(
    somalign_fit(qry, ref, solver = "log_domain", epsilon = 0.02,
                 label_guided = TRUE, max_iter = 50L, tol = 1e-7),
    regexp = "did not converge|delta"
  )

  # Annealing should reach full convergence without warning
  expect_no_warning(
    somalign_fit(qry, ref, solver = "annealing", epsilon = 0.02,
                 label_guided = TRUE, max_iter = 500L, tol = 1e-7,
                 anneal_start = 100, anneal_stages = 10L)
  )
})
```

### 8.4 `anneal_stages = 1` behaves like `log_domain`

```r
test_that("anneal_stages = 1 produces same result as log_domain cold start", {
  ref <- tiny_reference()
  query <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  colnames(query) <- ref$features
  qry <- somalign_query(query, ref, som_query = make_som(rbind(c(-1, 0), c(1, 0))))

  fit_log <- somalign_fit(qry, ref, solver = "log_domain", epsilon = 0.1)
  fit_ann <- somalign_fit(qry, ref, solver = "annealing",  epsilon = 0.1,
                          anneal_stages = 1L)
  expect_equal(fit_ann$transport_plan, fit_log$transport_plan, tolerance = 1e-10)
})
```

### 8.5 Invalid `anneal_start` is caught

```r
test_that("anneal_start < 1 is rejected", {
  ref <- tiny_reference()
  qry <- somalign_query(ref$codebook, ref, som_query = ref$som_ref,
                        codebook_space = "reference_scaled")
  expect_error(
    somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1, anneal_start = 0.5),
    "`anneal_start` must be >= 1"
  )
})
```

### 8.6 `diagnostics$solver$used` is "annealing"

```r
test_that("diagnostics$solver$used is 'annealing' when annealing solver is requested", {
  ref <- tiny_reference()
  qry <- somalign_query(ref$codebook, ref, som_query = ref$som_ref,
                        codebook_space = "reference_scaled")
  fit <- somalign_fit(qry, ref, solver = "annealing", epsilon = 0.1)
  expect_equal(somalign_diagnostics(fit)$solver$used, "annealing")
})
```

---

## 9. Docs and NAMESPACE

### 9.1 No new exports

All new functions (`.somalign_solve_annealing`, `.somalign_cooling_schedule`)
are internal (dot-prefixed). No additions to `NAMESPACE` are needed.

### 9.2 roxygen `@param solver` update

Update in `R/fit.R` (the shared `@param solver` block at line 18, and the copy
in `somalign_fit_two_pass` at line 499), `R/anchored.R`, and `R/diagnostics.R`:

```
#' @param solver Sinkhorn solver variant. `"internal"` (default) and `"auto"`
#'   both use the primal-domain scaling iteration. `"log_domain"` uses a
#'   numerically stable log-potential variant that avoids kernel underflow.
#'   `"annealing"` runs the log-domain solver across a geometric epsilon
#'   cooling schedule (starting at `anneal_start * epsilon`, cooling to
#'   `epsilon` over `anneal_stages` stages), warm-starting each stage from
#'   the previous stage's dual potentials. Recommended for label-guided fits
#'   or any fit with small `epsilon` (< 0.05) where cold-start Sinkhorn
#'   is slow or non-convergent.
```

### 9.3 New `@param` entries for tuning args

Add below `@param solver` in `somalign_fit`, `somalign_fit_anchored`,
`somalign_fit_two_pass`, and `somalign_sensitivity_grid`:

```
#' @param anneal_start Positive scalar >= 1. When `solver = "annealing"`,
#'   the starting epsilon is `anneal_start * epsilon`. Default `10`.
#'   Ignored when `solver != "annealing"`.
#' @param anneal_stages Positive integer. Number of cooling stages in the
#'   annealing schedule, including the final stage at the target `epsilon`.
#'   Default `10L`. A value of `1` degenerates to a cold-start log-domain
#'   solve. Ignored when `solver != "annealing"`.
#' @param anneal_factor Positive scalar < 1, or `NULL` (default). When not
#'   `NULL`, overrides the auto-computed per-stage cooling ratio. Use only
#'   when you need a non-geometric schedule. Ignored when
#'   `solver != "annealing"`.
```

---

## 10. BiocCheck Constraints (≤50-line exported bodies)

The three exported functions that gain new args are `somalign_fit`,
`somalign_fit_anchored`, and `somalign_fit_two_pass`. The new `anneal_*`
parameters add 3 lines to each function signature and 0 lines to the function
body (all dispatch logic is in internal helpers). The bodies remain well within
50 lines.

`somalign_sensitivity_grid` currently has a ~55-line body (lines 67–121); the
new args add only 3 formals and the `match.arg` line gains 2 extra choices. No
BiocCheck violation risk.

New internal helpers `.somalign_solve_annealing` (~45 lines) and
`.somalign_cooling_schedule` (~15 lines) are not exported and are not subject to
the 50-line rule.

---

## 11. Effort, Risks, and Dependencies

### 11.1 Effort

**Total: ~2 developer-days.**

| Task | Lines changed/added | Estimated time |
|------|---------------------|----------------|
| Warm-start interface in `.somalign_solve_internal_log` | +5 lines | 0.5 h |
| `.somalign_cooling_schedule` helper | +15 lines | 0.5 h |
| `.somalign_solve_annealing` driver | +45 lines | 2 h |
| `.somalign_solve_ot` dispatch branch | +15 lines | 0.5 h |
| Signature plumbing (4 exported fns + 2 internals) | ~30 lines | 1.5 h |
| Diagnostics wiring | +5 lines | 0.5 h |
| Tests (6 test cases) | ~100 lines | 3 h |
| Documentation | ~30 lines | 1 h |
| **Total** | **~250 lines** | **~9–10 h** |

### 11.2 Risks

| Risk | Likelihood | Mitigation |
|------|-----------|------------|
| Interior-stage non-convergence produces bad warm start | Medium | Use `inner_tol = 10 * tol`; enough stages so each stage is small; test 8.3 catches regression |
| `tau_a`/`tau_b` cached incorrectly in future refactor | Low | Document prominently in code comment; test 8.4 (single-stage equality) would catch it |
| Total iteration budget split too thinly for small codebooks | Low | `inner_iter = max(ceil(max_iter/n), 20L)` floor; test 8.1 verifies accuracy |
| Performance regression on standard solver paths | None | New code only reached via `solver == "annealing"` branch |

### 11.3 Dependencies shared with Ideas #1 and #3

- **Idea #1 (epsilon-sweep diagnostic)**: Idea #1 proposes exposing `log Z =
  sum(lse_g)` from `.somalign_solve_internal_log`. If implemented first, it
  also requires adding to the return list — which is compatible with the `f`/`g`
  addition planned here. Coordinate to do a single return-list expansion.
- **Idea #3 (if any per-epsilon machinery)**: Any epsilon-schedule infrastructure
  shared with Idea #3 should live in `.somalign_cooling_schedule` to avoid
  duplication. The geometric schedule function is generic.

### 11.4 Implementation order

1. Add `f`/`g` to `.somalign_solve_internal_log` return + `f_init`/`g_init`
   args (also needed by Idea #1 if it exposes log Z — coordinate).
2. Add `.somalign_cooling_schedule`.
3. Add `.somalign_solve_annealing`.
4. Update `.somalign_solve_ot` dispatch.
5. Plumb `anneal_*` args through all exported fns and internals.
6. Update `.somalign_build_diagnostics`.
7. Write tests; run `devtools::check()`.
8. Update roxygen.
