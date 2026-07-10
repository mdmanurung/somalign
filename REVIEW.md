# somalign — Thorough Code Review

**Reviewer:** Claude Code (claude-opus-4-8)  
**Date:** 2026-07-10  
**Package version:** 0.0.0.9000  
**R version:** 4.5.1 (R4_51 conda env)  

All non-trivial claims in this review were **empirically demonstrated**, not merely inferred
from static analysis. Evidence is referenced inline.

---

## 1. Package Health (R CMD check + test suite)

### 1.1 R CMD check

```
0 errors | 0 warnings | 1 NOTE
```

The single NOTE is:
```
* checking for future file timestamps ... NOTE
  unable to verify current time
```
This is an HPC sandboxed-clock artifact (unrelated to the package).

**Conclusion:** The package is structurally sound. `R CMD check --as-cran` passes cleanly.

### 1.2 Test suite

```
FAIL 0 | WARN 0 | SKIP 1 | PASS 49
```

The 1 skip is the POT comparison test, correctly gated by `SOMALIGN_RUN_POT_TESTS=true`.
All 49 tests pass. The suite covers:
- Reference/query validation and error messages (`test-validation.R`, `test-reference.R`)
- OT solver dispatch, label transfer gating, results shape (`test-ot-labels-results.R`)
- End-to-end pipeline with real `kohonen::som()` training and sensitivity grid
  (`test-training-integration.R`)

---

## 2. Severity-Ranked Findings

### 🔴 HIGH — Silent non-convergence (`R/ot.R:86–100`)

**What:** The generalized Sinkhorn loop (`for (iter in seq_len(max_iter))`) exits silently
when `max_iter` iterations are reached without converging. The only signal is
`iterations == max_iter` stored in `fit$diagnostics$solver$iterations`, which users are
unlikely to inspect.

**Demonstrated (`benchmarks/numerical_experiments.R`, Experiment 3 & 4):**
```
epsilon=0.01: converged=FALSE after 1000 iters  (default eps is 0.05 — close!)
epsilon=0.001: delta=0.74 after 1000 iters  (very far from convergence)

somalign_fit with max_iter=2: iterations=2, warnings/messages='(none)'
```
Non-convergence is common at user-accessible parameter values and is **never communicated**.

**Why it matters:** A non-converged transport plan can have incorrect mass distribution,
leading to systematically wrong label transfers and node shifts — the core outputs of the
package. There is no indication anything went wrong.

**Suggested fix:**
```r
# In .somalign_solve_internal, after the loop:
if (iterations == max_iter) {
  warning(
    "Sinkhorn solver did not converge in ", max_iter, " iterations ",
    "(delta = ", signif(delta, 3), "). ",
    "Try increasing `max_iter`, `epsilon`, or `rho_*`.",
    call. = FALSE
  )
}
```

---

### 🔴 HIGH — Kernel underflow silently destroys cost information (`R/ot.R:78–80`)

**What:** The Sinkhorn kernel `k <- exp(-cost / epsilon)` is floored at `.Machine$double.xmin`
(`~2.2e-308`) to keep iterations alive. At sufficiently small `epsilon`, **all** kernel entries
underflow below this floor and become identical. The algorithm then runs 1000 iterations
on a cost-independent kernel, producing a uniform plan with no transport structure.

**Demonstrated (Experiment 2):**
```
epsilon    iters    K_min         plan_sum   max_plan_row   degenerate?
1.0e-04    1000     0.000e+00     0.994      0.497          YES (floor hit)
1.0e-05    1000     0.000e+00     1.000      0.500          YES (floor hit)
1.0e-08    1000     0.000e+00     1.000      0.500          YES (floor hit)
```
At `eps=1e-4` with costs O(1), the plan is essentially uniform regardless of the actual
cost matrix. No warning is emitted.

**Why it matters:** The cost matrix encodes the meaningful structure (codebook distances).
A uniform plan means every query node maps equally to all reference nodes, destroying label
transfer and correction accuracy.

**Suggested fix:** Warn when the raw kernel has entries below the floor:
```r
k_raw <- exp(-cost / epsilon)
dominated <- sum(k_raw < .Machine$double.xmin) / length(k_raw)
if (dominated > 0.01) {
  warning(
    round(100 * dominated), "% of kernel entries underflowed at epsilon=", epsilon,
    ". Consider increasing epsilon. Current safe lower bound: ",
    signif(-max(cost) / log(.Machine$double.xmin), 2),
    call. = FALSE
  )
}
k <- pmax(k_raw, .Machine$double.xmin)
```

---

### 🟡 MEDIUM — Dense O(n_samples × n_nodes) distance matrix is the memory hotspot
(`R/utils.R:200–201`, `R/fit.R:73–75`)

**What:** `.somalign_nearest_code()` allocates a full `n_samples × n_nodes` matrix via
`outer()`:
```r
d2 <- outer(rowSums(x * x), rowSums(codebook * codebook), "+") -
      2 * tcrossprod(x, codebook)
```
This is called:
1. Once in `somalign_query()` to map samples to query nodes.
2. **Twice** in `somalign_fit()` (lines 73 and 75, for direct and corrected projection).

**Demonstrated (benchmark Section 2):** At n=100,000 and a 10×10 grid (100 nodes), each
call allocates ~80 MB. The three calls per pipeline consume ~240 MB peak for the distance
matrices alone, not counting other temporaries.

**Memory formula:** `n_samples × n_nodes × 8 bytes` per call.

| n_samples | n_nodes | per call |
|----------:|--------:|---------:|
| 10,000    | 100     | 8 MB     |
| 100,000   | 100     | 80 MB    |
| 1,000,000 | 100     | 800 MB   |
| 100,000   | 400     | 320 MB   |

**Why it matters:** At cytometry scales (100k–1M cells, larger grids) this causes OOM or
heavy swapping, degrading all other in-memory computation.

**Suggested fix (short-term):** Chunk the rows of `x` and compute nearest-code per chunk,
reducing peak allocation to `chunk_size × n_nodes`. Example:
```r
.somalign_nearest_code_chunked <- function(x, codebook, chunk_size = 10000L) {
  n <- nrow(x)
  unit <- integer(n); distance <- numeric(n)
  for (i in ceiling(seq(1, n, by = chunk_size))) {
    idx <- i:min(i + chunk_size - 1L, n)
    res <- .somalign_nearest_code(x[idx, , drop=FALSE], codebook)
    unit[idx] <- res$unit; distance[idx] <- res$distance
  }
  list(unit = unit, distance = distance)
}
```
**Longer-term:** Rcpp or `FNN::knn.index` (which uses KD-trees and is much faster for
many features).

---

### 🟡 MEDIUM — `match_mass_ratio > 1` silently clamped (`R/fit.R:52–53`)

**What:**
```r
match_mass_ratio <- ifelse(query$node_masses > 0, row_mass / query$node_masses, 0)
match_fraction <- pmin(match_mass_ratio, 1)
```
`match_mass_ratio` can exceed 1.0 in unbalanced OT (the solver redistributes mass), yet
it's silently clamped to 1. This hides an anomaly that could indicate a numerical issue
or a very ill-conditioned problem.

**Suggested fix:** Emit a debug message or include `match_mass_ratio > 1` counts in
diagnostics so the user can inspect.

---

### 🟡 MEDIUM — No OT correctness tests in the test suite

**What:** None of the 49 passing tests verifies that:
- Plan row/col sums are approximately equal to target masses (unbalanced tolerance).
- The internal Sinkhorn plan agrees with the POT reference implementation.
- The plan is non-degenerate for a known well-conditioned input.

The tests only check that the plan is finite, non-negative, and has the right shape.

**Why it matters:** Structural checks pass even when the solver produces a degenerate plan
(as demonstrated in Experiment 2 for small `epsilon`).

**Suggested fix:** Add a test like:
```r
test_that("Sinkhorn marginals are approximately preserved for well-conditioned input", {
  set.seed(1)
  cost <- matrix(runif(6, 0.1, 2), 2, 3)
  a <- c(0.5, 0.5); b <- c(1/3, 1/3, 1/3)
  plan <- somalign:::.somalign_solve_internal(cost, a, b, 0.1, 1, 1, 1000, 1e-9)$plan
  expect_equal(rowSums(plan), a, tolerance = 0.15)  # unbalanced tolerance
  expect_equal(colSums(plan), b, tolerance = 0.10)
  expect_true(max(abs(plan)) > 0.01)  # non-degenerate
})
```

---

### 🟡 MEDIUM — Marginal deviation is large and undocumented

**What:** At the default `rho_query = rho_ref = 1, epsilon = 0.05`, the row-sum deviation
from target masses can be **13% or more** (Experiment 1: `max |row_sum - a_i| = 0.132` for
`a_i = 0.25`). This is correct behaviour for unbalanced OT — mass destruction is the point
of the `rho` parameters — but it is nowhere documented in `?somalign_fit` or the vignettes.

New users who interpret the transport plan as a row-stochastic assignment matrix will be
misled.

**Suggested fix:** Add a note to `?somalign_fit` explaining that row/col sums of the raw
`transport_plan` will not sum to the node masses (they sum to less, by design), and point
users to `fit$diagnostics$ot$max_row_mass_error` and `max_col_mass_error` for quantification.

---

### 🟡 MEDIUM — Per-node loops in several utility functions

Three functions iterate over `n_nodes` in R loops where vectorised operations would suffice:

1. **`.somalign_distance_quantiles` (`R/utils.R:246–255`):** Calls `stats::quantile()`
   inside a `for` loop over `n_nodes`. Equivalent to `apply(matrix_of_distances, 1, quantile)`.
   Fine at small grids; measurably slow at 400 nodes.

2. **`.somalign_reference_top_labels` (`R/results.R:73–81`):** Iterates over `n_nodes`
   calling `which.max(row)`. Vectorisable with `max.col(label_prob)`.

3. **`.somalign_label_probabilities` (`R/utils.R:288–297`):** Per-node loop to compute
   label counts. Could use sparse matrix operations or `tapply`.

**Suggested fixes:**
```r
# .somalign_reference_top_labels — replace loop with:
idx       <- max.col(label_prob, ties.method = "first")
label     <- colnames(label_prob)[idx]
confidence <- label_prob[cbind(seq_len(nrow(label_prob)), idx)]
label[rowSums(label_prob) == 0] <- NA_character_
confidence[rowSums(label_prob) == 0] <- NA_real_
```

---

### 🟡 MEDIUM — POT test skipped by default; no CI integration

**What:** The only test of the POT solver path requires `SOMALIGN_RUN_POT_TESTS=true`.
It only checks plan dimensions and finiteness, not numerical agreement with the internal
solver. The test file does not test the label transfer output under POT, only the transport
plan shape.

Now that POT is installed, this test can be activated. The cross-comparison in the
benchmark (Section 5) shows the solvers agree to within `~1e-4` max absolute difference —
this should become a formal test.

---

### 🟢 LOW — Noisy `UserWarning` from POT (`R/ot.R`, `.somalign_solve_pot`)

**What:** When calling POT's `sinkhorn_unbalanced` with `reg_type = "entropy"`, POT emits:
```
UserWarning: If reg_type = entropy, then the matrix c is overwritten by the one matrix.
```
`c` here is POT's **internal regularization reference measure** (an auxiliary parameter,
defaulting to `a bᵀ`), **not** the cost matrix. The cost is passed as `M = cost` (line 57
of `R/ot.R`). The warning is benign — POT locally resets `c` to a ones-matrix under
entropic regularization, which is the mathematically correct behaviour for the KL-unbalanced
variant. It has no effect on the returned transport plan.

**Why it matters (minor):** The warning leaks through reticulate to stderr, which can alarm
users. It appears on every call (once per sensitivity-grid iteration if solver = "pot").

**Suggested fix:** Suppress via Python's `warnings` module before the call:
```r
# In .somalign_solve_pot, before fn(...)
reticulate::py_run_string("import warnings; warnings.filterwarnings('ignore', category=UserWarning, module='ot')")
```
or upgrade POT (the warning may be removed in a future release).

---

### 🟢 LOW — Missing `\examples{}` in all 9 man pages

**What:** Every exported function's `.Rd` file has no `\examples{}` section.
`R CMD check` shows `checking examples ... NONE` without a NOTE in this environment,
but CRAN's submission portal emits a NOTE for documented functions with no examples.

**Affected files:** `man/somalign_diagnostics.Rd`, `man/somalign_fit.Rd`,
`man/somalign_query.Rd`, `man/somalign_reference.Rd`,
`man/somalign_reference_from_nodes.Rd`, `man/somalign_results.Rd`,
`man/somalign_sensitivity_grid.Rd`, `man/somalign_train_reference.Rd`,
`man/somalign-package.Rd`.

**Suggested fix:** Add at least a `\dontrun{}` or `\donttest{}` block to each function.
The vignette code is a natural source.

---

### 🟢 LOW — Placeholder author/email in DESCRIPTION

**What:**
```
Authors@R: person("somalign contributors", role = c("aut", "cre"), email = "noreply@example.com")
```
This must be replaced with a real name and contact before any CRAN submission.

---

### 🟢 LOW — `bench` and `microbenchmark` not in Suggests

The benchmark scripts in `benchmarks/` depend on these packages but they are not declared
in `DESCRIPTION`. Not an R CMD check error (benchmarks are not in `tests/`), but worth
noting for reproducibility.

---

### 🟢 LOW — `somalign_sensitivity_grid` is unparallelized (`R/diagnostics.R:44–54`)

The `for` loop over the parameter grid runs all fits sequentially. For large grids
(many epsilon × rho combinations) this is the bottleneck. `parallel::mclapply` or
`future.apply::future_lapply` would give easy multi-core speedup.

---

### 🟢 LOW — `NEWS.md` is a stub

`NEWS.md` exists but contains only a placeholder. It should track changes as the package
evolves toward 0.1.0.

---

## 3. Design Assessment

### 3.1 Conservative-primary / auxiliary-corrected split

The package makes the right architectural choice: **direct nearest-node projection is the
primary result** (`old_som_unit`, `old_som_label`, `final_status`, `outside_reference_distance`),
and OT-corrected columns are auxiliary. The docstring in `R/results.R:5–8` makes this
explicit. This avoids the common pitfall of presenting a statistically regularised estimate
as the ground truth.

### 3.2 API surface

The exported API (`somalign_train_reference`, `somalign_reference`, `somalign_reference_from_nodes`,
`somalign_query`, `somalign_fit`, `somalign_results`, `somalign_diagnostics`,
`somalign_sensitivity_grid`) is well-sized: small enough to be learnable, complete enough
to support the main use cases (fresh SOM, pre-trained SOM, node-level artifacts).

The handling of pre-trained query SOMs requires the user to supply data in the
`reference_scaled` coordinate space and declare this explicitly via `codebook_space`
(enforced with `stop()`). The vignette demonstrates this correctly. The enforcement is
correct and the error message is informative.

### 3.3 Parameter defaults

| Parameter | Default | Assessment |
|-----------|---------|------------|
| `epsilon` | 0.05 | **Reasonable** for unit-scale codebook distances. Warn-range: < 1e-3. |
| `rho_query / rho_ref` | 1 | **May be too relaxed.** At rho=1, 13%+ of mass can be destroyed per node (Experiment 1). Users with balanced datasets may want rho >> 1. Document the trade-off. |
| `min_match_fraction` | 0.05 | **Conservative but defensible.** Allows label transfer even when only 5% of a query node's mass is matched. Document what this means. |
| `confidence_threshold` | 0.6 | **Reasonable** for top-1 label probability. |
| `max_iter` | 1000 | **Too high as a silent limit.** Should warn at convergence failure. |
| `tol` | 1e-7 | Fine for converged cases; never reached at low epsilon. |
| novelty quantile | 0.95 | Sensible default; documented in vignette. |

### 3.4 Numerical architecture

The squared-distance identity `‖x - y‖² = ‖x‖² + ‖y‖² - 2xᵀy` is used correctly in both
`.somalign_nearest_code` (R/utils.R:200–201) and `.somalign_pairwise_distance` (R/utils.R:209).
The `pmax(d2, 0)` guard against floating-point negatives is correct. Accuracy is machine
precision (Experiment 7: max diff vs. naive = 8.9e-16).

The generalized Sinkhorn iteration (`R/ot.R:87–104`) correctly implements the KL-unbalanced
variant with scaling factors `tau_a = ρ_q/(ρ_q + ε)` and `tau_b = ρ_r/(ρ_r + ε)`. The
`u[!is.finite(u)] <- 0` guard is correct for avoiding NaN propagation. The convergence
criterion (relative change in `u, v`) is standard.

### 3.5 Scalability summary (empirical — see `benchmarks/RESULTS.md`)

**Stage decomposition at n=10k, p=20, 10×10 grid:**

| Stage | Median ms | Alloc MB | Who owns it |
|-------|----------:|---------:|-------------|
| `somalign_train_reference` | 1,165 | 119 | `kohonen::som()` |
| `somalign_query` | 918 | 90 | `kohonen::som()` |
| `fit: cost matrix build` | 0.3 | 0.6 | somalign |
| **`fit: OT solve (internal)`** | **7** | **2.8** | **somalign core** |
| **`fit: project_samples`** | **104** | **55** | **somalign HOTSPOT** |
| `somalign_results` | 3 | 1.7 | somalign |

**`somalign_fit` end-to-end: 164 ms** (measured directly). The dominant contributions
are the two `project_samples` calls; OT is only 7 ms. The `kohonen` dependency contributes
~2 s — roughly 12× somalign's own fit cost. Note: standalone `project_samples` benchmarks
(~104 ms) reflect single-iteration cold-cache timing and may exceed the in-fit cost due to
GC pressure with `bench::mark()`'s warm-up behaviour.

**n_samples scaling (10×10 grid, single `project_samples` call):**

| n_samples | ms | alloc MB |
|----------:|---:|---------:|
| 1,000 | 44 | 6 |
| 10,000 | 98 | 55 |
| 100,000 | 396 | 555 |
| 1,000,000 | ~5,000 | ~800* |

*`proc.time()` estimate; full `bench::mark()` OOMed. Two calls per `somalign_fit` + one in `somalign_query` = 3 total per pipeline. At n=100k: ~1.7 GB peak for projections alone.

**Grid size scaling (n=10k, `somalign_fit` stages):**

| Grid | n_nodes | OT solve (R) | OT solve (POT) | project_samples |
|------|--------:|-------------:|---------------:|----------------:|
| 2×2 | 4 | 4 ms | 3.4 ms | 2 ms |
| 5×5 | 25 | 5 ms | 3.6 ms | 63 ms |
| 10×10 | 100 | 8 ms | 4.2 ms | 131 ms |
| 15×15 | 225 | 25 ms | 6.6 ms | 68 ms |
| 20×20 | 400 | 64 ms | 13 ms | 157 ms |

OT solve grows quadratically in n_nodes. POT's C backend is **5× faster** than pure-R
at 20×20 (13 vs 64 ms). Plan agreement: max|Δ| < 1e-7 at all grid sizes.
*(The `project_samples` column is non-monotonic at 10×10 → 15×15 → 20×20 due to GC
pressure in single-iteration `bench::mark()` runs; treat it as approximate. OT-solve and
cost-build trends are robust.)*

The bottleneck for large-scale use is **not the OT solver** but the dense distance
matrix in `project_samples`. Chunked projection or an approximate nearest-neighbour
backend (e.g., `FNN::knn.index`) would unlock true cytometry scales (1M+ cells).
At grids > 15×15, switching `solver = "pot"` gives a meaningful OT speedup.

---

## 4. Test Coverage Assessment

| Area | Covered? | Quality |
|------|----------|---------|
| Input validation (errors) | ✅ | Good — covers feature names, NA, Inf, zero variance |
| Object structure / S3 classes | ✅ | Good |
| Label transfer gating | ✅ | Tested |
| OT solver dispatch | ✅ | Tested (auto fallback, parameter validation) |
| OT plan correctness (marginals) | ❌ | **Gap** — only checks finite/non-negative |
| Convergence failure / non-convergence | ❌ | **Gap** — no test for max_iter hit |
| Underflow at small epsilon | ❌ | **Gap** — no regression test |
| Internal vs POT numerical agreement | ❌ | **Gap** — skipped by default |
| Pre-trained query SOM (coordinate system) | ✅ | Tested in integration suite |
| Sensitivity grid | ✅ | Shape and content tested |
| `somalign_reference_from_nodes` | ✅ | Tested (unknown_reference_distance case) |
| `somalign_diagnostics` structure | ✅ | Field names checked |

---

## 5. Documentation Assessment

| Item | Status |
|------|--------|
| All exported functions documented | ✅ |
| `\examples{}` in man pages | ❌ All 9 missing |
| Vignette 1 (`somalign.Rmd`) | ✅ Runs end-to-end correctly |
| Vignette 2 (`pretrained-old-and-new-soms.Rmd`) | ✅ Runs end-to-end correctly |
| `@return` tags complete | ✅ |
| Parameter defaults documented | ⚠️ Partial — `rho` trade-offs, marginal behavior undocumented |
| `NEWS.md` | ⚠️ Stub only |
| `DESCRIPTION` author | ❌ Placeholder |

---

## 6. Summary: Priority Action List

| # | Severity | Action |
|---|----------|--------|
| 1 | 🔴 HIGH | Emit `warning()` when Sinkhorn loop hits `max_iter` without convergence (`R/ot.R:100`) |
| 2 | 🔴 HIGH | Warn when kernel underflows at small epsilon (`R/ot.R:78–80`) |
| 3 | 🟡 MED | Add chunked or approximate nearest-code for large `n_samples` (`R/utils.R:197–205`) |
| 4 | 🟡 MED | Add OT marginal correctness test and non-convergence regression test |
| 5 | 🟡 MED | Document expected marginal deviation and `rho` mass-destruction trade-off |
| 6 | 🟡 MED | Vectorise per-node loops in `results.R`, `utils.R` (label prob, distance quantiles) |
| 7 | 🟡 MED | Enable and extend POT comparison test; add numerical agreement assertion |
| 8 | 🟢 LOW | Suppress noisy-but-benign POT `UserWarning` about ref-measure `c` (`R/ot.R`, `.somalign_solve_pot`) |
| 9 | 🟢 LOW | Add `\examples{}` to all 9 man pages |
| 10 | 🟢 LOW | Replace placeholder author/email in `DESCRIPTION` |
| 11 | 🟢 LOW | Add `parallel::mclapply` support to `somalign_sensitivity_grid` |

Fixes 1 and 2 are one-liners and should ship in the next commit. Fix 3 is the prerequisite
for cytometry-scale use. Fix 8 is a safety fix for the POT path — one line (`cost_copy <- cost`).
Fixes 4–7 strengthen correctness guarantees before a 0.1.0 release.
