# somalign — Refactor long functions into internal helpers

## Context

BiocCheck warns on every function longer than 50 lines. Seven functions
in `somalign` currently trip that check. The goal is to shorten them by
extracting cohesive blocks into **internal, non-exported helpers**
(dot-prefixed `.somalign_*`, no `@export`) so the warnings clear
**without enlarging the user-facing API** — users still see the same
exported functions; the new helpers are invisible in the namespace and
generate no help pages.

Key constraint: **BiocCheck flags *all* functions over 50 lines,
exported or not.** So the goal is not “move one big block into one
helper” but “decompose each long function into several balanced pieces,
each comfortably under 50 lines.” A 90-line internal helper would just
relocate the warning.

Second constraint: **R cannot be executed in this environment**, so the
refactor is verified in-session by read-through only (free-variable and
return-contract audit). Running the test suite and BiocCheck is a
user-side handoff.

## Targets (current body span, inclusive)

| Function | <File:line> | Lines | Action |
|----|----|----|----|
| `somalign_fit` | R/fit.R:56 | ~149 | Thin orchestrator + 5 helpers |
| `.somalign_solve_internal` | R/ot.R:33 | ~80 | Extract kernel setup + convergence warnings (keep loop) |
| `somalign_reference_from_nodes` | R/reference.R:153 | ~69 | Extract validation, messages, quantile resolution |
| `.somalign_transfer_labels` | R/fit.R:206 | ~66 | Extract empty-case builder + second-label ranking |
| `somalign_sensitivity_grid` | R/diagnostics.R:56 | ~65 | Extract row-summary + parallel/sequential dispatch |
| `somalign_reference` | R/reference.R:74 | ~54 | Extract center/scale resolution (marginal) |
| `.somalign_get_codebook` | R/utils.R:141 | ~50 | Extract codebook extraction (borderline) |

`somalign_train_reference` (R/reference.R:21) is now ~24 lines after
earlier refactoring and is **no longer a target**.

## Conventions (match existing code)

- New helpers are `.somalign_`-prefixed, defined at the top level of the
  same file as their caller, with **no roxygen block** (like
  `.somalign_row_normalize`, `.somalign_solve_ot`). No `@export`, no
  `@keywords internal`.
- Keep the existing `pkg::fun` call style inside moved code
  (e.g. [`stats::median`](https://rdrr.io/r/stats/median.html) at
  fit.R:73). Unexported helpers using `pkg::fun` generate **no NAMESPACE
  entry, no `.Rd`, and no new `importFrom`** — this is what makes the
  refactor safe to do without the roxygen tooling we can’t run here.

## Extraction seams (per function)

### `somalign_fit` (R/fit.R) — thin orchestrator

Keep the signature, the `.somalign_check_*` guards, `match.arg`, and the
final call chain. Extract five helpers: -
`.somalign_align_transport(...)` — lines 72–100: cost, `cost_scale`,
`.somalign_solve_ot`, `correspondence`, row/col mass,
`match_mass_ratio`/`match_fraction`, and the `match_mass_ratio > 1`
**[`message()`](https://rdrr.io/r/base/message.html) (verbatim)**.
Returns a list of all of these. -
`.somalign_project_pair(query, reference, node_shifts, chunk_size)` —
lines 120–123: `direct`, `corrected`, `correction_norm`. -
`.somalign_build_diagnostics(transport, query, reference, node_shifts, projection, epsilon, rho_query, rho_ref)`
— lines 125–161: pure list assembly. **Field names must survive
verbatim** — `diagnostics$solver$converged/ final_delta/cost_scale` and
the `ot`/`nodes`/`projection` sublists are asserted by tests. If this
lands near 50, split the `ot`/`nodes` sublists into a nested helper. -
`.somalign_fit_warnings(diagnostics)` — lines 163–184: the
mass-destruction and outside-fraction
**[`warning()`](https://rdrr.io/r/base/warning.html)s, text verbatim**
(untested but behavioral). - `.somalign_new_fit(...)` — lines 186–203:
the `structure(list(...), class = "somalign_fit")` assembly.

Resulting `somalign_fit` body ≈ signature + ~8 call statements → under
50.

### `.somalign_solve_internal` (R/ot.R) — extract edges, keep the loop

**Do not carve the `u`/`v`/`delta` iteration (lines 66–83)** — that is
where the mutable state lives. Extract only the edges: -
`.somalign_sinkhorn_kernel(cost, epsilon)` — lines 41–58: `k_raw`, the
underflow **[`warning()`](https://rdrr.io/r/base/warning.html)
(verbatim, incl. the `...underflowed` text)**, returns floored `k`. -
`.somalign_warn_convergence(final_delta, iterations, max_iter, tol)` —
lines 85–105: the non-finite-delta and non-convergence
**[`warning()`](https://rdrr.io/r/base/warning.html)s (verbatim)**. Keep
`tau_a`/`tau_b`, `u`/`v` init, the loop, `converged`, and plan assembly
in the body.

### `.somalign_transfer_labels` (R/fit.R)

- `.somalign_empty_label_transfer(n_nodes, match_fraction)` — lines
  212–223: the no-label early-return `data.frame` (same column
  names/types).
- `.somalign_second_labels(...)` — lines 236–248: second-label /
  second-confidence computation. Keep the main probs → top-label →
  entropy → acceptance flow and the final `data.frame`.

### `somalign_reference_from_nodes` (R/reference.R)

- `.somalign_validate_node_codebook(codebook, features)` — lines
  161–176: matrix coercion, feature-vector checks, feature selection,
  finite check; returns the prepared codebook.
- `.somalign_warn_from_nodes(label_prob, distance_quantiles)` — lines
  183–194: the two [`message()`](https://rdrr.io/r/base/message.html)s
  (**verbatim**).
- `.somalign_resolve_global_quantiles(distance_quantiles, global_distance_quantiles)`
  — lines 196–202.

### `somalign_sensitivity_grid` (R/diagnostics.R)

`.run_one` captures `...`, `query`, `reference`, `solver`, `grid`
lexically — **keep it as a local closure** (hoisting it loses free `...`
forwarding). Extract: -
`.somalign_grid_row_summary(fit, epsilon, rho_query, rho_ref)` — lines
86–98: the pure per-row `data.frame` builder that `.run_one` calls. -
`.somalign_run_grid(n, run_one, parallel)` — lines 101–119: the
BiocParallel-vs-sequential dispatch, taking the closure as an argument.

### `somalign_reference` (R/reference.R) — marginal (54 → ~42)

- `.somalign_resolve_center_scale(center, scale, data)` — lines 85–96:
  fill `center`/`scale` when either is `NULL`. Returns a list.

### `.somalign_get_codebook` (R/utils.R) — borderline (50)

- `.somalign_extract_codes(som, what)` — lines 142–168: the “pull a
  codebook out of a matrix / kohonen object / `$codes` list” logic. Keep
  feature-selection/ validation (169–189) in the caller.

## Safety discipline (because R is not run here)

For **each** extraction, before moving the block: 1. **Free-variable
audit.** List every symbol the block *reads*; confirm each is either a
parameter passed into the new helper or defined *inside* the moved
block. R’s lexical scoping means a top-level helper sees only its args +
the package namespace — a forgotten free variable errors (“object not
found”) or, worse, silently binds a wrong same-named symbol. This is the
primary risk. 2. **Return-contract audit.** Confirm the helper returns
exactly what the caller consumes, including list field names and the
`attr(shifts, "correction_allowed")` attribute contract from
`.somalign_node_shifts` that `.somalign_build_diagnostics` reads. 3.
**Byte-identical behavior.** Every
[`warning()`](https://rdrr.io/r/base/warning.html),
[`message()`](https://rdrr.io/r/base/message.html),
[`stop()`](https://rdrr.io/r/base/stop.html) string, every diagnostic
field name, and every returned attribute is copied verbatim. Tests
assert `expect_warning(..., "zero")`, the
`converged/final_delta/ cost_scale` fields, and object structure;
untested messages are still behavior.

## Not touched

- `NAMESPACE` and `man/*.Rd` — no exports change, so no regeneration
  needed.
- No new `importFrom` — moved code keeps `pkg::fun` calls.
- No behavior, signature, default, or numeric change to any exported
  function.

## Suggested implementation order (safest first)

1.  `somalign_fit` diagnostics/return/warning helpers (pure
    construction, near-zero risk).
2.  `somalign_fit` `.somalign_align_transport` +
    `.somalign_project_pair`.
3.  `.somalign_transfer_labels`, `somalign_reference_from_nodes`,
    `somalign_sensitivity_grid`.
4.  `somalign_reference`, `.somalign_get_codebook` (marginal — done for
    cleanliness).
5.  `.somalign_solve_internal` edges **last** (most delicate; loop stays
    intact).

After each function, re-check the caller’s span drops under 50 with
margin and that no new helper is itself over 50.

## Verification (user-side handoff — cannot run in this environment)

In-session verification is a free-variable / return-contract
read-through only. Run locally:

1.  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    — confirm **no** `NAMESPACE` diff (proves nothing became exported
    and no `importFrom` is missing).
2.  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
    — full `tests/testthat/` suite must pass unchanged; watch
    `test-ot-correctness.R` (warning strings, `diagnostics$solver`
    fields), `test-training-integration.R`, and
    `test-chunked-projection.R`.
3.  `BiocCheck::BiocCheck(".")` — confirm the “\> 50 lines” warnings for
    the seven targets are gone and no new function trips the check.
4.  `rcmdcheck::rcmdcheck(args = "--as-cran")` — confirm clean.
