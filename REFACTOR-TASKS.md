# somalign long-function refactor task list

This checklist implements `REFACTOR-PLAN.md` by extracting internal helpers only.
The goal is to clear BiocCheck function-length warnings without changing exported
functions, signatures, defaults, object structures, diagnostics fields, result
columns, warnings, messages, or numeric behavior.

## Global rules

- [ ] Work on one target function at a time and inspect the diff before moving to the next target.
- [ ] Add new helpers as top-level `.somalign_*` functions in the same `R/*.R` file as their caller.
- [ ] Do not add roxygen blocks, `@export`, `@keywords internal`, or any NAMESPACE entry for new helpers.
- [ ] Keep existing `pkg::fun` call style inside moved code, such as `stats::median` and `BiocParallel::bplapply`.
- [ ] Preserve every exported function signature and default value exactly.
- [ ] Preserve every `structure(..., class = ...)` class name exactly.
- [ ] Preserve every public list field name in `somalign_reference`, `somalign_query`, and `somalign_fit` objects exactly.
- [ ] Preserve every `somalign_results()` output column name and column order exactly.
- [ ] Preserve every diagnostic field name and sublist name exactly.
- [ ] Preserve every existing `stop()`, `warning()`, and `message()` string exactly.
- [ ] Before each extraction, write down the moved block's free variables and make each one either a helper argument or a variable defined inside the helper.
- [ ] Before each extraction, write down the helper return contract and confirm the caller consumes exactly those names or values.
- [ ] After each extraction, re-read the caller and helper for missing variables, changed evaluation order, changed warnings, and changed return objects.
- [ ] After each target, check that the original long function is below 50 lines with margin.
- [ ] After each target, check that no new helper is itself above 50 lines.

## Preflight

- [ ] Confirm the working tree state with `git status --short`.
- [ ] Confirm `REFACTOR-PLAN.md` is present and read it fully.
- [ ] Re-open these files before editing: `R/fit.R`, `R/ot.R`, `R/reference.R`, `R/diagnostics.R`, `R/utils.R`.
- [ ] Re-open API consumers before editing `somalign_fit`: `R/results.R` and `R/print.R`.
- [ ] Re-open relevant tests before editing: `tests/testthat/test-ot-labels-results.R`, `tests/testthat/test-ot-correctness.R`, `tests/testthat/test-ot-warnings.R`, `tests/testthat/test-reference.R`, `tests/testthat/test-training-integration.R`, and `tests/testthat/test-chunked-projection.R`.

## Phase 1: Refactor `somalign_fit` pure assembly helpers in `R/fit.R`

- [ ] In `R/fit.R`, identify the current `diagnostics <- list(...)` block in `somalign_fit`.
- [ ] Audit free variables for the diagnostics block: `ot`, `epsilon`, `rho_query`, `rho_ref`, `cost_scale`, `plan`, `row_mass`, `col_mass`, `query`, `reference`, `match_fraction`, `match_mass_ratio`, `node_shifts`, `direct`, and `corrected`.
- [ ] Define `.somalign_build_diagnostics(transport, query, reference, node_shifts, projection, epsilon, rho_query, rho_ref)` below `somalign_fit`.
- [ ] Make `.somalign_build_diagnostics()` read `ot`, `cost_scale`, `plan`, `row_mass`, `col_mass`, `match_fraction`, and `match_mass_ratio` from the `transport` list.
- [ ] Make `.somalign_build_diagnostics()` read `direct` and `corrected` from the `projection` list.
- [ ] Preserve the top-level diagnostics names exactly: `solver`, `ot`, `nodes`, `projection`.
- [ ] Preserve `diagnostics$solver` names exactly: `requested`, `used`, `notes`, `iterations`, `converged`, `final_delta`, `epsilon`, `rho_query`, `rho_ref`, `cost_scale`.
- [ ] Preserve `diagnostics$ot` names exactly: `transport_mass`, `row_mass`, `col_mass`, `query_mass`, `reference_mass`, `match_fraction`, `match_mass_ratio`, `max_row_mass_error`, `max_col_mass_error`.
- [ ] Preserve `diagnostics$nodes` data frame columns exactly: `query_node`, `query_mass`, `transported_mass`, `match_fraction`, `correction_allowed`, `correction_norm`.
- [ ] Preserve `diagnostics$projection` names exactly: `outside_direct_fraction`, `outside_corrected_fraction`.
- [ ] Confirm `.somalign_build_diagnostics()` uses `attr(node_shifts, "correction_allowed")` exactly as before.
- [ ] Replace the inline diagnostics block in `somalign_fit` with a call to `.somalign_build_diagnostics(...)`.
- [ ] Inspect the diff for any changed diagnostic name, changed order, or changed expression.

- [ ] Identify the mass-destruction and outside-fraction warning block in `somalign_fit`.
- [ ] Audit free variables for the warning block: `diagnostics` and `query$node_masses`.
- [ ] Define `.somalign_fit_warnings(diagnostics)` below `somalign_fit`.
- [ ] Inside `.somalign_fit_warnings()`, compute `query_total_mass <- sum(diagnostics$ot$query_mass)` so the helper does not need `query`.
- [ ] Move both warning calls into `.somalign_fit_warnings()` with text and `call. = FALSE` unchanged.
- [ ] Replace the inline warning block in `somalign_fit` with `.somalign_fit_warnings(diagnostics)`.
- [ ] Confirm `.somalign_fit_warnings()` returns invisibly or relies on the last expression only; the caller must not use its return value.
- [ ] Inspect the diff for warning text changes.

- [ ] Identify the final `structure(list(...), class = "somalign_fit")` block in `somalign_fit`.
- [ ] Audit free variables for the fit object block: `query`, `reference`, `cost`, `plan`, `correspondence`, `label_transfer`, `node_shifts`, `direct`, `corrected`, `correction_norm`, and `diagnostics`.
- [ ] Define `.somalign_new_fit(query, reference, transport, label_transfer, node_shifts, projection, diagnostics)` below `somalign_fit`.
- [ ] Make `.somalign_new_fit()` read `cost`, `plan`, and `correspondence` from `transport`.
- [ ] Make `.somalign_new_fit()` read `direct`, `corrected`, and `correction_norm` from `projection`.
- [ ] Preserve the `somalign_fit` list names exactly: `query`, `reference`, `cost`, `transport_plan`, `correspondence`, `label_transfer`, `node_shifts`, `projection`, `diagnostics`.
- [ ] Preserve the nested `projection` names exactly: `direct`, `corrected`, `correction_norm`.
- [ ] Preserve `class = "somalign_fit"` exactly.
- [ ] Replace the inline final structure block with `.somalign_new_fit(...)`.
- [ ] Inspect `R/results.R` mentally against the new object to confirm `fit$projection$direct`, `fit$projection$corrected`, `fit$projection$correction_norm`, and `fit$label_transfer` still exist.
- [ ] Inspect `R/print.R` mentally against the new object to confirm `fit$diagnostics$solver$used` and `fit$diagnostics$ot$transport_mass` still exist.

## Phase 2: Refactor `somalign_fit` transport and projection helpers in `R/fit.R`

- [ ] Identify the cost, scaling, OT solve, correspondence, mass-ratio, and `match_mass_ratio > 1` message block in `somalign_fit`.
- [ ] Audit free variables for the transport block: `query`, `reference`, `epsilon`, `rho_query`, `rho_ref`, `solver`, `max_iter`, and `tol`.
- [ ] Define `.somalign_align_transport(query, reference, epsilon, rho_query, rho_ref, solver, max_iter, tol)` below `somalign_fit`.
- [ ] Move the cost matrix calculation into `.somalign_align_transport()`.
- [ ] Move the `cost_scale <- stats::median(cost[cost > 0])` calculation into `.somalign_align_transport()`.
- [ ] Preserve the non-finite or zero `cost_scale` fallback to `1` exactly.
- [ ] Move the `.somalign_solve_ot(...)` call into `.somalign_align_transport()` with every argument unchanged.
- [ ] Move `plan`, `correspondence`, `row_mass`, `col_mass`, `match_mass_ratio`, and `match_fraction` calculation into `.somalign_align_transport()`.
- [ ] Move the `match_mass_ratio > 1` message into `.somalign_align_transport()` with text unchanged.
- [ ] Return a named list from `.somalign_align_transport()` containing exactly the fields consumed later: `cost`, `cost_scale`, `ot`, `plan`, `correspondence`, `row_mass`, `col_mass`, `match_mass_ratio`, and `match_fraction`.
- [ ] Replace the inline transport block in `somalign_fit` with `transport <- .somalign_align_transport(...)`.
- [ ] Update later `somalign_fit` references to use `transport$correspondence`, `transport$row_mass`, and `transport$match_fraction`.
- [ ] Inspect the diff for changed numeric expressions or changed message text.

- [ ] Identify the direct projection, corrected matrix, corrected projection, and per-sample correction norm block in `somalign_fit`.
- [ ] Audit free variables for the projection block: `query`, `reference`, `node_shifts`, and `chunk_size`.
- [ ] Define `.somalign_project_pair(query, reference, node_shifts, chunk_size)` below `somalign_fit`.
- [ ] Move `direct <- .somalign_project_samples(...)` into `.somalign_project_pair()`.
- [ ] Move `corrected_matrix <- query$scaled_data + node_shifts[query$sample_unit, , drop = FALSE]` into `.somalign_project_pair()`.
- [ ] Move `corrected <- .somalign_project_samples(...)` into `.somalign_project_pair()`.
- [ ] Move `correction_norm <- sqrt(rowSums(node_shifts[query$sample_unit, , drop = FALSE]^2))` into `.somalign_project_pair()`.
- [ ] Return a named list with exactly `direct`, `corrected`, and `correction_norm`.
- [ ] Replace the inline projection block in `somalign_fit` with `projection <- .somalign_project_pair(...)`.
- [ ] Confirm `somalign_fit` now calls, in order, checks, `match.arg`, `.somalign_align_transport()`, `.somalign_transfer_labels()`, `.somalign_node_shifts()`, `.somalign_project_pair()`, `.somalign_build_diagnostics()`, `.somalign_fit_warnings()`, and `.somalign_new_fit()`.
- [ ] Confirm the resulting `somalign_fit` body is comfortably under 50 lines.

## Phase 3: Refactor `.somalign_transfer_labels` in `R/fit.R`

- [ ] Identify the no-label early-return data frame in `.somalign_transfer_labels`.
- [ ] Audit free variables for the no-label block: `n_nodes` and `match_fraction`.
- [ ] Define `.somalign_empty_label_transfer(n_nodes, match_fraction)` below `.somalign_transfer_labels`.
- [ ] Move the no-label data frame into `.somalign_empty_label_transfer()` unchanged.
- [ ] Preserve column names exactly: `query_node`, `label`, `confidence`, `second_label`, `second_confidence`, `entropy`, `match_fraction`, `accepted`.
- [ ] Replace the inline no-label early return with `return(.somalign_empty_label_transfer(n_nodes, match_fraction))`.

- [ ] Identify the second-label and second-confidence block in `.somalign_transfer_labels`.
- [ ] Audit free variables for the second-label block: `probs_norm`, `top_idx`, `n_nodes`, `probs`, `has_mass`, and `label_names`.
- [ ] Define `.somalign_second_labels(probs_norm, top_idx, has_mass, label_names)` below `.somalign_transfer_labels`.
- [ ] Inside `.somalign_second_labels()`, compute `n_nodes <- nrow(probs_norm)`.
- [ ] Preserve the existing `ncol(probs) == 1L` behavior by checking `ncol(probs_norm) == 1L`.
- [ ] Preserve the zeroing of top probabilities before computing second labels.
- [ ] Preserve the rule that second labels with missing or zero confidence become `NA_character_`.
- [ ] Return a named list with exactly `second_label` and `second_confidence`.
- [ ] Replace the inline second-label block with a call to `.somalign_second_labels(...)`.
- [ ] Preserve the final `.somalign_transfer_labels()` data frame column names and order exactly.
- [ ] Confirm `.somalign_transfer_labels` and both new helpers are each under 50 lines.

## Phase 4: Refactor `somalign_reference_from_nodes` in `R/reference.R`

- [ ] Identify only the codebook/features validation block in `somalign_reference_from_nodes`: matrix coercion, feature-vector checks, duplicated feature checks, default colnames, feature selection, and finite codebook validation.
- [ ] Do not move `center <- .somalign_named_numeric(...)`, `scale <- .somalign_named_numeric(...)`, or `.somalign_validate_scale(scale)` into the codebook helper.
- [ ] Audit free variables for the codebook/features validation block: `codebook` and `features`.
- [ ] Define `.somalign_validate_node_codebook(codebook, features)` below `somalign_reference_from_nodes`.
- [ ] Move only the codebook/features validation block into `.somalign_validate_node_codebook()`.
- [ ] Preserve both `features` stop messages exactly.
- [ ] Preserve the behavior that missing codebook colnames are set to `features`.
- [ ] Return the prepared `codebook`.
- [ ] Replace the inline codebook/features validation block with `codebook <- .somalign_validate_node_codebook(codebook, features)`.
- [ ] Leave center/scale validation immediately after the helper call in the caller.

- [ ] Identify the two informational `message()` calls in `somalign_reference_from_nodes`.
- [ ] Audit free variables for the message block: `label_prob` and `distance_quantiles`.
- [ ] Define `.somalign_warn_from_nodes(label_prob, distance_quantiles)` below `somalign_reference_from_nodes`.
- [ ] Move both `message()` calls into `.somalign_warn_from_nodes()` with text unchanged.
- [ ] Replace the inline message block with `.somalign_warn_from_nodes(label_prob, distance_quantiles)`.

- [ ] Identify the global distance quantile resolution block in `somalign_reference_from_nodes`.
- [ ] Audit free variables for the quantile block: `distance_quantiles` and `global_distance_quantiles`.
- [ ] Define `.somalign_resolve_global_quantiles(distance_quantiles, global_distance_quantiles)` below `somalign_reference_from_nodes`.
- [ ] Move the `is.null(global_distance_quantiles)` fallback into `.somalign_resolve_global_quantiles()`.
- [ ] Preserve `global_distance_quantiles <- as.numeric(global_distance_quantiles)` exactly.
- [ ] Preserve the existing behavior of assigning names only when `names(global_distance_quantiles)` is `NULL`.
- [ ] Return the resolved `global_distance_quantiles`.
- [ ] Replace the inline quantile block with `global_distance_quantiles <- .somalign_resolve_global_quantiles(...)`.
- [ ] Preserve the final `somalign_reference` object list names and class exactly.
- [ ] Confirm `somalign_reference_from_nodes` and all new helpers are each under 50 lines.

## Phase 5: Refactor `somalign_sensitivity_grid` in `R/diagnostics.R`

- [ ] Identify the per-grid-row data frame construction inside local `.run_one`.
- [ ] Audit free variables for the row-summary block: `fit`, `grid$epsilon[i]`, `grid$rho_query[i]`, and `grid$rho_ref[i]`.
- [ ] Define `.somalign_grid_row_summary(fit, epsilon, rho_query, rho_ref)` below `somalign_sensitivity_grid`.
- [ ] Move `diag <- somalign_diagnostics(fit)` into `.somalign_grid_row_summary()`.
- [ ] Move the row data frame into `.somalign_grid_row_summary()` unchanged.
- [ ] Preserve row-summary column names and order exactly: `epsilon`, `rho_query`, `rho_ref`, `solver`, `transport_mass`, `mean_match_fraction`, `max_row_mass_error`, `max_col_mass_error`, `accepted_label_fraction`, `outside_direct_fraction`, `outside_corrected_fraction`.
- [ ] In local `.run_one`, keep the `somalign_fit(...)` call local so `...`, `query`, `reference`, `solver`, and `grid` remain lexically captured.
- [ ] Replace the inline row data frame in `.run_one` with `.somalign_grid_row_summary(fit, grid$epsilon[i], grid$rho_query[i], grid$rho_ref[i])`.

- [ ] Identify the parallel/sequential dispatch block in `somalign_sensitivity_grid`.
- [ ] Audit free variables for the dispatch block: `parallel`, `.run_one`, and `grid`.
- [ ] Define `.somalign_run_grid(n, run_one, parallel)` below `somalign_sensitivity_grid`.
- [ ] Move the `parallel` branch into `.somalign_run_grid()`.
- [ ] Preserve the `BiocParallel` availability check and error text exactly.
- [ ] Preserve the `BiocParallel::bplapply(seq_len(n), run_one)` behavior.
- [ ] Preserve the sequential `vector("list", n)` plus `for` loop behavior.
- [ ] Return `rows` from `.somalign_run_grid()`.
- [ ] Replace the inline dispatch block with `rows <- .somalign_run_grid(nrow(grid), .run_one, parallel)`.
- [ ] Preserve the final `do.call(rbind, rows)` in `somalign_sensitivity_grid`.
- [ ] Confirm `somalign_sensitivity_grid` and both new helpers are each under 50 lines.

## Phase 6: Refactor `somalign_reference` in `R/reference.R`

- [ ] Identify the center/scale resolution block in `somalign_reference`.
- [ ] Audit free variables for the center/scale block: `center`, `scale`, and `data`.
- [ ] Define `.somalign_resolve_center_scale(center, scale, data)` below `somalign_reference`.
- [ ] Move only the `if (is.null(center) || is.null(scale))` block into `.somalign_resolve_center_scale()`.
- [ ] Preserve the behavior that provided `center` or `scale` values are not recomputed.
- [ ] Return a named list with exactly `center` and `scale`.
- [ ] Replace the inline block with `resolved <- .somalign_resolve_center_scale(center, scale, data)`.
- [ ] Assign `center <- resolved$center` and `scale <- resolved$scale`.
- [ ] Leave `.somalign_named_numeric()` and `.somalign_validate_scale()` in the caller.
- [ ] Preserve the final `somalign_reference` object list names and class exactly.
- [ ] Confirm `somalign_reference` and the new helper are each under 50 lines.

## Phase 7: Refactor `.somalign_get_codebook` in `R/utils.R`

- [ ] Identify the code extraction block in `.somalign_get_codebook`: initialize `codes`, handle matrix/data frame input, try `kohonen::getCodes()`, inspect `$codes`, choose list element, error if still null, and coerce via `.somalign_as_matrix()`.
- [ ] Audit free variables for the extraction block: `som` and `what`.
- [ ] Define `.somalign_extract_codes(som, what)` below `.somalign_get_codebook`.
- [ ] Move the extraction block into `.somalign_extract_codes()` unchanged.
- [ ] Preserve the `requireNamespace("kohonen", quietly = TRUE)` behavior.
- [ ] Preserve the `tryCatch(kohonen::getCodes(som), error = function(e) NULL)` behavior.
- [ ] Preserve the `$codes` fallback behavior.
- [ ] Preserve the list-codebook selection behavior: prefer `"data"` when present, otherwise first element.
- [ ] Preserve the `Could not extract a SOM codebook` stop text exactly.
- [ ] Preserve `.somalign_as_matrix(codes, what = paste0(what, " codebook"))`.
- [ ] Replace the inline extraction block with `codes <- .somalign_extract_codes(som, what)`.
- [ ] Leave all feature-selection and finite-validation logic in `.somalign_get_codebook()`.
- [ ] Confirm `.somalign_get_codebook` and `.somalign_extract_codes` are each under 50 lines.

## Phase 8: Refactor `.somalign_solve_internal` in `R/ot.R`

- [ ] Identify the Sinkhorn kernel setup and underflow warning block in `.somalign_solve_internal`.
- [ ] Audit free variables for the kernel block: `cost`, `epsilon`, and `tiny`.
- [ ] Define `.somalign_sinkhorn_kernel(cost, epsilon, tiny)` below `.somalign_solve_internal`.
- [ ] Keep `tiny <- .Machine$double.xmin` in `.somalign_solve_internal` before the helper call because the iteration loop still uses `tiny`.
- [ ] Move `k_raw`, `underflow_fraction`, the underflow warning, and `pmax(k_raw, tiny)` into `.somalign_sinkhorn_kernel()`.
- [ ] Preserve the underflow warning text exactly, including `underflowed`.
- [ ] Preserve `safe_eps <- signif(-max(cost) / log(.Machine$double.xmin), 3)` exactly.
- [ ] Return only the floored kernel matrix `k`.
- [ ] Replace the inline kernel block with `k <- .somalign_sinkhorn_kernel(cost, epsilon, tiny)`.
- [ ] Confirm the `u`/`v`/`delta` iteration loop remains in `.somalign_solve_internal`.
- [ ] Confirm `tau_a`, `tau_b`, `u`, `v`, `delta`, and `iterations` remain in `.somalign_solve_internal`.

- [ ] Identify the non-finite-delta and non-convergence warning block in `.somalign_solve_internal`.
- [ ] Audit free variables for the convergence-warning block: `final_delta`, `iterations`, `max_iter`, and `tol`.
- [ ] Define `.somalign_warn_convergence(final_delta, iterations, max_iter, tol)` below `.somalign_solve_internal`.
- [ ] Move only the two convergence warning branches into `.somalign_warn_convergence()`.
- [ ] Preserve both convergence warning texts exactly.
- [ ] Replace the inline convergence warning block with `.somalign_warn_convergence(final_delta, iterations, max_iter, tol)`.
- [ ] Leave `converged <- is.finite(final_delta) && final_delta < tol` in `.somalign_solve_internal`.
- [ ] Leave plan assembly in `.somalign_solve_internal`.
- [ ] Confirm `.somalign_solve_internal`, `.somalign_sinkhorn_kernel`, and `.somalign_warn_convergence` are each under 50 lines.

## Static review after all edits

- [ ] Re-read every new helper and confirm it has no accidental roxygen block.
- [ ] Search for accidental exports with `rg -n "@export|export\\(" R NAMESPACE`.
- [ ] Confirm no new helper name appears in `NAMESPACE`.
- [ ] Search all edited files for changed public names with `rg -n "transport_plan|correspondence|label_transfer|node_shifts|projection|diagnostics|outside_direct_fraction|outside_corrected_fraction|final_delta|cost_scale" R tests`.
- [ ] Inspect `git diff -- R/fit.R R/ot.R R/reference.R R/diagnostics.R R/utils.R`.
- [ ] Confirm `R/results.R` did not need changes.
- [ ] Confirm `R/print.R` did not need changes.
- [ ] Confirm no `man/*.Rd` files changed unless roxygen was intentionally run later.
- [ ] Confirm no generated `docs/` pkgdown files changed.

## Functional verification

- [ ] Run `devtools::document()` in an R environment.
- [ ] Confirm `git diff -- NAMESPACE man` shows no unintended changes.
- [ ] Run `devtools::test()` in an R environment.
- [ ] If a test fails, inspect whether the failure is from changed object structure, changed warning/message text, changed result columns, or changed numeric output.
- [ ] Run `BiocCheck::BiocCheck(".")` in an R environment.
- [ ] Confirm the seven original function-length warnings are gone.
- [ ] Confirm BiocCheck does not report any new helper over 50 lines.
- [ ] Run `rcmdcheck::rcmdcheck(args = "--as-cran")` in an R environment.
- [ ] Capture the final test, BiocCheck, and R CMD check outcomes in the final report.

## Stop conditions

- [ ] Stop if an extraction requires changing an exported function signature or default value.
- [ ] Stop if an extraction requires changing the `somalign_fit`, `somalign_reference`, or `somalign_query` object field contract.
- [ ] Stop if an extraction changes `somalign_results()` columns or column order.
- [ ] Stop if `devtools::document()` creates new exports or unexpected `man/*.Rd` files.
- [ ] Stop if a helper remains above 50 lines after extraction and cannot be split without changing behavior.
- [ ] Stop if R is unavailable for final validation and report the read-through checks completed plus the exact commands still needing user-side execution.
