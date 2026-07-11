# somalign 0.99.0

* Added `somalign_train_reference()` and `somalign_reference()` to train or
  wrap an existing `kohonen` SOM as a fixed reference object, including node
  mass, per-node label probability, and per-node distance quantile summaries.

* Added `somalign_reference_from_nodes()` to reconstruct a reference object
  from pre-computed node-level artifacts (codebook, center, scale, masses,
  label probabilities, distance quantiles) without needing the original data.

* Added `somalign_query()` to scale query data against reference parameters,
  optionally train a query SOM, and return a `somalign_query` object
  containing per-sample node assignments and distances.

* Added `somalign_fit()` to align a query SOM to the fixed reference SOM via
  codebook-level unbalanced entropic optimal transport. Uses an internal
  pure-R generalized Sinkhorn solver; `"auto"` is retained as a compatibility
  alias for the internal solver.

* Added `somalign_results()` to extract a per-sample data frame with direct
  reference projection columns (canonical) and transport-corrected projection
  columns (auxiliary, for annotation and visualisation).

* Added `somalign_diagnostics()` to return solver, OT plan, node-level, and
  projection diagnostics from a `somalign_fit` object.

* Added `somalign_sensitivity_grid()` to sweep a grid of OT hyperparameters
  (`epsilon`, `rho_query`, `rho_ref`) and return a summary data frame for
  tuning and stability assessment. Supports optional parallelisation via
  `BiocParallel::bplapply()` when `parallel = TRUE`.

* **Breaking change** — `somalign_fit()` default `epsilon` changed from `0.05`
  to `0.5`. The cost matrix is now normalised by its median positive entry
  before the Sinkhorn kernel is computed, making `epsilon` scale- and
  dimension-invariant. Users who previously obtained meaningful regularisation
  with custom `epsilon` values should re-tune `epsilon` against the new
  normalised cost scale (inspect `diagnostics$solver$cost_scale`).

* `somalign_fit()` diagnostics now include `diagnostics$solver$converged`
  (logical), `diagnostics$solver$final_delta` (final iterate change), and
  `diagnostics$solver$cost_scale` (cost normalisation factor). A warning is
  now emitted when more than 50 % of query mass is destroyed during transport
  or when more than 50 % of query samples fall outside reference distance
  thresholds. A non-finite iterate delta also triggers an explicit warning
  (previously such degenerate solves could finish silently).

* `somalign_query()` gains a `codebook_space` argument (default
  `"reference_scaled"`). When a user-supplied `som_query` is provided, setting
  `codebook_space = "raw"` instructs `somalign_query()` to re-scale the
  codebook into reference-scaled space, mirroring the existing behaviour of
  `somalign_reference()`.

* `somalign_reference_from_nodes()` now emits informational messages when no
  label probabilities or distance quantiles are supplied, so users are aware
  that label transfer and outside-reference detection are disabled.

* `somalign_validate_ot_inputs()` (internal) now warns when all query or
  reference node masses are zero, preventing a silent all-zero transport plan.

* Fixed a silent fallback in `.somalign_thresholds()`: when the requested
  distance-quantile column is absent, a warning is now emitted instead of
  silently substituting the third column by position.
