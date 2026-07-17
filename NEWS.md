# somalign 0.99.4

## New features

* `somalign_correct_expression()` returns a cell-level (cells by markers)
  batch-corrected marker expression matrix for downstream visualisation and
  differential expression. The correction is confined to the anchor-estimated
  batch subspace and smoothed across each cell's nearest SOM nodes with a
  Gaussian kernel, preserving variation orthogonal to the batch direction.
  Requires a subspace-aware fit from `somalign_fit_anchored(correction =
  "subspace"/"both")` or `somalign_fit_two_pass()`. It is an auxiliary
  correction aid, not the primary label-transfer product; see
  `?somalign_correct_expression` and `vignette("anchor-samples")`.

* `somalign_reference()`, `somalign_train_reference()`, and
  `somalign_reference_from_som()` now compute and store `reference$node_var`
  (per-node per-marker variance in reference-scaled space) by default
  (`compute_node_var = TRUE`; set `FALSE` to opt out).
  `somalign_reference_from_nodes()` accepts a pre-computed `node_var` matrix
  directly. `somalign_results()` uses it to add a calibrated chi-squared
  alternative to the distance-quantile `outside_reference_distance` flag:
  `outside_reference_surprisal`, `outside_reference_pvalue`, and
  `outside_reference_top_marker` (the single worst-contributing marker per
  cell -- useful for pinpointing marker-specific batch artifacts). An
  optional `outside_pvalue_threshold` argument adds a boolean
  `outside_reference_pvalue_flag` column. `NA` when the reference lacks
  `node_var`.

* `somalign_results()` now exposes `transferred_label_second`,
  `transferred_label_second_confidence`, and `transferred_label_margin`
  (top confidence minus second confidence, from the same query node's
  label transfer). `transferred_label` is an argmax over label transport
  probabilities and can be brittle when two labels receive close mass;
  the margin lets callers triage low-confidence label transfers
  (e.g. cells on nodes where the top and runner-up label are nearly tied)
  without re-deriving them from `fit$label_transfer`.

* `somalign_query_from_som()` — zero-reprojection query constructor that reuses
  `som$unit.classif` directly for per-cell node assignments, bypassing the
  O(N × nodes) nearest-code search that `somalign_query()` performs.  Accepts
  an optional pre-transformed `codebook` (e.g. after winsorisation and rescaling
  into reference-scaled space) and a `codebook_space` argument identical to
  `somalign_query()`.  The returned `somalign_query` object is fully compatible
  with `somalign_fit()`.

  `sample_distance` is set to `NA` (the field is not read by `somalign_fit()`).
  When the codebook has been non-linearly post-processed, the reused assignments
  are approximate: cells near node boundaries may flip, but this affects only a
  small fraction of tail cells and does not measurably impact OT alignment quality.

* `somalign_sensitivity_grid()` — direct projection cached across grid points.
  The O(N × reference\_nodes) projection of query cells onto the reference
  codebook is identical for every `(epsilon, rho_query, rho_ref)` combination,
  so it is now computed once before the grid loop and reused.  For a K-point
  grid this saves (K − 1) projection passes.  The function signature is otherwise
  unchanged (all `somalign_fit()` parameters are accepted explicitly; `...`
  forwarding has been removed).

# somalign 0.99.3

## New features

* `somalign_reference_from_som()` — zero-reprojection reference constructor
  for trained kohonen SOM objects.  Reuses all information already stored in
  the trained object: node masses from `som$unit.classif` (exact counts over
  every training cell at no computational cost), label probabilities from the
  supervised Y-layer codebook `codes[[2]]` (enables label transfer without
  passing any `labels` vector), and per-node distance thresholds recomputed
  in reference-scaled X-space from the embedded `som$data[[1]]`.  The
  distance computation is O(N × p) — a single subtraction per cell against
  its already-known assigned node — with no argmax and no O(N × nodes) memory
  peak.  On a 44.6 M-cell pilot cohort this eliminates the bottleneck that
  previously required subsampling to 100 k cells.

  The existing `somalign_reference()` and `somalign_reference_from_nodes()`
  APIs are unchanged.

  Note: on an `xyf`/`supersom` SOM trained with equal X+Y layer weights,
  `unit.classif` reflects a joint supervised assignment.  Node masses
  therefore match the SOM's own partition rather than a pure-X nearest-node
  assignment.  Distance quantiles are computed in X-only space so
  outside-reference thresholds remain on the same scale as
  `somalign_fit()`'s query distances.

# somalign 0.99.2

## Breaking changes

* Default `epsilon` lowered from `0.5` to `0.1` in `somalign_fit()`,
  `somalign_fit_anchored()`, and `somalign_som_stability()`. Once cost
  normalisation was added in v0.99.1, `epsilon = 0.5` left the transport plan
  too diffuse: on the Nuñez 2023 full-spectrum flow and CyTOF benchmarks, label
  posteriors collapsed onto the reference class prior and CyTOF acceptance
  dropped to zero. At `epsilon = 0.1` minority classes transfer correctly and
  the corrected JS divergence is 4× lower. Code that sets `epsilon` explicitly
  is unaffected; if you relied on the default, re-check your acceptance rates
  and correction norms.

* Default `epsilon_global` in `somalign_fit_two_pass()` lowered from `0.5` to
  `0.3`. The global pass aligns coarse structure and tolerates more smoothing
  than the local pass, so it keeps a larger `epsilon` than the single-pass
  default.

## New features

* `somalign_fit()` gains a `solver = "log_domain"` option. The log-domain
  Sinkhorn variant works with log-potentials $f$ (query) and $g$ (reference)
  rather than the primal scaling vectors, avoiding `exp(-C/epsilon)` kernel
  underflow entirely. It is slower per iteration but tolerates cost/epsilon
  ratios that cause the default `"internal"` solver to warn.

* `somalign_som_stability()` — new function that sweeps a vector of random
  seeds for query SOM training, holds the reference fixed, and returns a
  per-seed summary data frame. Reports `transport_mass`,
  `mean_match_fraction`, `max_row_mass_error`, `accepted_label_fraction`,
  `outside_direct_fraction`, `outside_corrected_fraction`,
  `mean_correction_norm`, and `converged` for each seed, quantifying
  run-to-run variance from SOM training randomness — the largest uncontrolled
  variance source in the `somalign` workflow.

## Improvements

* `somalign_fit()` diagnostics now include two scale-invariant marginal-error
  fields: `diagnostics$solver$rel_marginal_row_error` and
  `diagnostics$solver$rel_marginal_col_error`. These normalise the absolute
  marginal violations by total query and reference mass respectively, making
  it easier to distinguish expected mass destruction from solver
  non-convergence.

## Documentation

* `vignette("algorithm")` gains a new **Limitations of the OT correction**
  section covering barycentric contraction (the correction target is a
  conditional mean that contracts toward dense reference nodes and grows with
  `epsilon`), the unpaired OT design (mass distributions, not individual
  samples, are matched — anchor-based methods exploit sample-level pairing),
  and the uncalibrated novelty score (`match_fraction` reflects the
  epsilon/rho regime as much as any true novel population).

# somalign 0.99.1

## Breaking changes

* `somalign_fit()` default `epsilon` changed from `0.05` to `0.5`. The cost
  matrix is now normalised by its median positive entry before the Sinkhorn
  kernel is computed, making `epsilon` scale- and dimension-invariant across
  feature sets of different size and spread. Users who set `epsilon` explicitly
  relative to raw codebook distances should re-tune against the normalised scale
  (inspect `diagnostics$solver$cost_scale`).

## New features and improvements

* `somalign_query()` gains a `codebook_space` argument (default
  `"reference_scaled"`). When a user-supplied `som_query` is provided, setting
  `codebook_space = "raw"` instructs `somalign_query()` to re-scale the
  codebook into reference-scaled space, mirroring the existing behaviour of
  `somalign_reference()`.

* `somalign_fit()` diagnostics now include `diagnostics$solver$converged`
  (logical), `diagnostics$solver$final_delta` (final iterate change), and
  `diagnostics$solver$cost_scale` (cost normalisation factor).

## Bug fixes and hardening

* `somalign_fit()` now emits a warning when more than 50 % of query mass is
  destroyed during transport, or when more than 50 % of query samples fall
  outside reference distance thresholds.

* A non-finite Sinkhorn iterate delta now triggers an explicit warning instead
  of finishing silently (previously, a NaN delta suppressed both the
  early-break condition and the non-convergence warning).

* `somalign_reference_from_nodes()` now emits informational messages when no
  label probabilities or distance quantiles are supplied, so users are aware
  that label transfer and outside-reference detection are disabled.

* The internal OT input validator now warns when all query or reference node
  masses are zero, preventing a silent all-zero transport plan.

* Fixed a silent fallback in `.somalign_thresholds()`: when the requested
  distance-quantile column is absent, a warning is now emitted rather than
  silently substituting the third column by position.

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
