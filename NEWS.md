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
