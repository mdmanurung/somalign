# Changelog

## somalign 0.99.2

### New features

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  gains a `solver = "log_domain"` option. The log-domain Sinkhorn
  variant works with log-potentials $`f`$ (query) and $`g`$ (reference)
  rather than the primal scaling vectors, avoiding `exp(-C/epsilon)`
  kernel underflow entirely. It is slower per iteration but tolerates
  cost/epsilon ratios that cause the default `"internal"` solver to
  warn.

- [`somalign_som_stability()`](https://mdmanurung.github.io/somalign/reference/somalign_som_stability.md)
  — new function that sweeps a vector of random seeds for query SOM
  training, holds the reference fixed, and returns a per-seed summary
  data frame. Reports `transport_mass`, `mean_match_fraction`,
  `max_row_mass_error`, `accepted_label_fraction`,
  `outside_direct_fraction`, `outside_corrected_fraction`,
  `mean_correction_norm`, and `converged` for each seed, quantifying
  run-to-run variance from SOM training randomness — the largest
  uncontrolled variance source in the `somalign` workflow.

### Improvements

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  diagnostics now include two scale-invariant marginal-error fields:
  `diagnostics$solver$rel_marginal_row_error` and
  `diagnostics$solver$rel_marginal_col_error`. These normalise the
  absolute marginal violations by total query and reference mass
  respectively, making it easier to distinguish expected mass
  destruction from solver non-convergence.

### Documentation

- [`vignette("algorithm")`](https://mdmanurung.github.io/somalign/articles/algorithm.md)
  gains a new **Limitations of the OT correction** section covering
  barycentric contraction (the correction target is a conditional mean
  that contracts toward dense reference nodes and grows with `epsilon`),
  the unpaired OT design (mass distributions, not individual samples,
  are matched — anchor-based methods exploit sample-level pairing), and
  the uncalibrated novelty score (`match_fraction` reflects the
  epsilon/rho regime as much as any true novel population).

## somalign 0.99.1

### Breaking changes

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  default `epsilon` changed from `0.05` to `0.5`. The cost matrix is now
  normalised by its median positive entry before the Sinkhorn kernel is
  computed, making `epsilon` scale- and dimension-invariant across
  feature sets of different size and spread. Users who set `epsilon`
  explicitly relative to raw codebook distances should re-tune against
  the normalised scale (inspect `diagnostics$solver$cost_scale`).

### New features and improvements

- [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  gains a `codebook_space` argument (default `"reference_scaled"`). When
  a user-supplied `som_query` is provided, setting
  `codebook_space = "raw"` instructs
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  to re-scale the codebook into reference-scaled space, mirroring the
  existing behaviour of
  [`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md).

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  diagnostics now include `diagnostics$solver$converged` (logical),
  `diagnostics$solver$final_delta` (final iterate change), and
  `diagnostics$solver$cost_scale` (cost normalisation factor).

### Bug fixes and hardening

- [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  now emits a warning when more than 50 % of query mass is destroyed
  during transport, or when more than 50 % of query samples fall outside
  reference distance thresholds.

- A non-finite Sinkhorn iterate delta now triggers an explicit warning
  instead of finishing silently (previously, a NaN delta suppressed both
  the early-break condition and the non-convergence warning).

- [`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md)
  now emits informational messages when no label probabilities or
  distance quantiles are supplied, so users are aware that label
  transfer and outside-reference detection are disabled.

- The internal OT input validator now warns when all query or reference
  node masses are zero, preventing a silent all-zero transport plan.

- Fixed a silent fallback in `.somalign_thresholds()`: when the
  requested distance-quantile column is absent, a warning is now emitted
  rather than silently substituting the third column by position.

## somalign 0.99.0

- Added
  [`somalign_train_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_train_reference.md)
  and
  [`somalign_reference()`](https://mdmanurung.github.io/somalign/reference/somalign_reference.md)
  to train or wrap an existing `kohonen` SOM as a fixed reference
  object, including node mass, per-node label probability, and per-node
  distance quantile summaries.

- Added
  [`somalign_reference_from_nodes()`](https://mdmanurung.github.io/somalign/reference/somalign_reference_from_nodes.md)
  to reconstruct a reference object from pre-computed node-level
  artifacts (codebook, center, scale, masses, label probabilities,
  distance quantiles) without needing the original data.

- Added
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md)
  to scale query data against reference parameters, optionally train a
  query SOM, and return a `somalign_query` object containing per-sample
  node assignments and distances.

- Added
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  to align a query SOM to the fixed reference SOM via codebook-level
  unbalanced entropic optimal transport. Uses an internal pure-R
  generalized Sinkhorn solver; `"auto"` is retained as a compatibility
  alias for the internal solver.

- Added
  [`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md)
  to extract a per-sample data frame with direct reference projection
  columns (canonical) and transport-corrected projection columns
  (auxiliary, for annotation and visualisation).

- Added
  [`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)
  to return solver, OT plan, node-level, and projection diagnostics from
  a `somalign_fit` object.

- Added
  [`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
  to sweep a grid of OT hyperparameters (`epsilon`, `rho_query`,
  `rho_ref`) and return a summary data frame for tuning and stability
  assessment. Supports optional parallelisation via
  [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
  when `parallel = TRUE`.
