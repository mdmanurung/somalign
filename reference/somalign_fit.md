# Align a query SOM to a reference SOM

Align a query SOM to a reference SOM

## Usage

``` r
somalign_fit(
  query,
  reference,
  epsilon = 0.1,
  rho_query = 1,
  rho_ref = 1,
  solver = c("internal", "log_domain", "auto", "annealing"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L,
  diagonal_boost = 0,
  label_guided = FALSE,
  anneal_start = 10,
  anneal_stages = 10L,
  anneal_factor = NULL,
  feature_weights = NULL,
  laplacian_lambda = 0
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Entropic regularisation strength. The cost matrix is normalised by its
  median positive entry before computing the Sinkhorn kernel, so
  `epsilon` is approximately scale- and dimension-invariant. The default
  `0.1` gives a sharp transport plan that preserves cell-type
  specificity for typical z-scored SOM codebooks. Larger values
  (0.3–0.5) produce smoother, more diffuse plans that can help
  convergence on noisy or high-dimensional data but dilute label
  posteriors and increase barycentric shrinkage in the corrected
  projection. Very small values (\< 0.05) make the transport
  increasingly discrete and may require `solver = "log_domain"` for
  numerical stability. The normalisation scale is stored in
  `diagnostics$solver$cost_scale`.

- rho_query:

  Query-side unbalanced mass relaxation.

- rho_ref:

  Reference-side unbalanced mass relaxation.

- solver:

  Sinkhorn solver variant. `"internal"` (default) and `"auto"` both use
  the primal-domain scaling iteration. `"log_domain"` uses a numerically
  stable log-potential variant that avoids kernel underflow for small
  `epsilon` or high-dimensional codebooks; it is slower per iteration
  but tolerates cost/epsilon ratios that cause `"internal"` to warn.
  `"annealing"` runs the log-domain solver across a geometric epsilon
  cooling schedule (starting at `anneal_start * epsilon`, cooling to
  `epsilon` over `anneal_stages` stages), warm-starting each stage from
  the previous stage's dual potentials. Recommended for `label_guided`
  fits or any fit with small `epsilon` (\< 0.05) where cold-start
  Sinkhorn is slow or non-convergent; never underflows, since it never
  exponentiates the kernel.

- min_match_fraction:

  Minimum transported fraction required before a query node label
  transfer is accepted.

- confidence_threshold:

  Minimum top-label probability required before a query node label
  transfer is accepted.

- correction_min_mass:

  Minimum transported node mass required before a correction shift is
  applied. Corrections also require the node match fraction to pass
  `min_match_fraction`.

- max_iter:

  Maximum internal Sinkhorn iterations.

- tol:

  Internal Sinkhorn convergence tolerance.

- chunk_size:

  Integer. Number of samples to project per chunk when computing nearest
  reference node. Use `Inf` or `NULL` for no chunking (allocates a full
  n_samples x n_nodes matrix). Default `10000L`.

- diagonal_boost:

  Non-negative scalar. Amount by which to reduce the normalised OT cost
  for each query node's nearest reference node. A positive value makes
  the transport plan prefer identity-like mappings, shrinking
  over-correction when the two codebooks are already close. Zero
  (default) leaves the cost unchanged. Values around 0.1–0.5 are a
  reasonable starting point; very large values concentrate all mass on
  the diagonal and the plan degrades toward simple nearest-neighbour
  assignment.

- label_guided:

  Logical. When `TRUE`, uses `query$label_prob` and
  `reference$label_prob` to add a large cost penalty for node pairs
  whose dominant labels disagree, constraining OT to transport mass
  predominantly between concordant cell-type nodes. Nodes where the
  maximum label probability is below 0.5 are treated as unlabeled and
  are never penalized. Errors if `label_guided = TRUE` but either
  `label_prob` is `NULL`.

- anneal_start:

  Positive scalar \>= 1. When `solver = "annealing"`, the starting
  epsilon is `anneal_start * epsilon`. Default `10`. Ignored when
  `solver != "annealing"`.

- anneal_stages:

  Positive integer. Number of cooling stages in the annealing schedule,
  including the final stage at the target `epsilon`. Default `10L`. A
  value of `1` degenerates to a cold-start log-domain solve. Ignored
  when `solver != "annealing"`.

- anneal_factor:

  Positive scalar \< 1, or `NULL` (default). When not `NULL`, overrides
  the auto-computed per-stage cooling ratio. Ignored when
  `solver != "annealing"`.

- feature_weights:

  Either `NULL` (default, squared-Euclidean cost) or a named
  non-negative numeric vector with one entry per feature (explicit
  diagonal Mahalanobis weights on the OT cost). Weights are applied as
  `sqrt(w_f)` per-column scaling of both codebooks before the squared
  Euclidean distance is computed, yielding cost \\\sum_f w_f (q\_{if} -
  r\_{jf})^2\\. The resolved vector is stored in
  `fit$diagnostics$cost_metric$feature_weights`. Projection and
  threshold distances
  ([`somalign_results()`](https://mdmanurung.github.io/somalign/reference/somalign_results.md))
  are unaffected – weighting applies only to the OT cost. See
  [`somalign_fit_anchored()`](https://mdmanurung.github.io/somalign/reference/somalign_fit_anchored.md)
  for `"anchor"`, which auto-estimates weights from anchor
  displacements.

- laplacian_lambda:

  Non-negative scalar. Graph-Laplacian regularisation strength for the
  node-shift field. When greater than zero, the M x p raw node shifts
  are smoothed by solving \\(W + \lambda L)\\s^\* = W\\s\\, where \\W =
  \mathrm{diag}(\text{node\\masses})\\ (with
  `correction_allowed == FALSE` nodes zeroed out) and \\L\\ is the graph
  Laplacian of the query SOM's hexagonal or rectangular neighbor graph.
  This penalises squared differences between adjacent-node shifts,
  producing a spatially coherent correction field instead of one where
  neighboring nodes can receive wildly different shifts from
  finite-sample OT noise. Default `0` (no smoothing, exact current
  behaviour). A natural starting range is `0.1`–`1.0` (same
  cost/squared-distance scale as `epsilon`); larger values increasingly
  collapse the field toward its mass-weighted mean. Requires the query
  SOM to carry 2-D grid coordinates (`query$som_query$grid$pts`, present
  for any
  [`kohonen::som()`](https://rdrr.io/pkg/kohonen/man/supersom.html)- or
  [`kohonen::supersom()`](https://rdrr.io/pkg/kohonen/man/supersom.html)-trained
  SOM); errors otherwise.

## Value

A `somalign_fit` object.

## Details

The transport plan row sums will not equal `query$node_masses` exactly –
this is by design. Unbalanced optimal transport allows mass destruction,
so some query mass may be absorbed rather than transported. Deviation
grows with lower `rho_query` / `rho_ref` values and higher `epsilon`.
Use `diagnostics$ot$max_row_mass_error` to quantify the deviation in a
given fit; for near-balanced data, increase `rho_query` (e.g.
`rho_query = 10`) to enforce tighter marginal constraints. A warning is
emitted automatically when more than 50% of query mass is destroyed.

The cost matrix is normalised by its median positive entry before the
Sinkhorn kernel is computed. This makes `epsilon` scale- and
dimension-invariant: the same value produces the same degree of
regularisation regardless of the number of features or the spread of
codebook coordinates. The raw normalisation factor is stored as
`diagnostics$solver$cost_scale`.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
fit <- somalign_fit(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
```
