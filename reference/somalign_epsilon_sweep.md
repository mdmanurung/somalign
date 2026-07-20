# Epsilon phase-transition sweep for principled epsilon selection

Runs the Sinkhorn OT solve across a log-spaced epsilon grid without the
per-cell projection step, computing the transport-plan order parameter
(mean fractional reference-node usage, a perplexity-based measure), its
susceptibility (rate of change with log epsilon), the dual free energy,
and the mutual information between query and reference nodes at each
grid point.

## Usage

``` r
somalign_epsilon_sweep(
  query,
  reference,
  epsilon_grid = NULL,
  n_grid = 25L,
  rho_query = 1,
  rho_ref = 1,
  solver = c("log_domain", "internal", "auto", "annealing"),
  max_iter = 1000,
  tol = 1e-07,
  diagonal_boost = 0,
  label_guided = FALSE,
  parallel = FALSE,
  anneal_start = 10,
  anneal_factor = NULL,
  anneal_stages = 10L,
  topology = FALSE
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon_grid:

  Numeric vector of epsilon values. If `NULL` (default), a log-spaced
  grid of `n_grid` values from 1e-3 to 5 is used.

- n_grid:

  Integer. Grid size when `epsilon_grid = NULL`. Default `25`.

- rho_query, rho_ref:

  Mass relaxation parameters passed to the OT solver.

- solver:

  Sinkhorn solver. Default `"log_domain"` (required for `log_Z`;
  `"internal"`/`"auto"` fill `log_Z` with `NA`). `"annealing"` runs the
  geometric epsilon-cooling schedule at each grid point (using
  `anneal_start`/`anneal_factor`/`anneal_stages`) and also reports
  `log_Z`.

- max_iter, tol:

  Sinkhorn convergence parameters.

- diagonal_boost:

  Non-negative cost reduction on nearest-reference-node entries. Default
  `0`.

- label_guided:

  Logical; see
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- parallel:

  Logical; see
  [`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md).

- anneal_start, anneal_factor, anneal_stages:

  Annealing-schedule tuning parameters, used only when
  `solver = "annealing"`. See
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- topology:

  Logical. When `TRUE`, append three topology columns to `table`:
  `n_components_query` (H0 component count of the unchanged query
  codebook, computed once), `n_components_corrected` (H0 component count
  of the corrected query codebook at each epsilon), and
  `biggest_merge_mass_frac` (fraction of total query node mass in the
  single largest H0 component of the corrected codebook at the
  auto-selected threshold – the key over-merging signal). The threshold
  is derived via
  [`somalign_topology_audit()`](https://mdmanurung.github.io/somalign/reference/somalign_topology_audit.md)
  conventions (median of `reference$distance_quantiles` 95th-percentile
  column). Default `FALSE`. Note: these columns always use *all* query
  nodes, equivalent to `somalign_topology_audit(fit, nodes = "all")`,
  because the correction-allowed node set is itself epsilon-dependent
  and subsetting by it would confound the component-count trend across
  the grid. The default `somalign_topology_audit(fit)` uses
  `nodes = "correction_allowed"`, so its `n_components_corrected` may
  differ from the sweep's; pass `nodes = "all"` to the audit for
  directly comparable numbers.

## Value

A list of class `"somalign_epsilon_sweep"` with components `table`,
`epsilon_c`, `epsilon_rec`, `cost_scale`. The `table` data frame has one
row per epsilon value with columns `epsilon`, `log_epsilon`, `Phi`,
`susceptibility`, `log_Z`, `mutual_information`,
`conditional_entropy_mean`, `expected_cost`, `transport_mass`,
`iterations`, `converged`. When `topology = TRUE`, three additional
columns are appended: `n_components_query`, `n_components_corrected`,
`biggest_merge_mass_frac`.

## Details

The order parameter \\\Phi(\epsilon)\\ is the mean effective fraction of
reference nodes used per query node: \\\Phi = \mathrm{mean}\_i(2^{H_i})
/ K\\, where \\H_i\\ is the (bit) conditional entropy of the
row-normalised transport plan for query node \\i\\ and \\K\\ is the
number of reference nodes. As \\\epsilon \to 0\\, \\\Phi \to 1/K\\; as
\\\epsilon \to \infty\\, \\\Phi \to 1\\. The susceptibility \\\chi =
d\Phi/d\log\epsilon\\ peaks at the critical epsilon \\\epsilon_c\\,
marking the crossover between a localised (transport-cost-dominated) and
delocalised (entropy-dominated) plan. `epsilon_rec` (\\0.3 \\
\epsilon_c\\) keeps the plan in the ordered phase with a safety margin.

The sweep avoids the per-cell projection step, so it runs approximately
one OT solve per epsilon value – much faster than
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md)
for the same epsilon range.

## See also

[`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md),
[`somalign_select_epsilon()`](https://mdmanurung.github.io/somalign/reference/somalign_select_epsilon.md),
[`somalign_sensitivity_grid()`](https://mdmanurung.github.io/somalign/reference/somalign_sensitivity_grid.md),
[`somalign_diagnostics()`](https://mdmanurung.github.io/somalign/reference/somalign_diagnostics.md)

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat + 0.5, ref,
                      grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
sw <- somalign_epsilon_sweep(qry, ref, n_grid = 8)
#> Warning: Sinkhorn solver did not converge after 1000 iterations (final delta = 4.214e-05). Consider increasing max_iter, raising epsilon, or reducing rho_query / rho_ref.
#> Warning: Sinkhorn solver did not converge after 1000 iterations (final delta = 1.220e-06). Consider increasing max_iter, raising epsilon, or reducing rho_query / rho_ref.
sw$epsilon_c
#> [1] 0.4386533
sw$epsilon_rec
#> [1] 0.131596
```
