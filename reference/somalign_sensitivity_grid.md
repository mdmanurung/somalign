# Run an OT sensitivity grid

Run an OT sensitivity grid

## Usage

``` r
somalign_sensitivity_grid(
  query,
  reference,
  epsilon,
  rho_query,
  rho_ref,
  solver = c("internal", "auto"),
  min_match_fraction = 0.05,
  confidence_threshold = 0.6,
  correction_min_mass = 1e-08,
  max_iter = 1000,
  tol = 1e-07,
  chunk_size = 10000L,
  diagonal_boost = 0,
  parallel = FALSE
)
```

## Arguments

- query:

  A `somalign_query` object.

- reference:

  A `somalign_reference` object.

- epsilon:

  Numeric vector of entropic regularisation values.

- rho_query:

  Numeric vector of query-side mass relaxation values.

- rho_ref:

  Numeric vector of reference-side mass relaxation values.

- solver:

  Solver passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).
  `"auto"` is accepted as a compatibility alias for the internal pure-R
  solver.

- min_match_fraction:

  Minimum match fraction threshold passed to each
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md)
  call. Default `0.05`.

- confidence_threshold:

  Minimum label confidence for accepted label transfer. Default `0.6`.

- correction_min_mass:

  Minimum OT mass for a node shift to be applied. Default `1e-8`.

- max_iter:

  Maximum Sinkhorn iterations. Default `1000`.

- tol:

  Sinkhorn convergence tolerance. Default `1e-7`.

- chunk_size:

  Integer. Number of samples per projection chunk. `NULL` processes all
  samples at once. Default `10000L`.

- diagonal_boost:

  Non-negative scalar added to same-node OT costs to discourage
  self-transport. Default `0`.

- parallel:

  Logical. When `TRUE`, grid rows are evaluated in parallel using
  [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
  with the registered `BiocParallel` back-end (see
  [`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)).
  Configure the back-end before calling this function, e.g.
  `BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))`.
  When `FALSE` (default) a sequential for-loop is used, which is fully
  reproducible across platforms.

## Value

A data frame with one row per parameter combination.

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
somalign_sensitivity_grid(qry, ref,
                          epsilon = c(0.05, 0.1),
                          rho_query = c(0.5, 1),
                          rho_ref = 1)
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.34); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.44); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.17); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#>   epsilon rho_query rho_ref   solver transport_mass mean_match_fraction
#> 1    0.05       0.5       1 internal      0.9778797           0.9517853
#> 2    0.10       0.5       1 internal      1.0304219           0.9788873
#> 3    0.05       1.0       1 internal      0.9818532           0.9688494
#> 4    0.10       1.0       1 internal      1.0217717           0.9890859
#>   max_row_mass_error max_col_mass_error accepted_label_fraction
#> 1         0.04917404         0.03447202                       0
#> 2         0.04375593         0.03647086                       0
#> 3         0.03268400         0.03805557                       0
#> 4         0.02342453         0.03828297                       0
#>   outside_direct_fraction outside_corrected_fraction
#> 1                     0.3                        0.2
#> 2                     0.3                        0.1
#> 3                     0.3                        0.1
#> 4                     0.3                        0.1
```
