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
  parallel = FALSE,
  ...
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

- parallel:

  Logical. When `TRUE`, grid rows are evaluated in parallel using
  [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
  with the registered `BiocParallel` back-end (see
  [`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)).
  Configure the back-end before calling this function, e.g.
  `BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))`.
  When `FALSE` (default) a sequential for-loop is used, which is fully
  reproducible across platforms.

- ...:

  Additional arguments passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

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
somalign_sensitivity_grid(qry, ref,
                          epsilon = c(0.05, 0.1),
                          rho_query = c(0.5, 1),
                          rho_ref = 1)
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.31); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.38); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.12); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.18); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#>   epsilon rho_query rho_ref   solver transport_mass mean_match_fraction
#> 1    0.05       0.5       1 internal      0.8711764           0.8958259
#> 2    0.10       0.5       1 internal      0.9210089           0.9214855
#> 3    0.05       1.0       1 internal      0.8991056           0.9252915
#> 4    0.10       1.0       1 internal      0.9372579           0.9449376
#>   max_row_mass_error max_col_mass_error accepted_label_fraction
#> 1         0.09013118         0.08489199                       0
#> 2         0.08120687         0.06314924                       0
#> 3         0.06178872         0.07665801                       0
#> 4         0.05393120         0.05850111                       0
#>   outside_direct_fraction outside_corrected_fraction
#> 1                     0.3                        0.3
#> 2                     0.3                        0.1
#> 3                     0.3                        0.2
#> 4                     0.3                        0.1
```
