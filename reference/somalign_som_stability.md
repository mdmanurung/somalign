# Assess alignment stability across query SOM random seeds

Trains a new query SOM for each seed and runs the full alignment
pipeline, holding the reference SOM fixed. The returned summary
quantifies how much OT alignment statistics vary with query SOM training
randomness — the largest uncontrolled variance source in the `somalign`
workflow.

## Usage

``` r
somalign_som_stability(
  query_data,
  reference,
  som_seeds = seq_len(5L),
  epsilon = 0.1,
  rho_query = 1,
  rho_ref = 1,
  grid = NULL,
  rlen = 100,
  alpha = c(0.05, 0.01),
  parallel = FALSE,
  ...
)
```

## Arguments

- query_data:

  Numeric query data.

- reference:

  A `somalign_reference` object.

- som_seeds:

  Integer vector of random seeds used to train query SOMs.

- epsilon:

  Entropic regularisation passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- rho_query, rho_ref:

  Mass-relaxation parameters passed to
  [`somalign_fit()`](https://mdmanurung.github.io/somalign/reference/somalign_fit.md).

- grid:

  Optional
  [`kohonen::somgrid()`](https://rdrr.io/pkg/kohonen/man/unit.distances.html)
  for query SOM training.

- rlen:

  Number of SOM training iterations.

- alpha:

  Learning-rate schedule.

- parallel:

  Logical. Use
  [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html)
  when `TRUE`.

- ...:

  Additional arguments passed to
  [`somalign_query()`](https://mdmanurung.github.io/somalign/reference/somalign_query.md).

## Value

A data frame with one row per seed containing key alignment summary
statistics: `som_seed`, `transport_mass`, `mean_match_fraction`,
`max_row_mass_error`, `accepted_label_fraction`,
`outside_direct_fraction`, `outside_corrected_fraction`,
`mean_correction_norm`, `converged`.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(40), nrow = 20, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
somalign_som_stability(mat, ref, som_seeds = 1:3,
                       grid = kohonen::somgrid(2, 2, "hexagonal"), rlen = 5)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.06); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.02); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#> somalign_fit: 2 query node(s) have match_mass_ratio > 1 (max 1.11); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
#>   som_seed transport_mass mean_match_fraction max_row_mass_error
#> 1        1       0.992981           0.9665548        0.033445192
#> 2        2       1.003795           0.9901143        0.006636695
#> 3        3       0.977416           0.9619554        0.044318134
#>   accepted_label_fraction outside_direct_fraction outside_corrected_fraction
#> 1                       0                     0.2                       0.05
#> 2                       0                     0.2                       0.00
#> 3                       0                     0.2                       0.10
#>   mean_correction_norm converged
#> 1            0.4730823      TRUE
#> 2            0.3180551      TRUE
#> 3            0.4679876      TRUE
```
