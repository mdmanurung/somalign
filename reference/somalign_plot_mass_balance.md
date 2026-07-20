# Plot node mass balance

Scatter plot of query node mass vs transported mass, coloured by match
fraction. Points lying on the diagonal received all their mass; points
below it had mass destroyed by the unbalanced OT solver.

## Usage

``` r
somalign_plot_mass_balance(fit)
```

## Arguments

- fit:

  A `somalign_fit` object.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
              dimnames = list(NULL, c("F1", "F2")))
ref <- somalign_train_reference(mat, grid = kohonen::somgrid(2, 2, "hexagonal"),
                                rlen = 5)
qry <- somalign_query(mat, ref, grid = kohonen::somgrid(2, 2, "hexagonal"),
                      rlen = 5)
fit <- somalign_fit(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_plot_mass_balance(fit)
```
