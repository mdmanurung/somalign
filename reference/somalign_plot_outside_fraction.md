# Plot fraction of cells outside reference thresholds

Compares the percentage of cells that fall outside the reference
distance threshold before (Direct) and after (Corrected) applying the OT
correction vectors.

## Usage

``` r
somalign_plot_outside_fraction(fit)
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
#> somalign_reference_from_som: SOM has no second code layer; label transfer will be disabled.
fit <- somalign_fit(qry, ref)
#> somalign_fit: 3 query node(s) have match_mass_ratio > 1 (max 1.23); this is expected in unbalanced OT. See diagnostics$ot$match_mass_ratio for details.
somalign_plot_outside_fraction(fit)
```
